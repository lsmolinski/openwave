"""
LAGRANGIAN-FIELD ENGINE (M5)

PDE-driven evolution of the field ψ on a voxel grid. The action is

    S[ψ] = ∫ dt d³x [ ½|ψ̇|² − ½c²|∇ψ|² − V(ψ) ]

Euler-Lagrange yields  ∂²ψ/∂t² = c²·∇²ψ − ∂V/∂ψ.  In M5.0d the potential is
zero (free wave); M5.2 plugs in Klein-Gordon mass + Close Eq. 23 + LdG.

Module layout (top-to-bottom):
    INITIAL-CONDITION SEEDING — seed_gaussian, seed_dispersion_modes,
                                seed_vacuum, seed_hedgehog
    DIFFERENTIAL OPERATORS    — Laplacian, divergence, curl, curl-curl
    ψ EVOLVING ENGINE         — evolve_psi (leapfrog/Verlet)
    WAVE TRACKERS             — update_trackers_psi (amp / freq EMA — writes
                                to Trackers; stateful, dynamics-gated)
    FIELD OBSERVABLES         — compute_energyH_density (Hamiltonian),
                                compute_energyF_density (Frank elastic);
                                both write to FieldObservables (stateless,
                                always-on per-frame measurements)
    GRADIENT-DESCENT RELAX    — relax_director_step (one step of gradient
                                descent on Frank energy w/ tangent projection
                                + unit-length renorm + soft core pin; M5.1
                                task 6 — port from Exp 2)
    WINDING-NUMBER DIAGNOSTIC — compute_winding_number (Q on a sphere
                                around a defect; CPU numpy; M5.1 task 8 —
                                port from Exp 3; diagnostic only)
    POSITION RENDER           — sample_position_to_render (granule viz)
    3-PLANE SAMPLING          — sample_avg_trackers   (amp / freq globals)
                                sample_avg_observables (energyH / energyF
                                globals; separate pass for SoC, 2026-05-11)
    FLUX MESH UPDATING        — update_flux_mesh_values (color mapping;
                                reads both Trackers and FieldObservables)
    DIRECTOR GLYPHS           — update_director_glyphs (M5.1 director viz)
"""

import math

import taichi as ti
import numpy as np

from openwave.common import colormap

# ================================================================
# INITIAL-CONDITION SEEDING
# ================================================================
# Seed kernels prime the triple-buffer (psi_prev, psi at t−dt and t) with a
# known analytical solution so evolve_psi can step it forward. Used by:
#   - test/UI xperiments to verify the leapfrog kernel visually
#   - M5.0h dispersion-relation regression test


@ti.kernel
def seed_gaussian(
    wave_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    amplitude_am: ti.f32,  # type: ignore
    voxels_per_wavelength: ti.f32,  # type: ignore
    polarization: ti.template(),  # type: ignore
    direction: ti.template(),  # type: ignore
):
    """
    Seed psi_am (t=0) and psi_prev_am (t=−dt) with a Gaussian-windowed
    wave packet that respects the Dirichlet BC by construction.

        ψ(x, t) = A · sin(k·x − ω·t) · exp(−((x − x_c) / σ)²) · ê_pol

    The Gaussian envelope smoothly tapers ψ to ~0 at the boundaries, so we
    don't inject a step discontinuity at the edges that would radiate garbage
    inward and blow up amplitudes (an unwindowed sin(kx) has ψ ≠ 0 wherever
    the boundary doesn't land on a zero of sin — i.e. almost always).

    σ is chosen as N/6 voxels (where N is the propagation-axis grid size) so
    the envelope ≈ exp(−9) ≈ 1.2e−4 at the boundaries — effectively zero.
    Setting psi_prev at t=−dt (with both phase advance and shifted envelope
    center) gives the leapfrog forward-traveling velocity.

    Args:
        wave_field: WaveField instance (writes psi_am, psi_prev_am)
        c_amrs: wave speed (am/rs) — from launcher.compute_timestep
        dt_rs: timestep (rs) — from launcher.compute_timestep
        amplitude_am: peak amplitude in attometers
        voxels_per_wavelength: spatial period in voxel units (≥12 for stable f32)
        polarization: ti.types.vector(3, ti.f32) — displacement direction
        direction: ti.i32 — propagation axis: 0=x, 1=y, 2=z
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    k_grid = 2.0 * ti.math.pi / voxels_per_wavelength  # rad / voxel
    omega_rs = c_amrs * k_grid / wave_field.dx_am  # rad / rs (continuum ω = c·k)

    # 3-axis Gaussian envelope so ψ vanishes on EVERY boundary (Dirichlet BC),
    # not just the boundaries normal to the propagation direction. Without the
    # transverse taper, ψ would be translation-invariant in y/z and the seed
    # would have a step discontinuity at the y/z faces — each step those faces
    # radiate high-frequency garbage inward.
    x_c = ti.cast(nx, ti.f32) * 0.5
    y_c = ti.cast(ny, ti.f32) * 0.5
    z_c = ti.cast(nz, ti.f32) * 0.5
    sigma_x = ti.cast(nx, ti.f32) / 6.0
    sigma_y = ti.cast(ny, ti.f32) / 6.0
    sigma_z = ti.cast(nz, ti.f32) / 6.0

    # Phase shift from one timestep advancing the wave by c·dt (in voxel units)
    phase_dt = omega_rs * dt_rs
    # Envelope center moves backward by c·dt for the t=−dt buffer (group velocity = phase velocity for c·k dispersion)
    voxel_step = c_amrs * dt_rs / wave_field.dx_am

    for i, j, k in ti.ndrange(nx, ny, nz):
        x = ti.cast(i, ti.f32)
        y = ti.cast(j, ti.f32)
        z = ti.cast(k, ti.f32)

        # Spatial coordinate along propagation axis (drives the carrier phase)
        x_axis = x
        if direction == 1:
            x_axis = y
        elif direction == 2:
            x_axis = z

        phase = k_grid * x_axis
        envelope_3d = (
            ti.exp(-(((x - x_c) / sigma_x) ** 2))
            * ti.exp(-(((y - y_c) / sigma_y) ** 2))
            * ti.exp(-(((z - z_c) / sigma_z) ** 2))
        )

        # ψ at t=0: A · sin(k·x_axis) · envelope_3d(r)
        wave_field.psi_am[i, j, k] = amplitude_am * ti.sin(phase) * envelope_3d * polarization

        # ψ at t=−dt: carrier phase advances by ω·dt; envelope center shifts
        # backward by voxel_step along the propagation axis (transverse axes
        # don't move — the wave isn't translating in y/z).
        x_shift = x
        y_shift = y
        z_shift = z
        if direction == 0:
            x_shift = x + voxel_step
        elif direction == 1:
            y_shift = y + voxel_step
        else:
            z_shift = z + voxel_step

        envelope_prev = (
            ti.exp(-(((x_shift - x_c) / sigma_x) ** 2))
            * ti.exp(-(((y_shift - y_c) / sigma_y) ** 2))
            * ti.exp(-(((z_shift - z_c) / sigma_z) ** 2))
        )
        wave_field.psi_prev_am[i, j, k] = (
            amplitude_am * ti.sin(phase + phase_dt) * envelope_prev * polarization
        )

        # Initialize psi_new_am to the t=0 value too — required for BC consistency
        # so the first swap_buffers doesn't clobber boundary voxels (psi_new_am's
        # default is zero; without this write, psi_am's boundary becomes 0 after
        # the first swap regardless of what the seed put there). See evolve_psi
        # docstring "BC consistency requirement" for the full rationale.
        wave_field.psi_new_am[i, j, k] = amplitude_am * ti.sin(phase) * envelope_3d * polarization


@ti.kernel
def seed_dispersion_modes(
    wave_field: ti.template(),  # type: ignore
    mode_indices: ti.template(),  # type: ignore
    mode_amplitudes: ti.template(),  # type: ignore
    mode_cos_omega_dt: ti.template(),  # type: ignore
    polarization: ti.template(),  # type: ignore
):
    """
    Seed psi_am (t=0) and psi_prev_am (t=−dt) with a superposition of
    standing-wave Dirichlet eigenmodes — the M5.0h dispersion-relation
    regression test.

    Each mode m ∈ [0, n_modes) is a 3D Dirichlet eigenmode of the discrete
    Laplacian on the cell-centered grid:

        shape_m(i,j,k) = sin(n_x π i/(N−1)) · sin(n_y π j/(N−1)) · sin(n_z π k/(N−1))

    where (n_x, n_y, n_z) = mode_indices[m]. With ψ=0 at i ∈ {0, N−1}
    (and same for j, k), these are exact eigenmodes of the 6-point
    Laplacian with eigenvalue
        λ_m = −(c²/dx²) · 2 · Σ_α (1 − cos(n_α π / (N−1)))
    so each mode oscillates at its own ω_m given by the discrete leapfrog
    dispersion relation (validated by this very test).

    For a *standing* wave with zero initial velocity:
        ψ(x, t) = Σ_m A_m · cos(ω_m t) · shape_m(x) · ê_pol
        ψ(x, 0)   = Σ_m A_m · shape_m · ê_pol
        ψ(x, −dt) = Σ_m A_m · cos(ω_m·dt) · shape_m · ê_pol

    Caller pre-computes cos(ω_m·dt) per mode and passes via mode_cos_omega_dt.

    Args:
        wave_field: WaveField (writes psi_am, psi_prev_am)
        mode_indices: ti.field shape (n_modes, 3) of i32 — (n_x, n_y, n_z) per mode
        mode_amplitudes: ti.field shape (n_modes,) of f32 — A_m peak amplitude (am)
        mode_cos_omega_dt: ti.field shape (n_modes,) of f32 — cos(ω_m·dt) for ψ_prev
        polarization: ti.types.vector(3, ti.f32) — displacement direction ê
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    nm1_x = ti.cast(nx - 1, ti.f32)
    nm1_y = ti.cast(ny - 1, ti.f32)
    nm1_z = ti.cast(nz - 1, ti.f32)
    n_modes = mode_amplitudes.shape[0]

    for i, j, k in ti.ndrange(nx, ny, nz):
        sum_t0 = ti.cast(0.0, ti.f32)
        sum_tprev = ti.cast(0.0, ti.f32)
        x = ti.cast(i, ti.f32)
        y = ti.cast(j, ti.f32)
        z = ti.cast(k, ti.f32)
        for m in range(n_modes):
            n_x_m = ti.cast(mode_indices[m, 0], ti.f32)
            n_y_m = ti.cast(mode_indices[m, 1], ti.f32)
            n_z_m = ti.cast(mode_indices[m, 2], ti.f32)
            shape = (
                ti.sin(n_x_m * ti.math.pi * x / nm1_x)
                * ti.sin(n_y_m * ti.math.pi * y / nm1_y)
                * ti.sin(n_z_m * ti.math.pi * z / nm1_z)
            )
            amp = mode_amplitudes[m]
            sum_t0 += amp * shape
            sum_tprev += amp * mode_cos_omega_dt[m] * shape
        wave_field.psi_am[i, j, k] = sum_t0 * polarization
        wave_field.psi_prev_am[i, j, k] = sum_tprev * polarization
        # BC consistency: write psi_new_am too. See evolve_psi docstring.
        wave_field.psi_new_am[i, j, k] = sum_t0 * polarization


@ti.kernel
def seed_vacuum(
    wave_field: ti.template(),  # type: ignore
):
    """
    Fill ψ with the topological-vacuum ground state n = ẑ at every voxel.

    In M5.1+ the same ψ field plays a director role: a unit-length Vector(3)
    whose direction at each voxel encodes the local orientation of the medium.
    Vacuum = uniform orientation = no defects. This kernel writes the
    canonical vacuum (n = ẑ) into both the current and previous buffers,
    giving a static initial condition with zero implied velocity if any
    subsequent leapfrog step is ever taken.

    The vacuum is the substrate on which topological defects (particles)
    will be seeded by `seed_hedgehog` — far from any defect, the field returns
    to this configuration via the w_vac blend in that kernel.

    Args:
        wave_field: WaveField (writes psi_am, psi_prev_am)
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    z_hat = ti.Vector([0.0, 0.0, 1.0])
    for i, j, k in ti.ndrange(nx, ny, nz):
        wave_field.psi_am[i, j, k] = z_hat
        wave_field.psi_prev_am[i, j, k] = z_hat
        # BC consistency: write psi_new_am too so the first swap_buffers
        # doesn't clobber the boundary. See evolve_psi docstring.
        wave_field.psi_new_am[i, j, k] = z_hat


@ti.kernel
def seed_hedgehog(
    wave_field: ti.template(),  # type: ignore
    centers_voxel: ti.template(),  # type: ignore
    signs: ti.template(),  # type: ignore
    domain_quarter_voxels: ti.f32,  # type: ignore
    n_defects: ti.i32,  # type: ignore
):
    """
    Seed N topological hedgehog defects via weighted superposition + renormalization.

    Direct port of Exp 2's `seed_hedgehog_pair` (validated in
    `research_hub/sandbox_phase3_lagrangian/exp2_hedgehog_energy.py:71-108`),
    generalized from a hardcoded pair to an N-defect array. Same math; same
    weighting choices; same soft-core handling.

    Per voxel (x = (i, j, k)):
        1. For each defect d at center c[d] with sign s[d] ∈ {±1}:
              r[d]      = |x − c[d]|                            distance
              r_safe    = max(r[d], 0.2 · dx)                   soft-core
              n_radial  = s[d] · (x − c[d]) / r_safe            unit radial
              w_defect  = 1 / (r[d] + 0.5)                      proximity weight (~1/r tail)
              n_combined += w_defect · n_radial
              r_nearest = min(r_nearest, r[d])
        2. Vacuum-blend weight (uses nearest-defect distance):
              w_vac = 1 / (1 + (r_nearest / D/4)⁴)
            ≈ 1 inside ~D/4 of any defect; → 0 far away (returns to vacuum)
        3. Final blend: n = w_vac · n_combined + (1 − w_vac) · ẑ
        4. Renormalize: n / |n|  (with epsilon guard for the singular-defect
           voxel where the radial directions cancel)

    The seeded field is written to BOTH psi_am and psi_prev_am so a leapfrog
    step run on this initial condition starts with zero velocity. M5.1's
    primary use is static gradient-descent relaxation (no leapfrog); this
    is for forward-compatibility with M5.2+ dynamic tests.

    Args:
        wave_field: WaveField (writes psi_am, psi_prev_am)
        centers_voxel: ti.field shape (n_defects, 3) of i32 — centers in voxel coords
        signs: ti.field shape (n_defects,) of i32 — +1 outward / −1 inward per defect
        domain_quarter_voxels: f32 — D/4 in voxel units (sets w_vac falloff radius)
        n_defects: i32 — number of active defects in the centers/signs arrays

    Caveats:
        - The 0.2·dx soft-core is a numerical workaround for the topological
          singularity at the defect center; M5.5 (Skyrme stabilizer) and
          M5.6 (LdG biaxial) replace this with proper field-theoretic
          regularization. Don't trust the field's behavior at the
          single-voxel core for any physics test before M5.5.
        - The weighted superposition is NOT a topologically-clean seed —
          two opposite-sign hedgehogs blended this way may have winding
          number that integrates to a non-integer until gradient-descent
          relaxation (M5.1 tasks 4) settles them onto the topological
          minimum. The integer winding is RECOVERED post-relaxation.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    z_hat = ti.Vector([0.0, 0.0, 1.0])
    soft_core_voxels = ti.cast(0.2, ti.f32)  # match Exp 2's 0.2 · DX softening
    proximity_floor = ti.cast(0.5, ti.f32)  # match Exp 2's "r + 0.5"

    for i, j, k in ti.ndrange(nx, ny, nz):
        x = ti.cast(i, ti.f32)
        y = ti.cast(j, ti.f32)
        z = ti.cast(k, ti.f32)

        n_combined = ti.Vector([0.0, 0.0, 0.0])
        r_nearest = ti.cast(1e10, ti.f32)

        for d in range(n_defects):
            cx = ti.cast(centers_voxel[d, 0], ti.f32)
            cy = ti.cast(centers_voxel[d, 1], ti.f32)
            cz = ti.cast(centers_voxel[d, 2], ti.f32)
            sgn = ti.cast(signs[d], ti.f32)

            ddx = x - cx
            ddy = y - cy
            ddz = z - cz
            r = ti.sqrt(ddx * ddx + ddy * ddy + ddz * ddz)
            r_safe = ti.max(r, soft_core_voxels)

            n_radial = ti.Vector(
                [
                    sgn * ddx / r_safe,
                    sgn * ddy / r_safe,
                    sgn * ddz / r_safe,
                ]
            )

            w_defect = 1.0 / (r + proximity_floor)
            n_combined += w_defect * n_radial

            r_nearest = ti.min(r_nearest, r)

        # Vacuum blend: ~1 within D/4, falls as r⁻⁴ far away
        ratio = r_nearest / domain_quarter_voxels
        w_vac = 1.0 / (1.0 + ratio * ratio * ratio * ratio)

        n_blended = w_vac * n_combined + (1.0 - w_vac) * z_hat

        # Renormalize to unit length (ε in denominator avoids /0 at the core)
        norm = n_blended.norm() + 1e-12
        n_unit = n_blended / norm

        wave_field.psi_am[i, j, k] = n_unit
        wave_field.psi_prev_am[i, j, k] = n_unit
        # BC consistency: write psi_new_am too. See evolve_psi docstring.
        wave_field.psi_new_am[i, j, k] = n_unit


# ================================================================
# DIFFERENTIAL OPERATORS
# ================================================================
# Vector calculus on the Vector(3) field ψ stored in psi_am. All operators
# are 2nd-order central-difference @ti.func kernels — caller is responsible
# for skipping boundary voxels (each operator's halo requirement is noted in
# its docstring).
#
#   compute_laplacian   — ∇²ψ          → Vector(3); 1-cell halo
#   compute_divergence  — ∇·ψ          → scalar;    1-cell halo
#   compute_curl        — ∇×ψ          → Vector(3); 1-cell halo
#   compute_curl_curl   — ∇×(∇×ψ)      → Vector(3); 2-cell halo (uses identity)
#
# The curl-curl is implemented via the vector identity
#       ∇×(∇×ψ) = ∇(∇·ψ) − ∇²ψ
# rather than applying compute_curl twice. Two reasons:
#   1. Reuses validated Laplacian — only the gradient-of-divergence part is
#      new code, easier to verify
#   2. 2-cell halo (vs. 4-cell for nested curl) keeps the loop range tight
#      and matches Exp 7 v2's implementation (the reference)
#
# Why @ti.func not @ti.kernel: these are inner-loop primitives meant to be
# called from kernels that iterate over voxels (evolve_psi, future force
# kernels, divergence-cleaning projections, etc.), not standalone passes.


@ti.func
def compute_laplacian(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """
    Compute the vector Laplacian ∇²ψ at voxel [i, j, k] using a 6-point stencil.

    Reads from `wave_field.psi_am` (the current-time buffer, ψ at t).
    Returns a Vector(3) — the discrete Laplacian acts component-wise:
        ∇²ψ = (∇²ψ_x, ∇²ψ_y, ∇²ψ_z)

    Standard 2nd-order finite-difference stencil:
        ∇²ψ ≈ (face_sum − 6·center) / dx²
    where face_sum = ψ[i±1] + ψ[j±1] + ψ[k±1] (six face neighbors).

    Note: the L/T (longitudinal vs transverse) decomposition of the field is a
    property of the field configuration (radial vs tangential displacement),
    NOT of the operator. M5's vector field handles both together; the L/T view
    is computed post hoc by projecting ψ onto radial/tangential directions.

    Args:
        wave_field: WaveField instance (reads psi_am)
        i, j, k: Voxel indices (must be interior: 0 < i < nx-1, etc. — caller
            handles boundary skip)

    Returns:
        ∇²ψ as Vector(3) in units [am / am²] = [1/am]
    """
    face_sum = (
        wave_field.psi_am[i + 1, j, k]
        + wave_field.psi_am[i - 1, j, k]
        + wave_field.psi_am[i, j + 1, k]
        + wave_field.psi_am[i, j - 1, k]
        + wave_field.psi_am[i, j, k + 1]
        + wave_field.psi_am[i, j, k - 1]
    )
    center = wave_field.psi_am[i, j, k]
    laplacian_psi_am = (face_sum - 6.0 * center) / (wave_field.dx_am**2)
    return laplacian_psi_am


@ti.func
def compute_divergence(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """
    Compute the scalar divergence ∇·ψ at voxel [i, j, k].

    Standard 2nd-order central-difference stencil:
        ∇·ψ ≈ (∂_x ψ_x + ∂_y ψ_y + ∂_z ψ_z)
            = [ ψ_x[i+1] − ψ_x[i−1]
              + ψ_y[j+1] − ψ_y[j−1]
              + ψ_z[k+1] − ψ_z[k−1] ] / (2·dx)

    Used by:
      - M5.1 winding-number tracker (∇·n on the director field is an
        observable near hedgehog defects)
      - M5.2 Close Eq. 23 — the constraint ∇·s = 0 must be enforced;
        divergence-cleaning projection reads this each step
      - M5.0e curl_curl identity ∇×(∇×ψ) = ∇(∇·ψ) − ∇²ψ

    Args:
        wave_field: WaveField instance (reads psi_am)
        i, j, k: Voxel indices (must be interior: 0 < i < nx-1, etc. —
            caller handles boundary skip; 1-cell halo)

    Returns:
        ∇·ψ as ti.f32 scalar in units [am / am] = [dimensionless]
    """
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)
    dpsi_x_dx = (wave_field.psi_am[i + 1, j, k][0] - wave_field.psi_am[i - 1, j, k][0]) * inv_2dx
    dpsi_y_dy = (wave_field.psi_am[i, j + 1, k][1] - wave_field.psi_am[i, j - 1, k][1]) * inv_2dx
    dpsi_z_dz = (wave_field.psi_am[i, j, k + 1][2] - wave_field.psi_am[i, j, k - 1][2]) * inv_2dx
    return dpsi_x_dx + dpsi_y_dy + dpsi_z_dz


@ti.func
def compute_curl(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """
    Compute the vector curl ∇×ψ at voxel [i, j, k].

    Standard 2nd-order central-difference stencil:
        (∇×ψ)_x = ∂_y ψ_z − ∂_z ψ_y
        (∇×ψ)_y = ∂_z ψ_x − ∂_x ψ_z
        (∇×ψ)_z = ∂_x ψ_y − ∂_y ψ_x

    Used by:
      - M5.2 Close Eq. 19 linear limit  ∂²_t Q = −c² · ∇×(∇×Q)
      - Spin-density observables (s = curl of velocity field analog)
      - Magnetic-channel diagnostics in M7+ (transverse component extraction)

    Args:
        wave_field: WaveField instance (reads psi_am)
        i, j, k: Voxel indices (must be interior: 0 < i < nx-1, etc. —
            caller handles boundary skip; 1-cell halo)

    Returns:
        ∇×ψ as Vector(3) in units [am / am] = [dimensionless]
    """
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)
    # ∂_y ψ_z, ∂_z ψ_y
    dpsi_z_dy = (wave_field.psi_am[i, j + 1, k][2] - wave_field.psi_am[i, j - 1, k][2]) * inv_2dx
    dpsi_y_dz = (wave_field.psi_am[i, j, k + 1][1] - wave_field.psi_am[i, j, k - 1][1]) * inv_2dx
    # ∂_z ψ_x, ∂_x ψ_z
    dpsi_x_dz = (wave_field.psi_am[i, j, k + 1][0] - wave_field.psi_am[i, j, k - 1][0]) * inv_2dx
    dpsi_z_dx = (wave_field.psi_am[i + 1, j, k][2] - wave_field.psi_am[i - 1, j, k][2]) * inv_2dx
    # ∂_x ψ_y, ∂_y ψ_x
    dpsi_y_dx = (wave_field.psi_am[i + 1, j, k][1] - wave_field.psi_am[i - 1, j, k][1]) * inv_2dx
    dpsi_x_dy = (wave_field.psi_am[i, j + 1, k][0] - wave_field.psi_am[i, j - 1, k][0]) * inv_2dx
    return ti.Vector(
        [
            dpsi_z_dy - dpsi_y_dz,
            dpsi_x_dz - dpsi_z_dx,
            dpsi_y_dx - dpsi_x_dy,
        ]
    )


@ti.func
def compute_curl_curl(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """
    Compute ∇×(∇×ψ) at voxel [i, j, k] via the vector identity
        ∇×(∇×ψ) = ∇(∇·ψ) − ∇²ψ

    The gradient of divergence is computed by central-differencing
    compute_divergence at the six face neighbors — that's what makes
    this a 2-cell halo operator (each face neighbor itself needs a 1-cell
    halo for its own divergence stencil, so [i±2, j, k] etc. are read).

    The Laplacian piece reuses the validated compute_laplacian.

    Used by M5.2 Close Eq. 19 linear limit  ∂²_t Q = −c² · ∇×(∇×Q).

    Args:
        wave_field: WaveField instance (reads psi_am)
        i, j, k: Voxel indices (must be interior with 2-cell halo:
            1 < i < nx-2, etc. — caller handles boundary skip)

    Returns:
        ∇×(∇×ψ) as Vector(3) in units [am / am²] = [1/am]
    """
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)

    # ∇(∇·ψ): central-difference the divergence at the six face neighbors
    div_xp = compute_divergence(wave_field, i + 1, j, k)
    div_xm = compute_divergence(wave_field, i - 1, j, k)
    div_yp = compute_divergence(wave_field, i, j + 1, k)
    div_ym = compute_divergence(wave_field, i, j - 1, k)
    div_zp = compute_divergence(wave_field, i, j, k + 1)
    div_zm = compute_divergence(wave_field, i, j, k - 1)

    grad_div = ti.Vector(
        [
            (div_xp - div_xm) * inv_2dx,
            (div_yp - div_ym) * inv_2dx,
            (div_zp - div_zm) * inv_2dx,
        ]
    )

    laplacian = compute_laplacian(wave_field, i, j, k)

    return grad_div - laplacian


# ================================================================
# ψ EVOLVING ENGINE — LAGRANGIAN FIELD EVOLUTION
# ================================================================
# evolve_psi evolves the field ψ via leapfrog/Verlet integration of the
# wave equation:
#     ∂²ψ/∂t² = c²·∇²ψ − ∂V/∂ψ
#
# In M5.0d this lands with V(ψ) = 0 (free wave). The nonlinear potential term
# −∂V/∂ψ is added in M5.2 (Close's Eq. 23 + Klein-Gordon mass term + LdG).
#
# Why "evolve_psi" not "propagate_wave": the operation is field evolution,
# not wave-specific. Leapfrog also evolves topology (defect drift), gradient
# descent for relaxation reuses the same buffers, etc. Naming follows what the
# function operates on (the field ψ), not the method (leapfrog) or one of the
# behaviors (wave propagation).
#
# Caller MUST call wave_field.swap_buffers() after this kernel returns, to
# cycle the triple buffer: psi_prev ← psi, psi ← psi_new.


@ti.kernel
def evolve_psi(
    wave_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    m_freq_rs: ti.f32,  # type: ignore
):
    """
    Evolve ψ one timestep via leapfrog/Verlet integration.

    Klein-Gordon Equation (V = ½·m²·|ψ|²):
        ∂²ψ/∂t² = c²·∇²ψ − m²·ψ

    Discrete (leapfrog/Verlet):
        ψ_new = 2·ψ − ψ_prev + (c·dt)²·∇²ψ − (m·dt)²·ψ

    Setting m_freq_rs = 0 recovers the free wave equation (V = 0).

    Reads:  wave_field.psi_am, wave_field.psi_prev_am
    Writes: wave_field.psi_new_am

    Boundary: fixed-value Dirichlet — boundary voxels are never updated;
    they stay at whatever value the active seeder wrote. The interpretation
    of that value depends on the seeder:
      - M5.0 wave seeders (seed_gaussian, seed_dispersion_modes) leave
        boundaries at ψ ≈ 0 (zero-displacement BC; quiet wave boundary)
      - M5.1+ topology seeders (seed_vacuum, seed_hedgehog) leave boundaries
        at ψ = ẑ (vacuum-director BC; far-field defect reference)
    Both are valid Dirichlet BCs at different fixed values. The propagator
    is agnostic — it just iterates (1, n−1) on each axis (the 6-point
    Laplacian's 1-cell halo requirement) and never touches the boundary.

    BC consistency requirement: ALL THREE BUFFERS (psi_prev_am, psi_am,
    psi_new_am) must hold the same boundary value before the first call.
    psi_new_am's boundary is never written by the propagator, so swap_buffers'
    `psi_am.copy_from(psi_new_am)` would otherwise clobber the seeded
    boundary value with psi_new_am's default zero. Seeders fix this by
    writing all three buffers; verified bug 2026-05-09 with seed_hedgehog.

    CFL stability: dt ≤ dx / (c·√3) for 3D wave equation. Caller is
    responsible for sizing dt below this bound (see _launcher CFL eval).

    Args:
        wave_field: WaveField instance (reads psi_am, psi_prev_am; writes psi_new_am)
        c_amrs: wave speed in scaled units (am/rs); already includes any
            slow-motion factor for visualization
        dt_rs: timestep in rontoseconds, sized below the CFL bound
        m_freq_rs: Klein-Gordon mass-frequency m·c²/ℏ in rad/rs storage units
            (electron: ~7.76e-7 rad/rs at SIM_SPEED=1; 0 disables the mass term)
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    c_dt_squared = (c_amrs * dt_rs) ** 2
    m_dt_squared = (m_freq_rs * dt_rs) ** 2

    # Interior voxels only (Dirichlet BC; 6-point Laplacian needs 1-cell halo)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        laplacian_psi_am = compute_laplacian(wave_field, i, j, k)

        # Leapfrog/Verlet update: ψ_new = 2·ψ − ψ_prev + (c·dt)²·∇²ψ − (m·dt)²·ψ
        wave_field.psi_new_am[i, j, k] = (
            2.0 * wave_field.psi_am[i, j, k]
            - wave_field.psi_prev_am[i, j, k]
            + c_dt_squared * laplacian_psi_am
            - m_dt_squared * wave_field.psi_am[i, j, k]
        )


# ================================================================
# TRACKERS — PER-VOXEL AMPLITUDE & FREQUENCY
# ================================================================
# Trackers run after evolve_psi each step. They derive observable scalars
# from ψ (the field) without touching it — read-only on psi_am, write to the
# tracker fields in-place. EMA / zero-crossing logic preserved from M4 but
# adapted to vector ψ:
#   - amp_local_emarms_am  : RMS of |ψ| over recent history (EMA on |ψ|²)
#   - freq_local_cross_rHz : zero-crossing rate of ψ_z (chosen polarization axis)
#
# PERFORMANCE NOTE (deferred to M5.0i): this is a SECOND full-grid pass over
# ψ — evolve_psi already touched every voxel; we re-stream ψ and ψ_prev
# from memory just to compute trackers. At 100M voxels this costs ~5 GB of
# extra memory bandwidth per step (~25–30% of total per-step traffic).
#
# M2 avoided this by computing trackers INSIDE the same loop as the wave
# update, but M2's wave kernel was analytical (no neighbor dependency), so
# trackers could read the just-written ψ trivially. M5's leapfrog can do the
# same: in a merged kernel, write ψ_new then immediately use it (with
# ψ_prev already in registers) for tracker EMA / zero-crossing — gives true
# central-difference ψ̇ for free, more accurate than the forward diff used
# below (psi_new gets consumed by swap_buffers before update_trackers
# can see it, so we fall back to (ψ − ψ_prev) / dt here).
#
# Splitting was deliberate for M5.0d: cleaner kernel boundaries, AMR-ready
# (octree retrofit in M5.6/M5.8 prefers single-purpose kernels), and trackers
# are skippable in tests (e.g., M5.0h dispersion regression). Profile-first
# in M5.0i — alongside swap_buffers' rotating-pointer optimization (see
# medium.py:swap_buffers docstring), this is the second-largest known
# bandwidth target on the per-step path.


@ti.kernel
def update_trackers(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
):
    """
    Update per-voxel amplitude and frequency trackers from ψ.

    Reads: wave_field.psi_am, wave_field.psi_prev_am
    Writes: trackers.amp_local_emarms_am,
            trackers.freq_local_cross_rHz,
            trackers.last_crossing

    Amplitude:
        |ψ|² is the squared vector magnitude (scalar). EMA:
            rms² ← α·|ψ|² + (1−α)·rms²_old
        with a small unconditional decay so stale regions fade.

    Frequency:
        Zero crossing tracked on ψ_z (polarization axis chosen by convention).
        Detect sign change between psi_prev and psi: positive-going crossing
        marks one period boundary. Period_rs = elapsed − last_crossing,
        smoothed via EMA. If no crossing: untouched, decays naturally.

    Note: vector polarization is unknown a priori — using ψ_z as the canonical
    crossing channel. For radial wave-center fields (annihilation1) this is
    the longitudinal component near the equator; for wave seeds along x
    with ŷ polarization, ψ_z stays at 0, so freq tracker stays at base value.
    M5.0e will replace this with magnitude-zero-crossing or per-component tracking.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    decay = ti.cast(0.999, ti.f32)  # ~1000 frames to ~37%
    alpha_rms = ti.cast(0.005, ti.f32)
    alpha_freq = ti.cast(0.05, ti.f32)
    min_period_rs = ti.cast(2.0, ti.f32) * dt_rs

    for i, j, k in ti.ndrange(nx, ny, nz):
        # RMS amplitude via EMA on |ψ|²
        psi_curr = wave_field.psi_am[i, j, k]
        psi_prev = wave_field.psi_prev_am[i, j, k]
        psi2 = psi_curr.norm_sqr()
        rms2_old = trackers.amp_local_emarms_am[i, j, k] ** 2
        rms2_new = alpha_rms * psi2 + (1.0 - alpha_rms) * rms2_old
        trackers.amp_local_emarms_am[i, j, k] = ti.sqrt(rms2_new) * decay

        # Zero-crossing frequency on ψ_z (canonical channel)
        prev_z = psi_prev[2]
        curr_z = psi_curr[2]
        if prev_z < 0.0 and curr_z >= 0.0:
            period_rs = elapsed_t_rs - trackers.last_crossing[i, j, k]
            if period_rs > min_period_rs:
                measured_freq = 1.0 / period_rs
                old_freq = trackers.freq_local_cross_rHz[i, j, k]
                trackers.freq_local_cross_rHz[i, j, k] = (
                    alpha_freq * measured_freq + (1.0 - alpha_freq) * old_freq
                )
            trackers.last_crossing[i, j, k] = elapsed_t_rs

        # Unconditional frequency decay (counteracted by zero-crossing updates)
        trackers.freq_local_cross_rHz[i, j, k] = trackers.freq_local_cross_rHz[i, j, k] * decay


# ================================================================
# POTENTIAL ENERGY HOOK — V(ψ)
# ================================================================
# Default V(ψ) = 0 (free wave) in M5.0g. M5.2 swaps this implementation in to
# add the nonlinear physics:
#   - Klein-Gordon mass:    V = ½ m² |ψ|²
#   - Close Eq. 23 nonlin:  V = … (from −u·∇s + w×s couplings)
#   - LdG biaxial:          V = polynomial in (Q, axes δ/1/g)
#
# --- NATURAL-UNIT BASIS — applied here in V(ψ) -----------------------------
# Basis convention: ℏ = c = λ̄_C(electron) = 1, where λ̄_C = ℏ/(m_e·c) is the
# reduced Compton wavelength of the electron (constants.COMPTON_WAVELENGTH_
# REDUCED_ELECTRON_AM). This makes m_e = 1 in natural units and aligns
# dimensional couplings with textbook field-theory form.
#
# Single bridge value: the Klein-Gordon mass-frequency m·c²/ℏ in rad/rs.
# Computed once in launcher as
#     m_freq_kg_rs = c_amrs / COMPTON_WAVELENGTH_REDUCED_ELECTRON_AM
# (electron rest-frame angular frequency, ~7.76e-7 rad/rs at SIM_SPEED=1).
# Threaded into V_psi via compute_energyH_density and (Step 2+) evolve_psi.
# For ψ in am, V_storage = ½·(m_freq_rs)²·|ψ|² produces an energy density in
# (am/rs)² — same units as the kinetic and gradient terms.
#
# M5.2 Step 2 (this version): plain Klein-Gordon mass term active —
# V_psi returns ½·m²·|ψ|² and evolve_psi adds −(m·dt)²·ψ to the leapfrog.
# Setting m_freq_rs = 0 from a caller recovers free-wave (V = 0) behavior;
# the M5.0h dispersion tests do this. Step 4 (Close Eq. 23) layers nonlinear
# terms on top so the potential minimum sits on the unit sphere (|ψ| = 1)
# instead of at ψ = 0 — required for hedgehog stability.
#
# Linear kernels (Laplacian, divergence, curl, curl-curl) stay in storage
# units (am, rs, rHz) throughout — dimensionally self-balancing, don't
# benefit from natural units. See M5.0f decision-record in 3d_path_to_m5.md.


@ti.func
def V_psi(
    psi: ti.template(),  # type: ignore
    m_freq_rs: ti.f32,  # type: ignore
):
    """
    Scalar potential V(ψ) at one voxel — Klein-Gordon mass term (M5.2 Step 2).

    V(ψ) = ½·m²·|ψ|²
    EL contribution: ∂V/∂ψ = m²·ψ → evolve_psi adds −(m·dt)²·ψ to the leapfrog.
    Setting m_freq_rs = 0 reduces this to the free wave (V = 0).

    NOTE: a plain KG term has its minimum at ψ = 0, which is the WRONG shape
    for stabilizing a |ψ|=1 hedgehog. M5.2 Step 4 (Close Eq. 23) and Step 6
    (LdG biaxial) layer Mexican-hat-style nonlinear terms on top so the
    minimum sits on the unit sphere. Plain KG by itself is mainly useful for
    verifying the V(ψ) hook + Klein-Gordon dispersion ω² = c²k² + m².

    Args:
        psi: Vector(3) field value at the voxel
        m_freq_rs: Klein-Gordon mass-frequency m·c²/ℏ in rad/rs storage units
            (electron: ~7.76e-7 rad/rs at SIM_SPEED=1; 0 disables the term)

    Returns:
        scalar V(ψ) in the same units as kinetic and gradient terms (am²/rs²)
    """
    return 0.5 * m_freq_rs * m_freq_rs * psi.norm_sqr()


# ================================================================
# ENERGY DENSITY (HAMILTONIAN) — PER-VOXEL ENERGY FIELD
# ================================================================
# H(x) = ½|ψ̇|² + ½c²|∇ψ|² + V(ψ)
#
# Naming convention: the quantity stored is ENERGY density (aJ per voxel).
# The "_H" suffix tags which formula was used to compute it (Hamiltonian).
# Future formulas would parallel: `_L` for Lagrangian density, `_K` for
# kinetic-only, etc.
#
# Populated by compute_energyH_density into observables.energyH_density_aJ
# each step; consumed by:
#   - xforce_motion.compute_force_vector  → F = −∇E per wave-center
#   - sample_avg_trackers (below)         → 3-plane mean → dashboard total
#   - launcher WAVE_MENU=4 flux mesh      → visualization
#
# This is a third full-grid pass per step (alongside evolve_psi and
# update_trackers). Performance optimization (merge into evolve_psi)
# is M5.0i territory; we land it as a separate kernel for clarity in M5.0g.


@ti.kernel
def compute_energyH_density(
    wave_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    m_freq_rs: ti.f32,  # type: ignore
):
    """
    Compute per-voxel energy density (Hamiltonian formula) into
    observables.energyH_density_aJ.
    H = ½|ψ̇|² + ½c²|∇ψ|² + V(ψ)

    Reads:  wave_field.psi_am, wave_field.psi_prev_am
    Writes: observables.energyH_density_aJ

    Boundary handling: 1-cell halo (gradient stencil); boundary voxels left
    at zero (Dirichlet BC in ψ ⇒ ψ ≈ 0 at edges ⇒ H ≈ 0 trivially).

    Kinetic ψ̇ uses forward difference (ψ − ψ_prev)/dt because psi_new is
    consumed by swap_buffers before this kernel runs — same compromise as
    update_trackers. Central-difference ψ̇ becomes available once those
    two kernels are merged into evolve_psi (M5.0i optimization).

    PHYSICAL-SCALING TODO — M5.2: the field is currently named `_aJ` but
    actually stores `(am/rs)²` per voxel — the kernel writes raw kinematic
    Hamiltonian density without the physical-energy conversion factor
    (ρ_medium × voxel_volume_am³ × INTERNAL_ENERGY_TO_AJ ≈ matches M4's
    `E = ρ·V·(f·A)²` formula). Sufficient for M5.0g–M5.0i because:
      (a) F = −∇E only depends on the gradient (ratios survive scaling)
      (b) Tests against Exp 4 KG dispersion check ω(k), not absolute E
    But once M5.2 introduces V(ψ) terms with explicit dimensional
    couplings (Klein-Gordon mass, Close Eq. 23, LdG potential), we MUST
    apply the physical-scaling factor here so the V term and the
    kinetic+gradient terms add in the same units. See dashboard "(rel.)"
    labels in _launcher.py for the corresponding display caveat.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    inv_dt = 1.0 / dt_rs
    half_c2 = 0.5 * c_amrs**2
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)

    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        psi = wave_field.psi_am[i, j, k]
        psi_prev = wave_field.psi_prev_am[i, j, k]

        # Kinetic: ½ |ψ̇|²
        psi_dot = (psi - psi_prev) * inv_dt
        kinetic = 0.5 * psi_dot.norm_sqr()

        # Gradient: ½ c² |∇ψ|²  (central-difference on each axis, component-wise)
        d_dx = (wave_field.psi_am[i + 1, j, k] - wave_field.psi_am[i - 1, j, k]) * inv_2dx
        d_dy = (wave_field.psi_am[i, j + 1, k] - wave_field.psi_am[i, j - 1, k]) * inv_2dx
        d_dz = (wave_field.psi_am[i, j, k + 1] - wave_field.psi_am[i, j, k - 1]) * inv_2dx
        gradient = half_c2 * (d_dx.norm_sqr() + d_dy.norm_sqr() + d_dz.norm_sqr())

        # Potential: V(ψ) — 0 in M5.0g + M5.2 Step 1; KG mass term lands in Step 2
        potential = V_psi(psi, m_freq_rs)

        observables.energyH_density_aJ[i, j, k] = kinetic + gradient + potential


# Note: the global energy aggregate (mean per voxel + grid total) is computed
# by the 3-plane sample_avg_trackers below alongside amp / freq — single pass
# over the slice planes, single GPU↔CPU sync per dashboard cadence. The launcher
# multiplies the returned mean by voxel_count to get the grid-total scalar.


# ================================================================
# FRANK ELASTIC ENERGY DENSITY (M5.1)
# ================================================================
# H_F = (K/2) · |∇n̂|²   per voxel scalar
# |∇n̂|² = Σᵢ Σ_α (∂ᵢ n̂_α)²   (9 central-difference terms: 3 axes × 3 components)
#
# This is the same gradient term that's already implicit in the wave equation's
# c²·∇²ψ Laplacian. The kernel doesn't add new physics — it provides the SCALAR
# E we need for the M5.1 downstream tests:
#   - task 6 (gradient descent):     monitor E falling monotonically
#   - task 7 (Coulomb 1/d fit):       E(d) data points across separations
#   - WAVE_MENU=5 visualization:      color the flux mesh by elastic stress
#
# Naming convention: parallels energyH_density_aJ — quantity is energy density,
# _F suffix tags the Frank formula (vs _H Hamiltonian, future _L Lagrangian).
# K_frank coupling currently hard-coded to 1.0 (Exp 2 baseline); M5.6 plumbs in
# the physical elastic constants alongside the LdG axis hierarchy.


K_FRANK = 1.0  # Frank elastic coupling — Exp 2 baseline; physical scaling in M5.6


@ti.kernel
def compute_energyF_density(
    wave_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    K_frank: ti.f32,  # type: ignore
):
    """
    Compute per-voxel Frank elastic energy density H_F = (K/2)·|∇n̂|² into
    observables.energyF_density_aJ.

    Reads:  wave_field.psi_am  (director field n̂; |n̂|=1 enforced by seeders
            and by gradient-descent relaxation step in M5.1 task 6)
    Writes: observables.energyF_density_aJ

    Same central-difference gradient stencil as compute_energyH_density's
    gradient term, sans the c² factor and sans the kinetic/potential terms.

    Boundary handling: 1-cell halo (gradient stencil); boundary voxels left
    at zero. This is consistent with energyH_density_aJ — both kernels skip
    the outermost voxels.

    PHYSICAL-SCALING TODO — M5.6: the field is named `_aJ` aspirationally but
    actually stores `(K_frank · (am/am)²)` per voxel (dimensionless until the
    LdG elastic constants land in physical units). Sufficient for M5.1 because:
      (a) gradient descent on F depends only on |∇n̂|² gradient (relative)
      (b) Coulomb 1/d fit cares about the FUNCTIONAL form, not absolute scale
    Once M5.6 introduces the biaxial LdG potential with physical elastic
    constants (k1, k2, k3 — splay/twist/bend), apply the appropriate factor
    here so F adds to the H total in the same physical units.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    half_K = 0.5 * K_frank
    inv_2dx = 1.0 / (2.0 * wave_field.dx_am)

    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        # Gradient: 9 central-difference terms — ∂_x/y/z applied to each of 3 components
        d_dx = (wave_field.psi_am[i + 1, j, k] - wave_field.psi_am[i - 1, j, k]) * inv_2dx
        d_dy = (wave_field.psi_am[i, j + 1, k] - wave_field.psi_am[i, j - 1, k]) * inv_2dx
        d_dz = (wave_field.psi_am[i, j, k + 1] - wave_field.psi_am[i, j, k - 1]) * inv_2dx
        grad_n_sqr = d_dx.norm_sqr() + d_dy.norm_sqr() + d_dz.norm_sqr()
        observables.energyF_density_aJ[i, j, k] = half_K * grad_n_sqr


# ================================================================
# GRADIENT-DESCENT RELAXATION (M5.1 task 6)
# ================================================================
# Direct port of Exp 2's relax() inner loop to Taichi:
#     dn = ∇²n − (n·∇²n) n         (tangent projection onto unit sphere)
#     n_new = n + τ · dn           (gradient-descent step)
#     n_new ← n_new / |n_new|       (unit-length constraint)
#     n_new[defect_core] = ±ẑ      (soft Dirichlet pin)
#
# Stability: τ < dx²/(2·dim·K) by the standard heat-equation CFL condition.
# For dim=3 in 3D the bound is τ < dx²/(6·K). Use ~50% of the bound for
# headroom; the launcher computes this from wave_field.dx_am at runtime.
#
# Output buffer convention: writes to psi_new_am (same buffer leapfrog uses).
# The caller copies psi_new_am to BOTH psi_am AND psi_prev_am so that
# subsequent leapfrog steps see ψ̇ = 0 (no spurious time-derivative from the
# relaxation update — relaxation is a STATIC operation, not a dynamics step).
# This is different from swap_buffers, which rotates buffers for time evolution.
#
# Boundary handling: psi_new_am boundary voxels are copied from psi_am (no
# update) — preserves the fixed-value Dirichlet BC (e.g. ẑ at the universe
# edge from seed_vacuum).
#
# Core pin: soft Dirichlet at the SINGLE closest voxel to each defect center.
# Without this, gradient descent can numerically dissolve the defect on a
# discrete grid (topology is not strictly preserved by discretization).
#
# M7 FORWARD-LINK: this kernel is mathematically the γ → ∞ limit of a damped
# wave equation `∂²_t ψ = c²∇²ψ − γ·∂_t ψ`. Although the relaxation use is
# purely numerical (no physics interpretation in M5.1), the same primitives
# — Laplacian stencil, tangent projection, soft-core pin — will be reused
# in M7 thermal-modulation kernels where γ becomes a PHYSICAL damping
# coefficient (radiation loss, phonon coupling, EM-load impedance). See
# `research_hub/3d_path_to_m5.md § Beyond M6 — thermal mechanics pathway`
# for the infrastructure-foundation discussion.


@ti.kernel
def relax_director_step(
    wave_field: ti.template(),  # type: ignore
    tau: ti.f32,  # type: ignore
    pin_centers: ti.template(),  # type: ignore
    pin_signs: ti.template(),  # type: ignore
    n_defects: ti.i32,  # type: ignore
):
    """
    One gradient-descent step on the Frank elastic energy `H_F = (K/2)·|∇n̂|²`.

    Reads:  wave_field.psi_am
    Writes: wave_field.psi_new_am  (caller copies → psi_am AND psi_prev_am)

    Args:
        wave_field: WaveField with the director field in psi_am
        tau: step size (must satisfy τ < dx²/(6·K) for stability)
        pin_centers: ti.field shape (n_defects, 3) i32 — defect centers in voxel coords
        pin_signs: ti.field shape (n_defects,) i32 — ±1 per defect (pin direction)
        n_defects: number of active defects in pin_centers/signs

    Note: K_frank is implicit in `tau`'s CFL derivation — the kernel itself
    doesn't multiply by K because gradient descent on `(K/2)·|∇n̂|²` gives
    `∂n/∂τ = K·∇²n` and the K can be absorbed into the step size (τ_eff = K·τ).
    For K=1 (M5.1 default) this is τ_eff = τ.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    # ---- Interior: tangent-projected gradient step + renormalization ----
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        n = wave_field.psi_am[i, j, k]
        lap = compute_laplacian(wave_field, i, j, k)
        n_dot_lap = n.dot(lap)
        dn = lap - n_dot_lap * n  # tangent projection (remove component along n)
        n_new = n + tau * dn
        norm = n_new.norm() + 1e-12
        wave_field.psi_new_am[i, j, k] = n_new / norm

    # ---- Boundary: preserve Dirichlet BC by copying from current ψ ----
    for j, k in ti.ndrange(ny, nz):
        wave_field.psi_new_am[0, j, k] = wave_field.psi_am[0, j, k]
        wave_field.psi_new_am[nx - 1, j, k] = wave_field.psi_am[nx - 1, j, k]
    for i, k in ti.ndrange(nx, nz):
        wave_field.psi_new_am[i, 0, k] = wave_field.psi_am[i, 0, k]
        wave_field.psi_new_am[i, ny - 1, k] = wave_field.psi_am[i, ny - 1, k]
    for i, j in ti.ndrange(nx, ny):
        wave_field.psi_new_am[i, j, 0] = wave_field.psi_am[i, j, 0]
        wave_field.psi_new_am[i, j, nz - 1] = wave_field.psi_am[i, j, nz - 1]

    # ---- Pin defect cores (soft Dirichlet at single closest voxel) ----
    for d in range(n_defects):
        ci = pin_centers[d, 0]
        cj = pin_centers[d, 1]
        ck = pin_centers[d, 2]
        sgn = ti.cast(pin_signs[d], ti.f32)
        wave_field.psi_new_am[ci, cj, ck] = ti.Vector([0.0, 0.0, sgn])


# ================================================================
# WINDING-NUMBER DIAGNOSTIC (M5.1 task 8)
# ================================================================
# Compute the topological charge Q of a director field on a sphere around
# a given center. For a +1 hedgehog Q ≈ +1, anti-hedgehog Q ≈ −1, vacuum
# Q ≈ 0. Pure CPU numpy; called occasionally (once per dashboard frame) as
# a diagnostic, not per-step. Port of Exp 3's `winding_number()`.
#
# Math:
#     Q = (1/4π) ∮_S n̂ · (∂_θ n̂ × ∂_φ n̂) dθ dφ
# Where S is a sphere of given radius around the defect center.
# Method: trilinear-sample ψ on a regular (θ, φ) grid; finite-difference
# the angular derivatives; integrate.


def compute_winding_number(
    psi_np,
    center_vox,
    radius_vox,
    n_theta=32,
    n_phi=64,
):
    """Compute topological winding number Q on a sphere of radius `radius_vox`
    around `center_vox` in voxel coordinates.

    Args:
        psi_np: (nx, ny, nz, 3) numpy array of director-field values
                (typically `wave_field.psi_am.to_numpy()`).
        center_vox: (cx, cy, cz) tuple of voxel-coord defect center.
        radius_vox: sphere radius in voxel units. Should be small enough that
                    the sphere fits inside the grid but large enough to be
                    away from the defect core (e.g., 4-10 voxels).
        n_theta: number of polar samples (default 32).
        n_phi: number of azimuthal samples (default 64).

    Returns:
        Q (float): the winding number. ≈ +1 for hedgehog, −1 for anti-hedgehog,
                   0 for vacuum-only enclosure. Discrete grids give Q within
                   ~5% of the integer for a typical relaxed defect.
    """
    nx, ny, nz, _ = psi_np.shape
    cx, cy, cz = center_vox

    # Spherical sample grid
    theta = np.linspace(1e-3, np.pi - 1e-3, n_theta)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    dth = float(theta[1] - theta[0])
    dph = float(phi[1] - phi[0])
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    sx = cx + radius_vox * np.sin(TH) * np.cos(PH)
    sy = cy + radius_vox * np.sin(TH) * np.sin(PH)
    sz = cz + radius_vox * np.cos(TH)

    # Trilinear sample of ψ at the sphere points
    i0 = np.clip(np.floor(sx).astype(int), 0, nx - 2)
    j0 = np.clip(np.floor(sy).astype(int), 0, ny - 2)
    k0 = np.clip(np.floor(sz).astype(int), 0, nz - 2)
    fx = np.clip(sx - i0, 0.0, 1.0)
    fy = np.clip(sy - j0, 0.0, 1.0)
    fz = np.clip(sz - k0, 0.0, 1.0)

    n_s = np.zeros(sx.shape + (3,), dtype=np.float64)
    for dxi in (0, 1):
        for dyi in (0, 1):
            for dzi in (0, 1):
                w = (fx if dxi else 1 - fx) * (fy if dyi else 1 - fy) * (fz if dzi else 1 - fz)
                n_s += w[..., None] * psi_np[i0 + dxi, j0 + dyi, k0 + dzi, :]

    # Re-normalize the sampled directors (interpolation slightly shrinks |n̂|)
    norm = np.linalg.norm(n_s, axis=-1, keepdims=True)
    n_s = n_s / np.maximum(norm, 1e-12)

    # Angular derivatives (central differences; φ wraps periodic)
    dn_dtheta = np.zeros_like(n_s)
    dn_dphi = np.zeros_like(n_s)
    dn_dtheta[1:-1] = (n_s[2:] - n_s[:-2]) / (2 * dth)
    dn_dtheta[0] = (n_s[1] - n_s[0]) / dth
    dn_dtheta[-1] = (n_s[-1] - n_s[-2]) / dth
    dn_dphi[:, :] = (np.roll(n_s, -1, axis=1) - np.roll(n_s, 1, axis=1)) / (2 * dph)

    # Integrand: n̂ · (∂_θ n̂ × ∂_φ n̂). The sinθ factor is implicit in the
    # cross product for a unit-vector field on a parameterized sphere.
    cross = np.cross(dn_dtheta, dn_dphi, axis=-1)
    integrand = (n_s * cross).sum(axis=-1)
    return float(integrand.sum() * dth * dph / (4.0 * np.pi))


# ================================================================
# POSITION RENDER
# ================================================================


@ti.kernel
def sample_position_to_render(
    wave_field: ti.template(),  # type: ignore
    amp_boost: ti.f32,  # type: ignore
    stride: ti.i32,  # type: ignore
    num_render: ti.i32,  # type: ignore
):
    """Sample granule positions from z flux_mesh plane with stride, writing to 1D position_render.

    Samples every `stride`-th voxel from the XY plane at the z flux mesh index,
    capping output at `num_render` particles for performance.
    """
    k = int(wave_field.flux_mesh_planes[2] * wave_field.grid_size[2])
    max_dim = ti.cast(wave_field.max_grid_size, ti.f32)
    sampled_ny = (wave_field.ny + stride - 1) // stride  # cols per sampled row

    for render_idx in range(num_render):
        si = render_idx // sampled_ny  # sampled row index
        sj = render_idx % sampled_ny  # sampled col index
        i = si * stride
        j = sj * stride
        displaced = amp_boost * wave_field.psi_am[i, j, k] / wave_field.dx_am + ti.Vector(
            [ti.cast(i, ti.f32), ti.cast(j, ti.f32), ti.cast(k, ti.f32)]
        )
        wave_field.position_render[render_idx] = displaced / max_dim


# ================================================================
# 3-PLANE SAMPLING FOR AVERAGE TRACKERS
# ================================================================
# PERFORMANCE NOTE: Full GPU reduction (atomic_add over all voxels) causes
# severe performance issues due to atomic contention with millions of voxels.
# The 3-plane sampling is a deliberate compromise:
# - Samples ~3N² voxels instead of N³ (e.g., 3% for 100³ grid)
# - Assumes isotropic field distribution (valid for most wave scenarios)
# - Acceptable accuracy vs massive performance gain
#
# AUTO-REDUCE CAVEAT (M5.0h, 2026-05-08): Taichi's "implicit sum reduction"
# pattern (`s = 0; for ... in ti.ndrange: s += ...`) lowers to atomic_add on
# the Metal backend and hits the same contention wall as explicit atomic_add.
# At 63³ × ~5 reductions/step the test pinned the GPU at 100% and stalled
# at ~25 steps/min. Workaround used in m5_0h_dispersion: replace the global
# inner product with sparse point sampling (one voxel read per mode of
# interest, FFT recovers ω from the mixed time series). Don't assume "auto-
# reduction" sidesteps this on Metal — it doesn't.
# ================================================================

# Cached slice buffers — separate caches per domain (Trackers vs FieldObservables)
# so each domain owns its own 3-plane sampling pass. Per `feedback_visual_rendering_priority`,
# each class is self-contained: Trackers samples amp/freq, FieldObservables samples H/F.
_slice_xy_amp = None
_slice_xy_freq = None
_slice_xz_amp = None
_slice_xz_freq = None
_slice_yz_amp = None
_slice_yz_freq = None
_slice_xy_energyH = None
_slice_xy_energyF = None
_slice_xz_energyH = None
_slice_xz_energyF = None
_slice_yz_energyH = None
_slice_yz_energyF = None


# --- TRACKERS slice-copy kernels (amp + freq) -----------------------------


@ti.kernel
def _copy_slice_xy_trackers(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """Copy amp / freq at the XY center plane (z=mid_z) to 2D buffers."""
    for i, j in slice_amp:
        slice_amp[i, j] = trackers.amp_local_emarms_am[i, j, mid_z]
        slice_freq[i, j] = trackers.freq_local_cross_rHz[i, j, mid_z]


@ti.kernel
def _copy_slice_xz_trackers(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """Copy amp / freq at the XZ center plane (y=mid_y) to 2D buffers."""
    for i, k in slice_amp:
        slice_amp[i, k] = trackers.amp_local_emarms_am[i, mid_y, k]
        slice_freq[i, k] = trackers.freq_local_cross_rHz[i, mid_y, k]


@ti.kernel
def _copy_slice_yz_trackers(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """Copy amp / freq at the YZ center plane (x=mid_x) to 2D buffers."""
    for j, k in slice_amp:
        slice_amp[j, k] = trackers.amp_local_emarms_am[mid_x, j, k]
        slice_freq[j, k] = trackers.freq_local_cross_rHz[mid_x, j, k]


# --- FIELD-OBSERVABLES slice-copy kernels (energyH + energyF) -------------


@ti.kernel
def _copy_slice_xy_observables(
    observables: ti.template(),  # type: ignore
    slice_energyH: ti.template(),  # type: ignore
    slice_energyF: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """Copy energyH / energyF at the XY center plane (z=mid_z) to 2D buffers."""
    for i, j in slice_energyH:
        slice_energyH[i, j] = observables.energyH_density_aJ[i, j, mid_z]
        slice_energyF[i, j] = observables.energyF_density_aJ[i, j, mid_z]


@ti.kernel
def _copy_slice_xz_observables(
    observables: ti.template(),  # type: ignore
    slice_energyH: ti.template(),  # type: ignore
    slice_energyF: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """Copy energyH / energyF at the XZ center plane (y=mid_y) to 2D buffers."""
    for i, k in slice_energyH:
        slice_energyH[i, k] = observables.energyH_density_aJ[i, mid_y, k]
        slice_energyF[i, k] = observables.energyF_density_aJ[i, mid_y, k]


@ti.kernel
def _copy_slice_yz_observables(
    observables: ti.template(),  # type: ignore
    slice_energyH: ti.template(),  # type: ignore
    slice_energyF: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """Copy energyH / energyF at the YZ center plane (x=mid_x) to 2D buffers."""
    for j, k in slice_energyH:
        slice_energyH[j, k] = observables.energyH_density_aJ[mid_x, j, k]
        slice_energyF[j, k] = observables.energyF_density_aJ[mid_x, j, k]


# --- Sampler functions ----------------------------------------------------


def sample_avg_trackers(
    wave_field,
    trackers,
):
    """
    Estimate global Trackers aggregates by sampling 3 orthogonal center planes.

    Computes:
        trackers.amp_global_emarms_am   ← √⟨amp²⟩  (RMS over the 3 planes)
        trackers.freq_global_avg_rHz    ← ⟨freq⟩   (mean over the 3 planes)

    Why 3-plane sampling: full 3D reductions with GPU atomic_add cause severe
    atomic contention at million-voxel scale and stalled the GUI for tens of
    seconds in early M5.0d.2 testing. 3-plane sampling reads ~3·N² voxels
    (3% of a 100³ grid), is exact for spatially-uniform fields, and is a sound
    estimator for smooth wave packets.

    Companion sampler: `sample_avg_observables` (handles energyH / energyF
    via the FieldObservables class). Each domain owns its own pass for clean
    separation of concerns.

    Args:
        wave_field: WaveField (grid dimensions)
        trackers: Trackers (reads amp/freq per-voxel fields; writes the two
            global aggregates)
    """
    global _slice_xy_amp, _slice_xy_freq
    global _slice_xz_amp, _slice_xz_freq
    global _slice_yz_amp, _slice_yz_freq

    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    # Initialize slice buffers once
    if _slice_xy_amp is None:
        _slice_xy_amp = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_freq = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xz_amp = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_freq = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_yz_amp = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_freq = ti.field(dtype=ti.f32, shape=(ny, nz))

    # Copy 3 center slices to 2D buffers (parallel kernels)
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    _copy_slice_xy_trackers(trackers, _slice_xy_amp, _slice_xy_freq, mid_z)
    _copy_slice_xz_trackers(trackers, _slice_xz_amp, _slice_xz_freq, mid_y)
    _copy_slice_yz_trackers(trackers, _slice_yz_amp, _slice_yz_freq, mid_x)

    # Transfer 2D slices to CPU for numpy operations (exclude boundary voxels)
    xy_amp = _slice_xy_amp.to_numpy()[1:-1, 1:-1]
    xy_freq = _slice_xy_freq.to_numpy()[1:-1, 1:-1]
    xz_amp = _slice_xz_amp.to_numpy()[1:-1, 1:-1]
    xz_freq = _slice_xz_freq.to_numpy()[1:-1, 1:-1]
    yz_amp = _slice_yz_amp.to_numpy()[1:-1, 1:-1]
    yz_freq = _slice_yz_freq.to_numpy()[1:-1, 1:-1]

    n_samples = xy_amp.size + xz_amp.size + yz_amp.size

    # RMS amplitude: √(⟨A²⟩) — energy-weighting convention
    total_amp_squared = (xy_amp**2).sum() + (xz_amp**2).sum() + (yz_amp**2).sum()
    trackers.amp_global_emarms_am[None] = float(np.sqrt(total_amp_squared / n_samples))

    # Plain mean for frequency
    total_freq = xy_freq.sum() + xz_freq.sum() + yz_freq.sum()
    trackers.freq_global_avg_rHz[None] = float(total_freq / n_samples)


def sample_avg_observables(
    wave_field,
    observables,
):
    """
    Estimate global FieldObservables aggregates by sampling 3 orthogonal planes.

    Computes:
        observables.energyH_global_avg_aJ ← ⟨E_H⟩  (mean Hamiltonian density)
        observables.energyF_global_avg_aJ ← ⟨E_F⟩  (mean Frank elastic density)

    Caller (launcher) derives the grid total trivially:
        energyH_total_aJ = energyH_global_avg_aJ × voxel_count

    Same 3-plane sampling rationale as `sample_avg_trackers`. Separate pass so
    FieldObservables is self-contained.

    Args:
        wave_field: WaveField (grid dimensions)
        observables: FieldObservables (reads energyH/energyF per-voxel fields;
            writes the two global aggregates)
    """
    global _slice_xy_energyH, _slice_xy_energyF
    global _slice_xz_energyH, _slice_xz_energyF
    global _slice_yz_energyH, _slice_yz_energyF

    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    # Initialize slice buffers once
    if _slice_xy_energyH is None:
        _slice_xy_energyH = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_energyF = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xz_energyH = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_energyF = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_yz_energyH = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_energyF = ti.field(dtype=ti.f32, shape=(ny, nz))

    # Copy 3 center slices to 2D buffers (parallel kernels)
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    _copy_slice_xy_observables(observables, _slice_xy_energyH, _slice_xy_energyF, mid_z)
    _copy_slice_xz_observables(observables, _slice_xz_energyH, _slice_xz_energyF, mid_y)
    _copy_slice_yz_observables(observables, _slice_yz_energyH, _slice_yz_energyF, mid_x)

    # Transfer 2D slices to CPU for numpy operations (exclude boundary voxels)
    xy_energyH = _slice_xy_energyH.to_numpy()[1:-1, 1:-1]
    xy_energyF = _slice_xy_energyF.to_numpy()[1:-1, 1:-1]
    xz_energyH = _slice_xz_energyH.to_numpy()[1:-1, 1:-1]
    xz_energyF = _slice_xz_energyF.to_numpy()[1:-1, 1:-1]
    yz_energyH = _slice_yz_energyH.to_numpy()[1:-1, 1:-1]
    yz_energyF = _slice_yz_energyF.to_numpy()[1:-1, 1:-1]

    n_samples = xy_energyH.size + xz_energyH.size + yz_energyH.size

    # Plain mean for energy density H (per voxel); launcher multiplies by
    # voxel_count to get the grid total.
    total_energyH = xy_energyH.sum() + xz_energyH.sum() + yz_energyH.sum()
    observables.energyH_global_avg_aJ[None] = float(total_energyH / n_samples)

    # Plain mean for Frank elastic energy density F (per voxel).
    total_energyF = xy_energyF.sum() + xz_energyF.sum() + yz_energyF.sum()
    observables.energyF_global_avg_aJ[None] = float(total_energyF / n_samples)


# ================================================================
# FLUX MESH VALUES UPDATING
# ================================================================
# The flux mesh is a VISUALIZATION layer — it converts simulation-side
# scalars/vectors to per-vertex colors and Z-axis warps so the user can "see"
# what the field is doing. It is not physics; nothing here feeds back into
# evolve_psi or any tracker. Treat it as a display driver.
#
# VECTOR-FIELD RENDERING NUANCE (worth remembering when M5.0g+ adds new
# wave_menu options for energy density, curl, divergence, energy flux,
# etc.):
#
# ψ is a Vector(3) field, but the flux mesh maps every voxel to ONE scalar
# (one color, one Z-warp height). The current "Displacement (Magnitude)"
# mode renders |ψ| = ψ.norm(). For a y-polarized seed ψ_y = A·sin(k·x):
#
#     |ψ| = A·|sin(k·x)|         period λ/2 — looks like 2× as many bumps
#     ψ_y = A· sin(k·x)          period λ — the "true" sinusoid (signed)
#
# Magnitude is the right scalar for energy-density-style observables (always
# positive, ∝ |ψ|²) but rectifies signed waves and visually doubles their
# spatial frequency. When designing future wave_menu modes, decide per
# observable:
#
#   - magnitude: |v|              → energy/intensity views (always ≥ 0)
#   - signed component: v · ê     → propagation/polarization views (signed)
#   - radial component: v · r̂    → longitudinal vs transverse decomposition
#   - dominant axis: max(|v_i|)·sign  → quick-and-dirty signed view
#
# For the flux mesh specifically (one scalar per vertex), each new wave_menu
# entry should pick one explicitly and not silently fall back to .norm().
# This becomes important once we have curl(ψ), force vectors, energy flux,
# etc. — all of which are vector quantities that need a deliberate scalar
# projection to be rendered.


@ti.kernel
def update_flux_mesh_values(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    wave_menu: ti.i32,  # type: ignore
    warp_mesh: ti.i32,  # type: ignore
):
    """
    Update flux mesh colors and vertices by sampling wave properties from voxel grid.

    Samples wave displacement at each plane vertex position and maps it to a color.
    Should be called every frame after wave propagation to update visualization.

    Args:
        wave_field: WaveField instance with flux mesh fields and ψ data
        trackers: Trackers — reads amp/freq (WAVE_MENU 2 / 3 + colormap range)
        observables: FieldObservables — reads energyH/energyF (WAVE_MENU 4 / 5)
        wave_menu: Selected Wave displayed with color palette
    """

    # ================================================================
    # XY Plane: Sample at z = fm_plane_z_idx
    # ================================================================
    # Always update all planes (conditionals cause GPU branch divergence)
    # wave_menu == 4 renders the energy density (Hamiltonian) field
    # `observables.energyH_density_aJ`, populated each step by M5.0g.
    # wave_menu == 5 renders the Frank elastic density `observables.energyF_density_aJ`,
    # populated by compute_energyF_density (M5.1 task 5).
    for i, j in ti.ndrange(wave_field.nx, wave_field.ny):
        # Sample displacement at this voxel
        disp_value = wave_field.psi_am[i, j, wave_field.fm_plane_z_idx]
        amp_value = trackers.amp_local_emarms_am[i, j, wave_field.fm_plane_z_idx]
        freq_value = trackers.freq_local_cross_rHz[i, j, wave_field.fm_plane_z_idx]
        energyH_value = observables.energyH_density_aJ[i, j, wave_field.fm_plane_z_idx]
        energyF_value = observables.energyF_density_aJ[i, j, wave_field.fm_plane_z_idx]
        univ_edge_z = wave_field.universe_size_am[2]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 5:  # Frank elastic density on ironbow (defect-focused palette)
            F_max = observables.energyF_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                energyF_value, 0.0, F_max
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                energyF_value / F_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size))
            )
        elif wave_menu == 4:  # Energy density (Hamiltonian) on ironbow
            H_max = observables.energyH_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                energyH_value, 0.0, H_max
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                energyH_value / H_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size))
            )
        elif wave_menu == 3:  # blueprint
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + wave_field.flux_mesh_planes[2] * (
                wave_field.nz / wave_field.max_grid_size
            )
        elif wave_menu == 2:  # ironbow
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                amp_value / univ_edge_z * warp_mesh
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )
        else:  # default to orange (wave_menu == 1 or 4)
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_orange_color(
                disp_value.norm(), 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_xy_vertices[i, j] = ti.Vector(
                [
                    (ti.cast(i, ti.f32) + 0.5) / wave_field.max_grid_size
                    + disp_value[0] / wave_field.universe_size_am[0] * warp_mesh,
                    (ti.cast(j, ti.f32) + 0.5) / wave_field.max_grid_size
                    + disp_value[1] / wave_field.universe_size_am[1] * warp_mesh,
                    wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
                    + disp_value[2] / wave_field.universe_size_am[2] * warp_mesh,
                ]
            )

    # ================================================================
    # XZ Plane: Sample at y = fm_plane_y_idx
    # ================================================================
    for i, k in ti.ndrange(wave_field.nx, wave_field.nz):
        # Sample displacement at this voxel
        disp_value = wave_field.psi_am[i, wave_field.fm_plane_y_idx, k]
        amp_value = trackers.amp_local_emarms_am[i, wave_field.fm_plane_y_idx, k]
        freq_value = trackers.freq_local_cross_rHz[i, wave_field.fm_plane_y_idx, k]
        energyH_value = observables.energyH_density_aJ[i, wave_field.fm_plane_y_idx, k]
        energyF_value = observables.energyF_density_aJ[i, wave_field.fm_plane_y_idx, k]
        univ_edge_y = wave_field.universe_size_am[1]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 5:  # Frank elastic density on ironbow
            F_max = observables.energyF_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                energyF_value, 0.0, F_max
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                energyF_value / F_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size))
            )
        elif wave_menu == 4:  # Energy density (Hamiltonian) on ironbow
            H_max = observables.energyH_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                energyH_value, 0.0, H_max
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                energyH_value / H_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size))
            )
        elif wave_menu == 3:  # blueprint
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + wave_field.flux_mesh_planes[1] * (
                wave_field.ny / wave_field.max_grid_size
            )
        elif wave_menu == 2:  # ironbow
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                amp_value / univ_edge_y * warp_mesh
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )
        else:  # default to orange (wave_menu == 1 or 4)
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_orange_color(
                disp_value.norm(), 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_xz_vertices[i, k] = ti.Vector(
                [
                    (ti.cast(i, ti.f32) + 0.5) / wave_field.max_grid_size
                    + disp_value[0] / wave_field.universe_size_am[0] * warp_mesh,
                    wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
                    + disp_value[1] / wave_field.universe_size_am[1] * warp_mesh,
                    (ti.cast(k, ti.f32) + 0.5) / wave_field.max_grid_size
                    + disp_value[2] / wave_field.universe_size_am[2] * warp_mesh,
                ]
            )

    # ================================================================
    # YZ Plane: Sample at x = fm_plane_x_idx
    # ================================================================
    for j, k in ti.ndrange(wave_field.ny, wave_field.nz):
        # Sample displacement at this voxel
        disp_value = wave_field.psi_am[wave_field.fm_plane_x_idx, j, k]
        amp_value = trackers.amp_local_emarms_am[wave_field.fm_plane_x_idx, j, k]
        freq_value = trackers.freq_local_cross_rHz[wave_field.fm_plane_x_idx, j, k]
        energyH_value = observables.energyH_density_aJ[wave_field.fm_plane_x_idx, j, k]
        energyF_value = observables.energyF_density_aJ[wave_field.fm_plane_x_idx, j, k]
        univ_edge_x = wave_field.universe_size_am[0]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 5:  # Frank elastic density on ironbow
            F_max = observables.energyF_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                energyF_value, 0.0, F_max
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                energyF_value / F_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size))
            )
        elif wave_menu == 4:  # Energy density (Hamiltonian) on ironbow
            H_max = observables.energyH_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                energyH_value, 0.0, H_max
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                energyH_value / H_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size))
            )
        elif wave_menu == 3:  # blueprint
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + wave_field.flux_mesh_planes[0] * (
                wave_field.nx / wave_field.max_grid_size
            )
        elif wave_menu == 2:  # ironbow
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                amp_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
        else:  # default to orange (wave_menu == 1 or 4)
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_orange_color(
                disp_value.norm(), 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_yz_vertices[j, k] = ti.Vector(
                [
                    wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
                    + disp_value[0] / wave_field.universe_size_am[0] * warp_mesh,
                    (ti.cast(j, ti.f32) + 0.5) / wave_field.max_grid_size
                    + disp_value[1] / wave_field.universe_size_am[1] * warp_mesh,
                    (ti.cast(k, ti.f32) + 0.5) / wave_field.max_grid_size
                    + disp_value[2] / wave_field.universe_size_am[2] * warp_mesh,
                ]
            )


# ================================================================
# DIRECTOR-GLYPH RENDERING (M5.1)
# ================================================================
# The flux mesh maps every voxel to ONE scalar (color + Z-warp). For a
# director field where |ψ|=1 by construction, scalar magnitude is uninformative
# (always 1 except at singular cores). To see TOPOLOGY, we need to render
# the *direction* — which is what director glyphs do.
#
# Each glyph is a line segment from voxel position to voxel + L · n̂, with
# both endpoints sharing a color computed from the signed components:
#     color = (n̂ + 1) / 2     ∈ [0, 1]³
#
# This signed-component RGB has a useful property: opposite directions get
# RGB-complementary colors (red↔cyan, green↔magenta, blue↔yellow). A +1
# hedgehog shows red on the +x face and cyan on the −x face; a −1 hedgehog
# inverts the pattern. Combined with the line offset (glyphs reach further
# in the n̂ direction), polarity is unambiguous without arrowheads.
#
# Sampling pattern: same 3-plane convention as flux_mesh — XY at fm_plane_z,
# XZ at fm_plane_y, YZ at fm_plane_x. At GLYPH_STRIDE=4 with a 64³ grid,
# total glyph count is ~768 segments; renders in well under 1 ms.
#
# Design doc: research_hub/3e_director_glyph_rendering.md.

# Half-arrowhead toggle (hardcoded for testing — flip to False to compare):
# When True, each glyph gets ONE extra barb at the tip end so head vs tail is
# unambiguous without arrowhead V-shapes. Cost: doubles the line count rendered
# when on (e.g. 768 → 1536 lines at 64³ / stride=4 — still trivial on GPU).
# Buffer is always allocated; the launcher gates the second scene.lines call.
SHOW_DIRECTOR_ARROWHEAD = True
ARROWHEAD_LENGTH_FRAC = 0.35  # barb length relative to shaft length
# Barb angle measured from −n̂ (the backward shaft direction), rotated toward
# the perpendicular axis. NOTE: smaller values give THINNER arrow-head profiles
# (the variable runs opposite to a "sweep from perpendicular" intuition):
#   90° → perpendicular barb (widest, no backward sweep)
#   36° → moderately swept back (thin arrow profile — current default)
#   <30° → very thin (barb collapses toward −n̂)
ARROWHEAD_ANGLE_DEG = 36.0
_arrow_rad = math.radians(ARROWHEAD_ANGLE_DEG)
ARROWHEAD_BACK_COMP = -math.cos(_arrow_rad)  # component along n̂ (negative = backward)
ARROWHEAD_PERP_COMP = math.sin(_arrow_rad)  # component along perpendicular axis


@ti.kernel
def update_director_glyphs(
    wave_field: ti.template(),  # type: ignore
    length: ti.f32,  # type: ignore
    arrow_length: ti.f32,  # type: ignore
    show_level: ti.i32,  # type: ignore
):
    """
    Update the director-glyph line-segment field by sampling ψ on the
    three flux-mesh planes (XY, XZ, YZ) at GLYPH_STRIDE.

    Per glyph, writes:
      - vertices[2k+0] = voxel position (normalized [0,1] camera space)
      - vertices[2k+1] = voxel position + length · n̂
      - colors[2k+0]   = colors[2k+1] = colormap palette of (1 − n̂[2]) ∈ [0, 2]
        scalar represents "twist away from vacuum": 0 at n=ẑ (vacuum), 1 at
        equatorial directors, 2 at the −ẑ pole. With dark-to-bright palettes
        (blueprint default; ironbow alt) the vacuum BLENDS into the black GUI
        background — only the defect-twisted region stands out visually.

    Half-arrowhead barb (when arrow_length > 0):
      - arrow_vertices[2k+0] = tip = pos + length · n̂
      - arrow_vertices[2k+1] = tip + arrow_length · barb_dir
      - barb_dir = ARROWHEAD_BACK_COMP · n̂ + ARROWHEAD_PERP_COMP · perp,
        where perp is a stable perpendicular to n̂ (cross with ẑ, falling
        back to x̂ when n̂ is nearly parallel to ẑ) and the back/perp
        components encode ARROWHEAD_ANGLE_DEG (54° default → swept-back
        arrowhead, ≈36° between barb and shaft). Same color as the shaft.
        The launcher decides whether to render this buffer.

    `show_level` mirrors SHOW_FLUX_MESH semantics for progressive plane reveal:
        0 → all planes off (degenerate glyphs; nothing renders)
        1 → XY plane only
        2 → XY + XZ planes
        3 → XY + XZ + YZ planes (all)
    Off-planes get degenerate glyphs (start == end, color = 0) which GGUI
    renders as invisible 0-length line segments.

    Boundary handling: 1e-12 epsilon in n̂ denominator avoids /0 at vacuum
    voxels where ψ ≈ 0 (none in M5.1's seeded fields, but guarding anyway).

    Args:
        wave_field: WaveField (reads psi_am; writes director_glyph_vertices /
            director_glyph_colors / director_glyph_arrow_vertices /
            director_glyph_arrow_colors)
        length: glyph shaft length in normalized camera coords (e.g. 0.02 for
            ~2% of the universe edge)
        arrow_length: barb length in the same normalized coords. Pass 0.0 to
            skip computing the barb (writes zero-length lines instead).
        show_level: 0..3, parallel to SHOW_FLUX_MESH (0 off, 1 XY, 2 +XZ, 3 all)
    """
    nx = wave_field.nx
    ny = wave_field.ny
    nz = wave_field.nz
    stride = wave_field.GLYPH_STRIDE
    max_dim = ti.cast(wave_field.max_grid_size, ti.f32)

    # Use the pre-computed round-up sampled counts (matches granule sampling
    # pattern in _launcher.py — last partial row is included).
    nx_s = wave_field.glyph_nx_s
    ny_s = wave_field.glyph_ny_s
    nz_s = wave_field.glyph_nz_s

    z_plane_idx = wave_field.fm_plane_z_idx
    y_plane_idx = wave_field.fm_plane_y_idx
    x_plane_idx = wave_field.fm_plane_x_idx

    zero_v = ti.Vector([0.0, 0.0, 0.0])

    # ----- XY plane at z = fm_plane_z_idx -----
    for si, sj in ti.ndrange(nx_s, ny_s):
        idx = (wave_field.glyph_offset_xy + si * ny_s + sj) * 2
        if show_level >= 1:
            # Clamp to grid extent — round-up sampling can give
            # (nx_s − 1) * stride > nx − 1 for some grid/stride combos.
            i = ti.min(si * stride, nx - 1)
            j = ti.min(sj * stride, ny - 1)
            psi = wave_field.psi_am[i, j, z_plane_idx]
            norm = psi.norm() + 1e-12
            n_hat = psi / norm

            pos = ti.Vector(
                [
                    ti.cast(i, ti.f32) / max_dim,
                    ti.cast(j, ti.f32) / max_dim,
                    ti.cast(z_plane_idx, ti.f32) / max_dim,
                ]
            )
            tip = pos + length * n_hat
            # Color = palette mapping of (1 − n_hat[2]) ∈ [0, 2]:
            # vacuum (n=ẑ) → 0 (DARK, blends into black GUI background, defect
            # is what stands out); equator (n in xy-plane) → 1 (mid); inward
            # south-pole (n=-ẑ) → 2 (BRIGHT, peak twist away from vacuum).
            # CHANGE THIS LINE to test palettes (both start dark, end bright):
            #   colormap.get_ironbow_color(...)    black → magenta → red → yellow-white
            #   colormap.get_blueprint_color(...)  dark blue → light blue (current)
            color = colormap.get_blueprint_color(1.0 - n_hat[2], 0.0, 2.0)
            wave_field.director_glyph_vertices[idx + 0] = pos
            wave_field.director_glyph_vertices[idx + 1] = tip
            wave_field.director_glyph_colors[idx + 0] = color
            wave_field.director_glyph_colors[idx + 1] = color

            # Half-arrowhead: barb perpendicular to n̂ at the tip end.
            # Stable perp = n̂ × ref where ref flips to x̂ when n̂ ≈ ẑ
            # to avoid the degenerate cross product at vacuum directors.
            ref = ti.Vector([0.0, 0.0, 1.0])
            if ti.abs(n_hat[2]) > 0.9:
                ref = ti.Vector([1.0, 0.0, 0.0])
            perp = n_hat.cross(ref).normalized()
            barb_dir = ARROWHEAD_BACK_COMP * n_hat + ARROWHEAD_PERP_COMP * perp
            barb_end = tip + arrow_length * barb_dir
            wave_field.director_glyph_arrow_vertices[idx + 0] = tip
            wave_field.director_glyph_arrow_vertices[idx + 1] = barb_end
            wave_field.director_glyph_arrow_colors[idx + 0] = color
            wave_field.director_glyph_arrow_colors[idx + 1] = color
        else:
            wave_field.director_glyph_vertices[idx + 0] = zero_v
            wave_field.director_glyph_vertices[idx + 1] = zero_v
            wave_field.director_glyph_colors[idx + 0] = zero_v
            wave_field.director_glyph_colors[idx + 1] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 0] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 1] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 0] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 1] = zero_v

    # ----- XZ plane at y = fm_plane_y_idx -----
    for si, sk in ti.ndrange(nx_s, nz_s):
        idx = (wave_field.glyph_offset_xz + si * nz_s + sk) * 2
        if show_level >= 2:
            i = ti.min(si * stride, nx - 1)
            k = ti.min(sk * stride, nz - 1)
            psi = wave_field.psi_am[i, y_plane_idx, k]
            norm = psi.norm() + 1e-12
            n_hat = psi / norm

            pos = ti.Vector(
                [
                    ti.cast(i, ti.f32) / max_dim,
                    ti.cast(y_plane_idx, ti.f32) / max_dim,
                    ti.cast(k, ti.f32) / max_dim,
                ]
            )
            tip = pos + length * n_hat
            # Color = palette mapping of (1 − n_hat[2]) ∈ [0, 2]:
            # vacuum (n=ẑ) → 0 (DARK, blends into black GUI background, defect
            # is what stands out); equator (n in xy-plane) → 1 (mid); inward
            # south-pole (n=-ẑ) → 2 (BRIGHT, peak twist away from vacuum).
            # CHANGE THIS LINE to test palettes (both start dark, end bright):
            #   colormap.get_ironbow_color(...)    black → magenta → red → yellow-white
            #   colormap.get_blueprint_color(...)  dark blue → light blue (current)
            color = colormap.get_blueprint_color(1.0 - n_hat[2], 0.0, 2.0)
            wave_field.director_glyph_vertices[idx + 0] = pos
            wave_field.director_glyph_vertices[idx + 1] = tip
            wave_field.director_glyph_colors[idx + 0] = color
            wave_field.director_glyph_colors[idx + 1] = color

            ref = ti.Vector([0.0, 0.0, 1.0])
            if ti.abs(n_hat[2]) > 0.9:
                ref = ti.Vector([1.0, 0.0, 0.0])
            perp = n_hat.cross(ref).normalized()
            barb_dir = ARROWHEAD_BACK_COMP * n_hat + ARROWHEAD_PERP_COMP * perp
            barb_end = tip + arrow_length * barb_dir
            wave_field.director_glyph_arrow_vertices[idx + 0] = tip
            wave_field.director_glyph_arrow_vertices[idx + 1] = barb_end
            wave_field.director_glyph_arrow_colors[idx + 0] = color
            wave_field.director_glyph_arrow_colors[idx + 1] = color
        else:
            wave_field.director_glyph_vertices[idx + 0] = zero_v
            wave_field.director_glyph_vertices[idx + 1] = zero_v
            wave_field.director_glyph_colors[idx + 0] = zero_v
            wave_field.director_glyph_colors[idx + 1] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 0] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 1] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 0] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 1] = zero_v

    # ----- YZ plane at x = fm_plane_x_idx -----
    for sj, sk in ti.ndrange(ny_s, nz_s):
        idx = (wave_field.glyph_offset_yz + sj * nz_s + sk) * 2
        if show_level >= 3:
            j = ti.min(sj * stride, ny - 1)
            k = ti.min(sk * stride, nz - 1)
            psi = wave_field.psi_am[x_plane_idx, j, k]
            norm = psi.norm() + 1e-12
            n_hat = psi / norm

            pos = ti.Vector(
                [
                    ti.cast(x_plane_idx, ti.f32) / max_dim,
                    ti.cast(j, ti.f32) / max_dim,
                    ti.cast(k, ti.f32) / max_dim,
                ]
            )
            tip = pos + length * n_hat
            # Color = palette mapping of (1 − n_hat[2]) ∈ [0, 2]:
            # vacuum (n=ẑ) → 0 (DARK, blends into black GUI background, defect
            # is what stands out); equator (n in xy-plane) → 1 (mid); inward
            # south-pole (n=-ẑ) → 2 (BRIGHT, peak twist away from vacuum).
            # CHANGE THIS LINE to test palettes (both start dark, end bright):
            #   colormap.get_ironbow_color(...)    black → magenta → red → yellow-white
            #   colormap.get_blueprint_color(...)  dark blue → light blue (current)
            color = colormap.get_blueprint_color(1.0 - n_hat[2], 0.0, 2.0)
            wave_field.director_glyph_vertices[idx + 0] = pos
            wave_field.director_glyph_vertices[idx + 1] = tip
            wave_field.director_glyph_colors[idx + 0] = color
            wave_field.director_glyph_colors[idx + 1] = color

            ref = ti.Vector([0.0, 0.0, 1.0])
            if ti.abs(n_hat[2]) > 0.9:
                ref = ti.Vector([1.0, 0.0, 0.0])
            perp = n_hat.cross(ref).normalized()
            barb_dir = ARROWHEAD_BACK_COMP * n_hat + ARROWHEAD_PERP_COMP * perp
            barb_end = tip + arrow_length * barb_dir
            wave_field.director_glyph_arrow_vertices[idx + 0] = tip
            wave_field.director_glyph_arrow_vertices[idx + 1] = barb_end
            wave_field.director_glyph_arrow_colors[idx + 0] = color
            wave_field.director_glyph_arrow_colors[idx + 1] = color
        else:
            wave_field.director_glyph_vertices[idx + 0] = zero_v
            wave_field.director_glyph_vertices[idx + 1] = zero_v
            wave_field.director_glyph_colors[idx + 0] = zero_v
            wave_field.director_glyph_colors[idx + 1] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 0] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 1] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 0] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 1] = zero_v
