"""
M5 ENGINE — INITIAL-CONDITION SEEDING (engine1_seeds)

Seed kernels that prime the leapfrog triple-buffer (psi_prev, psi, psi_new)
with a known field configuration. Split from the former monolithic
lagrangian_engine.py (2026-05-25 SoC refactor). Concern order:
  engine1_seeds → engine2_pde → engine2_pde →
  engine3_observables → engine4_render
"""

import taichi as ti

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
    `research/scripts/sandbox_lagrangian/exp2_hedgehog_energy.py:71-108`),
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
# MATRIX-FIELD SEEDERS (M5.4 — order parameter M = O·D·O^T)
# ================================================================
# Matrix-substrate analogs of seed_vacuum / seed_hedgehog. The director n̂ is
# computed exactly as in the ψ versions (same blend, same soft-core), then
# embedded in the uniaxial order parameter
#     M = δ·I + (1−δ)·n̂⊗n̂        (eigenvalues 1, δ, δ; principal eigvec = n̂)
# Each seeder writes BOTH the matrix triple buffer (M_prev/M/M_new — triple-buffer
# BC rule, so the matrix leapfrog landing in M5.5 starts BC-consistent) AND the
# derived director_nhat with the correct PHYSICAL sign. Seeding director_nhat
# directly keeps eigen_decompose's sign-continuity branch smooth from t=0 and lets
# the M5.4 director-equivalent Coulomb relaxation operate on a clean director.


@ti.func
def uniaxial_M(n, delta):  # type: ignore
    """Uniaxial order parameter M = δ·I + (1−δ)·n̂⊗n̂ from a unit director n̂.

    Eigenvalues (1, δ, δ), principal eigenvector = n̂ — the M5.4 embedding of a
    director into the matrix substrate. Verified by the eigen_decompose round-trip
    (spectrum recovers (1, δ, δ); director recovers n̂ at 0.9995).
    """
    eye = ti.Matrix.identity(ti.f32, 3)
    return delta * eye + (1.0 - delta) * n.outer_product(n)


@ti.kernel
def seed_vacuum_M(
    wave_field: ti.template(),  # type: ignore
    delta: ti.f32,  # type: ignore
):
    """Fill the matrix field with the topological-vacuum ground state (n̂ = ẑ).

    M_vac = δ·I + (1−δ)·ẑ⊗ẑ everywhere (uniform orientation = no defects). Writes
    all three matrix buffers + director_nhat (= ẑ). Matrix analog of seed_vacuum.

    Args:
        wave_field: WaveField (writes M_am, M_prev_am, M_new_am, director_nhat)
        delta: uniaxial minor-axis eigenvalue (wave_field.lc_delta)
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    z_hat = ti.Vector([0.0, 0.0, 1.0])
    m_vac = uniaxial_M(z_hat, delta)
    for i, j, k in ti.ndrange(nx, ny, nz):
        wave_field.M_am[i, j, k] = m_vac
        wave_field.M_prev_am[i, j, k] = m_vac
        wave_field.M_new_am[i, j, k] = m_vac
        wave_field.director_nhat[i, j, k] = z_hat
        wave_field.director_nhat_new[i, j, k] = z_hat


@ti.kernel
def seed_hedgehog_M(
    wave_field: ti.template(),  # type: ignore
    centers_voxel: ti.template(),  # type: ignore
    signs: ti.template(),  # type: ignore
    domain_quarter_voxels: ti.f32,  # type: ignore
    n_defects: ti.i32,  # type: ignore
    delta: ti.f32,  # type: ignore
):
    """Seed N hedgehog defects on the matrix substrate.

    Identical director construction to seed_hedgehog (weighted radial superposition
    + vacuum blend + renormalization — see that kernel's docstring for the math and
    caveats), then embeds the unit director via uniaxial_M and writes the matrix
    triple buffer + director_nhat. The multi-defect API drives the M5.4 Coulomb-pair
    gate directly (centers ±d/2 along x, signs ±1).

    Args:
        wave_field: WaveField (writes M_am, M_prev_am, M_new_am, director_nhat)
        centers_voxel: ti.field (n_defects, 3) i32 — centers in voxel coords
        signs: ti.field (n_defects,) i32 — +1 outward / −1 inward per defect
        domain_quarter_voxels: f32 — D/4 in voxel units (w_vac falloff radius)
        n_defects: i32 — number of active defects
        delta: uniaxial minor-axis eigenvalue (wave_field.lc_delta)
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    z_hat = ti.Vector([0.0, 0.0, 1.0])
    soft_core_voxels = ti.cast(0.2, ti.f32)
    proximity_floor = ti.cast(0.5, ti.f32)

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
                [sgn * ddx / r_safe, sgn * ddy / r_safe, sgn * ddz / r_safe]
            )
            w_defect = 1.0 / (r + proximity_floor)
            n_combined += w_defect * n_radial
            r_nearest = ti.min(r_nearest, r)

        ratio = r_nearest / domain_quarter_voxels
        w_vac = 1.0 / (1.0 + ratio * ratio * ratio * ratio)
        n_blended = w_vac * n_combined + (1.0 - w_vac) * z_hat
        n_unit = n_blended / (n_blended.norm() + 1e-12)

        m = uniaxial_M(n_unit, delta)
        wave_field.M_am[i, j, k] = m
        wave_field.M_prev_am[i, j, k] = m
        wave_field.M_new_am[i, j, k] = m
        wave_field.director_nhat[i, j, k] = n_unit
        wave_field.director_nhat_new[i, j, k] = n_unit


