"""
LAGRANGIAN-FIELD ENGINE (M5)

PDE-driven evolution of the field ψ on a voxel grid. The action is

    S[ψ] = ∫ dt d³x [ ½|ψ̇|² − ½c²|∇ψ|² − V(ψ) ]

Euler-Lagrange yields  ∂²ψ/∂t² = c²·∇²ψ − ∂V/∂ψ.  In M5.0d the potential is
zero (free wave); M5.2 plugs in Klein-Gordon mass + Close Eq. 23 + LdG.

Module layout (top-to-bottom):
    DIFFERENTIAL OPERATORS    — Laplacian; M5.0e adds curl/div/curl-curl
    INITIAL-CONDITION SEEDING — seed_wave (test/UI verification)
    ψ PROPAGATION ENGINE      — propagate_psi (leapfrog/Verlet)
    FIELD OBSERVABLES         — update_trackers_psi (amp/freq EMA)
    TOTAL HAMILTONIAN         — compute_total_hamiltonian (scalar reduction)
    POSITION RENDER           — sample_position_to_render (granule viz)
    3-PLANE SAMPLING          — sample_avg_trackers (cheap global means)
    FLUX MESH UPDATING        — update_flux_mesh_values (color mapping)
"""

import taichi as ti
import numpy as np

from openwave.common import colormap

# ================================================================
# INITIAL-CONDITION SEEDING
# ================================================================
# Seed kernels prime the triple-buffer (psi_prev, psi at t−dt and t) with a
# known analytical solution so propagate_psi can step it forward. Used by:
#   - test/UI xperiments to verify the leapfrog kernel visually
#   - M5.0h dispersion-relation regression test


@ti.kernel
def seed_wave(
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


# ================================================================
# DIFFERENTIAL OPERATORS
# ================================================================
# 6-point Laplacian stencil for the vector field ψ. Acts component-wise on
# psi_am (ti.Vector.field(3, ...)); Taichi handles Vector(3) arithmetic
# natively, so the stencil applied to a Vector field returns a Vector field.
#
# Curl, divergence, and curl(curl) operators land in this section in M5.0e.


@ti.func
def compute_laplacian_psi(
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


# ================================================================
# ψ PROPAGATION ENGINE — LAGRANGIAN FIELD EVOLUTION
# ================================================================
# propagate_psi evolves the field ψ via leapfrog/Verlet integration of the
# wave equation:
#     ∂²ψ/∂t² = c²·∇²ψ − ∂V/∂ψ
#
# In M5.0d this lands with V(ψ) = 0 (free wave). The nonlinear potential term
# −∂V/∂ψ is added in M5.2 (Close's Eq. 23 + Klein-Gordon mass term + LdG).
#
# Why "propagate_psi" not "propagate_wave": the operation is field evolution,
# not wave-specific. Leapfrog also evolves topology (defect drift), gradient
# descent for relaxation reuses the same buffers, etc. Naming follows what the
# function operates on (the field ψ), not the method (leapfrog) or one of the
# behaviors (wave propagation).
#
# Caller MUST call wave_field.swap_buffers() after this kernel returns, to
# cycle the triple buffer: psi_prev ← psi, psi ← psi_new.


@ti.kernel
def propagate_psi(
    wave_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
):
    """
    Evolve ψ one timestep via leapfrog/Verlet integration.

    Wave Equation (free wave, V=0):
        ∂²ψ/∂t² = c²·∇²ψ

    Discrete (leapfrog/Verlet):
        ψ_new = 2·ψ − ψ_prev + (c·dt)²·∇²ψ

    Reads:  wave_field.psi_am, wave_field.psi_prev_am
    Writes: wave_field.psi_new_am

    Boundary: Dirichlet (ψ = 0 at edges) — interior voxels only.
    The 6-point Laplacian needs a 1-cell halo, so the loop range is
    (1, n−1) on each axis. Boundary voxels are never updated; they stay
    at whatever value they were initialized with (0 by default).

    CFL stability: dt ≤ dx / (c·√3) for 3D wave equation. Caller is
    responsible for sizing dt below this bound (see _launcher CFL eval).

    Args:
        wave_field: WaveField instance (reads psi_am, psi_prev_am; writes psi_new_am)
        c_amrs: wave speed in scaled units (am/rs); already includes any
            slow-motion factor for visualization
        dt_rs: timestep in rontoseconds, sized below the CFL bound
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    c_dt_squared = (c_amrs * dt_rs) ** 2

    # Interior voxels only (Dirichlet BC; 6-point Laplacian needs 1-cell halo)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        laplacian_psi_am = compute_laplacian_psi(wave_field, i, j, k)

        # Leapfrog/Verlet update: ψ_new = 2·ψ − ψ_prev + (c·dt)²·∇²ψ
        wave_field.psi_new_am[i, j, k] = (
            2.0 * wave_field.psi_am[i, j, k]
            - wave_field.psi_prev_am[i, j, k]
            + c_dt_squared * laplacian_psi_am
        )


# ================================================================
# FIELD OBSERVABLES — PER-VOXEL AMPLITUDE & FREQUENCY TRACKERS
# ================================================================
# Trackers run after propagate_psi each step. They derive observable scalars
# from ψ (the field) without touching it — read-only on psi_am, write to the
# tracker fields in-place. EMA / zero-crossing logic preserved from M4 but
# adapted to vector ψ:
#   - amp_local_emarms_am  : RMS of |ψ| over recent history (EMA on |ψ|²)
#   - freq_local_cross_rHz : zero-crossing rate of ψ_z (chosen polarization axis)
#
# PERFORMANCE NOTE (deferred to M5.0i): this is a SECOND full-grid pass over
# ψ — propagate_psi already touched every voxel; we re-stream ψ and ψ_prev
# from memory just to compute trackers. At 100M voxels this costs ~5 GB of
# extra memory bandwidth per step (~25–30% of total per-step traffic).
#
# M2 avoided this by computing trackers INSIDE the same loop as the wave
# update, but M2's wave kernel was analytical (no neighbor dependency), so
# trackers could read the just-written ψ trivially. M5's leapfrog can do the
# same: in a merged kernel, write ψ_new then immediately use it (with
# ψ_prev already in registers) for tracker EMA / zero-crossing — gives true
# central-difference ψ̇ for free, more accurate than the forward diff used
# below (psi_new gets consumed by swap_buffers before update_trackers_psi
# can see it, so we fall back to (ψ − ψ_prev) / dt here).
#
# Splitting was deliberate for M5.0d: cleaner kernel boundaries, AMR-ready
# (octree retrofit in M5.6/M5.8 prefers single-purpose kernels), and trackers
# are skippable in tests (e.g., M5.0h dispersion regression). Profile-first
# in M5.0i — alongside swap_buffers' rotating-pointer optimization (see
# medium.py:swap_buffers docstring), this is the second-largest known
# bandwidth target on the per-step path.


@ti.kernel
def update_trackers_psi(
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
# TOTAL HAMILTONIAN — DASHBOARD ENERGY OBSERVABLE
# ================================================================
# H = ½|ψ̇|² + ½c²|∇ψ|² + V(ψ)  ;  V=0 in M5.0d.
# Estimated via 3-plane center-slice sampling — same compromise as
# sample_avg_trackers, for the same reason: a full 3D atomic reduction on a
# 100M-voxel grid contends so badly it can stall the GUI for tens of seconds
# (looks like a freeze). 3-plane sampling reads ~3·N² voxels and is exact for
# spatially-uniform fields and a good estimator for smooth wave packets.
# Per-voxel H density field + force coupling land in M5.0g.


@ti.kernel
def _hamiltonian_slice_xy(
    wave_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    H_slice: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """H per voxel on the XY center plane (fixed z = mid_z). Excludes boundary."""
    nx, ny = wave_field.nx, wave_field.ny
    inv_dt = 1.0 / dt_rs
    half_c2 = 0.5 * c_amrs**2

    for i, j in ti.ndrange((1, nx - 1), (1, ny - 1)):
        psi_dot = (wave_field.psi_am[i, j, mid_z] - wave_field.psi_prev_am[i, j, mid_z]) * inv_dt
        kinetic = 0.5 * psi_dot.norm_sqr()
        d_dx = (wave_field.psi_am[i + 1, j, mid_z] - wave_field.psi_am[i - 1, j, mid_z]) / (
            2.0 * wave_field.dx_am
        )
        d_dy = (wave_field.psi_am[i, j + 1, mid_z] - wave_field.psi_am[i, j - 1, mid_z]) / (
            2.0 * wave_field.dx_am
        )
        d_dz = (wave_field.psi_am[i, j, mid_z + 1] - wave_field.psi_am[i, j, mid_z - 1]) / (
            2.0 * wave_field.dx_am
        )
        grad_sq = d_dx.norm_sqr() + d_dy.norm_sqr() + d_dz.norm_sqr()
        H_slice[i, j] = kinetic + half_c2 * grad_sq


@ti.kernel
def _hamiltonian_slice_xz(
    wave_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    H_slice: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """H per voxel on the XZ center plane (fixed y = mid_y). Excludes boundary."""
    nx, nz = wave_field.nx, wave_field.nz
    inv_dt = 1.0 / dt_rs
    half_c2 = 0.5 * c_amrs**2

    for i, k in ti.ndrange((1, nx - 1), (1, nz - 1)):
        psi_dot = (wave_field.psi_am[i, mid_y, k] - wave_field.psi_prev_am[i, mid_y, k]) * inv_dt
        kinetic = 0.5 * psi_dot.norm_sqr()
        d_dx = (wave_field.psi_am[i + 1, mid_y, k] - wave_field.psi_am[i - 1, mid_y, k]) / (
            2.0 * wave_field.dx_am
        )
        d_dy = (wave_field.psi_am[i, mid_y + 1, k] - wave_field.psi_am[i, mid_y - 1, k]) / (
            2.0 * wave_field.dx_am
        )
        d_dz = (wave_field.psi_am[i, mid_y, k + 1] - wave_field.psi_am[i, mid_y, k - 1]) / (
            2.0 * wave_field.dx_am
        )
        grad_sq = d_dx.norm_sqr() + d_dy.norm_sqr() + d_dz.norm_sqr()
        H_slice[i, k] = kinetic + half_c2 * grad_sq


@ti.kernel
def _hamiltonian_slice_yz(
    wave_field: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    H_slice: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """H per voxel on the YZ center plane (fixed x = mid_x). Excludes boundary."""
    ny, nz = wave_field.ny, wave_field.nz
    inv_dt = 1.0 / dt_rs
    half_c2 = 0.5 * c_amrs**2

    for j, k in ti.ndrange((1, ny - 1), (1, nz - 1)):
        psi_dot = (wave_field.psi_am[mid_x, j, k] - wave_field.psi_prev_am[mid_x, j, k]) * inv_dt
        kinetic = 0.5 * psi_dot.norm_sqr()
        d_dx = (wave_field.psi_am[mid_x + 1, j, k] - wave_field.psi_am[mid_x - 1, j, k]) / (
            2.0 * wave_field.dx_am
        )
        d_dy = (wave_field.psi_am[mid_x, j + 1, k] - wave_field.psi_am[mid_x, j - 1, k]) / (
            2.0 * wave_field.dx_am
        )
        d_dz = (wave_field.psi_am[mid_x, j, k + 1] - wave_field.psi_am[mid_x, j, k - 1]) / (
            2.0 * wave_field.dx_am
        )
        grad_sq = d_dx.norm_sqr() + d_dy.norm_sqr() + d_dz.norm_sqr()
        H_slice[j, k] = kinetic + half_c2 * grad_sq


# Module-level slice buffers, lazily allocated on first call
_h_slice_xy = None
_h_slice_xz = None
_h_slice_yz = None


def compute_total_hamiltonian(wave_field, c_amrs, dt_rs):
    """
    Estimated total Hamiltonian H = Σ (½|ψ̇|² + ½c²|∇ψ|²) via 3-plane sampling.

    Returns ⟨H⟩ × voxel_count where ⟨H⟩ is the mean per-voxel H over the three
    orthogonal center planes — exact for spatially-uniform fields, a sound
    estimator for smooth wave packets, and ~10⁴× cheaper than full reduction
    on million-voxel grids.

    Units: am²/rs² per voxel; physical scaling lands with M5.0f.
    """
    global _h_slice_xy, _h_slice_xz, _h_slice_yz

    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    if _h_slice_xy is None:
        _h_slice_xy = ti.field(dtype=ti.f32, shape=(nx, ny))
        _h_slice_xz = ti.field(dtype=ti.f32, shape=(nx, nz))
        _h_slice_yz = ti.field(dtype=ti.f32, shape=(ny, nz))

    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    _hamiltonian_slice_xy(wave_field, c_amrs, dt_rs, _h_slice_xy, mid_z)
    _hamiltonian_slice_xz(wave_field, c_amrs, dt_rs, _h_slice_xz, mid_y)
    _hamiltonian_slice_yz(wave_field, c_amrs, dt_rs, _h_slice_yz, mid_x)

    xy = _h_slice_xy.to_numpy()[1:-1, 1:-1]
    xz = _h_slice_xz.to_numpy()[1:-1, 1:-1]
    yz = _h_slice_yz.to_numpy()[1:-1, 1:-1]
    mean_H_per_voxel = float((xy.sum() + xz.sum() + yz.sum()) / (xy.size + xz.size + yz.size))
    return mean_H_per_voxel * wave_field.voxel_count


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
# ================================================================

# Cached slice buffers (initialized on first call)
_slice_xy_amp = None
_slice_xy_freq = None
_slice_xz_amp = None
_slice_xz_freq = None
_slice_yz_amp = None
_slice_yz_freq = None


@ti.kernel
def _copy_slice_xy(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """Copy center XY slice (fixed z) to 2D buffer."""
    for i, j in slice_amp:
        slice_amp[i, j] = trackers.amp_local_emarms_am[i, j, mid_z]
        slice_freq[i, j] = trackers.freq_local_cross_rHz[i, j, mid_z]


@ti.kernel
def _copy_slice_xz(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """Copy center XZ slice (fixed y) to 2D buffer."""
    for i, k in slice_amp:
        slice_amp[i, k] = trackers.amp_local_emarms_am[i, mid_y, k]
        slice_freq[i, k] = trackers.freq_local_cross_rHz[i, mid_y, k]


@ti.kernel
def _copy_slice_yz(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """Copy center YZ slice (fixed x) to 2D buffer."""
    for j, k in slice_amp:
        slice_amp[j, k] = trackers.amp_local_emarms_am[mid_x, j, k]
        slice_freq[j, k] = trackers.freq_local_cross_rHz[mid_x, j, k]


def sample_avg_trackers(
    wave_field,
    trackers,
):
    """
    Estimate RMS amplitude and average frequency by sampling 3 orthogonal planes.

    Samples XY, XZ, and YZ center slices to avoid full 3D reduction.
    This is a deliberate performance compromise - full GPU reduction with
    atomic operations causes severe contention with millions of voxels.

    For isotropic fields, center-plane sampling provides representative estimates.

    Args:
        wave_field: WaveField instance containing grid dimensions
        trackers: WaveTrackers instance with per-voxel and average fields
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
    _copy_slice_xy(trackers, _slice_xy_amp, _slice_xy_freq, mid_z)
    _copy_slice_xz(trackers, _slice_xz_amp, _slice_xz_freq, mid_y)
    _copy_slice_yz(trackers, _slice_yz_amp, _slice_yz_freq, mid_x)

    # Transfer 2D slices to CPU for numpy operations (exclude boundary voxels)
    xy_amp = _slice_xy_amp.to_numpy()[1:-1, 1:-1]
    xy_freq = _slice_xy_freq.to_numpy()[1:-1, 1:-1]
    xz_amp = _slice_xz_amp.to_numpy()[1:-1, 1:-1]
    xz_freq = _slice_xz_freq.to_numpy()[1:-1, 1:-1]
    yz_amp = _slice_yz_amp.to_numpy()[1:-1, 1:-1]
    yz_freq = _slice_yz_freq.to_numpy()[1:-1, 1:-1]

    # Compute RMS amplitude: √(⟨A²⟩) for correct energy weighting
    total_amp_squared = (xy_amp**2).sum() + (xz_amp**2).sum() + (yz_amp**2).sum()
    total_freq = xy_freq.sum() + xz_freq.sum() + yz_freq.sum()
    n_samples = xy_amp.size + xz_amp.size + yz_amp.size

    trackers.amp_global_emarms_am[None] = float(np.sqrt(total_amp_squared / n_samples))
    trackers.freq_global_avg_rHz[None] = float(total_freq / n_samples)


# ================================================================
# FLUX MESH VALUES UPDATING
# ================================================================
# The flux mesh is a VISUALIZATION layer — it converts simulation-side
# scalars/vectors to per-vertex colors and Z-axis warps so the user can "see"
# what the field is doing. It is not physics; nothing here feeds back into
# propagate_psi or any tracker. Treat it as a display driver.
#
# VECTOR-FIELD RENDERING NUANCE (worth remembering when M5.0g+ adds new
# wave_menu options for Hamiltonian density, curl, divergence, energy flux,
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
    wave_menu: ti.i32,  # type: ignore
    warp_mesh: ti.i32,  # type: ignore
):
    """
    Update flux mesh colors and vertices by sampling wave properties from voxel grid.

    Samples wave displacement at each plane vertex position and maps it to a color.
    Should be called every frame after wave propagation to update visualization.

    Args:
        wave_field: WaveField instance containing flux mesh fields and displacement data
        trackers: WaveTrackers instance with amplitude/frequency data for color scaling
        wave_menu: Selected Wave displayed with color palette
    """

    # ================================================================
    # XY Plane: Sample at z = fm_plane_z_idx
    # ================================================================
    # Always update all planes (conditionals cause GPU branch divergence)
    # wave_menu == 4 (Hamiltonian density) deferred to M5.0g; falls through
    # to displacement view in the meantime.
    for i, j in ti.ndrange(wave_field.nx, wave_field.ny):
        # Sample displacement at this voxel
        disp_value = wave_field.psi_am[i, j, wave_field.fm_plane_z_idx].norm()
        amp_value = trackers.amp_local_emarms_am[i, j, wave_field.fm_plane_z_idx]
        freq_value = trackers.freq_local_cross_rHz[i, j, wave_field.fm_plane_z_idx]
        univ_edge_z = wave_field.universe_size_am[2]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 4:  # ironbow
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                disp_value / univ_edge_z * warp_mesh
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
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
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                disp_value / univ_edge_z * warp_mesh
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )

    # ================================================================
    # XZ Plane: Sample at y = fm_plane_y_idx
    # ================================================================
    for i, k in ti.ndrange(wave_field.nx, wave_field.nz):
        # Sample displacement at this voxel
        disp_value = wave_field.psi_am[i, wave_field.fm_plane_y_idx, k].norm()
        amp_value = trackers.amp_local_emarms_am[i, wave_field.fm_plane_y_idx, k]
        freq_value = trackers.freq_local_cross_rHz[i, wave_field.fm_plane_y_idx, k]
        univ_edge_y = wave_field.universe_size_am[1]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 4:  # ironbow
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                disp_value / univ_edge_y * warp_mesh
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
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
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                disp_value / univ_edge_y * warp_mesh
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )

    # ================================================================
    # YZ Plane: Sample at x = fm_plane_x_idx
    # ================================================================
    for j, k in ti.ndrange(wave_field.ny, wave_field.nz):
        # Sample displacement at this voxel
        disp_value = wave_field.psi_am[wave_field.fm_plane_x_idx, j, k].norm()
        amp_value = trackers.amp_local_emarms_am[wave_field.fm_plane_x_idx, j, k]
        freq_value = trackers.freq_local_cross_rHz[wave_field.fm_plane_x_idx, j, k]
        univ_edge_x = wave_field.universe_size_am[0]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 4:  # ironbow
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                disp_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
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
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                disp_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
