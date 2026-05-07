"""
ENERGY-WAVE ENGINE

ON VECTOR-WAVE METHOD

Wave Physics Engine @spacetime module. Wave dynamics and motion.
"""

import taichi as ti
import numpy as np

from openwave.common import colormap, constants

# ================================================================
# Energy-Wave Oscillation Parameters
# ================================================================
base_amplitude_am = constants.EWAVE_AMPLITUDE / constants.ATTOMETER  # am, oscillation amplitude
base_wavelength = constants.EWAVE_LENGTH  # in meters
base_wavelength_am = constants.EWAVE_LENGTH / constants.ATTOMETER  # in attometers
base_frequency = constants.EWAVE_FREQUENCY  # in Hz
base_frequency_rHz = constants.EWAVE_FREQUENCY * constants.RONTOSECOND  # in rHz (1/rontosecond)
rho = constants.MEDIUM_DENSITY  # medium density for Gaussian energy calc (kg/m³)
rho_qgam = constants.MEDIUM_DENSITY_QGAM  # qg/am³, for energy computation in scaled units


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
# WAVE PROPAGATION ENGINE — LEGACY (M4 ANALYTICAL, TO BE DELETED IN M5.0d.2)
# ================================================================


@ti.kernel
def propagate_wave(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
):
    """
    Compute wave displacement using Wolff-LaFreniere analytical wave equation.

    Physical Properties:
        - Near center: standing wave behavior (finite amplitude at r=0)
        - Far from center: transitions to traveling wave
        - 1/r amplitude falloff (energy conserving, from Wolff)
        - Electron core diameter = λ (full wavelength)
        - Superposition of multiple wave centers supported

    Wave Trackers Updated:
        - amp_local_emarms_am: RMS amplitude via EMA on ψ² (for visualization)
        - amp_local_phasorrms_am: RMS amplitude via phasor components (for energy/force gradients)
        - freq_local_cross_rHz: Frequency via zero-crossing detection with EMA smoothing

    See research/01_wolff_lafreniere.md for full derivation and theory.

    Args:
        wave_field: WaveField instance containing displacement arrays and grid info
        trackers: WaveTrackers instance for tracking wave properties
        wave_center: WaveCenter instance with source positions and phase offsets
        dt_rs: Timestep size (rs)
        elapsed_t_rs: Elapsed simulation time (rs)
    """
    # Grid dimensions for boundary handling
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    # Compute angular frequency (ω = 2πf) for temporal phase variation
    omega_rs = (
        2.0 * ti.math.pi * base_frequency_rHz / wave_field.scale_factor
    )  # angular frequency (rad/rs)

    # Compute angular wave number (k = 2π/λ) for spatial phase variation
    wavelength_grid = base_wavelength * wave_field.scale_factor / wave_field.dx
    k_grid = 2.0 * ti.math.pi / wavelength_grid  # radians per grid index

    # Temporal phase: φ = ω·t, oscillatory in time
    temporal_phase = omega_rs * elapsed_t_rs

    # ================================================================
    # WAVE PROPAGATION: Update voxels using wave functions
    # ================================================================
    # Update all voxels
    for i, j, k in ti.ndrange(nx, ny, nz):
        # Reset before accumulation
        prev_disp = wave_field.psi_am[i, j, k]
        wave_field.psi_am[i, j, k] = ti.Vector([0.0, 0.0, 0.0])

        # Loop over wave-centers (wave superposition principle causing interference)
        for wc_idx in range(wave_center.num_sources):
            # Skip inactive (annihilated) WCs
            if wave_center.active[wc_idx] == 0:
                continue

            # Get direction from voxel to wave center (normalized vector)
            dir_vec = [i, j, k] - wave_center.position_grid[wc_idx]
            r_grid = dir_vec.norm() + 1e-10
            direction = dir_vec / r_grid
            # Center voxel: direction is zero (undefined for point source),
            # use unit-z fallback so scalar oscillator reaches amplitude tracker
            direction += ti.cast(r_grid < 0.5, ti.f32) * ti.Vector([0.0, 0.0, 1.0])

            # Spatial phase: φ = k·r, creates spherical wave fronts, dimensionless, in radians
            spatial_phase = k_grid * r_grid

            # Cache source-specific phase offset
            source_offset = wave_center.offset[wc_idx]

            # ================================================================
            # WEIGHTED PARTIAL STANDING WAVE
            # ================================================================
            # Superposition of in-wave + out-wave, where the in-wave fades with distance
            # Near the source: standing wave (in-wave + out-wave interfere)
            # Far from source: pure traveling wave (in-wave fades)
            #
            # ψ(r,t) = A · [w·sin(kr+ωt+φ) + sin(kr-ωt-φ)] / kr
            #
            # Weight(r,λ): smooth decay from 1 → 0, controls standing → traveling transition
            # Standing limit (w=1, r→0): 2·cos(ωt + φ) via sinc normalization
            # Traveling limit (w=0): sin(kr - ωt - φ) / kr
            # ================================================================
            # In-wave weight: sharp Lorentzian rolloff at 1.25λ
            transition = 1 + 1 / 4  # number of wavelengths (λ)
            weight = 1.0 / (1.0 + (r_grid / (transition * wavelength_grid)) ** 8)

            # Weighted partial standing wave oscillator
            oscillator = ti.select(
                r_grid < 0.5,  # center voxel: analytical limit
                2.0 * ti.cos(temporal_phase + source_offset),  # standing wave limit: 2·cos(ωt + φ)
                (
                    weight * ti.sin(spatial_phase + (temporal_phase + source_offset))  # in-wave
                    + ti.sin(spatial_phase - (temporal_phase + source_offset))  # out-wave
                )
                / spatial_phase,  # 1/kr decay (sinc envelope)
            )

            # Accumulate this source's contribution (wave superposition)
            wave_field.psi_am[i, j, k] += (
                base_amplitude_am * wave_field.scale_factor * oscillator
            ) * direction

        # TODO: consider precision rounding to ensure perfect cancellation
        # Precision rounding to ensure wave cancellation
        # Critical for opposing phase sources (180°) that should annihilate
        # Floating-point: (+1.250001) + (-1.249999) = 0.000002 (imperfect cancel)
        # With rounding: (+1.2500) + (-1.2500) = 0.0 (perfect cancel)
        # precision = ti.cast(1e4, ti.f32)  # round to 4 decimal places
        # wave_field.psi_am[i, j, k] = (
        #     ti.round(wave_field.psi_am[i, j, k] * precision) / precision
        # )

        curr_disp = wave_field.psi_am[i, j, k]

        # ================================================================
        # WAVE-TRACKERS: Update amplitude and frequency trackers for visualization and forces
        # ================================================================
        # # PEAK AMPLITUDE tracking
        # ti.atomic_max(trackers.ampL_local_peak_am[i, j, k], curr_disp)
        # decay_factor_peak = ti.cast(0.999, ti.f32)  # ~100 frames to ~37%, ~230 to ~10%
        # trackers.ampL_local_peak_am[i, j, k] = (
        #     trackers.ampL_local_peak_am[i, j, k] * decay_factor_peak
        # )

        # TODO: review EMS method for amplitude tracking (values compared to peak amp for force)
        # TODO: review how vector displacement (has direction component) is handled in RMS calculation
        # RMS AMPLITUDE tracking via EMA on ψ² (squared displacement)
        # Running RMS: tracks √⟨ψ²⟩ - the energy-equivalent amplitude (Energy ∝ ψ²)
        # Used for: energy calculation, force gradients, visualization scaling
        # Physics: particles respond to time-averaged energy density, not
        # instantaneous displacement (inertia acts as low-pass filter at ~10²⁵ Hz)
        # EMA on ψ²: rms² = α * ψ² + (1 - α) * rms²_old, then rms = √(rms²)
        # α controls adaptation speed: higher = faster response, lower = smoother
        # RMS amplitude
        disp2 = wave_field.psi_am[i, j, k].norm() ** 2
        current_rms2 = trackers.amp_local_emarms_am[i, j, k] ** 2
        alpha_rms = 0.005  # EMA smoothing factor for RMS tracking
        new_rms2 = alpha_rms * disp2 + (1.0 - alpha_rms) * current_rms2
        new_amp = ti.sqrt(new_rms2)

        # Unconditional decay clears trails from moving sources
        # Active regions counteract decay via EMA update from strong displacement
        # Stale regions (waves propagated away) decay to zero over time
        decay_factor = ti.cast(0.99, ti.f32)  # ~100 frames to ~37%, ~230 to ~10%
        trackers.amp_local_emarms_am[i, j, k] = new_amp * decay_factor

        # TODO: review new frequency tracking method
        # FREQUENCY tracking, via zero-crossing detection with EMA smoothing
        # Detect positive-going zero crossing (negative → positive transition)
        # Period = time between consecutive positive zero crossings
        # More robust than peak detection since it's amplitude-independent
        # EMA smoothing: f_new = α * f_measured + (1 - α) * f_old
        # α controls adaptation speed: higher = faster response, lower = smoother

        # if prev_disp < 0.0 and curr_disp >= 0.0:  # Zero crossing detected
        #     period_rs = elapsed_t_rs - trackers.last_crossing[i, j, k]
        #     if period_rs > dt_rs * 2:  # Filter out spurious crossings
        #         measured_freq = 1.0 / period_rs  # in rHz
        #         current_freq = trackers.freq_local_cross_rHz[i, j, k]
        #         alpha_freq = 0.05  # EMA smoothing factor for frequency
        #         trackers.freq_local_cross_rHz[i, j, k] = (
        #             alpha_freq * measured_freq + (1.0 - alpha_freq) * current_freq
        #         )
        #     trackers.last_crossing[i, j, k] = elapsed_t_rs

        # # Unconditional frequency decay (counteracted by zero-crossing updates in active regions)
        # trackers.freq_local_cross_rHz[i, j, k] *= decay_factor

        trackers.freq_local_cross_rHz[i, j, k] = (
            base_frequency_rHz  # TODO: fixed frequency for testing, replace with above method for dynamic frequency
        )

        # ================================================================
        # LOCAL ENERGY PER VOXEL: E = ρ · V · (f · A)² in aJ (attojoules)
        # F = -∇E (force = negative energy gradient, computed in force_motion)
        # ================================================================
        # rho_qgam (qg/am³), dx_am³ (am³), f_rHz (1/rs), rms_am (am)
        # Internal units: qg·am²/rs² × 1000 → aJ
        dx_am = wave_field.dx_am
        amp_am = trackers.amp_local_emarms_am[i, j, k]
        trackers.energy_local_aJ[i, j, k] = (
            rho_qgam
            * dx_am**3
            * (base_frequency_rHz * amp_am) ** 2
            * constants.INTERNAL_ENERGY_TO_AJ
        )


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
_slice_xy_energy = None
_slice_xz_amp = None
_slice_xz_freq = None
_slice_xz_energy = None
_slice_yz_amp = None
_slice_yz_freq = None
_slice_yz_energy = None


@ti.kernel
def _copy_slice_xy(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    slice_energy: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """Copy center XY slice (fixed z) to 2D buffer."""
    for i, j in slice_amp:
        slice_amp[i, j] = trackers.amp_local_emarms_am[i, j, mid_z]
        slice_freq[i, j] = trackers.freq_local_cross_rHz[i, j, mid_z]
        slice_energy[i, j] = trackers.energy_local_aJ[i, j, mid_z]


@ti.kernel
def _copy_slice_xz(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    slice_energy: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """Copy center XZ slice (fixed y) to 2D buffer."""
    for i, k in slice_amp:
        slice_amp[i, k] = trackers.amp_local_emarms_am[i, mid_y, k]
        slice_freq[i, k] = trackers.freq_local_cross_rHz[i, mid_y, k]
        slice_energy[i, k] = trackers.energy_local_aJ[i, mid_y, k]


@ti.kernel
def _copy_slice_yz(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    slice_energy: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """Copy center YZ slice (fixed x) to 2D buffer."""
    for j, k in slice_amp:
        slice_amp[j, k] = trackers.amp_local_emarms_am[mid_x, j, k]
        slice_freq[j, k] = trackers.freq_local_cross_rHz[mid_x, j, k]
        slice_energy[j, k] = trackers.energy_local_aJ[mid_x, j, k]


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
    global _slice_xy_amp, _slice_xy_freq, _slice_xy_energy
    global _slice_xz_amp, _slice_xz_freq, _slice_xz_energy
    global _slice_yz_amp, _slice_yz_freq, _slice_yz_energy

    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    # Initialize slice buffers once
    if _slice_xy_amp is None:
        _slice_xy_amp = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_freq = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_energy = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xz_amp = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_freq = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_energy = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_yz_amp = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_freq = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_energy = ti.field(dtype=ti.f32, shape=(ny, nz))

    # Copy 3 center slices to 2D buffers (parallel kernels)
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    _copy_slice_xy(trackers, _slice_xy_amp, _slice_xy_freq, _slice_xy_energy, mid_z)
    _copy_slice_xz(trackers, _slice_xz_amp, _slice_xz_freq, _slice_xz_energy, mid_y)
    _copy_slice_yz(trackers, _slice_yz_amp, _slice_yz_freq, _slice_yz_energy, mid_x)

    # Transfer 2D slices to CPU for numpy operations
    # Exclude boundary voxels
    xy_amp = _slice_xy_amp.to_numpy()[1:-1, 1:-1]
    xy_freq = _slice_xy_freq.to_numpy()[1:-1, 1:-1]
    xy_energy = _slice_xy_energy.to_numpy()[1:-1, 1:-1]
    xz_amp = _slice_xz_amp.to_numpy()[1:-1, 1:-1]
    xz_freq = _slice_xz_freq.to_numpy()[1:-1, 1:-1]
    xz_energy = _slice_xz_energy.to_numpy()[1:-1, 1:-1]
    yz_amp = _slice_yz_amp.to_numpy()[1:-1, 1:-1]
    yz_freq = _slice_yz_freq.to_numpy()[1:-1, 1:-1]
    yz_energy = _slice_yz_energy.to_numpy()[1:-1, 1:-1]

    # Compute RMS amplitude: √(⟨A²⟩) for correct energy weighting
    # amp_local_emarms_am contains per-voxel RMS values, square them for energy
    total_amp_squared = (xy_amp**2).sum() + (xz_amp**2).sum() + (yz_amp**2).sum()
    total_freq = xy_freq.sum() + xz_freq.sum() + yz_freq.sum()
    total_energy = xy_energy.sum() + xz_energy.sum() + yz_energy.sum()
    n_samples = xy_amp.size + xz_amp.size + yz_amp.size

    trackers.amp_global_emarms_am[None] = float(np.sqrt(total_amp_squared / n_samples))
    trackers.freq_global_avg_rHz[None] = float(total_freq / n_samples)
    trackers.energy_global_avg_aJ[None] = float(total_energy / n_samples)


# ================================================================
# FLUX MESH VALUES UPDATING
# ================================================================


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
    for i, j in ti.ndrange(wave_field.nx, wave_field.ny):
        # Sample longitudinal displacement at this voxel
        disp_value = wave_field.psi_am[i, j, wave_field.fm_plane_z_idx].norm()
        amp_value = trackers.amp_local_emarms_am[i, j, wave_field.fm_plane_z_idx]
        freq_value = trackers.freq_local_cross_rHz[i, j, wave_field.fm_plane_z_idx]
        energy_value = trackers.energy_local_aJ[i, j, wave_field.fm_plane_z_idx]
        univ_edge_z = wave_field.universe_size_am[2]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 4:  # ironbow
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                energy_value, 0.0, trackers.energy_global_avg_aJ[None] * 2
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                energy_value / univ_edge_z * warp_mesh
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
        else:  # default to orange (wave_menu == 1)
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
        # Sample longitudinal displacement at this voxel
        disp_value = wave_field.psi_am[i, wave_field.fm_plane_y_idx, k].norm()
        amp_value = trackers.amp_local_emarms_am[i, wave_field.fm_plane_y_idx, k]
        freq_value = trackers.freq_local_cross_rHz[i, wave_field.fm_plane_y_idx, k]
        energy_value = trackers.energy_local_aJ[i, wave_field.fm_plane_y_idx, k]
        univ_edge_y = wave_field.universe_size_am[1]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 4:  # ironbow
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                energy_value, 0.0, trackers.energy_global_avg_aJ[None] * 2
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                energy_value / univ_edge_y * warp_mesh
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
        else:  # default to orange (wave_menu == 1)
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
        # Sample longitudinal displacement at this voxel
        disp_value = wave_field.psi_am[wave_field.fm_plane_x_idx, j, k].norm()
        amp_value = trackers.amp_local_emarms_am[wave_field.fm_plane_x_idx, j, k]
        freq_value = trackers.freq_local_cross_rHz[wave_field.fm_plane_x_idx, j, k]
        energy_value = trackers.energy_local_aJ[wave_field.fm_plane_x_idx, j, k]
        univ_edge_x = wave_field.universe_size_am[0]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 4:  # ironbow
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                energy_value, 0.0, trackers.energy_global_avg_aJ[None] * 2
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                energy_value / univ_edge_x * warp_mesh
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
        else:  # default to orange (wave_menu == 1)
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_orange_color(
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                disp_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
