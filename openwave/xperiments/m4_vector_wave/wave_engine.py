"""
ENERGY-WAVE ENGINE

ON WOLFF-LAFRENIERE METHOD

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


# ================================================================
# WAVE PROPAGATION ENGINE
# ================================================================


@ti.kernel
def select_voxels(
    wave_field: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
):
    """
    Select voxels for optimized wave computation (neighbor-only mode).

    Populates wave_field.selected_voxels with the center + 6 neighbors
    of each active wave center. Used for force gradient calculations
    when flux mesh visualization is disabled.

    Args:
        wave_field: WaveField instance with selected_voxels buffer
        wave_center: WaveCenter instance with source positions
    """
    # Reset counter
    wave_field.num_selected_voxels[None] = 0

    # 6-neighbor offsets: +x, -x, +y, -y, +z, -z (plus center = 0,0,0)
    # Using static unrolling for performance
    for wc_idx in range(wave_center.num_sources):
        if wave_center.active[wc_idx] == 0:
            continue

        # Get wave center grid position
        cx = wave_center.position_grid[wc_idx][0]
        cy = wave_center.position_grid[wc_idx][1]
        cz = wave_center.position_grid[wc_idx][2]

        # Add center voxel (offset 0)
        idx = ti.atomic_add(wave_field.num_selected_voxels[None], 1)
        if idx < wave_field.max_selected_voxels:
            wave_field.selected_voxels[idx] = ti.Vector([cx, cy, cz])

        # Add 6 neighbors with boundary clamping
        # +x neighbor
        idx = ti.atomic_add(wave_field.num_selected_voxels[None], 1)
        if idx < wave_field.max_selected_voxels:
            wave_field.selected_voxels[idx] = ti.Vector(
                [ti.min(cx + 1, wave_field.nx - 1), cy, cz]
            )

        # -x neighbor
        idx = ti.atomic_add(wave_field.num_selected_voxels[None], 1)
        if idx < wave_field.max_selected_voxels:
            wave_field.selected_voxels[idx] = ti.Vector([ti.max(cx - 1, 0), cy, cz])

        # +y neighbor
        idx = ti.atomic_add(wave_field.num_selected_voxels[None], 1)
        if idx < wave_field.max_selected_voxels:
            wave_field.selected_voxels[idx] = ti.Vector(
                [cx, ti.min(cy + 1, wave_field.ny - 1), cz]
            )

        # -y neighbor
        idx = ti.atomic_add(wave_field.num_selected_voxels[None], 1)
        if idx < wave_field.max_selected_voxels:
            wave_field.selected_voxels[idx] = ti.Vector([cx, ti.max(cy - 1, 0), cz])

        # +z neighbor
        idx = ti.atomic_add(wave_field.num_selected_voxels[None], 1)
        if idx < wave_field.max_selected_voxels:
            wave_field.selected_voxels[idx] = ti.Vector(
                [cx, cy, ti.min(cz + 1, wave_field.nz - 1)]
            )

        # -z neighbor
        idx = ti.atomic_add(wave_field.num_selected_voxels[None], 1)
        if idx < wave_field.max_selected_voxels:
            wave_field.selected_voxels[idx] = ti.Vector([cx, cy, ti.max(cz - 1, 0)])


@ti.func
def compute_voxel_wave(
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
    k_grid: ti.f32,  # type: ignore
    temporal_phase: ti.f32,  # type: ignore
):
    """
    Compute wave displacement and envelope for a single voxel.

    This function contains the core wave physics computation shared by both
    full-grid and neighbor-only iteration modes.

    Wolff-LaFreniere Combined Form:
        ψ(r,t) = A · [sin(ωt - kr) - sin(ωt)] / r

    Expanded Form (used in implementation):
        ψ(r,t) = A · [-cos(ωt)·sin(kr)/r - sin(ωt)·(1-cos(kr))/r]
               = A · [-cos(ωt)·Phase(r) - sin(ωt)·Quadrature(r)]

    Components:
        Phase:      sin(kr)/r  → k as r→0  (standing wave envelope)
        Quadrature: (1-cos(kr))/r → 0 as r→0  (traveling wave component)

    Physical Properties:
        - Near center: standing wave behavior (finite amplitude at r=0)
        - Far from center: transitions to traveling wave
        - 1/r amplitude falloff (energy conserving, from Wolff)
        - Electron core diameter = λ (full wavelength, from LaFreniere)
        - Superposition of multiple wave centers supported

    Wave Trackers Updated:
        - amp_local_rms_am: RMS amplitude via EMA on ψ² (for energy/force gradients)
        - freq_local_cross_rHz: Frequency via zero-crossing detection with EMA smoothing

    See research/01_wolff_lafreniere.md for full derivation and theory.

    Args:
        i, j, k: Voxel grid indices
        wave_field: WaveField instance
        trackers: WaveTrackers instance
        wave_center: WaveCenter instance
        dt_rs: Timestep size (rs)
        elapsed_t_rs: Elapsed simulation time (rs)
        k_grid: Angular wave number (radians per grid index)
        temporal_phase: Current temporal phase (ω·t)
    """
    # Reset before accumulation
    prev_disp = wave_field.displacement_am[i, j, k]
    wave_field.displacement_am[i, j, k] = ti.Vector([0.0, 0.0, 0.0])

    # ================================================================
    # WAVE PROPAGATION: Update voxels using wave functions
    # ================================================================
    # Loop over wave-centers (wave superposition principle causing interference)
    for wc_idx in range(wave_center.num_sources):
        # Skip inactive (annihilated) WCs
        if wave_center.active[wc_idx] == 0:
            continue

        # Compute radial distance from wave source (in grid indices)
        r_grid = ti.sqrt(
            (i - wave_center.position_grid[wc_idx][0]) ** 2
            + (j - wave_center.position_grid[wc_idx][1]) ** 2
            + (k - wave_center.position_grid[wc_idx][2]) ** 2
        )

        # Cache source-specific phase offset
        source_offset = wave_center.offset[wc_idx]

        # Spatial phase: φ = k·r, creates spherical wave fronts, dimensionless, in radians
        spatial_phase = k_grid * r_grid

        # ================================================================
        # Combined and Adjusted WOLFF-LAFRENIERE canonical form:
        #   ψ(r,t) = A · [sin(ωt - kr) - sin(ωt)] / r
        # Expanded form:
        #   ψ(r,t) = A · [-cos(ωt) · sin(kr)/r - sin(ωt) · (1 - cos(kr))/r]
        #   ψ(r,t) = A · [-cos(ωt) · Phase(r) - sin(ωt) · Quadrature(r)]
        # ================================================================
        # Phase term: sin(kr)/r → k as r→0 (physical units)
        phase_term = ti.select(
            r_grid < 0.5,  # threshold in grid units (catches center voxel only)
            k_grid,  # analytical limit
            ti.sin(spatial_phase) / r_grid,
        )

        # Quadrature term: (1-cos(kr))/r → 0 as r→0
        quadrature_term = ti.select(
            r_grid < 0.5,  # threshold in grid units (catches center voxel only)
            0.0,  # analytical limit
            (1 - ti.cos(spatial_phase)) / r_grid,
        )

        # Oscillator with source_offset phase shift
        oscillator = (
            -ti.cos(temporal_phase + source_offset) * phase_term
            - ti.sin(temporal_phase + source_offset) * quadrature_term
        )

        wave_field.displacement_am[i, j, k] += (
            base_amplitude_am * wave_field.scale_factor * oscillator
        )

    # Precision rounding to ensure wave cancellation
    # Critical for opposing phase sources (180°) that should annihilate
    # Floating-point: (+1.250001) + (-1.249999) = 0.000002 (imperfect cancel)
    # With rounding: (+1.2500) + (-1.2500) = 0.0 (perfect cancel)
    precision = ti.cast(1e4, ti.f32)  # round to 4 decimal places
    wave_field.displacement_am[i, j, k] = (
        ti.round(wave_field.displacement_am[i, j, k] * precision) / precision
    )

    curr_disp = wave_field.displacement_am[i, j, k]

    # ================================================================
    # WAVE-TRACKERS: Update amplitude and frequency trackers for visualization and forces
    # ================================================================
    # # PEAK AMPLITUDE tracking
    # ti.atomic_max(trackers.ampL_local_peak_am[i, j, k], curr_disp)
    # decay_factor_peak = ti.cast(0.999, ti.f32)  # ~100 frames to ~37%, ~230 to ~10%
    # trackers.ampL_local_peak_am[i, j, k] = (
    #     trackers.ampL_local_peak_am[i, j, k] * decay_factor_peak
    # )

    # RMS AMPLITUDE tracking via EMA on ψ² (squared displacement)
    # Running RMS: tracks √⟨ψ²⟩ - the energy-equivalent amplitude (Energy ∝ ψ²)
    # Used for: energy calculation, force gradients, visualization scaling
    # Physics: particles respond to time-averaged energy density, not
    # instantaneous displacement (inertia acts as low-pass filter at ~10²⁵ Hz)
    # EMA on ψ²: rms² = α * ψ² + (1 - α) * rms²_old, then rms = √(rms²)
    # α controls adaptation speed: higher = faster response, lower = smoother
    # RMS amplitude
    disp2 = wave_field.displacement_am[i, j, k] ** 2
    current_rms2 = trackers.amp_local_rms_am[i, j, k] ** 2
    alpha_rms = 0.005  # EMA smoothing factor for RMS tracking
    new_rms2 = alpha_rms * disp2 + (1.0 - alpha_rms) * current_rms2
    new_amp = ti.sqrt(new_rms2)

    # Unconditional decay clears trails from moving sources
    # Active regions counteract decay via EMA update from strong displacement
    # Stale regions (waves propagated away) decay to zero over time
    decay_factor = ti.cast(0.99, ti.f32)  # ~100 frames to ~37%, ~230 to ~10%
    trackers.amp_local_rms_am[i, j, k] = new_amp * decay_factor

    # TODO: review new frequency tracking method
    # FREQUENCY tracking, via zero-crossing detection with EMA smoothing
    # Detect positive-going zero crossing (negative → positive transition)
    # Period = time between consecutive positive zero crossings
    # More robust than peak detection since it's amplitude-independent
    # EMA smoothing: f_new = α * f_measured + (1 - α) * f_old
    # α controls adaptation speed: higher = faster response, lower = smoother
    if prev_disp < 0.0 and curr_disp >= 0.0:  # Zero crossing detected
        period_rs = elapsed_t_rs - trackers.last_crossing[i, j, k]
        if period_rs > dt_rs * 2:  # Filter out spurious crossings
            measured_freq = 1.0 / period_rs  # in rHz
            current_freq = trackers.freq_local_cross_rHz[i, j, k]
            alpha_freq = 0.05  # EMA smoothing factor for frequency
            trackers.freq_local_cross_rHz[i, j, k] = (
                alpha_freq * measured_freq + (1.0 - alpha_freq) * current_freq
            )
        trackers.last_crossing[i, j, k] = elapsed_t_rs

    # Unconditional frequency decay (counteracted by zero-crossing updates in active regions)
    trackers.freq_local_cross_rHz[i, j, k] *= decay_factor


@ti.kernel
def propagate_wave_neighbors(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
    num_selected: ti.i32,  # type: ignore
):
    """
    Compute wave displacement for selected neighbor voxels only (optimized mode).

    Processes only the voxels selected by select_voxels() - the center + 6 neighbors
    of each active wave center. Used for force gradient calculations when flux mesh
    visualization is disabled.

    Performance: ~7*N voxels instead of nx*ny*nz (massive gain for large grids).

    Args:
        wave_field: WaveField instance with selected_voxels populated
        trackers: WaveTrackers instance for amplitude envelope
        wave_center: WaveCenter instance with source positions
        dt_rs: Timestep size (rs)
        elapsed_t_rs: Elapsed simulation time (rs)
        num_selected: Number of selected voxels to process
    """
    # Compute wave parameters
    omega_rs = 2.0 * ti.math.pi * base_frequency_rHz / wave_field.scale_factor
    wavelength_grid = base_wavelength * wave_field.scale_factor / wave_field.dx
    k_grid = 2.0 * ti.math.pi / wavelength_grid
    temporal_phase = omega_rs * elapsed_t_rs

    # Iterate only over selected voxels (neighbors of wave centers)
    for sel_idx in range(num_selected):
        i = wave_field.selected_voxels[sel_idx][0]
        j = wave_field.selected_voxels[sel_idx][1]
        k = wave_field.selected_voxels[sel_idx][2]
        compute_voxel_wave(
            i, j, k, wave_field, trackers, wave_center, dt_rs, elapsed_t_rs, k_grid, temporal_phase
        )


@ti.kernel
def propagate_wave_full(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    wave_center: ti.template(),  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
):
    """
    Compute wave displacement using Wolff-LaFreniere analytical wave equation (full grid).

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

    # FULL GRID MODE: Compute all voxels (needed for flux mesh visualization)
    for i, j, k in ti.ndrange(nx, ny, nz):
        compute_voxel_wave(
            i, j, k, wave_field, trackers, wave_center, dt_rs, elapsed_t_rs, k_grid, temporal_phase
        )


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
_slice_xy_ampT = None
_slice_xy_freq = None
_slice_xz_amp = None
_slice_xz_ampT = None
_slice_xz_freq = None
_slice_yz_amp = None
_slice_yz_ampT = None
_slice_yz_freq = None


@ti.kernel
def _copy_slice_xy(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_ampP: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_z: ti.i32,  # type: ignore
):
    """Copy center XY slice (fixed z) to 2D buffer."""
    for i, j in slice_amp:
        slice_amp[i, j] = trackers.amp_local_rms_am[i, j, mid_z]
        slice_ampP[i, j] = trackers.amp_local_peak_am[i, j, mid_z]
        slice_freq[i, j] = trackers.freq_local_cross_rHz[i, j, mid_z]


@ti.kernel
def _copy_slice_xz(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_ampP: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_y: ti.i32,  # type: ignore
):
    """Copy center XZ slice (fixed y) to 2D buffer."""
    for i, k in slice_amp:
        slice_amp[i, k] = trackers.amp_local_rms_am[i, mid_y, k]
        slice_ampP[i, k] = trackers.amp_local_peak_am[i, mid_y, k]
        slice_freq[i, k] = trackers.freq_local_cross_rHz[i, mid_y, k]


@ti.kernel
def _copy_slice_yz(
    trackers: ti.template(),  # type: ignore
    slice_amp: ti.template(),  # type: ignore
    slice_ampP: ti.template(),  # type: ignore
    slice_freq: ti.template(),  # type: ignore
    mid_x: ti.i32,  # type: ignore
):
    """Copy center YZ slice (fixed x) to 2D buffer."""
    for j, k in slice_amp:
        slice_amp[j, k] = trackers.amp_local_rms_am[mid_x, j, k]
        slice_ampP[j, k] = trackers.amp_local_peak_am[mid_x, j, k]
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
    global _slice_xy_amp, _slice_xy_ampP, _slice_xy_freq
    global _slice_xz_amp, _slice_xz_ampP, _slice_xz_freq
    global _slice_yz_amp, _slice_yz_ampP, _slice_yz_freq

    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz

    # Initialize slice buffers once
    if _slice_xy_amp is None:
        _slice_xy_amp = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_ampP = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xy_freq = ti.field(dtype=ti.f32, shape=(nx, ny))
        _slice_xz_amp = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_ampP = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_xz_freq = ti.field(dtype=ti.f32, shape=(nx, nz))
        _slice_yz_amp = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_ampP = ti.field(dtype=ti.f32, shape=(ny, nz))
        _slice_yz_freq = ti.field(dtype=ti.f32, shape=(ny, nz))

    # Copy 3 center slices to 2D buffers (parallel kernels)
    mid_x, mid_y, mid_z = nx // 2, ny // 2, nz // 2
    _copy_slice_xy(trackers, _slice_xy_amp, _slice_xy_ampP, _slice_xy_freq, mid_z)
    _copy_slice_xz(trackers, _slice_xz_amp, _slice_xz_ampP, _slice_xz_freq, mid_y)
    _copy_slice_yz(trackers, _slice_yz_amp, _slice_yz_ampP, _slice_yz_freq, mid_x)

    # Transfer 2D slices to CPU for numpy operations
    # Exclude boundary voxels
    xy_amp = _slice_xy_amp.to_numpy()[1:-1, 1:-1]
    xy_ampP = _slice_xy_ampP.to_numpy()[1:-1, 1:-1]
    xy_freq = _slice_xy_freq.to_numpy()[1:-1, 1:-1]
    xz_amp = _slice_xz_amp.to_numpy()[1:-1, 1:-1]
    xz_ampP = _slice_xz_ampP.to_numpy()[1:-1, 1:-1]
    xz_freq = _slice_xz_freq.to_numpy()[1:-1, 1:-1]
    yz_amp = _slice_yz_amp.to_numpy()[1:-1, 1:-1]
    yz_ampP = _slice_yz_ampP.to_numpy()[1:-1, 1:-1]
    yz_freq = _slice_yz_freq.to_numpy()[1:-1, 1:-1]

    # Compute RMS amplitude: √(⟨A²⟩) for correct energy weighting
    # amp_local_rms_am contains per-voxel RMS values, square them for energy
    total_amp_squared = (xy_amp**2).sum() + (xz_amp**2).sum() + (yz_amp**2).sum()
    total_ampP_squared = (xy_ampP**2).sum() + (xz_ampP**2).sum() + (yz_ampP**2).sum()
    total_freq = xy_freq.sum() + xz_freq.sum() + yz_freq.sum()
    n_samples = xy_amp.size + xz_amp.size + yz_amp.size

    trackers.amp_global_rms_am[None] = float(np.sqrt(total_amp_squared / n_samples))
    trackers.ampP_global_rms_am[None] = float(np.sqrt(total_ampP_squared / n_samples))
    trackers.freq_global_avg_rHz[None] = float(total_freq / n_samples)


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
        psi_value = wave_field.displacement_am[i, j, wave_field.fm_plane_z_idx]
        vel_value = wave_field.velocity_am[i, j, wave_field.fm_plane_z_idx]
        amp_value = trackers.amp_local_rms_am[i, j, wave_field.fm_plane_z_idx]
        ampP_value = trackers.amp_local_peak_am[i, j, wave_field.fm_plane_z_idx]
        freq_value = trackers.freq_local_cross_rHz[i, j, wave_field.fm_plane_z_idx]
        univ_edge_z = wave_field.universe_size_am[2]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 5:  # blueprint
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + wave_field.flux_mesh_planes[2] * (
                wave_field.nz / wave_field.max_grid_size
            )
        elif wave_menu == 4:  # greenyellow
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_greenyellow_color(
                ampP_value,
                -trackers.amp_global_rms_am[None] * 2,
                trackers.amp_global_rms_am[None] * 2,
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                ampP_value / univ_edge_z * warp_mesh
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )
        elif wave_menu == 3:  # viridis
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_viridis_color(
                amp_value, 0, trackers.amp_global_rms_am[None] * 2
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                amp_value / univ_edge_z * warp_mesh
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )
        elif wave_menu == 2:  # bluered
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_bluered_color(
                vel_value,
                -trackers.amp_global_rms_am[None] * 2,
                trackers.amp_global_rms_am[None] * 2,
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                vel_value / univ_edge_z * warp_mesh
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )
        else:  # default to greenyellow (wave_menu == 1)
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_greenyellow_color(
                psi_value,
                -trackers.amp_global_rms_am[None] * 2,
                trackers.amp_global_rms_am[None] * 2,
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                psi_value / univ_edge_z * warp_mesh
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )

    # ================================================================
    # XZ Plane: Sample at y = fm_plane_y_idx
    # ================================================================
    for i, k in ti.ndrange(wave_field.nx, wave_field.nz):
        # Sample longitudinal displacement at this voxel
        psi_value = wave_field.displacement_am[i, wave_field.fm_plane_y_idx, k]
        vel_value = wave_field.velocity_am[i, wave_field.fm_plane_y_idx, k]
        amp_value = trackers.amp_local_rms_am[i, wave_field.fm_plane_y_idx, k]
        ampP_value = trackers.amp_local_peak_am[i, wave_field.fm_plane_y_idx, k]
        freq_value = trackers.freq_local_cross_rHz[i, wave_field.fm_plane_y_idx, k]
        univ_edge_y = wave_field.universe_size_am[1]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 5:  # blueprint
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + wave_field.flux_mesh_planes[1] * (
                wave_field.ny / wave_field.max_grid_size
            )
        elif wave_menu == 4:  # greenyellow
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_greenyellow_color(
                ampP_value,
                -trackers.amp_global_rms_am[None] * 2,
                trackers.amp_global_rms_am[None] * 2,
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                ampP_value / univ_edge_y * warp_mesh
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )
        elif wave_menu == 3:  # viridis
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_viridis_color(
                amp_value, 0, trackers.amp_global_rms_am[None] * 2
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                amp_value / univ_edge_y * warp_mesh
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )
        elif wave_menu == 2:  # bluered
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_bluered_color(
                vel_value,
                -trackers.amp_global_rms_am[None] * 2,
                trackers.amp_global_rms_am[None] * 2,
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                vel_value / univ_edge_y * warp_mesh
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )
        else:  # default to greenyellow (wave_menu == 1)
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_greenyellow_color(
                psi_value,
                -trackers.amp_global_rms_am[None] * 2,
                trackers.amp_global_rms_am[None] * 2,
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                psi_value / univ_edge_y * warp_mesh
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )

    # ================================================================
    # YZ Plane: Sample at x = fm_plane_x_idx
    # ================================================================
    for j, k in ti.ndrange(wave_field.ny, wave_field.nz):
        # Sample longitudinal displacement at this voxel
        psi_value = wave_field.displacement_am[wave_field.fm_plane_x_idx, j, k]
        vel_value = wave_field.velocity_am[wave_field.fm_plane_x_idx, j, k]
        amp_value = trackers.amp_local_rms_am[wave_field.fm_plane_x_idx, j, k]
        ampP_value = trackers.amp_local_peak_am[wave_field.fm_plane_x_idx, j, k]
        freq_value = trackers.freq_local_cross_rHz[wave_field.fm_plane_x_idx, j, k]
        univ_edge_x = wave_field.universe_size_am[0]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 5:  # blueprint
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_blueprint_color(
                freq_value, 0.0, trackers.freq_global_avg_rHz[None] * 2
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = freq_value / trackers.freq_global_avg_rHz[
                None
            ] / 3000 * warp_mesh + wave_field.flux_mesh_planes[0] * (
                wave_field.nx / wave_field.max_grid_size
            )
        elif wave_menu == 4:  # greenyellow
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_greenyellow_color(
                ampP_value,
                -trackers.amp_global_rms_am[None] * 2,
                trackers.amp_global_rms_am[None] * 2,
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                ampP_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
        elif wave_menu == 3:  # viridis
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_viridis_color(
                amp_value, 0, trackers.amp_global_rms_am[None] * 2
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                amp_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
        elif wave_menu == 2:  # bluered
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_bluered_color(
                vel_value,
                -trackers.amp_global_rms_am[None] * 2,
                trackers.amp_global_rms_am[None] * 2,
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                vel_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
        else:  # default to greenyellow (wave_menu == 1)
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_greenyellow_color(
                psi_value,
                -trackers.amp_global_rms_am[None] * 2,
                trackers.amp_global_rms_am[None] * 2,
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                psi_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
