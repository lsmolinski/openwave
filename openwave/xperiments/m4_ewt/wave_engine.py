"""
ENERGY-WAVE ENGINE

ON VECTOR-WAVE MODEL

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
# WAVE PROPAGATION ENGINE  (vector PDE solver)
# ================================================================
# Time-integrated leapfrog vector wave equation:  ∂²ψ/∂t² = c²∇²ψ
# (the nonlinear restoring term −dV(ψ) is added in P2).
#
# A single 3-component field ψ (psi_am) is evolved over ALL voxels — a PDE has
# no "selective voxel" shortcut. The structure mirrors M2's free-wave engine,
# vectorized componentwise; seeding takes the multi-center superposition insight
# from M5. NO vector calculus (div/curl) and NO director field (see #203 scope).
# ================================================================


@ti.kernel
def seed_wave(
    wave_field: ti.template(),  # type: ignore
    seed_mode: ti.i32,  # type: ignore
    boost: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
):
    """
    Seed the medium's BASE WAVE: the always-on background oscillation of the medium
    at its ground state (EWT). This is NOT sourced from the wave centers — in EWT the
    base wave fills the whole medium and the wave centers are localized disturbances on
    top of it, re-driven separately by the WC interaction modes (P3). The base wave
    radiates from the DOMAIN (universe) center voxel as a radial (longitudinal)
    displacement ψ = A·profile(r)·r̂. Interior only (Dirichlet ψ=0 at the outer boundary).

    seed_mode:
        0 = gaussian : A·exp(-r²/2σ²)·r̂ (σ = λ/2) localized pulse, released from rest
        1 = radial   : A·exp(-r²/2σ²)·cos(k·r)·r̂ cosine-modulated pulse, released from rest
        2 = full     : A·cos(ω·t − k·r)·r̂ domain-filling outgoing radial wave (à la M2
                       charge_full); ψ_prev is set at t = −dt to give an initial outgoing
                       velocity. Closest emulation of the always-on base wave.

    Args:
        wave_field: WaveField instance (psi_am / psi_prev_am / psi_new_am)
        seed_mode: base-wave profile selector
        boost: amplitude multiplier
        dt_rs: timestep (rs), used by mode 2 to set the t = −dt previous buffer
    """
    # Domain (universe) center, in grid indices — the base-wave source
    cx = wave_field.nx * 0.5
    cy = wave_field.ny * 0.5
    cz = wave_field.nz * 0.5

    # Wave numbers in grid/scaled units; Gaussian width σ = λ/2 (grid units)
    wavelength_grid = base_wavelength * wave_field.scale_factor / wave_field.dx
    k_grid = 2.0 * ti.math.pi / wavelength_grid  # radians per grid index
    omega_rs = 2.0 * ti.math.pi * base_frequency_rHz / wave_field.scale_factor  # rad/rs
    sigma_grid = wavelength_grid / 2.0
    amp = base_amplitude_am * wave_field.scale_factor * boost

    # Seed ψ (and ψ_prev) over interior voxels (Dirichlet ψ=0 at edges)
    for i, j, k in ti.ndrange(
        (1, wave_field.nx - 1), (1, wave_field.ny - 1), (1, wave_field.nz - 1)
    ):
        # Radial direction from the DOMAIN CENTER to this voxel
        dir_vec = ti.Vector(
            [ti.cast(i, ti.f32) - cx, ti.cast(j, ti.f32) - cy, ti.cast(k, ti.f32) - cz]
        )
        r_grid = dir_vec.norm()
        r_hat = dir_vec / (r_grid + 1e-6)
        node = ti.select(r_grid < 0.5, 0.0, 1.0)  # node at the exact center (r̂ undefined)
        envelope = ti.exp(-(r_grid * r_grid) / (2.0 * sigma_grid * sigma_grid))

        # Default: gaussian pulse released from rest (mode 0)
        psi = amp * envelope * r_hat * node
        psi_prev = psi
        if seed_mode == 1:  # cosine-modulated radial pulse, released from rest
            psi = amp * envelope * ti.cos(k_grid * r_grid) * r_hat * node
            psi_prev = psi
        elif seed_mode == 2:  # domain-filling outgoing base wave (à la charge_full)
            psi = amp * ti.cos(-k_grid * r_grid) * r_hat * node
            psi_prev = amp * ti.cos(omega_rs * (-dt_rs) - k_grid * r_grid) * r_hat * node

        wave_field.psi_am[i, j, k] = psi
        wave_field.psi_prev_am[i, j, k] = psi_prev

    # Clear the next buffer
    for i, j, k in ti.ndrange(wave_field.nx, wave_field.ny, wave_field.nz):
        wave_field.psi_new_am[i, j, k] = ti.Vector([0.0, 0.0, 0.0])


@ti.func
def compute_laplacian(
    wave_field: ti.template(),  # type: ignore
    i: ti.i32,  # type: ignore
    j: ti.i32,  # type: ignore
    k: ti.i32,  # type: ignore
):
    """
    Vector Laplacian ∇²ψ at voxel [i,j,k], 6-point face-neighbor stencil applied
    componentwise (a plain vector Laplacian, NO div/curl decomposition).
        ∇²ψ ≈ (face_sum - 6·center) / dx²

    Args:
        wave_field: WaveField instance (psi_am)
        i, j, k: interior voxel indices (0 < i,j,k < n-1)

    Returns:
        3-vector Laplacian in units [am / am²] = [1/am]
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
    return (face_sum - 6.0 * center) / (wave_field.dx_am**2)


@ti.kernel
def propagate_wave(
    wave_field: ti.template(),  # type: ignore
    trackers: ti.template(),  # type: ignore
    c_amrs: ti.f32,  # type: ignore
    dt_rs: ti.f32,  # type: ignore
    elapsed_t_rs: ti.f32,  # type: ignore
):
    """
    Evolve ψ one step by leapfrog over all interior voxels (vector wave PDE).

    Discrete form (Verlet / leapfrog), linear for P1:
        ψ_new = 2ψ - ψ_prev + (c·dt)²·∇²ψ
    The nonlinear restoring term −dt²·dV(ψ) is added in P2.

    Boundary: Dirichlet ψ=0 at edges (boundary voxels skipped).
    Trackers (RMS amplitude, zero-crossing frequency) and per-voxel energy are
    updated in-kernel using M4's existing fields/formulas, then buffers are swapped.

    Args:
        wave_field: WaveField (psi_am / psi_prev_am / psi_new_am)
        trackers: WaveTrackers (amp / freq / energy fields)
        c_amrs: wave speed (am/rs)
        dt_rs: timestep (rs)
        elapsed_t_rs: elapsed simulation time (rs)
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    c2dt2 = (c_amrs * dt_rs) ** 2

    # Domain center: reference for the radial (longitudinal) zero-crossing scalar
    cx = ti.cast(nx, ti.f32) * 0.5
    cy = ti.cast(ny, ti.f32) * 0.5
    cz = ti.cast(nz, ti.f32) * 0.5

    # Update interior voxels (Dirichlet BC: ψ=0 at edges; 6-pt stencil needs a 1-cell buffer)
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        laplacian = compute_laplacian(wave_field, i, j, k)

        # Leap-Frog update (linear; nonlinear V added in P2)
        # Standard form: ψ_new = 2ψ - ψ_prev + (c·dt)²·∇²ψ
        psi_new = (
            2.0 * wave_field.psi_am[i, j, k] - wave_field.psi_prev_am[i, j, k] + c2dt2 * laplacian
        )
        wave_field.psi_new_am[i, j, k] = psi_new

        # WAVE-TRACKERS ============================================
        # RMS AMPLITUDE via EMA on |ψ|², with M4's unconditional decay (clears stale trails)
        disp2 = psi_new.norm() ** 2
        current_rms2 = trackers.amp_local_emarms_am[i, j, k] ** 2
        alpha_rms = 0.005  # EMA smoothing factor for RMS tracking
        new_amp = ti.sqrt(alpha_rms * disp2 + (1.0 - alpha_rms) * current_rms2)
        decay_factor = ti.cast(0.99, ti.f32)  # ~100 frames to ~37%
        trackers.amp_local_emarms_am[i, j, k] = new_amp * decay_factor

        # FREQUENCY via zero-crossing of the radial (longitudinal) projection s = ψ·r̂
        # Reference direction = from domain center (a signed scalar that oscillates everywhere)
        rdir = ti.Vector(
            [ti.cast(i, ti.f32) - cx, ti.cast(j, ti.f32) - cy, ti.cast(k, ti.f32) - cz]
        )
        r_hat = rdir / (rdir.norm() + 1e-6)
        s_prev = wave_field.psi_am[i, j, k].dot(r_hat)
        s_curr = psi_new.dot(r_hat)
        if s_prev < 0.0 and s_curr >= 0.0:  # positive-going zero crossing
            period_rs = elapsed_t_rs - trackers.last_crossing[i, j, k]
            if period_rs > dt_rs * 2.0:  # reject spurious crossings
                measured_freq = 1.0 / period_rs  # in rHz
                current_freq = trackers.freq_local_cross_rHz[i, j, k]
                alpha_freq = 0.05  # EMA smoothing factor for frequency
                trackers.freq_local_cross_rHz[i, j, k] = (
                    alpha_freq * measured_freq + (1.0 - alpha_freq) * current_freq
                )
            trackers.last_crossing[i, j, k] = elapsed_t_rs

        # LOCAL ENERGY PER VOXEL: E = ρ · V · (f · A)² in aJ (M4 formula kept; f = base frequency)
        # F = -∇E (force = negative energy gradient, computed in force_motion)
        amp_am = trackers.amp_local_emarms_am[i, j, k]
        trackers.energy_local_aJ[i, j, k] = (
            rho_qgam
            * wave_field.dx_am**3
            * (base_frequency_rHz * amp_am) ** 2
            * constants.INTERNAL_ENERGY_TO_AJ
        )

    # Swap time levels: prev ← current, current ← new
    for i, j, k in ti.ndrange(nx, ny, nz):
        wave_field.psi_prev_am[i, j, k] = wave_field.psi_am[i, j, k]
        wave_field.psi_am[i, j, k] = wave_field.psi_new_am[i, j, k]


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
        if wave_menu == 1:  # orange
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_orange_color(
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                disp_value / univ_edge_z * warp_mesh
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )
        elif wave_menu == 2:  # ironbow
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                amp_value / univ_edge_z * warp_mesh
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
        elif wave_menu == 4:  # ironbow
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                energy_value, 0.0, trackers.energy_global_avg_aJ[None] * 2
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                energy_value / univ_edge_z * warp_mesh / 10
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
        if wave_menu == 1:  # orange
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_orange_color(
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                disp_value / univ_edge_y * warp_mesh
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )
        elif wave_menu == 2:  # ironbow
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                amp_value / univ_edge_y * warp_mesh
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
        elif wave_menu == 4:  # ironbow
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                energy_value, 0.0, trackers.energy_global_avg_aJ[None] * 2
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                energy_value / univ_edge_y * warp_mesh / 10
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
        if wave_menu == 1:  # orange
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_orange_color(
                disp_value, 0, trackers.amp_global_emarms_am[None] * 4
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                disp_value / univ_edge_x * warp_mesh
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
        elif wave_menu == 2:  # ironbow
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                amp_value, 0, trackers.amp_global_emarms_am[None] * 2
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                amp_value / univ_edge_x * warp_mesh
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
        elif wave_menu == 4:  # ironbow
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                energy_value, 0.0, trackers.energy_global_avg_aJ[None] * 2
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                energy_value / univ_edge_x * warp_mesh / 10
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
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
