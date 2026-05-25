"""
M5 ENGINE — RENDERING / VISUALIZATION (engine4_render)

The viz stack — OpenWave's biggest differential as a simulator:
  - sample_position_to_render  granule point-cloud
  - update_flux_mesh_values    3-plane scalar→color/warp mesh
  - update_director_glyphs     director-orientation line glyphs

Reads pre-computed fields off wave_field; no engine-func calls yet (M5.4 adds
engine2_pde.eigen_decompose). Full repurposing plan in 4b_rendering_features.md.
"""

import math

import taichi as ti

from openwave.common import colormap

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
        else:  # default to orange (wave_menu == 1)
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
        else:  # default to orange (wave_menu == 1)
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
        else:  # default to orange (wave_menu == 1)
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
# DIRECTOR-GLYPH RENDERING
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
# Design doc: research/2b_director_glyph_rendering.md.

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
