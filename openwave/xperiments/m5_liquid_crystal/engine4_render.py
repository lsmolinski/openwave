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

# Single-color glyph mode (M5.6.5b "Glyph Color = single"): a flat color so weak/far-field
# glyphs stay visible (the gradient palettes fade them toward black). Sourced from the
# colormap module (COLOR_MEDIUM = light blue) — single source of truth, no hardcoded RGB.
_GLYPH_SINGLE_COLOR = colormap.COLOR_MEDIUM[1]  # (r, g, b) tuple, resolved at kernel compile time

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

    M5.4 — DIRECTOR POINT-CLOUD (interim): the granule sphere sits at voxel + amp_boost·n̂,
    where n̂ = director_nhat (principal eigenvector of M). The old "granule pushed by the
    displacement wave" has no meaning on an orientation field. This is the cheap 4b interim;
    the full biaxial-ellipsoid showcase is deferred to M5.6 (in uniaxial M5.4 the two minor
    eigen-axes are degenerate, so only a surface-of-revolution reads correctly — premature
    until δ≠g). amp_boost = WARP_MESH (dimensionless director offset; no /dx_am — n̂ is unit).
    """
    k = int(wave_field.flux_mesh_planes[2] * wave_field.grid_size[2])
    max_dim = ti.cast(wave_field.max_grid_size, ti.f32)
    sampled_ny = (wave_field.ny + stride - 1) // stride  # cols per sampled row

    for render_idx in range(num_render):
        si = render_idx // sampled_ny  # sampled row index
        sj = render_idx % sampled_ny  # sampled col index
        i = si * stride
        j = sj * stride
        displaced = amp_boost * wave_field.director_nhat[i, j, k] + ti.Vector(
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
        # Sample director (M5.4: orientation, not displacement) + matrix observables
        dir_value = wave_field.director_nhat[i, j, wave_field.fm_plane_z_idx]
        amp_value = trackers.amp_local_emarms_am[i, j, wave_field.fm_plane_z_idx]
        freq_value = trackers.freq_local_cross_rHz[i, j, wave_field.fm_plane_z_idx]
        energyH_value = observables.energyH_density_aJ[i, j, wave_field.fm_plane_z_idx]
        energyF_value = observables.energyF_density_aJ[i, j, wave_field.fm_plane_z_idx]
        univ_edge_z = wave_field.universe_size_am[2]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 1:  # Orientation deviation ‖n̂−ẑ‖ on orange (wave_menu == 1)
            dev = (dir_value - ti.Vector([0.0, 0.0, 1.0])).norm()  # 0 at vacuum ẑ, →2 at −ẑ
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_orange_color(dev, 0.0, 2.0)
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                dev / 2.0 * 0.3 * warp_mesh / 300.0
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
        elif wave_menu == 4:  # Energy density (Hamiltonian) on ironbow
            H_max = observables.energyH_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                energyH_value, 0.0, H_max
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                energyH_value / H_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size))
            )
        elif wave_menu == 5:  # Frank elastic density on ironbow (defect-focused palette)
            F_max = observables.energyF_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_ironbow_color(
                energyF_value, 0.0, F_max
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                energyF_value / F_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size))
            )
        elif (
            wave_menu == 6
        ):  # EM divergence ∇·n̂ (splay = Coulomb charge) on greenyellow diverging
            div_s = observables.director_div_absmax[None] + 1e-12
            div_v = observables.director_div_field[i, j, wave_field.fm_plane_z_idx]
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_greenyellow_color(
                div_v, -div_s, div_s
            )
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                div_v / div_s * 0.3 * warp_mesh / 300.0
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )
        elif wave_menu == 7:  # EM curl ‖∇×n̂‖ (twist+bend = circulation/B) on orange
            # Shared director-distortion scale (M5.6.5b): scale curl against max(|∇·n̂|, ‖∇×n̂‖),
            # NOT its own max — so a static hedgehog's tiny discretization-noise curl reads
            # near-black (no circulation = no B for a static charge); real twist/dynamics lights it.
            curl_s = (
                ti.max(observables.director_div_absmax[None], observables.director_curl_max[None])
                + 1e-12
            )
            curl_v = observables.director_curl_mag_field[i, j, wave_field.fm_plane_z_idx]
            wave_field.fluxmesh_xy_colors[i, j] = colormap.get_orange_color(curl_v, 0.0, curl_s)
            wave_field.fluxmesh_xy_vertices[i, j][2] = (
                curl_v / curl_s * 0.3 * warp_mesh / 300.0
                + wave_field.flux_mesh_planes[2] * (wave_field.nz / wave_field.max_grid_size)
            )

    # ================================================================
    # XZ Plane: Sample at y = fm_plane_y_idx
    # ================================================================
    for i, k in ti.ndrange(wave_field.nx, wave_field.nz):
        # Sample director (M5.4: orientation, not displacement) + matrix observables
        dir_value = wave_field.director_nhat[i, wave_field.fm_plane_y_idx, k]
        amp_value = trackers.amp_local_emarms_am[i, wave_field.fm_plane_y_idx, k]
        freq_value = trackers.freq_local_cross_rHz[i, wave_field.fm_plane_y_idx, k]
        energyH_value = observables.energyH_density_aJ[i, wave_field.fm_plane_y_idx, k]
        energyF_value = observables.energyF_density_aJ[i, wave_field.fm_plane_y_idx, k]
        univ_edge_y = wave_field.universe_size_am[1]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 1:  # Orientation deviation ‖n̂−ẑ‖ on orange (wave_menu == 1)
            dev = (dir_value - ti.Vector([0.0, 0.0, 1.0])).norm()  # 0 at vacuum ẑ, →2 at −ẑ
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_orange_color(dev, 0.0, 2.0)
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                dev / 2.0 * 0.3 * warp_mesh / 300.0
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
        elif wave_menu == 4:  # Energy density (Hamiltonian) on ironbow
            H_max = observables.energyH_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                energyH_value, 0.0, H_max
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                energyH_value / H_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size))
            )
        elif wave_menu == 5:  # Frank elastic density on ironbow
            F_max = observables.energyF_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_ironbow_color(
                energyF_value, 0.0, F_max
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                energyF_value / F_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size))
            )
        elif (
            wave_menu == 6
        ):  # EM divergence ∇·n̂ (splay = Coulomb charge) on greenyellow diverging
            div_s = observables.director_div_absmax[None] + 1e-12
            div_v = observables.director_div_field[i, wave_field.fm_plane_y_idx, k]
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_greenyellow_color(
                div_v, -div_s, div_s
            )
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                div_v / div_s * 0.3 * warp_mesh / 300.0
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )
        elif wave_menu == 7:  # EM curl ‖∇×n̂‖ (twist+bend = circulation/B) on orange
            # Shared director-distortion scale (M5.6.5b): scale curl against max(|∇·n̂|, ‖∇×n̂‖),
            # NOT its own max — so a static hedgehog's tiny discretization-noise curl reads
            # near-black (no circulation = no B for a static charge); real twist/dynamics lights it.
            curl_s = (
                ti.max(observables.director_div_absmax[None], observables.director_curl_max[None])
                + 1e-12
            )
            curl_v = observables.director_curl_mag_field[i, wave_field.fm_plane_y_idx, k]
            wave_field.fluxmesh_xz_colors[i, k] = colormap.get_orange_color(curl_v, 0.0, curl_s)
            wave_field.fluxmesh_xz_vertices[i, k][1] = (
                curl_v / curl_s * 0.3 * warp_mesh / 300.0
                + wave_field.flux_mesh_planes[1] * (wave_field.ny / wave_field.max_grid_size)
            )

    # ================================================================
    # YZ Plane: Sample at x = fm_plane_x_idx
    # ================================================================
    for j, k in ti.ndrange(wave_field.ny, wave_field.nz):
        # Sample director (M5.4: orientation, not displacement) + matrix observables
        dir_value = wave_field.director_nhat[wave_field.fm_plane_x_idx, j, k]
        amp_value = trackers.amp_local_emarms_am[wave_field.fm_plane_x_idx, j, k]
        freq_value = trackers.freq_local_cross_rHz[wave_field.fm_plane_x_idx, j, k]
        energyH_value = observables.energyH_density_aJ[wave_field.fm_plane_x_idx, j, k]
        energyF_value = observables.energyF_density_aJ[wave_field.fm_plane_x_idx, j, k]
        univ_edge_x = wave_field.universe_size_am[0]

        # Map value to color/vertex using selected gradient
        # Scale range to 2× average for headroom without saturation (allows peak visualization)
        if wave_menu == 1:  # Orientation deviation ‖n̂−ẑ‖ on orange (wave_menu == 1)
            dev = (dir_value - ti.Vector([0.0, 0.0, 1.0])).norm()  # 0 at vacuum ẑ, →2 at −ẑ
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_orange_color(dev, 0.0, 2.0)
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                dev / 2.0 * 0.3 * warp_mesh / 300.0
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
        elif wave_menu == 4:  # Energy density (Hamiltonian) on ironbow
            H_max = observables.energyH_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                energyH_value, 0.0, H_max
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                energyH_value / H_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size))
            )
        elif wave_menu == 5:  # Frank elastic density on ironbow
            F_max = observables.energyF_global_avg_aJ[None] * 4.0 + 1e-10
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_ironbow_color(
                energyF_value, 0.0, F_max
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                energyF_value / F_max * 0.3 * warp_mesh / 300.0
                + (wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size))
            )
        elif (
            wave_menu == 6
        ):  # EM divergence ∇·n̂ (splay = Coulomb charge) on greenyellow diverging
            div_s = observables.director_div_absmax[None] + 1e-12
            div_v = observables.director_div_field[wave_field.fm_plane_x_idx, j, k]
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_greenyellow_color(
                div_v, -div_s, div_s
            )
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                div_v / div_s * 0.3 * warp_mesh / 300.0
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
            )
        elif wave_menu == 7:  # EM curl ‖∇×n̂‖ (twist+bend = circulation/B) on orange
            # Shared director-distortion scale (M5.6.5b): scale curl against max(|∇·n̂|, ‖∇×n̂‖),
            # NOT its own max — so a static hedgehog's tiny discretization-noise curl reads
            # near-black (no circulation = no B for a static charge); real twist/dynamics lights it.
            curl_s = (
                ti.max(observables.director_div_absmax[None], observables.director_curl_max[None])
                + 1e-12
            )
            curl_v = observables.director_curl_mag_field[wave_field.fm_plane_x_idx, j, k]
            wave_field.fluxmesh_yz_colors[j, k] = colormap.get_orange_color(curl_v, 0.0, curl_s)
            wave_field.fluxmesh_yz_vertices[j, k][0] = (
                curl_v / curl_s * 0.3 * warp_mesh / 300.0
                + wave_field.flux_mesh_planes[0] * (wave_field.nx / wave_field.max_grid_size)
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
    observables: ti.template(),  # type: ignore
    length: ti.f32,  # type: ignore
    show_level: ti.i32,  # type: ignore
    div_scale: ti.f32,  # type: ignore
    size_mode: ti.i32,  # type: ignore
    color_mode: ti.i32,  # type: ignore
):
    """
    Update the director-glyph line-segment field by sampling ψ on the
    three flux-mesh planes (XY, XZ, YZ) at GLYPH_STRIDE.

    Per glyph, writes a CENTERED segment (VIZ.1) — the director is an apolar
    nematic axis (n̂ ≡ −n̂, arbitrary eigenvector sign, no head/tail), so it is
    drawn centered on the voxel and carries NO arrowhead barb:
      - vertices[2k+0] = base = voxel − ½·shaft·n̂
      - vertices[2k+1] = tip  = voxel + ½·shaft·n̂
      - colors[2k+0]   = colors[2k+1] = single COLOR_MEDIUM, or (color_mode=1)
        the ∇·n̂ charge gradient (greenyellow, matches WM6 mesh).
      - arrow buffers are zeroed (no barb — apolar). Centering also makes the
        glyph gauge-stable: n̂→−n̂ only swaps the endpoints, so the apolar
        sign-flip is invisible (kills the 180° "slosh" artifact).

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
        show_level: 0..3, parallel to SHOW_FLUX_MESH (0 off, 1 XY, 2 +XZ, 3 all)
        div_scale: |∇·n̂| normalization for size_mode/color_mode
        size_mode: 0=unit shaft, 1=shaft ∝ |∇·n̂| charge density
        color_mode: 0=single COLOR_MEDIUM, 1=∇·n̂ greenyellow gradient
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
            n_hat = wave_field.director_nhat[i, j, z_plane_idx]  # M5.4: principal eigenvector of M
            div_v = observables.director_div_field[i, j, z_plane_idx]

            # +0.5 on the two in-plane axes aligns the glyph origin with the
            # flux-mesh cell center (create_flux_mesh: verts at (idx+0.5)/max_dim).
            # Perpendicular uses the mesh's own continuous plane coord (fm_plane_z_pos),
            # not z_plane_idx/max_dim — keeps the glyph welded to the sheet on any grid.
            pos = ti.Vector(
                [
                    (ti.cast(i, ti.f32) + 0.5) / max_dim,
                    (ti.cast(j, ti.f32) + 0.5) / max_dim,
                    wave_field.fm_plane_z_pos,
                ]
            )
            shaft = length
            if size_mode == 1:  # magnitude = |∇·n̂| charge density
                shaft = length * ti.min(ti.abs(div_v) / (div_scale + 1e-12), 1.0)
            # VIZ.1 (M5.6.5b): the director glyph is ALWAYS centered on the voxel
            # (base = pos − ½·shaft·n̂). The director is an apolar nematic axis
            # (n̂ ≡ −n̂, the eigenvector sign is arbitrary), so it has no head/tail
            # and its rendered center belongs at the voxel it samples. Centering also
            # makes it gauge-stable: n̂→−n̂ just swaps the two endpoints → the segment
            # is visually identical → no 180° "slosh" from the apolar sign-flip. Any
            # remaining motion is real (tilt + free-defect orientation dispersal).
            base = pos - 0.5 * shaft * n_hat
            tip = base + shaft * n_hat
            # Color = palette mapping of (1 − n_hat[2]) ∈ [0, 2]:
            # vacuum (n=ẑ) → 0 (DARK, blends into black GUI background, defect
            # is what stands out); equator (n in xy-plane) → 1 (mid); inward
            # south-pole (n=-ẑ) → 2 (BRIGHT, peak twist away from vacuum).
            # CHANGE THIS LINE to test palettes (both start dark, end bright):
            #   colormap.get_ironbow_color(...)    black → magenta → red → yellow-white
            #   colormap.get_blueprint_color(...)  dark blue → light blue (current)
            # M5.6.5b: color by ∇·n̂ (charge), greenyellow signed — matches WM6 mesh
            color = ti.Vector([_GLYPH_SINGLE_COLOR[0], _GLYPH_SINGLE_COLOR[1], _GLYPH_SINGLE_COLOR[2]])  # COLOR_MEDIUM single (see far field)
            if color_mode == 1:
                color = colormap.get_greenyellow_color(div_v, -div_scale, div_scale)
            wave_field.director_glyph_vertices[idx + 0] = base
            wave_field.director_glyph_vertices[idx + 1] = tip
            wave_field.director_glyph_colors[idx + 0] = color
            wave_field.director_glyph_colors[idx + 1] = color

            # No barb: the director is apolar (n̂ ≡ −n̂, no head/tail). Zero the
            # arrow buffers so the shared arrow-render pass draws nothing for it.
            wave_field.director_glyph_arrow_vertices[idx + 0] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 1] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 0] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 1] = zero_v
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
            n_hat = wave_field.director_nhat[i, y_plane_idx, k]  # M5.4: principal eigenvector of M
            div_v = observables.director_div_field[i, y_plane_idx, k]

            # +0.5 in-plane (i, k); perpendicular = mesh's continuous fm_plane_y_pos.
            pos = ti.Vector(
                [
                    (ti.cast(i, ti.f32) + 0.5) / max_dim,
                    wave_field.fm_plane_y_pos,
                    (ti.cast(k, ti.f32) + 0.5) / max_dim,
                ]
            )
            shaft = length
            if size_mode == 1:  # magnitude = |∇·n̂| charge density
                shaft = length * ti.min(ti.abs(div_v) / (div_scale + 1e-12), 1.0)
            # VIZ.1 (M5.6.5b): the director glyph is ALWAYS centered on the voxel
            # (base = pos − ½·shaft·n̂). The director is an apolar nematic axis
            # (n̂ ≡ −n̂, the eigenvector sign is arbitrary), so it has no head/tail
            # and its rendered center belongs at the voxel it samples. Centering also
            # makes it gauge-stable: n̂→−n̂ just swaps the two endpoints → the segment
            # is visually identical → no 180° "slosh" from the apolar sign-flip. Any
            # remaining motion is real (tilt + free-defect orientation dispersal).
            base = pos - 0.5 * shaft * n_hat
            tip = base + shaft * n_hat
            # Color = palette mapping of (1 − n_hat[2]) ∈ [0, 2]:
            # vacuum (n=ẑ) → 0 (DARK, blends into black GUI background, defect
            # is what stands out); equator (n in xy-plane) → 1 (mid); inward
            # south-pole (n=-ẑ) → 2 (BRIGHT, peak twist away from vacuum).
            # CHANGE THIS LINE to test palettes (both start dark, end bright):
            #   colormap.get_ironbow_color(...)    black → magenta → red → yellow-white
            #   colormap.get_blueprint_color(...)  dark blue → light blue (current)
            color = ti.Vector([_GLYPH_SINGLE_COLOR[0], _GLYPH_SINGLE_COLOR[1], _GLYPH_SINGLE_COLOR[2]])  # COLOR_MEDIUM single (see far field)
            if color_mode == 1:
                color = colormap.get_greenyellow_color(div_v, -div_scale, div_scale)
            wave_field.director_glyph_vertices[idx + 0] = base
            wave_field.director_glyph_vertices[idx + 1] = tip
            wave_field.director_glyph_colors[idx + 0] = color
            wave_field.director_glyph_colors[idx + 1] = color

            # No barb: the director is apolar (n̂ ≡ −n̂, no head/tail). Zero the
            # arrow buffers so the shared arrow-render pass draws nothing for it.
            wave_field.director_glyph_arrow_vertices[idx + 0] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 1] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 0] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 1] = zero_v
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
            n_hat = wave_field.director_nhat[x_plane_idx, j, k]  # M5.4: principal eigenvector of M
            div_v = observables.director_div_field[x_plane_idx, j, k]

            # +0.5 in-plane (j, k); perpendicular = mesh's continuous fm_plane_x_pos.
            pos = ti.Vector(
                [
                    wave_field.fm_plane_x_pos,
                    (ti.cast(j, ti.f32) + 0.5) / max_dim,
                    (ti.cast(k, ti.f32) + 0.5) / max_dim,
                ]
            )
            shaft = length
            if size_mode == 1:  # magnitude = |∇·n̂| charge density
                shaft = length * ti.min(ti.abs(div_v) / (div_scale + 1e-12), 1.0)
            # VIZ.1 (M5.6.5b): the director glyph is ALWAYS centered on the voxel
            # (base = pos − ½·shaft·n̂). The director is an apolar nematic axis
            # (n̂ ≡ −n̂, the eigenvector sign is arbitrary), so it has no head/tail
            # and its rendered center belongs at the voxel it samples. Centering also
            # makes it gauge-stable: n̂→−n̂ just swaps the two endpoints → the segment
            # is visually identical → no 180° "slosh" from the apolar sign-flip. Any
            # remaining motion is real (tilt + free-defect orientation dispersal).
            base = pos - 0.5 * shaft * n_hat
            tip = base + shaft * n_hat
            # Color = palette mapping of (1 − n_hat[2]) ∈ [0, 2]:
            # vacuum (n=ẑ) → 0 (DARK, blends into black GUI background, defect
            # is what stands out); equator (n in xy-plane) → 1 (mid); inward
            # south-pole (n=-ẑ) → 2 (BRIGHT, peak twist away from vacuum).
            # CHANGE THIS LINE to test palettes (both start dark, end bright):
            #   colormap.get_ironbow_color(...)    black → magenta → red → yellow-white
            #   colormap.get_blueprint_color(...)  dark blue → light blue (current)
            color = ti.Vector([_GLYPH_SINGLE_COLOR[0], _GLYPH_SINGLE_COLOR[1], _GLYPH_SINGLE_COLOR[2]])  # COLOR_MEDIUM single (see far field)
            if color_mode == 1:
                color = colormap.get_greenyellow_color(div_v, -div_scale, div_scale)
            wave_field.director_glyph_vertices[idx + 0] = base
            wave_field.director_glyph_vertices[idx + 1] = tip
            wave_field.director_glyph_colors[idx + 0] = color
            wave_field.director_glyph_colors[idx + 1] = color

            # No barb: the director is apolar (n̂ ≡ −n̂, no head/tail). Zero the
            # arrow buffers so the shared arrow-render pass draws nothing for it.
            wave_field.director_glyph_arrow_vertices[idx + 0] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 1] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 0] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 1] = zero_v
        else:
            wave_field.director_glyph_vertices[idx + 0] = zero_v
            wave_field.director_glyph_vertices[idx + 1] = zero_v
            wave_field.director_glyph_colors[idx + 0] = zero_v
            wave_field.director_glyph_colors[idx + 1] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 0] = zero_v
            wave_field.director_glyph_arrow_vertices[idx + 1] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 0] = zero_v
            wave_field.director_glyph_arrow_colors[idx + 1] = zero_v

@ti.kernel
def update_em_vector_glyphs(
    wave_field: ti.template(),  # type: ignore
    observables: ti.template(),  # type: ignore
    length: ti.f32,  # type: ignore
    scale: ti.f32,  # type: ignore
    show_level: ti.i32,  # type: ignore
    size_mode: ti.i32,  # type: ignore
    color_mode: ti.i32,  # type: ignore
):
    """M5.6.5b — B-direction glyphs: half-arrow segments along ∇×n̂ (the curl/circulation vector).

    Reuses the director-glyph BUFFERS (shaft + arrowhead) + 3-plane sampling, but the segment
    points along the curl vector observables.director_curl_field, colored by ‖∇×n̂‖ on the
    BLUERED gradient. Shaft length ∝ min(‖curl‖/scale, 1) so the view declutters where there is
    no circulation (static charge ⇒ ~zero-length, invisible; real twist ⇒ visible arrows with a
    half-barb tip showing the B direction). The barb scales with the shaft (also vanishes at
    curl≈0). `scale` = the shared distortion magnitude max(|∇·n̂|, ‖∇×n̂‖) (matches WAVE_MENU 7).
    show_level mirrors SHOW_DIRECTORS (0 off, 1 XY, 2 +XZ, 3 all). Writes director_glyph_vertices/
    colors + director_glyph_arrow_vertices/colors; the launcher renders both via scene.lines.
    """
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    stride = wave_field.GLYPH_STRIDE
    max_dim = ti.cast(wave_field.max_grid_size, ti.f32)
    nx_s, ny_s, nz_s = wave_field.glyph_nx_s, wave_field.glyph_ny_s, wave_field.glyph_nz_s
    z_idx, y_idx, x_idx = wave_field.fm_plane_z_idx, wave_field.fm_plane_y_idx, wave_field.fm_plane_x_idx
    inv_s = 1.0 / (scale + 1e-12)
    zero_v = ti.Vector([0.0, 0.0, 0.0])

    for plane in range(3):
        n_si = nx_s if plane < 2 else ny_s          # plane 0=XY, 1=XZ, 2=YZ
        n_sj = ny_s if plane == 0 else nz_s
        base = wave_field.glyph_offset_xy
        if plane == 1:
            base = wave_field.glyph_offset_xz
        if plane == 2:
            base = wave_field.glyph_offset_yz
        for sa, sb in ti.ndrange(n_si, n_sj):
            idx = (base + sa * n_sj + sb) * 2
            on = (plane == 0 and show_level >= 1) or (plane == 1 and show_level >= 2) \
                or (plane == 2 and show_level >= 3)
            if on:
                # map (sa, sb) → voxel (i, j, k) on the active plane
                i = ti.min(sa * stride, nx - 1) if plane < 2 else x_idx
                j = y_idx
                k = z_idx
                if plane == 0:
                    j = ti.min(sb * stride, ny - 1)
                if plane == 1:
                    k = ti.min(sb * stride, nz - 1)
                if plane == 2:
                    j = ti.min(sa * stride, ny - 1)
                    k = ti.min(sb * stride, nz - 1)
                curl = observables.director_curl_field[i, j, k]
                mag = curl.norm()
                dirv = curl / (mag + 1e-12)
                # +0.5 on the two in-plane axes (flux-mesh cell center); the
                # perpendicular axis is then snapped to the mesh's continuous plane
                # coord (fm_plane_*_pos) so glyphs sit exactly on the sheet.
                half = ti.Vector([0.5, 0.5, 0.5])
                if plane == 0:      # XY slice → k is perpendicular
                    half[2] = 0.0
                elif plane == 1:    # XZ slice → j is perpendicular
                    half[1] = 0.0
                else:               # YZ slice → i is perpendicular
                    half[0] = 0.0
                pos = (ti.Vector([ti.cast(i, ti.f32), ti.cast(j, ti.f32),
                                  ti.cast(k, ti.f32)]) + half) / max_dim
                if plane == 0:
                    pos[2] = wave_field.fm_plane_z_pos
                elif plane == 1:
                    pos[1] = wave_field.fm_plane_y_pos
                else:
                    pos[0] = wave_field.fm_plane_x_pos
                # shaft: unit (dirv≈unit wherever curl detectable) or field magnitude (declutters)
                shaft = length
                if size_mode == 1:
                    shaft = length * ti.min(mag * inv_s, 1.0)
                tip = pos + shaft * dirv
                color = ti.Vector([_GLYPH_SINGLE_COLOR[0], _GLYPH_SINGLE_COLOR[1], _GLYPH_SINGLE_COLOR[2]])  # COLOR_MEDIUM single (see far field)
                if color_mode == 1:
                    color = colormap.get_orange_color(mag, 0.0, scale + 1e-12)  # magnitude, matches WM7
                wave_field.director_glyph_vertices[idx + 0] = pos
                wave_field.director_glyph_vertices[idx + 1] = tip
                wave_field.director_glyph_colors[idx + 0] = color
                wave_field.director_glyph_colors[idx + 1] = color
                # half-arrow barb at the tip (stable perp; /(norm+eps) avoids NaN at curl≈0)
                ref = ti.Vector([0.0, 0.0, 1.0])
                if ti.abs(dirv[2]) > 0.9:
                    ref = ti.Vector([1.0, 0.0, 0.0])
                perp = dirv.cross(ref)
                perp = perp / (perp.norm() + 1e-12)
                barb_dir = ARROWHEAD_BACK_COMP * dirv + ARROWHEAD_PERP_COMP * perp
                barb_end = tip + (ARROWHEAD_LENGTH_FRAC * shaft) * barb_dir
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
