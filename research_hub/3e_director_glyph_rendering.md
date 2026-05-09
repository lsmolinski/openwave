# Director-Glyph Rendering — Design Note

**Status**: design-stage (2026-05-08), to be implemented as part of M5.1 alongside `seed_hedgehog`.

**Purpose**: visualize the **vector orientation** of ψ at each voxel — the structure that carries topology — distinct from `flux_mesh` which renders only scalar projections (magnitude, energy density). Without this, you can compute a winding number but you can't *see* a hedgehog. Critical for M5.1 validation (does the relaxed seed actually look like a hedgehog?) and persists through M5.4 (multi-defect dynamics) and M5.8 (Zitterbewegung at the defect core).

## Background — why scalar visualization is insufficient for topology

OpenWave's existing visualization (`update_flux_mesh_values` in `lagrangian_engine.py`) maps each voxel to ONE scalar (color + Z-warp height). For magnitude-style observables (`|ψ|`, `H_density`) this is the right choice. But topology lives in the **directional structure** of the field — a hedgehog is "directors radiating outward from a central point", which is a property of `ψ̂ = ψ / |ψ|`, not of `|ψ|`.

A magnitude rendering of a hedgehog looks visually identical to vacuum (uniform `|ψ|` everywhere). The only way to *see* topology is to render the vector orientation directly. That's what director-glyph rendering does.

This is also useful in the M5.0 wave regime — the same renderer shows vector polarization at each voxel for traveling waves. For an x-traveling y-polarized wave (the smoke test), every glyph oscillates between ±ŷ as the wave passes — useful even before topology lands.

## Approach: line-segment glyphs on the 3 flux-mesh planes

Three escalating options; M5.1 starts with A.

### A. 3-plane line glyphs (M5.1 first cut)

Sample ψ at every Nth voxel on the three flux-mesh planes (`fm_plane_x_idx`, `fm_plane_y_idx`, `fm_plane_z_idx` — same locations the existing flux mesh uses). For each sample:

- Compute `n̂ = ψ / |ψ|` (with epsilon guard against zero magnitude in vacuum)
- Render a line segment from voxel position to `voxel + L · n̂` (small offset, e.g. `L = 0.02` in normalized [0,1] camera space)
- Color it by direction (signed-component RGB — see "Color encoding" below)

Per Taichi GGUI: `render.scene.lines(vertex_field, per_vertex_color=color_field, width=2.0)`. Same primitive used by `axis_field` and `grid_lines` in `render.py`.

**Glyph density**: at 256³ with stride 4, you get `(64·64) × 3 ≈ 12k glyphs` total — fast on Metal. At 384³ stride 4 → ~28k glyphs, still cheap.

**User control**: a `SHOW_DIRECTORS` checkbox in `display_controls()` parallel to `SHOW_FLUX_MESH`.

### B. 3D glyph cloud everywhere

Same as A but iterate over the whole grid every Nth voxel. Drowns in noise quickly; only useful at low resolution or with aggressive stride. Probably skip unless we hit a specific need.

### C. Glyphs only near defects (M5.4+)

Once `winding_number(center, r)` is implemented (M5.1 task #6), use it to locate defects, then render glyphs only in a sphere around each defect core. Adaptive density (denser near core, sparser far). Best for M5.4 multi-defect scenarios where 3-plane sampling might miss off-axis defects. Add later as needed.

## Color encoding direction

You want **opposite directions to look visually opposite** so a +1 hedgehog (outward) is distinguishable from a −1 hedgehog (inward) at a glance. Two mapping choices, only one is correct:

| Mapping | Effect |
| --- | --- |
| `color = abs(n̂)` (absolute) | +1 and −1 hedgehogs look identical (both "rainbow starburst") — **WRONG** |
| `color = (n̂ + 1) / 2` (signed-to-unit) | +x is red, −x is cyan, +y green, −y magenta, +z blue, −z yellow — opposites are RGB-complementary — **CORRECT** |

Use the signed mapping. A +1 hedgehog will show red on the right side (n̂ points +x) and cyan on the left (n̂ points −x); flipping the sign inverts the color pattern. You can tell defect polarity without needing arrowheads.

## Why arrowheads are skippable for first pass

Pure line segments are bidirectional — a line from `(0,0,0)` to `(1,0,0)` looks identical to one from `(1,0,0)` to `(0,0,0)`. But because the kernel draws *from voxel center toward `voxel + L · n̂`*, direction is encoded in the offset:

- Right side of a +1 hedgehog: `n̂ = +x̂` → glyphs reach further right
- Left side of a +1 hedgehog: `n̂ = −x̂` → glyphs reach further left

Combined with signed-component RGB color, direction is unambiguous even without arrowheads.

If we want explicit arrowheads later, each glyph becomes 3 line segments (shaft + 2 head lines) = 6 vertices instead of 2. Tripled vertex count, slightly nicer aesthetics. Add only if first-pass rendering proves visually ambiguous for some configuration.

## Implementation sketch (~80 lines total)

```python
# medium.py — in WaveField.__init__:
GLYPH_STRIDE = 4
n_glyphs = (
    (self.nx // GLYPH_STRIDE) * (self.ny // GLYPH_STRIDE)        # XY plane
    + (self.nx // GLYPH_STRIDE) * (self.nz // GLYPH_STRIDE)      # XZ plane
    + (self.ny // GLYPH_STRIDE) * (self.nz // GLYPH_STRIDE)      # YZ plane
)
# 2 vertices per glyph for the line segment
self.director_glyph_vertices = ti.Vector.field(3, ti.f32, shape=(2 * n_glyphs,))
self.director_glyph_colors  = ti.Vector.field(3, ti.f32, shape=(2 * n_glyphs,))
```

```python
# lagrangian_engine.py — new kernel:
@ti.kernel
def update_director_glyphs(
    wave_field: ti.template(),
    stride: ti.i32,
    length: ti.f32,
):
    """3-plane sample of n̂ = ψ / |ψ|, written as line-segment glyphs
    with signed-component RGB colors. Mirror of update_flux_mesh_values
    but for vector orientation instead of scalar magnitude.

    Per-glyph encoding:
        vertices[2k+0] = voxel position (normalized [0,1] camera space)
        vertices[2k+1] = voxel position + length · n̂
        colors[2k+0]   = colors[2k+1] = (n̂ + 1) / 2   — signed-component RGB

    Boundary handling: 1e-10 epsilon in n̂ denominator avoids /0 in vacuum
    voxels where ψ ≈ 0; degenerate glyphs collapse to a single point and
    don't render visibly.
    """
    # Loop body: 3 nested ti.ndrange blocks, one per plane (XY, XZ, YZ)
    # Each writes its slice of glyph_vertices/colors at flat offset
    # offset_xy = 0
    # offset_xz = (nx//stride) * (ny//stride)
    # offset_yz = offset_xz + (nx//stride) * (nz//stride)
    ...
```

```python
# _launcher.py — in main render loop, alongside flux_mesh rendering:
if state.SHOW_DIRECTORS:
    lagrange.update_director_glyphs(state.wave_field, stride=4, length=0.02)
    render.scene.lines(
        state.wave_field.director_glyph_vertices,
        per_vertex_color=state.wave_field.director_glyph_colors,
        width=2.0,
    )
```

```python
# _launcher.py — in display_controls(), next to SHOW_FLUX_MESH:
state.SHOW_DIRECTORS = sub.checkbox("Show Directors", state.SHOW_DIRECTORS)
```

## What you'll see

### M5.0 wave regime (V=0, no topology)

Glyphs at every Nth voxel oscillate between ±ψ̂_polarization as a wave passes. For the smoke-test x-traveling y-polarized wave: rows of glyphs all swing together between green (+ŷ) and magenta (−ŷ). Traveling-wave fronts visible as bands of synchronized color.

### M5.1 hedgehog (after relaxation)

Concentrated radial pattern of color-varying glyphs in the central ~N/3 of the grid, smoothly transitioning to uniform-color (one direction = vacuum) toward the boundaries. A +1 hedgehog: red on the +x face, cyan on the −x face, green on +y, magenta on −y, blue on +z, yellow on −z — the RGB color cube unfolded onto the grid. A −1 hedgehog: same colors, mirrored across the center.

### M5.4 hedgehog + anti-hedgehog pair

Two color-cube starbursts of opposite polarity, sitting in vacuum with a "compromise zone" between them where the directors fight to satisfy both topologies — visibly distorted glyph orientations along the line connecting them.

## Reference: Exp 2's hedgehog seed extent

The hedgehog seed is **not** purely radial across the entire grid — it's a weighted blend with vacuum that concentrates the radial structure within ~D/4 of each defect (where D = domain edge). From `research_hub/sandbox_phase3_lagrangian/exp2_hedgehog_energy.py:71-108`:

```python
# Per-defect radial proximity weight (long tail, ~1/r):
w_defect = 1.0 / (r + 0.5)

# Vacuum-blend weight (stays ~1 near defects, falls as r⁻⁴ far away):
w_vac = 1.0 / (1.0 + (r_nearest / (DOMAIN / 4)) ** 4)

# Blend: radial near defects, ẑ-vacuum far away
n = w_vac · n_radial + (1 - w_vac) · ẑ

# Renormalize:
n = n / |n|
```

| `r / (D/4)` | `w_vac` | What the seed shows |
| --- | --- | --- |
| 0 (at core) | 1.00 | pure radial |
| 0.5 | 0.94 | still mostly radial |
| 1.0 (at D/4) | 0.50 | half radial / half vacuum |
| 1.5 | 0.16 | mostly vacuum, faint radial tail |
| 2.0 (at D/2 = boundary) | 0.06 | essentially vacuum |

After Frank-energy gradient descent (M5.1 task #4), the relaxed profile typically extends a bit further than the seed (curvature minimization spreads the rotation pattern out) — expect the visible radial structure to fill ~N/3 voxels in radius before fading to uniform-color vacuum at the boundary.

If, after relaxation, you see a *sharp* boundary between the radial zone and vacuum (rather than a smooth transition), the relaxation hasn't converged yet — keep iterating.

## Sequencing

This belongs in M5.1 (alongside `seed_hedgehog` + relaxation), NOT M5.0. Without `seed_hedgehog` there's no interesting director field to render — the smoke-test wave is fine for visual confirmation, but the real motivation for the renderer is "is this hedgehog actually a hedgehog?" which only matters once we have hedgehogs.

Order within M5.1:

1. Implement seeders + Frank energy + relaxation (tasks #1–#4)
2. Implement the renderer (this doc) — visually confirm relaxed hedgehog looks right
3. Run the 1/d Coulomb gating test (task #5) — numerical confirmation
4. Implement winding_number tracker (task #6) — topological confirmation

Renderer between #4 and #5 is ideal: visually catching a bad relaxation before running the (slower) Coulomb sweep saves time when the seed has bugs.
