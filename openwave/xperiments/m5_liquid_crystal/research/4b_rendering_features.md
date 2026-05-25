# Rendering / Visualization Features — Catalog + Matrix-Substrate Migration

**Purpose**: single source of truth for OpenWave M5's rendering/visualization stack — the **biggest differential OpenWave has as a simulator** (per `feedback_visual_rendering_priority`). Part 1 catalogs the as-shipped feature set (M5.1, Vector(3) ψ). Part 2 is the M5.4 repurposing plan when the substrate becomes the matrix field `M = ODO^T`. None of these features is dropped in the migration; each is re-sourced.

> Director-glyph **design + as-shipped specifics** live in [`2b_director_glyph_rendering.md`](2b_director_glyph_rendering.md). This doc is the broader stack catalog + migration plan. M5.4 implementers: Part 2 is the work list; cross-ref [`0b_M5_roadmap.md` § Phase M5.4](0b_M5_roadmap.md#phase-m54--matrix-field-substrate-migration).

---

## Part 1 — Current feature inventory (M5.1, Vector(3) ψ)

### The rendering features

| Feature | Code | Behavior |
| --- | --- | --- |
| flux_mesh | `update_flux_mesh_values`, `render_flux_mesh`, `fluxmesh_*_vertices/colors` | 3 plane meshes (XY/XZ/YZ); each voxel → 1 scalar → color + perpendicular warp |
| scalar warp_mesh | WAVE_MENU 2–5 | raises mesh perpendicular coord by a scalar (amp/freq/energy) |
| vector_warp | WAVE_MENU 1 | deforms vertex in all 3 axes by the ψ components (mesh ripples by the vector field) |
| granule render | `sample_position_to_render`, `position_render` | granule sphere sits at voxel + `amp_boost · ψ` (`amp_boost = WARP_MESH`) |
| vector glyphs | `update_director_glyphs`, `director_glyph_*` | line voxel → voxel + L·n̂, colored by `(1 − n̂_z)`, optional half-arrowhead |

### Color-coding modes (WAVE_MENU)

| Mode | Palette | Reads |
| --- | --- | --- |
| 1 Displacement (Magnitude) | orange | `ψ.norm()` (+ vector_warp) |
| 2 Amplitude (EMA RMS) | ironbow | `amp_local_emarms_am` |
| 3 Frequency (L&T) | blueprint | `freq_local_cross_rHz` |
| 4 ENERGY (Hamiltonian) | ironbow | `energyH_density_aJ` |
| 5 ENERGY (Frank Elastic) | ironbow | `energyF_density_aJ` |

Palettes (`get_orange_color`, `get_ironbow_color`, `get_blueprint_color`) live in the colormap module and are substrate-agnostic.

### Controls

| Control | Range | Drives |
| --- | --- | --- |
| `SHOW_FLUX_MESH` | 0..3 | flux-mesh plane reveal (0 off, 1 XY, 2 +XZ, 3 all) |
| `SHOW_DIRECTORS` | 0..3 | glyph plane reveal (same semantics) |
| `SHOW_GRANULES` | bool | granule point-cloud (only when flux mesh active) |
| `VIZ_STRIDE` | int | shared sampling stride for glyphs AND granules |
| `WARP_MESH` | float | warp magnitude (scalar warp + granule `amp_boost`) |
| `WAVE_MENU` | 1..5 | active color-coding mode |

---

## Part 2 — Matrix-substrate migration (M5.4+)

When M5.4 replaces the Vector(3) director ψ with the paper's real-symmetric matrix field `M(x) = O(x)·D·O^T` (`D = diag(g, 1, δ)` frozen global; `O(x)` the per-voxel dynamical rotation), **every** rendering feature above must be re-sourced. None is dropped.

### The two-role split — why each feature lands differently

The Vector(3) ψ does **double duty** today, and only one role survives:

| ψ role today | Consumers | Survives? |
| --- | --- | --- |
| Wave displacement (M5.0) | Gaussian/dispersion seeds, vector_warp, granule push, \|ψ\|-amplitude, ψ_z-frequency | retires with the wave regime |
| Director orientation (M5.1+) | glyphs, Frank energy, winding Q, Coulomb relax | survives, re-sourced from M |

In the matrix world there is **no displacement vector** — M stores orientation. The director becomes a *derived* quantity: the eigenvector of M, not the stored field.

### The lynchpin — one new kernel unlocks the whole stack

A single new primitive carries nearly all the repurposing:

```text
eigen_decompose(M)  →  (n̂₁, n̂₂, n̂₃, λ₁, λ₂, λ₃)   @ti.func, per voxel
```

It produces a `director_nhat` Vector(3) (principal eigenvector) + the three eigenvectors/eigenvalues. Once it exists: glyphs read `director_nhat`, vector_warp reads it, granule-director reads it, and the eigenvalue spread drives the ellipsoid option. The scalar observable kernels (`‖M−D‖_F`, energyH/F recomputed from M) feed flux_mesh, scalar warp, and the color modes. **Build this kernel early in M5.4 — it is the rendering and tracker dependency**, and it is folded into the M5.3 feasibility spike so it is proven in Taichi before the production migration.

### Feature-by-feature mapping

| Feature | Verdict | New source |
| --- | --- | --- |
| flux_mesh (3-plane sampler) | ✅ survives unchanged | any matrix-derived scalar field |
| scalar warp_mesh | ✅ survives | `‖M−D‖_F`, `Tr(M²)`, energyH/F, clock ω |
| vector_warp | ⚠️ repurpose | warp by `director_nhat` deviation from vacuum ẑ |
| granule render | ❌ dies as displacement → repurpose | director point-cloud, or biaxial ellipsoid per O(x) |
| vector glyphs | ✅ survives, gets richer | principal eigenvector(s) of M |

Details per feature:

**flux_mesh** — substrate-agnostic. A scalar→color/warp driver: samples one scalar per voxel on the 3 planes and maps it; indifferent to whether the scalar came from ψ or M. The 3-plane sampling, offset/index math, and palette calls all carry over verbatim. Only the *scalar fields it reads* change (recomputed from M).

**scalar warp_mesh** — survives. Raises each vertex's perpendicular coordinate by a scalar. Feed it any matrix scalar: orientation deviation `‖M−D‖_F`, `Tr(M²)`, recomputed energy densities, or clock frequency. `WARP_MESH` slider semantics unchanged.

**vector_warp** — repurpose. No displacement vector to ripple by. Re-source the warp 3-vector from the **director deviation** `n̂ − ẑ` (or the EM-tilt-gradient vector once M5.5 lands). The mesh bulges/ripples where directors tilt away from vacuum — preserving the "mesh deforms by a vector field" aesthetic, now showing *orientation* rather than displacement.

**granule render** — dies in current form ("granule pushed by the displacement wave" is a wave-model holdover; an orientation field has nothing to push by). Three repurpose options, increasing fidelity:

| Option | Shows | Cost |
| --- | --- | --- |
| (a) retire | granules off in matrix xperiments | trivial |
| (b) director point-cloud | granule at voxel + L·n̂ | low (reuses `director_nhat`) |
| (c) biaxial ellipsoid | per-granule ellipsoid, axes scaled by (g,1,δ), oriented by O(x) | medium |

Option (c) is the richest and most LC-faithful — the literal "biaxial top at each voxel" picture from the Duda thread study ([`4a_convo_2026.05.12.md` §5](4a_convo_2026.05.12.md)). Recommend (b) as cheap interim, (c) as the showcase.

**vector glyphs** — survives, gets richer. Source becomes the principal eigenvector of M instead of `ψ/‖ψ‖`. Biaxial M has up to **three** eigenvector families → option to render three glyph sets, color-coded by eigenvalue (g / 1 / δ axes). The arrowhead (`SHOW_DIRECTOR_ARROWHEAD`) should **retire**: a nematic director is apolar (n ≡ −n), so an eigenvector has no head/tail — correct LC physics, not a regression. The 2b Option C ("glyphs only near defects") becomes more valuable for M5.4 multi-defect scenes.

### Color-coding modes — source remap

Palettes survive untouched; only the scalar each mode reads changes:

| Mode | Current source | Matrix source |
| --- | --- | --- |
| Displacement | `ψ.norm()` | `‖M−D‖_F` (orientation deviation) |
| Amplitude (EMA RMS) | EMA of `\|ψ\|²` | EMA of `‖M−D‖_F` — thermal **A** |
| Frequency (L&T) | ψ_z zero-crossing | O(x) rotation rate at core — clock **ω**, thermal ω |
| Envelope (new) | n/a | slow-varying envelope of the `‖M−D‖_F` oscillation |
| ENERGY (Hamiltonian) | kinetic+grad+V on ψ | recomputed from M: `‖Ṁ‖² + ‖∇M‖² + V(M)` |
| ENERGY (Frank Elastic) | `(K/2)\|∇n̂\|²` | matrix elastic term from `∇M` / commutator |

Two notes:

- **Displacement → orientation-deviation**: the orange "displacement" mode loses its literal meaning but stays useful renamed "Orientation deviation" — how far the local frame is twisted from vacuum D.
- **Envelope is a new mode worth adding**: the slowly-varying amplitude profile of the orientation-excitation oscillation. It is the natural view for the **M5.7 resonance hunt** — an envelope that stays localized over many clock periods *is* a long-lived resonance, visually. Pairs with the redefined frequency tracker.

### Tracker redefinitions (feed the color modes)

Amplitude/frequency trackers are not retired — re-derived against M, and physically sharper:

| Tracker | New matrix definition | Doubles as |
| --- | --- | --- |
| Amplitude | Frobenius deviation `‖M−D‖_F` (EMA) | thermal **A** (SABER joint (A, ω)) |
| Frequency | O(x) rotation rate at the defect core | de Broglie clock **ω** (M5.8 headline), thermal ω |

Tracker *infrastructure* (per-voxel EMA fields, crossing detection, global aggregates) carries over unchanged; only the sampled quantity changes. The ψ_z zero-crossing was always a convention hack (the kernel docstring admits it) — the matrix clock-rate is the genuine article.

### Decisions that need Rodrigo's input (UX calls, not mechanical ports)

| Decision | Options |
| --- | --- |
| Granule repurpose | retire / director point-cloud / biaxial ellipsoid |
| Glyph multiplicity | principal eigenvector only / all 3 eigenvector families |
| Arrowheads | retire (apolar director) / keep as motion hint |
| "Displacement" mode | rename to "Orientation deviation" / drop |
| Envelope mode | add now (M5.7-ready) / defer |

### What does NOT change

- The 3-plane sampling convention (XY at fm_plane_z, XZ at fm_plane_y, YZ at fm_plane_x)
- All color palettes and the colormap module
- `SHOW_FLUX_MESH` / `SHOW_DIRECTORS` / `SHOW_GRANULES` / `VIZ_STRIDE` / `WARP_MESH` controls and their 0..3 plane-reveal semantics
- The flux-mesh vertex/offset/index math
- GGUI `scene.lines` / mesh primitives

Net: the rendering layer is **not a rewrite**. It is "add the eigen-decomposition kernel, recompute the scalar observables from M, re-point ~5 reads." The real M5.4 cost stays in the substrate + operators (matrix triple buffer in `medium.py`, commutator `[M_μ, M_ν]`, matrix Laplacian, matrix seeders, Frank/Coulomb on M).
