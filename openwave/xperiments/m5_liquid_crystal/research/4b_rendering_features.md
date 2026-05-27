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
| 1 Deviation (Magnitude) | orange | `ψ.norm()` (+ vector_warp) |
| 2 Thermal Amp (EMA RMS) | ironbow | `amp_local_emarms_am` |
| 3 Thermal Clock (omega) | blueprint | `freq_local_cross_rHz` |
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

### Decisions taken (Rodrigo, 2026-05-26 — M5.4 step 5)

| Decision | Chosen | Implementation |
| --- | --- | --- |
| Granule repurpose | biaxial ellipsoid (showcase) — **deferred to M5.6**; director point-cloud interim now | M5.4 uniaxial → the two minor eigen-axes are DEGENERATE (λ₂=λ₃=δ), so a biaxial ellipsoid's minor axes are arbitrary per voxel. Only a surface-of-revolution reads correctly, and it's a heavy GGUI mesh feature that shows its value only once `δ≠g` (M5.6). `sample_position_to_render` now does the cheap interim: granule at voxel + `amp_boost·n̂`. Full ellipsoid surface → M5.6. |
| Glyph multiplicity | principal-only (forced) | M5.4 is uniaxial — the 3-eigenvector (biaxial) glyph view is degenerate/meaningless until M5.6. Glyphs read `director_nhat`. |
| Arrowheads | **keep** as motion hint | `SHOW_DIRECTOR_ARROWHEAD` unchanged; the half-barb stays even though the director is formally apolar. |
| "Displacement" mode | rename "Orientation deviation" | WAVE_MENU=1 now colors by `‖n̂−ẑ‖` (0 at vacuum, →2 at −ẑ) + perpendicular warp (matches the modes 2–5 scalar pattern; the old 3-axis vector ripple dropped — restore later if wanted). |
| Envelope mode | **defer to M5.7** | M5.4 is static (no oscillation to envelope). Lands when M5.7 dynamics arrive. |

**Granule degeneracy finding (M5.4 step 5):** the biaxial ellipsoid is genuinely premature in uniaxial M5.4 — recommend building the prolate-spheroid → biaxial-ellipsoid *surface* mesh at M5.6 when the structure becomes truly biaxial. The director point-cloud is the faithful interim.

### What does NOT change

- The 3-plane sampling convention (XY at fm_plane_z, XZ at fm_plane_y, YZ at fm_plane_x)
- All color palettes and the colormap module
- `SHOW_FLUX_MESH` / `SHOW_DIRECTORS` / `SHOW_GRANULES` / `VIZ_STRIDE` / `WARP_MESH` controls and their 0..3 plane-reveal semantics
- The flux-mesh vertex/offset/index math
- GGUI `scene.lines` / mesh primitives

Net: the rendering layer is **not a rewrite**. It is "add the eigen-decomposition kernel, recompute the scalar observables from M, re-point ~5 reads." The real M5.4 cost stays in the substrate + operators (matrix triple buffer in `medium.py`, commutator `[M_μ, M_ν]`, matrix Laplacian, matrix seeders, Frank/Coulomb on M).

---

## Part 3 — Wiring viz to physical observables (M5.6.5b): "how do I *see* X?"

**Premise (Rodrigo, 2026-05-27).** Through M5.5 the viz channels were mostly placeholders — only the energy flux-meshes (`WAVE_MENU` 4/5) carried real signal. The channels are general-purpose *displays*; the value comes from **wiring each to the observable that builds physical intuition** for the phenomena OpenWave investigates (matter / forces / EM / heat). This part is the observable→channel map: pick a physics question, pick the channel + observable that answers it visually.

> **Repo-discipline note** (per `feedback_repo_discipline`, cardinal rule): this doc frames the wiring by **physics intuition** only. The thermal observables are the `(A, ω)` excess of `project_thermal_amplitude_hypothesis` (hypothesis *name* used as a physics cross-ref, OK in OpenWave). SABER Direct-Heat-Conversion device specs / engineering targets / product framing — the *applied motivation* — stay in the **private SABER repo**, never here.

### The map — physics question → channel + observable

| Phenomenon | "How do I see it?" | Channel | Observable | Status |
| --- | --- | --- | --- | --- |
| **EM — electric / charge** | where is the charge, `∇·E`? | flux_mesh color (signed greenyellow) + E/director glyph | `∇·n̂` (director splay) — diverges at defect cores like Coulomb charge | ✅ DONE (WAVE_MENU 6 + E glyph) |
| **EM — magnetic / circulation** | where is `B`, the circulation? | flux_mesh color (orange magnitude) + B/curl glyph (half-arrow direction) | `‖∇×n̂‖` (twist+bend) — magnitude colors; the curl-vector glyph shows direction | ✅ DONE (WAVE_MENU 7 + B glyph) |
| **Thermal energy `A`** | where is the heat *amplitude*? | flux_mesh color (ironbow) | `‖M−D‖_F` EMA = thermal **A** | ✅ wired (WAVE_MENU 2) |
| **Thermal clock `ω`** | how fast does it tick (heat↔time)? | flux_mesh color (blueprint) | `‖Ṁ‖_F` EMA = clock **ω** | ✅ wired (WAVE_MENU 3) |
| **Thermal energy (joint)** | the heat *energy* content in one view | flux_mesh color | **`(A·ω)²`** — NOT `A·ω`. For a defect oscillator `E_kin = ½m(Aω)²`, so `A·ω` alone is peak *velocity* (length/time); the energy-dimensional quantity is `(A·ω)²`. Combine the two trackers as a product-of-squares. | 🚧 DEFERRED → 9b |
| **Field energy** | total Hamiltonian / elastic energy | flux_mesh color + scalar warp | `energyH`, `energyF` | ✅ wired (4/5) |
| **Modulation response** | apply an EM-wave lever → see the thermal energy shift | time-trace / before-after of `(A·ω)²` under an EM-wave seed | `Δ(A·ω)²` vs an applied tilt-wave (the SABER-method "modulation" picture, physics framing only) | 🚧 DEFERRED → 9b/M5.7 (dynamic; needs an EM-wave seeder) |

### Channel assignments — what each viz primitive is *best* at

| Channel | Best physical use | Current | Target |
| --- | --- | --- | --- |
| flux_mesh **color** | any scalar field on the 3 planes | energy H/F + thermal A/ω + EM `∇·n̂` (WM6), `‖∇×n̂‖` (WM7) | + `(A·ω)²` joint-thermal (→ 9b) |
| flux_mesh **scalar warp** | lift the plane by a scalar (3D relief of a field) | energy/amp | pair with whichever color scalar is active |
| **vector_warp** | ripple the mesh by a *vector* field — best for showing a field *direction* | director deviation `n̂−ẑ` | superseded for B-direction by the **B/curl glyph** (cleaner than mesh ripple) |
| **director glyphs** | the orientation field itself = the LC "field lines" | principal `n̂` + E/B glyph toggle (shaft = magnitude, half-arrow = direction, color = field/`COLOR_MEDIUM`) | these ARE the EM-analog field lines |
| **granules** | point cloud at voxels; cheap scalar carrier | `voxel + L·n̂` | color each granule by local thermal `A` — a sparse "heat map" (→ 9b) |

### Why `∇·n̂` and `∇×n̂` are the right "see EM" observables

The director `n̂` is the LC analog of the field; its **distortion modes** map onto the EM sector verified in M5.6.4 (`5a §5d`): **splay** `∇·n̂` localizes at defect cores exactly like Coulomb charge density (M5.1 already showed the splay/dumbbell geometry around `±` defects, `3a_coulomb_visual_geometry.md`), and **twist+bend** `∇×n̂` is the circulation/`B`-like mode. Computing `∇·n̂` (signed scalar) and `∇×n̂` (vector + magnitude) as production observables turns the abstract "EM emerges from tilts" result into something you can *look at* on the flux mesh. Cheap: the `∇·`/`∇×` kernels already exist (`divergence`/`curl` on the dormant ψ); they re-point onto `director_nhat`.

### Implementation order (M5.6.5b)

1. ✅ **EM div/curl DONE** (freshest — M5.6.4 just validated the EM sector): `compute_director_em` → `director_div` (scalar), `director_curl` (vector) + `director_curl_mag` (scalar), wired as WAVE_MENU 6/7 + E/B glyphs with size/color toggles. Glyph↔flux-mesh cell-center alignment fixed (see as-built log). **This is the M5.6 PR boundary.**

The rest are **deferred past the M5.6 PR** (Rodrigo 2026-05-27):

1. **gauge-stable charge** — `|∇·n̂|` (unsigned) or topological-winding density, so the charge view doesn't flip sign under Evolve-PDE (kills the apolar-director sign artifact, see caveat below). → lands with **M5.7** (first sustained dynamic runs).
1. **`(A·ω)²` joint-thermal energy** — product-of-squares color mode from the two existing trackers; `(A·ω)²` (not `A·ω`) is the energy-dimensional quantity. → **9b**.
1. **granule heat-map** — color each granule by local thermal `A`. → **9b**.
1. **modulation-response** — apply an EM-wave lever, view `Δ(A·ω)²` (the SABER-method "modulation" picture, physics framing only). → **9b/M5.7** (dynamic; needs an EM-wave seeder).

### As-built log (update as features land)

**✅ EM div/curl — "see charge / see B" (2026-05-27).**

| Aspect | As shipped |
| --- | --- |
| Observable | `compute_director_em` (`engine3_observables.py`): `∇·n̂` → `director_div_field` (signed), `‖∇×n̂‖` → `director_curl_mag_field`. Central-diff stencil, 1-cell halo. |
| Color scale | `compute_director_em_scale`: single-plane `atomic_max` over the 3 render planes → `director_div_absmax`, `director_curl_max` (no full-grid reduction, per `feedback_taichi_metal_atomics`). |
| WAVE_MENU 6 | **EM div (charge/E)** — `∇·n̂` on greenyellow diverging, symmetric scale `[−div_absmax, +div_absmax]` (signed: ± charge). |
| WAVE_MENU 7 | **EM curl (circ/B)** — `‖∇×n̂‖` on orange, scaled against the **shared** distortion magnitude `max(div_absmax, curl_max)` — NOT its own max. |
| Validation | hedgehog `n̂=r̂`: `∇·n̂ = 2/r` (charge), matched to 0.02%; `‖∇×n̂‖ ≈ 0` (radial is curl-free). Headless: render branches compile, finite. |

**The shared-scale design (why curl isn't self-normalized).** A static hedgehog is a *pure charge* — its director is curl-free, so `‖∇×n̂‖` is zero up to ~5e-4 grid-discretization noise. If curl is colored against *its own* max, that noise stretches to full brightness → spurious "rings" around the core (the near-degenerate eigenvector zone amplifies it). Scaling curl against `max(div, curl)` instead makes the B view **honest**: for a static charge `curl/div ≈ 1.6%` → near-black (correctly "no B"); real circulation (a twisting director under Evolve-PDE, or a moving/biaxial-dynamic charge) grows curl to ~div and lights it up. This is the right physics reading — **a static charge has E but no B; B appears with motion/twist**.

**Reading the views.** WAVE_MENU 6 = where the charge is (greenyellow diverging ± at defect cores, neutral bulk). WAVE_MENU 7 = where circulation is (dark for a static defect; lights up where the director twists). Under **Evolve PDE** both evolve in real time — that is *not* a glitch, it is the dynamical field: the charge redistributes and the curl/B grows as the director twists (the circulation generated by the moving/oscillating defect). Pausing freezes them. **This holds with V on too:** `∇·n̂` and `∇×n̂` are *orientation* observables, and V (rotation-invariant, M5.6.5d) pins only the *amplitude* `Tr(M²)`, not the director orientation — so the charge/B views show the freely-evolving twist (the QM/clock sector) regardless of V. To *see* V's confinement on screen, use the amplitude-based views (ENERGY Hamiltonian / Thermal Amp), which stay localized with V on while the orientation views (charge/B/glyphs) still spread.

**⚠️ Sign caveat — `∇·n̂` is gauge-dependent for the apolar director.** The director is nematic (`n̂ ≡ −n̂`, an eigenvector has no head/tail), and `n̂ → −n̂` flips `∇·n̂ → −∇·n̂`. `eigen_decompose` tracks the eigenvector sign **temporally per voxel** (flip if `n̂·n̂_prev < 0`) from the seeded outward `r̂`, so a *static* defect shows a consistent, meaningful ± charge. But under Evolve-PDE the sign continuity is **not enforced spatially** — as neighbouring voxels rotate differently during the slosh, their sign conventions drift apart and `∇·n̂` picks up the mismatch as **spurious local +/− flips** (observed: a clean `+` hedgehog develops green/negative patches after ~35 s of evolution). So: the charge **magnitude / where splay concentrates** is meaningful under dynamics; the **local sign is not a robust observable**. The gauge-invariant, conserved charge is the **topological winding** (Brouwer degree, M5.1 `compute_winding_number`), which never flips. Options if a gauge-stable dynamic charge view is wanted: (a) show `|∇·n̂|` unsigned (loses ± distinction, fine for a single defect); (b) gauge-fix the sign by a defect-relative orientation (`n̂·r̂_defect > 0`); (c) point the "charge" view at the topological winding density. Interim: the signed view stays (best for static intuition), with this caveat.

**✅ B-direction glyph — "see B as a vector" (2026-05-27).** For reading a *direction*, a glyph beats vector_warp (a rippling mesh shows magnitude-ish relief, not direction). `update_em_vector_glyphs` (`engine4_render.py`) draws a line segment along `∇×n̂` (`director_curl_field`) at each 3-plane-sampled voxel, **reusing the director-glyph buffers + offsets** (no new buffer/render path). Three glyph channels, deliberately assigned: **shaft length = magnitude**, **arrow (half-barb) = direction**, **color = the field value** (matching the corresponding flux-mesh WAVE_MENU). A glyph has 3 independent channels, so once the shaft carries magnitude, color is *free* — no information is lost by also encoding a scalar in it.

| Glyph (UI: "EM field Glyph (off=E on=B)") | shaft | arrow barb | color | matches mesh |
| --- | --- | --- | --- | --- |
| **E / director** (`GLYPH_VECTOR=0`) | uniform (n̂ is a *unit* field — its "magnitude" is the charge, shown by color) | director direction | `∇·n̂` signed → **greenyellow** | WM6 |
| **B / curl** (`GLYPH_VECTOR=1`) | `∝ min(‖∇×n̂‖/scale,1)` (genuine vector) | curl direction | `‖∇×n̂‖` → **orange** | WM7 |

Both reuse the one director-glyph buffer (shaft + arrowhead) + 3-plane sampling, switched by the toggle. Shared scale `max(div_absmax, curl_max)`. Headless-verified: both compile, finite; the B glyph auto-declutters — a static charge (curl≈0) gives ~zero-length, near-black arrows (no circulation = no B); real twist grows visible orange arrows with barbs showing the circulation direction. Safe `/(norm+eps)` avoids NaN at curl≈0.

**Two decoupling toggles (2026-05-27) — separate "field strength/value" from "see-everywhere structure".** The shaft/color assignment above is the *default*; two checkboxes below "EM field Glyph" let either channel be overridden:

| Toggle | off (0) | on (1) |
| --- | --- | --- |
| **Glyph Size (unit, magnitude)** | **unit** — shaft = `length·dir` (every glyph visible; the field-*line/direction* structure everywhere, incl. far field) | **magnitude** — shaft `∝` field magnitude (E: `\|∇·n̂\|` charge density, B: `‖∇×n̂‖`); declutters to where the field is strong |
| **Glyph Color (single, gradient)** | **single** — flat `COLOR_MEDIUM` (sourced from `colormap.py`, light blue; module constant `_GLYPH_SINGLE_COLOR`, no hardcoded RGB) — every glyph the same color so weak/far-field glyphs stay visible (gradient fades them to ~black) | **gradient** — the field value (E: `∇·n̂` greenyellow, B: `‖∇×n̂‖` orange), matching the WM6/WM7 mesh |

The `unit`+`single` combination is the **far-field inspection** view (uniform, fully-visible field lines everywhere); `magnitude`+`gradient` is the **field-strength** view (where + how strong + which way). The barb scales with the shaft, so it shrinks with the glyph in magnitude mode. Headless-verified across all 4 combinations (finite).

**✅ Glyph ↔ flux-mesh alignment fix (2026-05-27, on-screen confirmed).** Glyph origins didn't sit on the flux-mesh cells: the glyph used the bare voxel node `i/max_dim`, while `create_flux_mesh` (`medium.py`) places every mesh vertex at the **cell center** `(i+0.5)/max_dim`, and the glyph's perpendicular coord used a lossy `int(plane·n)/max_dim` that diverges from the mesh's continuous `plane·n/max_dim` on odd-dimensioned axes (grids are auto-derived from universe-size ÷ target-voxels, so odd counts are common). Fix: glyphs now read the mesh's **own** coordinates — `+0.5` on the two in-plane axes + the new continuous `fm_plane_{x,y,z}_pos` (added to `medium.py`, same formula the mesh uses) for the perpendicular. The field is still *sampled* at the integer voxel; only the render *placement* moved (exactly what the mesh does). Result: glyph tails sit dead-center on each mesh cell and flat in the sheet, for any grid parity / plane fraction. One source of truth for "where voxel `i` is in normalized space."

**Why the B color is *magnitude* (orange), not bluered-N/S — and what we'd lose.** `‖∇×n̂‖` is a magnitude (≥0), so it *cannot* show N/S poles. Red=N/blue=S needs a **signed** scalar — a projection `B·axis` (`B·ẑ` or the dipole axis). That's deferred (not lost): (1) it's axis-relative — only reads as N/S for a dipole-aligned field; (2) **there is no magnetic dipole yet** — the static hedgehog is a pure electric charge (`B=∇×n̂≈0`), so bluered-N/S would color near-zero noise. It becomes the right choice once a configuration produces a real circulating B → see **M5.6.5f (magnetic-dipole viz xperiment)** in `0b_M5_roadmap.md`.

**🚧 Deferred past the M5.6 PR (Rodrigo 2026-05-27)** — none blocks the PR; each homed to the phase where it becomes meaningful:

| Feature | What | Home | Why deferred |
| --- | --- | --- | --- |
| **gauge-stable charge** | `\|∇·n̂\|` unsigned or topological-winding density, so the charge view doesn't sign-flip under Evolve-PDE (apolar caveat above) | M5.7 | Needs sustained dynamic runs to be worth it |
| **`(A·ω)²` joint-thermal** | energy-dimensional color mode. `(A·ω)²` not `A·ω` — `E_kin=½m(Aω)²`, so `A·ω` alone is peak *velocity*; the energy quantity is `(A·ω)²` | 9b | Thermal program; the (A,ω) hypothesis test |
| **granule heat-map** | color each granule by local thermal `A` (sparse heat map) | 9b | Thermal program |
| **modulation-response** | apply an EM-wave lever → view `Δ(A·ω)²` (SABER-method "modulation", physics framing only) | 9b/M5.7 | Dynamic; needs an EM-wave seeder |
| **B dipole N/S (bluered)** | signed `B·axis` red=N/blue=S | M5.6.5f@M5.8 | Needs a real circulating B (no dipole yet) |
