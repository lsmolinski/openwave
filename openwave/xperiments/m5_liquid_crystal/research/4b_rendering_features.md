# Rendering / Visualization Features — Catalog + Matrix-Substrate Migration

**Purpose**: single source of truth for OpenWave M5's rendering/visualization stack — the **biggest differential OpenWave has as a simulator** (per `feedback_visual_rendering_priority`). Part 1 catalogs the as-shipped feature set (M5.1, Vector(3) ψ). Part 2 is the M5.4 repurposing plan when the substrate becomes the matrix field `M = ODO^T`. None of these features is dropped in the migration; each is re-sourced.

> Director-glyph **design + as-shipped specifics** live in [`2b_director_glyph_rendering.md`](2b_director_glyph_rendering.md). This doc is the broader stack catalog + migration plan. M5.4 implementers: Part 2 is the work list; cross-ref [`0b_M5_roadmap.md` § Phase M5.4](0b_M5_roadmap.md#phase-m54--matrix-field-substrate-migration).
>
> **Reading order (2026-05-30):** Parts 1–3 are the as-built history (M5.1 inventory → M5.4 migration → M5.6.5b EM-wiring as-shipped). **Part 4 is the consolidated target** — the full "what we want to *see*" observable catalog for the thermal/EM research, organized from Rodrigo's 2026-05-30 vision. **Part 5 is the implementation timeline + the M5.8 4×4 safety contract + the placeholder-sample strategy** (hard-coded analytic fields to validate rendering *before* the physics produces it). Start at Part 4 for the plan; Parts 1–3 for why things are the way they are.

---

## Part 1 — Current feature inventory (M5.1, Vector(3) ψ)

### The rendering features

| Feature | Code | Behavior |
| --- | --- | --- |
| flux_mesh | `update_flux_mesh_values`, `render_flux_mesh`, `fluxmesh_*_vertices/colors` | 3 plane meshes (XY/XZ/YZ); each voxel → 1 scalar → color + perpendicular warp |
| scalar warp_mesh | WAVE_MENU 2–5 | raises mesh perpendicular coord by a scalar (amp/freq/energy) |
| vector_warp | WAVE_MENU 1 | deforms vertex in all 3 axes by the ψ components (mesh ripples by the vector field) |
| vector glyphs | `update_director_glyphs`, `director_glyph_*` | line voxel → voxel + L·n̂, colored by `(1 − n̂_z)`, optional half-arrowhead |
| granule render | `sample_position_to_render`, `position_render` | granule sphere sits at voxel + `amp_boost · ψ` (`amp_boost = WARP_MESH`) |

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
| **EM — magnetic / rotation** | where is `B`, the rotational? | flux_mesh color (orange magnitude) + B/curl glyph (half-arrow direction) | `‖∇×n̂‖` (twist+bend) — magnitude colors; the curl-vector glyph shows direction | ✅ DONE (WAVE_MENU 7 + B glyph) |
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

1. **gauge-stable charge** — UPDATE 2026-05-30 (VIZ.1): the **director-glyph** half is ✅ DONE (centered + barbless = gauge-stable). The **signed-charge WM6** half is **left honest-but-flipping** and deferred to **M5.8** via topological-winding density (`|∇·n̂|`-unsigned was tried + dropped as redundant). The charge-region *expansion* is real free-defect dispersal, not the artifact. See §4.4 + §5.2 #1b.
1. **`(A·ω)²` joint-thermal energy** — product-of-squares color mode from the two existing trackers; `(A·ω)²` (not `A·ω`) is the energy-dimensional quantity. → **9b**.
1. **granule heat-map** — color each granule by local thermal `A`. → **9b**.
1. **modulation-response** — apply an EM-wave lever, view `Δ(A·ω)²` (the SABER-method "modulation" picture, physics framing only). → **9b/M5.7** (dynamic; needs an EM-wave seeder).

### As-built log (update as features land)

**✅ EM div/curl — "see charge / see B" (2026-05-27).**

| Aspect | As shipped |
| --- | --- |
| Observable | `compute_director_em` (`engine3_observables.py`): `∇·n̂` → `director_div_field` (signed), `‖∇×n̂‖` → `director_curl_mag_field`. Central-diff stencil, 1-cell halo. |
| Color scale | `compute_director_em_scale`: single-plane `atomic_max` over the 3 render planes → `director_div_absmax`, `director_curl_max` (no full-grid reduction, per `feedback_taichi_metal_atomics`). |
| WAVE_MENU_6 | **EM div (charge/E)** — `∇·n̂` on greenyellow diverging, symmetric scale `[−div_absmax, +div_absmax]` (signed: ± charge). |
| WAVE_MENU_7 | **EM curl (circ/B)** — `‖∇×n̂‖` on orange, scaled against the **shared** distortion magnitude `max(div_absmax, curl_max)` — NOT its own max. |
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
| **gauge-stable charge** | director glyph (centered+barbless) ✅ DONE (VIZ.1); signed-charge WM6 deferred to M5.8 via winding-density | **glyph: ✅ done; signed charge: M5.8** | re-planned 2026-05-30: glyph centering shipped; WM6 left honest-but-flipping (expansion is real physics) — see §4.4 + §5.2 #1b |
| **`(A·ω)²` joint-thermal** | energy-dimensional color mode. `(A·ω)²` not `A·ω` — `E_kin=½m(Aω)²`, so `A·ω` alone is peak *velocity*; the energy quantity is `(A·ω)²` | 9b | Thermal program; the (A,ω) hypothesis test |
| **granule heat-map** | color each granule by local thermal `A` (sparse heat map) | 9b | Thermal program |
| **modulation-response** | apply an EM-wave lever → view `Δ(A·ω)²` (SABER-method "modulation", physics framing only) | 9b/M5.7 | Dynamic; needs an EM-wave seeder |
| **B dipole N/S (bluered)** | signed `B·axis` red=N/blue=S | **sample pre-M5.8 (Part 5 #4), real @M5.8 (#7)** | re-planned 2026-05-30: the *render* is buildable now against a placeholder dipole (§4.5/§5.3); only the real circulating-B *source* waits for M5.8 |

> **Re-plan note (2026-05-30):** Part 5 supersedes the homes above for the two EM carry-overs.
> **gauge-stable charge** and the **magnetic-dipole viz** both have a pre-M5.8 path now (§4.4 and
> §4.5 with the placeholder strategy §5.3). The curl-vector mesh-warp (§4.3) and the 3-way glyph
> select (§4.2) join them as the pre-M5.8 batch. See Part 5 §5.2 for the full ordered timeline.

---

## Part 4 — The consolidated viz target (Rodrigo 2026-05-30): "what I want to *see*"

This part organizes the full observable wishlist for the EM + thermal research into one coherent
catalog: **the physical quantities to display, what each *is* in the substrate, and the render
channel that conveys it.** It is the ground-work for the SABER/thermal experiments (9b) — the
channels are how we will *read* thermal reality and *watch* modulation methods work in the
simulator. Physics-framing only here (cardinal rule); device/engineering motivation stays in the
private SABER repo.

### 4.1 The observable catalog — what each quantity IS

Everything visible derives from the director `n̂` (principal eigenvector of `M`) and the matrix
deviation. Two source families: **EM = orientation distortion of `n̂`** (tilts/gradients), **thermal
= the clock's rotational state** (`(A, ω)` of the twist). The director itself is the origin of both
EM observables (`∇·n̂`, `∇×n̂`) — Rodrigo's read is correct.

| Quantity | What it IS in the substrate | Scalar / vector | Source field (live) | Status |
| --- | --- | --- | --- | --- |
| **Director `n̂`** | principal eigenvector of `M` — the LC "grain"; apolar (`n̂≡−n̂`), unit length | unit vector | `director_nhat` | ✅ live |
| **Charge density** | `∇·n̂` (director splay) — diverges at defect cores like Coulomb charge | signed scalar | `director_div_field` | ✅ live (WM6) |
| **Charge sign** | sign of `∇·n̂` (+ outward / − inward hedgehog) — color-coded green=−/yellow=+ | sign | `director_div_field` sign | ⚠️ flips under Evolve-PDE → §4.4 fix |
| **Electric field lines** | the director field itself (E ∝ the `n̂` texture); convention: lines point **+ → −** | vector (= `n̂`) | `director_nhat` | ✅ live (E-glyph); +→− convention = §4.2 |
| **Magnetic field strength** | `‖∇×n̂‖` (twist+bend circulation magnitude) — ≈0 for a static charge | scalar ≥0 | `director_curl_mag_field` | ✅ live (WM7) |
| **Magnetic field lines** | the curl vector `∇×n̂` — direction = circulation; handedness N→S | vector | `director_curl_field` (**raw vector already stored**) | ✅ vector glyph; mesh-warp = §4.3 |
| **Magnetic moment μ** | net `∮ B` of the defect / its dipole axis — a single vector per defect | vector (per defect) | derived (reduction over `director_curl_field`) | 🚧 §4.5 (needs a circulating B) |
| **Thermal amplitude `A`** | `‖M−D‖_F` EMA = the clock's rotational **radius** (`≈δ/2`), grows with heat | scalar | `amp_local_emarms_am` | ✅ live (WM2) |
| **Thermal clock `ω`** | `‖Ṁ‖_F` EMA = the Zitterbewegung rate (director angular speed) | scalar | `freq_local_cross_rHz` | ✅ live (WM3) |
| **Joint thermal `(A·ω)²`** | energy-dimensional heat content (`E_kin=½m(Aω)²`; `A·ω` alone is peak velocity) | scalar | product of the two trackers | 🚧 9b |

> **Zitterbewegung note (from `0c §L7`):** the thermal `A` is the *radius* of the spinning
> director (constant at ground state, `≈δ/2`), `ω` its rate — the two are the AM/FM channels. WM2 +
> WM3 already show them; the joint `(A·ω)²` view is the single-scalar heat-energy display (→ 9b).

#### 4.1.1 How to *see* the Zitterbewegung clock — the spinning director, its ω, its radius

This is the central thermal observable (and the M5.8 headline), so it gets an explicit
"how-do-I-see-it" recipe. The clock = the director frame **twisting about the director axis `a`**
(the twist generator `Gx` rotates the `b`–`c` plane, `0c §L7` + `5a §7a`); its `ω` is the spin
rate, its **radius** is the rotational amplitude `≈δ/2` (the thermal `A`).

> ⚠️ **Correction (Rodrigo 2026-05-30) — the director glyph does NOT show the clock spin.** Two
> independent reasons: **(1)** the clock twists the secondary axes *around* `a`, so the director
> `a` itself **stays put** (it's the axle, not the clock-hand — invariant under the twist).
> **(2)** Even if it rolled about its own axis, a single line segment looks *identical at every
> roll angle* — a line has no feature to track roll, and the **barb shows head-vs-tail (direction
> along the shaft), never roll-about-the-shaft**. So removing/keeping the arrowhead is irrelevant
> to the spin. What you see sloshing on the director glyph under Evolve-PDE is the director
> **tilting** (`a` changing *direction* — the EM/tilt sector), **not** the clock-twist. To *see*
> the spin you must render the **secondary (δ) axis** sweeping around `a`.

Ways to watch the clock, corrected:

| Want to see | Channel | What you watch | Status |
| --- | --- | --- | --- |
| **the spin** (the twist about `a`) | **the δ cross-bar in the "Director + Delta" glyph** — the **middle (δ) eigenvector** (`director_mid`), shown as the shorter CYAN arm of the ellipsoid-wireframe "+" | the `b`/δ axis *sweeps around* the director `a` — would be the visible clock-hand if coherent. **In free 3D it tilts/disperses (M5.7.2), not coherently spins** — coherent spin needs M5.8/9b. Non-degenerate only on biaxial `diag(1,δ,0)`. | ✅ **DONE (VIZ.3, glyph state 1 "Director + Delta")** |
| spin (alt, cheap) | **perpendicular tick** on the director glyph | add one short barb *perpendicular* to the shaft, oriented by the δ-axis — a mark whose angle tracks the roll. Cheapest "see the roll" hack. | 🚧 needs building |
| **the rate `ω`** (how fast it ticks) | flux_mesh **WM3 "Thermal Clock"** (blueprint) | color = `‖Ṁ‖_F` EMA = angular speed; brighter = faster. *Measures* the rate (doesn't show rotation). | ✅ live (WM3) |
| **the radius `A`** (rotational amplitude) | flux_mesh **WM2 "Thermal Amp"** (ironbow) | color = `‖M−D‖_F` EMA = how far the frame swings from vacuum = the `≈δ/2` radius; grows with heat. | ✅ live (WM2) |
| **the joint energy** `(A·ω)²` | flux_mesh color (new mode) | single scalar = heat-energy content (`E_kin=½m(Aω)²`) | 🚧 9b |
| **the orbit traced** (optional) | granule **Zitterbewegung tracer** (§4.6) | a granule traces the **δ-axis tip's** orbit around `a` over time — makes the spin *visible as a path* | 🚧 9b candidate |

**Key reading:** the **secondary-δ-axis glyph** is the only channel that *visually shows the spin*
(the clock-hand sweeping around the director-axle); the director glyph shows tilt, not twist. Under
**Evolve PDE**, WM2 (radius) and WM3 (rate) are the AM and FM channels — together they *are* the
joint `(A, ω)` thermal state (heating pumps WM2/AM and/or shifts WM3/FM; M5.7.3 saw both respond).
WM2/WM3 already exist; the δ cross-bar in the Director glyph (VIZ.3 ✅, the shorter CYAN arm of the
ellipsoid "+") is the geometric channel that shows the δ-axis rotation directly.

> ⚠️ **3D caveat (`0c §L7`):** the *free* clock disperses in 3D (M5.7.2) — what you see spinning
> then radiating away is the free mode losing coherence. The self-sustaining clock is the **4D**
> Zitterbewegung (M5.8); a **driven** one (9b) holds a steady `(A, ω)`. So pre-M5.8 these channels
> show the *transient* spin + dispersal, which is itself the informative picture.

### 4.2 Glyph (vector-field) displays — split the one EM toggle into four

Today: one `GLYPH_VECTOR` toggle (off=E/director, on=B/curl) + Glyph-Size + Glyph-Color decoupling
checkboxes (Part 3 as-built). **Target:** make the vector fields independently selectable — they
answer different questions and a researcher wants to compare them. The **δ-clock-hand** is the new
4th channel (the only glyph that shows the spin — `§4.1.1` correction).

**4-state select (as-built VIZ.3, refined 2026-05-30): 4 mutually-exclusive checkboxes —
`0 = Director Vector`, `1 = Director + Delta Vectors`, `2 = Electric Field`, `3 = Magnetic Field`.**
A slider (`0=n+d 1=E 2=B`) was the first cut; switched to checkboxes because the labels read
clearer than slider positions, and to add a **Director-only** state so the δ bar never shows alone
(when the frame gets cluttered, drop to `n̂` axis only). The `Glyph Size`/`Glyph Color` toggles apply
**only to the field states (E/B)** — the Director states are *orientation, not a field*, so they
ignore them (always unit + fixed colors). UI spells "delta" — GGUI cannot render the Greek `δ`.

| Glyph state | Direction (barb) | Size / Color | What it is |
| --- | --- | --- | --- |
| **0 Director Vector** | none (apolar axis) | **fixed** (ignores size/color) — `n̂` unit + COLOR_MEDIUM | principal axis `a`=`director_nhat` only — the clean orientation arrow (arrow buffer blanked) |
| **1 Director + Delta** | none (both axes apolar) | **fixed** (ignores size/color) — `n̂` unit + COLOR_MEDIUM; δ bar shorter (∝ λ₂/λ₁) + CYAN/COLOR_FIELD | the biaxial-frame **ellipsoid-wireframe "+"**: `a`=`director_nhat` + middle axis `b`=`director_mid`. Shows orientation (tilt); the δ bar is the would-be clock-hand (free 3D = tilts/disperses, M5.7.2 — coherent spin needs M5.8/9b) |
| **2 Electric Field** | `+ → −` half-barb (gauge-arbitrary sign until M5.8) | **honors** size (`∝\|∇·n̂\|` charge) + color (greenyellow charge) | `director_nhat` as a polar field line (E ∝ `n̂`) |
| **3 Magnetic Field** | `N → S` (along `∇×n̂`) | honors size (`∝‖∇×n̂‖`) + color | `director_curl_field`, via `update_em_vector_glyphs` |

**Why the Director states ignore size/color (Rodrigo 2026-05-30):** they represent the *ellipsoid
orientation*, not an E/B field — there is no charge/strength to encode, so magnitude-scaling and
charge-coloring are meaningless. State 0 is the bare director axis; state 1 adds the δ cross-bar to
make the biaxial frame legible as a wireframe "+". (Both share `mode=0` in `_write_glyph`; a
`show_delta` flag toggles the δ arm — blanked in state 0, drawn in state 1.)

The **Size/Color toggles** (`Glyph Size unit/magnitude`, `Glyph Color single/gradient`) act only on
the **E/B field states**: `unit` = every line same length (structure-everywhere view) vs `magnitude`
= shaft `∝` the field value (charge density `\|∇·n̂\|` for E, strength `‖∇×n̂‖` for B); `single`
COLOR_MEDIUM vs the value gradient. The Director state is fixed (orientation, not a field).

**Centering convention (decided 2026-05-30) — ALL glyphs centered on the voxel; the barb is the
only polar/apolar distinction.** A glyph is a short *tangent-segment of the field line through the
voxel*, so centered (`base = pos − ½·shaft·v̂ → tip = pos + ½·shaft·v̂`) is the faithful
"field-line" representation for *all* of them. The barb is what differs:

| Glyph | Centered? | Barb? | Why |
| --- | --- | --- | --- |
| **Director `n̂`**, **δ-clock-hand** | ✅ yes | ❌ no | apolar axes (`v̂ ≡ −v̂`) — pass *through* the point, no head/tail (centering is also a free gauge-fix, VIZ.1) |
| **E** (`+→−`), **B** (`N→S`) | ✅ yes | ✅ yes | polar field lines — pass *through* the voxel (so centered), but have a real direction (so barbed) |

Director centering is already **shipped** (VIZ.1); **E/B centering lands in VIZ.3** when their
glyphs are rebuilt (the current shipped B-glyph is base-at-voxel — switch to centered for field-line
faithfulness + visual consistency).

### 4.3 flux_mesh (scalar/warp) displays — keep, with one upgrade

The WAVE_MENU flux-meshes are good as-is. One change Rodrigo flagged:

| WAVE_MENU | Keep / change | Detail |
| --- | --- | --- |
| 6 **EM div** (`∇·n̂`) | **keep** | signed charge scalar on greenyellow diverging — correct as a scalar field |
| 7 **EM curl** | **upgrade (warp + color toggle)** | currently warps the mesh vertex by `‖∇×n̂‖` (a *scalar* → perpendicular-only lift) on orange. **Change (warp):** warp the vertex by the **raw `∇×n̂` vector** (`director_curl_field`, already stored) → the mesh deforms as a *twist in fabric*, showing the B-field rotation + handedness, not just its magnitude (the "vector_warp" idea from Part 2, now sourced from the curl vector). **Change (color, optional):** add a color toggle — `orange ‖∇×n̂‖` (magnitude, default/honest) **or** `bluered (∇×n̂)·axis` (**signed** axial projection = N/S poles, §4.5). |
| 2/3 **Thermal A/ω** | keep | the wired thermal channels (the Zitterbewegung clock — §4.1.1) |
| 4/5 **Energy H/F** | keep | |

> **Why the curl-vector warp is cheap:** `director_curl_field` (the raw `∇×n̂` vector) is *already
> computed and stored* (`engine3_observables.py:417`) — WM7 just throws away the direction and uses
> `.norm()`. The upgrade is a new render branch that displaces the vertex by all 3 curl components
> (like WAVE_MENU 1's old vector_warp did for ψ), no new physics kernel.

**Blue-red signed-color — warp and color are separable** (Rodrigo 2026-05-30):

> The *warp*
> (twist-in-fabric) shows B rotation/handedness from the **vector** `∇×n̂` and works regardless of
> magnitude. The *color* is an independent scalar channel: the current `‖∇×n̂‖` is a **magnitude
> (≥0) → cannot show N/S poles**; to color by pole you need a **signed** scalar — the axial
> projection `B·axis = (∇×n̂)·axis` (red = +axis/N, blue = −axis/S), where `|value|` is the field
> strength *along that axis* and the sign is the pole. **This signed-bluered color is exactly the
> §4.5 magnetic-dipole coloring (the 5f carry-over)** — same field, same requirements: (1) a chosen
> projection axis (the dipole/spin axis), (2) a *real circulating B* to be non-trivial. For a static
> hedgehog `∇×n̂≈0`, so bluered colors near-zero noise — which is why orange-magnitude is the
> *default*. **Plan:** ship VIZ.2's vector-warp with **both** color options as a toggle (orange
> magnitude = honest static default; bluered signed `B·axis` = the dipole N/S view, lights up with a
> real/placeholder circulating B, axis-selectable). This folds the §4.5 N/S coloring into the WM7
> upgrade — one mesh, two color modes, validated together against the §5.3 placeholder dipole.

### 4.4 The gauge-stable charge fix (M5.6.5b carry-over) — the sign that won't flip

**Problem (Part 3 caveat):** `∇·n̂` sign flips spuriously under Evolve-PDE because the apolar
director's sign (`n̂≡−n̂`) drifts between neighbouring voxels during the slosh. The charge
*magnitude / location* is meaningful; the *local sign* is not robust.

**Fix options** (ranked simplest → most-robust):

1. `|∇·n̂|` unsigned — trivial (abs + non-diverging colormap); loses ± distinction (fine for a
   single defect).
1. defect-relative sign-fix — gauge-fix `n̂` to point outward from the defect center
   (`n̂·r̂_defect > 0`) before taking `∇·`; preserves ± for known defect positions. **The pin
   machinery already exists** — `relax_director_step` (`engine2_pde.py:748`) takes `pin_centers` +
   `pin_signs` and enforces a defect-relative director sign, but it's wired as a *seeding/relaxation*
   step (the M5.1 Frank-energy descent), NOT run per-step under Evolve-PDE. So option 2 reuses that
   logic: a lightweight per-step (or every-N-step) spatial re-pin of `director_nhat` sign by
   `n̂·r̂_defect`. Caveat: needs the live defect center(s) (fine for a seeded single/known defect;
   harder once defects move).
1. **topological winding density** (Brouwer degree) — the gauge-invariant, conserved charge that
   *never* flips. Most robust, most work: only the CPU total-`Q` sphere diagnostic
   (`compute_winding_number`) exists today; a per-voxel winding-density field would be new.

**✅ As-built decision (VIZ.1, 2026-05-30) — what we actually shipped + why:**

- **Director glyph → CENTERED + BARBLESS, unconditionally** (not a toggle). A director is an
  *apolar nematic axis* (`n̂≡−n̂`, no head/tail) — an axis *through* a point, not an arrow *from* it
  — so the **physically-correct convention** is a segment centered on the voxel with no barb. The
  old base-at-voxel + barb was a Vector(3)-`ψ` displacement holdover. Centering is a free
  gauge-fix as a bonus (`n̂→−n̂` swaps endpoints = identical). Removed the `GLYPH_CENTERED` toggle,
  the `arrow_length`/`centered` params, and the barb math. Barbs are kept only for the **polar**
  E/B field-line glyphs (VIZ.3), which genuinely have a direction.
- **`|∇·n̂|` unsigned charge mode → NOT shipped (dropped the planned WM8).** Redundant with the
  signed WM6 (same field, just sign-stripped) → UI clutter (Rodrigo). Option 1 above is therefore
  *available in principle* but deliberately not exposed.
- **Signed charge WM6 → left honest-but-flipping (option b).** The cosmetic local sign-flip stays;
  the gauge-invariant fix (option 3, **topological winding density** — the conserved charge that
  *cannot* flip) is **deferred to M5.8** (sustained dynamic runs, where reliable ± charge between
  defects becomes load-bearing). Per-frame re-pin (option 2) is the cheaper interim if ever needed
  before then, but band-aiding the eigenvector sign is not the real answer — winding density is.

**Crucial distinction (confirmed on-screen by the H-contained / F-expanding GIFs, 2026-05-30):** the
charge-region **expansion** + director **tilt** seen under Evolve-PDE is **real physics, NOT a bug**
— it is the free-defect orientation **dispersal** (M5.7.2 / M5.6.5c: `V` confines amplitude
`Tr(M²)` so the Hamiltonian H stays gathered, but `V` is rotation-invariant so it does *not* confine
director orientation → the Frank elastic F spreads). Only the **local sign-flip** is a gauge
artifact. A gauge-fix should remove the flip and leave the expansion — the expansion is the result.

**⚠️ Scope expansion → ✅ FIXED (Rodrigo 2026-05-30) — the same apolar gauge also corrupted the
director GLYPH.** Under Evolve-PDE the director glyphs *sloshed direction* (sudden 180° flips), the
**same `n̂≡−n̂` sign-flip artifact** as the charge view. Root cause: the glyph was drawn
**asymmetric** (`pos → pos + L·n̂`), so `n̂→−n̂` flipped the segment 180°. **Fixed by rendering the
director glyph CENTERED** (`pos − ½·shaft·n̂ → pos + ½·shaft·n̂`): `n̂→−n̂` now merely swaps the two
endpoints → the segment is **visually identical** → the flip is invisible. Centered + barbless is
the correct apolar convention (above), so this is unconditional, not a toggle. After centering, the
only motion left on the director glyph is the **real** tilt (EM) + free-defect orientation
**dispersal** (M5.7.2) — both genuine physics a gauge-fix must NOT remove.

| What you see on the director glyph | What it is | Removed by centering? |
| --- | --- | --- |
| sudden **180° flips** | apolar sign artifact (gauge) | ✅ yes |
| smooth **direction drift** | real **tilt** (EM/`a`-axis motion) | no — genuine |
| gradual **scrambling** over many steps | free-defect orientation **dispersal** (M5.7.2) | no — genuine |

So **two render targets, one apolar root**: (a) charge field → `|∇·n̂|` unsigned (or defect-relative
gauge-fix / winding); (b) director glyph → **centered render**. A single deeper fix — spatially
gauge-fixing `n̂` (`n̂·r̂_defect > 0`) *before* both the `∇·` and the glyph read it — would stabilize
both at once but needs the defect center; the centered-glyph + `|∇·n̂|` pair is the cheaper
independent route. **And note (from `§4.1.1`):** the director's principal axis shows only tilt
(never the clock spin) — the **δ cross-bar** (the middle-eigenvector arm of the VIZ.3 ellipsoid-cross,
Director glyph state 0) is what shows the δ-axis rotation; centering the director cleans up its view
but the spin lives on the δ arm (and in free 3D that arm tilts/disperses rather than coherently spins
— M5.7.2; coherent spin = M5.8/9b).

### 4.5 Magnetic moment + the dipole N/S coloring (M5.6.5f) — sample-first

**The honest blocker:** a *static* hedgehog is a **pure electric charge** — `∇×n̂ ≈ 0`, no
circulating B, no poles. So `B·axis` (the signed projection that would give red=N/blue=S) colors
near-zero noise. A real magnetic moment needs a **circulating B**, which only appears with a
twisting/spinning defect — i.e. the Zitterbewegung clock (M5.8) or a seeded current-loop/vortex.

**Two-stage plan:**

| Stage | What | When |
| --- | --- | --- |
| **Sample (test the render)** ✅ **DONE (VIZ.4, 2026-05-30)** | `_viz_sample_dipole` xparameter writes a hard-coded analytic dipole `B(r)=amp·[3(m̂·r̂)r̂−m̂]/max(r,r0)³` about `DIPOLE_AXIS` (=`ẑ`) into `director_curl_field` each frame (`fill_dipole_sample_B`, gated on `DIPOLE_SAMPLE`); wires the **radial** bluered N/S coloring (WM7 + `CURL_COLOR=1`) + B glyphs (state 3) + a **YELLOW moment vector glyph** (`update_moment_glyph` → `moment_glyph_*` buffers). | **pre-M5.8** (§5 placeholder strategy) |
| **Real** | point the same render at the *actual* circulating B from a twisting biaxial defect under Evolve-PDE (+ auto-axis from the net circulation); **delete ALL the placeholder scaffolding** incl the YELLOW hard-coded moment (roadmap 5f stage-2) | **M5.8** (the clock generates the real moment) |

This is the M5.6.5f carry-over: **the render path is built + smoke-tested** (Sample ✅); the *physics
source* swaps in at M5.8 with no render change.

**Axial vs radial — which scalar gives N/S poles (Rodrigo 2026-05-30):** the bluered color projects
`∇×n̂` (=B) onto a chosen direction. The choice matters:

| Projection | Formula | Dipole pattern | Bar-magnet N/S? |
| --- | --- | --- | --- |
| **Axial** (VIZ.2 original) | `B·ẑ` (fixed axis) | `B∥m̂` along the *whole* axis ⇒ RED at **both** ±ẑ ends + BLUE equatorial belt; axial:equatorial = 2:1 | ❌ reads as "2 red lobes" |
| **Radial** (VIZ.4, dipole) | `B·r̂` from center, `∝ cosθ` | RED N-hemisphere (B flows OUT), BLUE S-hemisphere (B flows IN), white equator | ✅ N-red above / S-blue below (Duda's slide) |

The dipole sample uses the **radial** projection (`_curl_signed_proj` + `curl_radial`/`curl_center`)
so it matches the bar magnet. General WM7 runs keep the axial projection (no defined center yet —
M5.8 wires radial to the real defect center). Both are real physics; radial is the one that shows
*poles*. Math confirmed (numpy + headless): radial top `+0.0093` RED/N, bottom `−0.0093` BLUE/S,
equator `0`. The biaxial hedgehog is still seeded underneath, so the `∇·n̂` charge (WM6 / E glyphs)
shows real structure for context — only the curl/B channel is the placeholder.

### 4.6 Granule positions — open question, candidate uses

`position_render` (the granule point-cloud) currently sits at `voxel + WARP_MESH·n̂` (director
point-cloud interim, M5.4 §step5). Candidate repurposes for the thermal program:

| Candidate | Shows | Cost |
| --- | --- | --- |
| **thermal heat-map** (Part 3 deferred) | color each granule by local thermal `A` → a sparse 3D heat map (vs the planar flux-mesh) | low (color by `amp_local_emarms_am`) |
| **biaxial ellipsoid** (Part 2 option c) | per-granule ellipsoid, axes `(1,δ,0)`, oriented by `O(x)` — the literal "biaxial top" | medium (new GGUI mesh; only meaningful biaxial, M5.6+) |
| **Zitterbewegung tracer** | granule traces the director's rotational orbit over time (the spinning-arrow tip) — makes the clock *visible as motion* | medium (needs a short position history) |

**Recommendation:** the **thermal heat-map** is the natural 9b use (sparse volumetric `A`, complements
the planar mesh); the **Zitterbewegung tracer** is the most pedagogically valuable for *seeing* the
clock spin (ties to `0c §L7`). Leave as open until 9b picks the priority.

---

## Part 5 — Implementation timeline + the M5.8 4×4 safety contract

### 5.1 Does M5.8 (3×3 → 4×4) break any of this? — the safety contract

**Short answer: no, with one well-isolated exception.** Audited against the live code:

| Render dependency | 3D-spatial or matrix-size-bound? | M5.8 impact |
| --- | --- | --- |
| `director_div_field`, `director_curl_field`, `‖∇×n̂‖` (`compute_director_em`) | **spatial** (∂/∂x,y,z on `director_nhat`) — independent of matrix size | **none** — grid stays 3D (`4a §6`); time is the 0-eigenvalue, NOT a 4th grid axis |
| thermal `A` = `‖M−D‖_F`, `ω` = `‖Ṁ‖_F` | Frobenius norms — work on 3×3 *or* 4×4 | **none** (D gains a 4th entry `0`; norm formula unchanged) |
| flux_mesh, warp, glyphs, granules | read `director_nhat` + scalar fields | **none** — consume *derived* fields, blind to source rank |
| **`principal_director` (Cardano eigensolver, `engine2_pde.py:258`)** | **hardcoded 3×3** (3-root Cardano) | **⚠️ THE ONE TOUCH POINT** — 4×4 `M` needs a 4×4 eigendecomposition (or: extract the spatial-3×3 sub-block for `n̂`, treat the time-row separately). This is the *only* viz-relevant kernel the promotion changes. |

**The contract for M5.8:** the 4×4 promotion must keep producing a 3-vector **`director_nhat`** and
the scalar **thermal/EM observable fields** with the *same names and meaning*. As long as
`eigen_decompose` still emits a spatial `director_nhat` (the principal *spatial* director) + the
Frobenius-deviation scalars, **every Part 4 render channel keeps working untouched**. The promotion
changes *how* `director_nhat` is computed (4×4 eigensolve, or spatial-subblock), not *what the
renderer reads*. → recorded as an M5.8 acceptance check in `0b_M5_roadmap.md`.

**Will future EM-field viz be compromised by M5.8?** No — the opposite. M5.8 *adds* the missing
ingredient: a twisting/spinning defect (the clock) produces the **real circulating B** that §4.5's
magnetic-moment/dipole viz needs. M5.8 is the *enabler* of the magnetic viz, not a threat to it.
**Evolve-PDE steps** are likewise safe: the EM/thermal observables are recomputed each frame from
whatever `M` is current; a 4×4 `M` evolving under the 4D action still yields `director_nhat` + the
deviation scalars every frame.

### 5.2 Implementation timeline — order by complexity, gated on M5.8 or not

| # | Feature | Pre-M5.8? | Complexity | Status |
| --- | --- | --- | --- | --- |
| 1 | **VIZ.1 — centered + barbless director glyph** (§4.4) — apolar-correct convention; kills the 180° slosh | ✅ yes | **low** | **✅ DONE + TESTED (2026-05-30)** |
| 1b | **Gauge-stable SIGNED charge** (§4.4) — winding-density (Brouwer) so WM6's ± can't flip | ⏳ **M5.8** | medium | deferred — see "When we address it" below; WM6 stays honest-but-flipping until then |
| 2 | **VIZ.2 — curl-vector mesh-warp + bluered N/S toggle** (§4.3) — WM7 warps by raw `∇×n̂` vector | ✅ yes | **low** | **✅ DONE (2026-05-30)** |
| 3 | **VIZ.3 — 4-state glyph select: Director / Director+Delta / E / B** (§4.2) — *Director states = ellipsoid orientation (agnostic to size/color); E/B polar + barbed; 4 mutually-exclusive checkboxes, "delta" spelled out for GGUI* | ✅ yes | **low–med** | **✅ DONE (2026-05-30)** — `director_mid` + `eigvec_for`; `_write_glyph` ti.func w/ `mode` + `show_delta`; refined per Rodrigo (added bare Director-only state so δ never shows alone; slider → 4 checkboxes) |
| 4 | **VIZ.4 — magnetic-dipole viz SAMPLE** (§4.5 stage 1) — hard-coded analytic B → bluered N/S + moment glyph | ✅ yes (placeholder) | **medium** | **✅ DONE (2026-05-30)** — `_viz_sample_dipole` xparam + `fill_dipole_sample_B` + `update_moment_glyph` + `moment_glyph_*` buffers; loader reads `GLYPH_VECTOR`/`SIZE`/`COLOR`/`CURL_COLOR`+dipole keys. Math + headless verified; **on-screen confirmed** — radial `B·r̂` gives N-red-above / S-blue-below + YELLOW moment glyph (the original axial projection's "2 red lobes" was fixed via `_curl_signed_proj`/`curl_radial`) |
| 5 | **Joint `(A·ω)²` thermal** (Part 3) — product-of-squares color mode | ⏸ 9b | low | pending |
| 6 | **Granule thermal heat-map** (§4.6) | ⏸ 9b | low | pending |
| 7 | **Magnetic-dipole viz REAL** (§4.5 stage 2) — point at the clock's actual circulating B + auto-axis | ⏳ M5.8 | low (render exists from #4) | pending |
| 8 | **Two-defect interaction** (M5.6.5e) | ⏳ M5.8 | high | pending |
| 9 | **Biaxial-ellipsoid granule** (Part 2 opt c) | ⏳ M5.6+ | medium | pending |

**The pre-M5.8 batch = #1–#4.** **#1 is done + tested.** #4 ships as a placeholder sample so the
dipole rendering is *seen and approved* before M5.8 produces the real B; then #7 re-points #4's
finished render at the live field — near-zero extra work. Items #5–#6 are 9b, and #1b/#7–#9 wait by
nature.

**When we address the signed-charge gauge (#1b):** per the plan, **only when sustained dynamic
charge-sign tracking actually matters** — that's **M5.8** (long dynamic runs) or the **two-defect /
9b** work where reliable ± charge *between* defects is load-bearing. The honest answer is **option 3,
topological winding density** (Brouwer degree, the conserved charge that physically cannot flip), and
it lands there. Until then: **(b)** — WM6 stays *honest-but-flipping*, and the charge-region
**expansion (not the flip) is the real physics story** (free-defect orientation dispersal, M5.7.2 /
M5.6.5c). Cheaper interim if ever needed before M5.8: a per-step defect-relative re-pin
(`n̂·r̂_defect > 0`, the `relax_director_step` `pin_signs` logic, currently seed-only).

### 5.3 Placeholder-sample strategy — validate rendering before the physics produces it

**Principle:** rendering and physics are decoupled. We can build + visually verify a render channel
against a **hard-coded analytic field** *now*, then swap the data source to the real computed field
when the physics lands — with no render change. This de-risks the viz work and lets Rodrigo approve
the *look* early.

| Channel | Placeholder field (hard-coded, analytic) | Swap to (real) when |
| --- | --- | --- |
| **Magnetic dipole B** (§4.5) | ideal dipole `B(r) = (μ₀/4π)[3(m̂·r̂)r̂ − m̂]/r³` about a chosen axis, written into `director_curl_field` | the twisting-defect B at M5.8 |
| **Magnetic moment μ** | a fixed vector at the defect center (e.g. `+ẑ`) | reduction over the real `∇×n̂` |
| **Joint `(A·ω)²`** | `A(r)·ω(r)` from analytic Gaussians, before the real trackers settle | the live EMA trackers |

**Implementation form:** a debug/sample xparameter (e.g. `_viz_sample_dipole.py`) or a
`VIZ_PLACEHOLDER` flag that fills the observable field with the analytic function in a one-shot
kernel, sim **paused**. The render path reads the field exactly as it will in production. **Delete
the placeholder once the real source is wired** (or keep it behind the flag as a render unit-test).
This mirrors how `_topo_*` xparameters already isolate seed configs for visual gates.

> **Repo discipline:** placeholder *fields* are physics-shaped analytic functions (fine in
> OpenWave). Keep any SABER device framing out — these samples exist to validate the *renderer*,
> which is open-source infrastructure.

### 5.4 Recommended next action

The **pre-M5.8 batch (#1–#4) is COMPLETE** (2026-05-30) — each small, independent, improving the
EM/thermal "seeing" the 9b research will lean on. **#1 (centered + barbless director glyph) ✅ done +
tested; #2 (curl-vector warp + bluered N/S toggle) ✅ done + tested; #3 (4-state glyph select:
Director / Director+Delta / E / B — E/B centered+barbed, director states centered+barbless) ✅ done +
tested; #4 (magnetic-dipole sample — analytic B + N/S + moment glyph) ✅ done + on-screen confirmed
(radial `B·r̂` → N-red-above / S-blue-below + YELLOW moment glyph).** Sprint fully closed. Defer #5–#9 per the table (the signed-charge gauge, item 1b, lands at M5.8 via winding density). Update
Part 3's deferred-table homes as each lands. **Next:** resume the curriculum (Lesson 2) with the viz
wired, then the M5.8 4D build.
