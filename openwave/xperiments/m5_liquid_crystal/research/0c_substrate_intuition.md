# M5 Substrate Intuition — the 3×3 matrix infrastructure (prep for 4×4)

**Purpose:** build working intuition for OpenWave's M5 Liquid-Crystal substrate — what the
matrix field *is*, what its parts mean physically, how a particle (defect) is built from it, how
the field evolves, where its energy/mass live, why it oscillates (the clock), how forces emerge,
and how we visualize all of it — so the M5.8 promotion to 4×4 (`5a §10b`) lands with minimal
knowledge gaps.

**Format:** built
step-by-step during a teaching session. Each lesson distills an intuition-first Q&A (math second,
always anchored to the live engine: `medium.py`, `engine1_seeds.py`, `engine2_pde.py`,
`engine3_observables.py`, `engine4_render.py`, `_launcher.py`). Lesson bodies fill in as we cover
them.

**Status legend:** ✅ done · 🔶 in progress · 🚧 next · *(blank = pending)*.

> **Curriculum refactor (2026-05-31).** Restructured 12 → 10 lessons, sprouted *by* learning the
> curriculum. The old L1–L2–L3 arch is re-split into **L1 (the grid/medium/vacuum)** + **L2 (each
> voxel's "personality" — the matrix, eigenvalues & the physics map, merging old L2+L3)**; the 4×4
> bridge moved up to **L3** (right after the content lesson, since it's an eigenvalue-spectrum
> extension); spin-½ folded into the clock at **L7**; handedness + composites is the finale at
> **L10**. Old numbering: L1→split L1/L2, L2+L3→L2, L10→L3, L4–L9 unchanged, L7+L11→L7, L12→L10.

---

## Curriculum

| # | Lesson | Covers (questions + *added topics*) |
| --- | --- | --- |
| [1](#lesson-1--the-medium-the-grid--the-vacuum) | [The medium, the grid & the vacuum](#lesson-1--the-medium-the-grid--the-vacuum) ✅ | *the medium = an LdG tensor-field `M(x)` on a 3D space grid, time-evolved; the order-parameter / coarse-graining reading; why a matrix not an arrow (the Vector(3)→matrix story); the vacuum/ground state*; "biaxial top at each voxel" |
| [2](#lesson-2--each-voxels-personality-m--odoᵀ-eigenvalues--the-physics-map) | [Each voxel's personality: `M = O·D·Oᵀ`, eigenvalues & the physics map](#lesson-2--each-voxels-personality-m--odoᵀ-eigenvalues--the-physics-map) | the 9 numbers (6 independent), `D`=eigenvalues=ellipsoid shape, `O`=eigenvectors=director frame, the director `n̂`; the eigenvalue→physics map (tilt→EM, twist→QM(ℏ), null→clock); the curvature operators `A_μ=[M,∂M]`, `F_μν=[M_μ,M_ν]` (force = curvature of the frame) + grad/div/curl/laplacian; *the M4 6-phasor-ellipse → ellipsoid bridge; natural units & δ↔ℏ* |
| [3](#lesson-3--the-4th-dimension-gravity-g--the-time-axis) | [The 4th dimension: gravity (`g`) + the time axis](#lesson-3--the-4th-dimension-gravity-g--the-time-axis) | the time axis / 0-eigenvalue, `D=diag(g,1,δ,0)`, `O∈SO(1,3)`, *teleparallelism*; gravity = time-axis scale `g`; the clock = rotation-into-time; **the two "times" (`dt` vs the matrix time index)** + the physical analogies; defers the *engine why* (negative-energy mechanism) → L7 |
| [4](#lesson-4--building-a-particle-the-biaxial-hedgehog--topology) | [Building a particle: the biaxial hedgehog & topology](#lesson-4--building-a-particle-the-biaxial-hedgehog--topology) | `O=[r̂ \| e_Θ \| e_Φ]` (the three vectors), eigenvalue melt, disclination; *+ winding number = quantized charge, Derrick's theorem → no static soliton* |
| [5](#lesson-5--energy-mass--the-ground-state) | [Energy, mass & the ground state](#lesson-5--energy-mass--the-ground-state) | *the action principle (ℒ=T−U → EOM); the energy Hamiltonian vs the Frank elastic energy; mass = stored field energy above vacuum (E=mc²); F = −∇E; the ground state* |
| [6](#lesson-6--dynamics-how-the-field-actually-moves) | [Dynamics: how the field actually moves](#lesson-6--dynamics-how-the-field-actually-moves) | *the leapfrog `evolve_M`; faithful (`4Σ‖[M_μ,Ṁ]‖²`) vs simple (`½‖Ṁ‖²`) kinetic; `V(M)` confines amplitude not orientation (the M5.7 root cause); energy conservation as the validation* |
| [7](#lesson-7--the-de-broglie-clock-engine--spin-½-zitterbewegung) | [The de Broglie clock-engine & spin-½ (Zitterbewegung)](#lesson-7--the-de-broglie-clock-engine--spin-½-zitterbewegung) | *why a topological defect can't relax → oscillates (knotted-rubber-band); the spinning-arrow visual; spinning vs oscillating; ω=2mc²/ℏ; the **engine** (Minkowski negative-energy self-propulsion — depth here); **spin-½** (SO(3) double-cover, 2ω doubling, L=ℏ/2); de Broglie λ; time-crystal* |
| [8](#lesson-8--force-emergence-coulomb-maxwell-magnetism-gravity) | [Force emergence: Coulomb, Maxwell, magnetism, gravity](#lesson-8--force-emergence-coulomb-maxwell-magnetism-gravity) | Coulomb (static topology, 1/d) ↔ Maxwell (dynamic tilts); electric (`∇·n̂`) / magnetic (`∇×n̂`) / gravitational (boosts); *EM orthogonality E⊥B in the tensor field*; magnetic moment; *magnetism as a dynamical correction to Coulomb (Feynman) vs* permanent-magnet static B with no moving charge |
| [9](#lesson-9--seeing-it-the-visualization-map) | [Seeing it: the visualization map](#lesson-9--seeing-it-the-visualization-map) | glyphs (direction=`n̂`, size, color), `flux_mesh`, `warp_mesh` scalar vs vector, granule positions, WAVE_MENU channels; *+ apolar `n̂≡−n̂` gauge sign-flip caveat* |
| [10](#lesson-10--handedness-chirality--composite-particles) | [Handedness, chirality & composite particles](#lesson-10--handedness-chirality--composite-particles) | the finale: **handedness/chirality** (traversal sign CW/CCW = ±; matter/antimatter; neutrino helicity; biaxial `π₁=Q₈` quaternion classes) + **composite particles** (9d); *seeds in L2 (ellipse handedness) + L4 (topology charge sign)* |

---

## Lesson 1 — The medium, the grid & the vacuum

> **Covers:** what the *medium* actually is — an LdG (Landau–de Gennes) symmetric-tensor field
> `M(x)` living on a 3D space grid and evolved in time; what it *represents* (order parameter /
> coarse-grained granule orbit); why M5 evolved from a Vector(3) `ψ` field to a matrix `M` (M5.2
> failed → M5.4 fixed); the **vacuum / ground state** (uniform `M=D`, no defect); the "biaxial top
> at each voxel" picture. *(The full decode of `M` — eigenvalue shape, the physics map — is L2.)*

### The one-sentence version

The M5 universe is **one field**: at every voxel of a 3D grid sits a little **oriented shape**
(a symmetric 3×3 matrix `M`), and the simulation is just those shapes pushing on their neighbors
and changing over time. Particles, charge, forces, the clock are *patterns* in this field — not
separate ingredients.

### What lives at each voxel — a "biaxial top"

Most simulators put a scalar (number) or a 3-vector (arrow) at each cell. M5 puts a tiny
**ellipsoid** — a *biaxial top*. An ellipsoid carries two independent things:

| Piece | What it is | Encoded by |
| --- | --- | --- |
| **Shape** | the 3 axis-lengths of the ellipsoid | the eigenvalues `D = diag(λ₁,λ₂,λ₃)` |
| **Orientation** | how that ellipsoid is rotated in space | the rotation `O` (its eigenvectors) |

Combined into one symmetric matrix: **`M = O·D·Oᵀ`**. The physics name is the **Landau–de Gennes
(LdG) order parameter** — the standard description of a liquid crystal. Duda's framework: *the
vacuum is a liquid-crystal-like medium; particles are defects in its alignment.* (How the
eigenvalues set the exact shape, and what each axis *means* physically, is the content of L2.)

### What the medium *represents* — the order-parameter / coarse-graining reading

![A granule tracing a fast elliptical orbit; time-averaged, its position-cloud covariance is the ellipsoid M — the M4 6-phasor ellipse to M5 ellipsoid bridge](images/granule_ellipse_small.gif)

`M` is an **order parameter**: by definition a coarse-grained average of whatever finer degrees of
freedom live below the voxel (a real liquid crystal's `M = ⟨n⊗n⟩` is averaged over many molecules).
The clean link to a *granule* picture is the **covariance**: a tiny granule tracing a high-frequency
**orbit** (the animation) has a position-cloud whose covariance `⟨(x−x̄)(x−x̄)ᵀ⟩` *is* a symmetric
matrix `O·D·Oᵀ` — eigenvalues = the orbit's variances (its **shape**), eigenvectors = its
**orientation**. So fast sub-voxel orbital motion, time-averaged at the scale we can resolve,
**presents as a static ellipsoid `M`**. This is the lineage M4 6-phasor ellipse (a granule's orbit)
→ M5 ellipsoid (its second moment) — opened in L2.

Honest scope: `M` being coarse-grained is *what it is* (and why it's tractable); a finer
granule/Planck substrate is *aligned with* OpenWave's wave-structure-of-matter roots (M1 granule
model, Yee/M3) but is an **open ontological hypothesis**, not something M5 asserts — M5 evolves `M`
as its fundamental *dynamical* field (no sub-voxel granule sim). Tractability forces this level:
voxel `~10⁻¹⁸ m` vs Planck `~1.6×10⁻³⁵ m` → resolving Planck *inside* a voxel needs `~10⁵¹`
sub-cells, impossible. The order parameter is the coarsest description that still carries the
topology / charge / clock physics.

### Why a matrix and not an arrow (the M5.2 → M5.4 story)

The field used to be a Vector(3) `ψ` (one arrow/voxel, `psi_am`, now legacy `medium.py:153`).
M5.2 tried to host the paper's physics on that arrow and **failed** (closed as an informative
negative, `medium.py:16`). An arrow can't carry:

| Limitation of one arrow | Why it matters | Matrix fix |
| --- | --- | --- |
| marks only **one** axis | three lepton families = hedgehogs of **three distinguishable axes** | `M` has 3 eigen-axes with distinct eigenvalues |
| is **polar** (head ≠ tail) | a nematic director is **apolar**: `n̂ ≡ −n̂` | `M = n̂⊗n̂` is blind to the sign of `n̂` |
| can't carry an internal **twist** | the QM clock + KG mass come from a *twist of the frame* | `O` is a full frame → it can twist |

So in M5.4 the substrate became the matrix `M` — a real `ti.Matrix.field(3,3)` triple buffer
`M_prev_am, M_am, M_new_am` (`medium.py:169-171`). Stored as the full 9 numbers, but only **6 are
independent** (`M` is symmetric, `Mᵀ=M`).

### The vacuum — the boring ground state

The vacuum is **every ellipsoid identical and aligned the same way** — `M(x)=D` everywhere,
uniform, no twist, no defect: flat, calm, lowest energy. A particle is what you get when the
alignment is knotted so the medium's elasticity can't comb it flat (a topological defect, Lesson
4). *Intuition: a 3D block of jelly with a grain to it — vacuum = perfectly straight uniform
grain; a particle = a permanent swirl the elasticity can't smooth out.*

### Q&A / clarifications (2026-05-29)

| # | Question | Answer (short) | Full in |
| --- | --- | --- | --- |
| 1 | why is a symmetric matrix an ellipsoid? | spectral theorem: `M=O·D·Oᵀ` = rotate → stretch by `D` → rotate back; applied to a unit sphere it gives an ellipsoid with **perpendicular** axes (symmetry ⟺ real eigenvalues + orthogonal axes). A non-symmetric matrix shears it into an axis-less blob. | L2 |
| 2 | `(g,1,δ)` biaxial vs `(1,δ,δ)` uniaxial? | the 3 entries are the axis-lengths: `g`≫1 (gravity/boost), `1` (EM/tilt, reference), `δ`~ℏ (QM/twist). **Biaxial** = all 3 distinct (triaxial, like a brick — 3 distinguishable axes). **Uniaxial** = 2 equal → degenerate spheroid (rugby ball — one special axis). M5.4 used uniaxial `(1,δ,δ)` because Coulomb only needs one axis; biaxial is needed for the twist/clock. | L2 (labels), L4 |
| 3 | where are `director_nhat` / `eigenvalues`? | stored at `medium.py:179-180`; **computed each frame** from `M_am` by the Cardano `eigen_decompose` in `engine2_pde.py`. Derived caches — `M` is the truth, these are read back out. | L2 |
| 4 | does the medium have an elastic force toward a ground state? | yes — orientation gradients cost energy (**Frank elastic**), giving a restoring force `F=−∇E` that relaxes toward the uniform vacuum. The drama: topology can forbid full relaxation → it oscillates. | L5 (energy), L7 (clock) |
| 5 | what is a "biaxial top"? | a rigid body by its 3 principal axes: spherical (3 equal), uniaxial/symmetric (2 equal — pencil/football), **biaxial** (all 3 different — book/phone, a full 3-axis orientation). | L1 |
| 6 | eigenvalues = shape, eigenvectors = orientation? how many? | yes: **3 eigenvalues** = shape (`D`), **3 perpendicular eigenvectors** = orientation (the columns of `O`). 6 independent numbers = **3 shape + 3 rotation angles**. | L2 |
| 7 | apolar director, but we draw a half-arrow? | the half-arrow is a **rendering aid, not physics** — the field is apolar (`n̂≡−n̂`, eigensolver sign arbitrary). The barb helps the eye follow orientation but lets neighbors disagree on sign → the gauge sign-flip artifact. | L9 |
| 8 | `director_nhat` is only the major axis — there's a 3rd orientation DoF? | **exactly.** `O∈SO(3)` has 3 rotational DoF; a unit vector `n̂` has 2 — the **missing 1 is rotation about `n̂` = the twist**. Meaningless for uniaxial (degenerate), **physical for biaxial** — and that twist DoF is where the QM phase / KG mass / the clock live. The director throws away exactly the DoF that becomes the de Broglie clock. | L2 (twist), L7 (clock) |
| 9 | what reality do the ellipsoids represent — granules at Planck scale? | `M` is an **order parameter** = coarse-grained average. A granule's high-frequency **orbit** has covariance `⟨(x−x̄)(x−x̄)ᵀ⟩ = O·D·Oᵀ` (orbit variances → eigenvalues = shape; orbit orientation → eigenvectors), so sub-voxel motion presents as a static ellipsoid. A finer granule/Planck substrate is *aligned with* OpenWave's roots but an **open hypothesis** — M5 evolves `M` directly. Tractability: voxel `~10⁻¹⁸ m` vs Planck `~1.6×10⁻³⁵ m` (`~10⁵¹` sub-cells to resolve). | L1 above, L2 (M4 ellipse), L7 |
| 10 | how do the eigenvalues map to the ellipsoid axes? | eigenvalue = axis length. Live 3D biaxial `diag(1, δ, 0)`: `1`→longest `a` (**`director n̂`**, EM axis), `δ`→middle `b` (QM, `~ℏ`), `0`→flat `c` (null → 4D clock). `O(x)` rotates the whole ellipsoid per voxel. (`medium.py:19` writes the general `diag(g,1,δ)`; `g`=gravity is the 4D addition.) See the figure in L2. | L2 (physics labels) |
| 11 | if `D` is frozen/global, why store + recompute `WaveField.eigenvalues` each frame? why `g,1,δ` not `a,b,c`? | the global `D` is the *ideal vacuum* spectrum (a constant — **not** stored per voxel). `WaveField.eigenvalues` stores the **local** eigenvalues of the actual `M(x,t)`, which **deviate** from `D`: (a) cores **melt** to isotropic `(1+δ)/3` (Faber regularization, `engine1_seeds.py:500`), (b) dynamics breathe the amplitude `Tr(M²)`. Recomputed every frame because the **director** (eigenvectors, via `O(x)`) changes every step — `eigen_decompose` returns eigenvectors *and* eigenvalues together. Symbols `g/1/δ/0` encode physics (gravity / EM-unity / QM`~ℏ` / null), not generic geometry — `a/b/c` would lose the hierarchy `g≫1≫δ~ℏ>0`. | L2 (map), L5 (V min) |
| 12 | does `δ ~ ℏ` corroborate the Planck-granule orbit hypothesis? | **intriguing rhyme, not corroboration.** Solid part: `δ` is the twist eigenvalue carrying the QM phase `exp(iEt/ℏ)` → KG mass (M5.6, verified). Speculative part: reading eigenvalues as orbit-variances (Q9) puts `ℏ` in the smallest orbital extent — but `ℏ` is an *action* and `δ` a *dimensionless ratio* (a role-identification in natural units), and the granule/Planck layer is itself unproven. Hold as a research thread, not evidence. | L7 (clock), L1/Q9 (granule) |

### Anchor in the live engine

| In `medium.py` | What it is |
| --- | --- |
| `M_am`, `M_prev_am`, `M_new_am` (`:169-171`) | the matrix field `M` at t, t−dt, t+dt — the substrate |
| `ti.Matrix.field(3,3,...)` shape `grid_size` | one 3×3 matrix per voxel on the cell-centered cubic grid |
| `D = diag(g,1,δ,0)`, `LC_DELTA=0.5` (`:19,:41`) | frozen eigenvalue spectrum; M5.4 uses uniaxial placeholder `(1,δ,δ)` |
| `director_nhat`, `eigenvalues` (`:179-180`) | *derived* each frame from `M` via `eigen_decompose` (Lesson 2) |
| `psi_am` (`:153`) | the **old** Vector(3) arrow — legacy, being retired |

**Takeaway:** the medium is a 3D grid of oriented ellipsoids (`M`); the vacuum is all of them
aligned; we use a matrix instead of an arrow because the physics needs shape + a twistable frame +
apolarity, which an arrow can't carry. The leftover "twist about `n̂`" DoF (Q8) is the seed of the
clock. *(Next, L2: decode `M` fully — what the eigenvalues are, and what moving each axis means.)*

**Anchors:** `4a §3/§5`, `medium.py`, M5.4 migration history.

---

## Lesson 2 — Each voxel's personality: `M = O·D·Oᵀ`, eigenvalues & the physics map

> **Covers:** the full decode of one voxel's matrix. Q1 (the 9 numbers, 6 independent), Q2/Q6
> (eigenvalues = shape `D`, eigenvectors = orientation `O`, the director `n̂`), Q3 (eigenvalues ↔
> matrix numbers ↔ `n̂`); the **ellipsoid axes** (how `1/δ/0` set the shape); the **eigenvalue→physics
> map** (tilt→EM, twist→QM(ℏ), null→clock); the **curvature operators** `A_μ=[M,∂_μM]`,
> `F_μν=[M_μ,M_ν]` — *a force field is a curvature (gradient) of the frame, not the frame itself* —
> plus grad/div/curl/laplacian (div=splay/charge, curl=circulation/B, laplacian=diffusion/wave);
> *the M4 6-phasor-ellipse → ellipsoid bridge (major axis / orbital normal / handedness=chirality);
> natural units & δ↔ℏ*. *(Merges old L2 + L3: the object and what its parts mean are one arc.)*

### The ellipsoid axes — how the eigenvalues set the shape

![Triaxial biaxial top: semi-axes a (longest, x) > b (medium, y) > c (shortest/flat, z); the director n-hat lies along the longest axis a](images/ellipsoid.png)

In the live **3D** substrate the biaxial vacuum spectrum is **`D = diag(1, δ, 0)`** (the M5.6
seeder, `engine1_seeds.py:477`) — three distinct axis-lengths. (`medium.py:19`'s header writes the
general `diag(g,1,δ)`; the gravity eigenvalue `g` is the **4D** addition — see L3.) The figure
above is drawn axis-aligned (`O = identity`):

| Eigenvalue | Size | Semi-axis in the figure | Physics label (the "why" → below) |
| --- | --- | --- | --- |
| `1` | largest (unity) | **`a`** (long axis, x) — **`director n̂` points here** | EM / tilt |
| `δ` | middle (`~ℏ`) | `b` (medium axis, y) | QM / twist |
| `0` | smallest (null) | `c` (short / flat axis, z) | the null/time axis → the 4D clock |

- **Where is `director_nhat`?** It's the **principal eigenvector** — the eigenvector of the
  *largest* eigenvalue (`1`, the EM axis) — so on the figure it runs **along the longest axis `a`**.
  For a hedgehog that's the radial direction (`n̂ = r̂`, the classic charge texture). The director
  captures only this *one* axis; the orientation of `b`/`c` *around* `a` is the leftover twist DoF
  from L1 Q8. (Flattening `c → 0` literally visualizes the **null** axis.)
- Bigger eigenvalue → longer axis (convention: `M` on the unit sphere ⇒ semi-axis = eigenvalue; the
  `√λ` convention keeps the same ordering). `D` is the *vacuum* shape (V(M)'s minimum); near a
  defect core it **melts** toward isotropic `(1+δ)/3`, and only `O(x)` (the orientation) varies
  freely per voxel.
- **Spectra by phase:** uniaxial M5.4 placeholder `diag(1, δ, δ)` (`δ=LC_DELTA=0.5` — one director
  axis, enough for Coulomb) → biaxial M5.6 `diag(1, δ, 0)` → full 4D `diag(g, 1, δ, 0)` adds `g` =
  gravity/boost (L3). Hierarchy `g ≫ 1 ≫ δ ~ ℏ > 0`.

### The M4 ellipse → M5 ellipsoid bridge (+ chirality, δ↔ℏ)

> **Seed intuition (to develop).** `M = O·D·Oᵀ` is literally an *ellipsoid at each voxel* — `D =
> diag(λ₁,λ₂,λ₃)` are the semi-axis lengths (the **shape**), `O` is orthogonal (the
> **orientation/rotation**). This is the 3D matrix generalization of M4's **6-phasor ellipse** (`R,
> Φ` per axis). An ellipse/ellipsoid carries orientation *for free*: the **major axis** = one
> direction in space, the **normal to the orbital plane** = another, and the **handedness** (CW vs
> CCW traversal) = a ± sign = **chirality** (the seed for L10). So one symmetric matrix encodes
> direction + shape + chirality together. *Natural units:* `δ ↔ ℏ` is a role-identification (the
> QM/twist eigenvalue plays the action quantum in the dimensionless scaling), not a dimensional
> equality — see L1 Q12.

### The eigenvalue→physics map + the curvature operators

> **Covers (from the merged old L3):** Q4 (how force fields are encoded) — the eigenvalue→physics
> map + *the yaw/pitch/roll framing of the three rotation generators; the curvature operators
> `A_μ=[M,∂_μM]`, `F_μν=[M_μ,M_ν]` — a force field is a **curvature** (gradient) of the frame, not
> the frame itself; the vector operators (grad / divergence / curl / laplacian) and what each means
> physically (div = splay/charge, curl = circulation/B, laplacian = diffusion/wave)*.
>
> The eigenvalue→physics map. Each axis = a kind of local orientation change: tilt→EM, twist→QM(ℏ),
> boost→gravity (L3, 4D), null→clock. The key idea: a force field is a curvature (gradient) of the
> frame, not the frame itself.
>
> **Seed (from L1 Q&A — twist vs tilt on the director):** the generators act on the *same* director
> `n̂` (the principal axis, eigenvalue `1`) in orthogonal ways → different physics. **Twist** = rotate
> the δ-axis (`director_mid`, λ=`δ`) and null axis (λ=`0`) *about* `n̂` (`n̂` stays put; generator
> `Gx`) → **QM / the clock** (L7). **Tilt** = rotate `n̂` *itself* toward another axis (its direction
> changes; generators `Gy,Gz`) → **EM**, whose *field* is the spatial gradient of those tilts
> (`∇·n̂`=charge, `∇×n̂`=B — L8). **Boost** (4D) → gravity (L3); the **null** (`0`) axis → the clock
> direction. Mnemonic: *twist about it = QM, tilt of it = EM.*

(to be filled during the session)

**Anchors:** `medium.py` (M storage), `engine2_pde.py` (Cardano eigensolver + operators), `4a §5/§6/§8`, `5a §1-2` (Eq.18-20), M4 6-phasor model, `1b`.

---

## Lesson 3 — The 4th dimension: gravity (`g`) + the time axis

> **Covers:** what the M5.8 promotion to 4×4 adds — the time axis / 0-eigenvalue, `D=diag(g,1,δ,0)`,
> `O∈SO(1,3)`, *teleparallelism* (the 4D liquid-crystal extension); **gravity = the time-axis scale
> `g`**; **the clock = a rotation into the time axis**; **the two "times"** (`dt` vs the matrix time
> index) + the physical analogies. Placed right after the content lesson because it's an
> eigenvalue-spectrum + frame-group extension of L2. The *engine why* (why the clock self-sustains —
> the Minkowski negative-energy mechanism) is developed at **L7**; here we set up the structure.
>
> 🔭 **What the 4th dimension adds to the 3×3 medium.** Everything in L1–L2 — the biaxial frame, the
> EM/tilt axis, the QM/twist axis, charge, the null axis — already lives in the **3×3 (spatial)**
> matrix. The 4×4 promotion doesn't replace any of it; it adds **one thing — a time axis to the
> order parameter** — and because *both* gravity and the clock are time-sector phenomena, both fall
> out of that single addition:
>
> | What 4D adds (one thing: the time axis) | Mechanism | Why it matters |
> | --- | --- | --- |
> | **Gravity** = the time axis's **scale** `g` | spectrum `diag(1, δ, 0)` → `diag(g, 1, δ, 0)`; `g` scales the *time* axis, and `∇g` (how it varies in space) = the gravitational field | mirrors GR (`Φ` lives in the time-time metric `g₀₀`; gravity = time-rate varying in space). `g` is **absent in the live 3D run** — no time slot → no gravity. (Render spec: L8 / `4b §4.7`.) ⚠️ framework's least-developed sector — design expectation, not yet verified. |
> | **The clock / engine** = a **rotation into** the time axis | the frame `O` goes `SO(3)` (spatial rotations) → `SO(1,3)` (Lorentz: rotations **and boosts**); the `(−+++)` sign flip makes the oscillating state **lower-energy than static** (M5.8.0a: `E(ω=0)=2.87 > E*=2.02`) → spinning is the ground state | the Zitterbewegung clock **propels itself** — fuel = the rest mass (`ℏω = mc²`). In free 3D (all-`+`) spinning only *costs* energy → it disperses (M5.7.2); **4D = the stable time-crystal** (`5a §10b`, L7). The *why* (negative-energy mechanism) → L7. |
>
> **The two "times" — `dt` vs the matrix time axis (Rodrigo's Q, 2026-05-31; the conceptual key).**
> There is **one** physical time playing **two roles**, and the 3×3 sim uses only one:
>
> | Role | Name | What it is | In 3×3? |
> | --- | --- | --- | --- |
> | **The projector** | `dt` | the parameter you *advance* to get the next frame (`M(t+dt)=M(t)+…`) — how any movie is computed | ✅ always |
> | **A direction to lean into** | the matrix's 4th slot | a spacetime *direction* the field's orientation `O` can point along — like x/y/z, but time | ❌ only at 4×4 |
>
> *Physical analogies (keep building these — this is the unlock for SABER engineering on the clock):*
>
> - **Clock-hand analogy (ties to L7's own language).** In L7 the δ-axis is the "clock-hand"
>   sweeping around the director `n̂`. The question is *what plane does the hand sweep in?* In **3×3**
>   the hand can only sweep a **spatial** plane (δ–0) — `dt` advances but the hand has no *time*
>   direction to wind into, so it's a spatial proxy that **disperses** (M5.7.2). In **4×4** the hand
>   sweeps a plane that **includes the time axis** — a space-time rotation — and *that* winding is
>   the genuine de Broglie phase `e^{-iEt/ℏ}`.
> - **Wristwatch-per-voxel (Rodrigo's image — correct):** each voxel carries its own clock phase,
>   and the phase *differing voxel-to-voxel* is literally what **momentum / de Broglie wavelength**
>   is. ✅ keeper. **Caveat (correct the "time-only" reading):** the clock is **not** a purely-time
>   rotation and **not** purely-spatial — it's a rotation in a plane that **mixes** a space axis with
>   the time axis. `n̂` stays the fixed axle; what winds is the frame in the δ-axis↔time-axis plane.
> - **The precise version:** a phase `e^{-iEt/ℏ}` needs *two* things — a `t` that **advances** (that's
>   `dt`) **and** a time direction to **wind into** (that's the 4th matrix axis). The 3×3 sim has the
>   first but not the second. So: **`dt` = time *passing*; the matrix time-axis = the particle
>   *having a hand that points into time*, so it can keep its own clock.** The clock is the hand
>   winding, not the passing. A 3×3 particle just sits while `dt` passes; a 4×4 particle carries a
>   hand that winds into time → the real clock.
> - 🚧 *more physical analogies to develop here at teach time* — Rodrigo flagged this as the key
>   intuition for enabling SABER engineering on the clock/time/gravity lever (`SABER 0_OVERVIEW §4`).
>   Take the time it needs.
>
> **Why gravity + clock arrive together.** Gravity = the time axis's *scale* (`g`); the clock = a
> *rotation into* the time axis. Both need the time axis to exist → both arrive at 4×4, neither
> sooner. That single structural addition is why the SABER unification (thermal/time/gravity as one
> lever) is even conceivable — `SABER 0_OVERVIEW §4`.
>
> **Clock vs engine — it's both.** "Clock" = the *measurement* (it ticks at a fixed `ω=2mc²/ℏ`, the
> de Broglie clock Catillon measured). "Engine" = the *mechanism* (a self-propelled rotation; Duda's
> "oscillation propelled by mass"). A **self-propelling clock** whose output is a precise frequency =
> the particle's mass. Not perpetual motion — the spin *is* the rest energy. (The energetic *why* is
> L7.)
>
> One-liner: **3×3 = the spatial particle (shape, EM, QM-twist, charge); the 4th dimension adds the
> time axis — and from it, gravity (`g` = time-scale) and the self-propelling clock-engine (a
> rotation into time the Minkowski sign keeps running).**

(to be filled during the session)

**Anchors:** `5a §10b`, `4a §6`, `theory/time_crystal.pdf`.

---

## Lesson 4 — Building a particle: the biaxial hedgehog & topology

> **Covers:** Q6 (the three vectors) — how `O(x)=[r̂ | e_Θ | e_Φ]` is laid out in space, the
> radial eigenvalue melt, the disclination line; *+ topological winding number = quantized charge,
> Derrick's theorem → why no stable static soliton exists (sets up the clock)*.
>
> Building a particle: the biaxial hedgehog. How O(x)=[r | e_0 | e_$] (the three vectors) is laid out in space, the eigenvalue melt + disclination, and why winding = quantized charge
>
> **★ Cover when we teach L4 (Rodrigo flagged 2026-05-30): "where is the charge sign?"** The
> *uniaxial* seed defines charge per-defect via `DEFECTS:[{"SIGN": ±1}]` → kernel `n̂ = SIGN·r̂`
> (outward `+1` / inward `−1`; `_topo_uniaxial2.py`, `engine1_seeds.py:439`). The *biaxial* seed
> (`_topo_biaxial1_v{on,off}.py`) has **NO sign knob** — `seed_biaxial_hedgehog_M` is single-center
> and hard-builds the radial frame `O=[r̂|e_Θ|e_Φ]` (one fixed winding). **Why:** a biaxial defect's
> charge is **not** a single `±1` — its order-parameter space is `SO(3)/D₂`, whose `π₁` is the
> **non-abelian quaternion group `Q₈`**, a richer classification than one sign bit. So "the charge
> sign" of a biaxial defect is a *quaternion-class label*, not a `±`. This is the discovery hook of
> the deferred two-defect demo (M5.6.5e → M5.8), and the seed for L10 (handedness). Teach: uniaxial
> `±` winding → biaxial `Q₈` classes.

(to be filled during the session)

**Anchors:** `engine1_seeds.py` (`seed_biaxial_hedgehog_M`), `5a §5b/§5e`, `1b`.

---

## Lesson 5 — Energy, mass & the ground state

> **Covers:** *the action principle (`ℒ = T − U`, least action → the Euler–Lagrange EOM); the
> energy **Hamiltonian** (the full conserved energy `Σ‖F_μν‖² + V`) vs the **Frank elastic
> energy** (the director-distortion piece, the classic LC energy); **mass = stored field energy
> above the vacuum** (`E = mc²`, the M5 `E ∝ K` lepton-mass result); **F = −∇E** (force is the
> gradient of energy); the ground state and why a defect is pinned above it*.

(to be filled during the session)

**Anchors:** `5a §1` (action) / `§6` (Hamiltonian), `1b` (E∝K mass), `5a §5c` (Faber mass scale), `3a` (F from E).

---

## Lesson 6 — Dynamics: how the field actually moves

> **Covers:** *the leapfrog time-stepper (`evolve_M`); the kinetic metric — faithful
> `4Σ‖[M_μ,Ṁ]‖²` vs the shipped simple `½‖Ṁ‖²`, the degeneracy, why the twist is dynamical only
> on a non-uniform (hedgehog) background; `V(M)` — confines amplitude `Tr(M²)` but NOT orientation
> (the root cause of the M5.7 free-dispersal nulls); energy conservation as the correctness test*.

(to be filled during the session)

**Anchors:** `engine2_pde.py`, `5a §5f/§5g/§9`.

---

## Lesson 7 — The de Broglie clock-engine & spin-½ (Zitterbewegung)

> **Covers:** *where the time-crystal / Zitterbewegung enters; how oscillation can be "propelled by
> mass"; whether the clock is a **spin** (ω only) or an **oscillation** (A & ω); the rotational
> axis (yaw/pitch/roll); `ω = 2mc²/ℏ`; **the engine — the Minkowski negative-energy self-propulsion
> mechanism (depth here)**; **spin-½** (SO(3) double-cover, the `2ω` doubling, `L=ℏ/2`); the de
> Broglie wavelength λ; the bridge to 4D / teleparallelism (structure set up in L3)*.
> *(Merges old L7 + the old spin-½ deep-dive: spin-½ is a property of this clock.)*
>
> **Seed intuition (to develop).**
>
> *Why it can't sit still — the knotted rubber band.* A rubber band stretched between two posts
> relaxes flat (a static ground state). Tie a topological **knot** in the middle: the tension still
> wants to relax, but topology forbids untying the knot. The band can neither relax fully nor sit
> statically at the knotted-stretched configuration — so the knot **vibrates** at a frequency set
> by the local elastic restoring force. The oscillation is the compromise between the *topological*
> constraint (cannot unwind) and the *energetic* constraint (wants minimum elastic energy): they
> can't both be satisfied in a static configuration, so the next-lowest-energy state is *moving*.
>
> *What to picture — a spinning arrow, not a bouncing ball.* Don't picture the defect as a point
> bouncing on a spring (translation). Picture a point with a **spinning arrow stuck through it**:
> the arrow is the local director orientation, and it **rotates** about an axis at `ω = 2mc²/ℏ`.
> The defect's *position* is fixed (or slowly drifting under external forces); the field's
> *orientation* at and around the defect rotates at the Zitterbewegung frequency.
>
> *The one-sentence version.* A topological defect in a Lagrangian field with the right potential
> (Duda's φ⁴ + curvature coupling, or M5's full LdG + Skyrme + KG) is **permanently displaced from
> the vacuum minimum but cannot fully relax due to topology**, so it oscillates around its
> constrained position at `ω = 2mc²/ℏ` — and the oscillation is **rotational** (the director winds
> around an axis), not translational.
>
> *Open questions to resolve in this lesson:* is the clock a steady **spin** (ω only, like a wheel / a spinning top
> / the Earth) or an **oscillation** (A & ω, like a pendulum)? And where do the **magnetic moment**,
> **spin-½**, and the **de Broglie λ** live in this picture?
>
> *Preliminary resolution (2026-05-30, from L1 Q&A):* it is a **SPIN**, not a pendulum — `ψ = ωt`,
> the phase angle winds *linearly* (toy model + Dirac `exp(iEt/ℏ)`). Why a spin and not a swing: the
> phase has **no preferred angle** (every phase has identical energy — the unconfined orientation
> direction, the same fact as "V confines amplitude not orientation", L6), so it rotates freely
> rather than being pulled back to a center. The rate isn't free, though — the time-crystal `−αω²`
> term energetically *selects* `ω = 2mc²/ℏ`. **Axis:** the **twist about the director `n̂`**
> (the principal axis, eigenvalue `1` = EM axis; generator `Gx`, a rotation in the `δ`–`0` plane,
> `5a §7a`) — `n̂` stays put while the δ-axis (`director_mid`, λ=`δ`) and null axis (λ=`0`) sweep
> *around* it (the "spinning arrow", and the leftover DoF from L1 Q8). Contrast with **EM = tilting
> `n̂` itself** (L1 Q2 / L8). Caveat: the clean
> steady spin is the *ideal/target* — in 3D the free defect disperses (M5.7.2), so confirming it at
> `ω=2mc²/ℏ` needs 4D (M5.8, structure in L3).
>
> *Refinement — "spin" does NOT mean "no amplitude" (2026-05-30):* a spin is **constant-amplitude**,
> not amplitude-free — like a point on a wheel, `(R·cos ωt, R·sin ωt)` has both a radius `R` and a
> rate `ω`. The clock's amplitude (radius/modulus, set by the rest mass) is *fixed* at the ground
> state, not absent; what makes it a spin vs a pendulum is the absence of a restoring force, not the
> absence of `A`. So the joint **`(A, ω)`** picture holds: `ω` = the rate (launcher WM3 "Thermal
> Clock"), `A` = the rotational radius/excitation amplitude (WM2 "Thermal Amp"). Ground state: both
> fixed. *Excitation/heat* modulates **both** — radius grows (AM) and/or rate shifts (FM) — exactly
> the dual channels M5.7.3 saw respond to driving (`A_core`→3×, director at `f_d`). (The
> heat-as-`(A,ω)`-excess *interpretation* is SABER/DHC, kept in the private repo per the cardinal
> rule; here it's just the substrate clock's two observables.)
>
> *What the "radius" is, and what spins (2026-05-30):* as the frame twists about the director `n̂`,
> the δ–`0` block `R(ψ)·diag(δ,0)·R(ψ)ᵀ` makes the field *components* oscillate — off-diagonal
> `(δ/2)sin(2ωt)`. So the **radius = the eigenvalue gap `(δ − 0)/2 = δ/2`** (the QM-`δ` axis sets
> its size), and the observable cycles at **`2ω`** — the apolar `n̂⊗n̂` doubling (a 180° turn looks
> identical), the likely origin of `ω_Zitt = 2mc²/ℏ` (confirm in M5.8). And *what* spins is the
> **orientation field as one collective phase-locked mode** localized on the defect (the director
> `n̂` = the fixed axle, the δ-axis `director_mid` = the clock-hand; all hands synchronized) — **not**
> the defect-as-a-point and **not** independent voxels.
>
> **★ "Collective mode" — flagged 2026-05-30 to unpack when we teach L7 (Rodrigo: idea not yet
> fully internalized).**
>
> *The mental image to hold:* a 3D field of **compass needles all precessing in lock-step**. Each
> needle (ellipsoid) turns, but the meaningful object is the **collective phase**. The director `n̂`
> at each point is the fixed **axle**; the δ-axis (`director_mid`) is the **clock-hand** sweeping
> around it — and all the clock-hands across the defect's neighborhood are **synchronized**.
>
> *Why it matters:* the **defect** is the stable topological knot that *hosts* the clock; the
> **clock** is the collective twist mode the knot carries. So "the particle" = defect (topology,
> permanent) **+** its intrinsic clock (the collective oscillation). And thermal `(A, ω)` is a
> **per-defect collective property**, not an individual-voxel one. *(3D caveat: the free collective
> mode disperses — M5.7.2 — so it needs 4D to self-sustain, M5.8, or a drive, 9b.)*
>
> *Questions to answer in L7 so this clicks:*
>
> 1. **What makes the needles lock-step?** — the elastic (Frank) coupling between neighbours + the
>    hedgehog's own `C_μν` source: the lowest-energy way to carry a twist is *coherently*, so one
>    collective phase `ψ(t)` wins over each voxel doing its own thing.
> 2. **How can a single `ψ(t)` describe a whole 3D region?** — it's a **collective coordinate /
>    normal mode**: like a drumhead's fundamental mode, *one* amplitude sets the whole membrane's
>    shape. `ψ` is that one number for the defect's twist mode.
> 3. **Why a "mode" and not a rigid point-spin?** — it has **spatial extent + a profile** (peaked at
>    the core, fading out, weighted by the texture) — a standing-wave/breather-like object, not a
>    structureless spinning dot.
> 4. **Precise vs idealized:** "lock-step" is the *coherent ideal*; really `ψ = ψ(x, t)` is a field
>    with a spatial profile. *How coherent/localized it stays* is exactly the M5.7.2 dispersal
>    result (free 3D radiates it away) → M5.8 (4D stabilizes) / 9b (drive sustains).
> 5. **Engine tie-in:** this is the `O(x) ∈ SO(3)` rotation DoF from `5a §9` / `m5_6_2b` — the
>    collective twist is a coherent excitation of that rotation field, not of `M`'s raw components.
>
> ---
>
> **★ THE ENGINE — depth (Rodrigo 2026-05-31): why the clock self-sustains.** The clock is not just
> a measurement — it is an **engine** that propels its own rotation. Develop here:
>
> - **The negative-energy mechanism (the fuel).** In 4D the kinetic term picks up the Minkowski
>   `(−+++)` sign: a rotation that leans into the time axis contributes *negatively* to the energy.
>   So the oscillating state is **lower-energy than the static one** — verified in the 1+1D toy model
>   (M5.8.0a: `E(ω=0)=2.87 > E*=2.02 = E(ω*=1.29)`). Spinning is the ground state; it costs nothing
>   to keep running because stopping would cost *more*. Fuel = the rest mass itself (`ℏω = mc²`).
> - **Why 3D can't (the contrast).** In the space-only 3×3 (all-`+` signature) a rotation only *adds*
>   energy → the free defect sheds it and disperses (M5.7.2). The engine needs the time axis's
>   negative-signature term, i.e. 4D (L3 sets up the structure). A *driven* 3D defect (9b) borrows
>   the energy externally — a motor with a power cord, not yet a self-contained engine.
> - **Not perpetual motion.** The spin *is* the rest energy; you can't extract it without destroying
>   the particle. What CAN be modulated is the *excess* above ground state (`(A,ω)` — the thermal
>   channel) and, speculatively, the rate/scale (time/gravity) — the SABER unification thesis
>   (`SABER 0_OVERVIEW §4`), an engineering bet on this mechanism, not a claim it's free energy.
>
> **★ SPIN-½ — folded in here (was the old L11 deep-dive).** Spin-½ is a property *of this clock*, so
> it lives in this lesson. Three threads:
>
> | Thread | One-line | Anchor |
> | --- | --- | --- |
> | **The `2ω` apolar doubling** | the clock's observable cycles at `2ω` because the order parameter is apolar (`n̂⊗n̂`, a 180° turn looks identical) → the origin of `ω_Zitt = 2mc²/ℏ` | the "radius/what-spins" note above |
> | **The SO(3) double-cover** | a 2π rotation does **not** restore the state (you need 4π) — the topological signature of spin-½; `O∈SO(3)` lifts to `SU(2)` | `5a §10` |
> | **Spin = the clock's angular momentum** | `L = ℏ/2` is the conserved angular momentum of the self-propelled rotation (the engine's "flywheel") | M5.8 |
>
> Full concreteness needs the 4D clock (M5.8, structure in L3). Tie-in: Duda's slide gyroscope/`L=ℏ/2`
> inset (the L8 magnetic-moment figure).

(to be filled during the session)

**Anchors:** `5a §10` (toy model), `theory/time_crystal.pdf`, `1b` (Derrick/time-crystal), `4a §6`.

---

## Lesson 8 — Force emergence: Coulomb, Maxwell, magnetism, gravity

> **Covers:** Q5 (Coulomb↔Maxwell, electric/magnetic/gravitational emergence), Q7 (magnetic
> moment — where/how to view), Q8 (permanent magnet static field with no moving charge) — static
> topology→Coulomb 1/d; dynamic tilts→Maxwell (both routes); electric=`∇·n̂` splay,
> magnetic=`∇×n̂` curl, gravitational=boosts; *EM orthogonality E⊥B in the tensor field*;
> *magnetism as a dynamical (relativistic) correction to Coulomb between moving charges (Feynman
> framing) vs the permanent magnet's static B from aligned spin-topology (no moving charge needed)*.

(to be filled during the session)

> **★ MAGNETIC MOMENT — dedicated unpack (Rodrigo flagged 2026-05-30: "I still don't fully grasp
> it").** Teach this against **Duda's electron slide** (`Screenshot 2026-05-28 at 3.15.48 PM` — the
> one with the bar magnet `m`, the `B` field lines, the spin axis + `ω` Larmor inset). That slide is
> literally the L8 target picture: **electric charge** `E ∝ 1/r²` (left), **magnetic dipole** `B ∝
> 1/r³` (center bar magnet), and the **gyroscope / spin** `L = ℏ/2` ticking at the Zitterbewegung
> `ω ≈ 2mc²/ℏ` (right). All three are *the same defect* seen through three observables.
>
> *Questions to answer in L8 so the moment clicks (relate each to what we built in the VIZ.4 session
> 2026-05-30):*
>
> | Question | Tie to what we built / where it lives |
> | --- | --- |
> | **What IS the magnetic moment `m`?** A vector: the axis + strength of a current loop / spinning charge. For our defect it is the **clock's spin axis** (the δ-twist rotation axis, L7) — *not* an independent thing. Spin ⇒ moment. | The **YELLOW moment glyph** (`update_moment_glyph`) is a literal arrow of `m̂`. In VIZ.4 it's a hard-coded `+ẑ` placeholder; at M5.8 it becomes `m̂ ∝ ∫∇×n̂` (the real net circulation) — roadmap 5f stage-2. |
> | **Why does a static hedgehog have NO moment?** It's a pure electric charge: `∇·n̂≠0` (splay) but `∇×n̂≈0` (no circulation) ⇒ no `B`, no poles. A moment needs a *circulating* `B`, which needs a *twisting/spinning* defect (the clock). | This is exactly why VIZ.4 needed a **placeholder** dipole — the static seed produces no real `B` to color yet (`4b §4.5`). The real `B` appears only at M5.8. |
> | **Where do the N/S poles come from?** `B = ∇×n̂` (a vector field). Color it by `B·r̂` (radial, from the defect center): red where `B` flows OUT (N hemisphere), blue where it flows IN (S) → `∝ cosθ` = the bar-magnet picture. *Axial* `B·ẑ` instead lights both ends red (the field's axial component) — real, but not "poles". | The **axial-vs-radial fix** we made this session (`_curl_signed_proj`, `curl_radial`). Rodrigo's "2 red spheres" observation IS the axial projection; radial gives Duda's N-red-above / S-blue-below. |
> | **Moment vs spin vs charge — how do they differ?** Charge = `∇·n̂` (scalar, monopole, `1/r²` field). Moment = `∇×n̂` integrated (vector, dipole, `1/r³` field). Spin = the *mechanical* rotation generating the moment (`L=ℏ/2`). The moment is the *magnetic shadow* of the spin. | Three WAVE_MENU/glyph channels: WM6 / E glyphs (charge), WM7 / B glyphs + moment glyph (moment), WM2/WM3 thermal A/ω (the spin rate, L7). |
> | **Why does B fall off so much faster than E?** (see the falloff headsup below — careful: field ≠ force ≠ our observable) | The `1/r³` (B) vs `1/r` (our `∇·n̂` E observable) gap is exactly why the B viz needed a different colormap calibration (γ-spread). |
> | **Permanent magnet (Q8) — static B, no moving charge?** Aligned spin-topology: many defects with their moments `m̂` locked parallel ⇒ macroscopic static `B`. No *translating* charge needed — the "current" is the frozen collective spin (the L7 collective mode). | Forward link to the L7 collective mode + the M5.8 multi-defect work (5e). |
>
> **★ FALLOFF HEADSUP (Rodrigo 2026-05-31) — don't conflate field / force / our observable.** The
> question "B ∝ 1/r³ vs E ∝ 1/r?" mixes three different things. Keep them separate:
>
> | | Electric | Magnetic | Gravitational |
> | --- | --- | --- | --- |
> | **What we RENDER** (glyph/mesh observable) | `∇·n̂` splay — hedgehog `n̂=r̂` ⇒ `2/r` → **1/r** | placeholder dipole `B` → **1/r³** | *not yet* — the boost-`g` field, M5.8 4D (`4b §4.7`) |
> | **Real FIELD of a point source** | charge (monopole): `E ∝ 1/r²` | dipole (no monopole): `B ∝ 1/r³` | mass (monopole): `g ∝ 1/r²` |
> | **Real FORCE law** | Coulomb: `F ∝ 1/r²` | dipole–dipole: `F ∝ 1/r⁴` | Newton: `F ∝ 1/r²` |
>
> So: ✅ our **B observable** is `1/r³` and our **E observable** is `1/r` (this IS why the B viz
> collapsed to black under a linear map — 9× steeper). ❌ but the real **electric field/force is
> `1/r²`** (Coulomb), NOT `1/r`; the `1/r` is specifically our `∇·n̂` *splay* observable, which
> tracks the Coulomb **potential** (∝1/r), not the field. In M5 the actual Coulomb behavior showed
> up as the **interaction energy** `E(d) ∝ 1/d` between two defects (M5.1/M5.4), so the force
> `−dE/dr ∝ 1/r²`. **Why magnetic *starts* steeper:** nature has electric monopoles (charges) →
> the E series starts at `1/r²`; it has **no magnetic monopole** → the B series starts one multipole
> higher, at the **dipole** = `1/r³`. Compare like-for-like and the asymmetry vanishes: an electric
> *dipole* field also falls as `1/r³`. (Quadrupole `1/r⁴`, etc.)
>
> **Gravity is the third case — and it's a monopole, like electricity.** Mass is a gravitational
> "charge" with only ONE sign (always positive ⇒ always attractive), so the gravitational field is a
> **monopole**: `g ∝ 1/r²`, force `F ∝ 1/r²` — the SAME falloff as the electric charge. So the family
> is: **monopole fields (electric charge, mass) → `1/r²`; dipole field (magnetism, no monopole) →
> `1/r³`.** When OpenWave renders gravity (M5.8 4D, the boost-`g` axis; viz spec in `4b §4.7`), expect
> it to spread on screen like E (gentle `1/r²`), not like the steep `1/r³` B that needed γ-compression
> — and use a single-sign sequential palette (no ± / bluered), since there's no "negative mass".
>
> Optional hands-on when teaching: launch **`_viz_sample_dipole`**, toggle WM7 bluered (radial N/S),
> flip to the Magnetic-Field glyphs (state 3, field lines), and watch the YELLOW `m̂` — then state
> plainly that this is a *placeholder shape*; the real moment is *generated by* the clock's spin at
> M5.8. The point of the session was to make the *picture* legible before the physics produces it.

**Anchors:** `engine3_observables.py` (`compute_director_em`), `5a §5d`, `3a`.

---

## Lesson 9 — Seeing it: the visualization map

> **Covers:** Q9 — how glyphs (direction=`n̂`, size=magnitude, color=observable), `flux_mesh`
> coloring, `warp_mesh` (scalar vs vector), and granule positions each render a piece of the
> physics; what every WAVE_MENU channel shows; *+ the apolar director `n̂≡−n̂` gauge sign-flip
> caveat*.

(to be filled during the session)

**Anchors:** `engine4_render.py`, `4b Part 3`, `_launcher.py`.

---

## Lesson 10 — Handedness, chirality & composite particles

> **The finale (M5.8 / M5.9-era).** Two intertwined threads: **handedness/chirality** (the ± that
> distinguishes matter from antimatter, and helicity for neutrinos) and **composite particles** (9d —
> how confirmed single defects combine). Both become load-bearing only once the clock (L7), the 4D
> structure (L3), and the lepton families (M5.9) are in place — hence last.
>
> 🚧 **Future finale (added 2026-05-31 from the old L12; not taught yet — slotted so it isn't
> lost).** **Seeds already exist in L2** (the M4-ellipse → ellipsoid bridge: traversal sign =
> chirality) **and L4** (topological charge sign / the biaxial `Q₈` classes).
>
> **Where handedness comes in (the threads to develop):**
>
> | Thread | One-line | Anchor |
> | --- | --- | --- |
> | **Orbit traversal sign** | the ellipse/ellipsoid is traversed CW or CCW — a ± that one symmetric matrix carries "for free" alongside direction + shape = **chirality** | L2 seed |
> | **Matter ↔ antimatter** | the charge-sign / winding-sign flip (the `±` hedgehog) — candidate for the matter/antimatter distinction | L4 (winding), `4b §4.4` |
> | **Neutrino helicity** | left/right-handed states; closed-vortex-loop candidates carry an intrinsic handedness | `1b`, M5.9 frontier |
> | **Biaxial subtlety** | a biaxial defect's "sign" is **not** a simple `±` — it's a quaternion class (`π₁(SO(3)/D₂)=Q₈`); handedness there is richer | L4, roadmap 5e |
>
> **Composite particles (9d).** Once single defects (leptons M5.9, quark-vortices) are confirmed,
> how they bind into composites (baryons, mesons) — the Cornell-potential quark string (M5.9) is the
> first composite hook. Deferred program: roadmap **9d**. Chirality is load-bearing here (a proton's
> quark content has definite handedness).
>
> Prereqs: L2 (the ellipsoid encodes chirality) + L4 (topology/winding sign) + L7 (clock) + M5.9
> (lepton families). Tie-in: the M5.6.5e two-defect demo (what "opposite handedness/charge" even
> means for a biaxial defect).

**Anchors:** `4a §5` (ellipse handedness), L2/L4 seeds, `1b`, 9d (composites).

---

## Appendix A — source questions (2026-05-29, voice-note batch 1)

The original questions this curriculum organizes:

1. What are the numerical representations in the 3×3/4×4 matrix? What do they represent physically?
1. What is the relationship between the eigenvalues, the matrix numbers, and the director vector?
1. How / what encodes force fields (EM, gravity, etc.)?
1. How do Coulomb and Maxwell relate? Where do electric / magnetic / gravitational forces emerge?
1. How is the biaxial hedgehog defined from the three vectors?
1. Where is / how to view the magnetic moment?
1. How do permanent magnets hold a permanent magnetic field with no moving charge?
1. How does all this translate to the M5 visualization (glyphs: direction, size, color + `flux_mesh` + granule positions)?

## Appendix B — added concepts (2026-05-29, batch 2) → lesson map

> Lesson numbers updated for the 2026-05-31 refactor (12 → 10 lessons).

| Added concept | Lands in |
| --- | --- |
| The medium (LdG tensor-field on a 3D grid, time-evolved); the vacuum state | L1 |
| The action principle | L5 |
| Particle mass / stored energy / ground state; Hamiltonian vs Frank elastic; F=−∇E | L5 |
| Time-crystal & Zitterbewegung; how oscillation is propelled by mass | L7 |
| Oscillation axes — yaw / pitch / roll | L2 (axes) + L7 (which axis is the clock) |
| charge/winding, spin, magnetic moment, de Broglie clock | L4 (winding) + L7 (spin, clock) + L8 (moment) |
| Vector operators: gradient, divergence, curl, laplacian | L2 |
| EM orthogonality (E⊥B) in the tensor field | L8 |
| Magnetism as dynamical correction to Coulomb (Feynman) vs static permanent magnets | L8 |
| Elliptical motion / 6-phasor ellipse → `M=O·D·Oᵀ` ellipsoid bridge | L2 |
| "Knotted rubber band" analogy (topology + energy → oscillation) | L7 (seed) |
| "Spinning arrow through a point" visual (rotational, not translational) | L7 (seed) |
| Spinning (ω) vs oscillating (A & ω); spin-½; de Broglie λ | L7 |
| The 4th dimension: gravity (`g`) + the time axis; the two "times" | L3 |
| 4D & teleparallelism | L3 |
| Handedness / chirality; matter-antimatter; composite particles | L10 |
