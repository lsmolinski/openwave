# ✅ PHASE 1c: Vector Wave Force — 3D Displacement Research

## Goal

Implement true vector displacement (independent L and T components) to test whether L→T spin conversion produces emergent charge-dependent force direction — the ONLY mechanism from Phase 1b that distinguished charges.

## Research Strategy

**Math-only exploration** — no visualization, no animation. Pure numpy scripts that:

1. Compute base wave field with L and T components at 3D grid points
1. Apply WC disturbances (L→T spin conversion, elastic)
1. Compute energy from independent components: `E = E_L + E_T`
1. Compute force: `F = -∇E_total`
1. Sweep separations, measure force direction, check charge sensitivity
1. Document findings at each step

Vector Scripts live in `research/FORCE_UNIFICATION/scripts_vector_wave/`. Not shipped product — raw math for exploration. Results go to this doc.

Old Scalar Scripts live in `research/FORCE_UNIFICATION/scripts_scalar_wave/`.

When the wave equations are validated, they get ported to the M4 Taichi engine (`wave_engine.py`) for 3D visualization.

---

## Phase 1b Carry-Forward (compact summary)

Phase 1b tested 10 WC disturbance models on 5 base wave modes in a 1D scalar sandbox:

- **Only L→T spin conversion (elastic Option G) distinguished charges** — opposite charge: 12/24 oscillates, same charge: 24/24 unclear. All other models: charge-blind, imposed ±1, or sinc-reintroduced
- **The 1D scalar sandbox has a fundamental limitation** — can't represent true independent L/T. Phasor channels (P, Q) combine into magnitude P²+Q², masking the L/T separation
- **Quadrature confirmed as strongest base wave** — two orthogonal channels = L/T duality, flat energy, encodes rotational degree of freedom (complex sinusoid connection)
- **Additive superposition ruled out** — base wave + WC waves produces same sinc oscillation regardless of base mode
- **Passive WC interactions fail** — reflections cancel in isotropic fields (M2 confirmed in 3D, re-confirmed in 1D)
- **Elastic disturbance is the mechanism** — wave character changes at WC (not reflects). Spin = wave transformation (Wolff: "the wave is spinning, not the particle")

Carry-forward tasks:

- [ ] Test energy redistribution: concentration near WC (r < K²λ), drainage in far field
- [ ] Determine how WC phase affects far-field drainage pattern → L→T spin is the candidate mechanism
- [ ] Test force emergence: drainage from WC1 disturbs WC2 → energy gradient → F = -∇E
- [ ] Validate against Coulomb reference (direction + 1/r² scaling)

---

## Force Mechanism: One Force, Two Directions

### F = -∇E Applied to Vector Components

In scalar (M3/1D), force comes from `F = -∇E` where `E = ρV(fA)²` and A is a single amplitude. In vector (M4/3D), the displacement has directional components, and energy decomposes:

```text
E_total(r) = E_L(r) + E_T(r) = ρV(f·A_L(r))² + ρV(f·A_T(r))²

F = -∇E_total = -∇E_L - ∇E_T
```

Force has TWO gradient contributions — one from each component. They're not separate force mechanisms — they're the same `F = -∇E` applied to different spatial components of the energy field:

- **∇E_L** → gradient of longitudinal energy → **electric force**
- **∇E_T** → gradient of transverse energy → **magnetic force**

One mechanism. Two directions. Two forces perceived at human scale.

### Why Energy is Not Scalar in 3D

**Energy requires a vector field representation in 3D space.** Scalar energy (`E = ρV(fA)²` with a single amplitude A) is a 1D simplification that works for a single wave source but breaks down when multiple sources exist — which is always the case in reality.

The reason: **3D spherical wave interference always promotes elliptical motion** in the medium elements. Even a single WC disturbing the isotropic base wave is enough — the base wave arrives from all directions, the WC disturbs it radially, and the superposition of omnidirectional incoming waves + radial disturbance at any off-axis point creates displacement in multiple directions → elliptical motion. Linear (scalar) motion only exists in the degenerate 1D case or exactly on-axis between two sources — not the physical default anywhere in 3D space.

This is why:

- **Complex sinusoids** (imaginary numbers) in quantum mechanics encode this vector nature — the real part is one axis of the ellipse, the imaginary part is the orthogonal axis. Complex numbers are required because energy is fundamentally vectorial in 3D
- **Spin** emerges from the handedness (CW vs CCW) of this elliptical motion — the direction the granule travels around its elliptical track. Spin up/down is the orientation of this rotation
- The **90° quadrature relationship** between electric and magnetic forces is built into the elliptical geometry — the two axes of the ellipse are perpendicular by definition, and the energy gradients along each axis produce forces in orthogonal directions

The scalar model (M3) collapsed all this directional information into magnitude — correct scaling (1/r²) but wrong direction (charge sign). The vector model (M4/Phase 1c) preserves it.

### What Defines Longitudinal vs Transverse

The decomposition is relative to the **radial direction from the WC experiencing the force**:

- **Longitudinal (L)**: granule displacement component along the radial direction (toward/away from WC) — compression/rarefaction wave
- **Transverse (T)**: granule displacement component perpendicular to radial (tangential to a sphere centered on WC) — shear wave

In 3D, the transverse plane at any point has **360° of freedom** — any direction perpendicular to the radius satisfies the definition. This explains key electromagnetic observations:

- **Electric force is always present** — the longitudinal component has exactly one direction (radial), no ambiguity, no cancellation. Every WC produces a radial energy gradient → electric force always exists
- **Magnetic force is conditional** — the transverse component can point in ANY direction in the perpendicular plane. In most matter, these directions are random → cancel to zero net field → no detectable magnetism. Only when transverse directions are **aligned/coherent** does a net transverse field emerge:
  - **Moving electrons**: current flow aligns electron spins → coherent transverse emission → magnetic field around wire
  - **Permanent magnets**: crystal lattice locks spin alignment → persistent coherent transverse field
  - **Spin coherence**: aligned spins produce net transverse wave field (magnetic); random spins cancel (non-magnetic)

### Gravitational Force as Residual

The **third force** (gravity) may emerge from the **total energy deficit** — the WC concentrates energy near itself (standing waves = mass), draining the far field. This density deficit is omnidirectional (not L or T specific) and manifests as the weak, always-attractive gravitational force:

- **Electric** = ∇E_L (longitudinal gradient, charge-dependent, strong)
- **Magnetic** = ∇E_T (transverse gradient, alignment-dependent, medium)
- **Gravitational** = ∇(E_L + E_T) deficit (total drainage, always attractive, weak — 10⁻⁴² of electric)

---

## Vector Wave Theory

### Elliptical Granule Motion in 3D

When each xyz component of displacement oscillates at the same ω:

```text
x(t) = R_x · cos(ωt + Φ_x)
y(t) = R_y · cos(ωt + Φ_y)
z(t) = R_z · cos(ωt + Φ_z)
```

The trajectory at any grid point is **always an ellipse** lying on some plane in 3D space — mathematically guaranteed (same principle as light polarization: linear, circular, elliptical are the only states for a monochromatic wave).

Degenerate cases:

- **All phases equal** (Φ_x = Φ_y = Φ_z): linear oscillation (ellipse collapsed to a line) — pure longitudinal
- **Phases differ by π/2, equal amplitudes**: circle — equal L and T
- **General case**: ellipse on a tilted plane — mixed L/T with specific handedness

### 6 Phasor Numbers Describe the Full State

The ellipse at each voxel is fully described by 6 numbers: `R_x, R_y, R_z, Φ_x, Φ_y, Φ_z`.

These encode:

- **Shape**: amplitude ratios → elongation (L/T balance)
- **Orientation**: phase differences → tilt of ellipse plane (force direction)
- **Size**: magnitudes → total energy at that point
- **Handedness**: rotation direction (CW vs CCW) → charge/spin sign

Decompose into real and imaginary parts:

```text
a = (R_x·cos(Φ_x),  R_y·cos(Φ_y),  R_z·cos(Φ_z))   — real part
b = (R_x·sin(Φ_x),  R_y·sin(Φ_y),  R_z·sin(Φ_z))   — imaginary part
```

These two vectors ARE the ellipse. Normal to the plane: `a × b`. Semi-axes from the 2×2 matrix of `a·a`, `a·b`, `b·b`.

### L/T Decomposition at Any Point

Given a WC at position `r_wc` and a field point at `r`:

```text
radial direction: r̂ = (r - r_wc) / |r - r_wc|

Longitudinal amplitude: A_L = |ψ · r̂|       (projection onto radial)
Transverse amplitude:   A_T = |ψ - (ψ · r̂)r̂|  (rejection from radial)

E_L = ρV(f·A_L)²
E_T = ρV(f·A_T)²
```

This decomposition is WC-relative — the same field point has different L/T decomposition relative to different WCs.

### Extracting the 6 Phasor Values

**Method 1: Phasor Superposition (Analytical)** — sum complex phasors from all sources in one pass. Instant and exact. Each source contributes a complex phasor per xyz component, rotated by source phase. Total phasor gives R and Φ per component directly.

**Method 2: Quadrature Sampling (from timestep simulation)** — sample displacement at two times separated by a quarter period (`t₀` and `t₀ + T/4`):

```text
s1 = displacement_x at t₀
s2 = displacement_x at t₀ + π/(2ω)

R_x = √(s1² + s2²)
Φ_x = atan2(-s2, s1) - ωt₀
```

Repeat for y and z. The phasor route is analytical (instant, exact). The quadrature route requires simulation time and is approximate if sources are moving.

### Force Emergence Mechanism (from Phase 1b)

The mechanism we're validating in Phase 1c:

1. WC1 disturbs the base wave → concentrates energy in standing wave core → creates far-field energy deficit (drainage)
2. Drainage reaches WC2 → **disturbs WC2's standing waves** — warps the energy field around WC2
3. Warped energy field creates an **energy gradient at WC2's position**
4. `F = -∇E` → WC2 moves toward lower energy density → **force and motion**

The force direction depends on HOW WC1's drainage interacts with WC2's standing waves — and this must depend on the L→T spin direction (charge). This is what Phase 1c tests with true vector displacement.

---

## Spin as L→T Conversion (Elastic Disturbance)

### The Mechanism

WCs don't just reflect waves — they **transform** them:

1. Incoming waves: pure longitudinal (compression/rarefaction from all directions)
2. At WC: longitudinal amplitude converts to transverse amplitude
3. Outgoing waves: mixed L + T (reduced longitudinal, increased transverse)
4. The transverse component is **NEW** — not in the incoming field
5. **in ≠ out** → breaks the isotropic cancellation that defeated all passive models
6. Standing waves can form from the L/T interference pattern

### Energy Conservation

```text
E ∝ A_L² + A_T² = constant
```

The conversion redistributes energy between modes without creating or destroying it. The conversion ratio may be the **fine-structure constant α** (L/T coupling strength).

### Charge from Spin Direction

- **CW rotation** (spin +½): L→+T conversion → electron (source_offset = π)
- **CCW rotation** (spin -½): L→-T conversion → positron (source_offset = 0)

The "charge" is the **direction of the L→T conversion** — not a label imposed on the particle, but an intrinsic property of HOW the wave transforms at the WC. This is Wolff's 720° spherical rotation — a property of 3D space.

### Why Complex Numbers in QM

The Schrödinger equation requires complex numbers: `ψ = ψ_real + i·ψ_imaginary`

This naturally encodes L/T duality:

- **Real part** = longitudinal amplitude
- **Imaginary part** = transverse amplitude (spin component)
- **|ψ|² = ψ_real² + ψ_imaginary²** = total energy from both components
- The **imaginary unit i** encodes the 90° phase relationship between L and T

Complex numbers in QM are not mathematical convenience — they're physics. The i IS the transverse component.

---

## On-Axis Vector Analysis (Known Limitation)

For two WCs on the x-axis, displacement vectors are anti-parallel (d̂₁ = +x̂, d̂₂ = -x̂). Only x-component is nonzero, so `|ψ_vec|` reduces exactly to `|ψ_scalar|`. No improvement for on-axis test.

```text
|ψ_vec|² = |ψ₁|² + |ψ₂|² + 2·|ψ₁|·|ψ₂|·cos(k·Δr + Δφ)·cos(θ_geo)
```

On-axis: cos(θ_geo) = -1 (constant). The vector advantage only manifests **off-axis** or with **transverse components** — which is exactly what L→T spin conversion produces.

---

## Spin, Magnetism, and Thermal Energy

**Spin** may be the unifying concept:

- **Magnetic moment from spin**: L→T conversion at WC creates transverse wave emission. Coherent spin alignment = net transverse field = magnetism
- **Thermal energy from spin**: incoherent spin fluctuations (random T directions) = heat. Higher temperature = more energetic L→T cycling
- **Gravitational shading from spin**: spin converts L→T, reducing longitudinal amplitude. This L reduction is the gravitational "shadow" — the mechanism behind shading/push-out gravity

The fine-structure constant α may be the L→T conversion ratio. If gravity emerges from repeated application of this coupling, the 10⁻⁴² EM-to-gravitational ratio emerges naturally.

### Toroidal Flow Geometry (Smoliński)

The internal "engine" of a particle operates in the Energy Domain with toroidal (doughnut-shaped) wave flows governed by r⁵ scaling:

- **Spin** = toroidal wave flow around the particle core
- **Magnetic moment** = natural dipole from circulating toroidal current
- **Thermal energy** = modulation of toroidal flow rate/precession
- The elliptical displacement trajectories may be cross-sections of toroidal flows

### Vortex Electron Model (Butto, 2021)

The electron as a superfluid irrotational vortex — supporting the toroidal flow picture:

- **Spin-½ = differential rotation**: core rotates at 2× boundary rate (720° core / 360° boundary)
- **Planck constant = vortex angular momentum**: Γ·m_e = h
- **Magnetic moment without g-factor**: μ = qcr from vortex hydrodynamics
- **Helmholtz theorems**: vortex filaments form closed loops → toroidal geometry

---

## Implementation Plan

Math scripts in `scripts_vector_wave/`. Detailed task checklist in [00_roadmap.md](00_roadmap.md).

## ✅ Step 1: 3D Vector Base Wave — COMPLETED

Built 3D isotropic base wave `ψ(r) = (ψ_x, ψ_y, ψ_z)`: N=200 longitudinal plane waves from Fibonacci sphere directions, complex vector phasor field on 64³ grid (8λ extent, 8 vox/λ). `scripts_vector_wave/step1_base_wave.py`

### Results (all 6 checks passed)

- **Mean energy**: matches theory `E = ρVf²A₀²/2` to 0.3% (ratio 1.0033)
- **Energy CV = 0.5705**: matches chi-squared(6) prediction `1/√3 ≈ 0.5774`. The "speckle" is fundamental to coherent wave interference — the energy has wavelength-scale fluctuations that average out at particle scales (electron = 100λ)
- **L/T split: E_L/E = 0.3307, E_T/E = 0.6693**: matches theory (1/3 longitudinal, 2/3 transverse). Holds at all 4 tested reference points. This is a fundamental property of isotropic vector fields — one radial direction carries 1/3 of the energy, two perpendicular directions carry 2/3
- **E_L + E_T = E_total**: exact to machine precision (err = 6.79e-16)
- **Force = speckle noise only**: |F|_mean matches CV·E·k₀ estimate (ratio 1.03). No large-scale gradient — confirmed null baseline
- **Convergence**: CV is independent of N_sources (confirmed 50 → 1000), converging to 1/√3 as predicted by chi-squared statistics

Key insight for Step 2: the undisturbed base wave has a 1/3 L, 2/3 T geometric projection at any reference point (see L/T Projection below). A WC that converts L→T via spin will shift this ratio locally, creating the energy gradient that produces force.

### Base Wave Equation

Each of N plane waves from the isotropic field:

```text
ψ_n(r,t) = (A₀/√N) · cos(k₀·k̂_n·r + ωt + φ_n) · k̂_n
```

| Term | Meaning |
| --- | --- |
| `ψ_n` | Vector displacement from the n-th plane wave — a 3D vector [am] |
| `A₀/√N` | Amplitude per source. A₀ = 0.92 am total. Divided by √N (not N) so total **energy** sums correctly: N × (A₀/√N)² = A₀². Uncorrelated waves add in quadrature (same principle as random walk: N steps of size 1 → distance √N) |
| `cos(...)` | The oscillation — a scalar that varies in space and time. Phase order: `cos(spatial + temporal + source_offset)` — OpenWave standard |
| `k₀` | Wavenumber = 2π/λ₀ — converts distance to phase [radians/am] |
| `k̂_n` | Unit direction vector — one of N Fibonacci sphere points. Represents the propagation direction of this plane wave |
| `k̂_n · r` | Dot product — projects grid point position onto propagation direction. A plane wave varies along k̂ and is constant on planes perpendicular to k̂ |
| `k₀·k̂_n·r` | Spatial phase — radians of oscillation between origin and point r along direction k̂_n |
| `+ωt` | Temporal phase, In-wave convention: `+ωt` = incoming (toward future WC), `-ωt` = outgoing (from WC). Base wave is all in-waves from the universe |
| `φ_n` | Random phase offset per source (uniform 0 to 2π, seeded for reproducibility). Represents waves from different distant matter arriving with uncorrelated phases |
| `· k̂_n` | Displacement direction — the wave displaces the medium **along** k̂_n (longitudinal). This is the critical vector part |

- **Wave Direction**: all directions from all matter of the universe, simulated using fibonacci sphere to center
- **Wave Mode**: pure longitudinal — every plane wave displaces along its propagation direction
- **Wave Geometry**: plane waves. The base wave comes from all matter in the universe at enormous distances (~10²⁶ am). At that distance, wavefront curvature across our grid (8λ = 228 am) is essentially zero — spherical waves become plane waves locally (same reason sunlight is treated as parallel rays). Spherical waves appear in Step 2 when we add a WC inside the grid
- **Sign Convention**: `cos(k₀·k̂·r + ωt + φ)` for in-waves, matching OpenWave standard throughout (weighted partial standing wave: `sin(kr + ωt)` = in-wave, `sin(kr - ωt)` = out-wave). In phasor form: in-wave uses `exp(-i·spatial_phase)`, out-wave uses `exp(+i·spatial_phase)`

### Fibonacci Sphere (Source Directions)

A method to place N points approximately uniformly on a sphere surface using the golden angle spiral. Each point is a unit vector k̂_n — one propagation direction. N=200 gives excellent isotropy (max directional bias = 0.00026). Same concept as M1's 8 cube vertices but with much better coverage (200 uniform directions vs 8 diagonal-biased directions).

The golden ratio φ = (1 + √5)/2 ≈ 1.618 produces the golden angle ≈ 137.5°. The algorithm distributes points uniformly in cos(θ) (compensating for less area near poles) and rotates each successive point by the golden angle in azimuth. Because φ is the most irrational number (hardest to approximate with fractions), no two points ever line up in regular rows or columns — same spiral pattern as sunflower seeds and pinecones.

Reference: <https://extremelearning.com.au/evenly-distributing-points-on-a-sphere/>

### L/T Projection (Geometric Baseline — No WC)

The L/T decomposition requires a **reference point** to define the radial direction — it is meaningless without one. In Step 1, the grid center is used as an arbitrary reference to validate the math; the physics enters in Step 2 when a WC defines the radial direction.

The base wave is **pure longitudinal** — every plane wave displaces along its propagation direction k̂. There are zero transverse waves in the field. The "2/3 transverse" result is a **geometric projection effect**: when waves arrive from all directions and you project their displacements onto a chosen radial direction r̂, only waves aligned with r̂ contribute fully to "longitudinal." Waves from other directions have components perpendicular to r̂ that register as "transverse" — even though each wave is individually 100% longitudinal along its own k̂.

Statistically: the average of cos²θ over a uniform sphere = 1/3. So 1/3 of displacement variance projects onto any chosen direction (L), and 2/3 onto the perpendicular plane (T). This is the **geometric baseline** — the ratio that exists before any WC creates actual transverse waves via spin conversion.

## ✅ Step 2: WC as L→T Converter (Spin) — COMPLETED

Single WC at grid origin emitting spherical out-wave with L + T components. `scripts_vector_wave/step2_single_wc.py`

### Results (all 5 checks passed)

- **Energy concentration**: 1.98x base wave energy at WC core (r ≈ 0), returns to 1.0x beyond r > 1λ
- **L/T ratio shifts dramatically at WC core**: from baseline 0.33 up to 0.64 (η=0, pure L out-wave) or down to 0.17 (η=1, pure T out-wave). The shift is local — beyond 1λ, baseline 1/3-2/3 restored
- **CW and CCW produce identical energy** for single WC — spin sign doesn't affect |P|² (|+t̂|² = |-t̂|²). Spin sign matters in Step 3 when two WCs interact
- **Physical η = α ≈ 1/137**: nearly identical to η=0 because α is so small (0.73%). The effect is real but tiny — electric force (L) overwhelmingly dominates magnetic (T)

### Effect of Conversion Fraction η on Energy at WC Core

| η | Description | Concentration | E_L/E at core | E_L/E baseline |
| --- | --- | --- | --- | --- |
| 0.0 | pure L (no spin) | 1.98x | 0.639 | 0.331 |
| α ≈ 0.0073 | physical (fine structure) | 1.98x | 0.636 | 0.331 |
| 0.1 | moderate conversion | 1.95x | 0.596 | 0.331 |
| 0.5 | equal L/T | 1.83x | 0.413 | 0.331 |
| 1.0 | pure T (complete conversion) | 1.71x | 0.166 | 0.331 |

### WC Out-Wave Equation

```text
P_wc(r) = A_wc · sinc(k₀r) · exp(+i·(k₀r + φ_wc)) · [√(1-η)·r̂ + √η·q·(ẑ×r̂)]
```

| Term | Meaning |
| --- | --- |
| `A_wc` | Out-wave amplitude = A₀ (same as base wave) |
| `sinc(k₀r)` | Spherical standing wave envelope: sin(k₀r)/(k₀r), peak = 1 at r=0, nodes at r = nλ/2 |
| `exp(+i·...)` | Out-wave phasor (positive exponent, vs in-wave negative) |
| `k₀r + φ_wc` | Spatial phase + source_offset (φ_wc = π for electron, 0 for positron) |
| `√(1-η)·r̂` | L component (radial direction, reduced by conversion fraction η) |
| `√η·q·(ẑ×r̂)` | T component (azimuthal, from spin). q = ±1 for CW/CCW. Magnitude `ẑ×r̂` = sin(θ) gives natural dipole pattern |

- **η = α ≈ 1/137**: the fine structure constant. Fraction of in-wave longitudinal energy converted to transverse by spin. From EWT: `α = e²/e₀²` = (elementary charge / Planck charge)². Reference: <https://energywavetheory.com/physics-constants/fine-structure-constant/>
- **Gravity connection**: the amplitude deficit (in-wave > out-wave) from L→T conversion accumulates to produce the gravitational coupling constant αG ≈ 2.4×10⁻⁴³. Reference: <https://energywavetheory.com/physics-constants/gravity-coupling-constants/>

## ✅ Step 2a: KEY FINDINGS — Spin Scale, Sinc Resolution, and New Direction

### The Sinc Oscillation Was Never a Bug — It's K=1 Physics

Per EWT (Jeff Yee), spin only occurs at electron scale K=10, NOT at the fundamental wave center K=1 (neutrino). A single WC has no spin, no charge, no L→T conversion. The sinc oscillation that's been our "main blocker" for 4 months is actually the **correct physics for K=1**:

| Scale | K | Charge | Spin | Far-field Coulomb? | Sinc oscillation? |
| --- | --- | --- | --- | --- | --- |
| Neutrino | 1 | neutral | none | NO | Correct — oscillation = lock-in |
| Electron | 10 | charged | yes (L→T) | YES | Must be broken by spin |

- **Neutrino (K=1)**: single WC, electrically neutral, no spin. "Not attracted to particles" — can pass through Earth without collision. The oscillatory sinc force at sub-wavelength distances IS the strong force / lock-in that holds neutrinos together into electrons
- **Electron (K=10)**: 10 WCs in 1-3-6 tetrahedral arrangement at standing wave nodes. Spin emerges when one WC goes "off node" → forced to reposition → continuous rotation through the structure. "The spin of the electron becomes the magnetic force — a new transverse wave." Spin requires energy → L→T conversion at rate α ≈ 1/137
- **Strong force**: ~137x stronger than electric (= 1/α). At standing wave nodes, full in-wave amplitude is used (no spin loss). Lock-in at sub-wavelength distances: quarks at 3λ_e, nucleons at 4λ_e
- **Force hierarchy**: all forces are ratios of the electric force. Strong = electric × (1/α). Electric = baseline (in-wave minus spin energy). Magnetic = transverse from spin. Gravitational = accumulated amplitude deficit

References: <https://energywavetheory.com/forces/electric/>, <https://energywavetheory.com/forces/unification-of-forces/>, <https://energywavetheory.com/forces/strong-force/>

### Same Phase vs Opposite Phase — Lock-in and Annihilation

The sinc interaction energy `E_int ∝ cos(k·Δr + Δφ)` naturally produces different well structures for same vs opposite phase, it's actually beautifully consistent:

```text
Same phase (Δφ = 0):     E_int ∝ +cos(k·Δr)    wells at λ/2, 3λ/2, 5λ/2 ...
Opposite phase (Δφ = π): E_int ∝ -cos(k·Δr)    wells at 0, λ, 2λ, 3λ ...
```

| | Same Phase | Opposite Phase |
| --- | --- | --- |
| r = 0 | Energy HILL (max repulsion) | Energy WELL (deepest — annihilation) |
| r = λ/2 | Energy WELL (1st lock-in) | Energy HILL (barrier) |
| r = λ | Energy HILL | Energy WELL (2nd, shallower) |

**Same phase**: lock-in at standing wave nodes → particle formation. Zero separation = maximum repulsion. This IS how 10 neutrinos form an electron (strong force binding at sub-wavelength distances).

**Opposite phase**: deepest well at r=0 → annihilation. Barriers at λ/2, 3λ/2 can temporarily trap (positronium-like bound states, lifetime ~125 ps to ~142 ns). High kinetic energy overcomes barriers → direct annihilation (particle accelerators).

### Annihilation: What We Can Compare to Observed Data

- **Electron-positron annihilation**: well-observed experimentally (positronium, accelerator collisions). K=10, has spin. This is the comparison target for our simulation
- **Neutrino-antineutrino annihilation**: NO observed data. EWT neutrino page doesn't mention it. Neutrinos are too weakly interacting to observe annihilation events
- **Neutrino observation** in general: extremely limited — neutral, passes through Earth, only detected via rare weak-force interactions. No electromagnetic or strong force data to compare against

**Implication for simulation**: we cannot validate K=1 (neutrino) force behavior against experiment because there's no observational data. The neutrino's only use case in our simulation is as a building block for the electron (K=10 tetrahedron). All observable force validation (Coulomb, annihilation, magnetic) must be at K≥10.

### Jeff Yee's Response (2026-03-23)

Key points from Jeff's reply to the Phase 1c research update:

1. **Near-field physics is the breakthrough**: "I think it's a good thing that we're getting the physics we want for near-field, because this is the key to explaining particle formation. Simulations based on real math/physics that show this is actually the real breakthrough."
2. **L→T conversion insight is novel**: "The interesting one here that I never considered is that L→T conversion causes the slight difference in symmetry."
3. **Single WC may not do Coulomb**: "two wave centers may be difficult/impossible to show with the Coulomb force because at least in my model, a single wave center doesn't do the transverse conversion. Perhaps it does and it's slight, but perhaps not — at the moment a neutrino (single wave center) is thought to be neutral anyway with no charge."
4. **Bringing in Dieter Hauger**: the co-author of the Yee & Hauger wavelength shells paper, who "came up with the idea/equation" for where standing waves convert to traveling waves. May have insights on the far-field Coulomb mechanism.
5. **Variable λ(r)** from the shell model may contribute — the ∇λ force term works independently of sinc oscillation (Phase 1d territory).

### M3 Electron Stability Test (repulsion4.py)

Tested the K=10 tetrahedral electron in M3 with fixes (leapfrog integrator, λ/2 lock spacing, gradient radius 3, damping 0.995). Still unstable because:

**15/45 WC pairs sit at non-node distances.** The tetrahedral geometry (1-3-6 arrangement) produces distances of √3×λ/2 ≈ 0.87λ and √2×λ/2 ≈ 0.71λ between non-adjacent WCs — between standing wave nodes, where the sinc force actively destabilizes.

Two paths to fix:

1. **Non-linear λ(r)** (Phase 1d): Yee & Hauger shells give variable node spacing near the WC core. The real electron tetrahedron may sit at nodes of a non-uniform λ(r) lattice, not the uniform λ/2 lattice
2. **M4 vector method** (this Phase 1c): L→T spin conversion creates transverse energy with different spatial structure than the sinc pattern. Equilibrium in a vector field (E = E_L + E_T) could accommodate geometry that scalar M3 cannot

### K Variable: Do We Need It in Phase 1c?

**For Steps 3-4 (two-WC force test + Coulomb)**: K is not needed. We test two WCs with spin (η=α) at variable separation. This tests the L→T mechanism for charge-dependent force. The WCs represent simplified "particles with spin" — not full K=10 electrons with tetrahedral structure.

**For electron stability**: K=10 with tetrahedral geometry requires modeling 10 WCs. The tetrahedral shape is defined by the `tetrahedral_10()` function in repulsion4.py. The challenge is not modeling the shape — it's that the shape doesn't fit the uniform sinc node lattice. This requires Phase 1d (variable λ) and/or vector forces.

**For validation against experiment**: all observable comparisons (Coulomb 1/r², annihilation, magnetic force) are at K≥10. K=1 physics (lock-in, neutrality) is internally consistent but has no external data to validate against.

**Recommendation**: proceed with Steps 3-4 using two WCs with spin (no K variable needed). The K=10 tetrahedral electron is a Phase 3+ challenge that requires both variable λ AND vector forces.

### EWT References for Step 2a

- <https://energywavetheory.com/subatomic-particles/neutrino/>
- <https://energywavetheory.com/subatomic-particles/electron/>
- <https://energywavetheory.com/subatomic-particles/particle-creation-and-decay/>
- <https://energywavetheory.com/forces/electric/>
- <https://energywavetheory.com/forces/unification-of-forces/>
- <https://energywavetheory.com/forces/strong-force/>

---

## ❌ Step 3: Two-WC Force Test — SPIN ALONE DOESN'T CREATE COULOMB

Tested two WCs with spin (η = α, 0.1, 0.5) at separations 2λ to 6λ. Four spin configurations (CW-CW, CCW-CCW, CW-CCW, CCW-CW). `scripts_vector_wave/step3_two_wc_force.py`

### Result: ALL configurations show MIXED force directions

No consistent charge-dependent pattern. Force direction flips with separation for all spin combinations. Same-spin does NOT always repel, opposite-spin does NOT always attract.

### Root Cause: On-Axis Limitation

The T component from spin is **perpendicular** to the axis connecting the two WCs (`t̂ = ẑ × r̂ = ŷ` on-axis). It creates transverse (magnetic-like) force, but does NOT contribute to the radial (electric) force along the connecting axis. The radial force is dominated by the L component, which has the same sinc oscillation as the scalar model.

This was already predicted in the "On-Axis Vector Analysis" section: "For two WCs on the x-axis, displacement vectors are anti-parallel. Only x-component is nonzero, so |ψ_vec| reduces exactly to |ψ_scalar|."

### Conclusion

- **Spin (L→T conversion) creates magnetic force** (transverse to WC axis) — to be validated in Phase 4
- **Spin does NOT create electric force** (radial/Coulomb) — the sinc oscillation persists in the L component
- **Electric force needs a different mechanism** — the remaining candidate is **variable λ(r)** from Phase 1d

### LaFreniere Phase Shift Discovery (from `sa_phaseshift.html`)

LaFreniere found that the electron core (one full λ diameter) creates a **λ/2 phase shift** in the wave passing through it. The medium is compressed inside the core (7x smaller volume than the first onion layer), waves accelerate through the core, and emerge phase-shifted by half a wavelength. Peaks become valleys. This phase shift — not spin — is what creates the **charge sign** (electron = shifted, positron = opposite shift).

Key findings from LaFreniere:

- The electron is a "pulsating wave center" — standing waves progressively transform to traveling waves
- The core is exactly 1λ in diameter — the first onion layer is 7x the core volume
- Medium compression inside the core causes wave acceleration (faster than c locally)
- The λ/2 phase shift is a physical consequence of core traversal geometry
- Two electrons close together: the phase shift determines constructive vs destructive interference → repulsion vs attraction
- The "capture phenomenon": equilibrium exists at λ/4 offset where attraction and repulsion balance → quark formation (gluonic field)
- Electron-positron pair: π/2 phase offset produces unidirectional radiation → magnetic fields. Direction reversed for opposite spin or λ/2 distance change → north/south poles

**Connection to Phase 1d**: the phase shift emerges from **variable λ(r) inside the core** — wavelength compression creates the phase advance. This is exactly what Phase 1d's `E = ρV(c·A/λ(r))²` with `∇λ` force term aims to capture. The variable-λ mechanism may produce the Coulomb force that spin alone cannot.

**Connection to spin**: spin adds the transverse component (magnetic force, π/2 phase offset for electron-positron radiation), but the charge sign comes from core geometry (λ/2 phase shift). Both mechanisms are needed: variable λ for electric, spin for magnetic.

Reference: `lafreniere/Gabriel_LaFreniere/sa_phaseshift.html`

## ✅ Step 4 → Phase 1d: Variable λ(r) for Coulomb Force

Steps 4-5 of Phase 1c (Coulomb validation, convergence with 1d) merge into Phase 1d. The vector displacement infrastructure from Steps 1-2 is ready. The missing mechanism for electric force is variable λ(r) — to be implemented in Phase 1d with:

- Yee & Hauger wavelength shells: `λ(n) = 2(K-n)λ` per shell
- WKB phase integral for phasor computation with variable k(r)
- Energy equation: `E = ρV(c·A/λ(r))²` where `∇λ` creates force from wavelength gradients
- LaFreniere core phase shift: λ/2 advance from compressed core (1λ diameter, 7x compression)
- Combined with vector displacement (L + T) for full force decomposition

---

## Convergence with Phase 1d (Variable λ)

Phases 1c and 1d converge into one solution:

- **Phase 1c** (this): vector displacement with independent L/T → charge-dependent force from spin direction
- **Phase 1d**: variable λ(r) in energy equation → force from wavelength gradients (`∇λ` term)

Both are needed: L→T alone may still oscillate if using sinc spatial structure. Variable λ breaks sinc periodicity. Together: charge-dependent force from spin + smooth spatial structure from variable λ = the unified answer.

---

## References

- M2 spin theory: `m2_laplace_propagation/research/14_spin_theory.md`
- M2 WC experiments: `m2_laplace_propagation/research/13_wave_center.md`
- Wolff: *Schrodinger's Universe and the Origin of the Natural Laws*
- Smoliński: *The Geometric Identity of Gravity* (soliton structure, r⁵, εM)
- Butto (2021): *A New Theory for the Essence and Origin of Electron Spin* (vortex model)
- Yee & Hauger: *The Geometry of Particle Standing Waves v1.1* (wavelength shells)
- LaFreniere: wave reflection/radiation pressure model
- Phase 1b full results: `01b_base_disturbance.md`
