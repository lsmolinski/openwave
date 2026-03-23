# 🔶 PHASE 1c: Vector Wave Force — 3D Displacement Research

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

## 🚧 Step 3: Two-WC Force Test

Two WCs at variable separation. Compute `F = -∇(E_L + E_T)` at WC positions. Sweep separations and test the critical question: does force direction depend on charge sign, emergent from wave physics?

## 🚧 Step 4: Coulomb Validation

Force magnitude vs Coulomb reference, 1/r² scaling, Newton's 3rd law.

## 🚧 Step 5: Convergence with Phase 1d

Add variable λ(r) to energy equation. Test combined vector displacement + variable λ.

---

## Other Force Computation Approaches

Beyond `F = -∇E_total`, vector displacement enables alternative force quantities:

- **Divergence** (∇·ψ): compression/rarefaction — scalar but signed
- **Curl** (∇×ψ): rotational displacement — vector, related to magnetic field
- **Energy flux** (ψ × ∂ψ/∂t): directional energy flow (radiation pressure, LaFreniere's mechanism)
- **Per-component force**: `F_x = -∂E_x/∂x`, `F_y = -∂E_y/∂y`, `F_z = -∂E_z/∂z`

Standing waves have zero net flux, traveling waves have nonzero flux — flux naturally separates near-field (standing, lock-in) from far-field (traveling, Coulomb).

All Maxwell's equations have divergence and curl components.

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
