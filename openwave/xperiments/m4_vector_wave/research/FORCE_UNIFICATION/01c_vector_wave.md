# 🔶 PHASE 1c: Vector Wave Force — 3D Displacement Research

## Goal

Implement true vector displacement (independent L and T components) to test whether L→T spin conversion produces emergent charge-dependent force direction — the ONLY mechanism from Phase 1b that distinguished charges.

## ⚠️ Key Concepts Under Investigation

- **BASE WAVE**:
  - The fundamental energy wave that fills all of space — the energy source for matter formation, force emergence, EM waves, and heat. Waves arriving from all directions create an isotropic field. In 1D, validated as a standing wave (`ψ = A₀·cos(kx)·cos(ωt)`) via Laplacian self-stabilization. In 3D, the isotropic superposition from all directions creates a uniform energy density field with vector displacement at every point. The base wave equation in 3D and its exact form are still to be determined — this is Step 1 of Phase 1c

- **WAVE CENTERS (WCs)**:
  - Locations in 3D space where the base wave is elastically disturbed — the wave passes through and comes out changed (not just reflected). Phase 1b tested 10 disturbance models: passive (reflect, absorb, clamp, scatter) all failed; elastic L→T spin conversion is the only charge-sensitive mechanism found. The WC equation and exact disturbance mechanism are still being explored — this is Step 2 of Phase 1c

- **3D SPACE & VECTOR WAVES**:
  - In 3D, spherical wave interference from the base wave + WC disturbance always promotes elliptical granule motion — even one WC + base wave is enough. Energy is not scalar in 3D — it requires a vector field with independent longitudinal (L) and transverse (T) components: `E = E_L + E_T`. Scalar models (M3/1D) collapse directional information into magnitude, producing correct 1/r² scaling but wrong force direction. Vector displacement is required to capture charge-dependent force — this is why Phase 1c must work in 3D

- **FORCE MECHANISM**:
  - One unified force `F = -∇E`, two directions. `∇E_L` = electric force (radial/longitudinal gradient), `∇E_T` = magnetic force (perpendicular/transverse gradient). L/T defined relative to radial direction from the WC experiencing the force. Electric always present (one radial direction). Magnetic conditional — transverse has 360° freedom, cancels unless aligned/coherent (moving electrons, permanent magnets, spin coherence). Gravity = residual total energy deficit (omnidirectional, always attractive, weak)

- **ELLIPTICAL MOTION & SPIN**:
  - The ellipse at each point encodes L amplitude (semi-major), T amplitude (semi-minor), handedness (CW/CCW rotation direction = spin up/down = electron/positron), and ellipse plane orientation. Described by 6 phasor numbers (R_x, R_y, R_z, Φ_x, Φ_y, Φ_z). Complex sinusoids in QM encode this naturally: real = L, imaginary = T, |ψ|² = E_L + E_T. The imaginary unit i IS the 90° quadrature relationship between L and T components. Spin is the L→T wave transformation at the WC (Wolff: "the wave is spinning, not the particle"), with 720° spherical rotation required by 3D space geometry

- **NON-LINEAR WAVES (Phase 1d)**:
  - Variable λ(r) near WCs — wavelength changes with distance from the WC (Yee & Hauger shells `r_wavelength = 2Κλ - 2nλ`, WKB phase integral). The current energy formula `E = ρV(fA)²` uses constant f and can't see λ variation. Phase 1d implements `E = ρV(c·A/λ(r))²` where the `∇λ` term creates force from wavelength gradients. Converges with Phase 1c — variable λ breaks sinc periodicity while vector displacement provides charge sensitivity

- **ENERGY CONSERVATION & WAVE STEEPNESS**:
  - Wave steepness A/λ = constant for isolated energy redistribution. External energy input (heating) increases steepness. WC spin converts L→T while conserving `E_L + E_T = const`. The conversion ratio may be the fine-structure constant α

- **VISUALIZATION & RENDERING**:
  - Phase 1c focuses on raw math (numpy scripts) for an accelerated discovery phase — compute, sweep, analyze, document. No matplotlib animation, no 3D rendering. When the wave equations are validated and force emergence is confirmed, the equations get ported to OpenWave's M4 Taichi-based 3D rendering infrastructure (launcher, medium, wave_engine, particles, force & motion) for visualization and public demonstration

- **VALIDATION TARGETS** (the complete force unification validation set):
  - **Far-field (Coulomb regime)**: force direction emergent (not imposed ±1 — passes emergence test). 1/r² magnitude scaling vs Coulomb reference. Newton's 3rd law (equal and opposite). Opposite charge attracts, same charge repels — at ALL separations (no sinc flips). Energy conserved
  - **Near-field (particle regime)**: same-phase lock-in — oscillatory force creates stable energy wells (quarks, orbital shells, bonding). Opposite-phase monotonic attraction → annihilation (wave cancellation at zero separation). Near-field to far-field transition boundary — clean crossover from lock-in to Coulomb
  - **Gravitational**: 10⁻⁴² EM-to-gravitational force ratio emerges from the model (shading / density deficit). Computed G matches Smoliński's Scilab reference values
  - **Magnetic**: transverse force component (∇E_T) reproduces magnetic field geometry. Spin alignment → coherent T field → detectable magnetism. Perpendicular to electric (∇E_L) at 90°
  - **Matter formation**: attraction + repulsion equilibrium holds particles together. Standing wave core forms at WC with correct radius (K²λ). Mass = energy in standing waves (E = mc²)

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

Math scripts to build, what to compute, what to measure.

### Step 1: 3D Vector Base Wave

Implement the 3D isotropic base wave with vector displacement:

- [ ] 3D grid with vector displacement at each point: `ψ(r) = (ψ_x, ψ_y, ψ_z)`
- [ ] Isotropic base wave: waves from all directions → uniform energy density
- [ ] L/T decomposition relative to any reference point
- [ ] Verify: base wave alone produces flat energy, zero force

### Step 2: WC as L→T Converter (Spin)

- [ ] Place WC at grid center
- [ ] Implement L→T conversion: reduce A_L, increase A_T at WC, with 1/r rolloff
- [ ] Charge-dependent direction: CW (+½) vs CCW (-½)
- [ ] Verify: single WC creates energy concentration (standing waves) near it, drainage far field
- [ ] Verify: energy conservation `E_L + E_T = constant`

### Step 3: Two-WC Force Test

- [ ] Place two WCs at variable separation
- [ ] Compute `E_total = E_L + E_T` at each grid point
- [ ] Compute `F = -∇E_total` at WC positions
- [ ] Sweep separation: does force direction depend on charge sign?
- [ ] The critical test: opposite charge attracts, same charge repels — emergent, not imposed?

### Step 4: Coulomb Validation

- [ ] Compare force magnitude vs Coulomb reference at multiple separations
- [ ] Check 1/r² scaling
- [ ] Check Newton's 3rd law (equal and opposite forces)

### Step 5: Convergence with Phase 1d (Variable λ)

- [ ] Add λ(r) to energy equation: `E = ρV(c·A/λ(r))²`
- [ ] Test combined: vector displacement + variable λ → does force direction + magnitude work?

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
