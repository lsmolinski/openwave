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

Scripts live in `research/FORCE_UNIFICATION/scripts_vector_wave/`. Not shipped product — raw math for exploration. Results go to this doc.

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
