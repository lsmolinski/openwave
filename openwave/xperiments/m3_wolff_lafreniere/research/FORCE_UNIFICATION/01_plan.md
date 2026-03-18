# PLAN: Force Unification from Wave Interference

## GOALS

PRIMARY (EMERGENT FORCES)

- Provide numerical evidence that all fundamental FORCES (electric, magnetic, gravitational, strong) emerge from Energy Wave interference patterns in a spacetime medium — validating Energy Wave Theory (EWT) through simulation.

SECONDARY (EMERGENT WAVES)

- Provide numerical evidence that ELECTROMAGNETIC WAVES and HEAT emerge from Energy Wave.

## ROADMAP

- ✅ [Build 1D wave engine (sandbox for rapid testing)](#-build-1d-wave-engine-sandbox-for-rapid-testing)
- [ ] Research wave equations that satisfy force unification
- [ ] [Build Physics invariant tests](#--physics-invariant-tests)
- [ ] Implement final equations on M3 Scalar Waves
- [ ] Build M4 Vector Waves
- [ ] Build Magnetic force simulation
- [ ] Develop Time Dynamics (M4 or new method)

---

## ✅ Build 1D wave engine (sandbox for rapid testing)

The full 3D Taichi simulator (m3) is powerful but slow to iterate on — each equation change requires running the GPU kernel, waiting for convergence, and visually inspecting 3D fields. A lightweight **1D wave engine** using matplotlib would dramatically accelerate equation development.

**Purpose**: A standalone Python script for rapid prototyping and validation of wave equations, phasor superposition, and force computation — without the overhead of the 3D simulator.

**Features**:

- 1D spatial domain (single axis through wave centers) with configurable resolution
- Arbitrary number of wave centers with controllable parameters:
  - Position along the axis
  - Phase offset (source_offset) for charge
  - Amplitude
- All 5 wave equation forms available as selectable options
- Real-time matplotlib animation showing:
  - Instantaneous displacement ψ(x, t) — oscillating wave
  - Phasor RMS envelope — steady-state amplitude profile
  - Energy density E(x) = ρ·V·(f·A)² profile
  - Force field F(x) = -∇E gradient
- Interactive parameter controls (sliders or keyboard):
  - WC separation distance
  - WC phase offsets
  - Weight function transition distance and sharpness
  - Wave equation selection
- Direct comparison overlays:
  - Phasor RMS vs EMA-RMS (validation)
  - Simulated force vs analytical Coulomb force (1/r² reference line)
  - Near-field vs far-field regime boundary marker

**Benefits**:

- **Fast iteration**: change equation, see result instantly — no GPU compile, no 3D rendering
- **Clean 1D profiles**: directly comparable to the 1D cross-sections in LaFreniere's reference animations
- **Isolated testing**: test each component (wave equation, phasor, force) independently without the full simulation chain
- **Parameter sweeps**: easily sweep separation distance to plot force vs distance curves for 1/r² validation
- **Debugging**: verify that phasor coefficients, weight functions, and force gradients produce correct values before deploying to the 3D engine

**Location**: `openwave/xperiments/m3_wolff_lafreniere/research/wave_engine_1D_v2.py`, inspired by the similar scrip v1, and old 1D wave engine test. Lets reuse variable naming, layout conventions and some other standards as a template, but start the equations logic from scratch, there is lot to be optimized.

## ✅ Decision: Weighted Partial Standing Wave as Primary Equation

Based on the analysis above, the **weighted partial standing wave** is selected as the primary wave equation for force unification research. The other four forms (Wolff-original, LaFreniere-Marcotte original, phase-warped Marcotte, combined Wolff-LaFreniere) remain in the code as commented reference implementations for comparison testing.

### Implications for Wave Equation Choice

The **weighted partial standing wave** form is the best candidate:

```text
ψ = A · [w(r)·sin(kr + ωt) + sin(kr - ωt)] / kr
```

Why it matches the animations:

- **Standing waves near center** (w ≈ 1): produces the fixed concentric rings visible around each particle in both animations
- **Traveling waves far out** (w → 0): the out-wave `sin(kr - ωt)/kr` propagates outward, exactly like the moving rings in the animations
- **Interference between particles**: the traveling out-waves from each source overlap in the gap. Their phase relationship (same offset = constructive, opposite offset = destructive) creates the amplitude modulation visible in the 1D cross-sections
- **Phasor captures the envelope**: the phasor RMS field should show higher values between same-charge particles and lower values between opposite-charge particles — matching the 1D profiles

The other wave equations have limitations for reproducing these patterns:

- **Wolff-original**: standing wave everywhere — no traveling component to create interference between distant particles
- **LaFreniere-Marcotte original**: transition is inherent, not sharp enough — standing character persists too far
- **Phase-warped Marcotte**: single traveling wave with core correction — no true standing wave nodes near center

## [ ] Physics invariant tests

Alongside the 1D engine, build a pytest test suite that validates each wave equation form against known physical properties before any formula changes are deployed to the 3D engine:

- Dimensional consistency of all force/field equations
- Boundary behavior: correct limits at r→0 (center voxel) and r→∞ (far-field decay)
- Energy conservation: total energy in the field should be constant for a static configuration
- Phasor equivalence: phasor RMS must match EMA-RMS after convergence
- Near-field/far-field transition: envelope must be continuous across the weight function boundary
- Charge symmetry: same-charge pair must produce repulsion, opposite-charge pair must produce attraction in the far-field

The 1D engine is ideal for these tests — fast execution, no GPU overhead, and deterministic results.

## Phase 1: Validate Phasor Amplitude Against Known Solutions

Before chasing force emergence, confirm the phasor gives correct amplitude patterns:

- [ ] Two same-charge (electron-electron) particles: verify destructive interference pattern between them
- [ ] Two opposite-charge (electron-positron) particles: verify constructive interference between them
- [ ] Single particle: verify sinc-like amplitude envelope matches wave equation choice
- [ ] Compare phasor RMS against EMA-RMS to confirm equivalence (they should match after EMA converges)

## Phase 2: Analyze Energy Density Gradients

With confirmed amplitude fields, study the energy density landscape:

- [ ] Plot E = ρV(fA)² along the axis connecting two particles at various separations
- [ ] Identify where constructive vs destructive interference occurs relative to particle positions
- [ ] Verify that the energy gradient (∇E) points in the expected direction (attractive for opposite charges, repulsive for same charges)
- [ ] Check whether the gradient magnitude follows 1/r² scaling with separation distance

## Phase 3: Isolate the Force Mechanism

Test which wave equation form produces the correct force behavior:

- [ ] Run each of the 5 wave equations with same test configuration (2 electrons, 2 positrons, electron-positron pair)
- [ ] Measure force magnitude vs separation for each
- [ ] Compare against Coulomb's law: F = ke·q₁q₂/r²
- [ ] Identify which equation (if any) produces the correct scaling
- [ ] Investigate whether standing wave nodes play a role (does the force depend on whether particles sit at nodes or antinodes?)
- [ ] Test non-linear wave equations where wavelength λ(r) varies with distance from wave center — the Yee & Hauger model predicts discrete wavelength shells: r_n = 2(K-n)λ, which changes the interference pattern and may correct the force scaling (see phasor_superposition.md, WKB/eikonal phase integral approach). Smoliński's r⁵ energy scaling inside the soliton's Energy Domain provides a specific non-linear relationship to test — this could define how λ(r) varies near the wave center core

## Phase 4: Separate Standing and Traveling Contributions

The phasor can be extended to track standing and traveling components independently:

- [ ] Decompose the weighted partial standing wave phasor into standing-only and traveling-only amplitudes
- [ ] Compute force from each component separately
- [ ] Determine which component (or combination) produces the correct force law
- [ ] Test if the standing wave component alone gives the electric force

## Phase 5: Multi-Particle Validation and Gravitational Force

Scale up to test gravitational shading and distinguish between two candidate gravity mechanisms:

**Two competing gravity models:**

1. **Shading / directional attenuation**: Wave centers absorb part of the in-wave energy and re-emit an attenuated out-wave (energy lost to magnetic/transverse wave conversion). The out-wave toward another particle is weaker than the out-wave in other directions → directional energy deficit → gravitational attraction. **Problem**: This requires attenuating the wave only in the direction of another particle while keeping all other directions unchanged — this is inherently directional and likely **requires m4 vector displacement** to implement, since scalar m3 has no concept of wave direction per voxel

2. **Buoyancy / medium density (Smoliński)**: The electrical wave equation itself does not change. Instead, wave centers modify the local medium density/pressure (via λ modulation → granule velocity → pressure). Gravity emerges from the energy calculation E = ρV(fA)² where **ρ and f become local variables** rather than constants. **Advantage for m3**: the wave equation stays the same, only the post-hoc energy and force calculation changes. This can be tested in the 1D sandbox by making ρ or f position-dependent

**Implications for 1D sandbox development:**

- For Phases 1-3 (electric force): ρ and f are constant everywhere — no changes needed
- For Phase 5 (gravity): add option to make ρ(x) and/or f(x) position-dependent, computed from the wave field itself (e.g., ρ depends on local amplitude or λ)
- This is a small extension to the energy calculation, not a wave equation rewrite

**Tests:**

- [ ] Cluster of same-charge particles: does the combined wave field show reduced amplitude behind the cluster (shading)?
- [ ] Two clusters separated by distance: does a net force emerge between them?
- [ ] Compare cluster-cluster force against single-particle-pair force — is there a residual (gravitational) component?
- [ ] Implement Smoliński's push-out/buoyancy mechanism: model gravity as a pressure deficit where solitons displace medium density, rather than pure wave attenuation
- [ ] Test buoyancy model in 1D sandbox: make ρ(x) = f(local wave field) and check if gravitational-scale force emerges from the modified energy equation
- [ ] Test density hierarchy scaling: does the 10⁻⁴² EM-to-gravitational force ratio emerge from the lattice geometry?
- [ ] Validate computed G against Smoliński's Scilab verification values (see references/MagnetismGravity.pdf)
