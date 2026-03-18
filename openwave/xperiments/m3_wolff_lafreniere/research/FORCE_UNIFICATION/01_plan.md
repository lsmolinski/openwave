# PLAN: Force Unification from Wave Interference

## GOALS

PRIMARY (EMERGENT FORCES)

- Provide numerical evidence that all fundamental FORCES (electric, magnetic, gravitational, strong) emerge from Energy Wave interference patterns in a spacetime medium — validating Energy Wave Theory (EWT) through simulation.

SECONDARY (EMERGENT WAVES)

- Provide numerical evidence that ELECTROMAGNETIC WAVES and HEAT emerge from Energy Wave.

## ROADMAP

### Phase 0: Tooling

- ✅ Build 1D wave sandbox (matplotlib, interactive controls)
- ✅ Phasor superposition (analytical amplitude, replaces EMA-RMS)
- ✅ Coulomb reference comparison in sandbox
- ✅ Decision: weighted partial standing wave as primary equation
- [ ] Build physics invariant tests (pytest, boundary limits, energy conservation)
- [ ] Validate phasor RMS: single WC sinc envelope, same/opposite charge interference, EMA-RMS equivalence

### Phase 1: Electric Force — Far-Field (1D Sandbox only)

- [ ] Resolve far-field oscillatory force (MAIN BLOCKER — sinc nodes in out-wave)
- [ ] Test energy gradient ∇E directly vs amplitude gradient A·∇A
- [ ] Test gradient sampling radius (1 grid vs 1λ vs 2λ window)
- [ ] Test pressure/velocity gradient (90° phase-shifted from displacement)
- [ ] Test standing vs traveling wave decomposition (force from each component separately)
- [ ] Investigate multi-variable energy gradient: ∇A + ∇f + ∇ρ contributions
- [ ] Validate 1/r² force scaling (Coulomb law match)
- [ ] Validate force direction: opposite charges attract, same charges repel (consistently)
- [ ] Parameter sweep: force vs separation from 2λ to 10λ
- [ ] Plot energy density landscape along axis at various separations

### Phase 2: Electric Force — Near-Field (1D Sandbox only)

- [ ] Same-phase lock-in: verify oscillatory force creates stable energy wells
- [ ] Opposite-phase monotonic attraction: verify consistent attraction to annihilation
- [ ] Characterize near-field → far-field transition boundary
- [ ] Validate dual-treatment: raw phasor (near-field) vs smoothed envelope (far-field)

### Phase 3: Electric Force — 3D Validation (M3 Taichi, port from 1D)

- [ ] Port validated 1D equations to M3 3D wave engine
- [ ] Reproduce 1D force results in 3D (2 WCs on axis)
- [ ] Test off-axis configurations (verify spherical symmetry)
- [ ] Validate force & motion integration (particles move correctly)
- [ ] Compare against Coulomb force at multiple separations

### Phase 4: Non-Linear Wave Equations (M3, 1D → 3D)

- [ ] Test variable λ(r) (Yee & Hauger discrete wavelength shells)
- [ ] Test Smoliński r⁵ energy scaling near WC core
- [ ] Evaluate impact on near-field lock-in and far-field force scaling
- [ ] Compare 5 wave equation forms under same test configuration

### Phase 5: Gravitational Force (M3, Multi-Particle)

- [ ] Test wave shading with particle clusters
- [ ] Test Smoliński buoyancy model: ρ(x) and f(x) as local variables
- [ ] Validate 10⁻⁴² EM-to-gravitational force ratio
- [ ] Validate computed G against Smoliński's Scilab reference values

### Phase 6: Magnetic Force (M4 Vector Waves)

- [ ] Build M4 vector wave engine with transverse displacement
- [ ] Validate elliptical displacement trajectories (6 phasor numbers)
- [ ] Model spin as toroidal wave flow
- [ ] Demonstrate magnetic force from transverse wave interference
- [ ] Separate longitudinal (electric) and transverse (magnetic) force components

### Phase 7: Time Dynamics (M4 or new method)

- [ ] Implement variable λ per voxel (local dt)
- [ ] Demonstrate time dilation from energy starvation mechanism
- [ ] Connect λ modulation → granule velocity → pressure → gravity
- [ ] Test force control via frequency/spin manipulation

### Phase 8: Emergent Waves

- [ ] Demonstrate photon-like traveling wave packets
- [ ] Test thermal energy as standing wave dynamics
- [ ] Validate electromagnetic wave emergence from medium disturbances

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

---

## PHASE DETAILS Details

Details for each roadmap phase. Checklist items live in the ROADMAP above — details here provide context only.

### Phase 0: Tooling — Details

**Physics invariant tests** (pytest): validate each wave equation form against known physical properties before deploying to 3D engine:

- Dimensional consistency of all force/field equations
- Boundary behavior: correct limits at r→0 (center voxel) and r→∞ (far-field decay)
- Energy conservation: total energy in the field should be constant for a static configuration
- Phasor equivalence: phasor RMS must match EMA-RMS after convergence
- Near-field/far-field transition: envelope must be continuous across the weight function boundary
- Charge symmetry: same-charge pair must produce repulsion, opposite-charge pair must produce attraction in the far-field

The 1D engine is ideal for these tests — fast execution, no GPU overhead, and deterministic results.

**Phasor validation** (prerequisite to force work): before chasing force emergence, confirm the phasor gives correct amplitude patterns — single particle sinc envelope, same-charge destructive interference, opposite-charge constructive interference, EMA-RMS equivalence.

### Phase 1: Electric Force — Far-Field — Details

The main blocker is the far-field oscillatory force: the sinc function sin(kr)/kr in the out-wave creates permanent nodes in the phasor RMS at all distances. Force direction flips every λ/2 of separation change, even where only smooth 1/r decay should exist. Confirmed in both 3D and 1D engines.

**Candidate solutions to test** (see 05_1D_wave.md for full analysis):

- Energy gradient ∇E directly vs amplitude gradient A·∇A — squaring A² changes oscillation structure
- Gradient sampling radius — 1 grid unit captures sinc detail a real particle wouldn't feel, try 1λ-2λ windows
- Pressure/velocity gradient — force may be 90° phase-shifted from displacement, following pressure nodes not displacement nodes
- Multi-variable energy gradient — E = ρV(fA)² has three variables (A, f, ρ) that can all create spatial gradients independently
- Standing vs traveling wave decomposition — decompose phasor into standing-only and traveling-only amplitudes, test force from each
- Compare all 5 wave equation forms under same test configuration

**Validation targets**: plot E = ρV(fA)² along the axis connecting two particles at various separations, identify constructive/destructive interference locations, verify gradient direction and 1/r² magnitude scaling against Coulomb reference.

### Phase 2: Electric Force — Near-Field — Details

Near-field has two expected behaviors based on charge phase:

- **Same phase (same charge)**: oscillatory lock-in — alternating attraction/repulsion every λ/2 creates stable energy wells where particles are trapped (quarks, orbital shells, bonding)
- **Opposite phase (opposite charge)**: monotonic attraction — consistent pull toward annihilation at zero separation where waves cancel completely

Currently both configurations show the same oscillatory behavior — the sinc node structure overrides the charge-phase signal. Resolving Phase 1 (far-field) may also resolve this.

The dual-treatment boundary (raw phasor for near-field, smoothed for far-field) needs implementation and tuning. The weight function transition parameter may serve double duty.

### Phase 3: Electric Force — 3D Validation — Details

Port the validated 1D equations and force computation to M3 Taichi 3D engine. Verify that 3D results match 1D on-axis results, then test off-axis configurations for spherical symmetry. Validate force & motion integration — particles should move correctly under computed forces.

### Phase 4: Non-Linear Wave Equations — Details

Test variable λ(r) from two sources:

- **Yee & Hauger**: discrete wavelength shells r_n = 2(K-n)λ — changes interference pattern and may correct force scaling (WKB/eikonal phase integral approach)
- **Smoliński**: r⁵ energy scaling inside soliton's Energy Domain — defines how λ(r) varies near the wave center core

Evaluate impact on both near-field lock-in stability and far-field force scaling.

### Phase 5: Gravitational Force — Details

Two competing gravity models:

1. **Shading / directional attenuation**: WCs absorb part of in-wave, re-emit attenuated out-wave toward other particles → directional energy deficit → gravity. Likely requires M4 vector displacement (directional attenuation can't be represented in scalar M3)

2. **Buoyancy / medium density (Smoliński)**: Wave equation unchanged. WCs modify local ρ and f via λ modulation → granule velocity → pressure. Gravity emerges from E = ρV(fA)² where ρ and f become local variables. Works in M3 — only the energy calculation changes

For 1D sandbox: Phases 1-3 use constant ρ and f. Phase 5 adds option for position-dependent ρ(x) and f(x) computed from the wave field.

### Phase 6: Magnetic Force — Details

Requires M4 vector wave method. Magnetic force arises from transverse wave interference. The elliptical displacement trajectories (6 phasor numbers per voxel) naturally encode both longitudinal (electric) and transverse (magnetic) components. Spin modeled as toroidal wave flow — see 06_m4_vector.md for full analysis.

### Phase 7: Time Dynamics — Details

Variable λ per voxel makes dt a local variable — simulation no longer uses uniform timesteps. Energy starvation mechanism: destructive interference drops amplitude → frequency increases (energy conservation) → local time speeds up. See 07_time_dynamics.md for full analysis.

### Phase 8: Emergent Waves — Details

Photons as traveling wave packets, thermal energy as standing wave dynamics. See 08_emergent_waves.md.
