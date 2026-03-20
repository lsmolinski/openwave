# PLAN: Force Unification from Wave Interference

## GOALS

### PRIMARY (EMERGENT FORCES)

- Provide numerical evidence that all fundamental FORCES (electric, magnetic, gravitational, strong) emerge from Energy Wave interference patterns in a spacetime medium — validating Energy Wave Theory (EWT) through simulation.

### SECONDARY (EMERGENT WAVES)

- Provide numerical evidence that ELECTROMAGNETIC WAVES and HEAT emerge from Energy Wave.

## ROADMAP

## [Phase 0: Tooling](#phase-0-tooling--details)

- ✅ Build 1D wave sandbox (matplotlib, interactive controls)
- ✅ Phasor superposition (analytical amplitude, replaces EMA-RMS)
- ✅ Coulomb reference comparison in sandbox
- ✅ Decision: weighted partial standing wave as primary equation
- ✅ Build physics invariant tests (pytest, boundary limits, energy conservation)
- ✅ Validate phasor RMS: single WC sinc envelope, same/opposite charge interference, EMA-RMS equivalence
- ✅ Parameter sweep: force vs separation from 2λ to 10λ (`sweep_force_vs_separation.py`)
- ✅ Validate 1/r² force scaling — interaction energy E_int ∝ |Z₁|·|Z₂| ∝ 1/r gives F ∝ 1/r² (confirmed numerically)

## [Phase 1: Electric Force — FAR-FIELD (1D Sandbox only)](#phase-1-electric-force--far-field--details)

- [ ] Resolve far-field oscillatory force (MAIN BLOCKER — sinc nodes in out-wave)
  - ✅ Implemented F = -∇E as standard force computation (replaces A·∇A chain rule expansion)
  - ✅ Tested gradient sampling radius / Gaussian smoothing (σ = 0.25λ to 2λ) — does not resolve; destroys charge signal
  - ✅ Tested smooth envelope interaction (|Z₁|·|Z₂| with imposed charge sign) — 17/17 direction + 1/r² scaling, but sign is not emergent
  - ✅ Ruled out numerical precision — 1D uses f64; M3/M4 f32 with sim-friendly units; oscillation is real math not artifact
  - ✅ Ruled out inertia filtering — phasor RMS already IS the time-averaged quantity; problem is spatial sinc, not temporal
  - ✅ Ruled out pressure/velocity gradient — 90° shift moves sinc zeros by λ/4 but preserves λ/2 oscillation period
  - ✅ Ruled out standing vs traveling wave decomposition — each component individually is smooth (1/kr), but coherent superposition of two sources still oscillates via cos(k·Δr) interference; oscillation is intrinsic to wave interference, not to standing/traveling character
  - ✅ Tested all 5 wave equations (Wolff, LaFreniere-Marcotte, Phase-warped, Combined, Weighted) — all produce force direction flips; confirms oscillation is intrinsic to coherent interference regardless of spatial function
- [ ] Validate force direction: opposite charges attract, same charges repel (emergent, not imposed)
- [ ] Plot energy density landscape along axis at various separations
- [ ] Compare 1D profiles against LaFreniere reference animations (constructive/destructive patterns)

> **All linear scalar candidates exhausted (9/9 ruled out).** Three remaining paths in order of implementation complexity. Paths B and C are deeply connected — non-linear toroidal dynamics naturally produce vector patterns whose directional properties may carry charge information. They may converge into a single solution.

### [Phase 1a: Base Wave + Disturbance Model](#phase-1a-base-wave--disturbance-model--details)

- [ ] Implement base wave model in 1D sandbox (equation #6): constant A₀ + smooth disturbances
- [ ] Choose disturbance function δ(r): Lorentzian 1/(1+(kr)²), exponential, or weight-function-based
- [ ] Test force direction: opposite charge attraction, same charge repulsion (emergent from q = cos(phase))
- [ ] Validate 1/r² force scaling against Coulomb reference
- [ ] Run sweep_force_vs_separation.py with equation #6
- [ ] If direction + scaling confirmed, compare energy landscape against LaFreniere reference

### [Phase 1b: Non-Linear Wave Equations (1D)](#phase-1b-non-linear-wave-equations-1d--details)

- [ ] Implement variable λ(r) in 1D sandbox (Yee & Hauger discrete wavelength shells, WKB phase integral)
- [ ] Implement variable ρ(x) in 1D sandbox (density from granule velocity / wave interference)
- [ ] Test Smoliński Ψ³ cubic non-linearity (NLS soliton stabilizer)
- [ ] Test F = -∇E with spatially variable ρ(x), f(x), A(x) — all three gradients contributing
- [ ] Evaluate whether non-linear spatial structure breaks the sinc periodicity and resolves force direction
- [ ] If successful, validate against Coulomb reference (direction + 1/r² scaling)

### [Phase 1c: Vector Wave Force (M4 displacement direction)](#phase-1c-vector-wave-force--details)

- [ ] Extend 1D sandbox with vector displacement (2D: x + y components)
- [ ] Compute divergence (∇·ψ), curl, or energy flux from vector field
- [ ] Test whether signed vector quantities (divergence, flux direction) recover charge-sign without oscillatory ambiguity
- [ ] If successful, validate against Coulomb reference (direction + 1/r² scaling)
- [ ] Connect to elliptical rotation handedness (M4 phasor: same-phase = same rotation, opposite = opposite)

## [Phase 2: Electric Force — NEAR-FIELD (1D Sandbox only)](#phase-2-electric-force--near-field--details)

- [ ] Same-phase lock-in: verify oscillatory force creates stable energy wells
- [ ] Opposite-phase monotonic attraction: verify consistent attraction to annihilation
- [ ] Characterize near-field → far-field transition boundary
- [ ] Validate dual-treatment: raw phasor (near-field) vs smoothed envelope (far-field)
- [ ] Test with Verlet/leapfrog integrator (energy conservation for lock-in stability)
- [ ] Test with f64 precision (check if numerical drift causes escape from wells)

## [Phase 3: Electric Force — 3D Validation (M3 Taichi, port from 1D)](#phase-3-electric-force--3d-validation--details)

- [ ] Port validated 1D equations to M3 3D wave engine
- [ ] Reproduce 1D force results in 3D (2 WCs on axis)
- [ ] Test off-axis configurations (verify spherical symmetry)
- [ ] Validate force & motion integration (particles move correctly)
- [ ] Compare against Coulomb force at multiple separations

## [Phase 4: Non-Linear + Vector Wave Equations (M3/M4, 1D → 3D)](#phase-4-non-linear--vector-wave-equations--details)

> **Conditional scheduling**: Ports validated 1D results from Phase 1b (non-linear) and/or Phase 1c (vector) to 3D engines. If Phase 1a solved the problem, this phase introduces non-linear/vector equations as enhancements. Phases 1b and 1c may converge here — non-linear toroidal dynamics naturally produce vector patterns that carry charge information.

- [ ] Port variable λ(r) to M3 3D engine (Yee & Hauger discrete wavelength shells)
- [ ] Port variable ρ(x) to M3 3D engine (density from granule velocity)
- [ ] Port Ψ³ cubic non-linearity (NLS soliton stabilizer) to 3D
- [ ] Test Smoliński r⁵ energy scaling near WC core
- [ ] Port vector force computation (divergence/curl/flux) to M4 3D engine
- [ ] Test elliptical rotation handedness as charge-sign indicator in M4
- [ ] Evaluate impact on near-field lock-in and far-field force scaling

## [Phase 5: Gravitational Force (M3, Multi-Particle)](#phase-5-gravitational-force--details)

- [ ] Test wave shading with particle clusters
- [ ] Test Smoliński buoyancy model: ρ(x) and f(x) as local variables
- [ ] Validate 10⁻⁴² EM-to-gravitational force ratio
- [ ] Validate computed G against Smoliński's Scilab reference values

## [Phase 6: Magnetic Force (M4 Vector Waves)](#phase-6-magnetic-force--details)

- [ ] Build M4 vector wave engine with transverse displacement
- [ ] Validate elliptical displacement trajectories (6 phasor numbers)
- [ ] Model spin as toroidal wave flow
- [ ] Demonstrate magnetic force from transverse wave interference
- [ ] Separate longitudinal (electric) and transverse (magnetic) force components

## [Phase 7: Time Dynamics (M4 or new method)](#phase-7-time-dynamics--details)

- [ ] Implement variable λ per voxel (local dt)
- [ ] Demonstrate time dilation from energy starvation mechanism
- [ ] Connect λ modulation → granule velocity → pressure → gravity
- [ ] Test force control via frequency/spin manipulation

## [Phase 8: Emergent Waves](#phase-8-emergent-waves--details)

- [ ] Demonstrate photon-like traveling wave packets
- [ ] Test thermal energy as standing wave dynamics
- [ ] Validate electromagnetic wave emergence from medium disturbances

---

## PHASE DETAILS

Details for each roadmap phase. Checklist items live in the ROADMAP above — details here provide context only.

## Phase 0: Tooling — Details

**1D wave sandbox** (`wave_engine_1D_v2.py`): The full 3D Taichi simulator (m3) is powerful but slow to iterate on. The 1D sandbox provides fast iteration with instant visual feedback — no GPU compile, no 3D rendering. Clean 1D profiles are directly comparable to LaFreniere's reference animations.

Features: configurable WC array (position, phase, amplitude), all 5 wave equation forms, real-time matplotlib animation (displacement + phasor RMS overlay, energy density, force field), interactive controls (separation slider, WC on/off toggles, phase offset toggle), Coulomb reference comparison with direction-match detection, force annotations at WC positions.

**Weighted partial standing wave** selected as primary equation:

```text
ψ = A · [w(r)·sin(kr + ωt) + sin(kr - ωt)] / kr
```

Why: standing waves near center (w ≈ 1), traveling waves far out (w → 0), interference between particles via traveling out-waves, phasor captures the envelope. Other forms have limitations: Wolff (no traveling), LaFreniere-Marcotte (transition not sharp enough), Phase-warped Marcotte (no true standing nodes near center).

**Physics invariant tests** (pytest): validate wave equations before deploying to 3D engine — dimensional consistency, boundary limits (r→0, r→∞), energy conservation, phasor/EMA-RMS equivalence, near-field/far-field continuity, charge symmetry.

**Phasor validation** (prerequisite to force work): confirm correct amplitude patterns — single particle sinc envelope, same-charge destructive interference, opposite-charge constructive interference, EMA-RMS equivalence.

## Phase 1: Electric Force — Far-Field — Details

The main blocker is the far-field oscillatory force: the sinc function sin(kr)/kr in the out-wave creates permanent nodes in the phasor RMS at all distances. Force direction flips every λ/2 of separation change, even where only smooth 1/r decay should exist. Confirmed in both 3D and 1D engines.

**Force computation**: F = -∇E where E(x) = ρ·V·(f·A(x))². Computing from ∇E directly (not chain-rule expansion) ensures future variable ρ(x), f(x), λ(x) are automatically captured.

**Tested and ruled out** (see 05_1D_wave.md for full analysis):

- ✅ Gradient sampling / Gaussian smoothing — destroys charge signal along with oscillation
- ✅ Smooth envelope with imposed charge sign — perfect 1/r² + direction, but sign is not emergent (equivalent to old analytical envelope)
- ✅ Numerical precision — 1D uses f64, M3/M4 use f32 with simulation-friendly units; oscillation is real math, not floating-point artifact
- ✅ Inertia filtering — phasor RMS already IS the time-averaged amplitude; the oscillatory problem is spatial (sinc structure), not temporal

- ✅ Pressure/velocity gradient — 90° phase shift moves sinc zeros by λ/4 but preserves λ/2 oscillation period; no benefit
- ✅ Standing vs traveling wave decomposition — each component individually is smooth (1/kr), but coherent superposition of two sources oscillates via cos(k·Δr); oscillation is intrinsic to wave interference, not standing/traveling character

- ✅ Alternative wave equations — all 5 forms (Wolff, LaFreniere-Marcotte, Phase-warped, Combined, Weighted) produce force direction flips. Confirms the oscillation is intrinsic to coherent wave interference regardless of spatial function choice

**All linear candidates exhausted.** The oscillatory force is a fundamental consequence of coherent wave interference — no linear wave equation or operation on the superposed field can eliminate it while preserving charge-dependent direction.

**Validation targets**: plot E(x) along the axis connecting two particles at various separations, identify constructive/destructive interference locations, verify gradient direction and 1/r² magnitude scaling against Coulomb reference.

## Phase 1a: Base Wave + Disturbance Model — Details

**Paradigm shift**: WCs are not wave SOURCES emitting into empty space. Space is filled with a pre-existing **base wave** (fundamental energy wave, isotropic standing waves, uniform A₀). WCs create DISTURBANCES in this field — redistributing energy, not injecting it.

**Why this is different from all 9 ruled-out approaches**: all previous candidates used additive coherent wave superposition (ψ = Σ waves), which always produces oscillatory interference via cos(k·Δr). The base wave model uses **smooth amplitude modulation** on a pre-existing field:

```text
A_total = A₀ + Σ q_n · A_peak · δ(r_n)
E = ρVf²·A_total²  →  dominant force term = 2A₀ · q₂ · ∇δ₂
```

The force is LINEAR in charge sign q = cos(phase), SMOOTH (no oscillation), and AMPLIFIED by A₀. See [05_1D_wave.md](05_1D_wave.md#possible-solution-base-wave--disturbance-model-phase-1a--new-idea) for full mathematical analysis, force direction derivation, and connection to EWT/LaFreniere/gravitational shading.

**Implementation**: add equation #6 to `wave_engine_1D_v2.py` with base amplitude A₀, smooth disturbance function δ(r), and charge sign q = cos(phase). Test with sweep script.

## Phase 1b: Non-Linear Wave Equations (1D) — Details

**Rationale**: All linear operations on the sinc function preserve its λ/2 periodicity. Only a non-linear wave equation — where the spatial structure itself is no longer a pure sinc — can break the oscillatory pattern. This is the last resort before accepting that the linear weighted partial standing wave cannot produce emergent charge-sign forces.

**Three energy gradient variables** — currently only A varies spatially; making ρ and f spatial variables turns E = ρV(fA)² into a multi-variable field where ∇E captures contributions from all three:

- **A** (amplitude): current phasor RMS — carries sinc oscillation
- **f / λ** (frequency): Yee & Hauger discrete wavelength shells give λ(r) = 2(K-n)λ per shell. WKB phase integral replaces kr with ∫k(r')dr'. This changes the sinc spatial structure — nodes are no longer equally spaced, breaking the periodic force oscillation
- **ρ** (density): granule velocity (∂ψ/∂t) determines local medium density/pressure. Wave interference changes granule velocities → changes local ρ. The ∇ρ contribution to ∇E has a different spatial structure than ∇A and may carry force information that amplitude alone cannot provide

**Implementation in 1D sandbox**: replace constant k in phasor coefficients with integrated phase from variable λ(r). Add ρ(x) field computed from local wave state. F = -∇E remains unchanged — it automatically captures all variable contributions.

See [05_1D_wave.md](05_1D_wave.md#possible-solution-non-linear-wave-equations-phase-1b-fallback--phase-4) for full analysis.

## Phase 2: Electric Force — Near-Field — Details

Near-field has two expected behaviors based on charge phase:

- **Same phase (same charge)**: oscillatory lock-in — alternating attraction/repulsion every λ/2 creates stable energy wells where particles are trapped (quarks, orbital shells, bonding)
- **Opposite phase (opposite charge)**: monotonic attraction — consistent pull toward annihilation at zero separation where waves cancel completely

Currently both configurations show the same oscillatory behavior — the sinc node structure overrides the charge-phase signal. Resolving Phase 1 (far-field) may also resolve this.

The dual-treatment boundary (raw phasor for near-field, smoothed for far-field) needs implementation and tuning. The weight function transition parameter may serve double duty.

## Phase 1c: Vector Wave Force — Details

**Problem**: F = -∇(|ψ|²) uses scalar magnitude, which discards vector direction information. On-axis, vector reduces to scalar (cos(θ_geo) = -1, constant). No help for the standard test case.

**Opportunity**: vector displacement carries information beyond magnitude — ellipse rotation direction (handedness), divergence (compression), curl (rotation), energy flux direction. These are **signed quantities** that could recover charge-phase information without the oscillatory cos(k·Δr) ambiguity.

**Why this might work**: all of classical EM is built on vector field operations (right-hand rule, cross products, curl). Spin, ellipses, spirals, toroids pervade quantum mechanics. The scalar |ψ|² approach may be fundamentally insufficient — charge-sign information may ONLY exist in vector field structure.

**Connection to Phase 1b (non-linear)**: non-linear internal dynamics (toroidal r⁵ flows, Ψ³ stabilization) naturally produce the elliptical/toroidal vector patterns whose directional properties carry charge info. Phases 1b and 1c may converge.

**Implementation**: extend 1D sandbox with 2D vector displacement (x + y components), compute signed vector quantities (∇·ψ, ∇×ψ, flux), test force from these. See [05_1D_wave.md](05_1D_wave.md#possible-solution-vector-wave-force-m4-displacement-direction) for full analysis.

## Phase 3: Electric Force — 3D Validation — Details

Port the validated 1D equations and force computation to M3 Taichi 3D engine. Verify that 3D results match 1D on-axis results, then test off-axis configurations for spherical symmetry. Validate force & motion integration — particles should move correctly under computed forces.

## Phase 4: Non-Linear + Vector Wave Equations — Details

Ports validated 1D results from Phase 1b (non-linear) and/or Phase 1c (vector) to 3D engines. Phases 1b and 1c are deeply connected and may converge here.

**Non-linear (from Phase 1b → M3 3D):**

- **Yee & Hauger**: discrete wavelength shells r_n = 2(K-n)λ — non-uniform node spacing breaks sinc periodicity (WKB/eikonal phase integral)
- **Smoliński**: r⁵ energy scaling inside soliton's Energy Domain — defines how λ(r) varies near the wave center core
- **Ψ³ cubic non-linearity**: NLS soliton stabilizer, modifies spatial function from pure sinc
- **Variable ρ(x)**: from granule velocity — ∇ρ contributes to ∇E with different spatial structure than ∇A

**Vector (from Phase 1c → M4 3D):**

- **Divergence/curl/flux** force computation in full 3D vector field
- **Elliptical rotation handedness** as charge-sign indicator (6-phasor model)
- **Toroidal flow geometry**: non-linear + vector converge in the Energy Domain

See [Phase 1b](#phase-1b-non-linear-wave-equations-1d--details), [Phase 1c](#phase-1c-vector-wave-force--details), and [05_1D_wave.md](05_1D_wave.md#possible-solution-non-linear-wave-equations-phase-1b-fallback--phase-4) for full analysis.

## Phase 5: Gravitational Force — Details

Two competing gravity models:

1. **Shading / directional attenuation**: WCs absorb part of in-wave, re-emit attenuated out-wave toward other particles → directional energy deficit → gravity. Likely requires M4 vector displacement (directional attenuation can't be represented in scalar M3)

2. **Buoyancy / medium density (Smoliński)**: Wave equation unchanged. WCs modify local ρ and f via λ modulation → granule velocity → pressure. Gravity emerges from E = ρV(fA)² where ρ and f become local variables. Works in M3 — only the energy calculation changes

For 1D sandbox: Phases 1-3 use constant ρ and f. Phase 5 adds option for position-dependent ρ(x) and f(x) computed from the wave field.

## Phase 6: Magnetic Force — Details

Requires M4 vector wave method. Magnetic force arises from transverse wave interference. The elliptical displacement trajectories (6 phasor numbers per voxel) naturally encode both longitudinal (electric) and transverse (magnetic) components. Spin modeled as toroidal wave flow — see 06_m4_vector.md for full analysis.

## Phase 7: Time Dynamics — Details

Variable λ per voxel makes dt a local variable — simulation no longer uses uniform timesteps. Energy starvation mechanism: destructive interference drops amplitude → frequency increases (energy conservation) → local time speeds up. See 07_time_dynamics.md for full analysis.

## Phase 8: Emergent Waves — Details

Photons as traveling wave packets, thermal energy as standing wave dynamics. See 08_emergent_waves.md.
