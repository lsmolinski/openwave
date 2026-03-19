# PLAN: Force Unification from Wave Interference

## GOALS

### PRIMARY (EMERGENT FORCES)

- Provide numerical evidence that all fundamental FORCES (electric, magnetic, gravitational, strong) emerge from Energy Wave interference patterns in a spacetime medium — validating Energy Wave Theory (EWT) through simulation.

### SECONDARY (EMERGENT WAVES)

- Provide numerical evidence that ELECTROMAGNETIC WAVES and HEAT emerge from Energy Wave.

## ROADMAP

### [Phase 0: Tooling](#phase-0-tooling--details)

- ✅ Build 1D wave sandbox (matplotlib, interactive controls)
- ✅ Phasor superposition (analytical amplitude, replaces EMA-RMS)
- ✅ Coulomb reference comparison in sandbox
- ✅ Decision: weighted partial standing wave as primary equation
- ✅ Build physics invariant tests (pytest, boundary limits, energy conservation)
- ✅ Validate phasor RMS: single WC sinc envelope, same/opposite charge interference, EMA-RMS equivalence
- ✅ Parameter sweep: force vs separation from 2λ to 10λ (`sweep_force_vs_separation.py`)
- ✅ Validate 1/r² force scaling — interaction energy E_int ∝ |Z₁|·|Z₂| ∝ 1/r gives F ∝ 1/r² (confirmed numerically)

### [Phase 1: Electric Force — FAR-FIELD (1D Sandbox only)](#phase-1-electric-force--far-field--details)

- [ ] Resolve far-field oscillatory force (MAIN BLOCKER — sinc nodes in out-wave)
  - ✅ Implemented F = -∇E as standard force computation (replaces A·∇A chain rule expansion)
  - ✅ Tested gradient sampling radius / Gaussian smoothing (σ = 0.25λ to 2λ) — does not resolve; destroys charge signal
  - ✅ Tested smooth envelope interaction (|Z₁|·|Z₂| with imposed charge sign) — 17/17 direction + 1/r² scaling, but sign is not emergent
  - ✅ Ruled out numerical precision — 1D uses f64; M3/M4 f32 with sim-friendly units; oscillation is real math not artifact
  - ✅ Ruled out inertia filtering — phasor RMS already IS the time-averaged quantity; problem is spatial sinc, not temporal
  - ✅ Ruled out pressure/velocity gradient — 90° shift moves sinc zeros by λ/4 but preserves λ/2 oscillation period
  - [ ] Test standing vs traveling wave decomposition (force from each component separately)
- [ ] Validate force direction: opposite charges attract, same charges repel (emergent, not imposed)
- [ ] Plot energy density landscape along axis at various separations
- [ ] Compare 1D profiles against LaFreniere reference animations (constructive/destructive patterns)

> **Fallback path**: If all Phase 1 linear approaches fail to resolve the oscillatory force, proceed to [Phase 1b: Non-Linear Wave Equations (1D)](#phase-1b-non-linear-wave-equations-1d--details) before Phase 2. If Phase 1 succeeds, non-linear equations remain as [Phase 4](#phase-4-non-linear-wave-equations-m3-1d--3d--details) (later optimization).

### [Phase 1b: Non-Linear Wave Equations (1D — FALLBACK)](#phase-1b-non-linear-wave-equations-1d--details)

- [ ] Implement variable λ(r) in 1D sandbox (Yee & Hauger discrete wavelength shells, WKB phase integral)
- [ ] Implement variable ρ(x) in 1D sandbox (density from granule velocity / wave interference)
- [ ] Test F = -∇E with spatially variable ρ(x), f(x), A(x) — all three gradients contributing
- [ ] Evaluate whether non-linear spatial structure breaks the sinc periodicity and resolves force direction
- [ ] If successful, validate against Coulomb reference (direction + 1/r² scaling)

### [Phase 2: Electric Force — NEAR-FIELD (1D Sandbox only)](#phase-2-electric-force--near-field--details)

- [ ] Same-phase lock-in: verify oscillatory force creates stable energy wells
- [ ] Opposite-phase monotonic attraction: verify consistent attraction to annihilation
- [ ] Characterize near-field → far-field transition boundary
- [ ] Validate dual-treatment: raw phasor (near-field) vs smoothed envelope (far-field)
- [ ] Test with Verlet/leapfrog integrator (energy conservation for lock-in stability)
- [ ] Test with f64 precision (check if numerical drift causes escape from wells)

### [Phase 3: Electric Force — 3D Validation (M3 Taichi, port from 1D)](#phase-3-electric-force--3d-validation--details)

- [ ] Port validated 1D equations to M3 3D wave engine
- [ ] Reproduce 1D force results in 3D (2 WCs on axis)
- [ ] Test off-axis configurations (verify spherical symmetry)
- [ ] Validate force & motion integration (particles move correctly)
- [ ] Compare against Coulomb force at multiple separations

### [Phase 4: Non-Linear Wave Equations (M3, 1D → 3D)](#phase-4-non-linear-wave-equations-m3-1d--3d--details)

> **Conditional scheduling**: If Phase 1b was triggered (fallback), this phase ports the validated non-linear 1D equations to 3D. If Phase 1 succeeded with linear equations, this phase introduces non-linear equations fresh in both 1D and 3D.

- [ ] Port variable λ(r) to M3 3D engine (Yee & Hauger discrete wavelength shells)
- [ ] Port variable ρ(x) to M3 3D engine (density from granule velocity)
- [ ] Test Smoliński r⁵ energy scaling near WC core
- [ ] Evaluate impact on near-field lock-in and far-field force scaling
- [ ] Compare 5 wave equation forms under same test configuration

### [Phase 5: Gravitational Force (M3, Multi-Particle)](#phase-5-gravitational-force--details)

- [ ] Test wave shading with particle clusters
- [ ] Test Smoliński buoyancy model: ρ(x) and f(x) as local variables
- [ ] Validate 10⁻⁴² EM-to-gravitational force ratio
- [ ] Validate computed G against Smoliński's Scilab reference values

### [Phase 6: Magnetic Force (M4 Vector Waves)](#phase-6-magnetic-force--details)

- [ ] Build M4 vector wave engine with transverse displacement
- [ ] Validate elliptical displacement trajectories (6 phasor numbers)
- [ ] Model spin as toroidal wave flow
- [ ] Demonstrate magnetic force from transverse wave interference
- [ ] Separate longitudinal (electric) and transverse (magnetic) force components

### [Phase 7: Time Dynamics (M4 or new method)](#phase-7-time-dynamics--details)

- [ ] Implement variable λ per voxel (local dt)
- [ ] Demonstrate time dilation from energy starvation mechanism
- [ ] Connect λ modulation → granule velocity → pressure → gravity
- [ ] Test force control via frequency/spin manipulation

### [Phase 8: Emergent Waves](#phase-8-emergent-waves--details)

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

**Remaining candidates to test (linear):**

- Standing vs traveling wave decomposition — decompose phasor into standing-only and traveling-only amplitudes, test force from each

**Fallback candidate (non-linear — requires implementation):**

- Non-linear wave equations — variable λ(r) and ρ(x) make the sinc spatial structure non-periodic, potentially resolving the force direction problem. See [Phase 1b](#phase-1b-non-linear-wave-equations-1d--details) and [05_1D_wave.md](05_1D_wave.md#possible-solution-non-linear-wave-equations-phase-1b-fallback--phase-4)

**Validation targets**: plot E(x) along the axis connecting two particles at various separations, identify constructive/destructive interference locations, verify gradient direction and 1/r² magnitude scaling against Coulomb reference.

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

## Phase 3: Electric Force — 3D Validation — Details

Port the validated 1D equations and force computation to M3 Taichi 3D engine. Verify that 3D results match 1D on-axis results, then test off-axis configurations for spherical symmetry. Validate force & motion integration — particles should move correctly under computed forces.

## Phase 4: Non-Linear Wave Equations (M3, 1D → 3D) — Details

If Phase 1b was triggered, this phase ports validated non-linear 1D equations to M3 3D. If Phase 1 succeeded linearly, this introduces non-linear equations fresh.

Variable λ(r) from two sources:

- **Yee & Hauger**: discrete wavelength shells r_n = 2(K-n)λ — changes interference pattern, non-uniform node spacing breaks sinc periodicity (WKB/eikonal phase integral approach)
- **Smoliński**: r⁵ energy scaling inside soliton's Energy Domain — defines how λ(r) varies near the wave center core

Variable ρ(x) from granule velocity — wave interference changes local cycling rate → local density/pressure. ∇ρ contributes to ∇E with different spatial structure than ∇A.

See [Phase 1b details](#phase-1b-non-linear-wave-equations-1d--details) for the theoretical foundation and [05_1D_wave.md](05_1D_wave.md#possible-solution-non-linear-wave-equations-phase-1b-fallback--phase-4) for full analysis.

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
