# ROADMAP

## [PHASE 0: Tooling](#phase-0-tooling--details)

- ✅ Build 1D wave sandbox (matplotlib, interactive controls)
- ✅ Phasor superposition (analytical amplitude, replaces EMA-RMS)
- ✅ Coulomb reference comparison in sandbox
- ✅ Decision: weighted partial standing wave as primary equation
- ✅ Build physics invariant tests (pytest, boundary limits, energy conservation)
- ✅ Validate phasor RMS: single WC sinc envelope, same/opposite charge interference, EMA-RMS equivalence
- ✅ Parameter sweep: force vs separation from 2λ to 10λ (`sweep_force_vs_separation.py`)
- ✅ Validate 1/r² force scaling — interaction energy E_int ∝ |Z₁|·|Z₂| ∝ 1/r gives F ∝ 1/r² (confirmed numerically)

## [PHASE 1: Electric Force — FAR-FIELD (1D Sandbox only)](01_electric_farfield.md)

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

> **All linear scalar candidates exhausted (9/9 ruled out). Signed disturbance model (Phase 1a) also ruled out — charge sign is not emergent from wave physics.** Three remaining paths. Paths C and D are deeply connected — non-linear toroidal dynamics naturally produce vector patterns whose directional properties may carry charge information. They may converge into a single solution.

### [PHASE 1a: Signed Disturbance (forced charge sign)](01_electric_farfield.md#phase-1a-signed-disturbance-forced-charge-sign--details) — RULED OUT

- ✅ Implemented signed disturbance model in 1D sandbox (equation #6): A₀ + q·δ(r) with `BASE_AMPLITUDE_RATIO`
- ✅ Tested δ(r) = 1/(1+kr) and 1/√(1+(kr)²) — smooth, 1/r far-field decay
- ✅ Same charge repulsion: 9/9 correct direction, near-constant Coulomb ratio
- ❌ Opposite charge attraction: 0/9 — asymmetric energy landscape, Newton's 3rd law violated
- ❌ **Charge sign is NOT emergent** — q = cos(phase) acts as a ±1 label on smooth potential, equivalent to previously ruled-out imposed-sign approach. Not genuine force emergence from wave interference

### [PHASE 1b: Base Wave + WC Energy Redistribution](01_electric_farfield.md#phase-1b-base-wave--wc-energy-redistribution--details)

- [ ] Model the base wave: isotropic standing wave field with known ρ, λ, A, f — granule displacement in multiple phases (not uniform cos(ωt)), uniform energy density
- [ ] Model WC disturbance: reflection/scattering of base wave → radial disturbance expansion, standing wave formation (nλ core, radius = K²λ), energy concentration near WC, energy deficit in far field
- [ ] Determine how WC phase affects the far-field energy drainage pattern (NOT via simple ±1 sign — the mechanism is unknown and must be discovered)
- [ ] Test force emergence: far-field drainage from WC1 disturbs WC2's standing waves → energy gradient at WC2 → force toward lower energy density
- [ ] Validate against Coulomb reference (direction + 1/r² scaling)

### [PHASE 1c: Non-Linear Wave Equations (1D)](01_electric_farfield.md#phase-1c-non-linear-wave-equations-1d--details)

- [ ] Implement variable λ(r) in 1D sandbox (Yee & Hauger discrete wavelength shells, WKB phase integral)
- [ ] Implement variable ρ(x) in 1D sandbox (density from granule velocity / wave interference)
- [ ] Test Smoliński Ψ³ cubic non-linearity (NLS soliton stabilizer)
- [ ] Test F = -∇E with spatially variable ρ(x), f(x), A(x) — all three gradients contributing
- [ ] Evaluate whether non-linear spatial structure breaks the sinc periodicity and resolves force direction
- [ ] If successful, validate against Coulomb reference (direction + 1/r² scaling)

### [PHASE 1d: Vector Wave Force (M4 displacement direction)](01_electric_farfield.md#phase-1d-vector-wave-force--details)

- [ ] Extend 1D sandbox with vector displacement (2D: x + y components)
- [ ] Compute divergence (∇·ψ), curl, or energy flux from vector field
- [ ] Test whether signed vector quantities (divergence, flux direction) recover charge-sign without oscillatory ambiguity
- [ ] If successful, validate against Coulomb reference (direction + 1/r² scaling)
- [ ] Connect to elliptical rotation handedness (M4 phasor: same-phase = same rotation, opposite = opposite)

## [PHASE 2: Electric Force — NEAR-FIELD (1D Sandbox only)](02_electric_nearfield.md)

- [ ] Same-phase lock-in: verify oscillatory force creates stable energy wells
- [ ] Opposite-phase monotonic attraction: verify consistent attraction to annihilation
- [ ] Characterize near-field → far-field transition boundary
- [ ] Validate dual-treatment: raw phasor (near-field) vs smoothed envelope (far-field)
- [ ] Test with Verlet/leapfrog integrator (energy conservation for lock-in stability)
- [ ] Test with f64 precision (check if numerical drift causes escape from wells)

## [PHASE 3: Electric Force — 3D Validation (Taichi, port from 1D)](03_electric_3D.md)

> **Conditional scheduling**: Ports validated 1D results from Phase 1c (non-linear) and/or Phase 1d (vector) to 3D engines. Phases 1c and 1d may converge here — non-linear toroidal dynamics naturally produce vector patterns that carry charge information.

M3 SCALAR:

- [ ] Port validated 1D equations to M3 3D wave engine
- [ ] Reproduce 1D force results in 3D (2 WCs on axis)
- [ ] Test off-axis configurations (verify spherical symmetry)
- [ ] Validate force & motion integration (particles move correctly)
- [ ] Compare against Coulomb force at multiple separations

M4 VECTOR:

- [ ] Port variable λ(r) to M3 3D engine (Yee & Hauger discrete wavelength shells)
- [ ] Port variable ρ(x) to M3 3D engine (density from granule velocity)
- [ ] Port Ψ³ cubic non-linearity (NLS soliton stabilizer) to 3D
- [ ] Test Smoliński r⁵ energy scaling near WC core
- [ ] Port vector force computation (divergence/curl/flux) to M4 3D engine
- [ ] Test elliptical rotation handedness as charge-sign indicator in M4
- [ ] Evaluate impact on near-field lock-in and far-field force scaling

## [PHASE 4: Magnetic Force (M4 Vector Waves)](04_magnetic_vector.md)

- [ ] Build M4 vector wave engine with transverse displacement
- [ ] Validate elliptical displacement trajectories (6 phasor numbers)
- [ ] Model spin as toroidal wave flow
- [ ] Demonstrate magnetic force from transverse wave interference
- [ ] Separate longitudinal (electric) and transverse (magnetic) force components

## [PHASE 5: Gravitational Force (M3, Multi-Particle)](05_gravitational.md)

- [ ] Test wave shading with particle clusters
- [ ] Test Smoliński buoyancy model: ρ(x) and f(x) as local variables
- [ ] Validate 10⁻⁴² EM-to-gravitational force ratio
- [ ] Validate computed G against Smoliński's Scilab reference values

## [PHASE 6: Time Dynamics (M4 or new method)](06_time_dynamics.md)

- [ ] Implement variable λ per voxel (local dt)
- [ ] Demonstrate time dilation from energy starvation mechanism
- [ ] Connect λ modulation → granule velocity → pressure → gravity
- [ ] Test force control via frequency/spin manipulation

## [PHASE 7: Emergent Waves](07_emergent_waves.md)

- [ ] Demonstrate photon-like traveling wave packets
- [ ] Test thermal energy as standing wave dynamics
- [ ] Validate electromagnetic wave emergence from medium disturbances

---

## PHASE 0: Tooling — Details

**1D wave sandbox** (`wave_engine_1D_v2.py`): The full 3D Taichi simulator (m3) is powerful but slow to iterate on. The 1D sandbox provides fast iteration with instant visual feedback — no GPU compile, no 3D rendering. Clean 1D profiles are directly comparable to LaFreniere's reference animations.

Features: configurable WC array (position, phase, amplitude), all 5 wave equation forms, real-time matplotlib animation (displacement + phasor RMS overlay, energy density, force field), interactive controls (separation slider, WC on/off toggles, phase offset toggle), Coulomb reference comparison with direction-match detection, force annotations at WC positions.

**Weighted partial standing wave** selected as primary equation:

```text
ψ = A · [w(r)·sin(kr + ωt) + sin(kr - ωt)] / kr
```

Why: standing waves near center (w ≈ 1), traveling waves far out (w → 0), interference between particles via traveling out-waves, phasor captures the envelope. Other forms have limitations: Wolff (no traveling), LaFreniere-Marcotte (transition not sharp enough), Phase-warped Marcotte (no true standing nodes near center).

**Physics invariant tests** (pytest): validate wave equations before deploying to 3D engine — dimensional consistency, boundary limits (r→0, r→∞), energy conservation, phasor/EMA-RMS equivalence, near-field/far-field continuity, charge symmetry.

**Phasor validation** (prerequisite to force work): confirm correct amplitude patterns — single particle sinc envelope, same-charge destructive interference, opposite-charge constructive interference, EMA-RMS equivalence.
