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

## [PHASE 1: FAR-FIELD Force (1D Sandbox only)](01_far_field.md)

- [ ] Resolve far-field oscillatory force (MAIN BLOCKER — sinc nodes in out-wave)
  - ✅ Implemented F = -∇E as standard force computation (replaces A·∇A chain rule expansion)
  - ✅ Tested gradient sampling radius / Gaussian smoothing (σ = 0.25λ to 2λ) — does not resolve; destroys charge signal
  - ✅ Tested smooth envelope interaction (|Z₁|·|Z₂| with imposed charge sign) — 17/17 direction + 1/r² scaling, but sign is not emergent
  - ✅ Ruled out numerical precision — 1D uses f64; M3/M4 f32 with sim-friendly units; oscillation is real math not artifact
  - ✅ Ruled out inertia filtering — phasor RMS already IS the time-averaged quantity; problem is spatial sinc, not temporal
  - ✅ Ruled out pressure/velocity gradient — 90° shift moves sinc zeros by λ/4 but preserves λ/2 oscillation period
  - ✅ Ruled out standing vs traveling wave decomposition — each component individually is smooth (1/kr), but coherent superposition of two sources still oscillates via cos(k·Δr) interference; oscillation is intrinsic to wave interference, not to standing/traveling character
  - ✅ Tested all 5 wave equations (Wolff, LaFreniere-Marcotte, Phase-warped, Combined, Weighted) — all produce force direction flips; confirms oscillation is intrinsic to coherent interference regardless of spatial function
- ✅ 1/r² force law scaling (RESOLVED) — interaction energy E_int ∝ |Z₁|·|Z₂| ∝ 1/r gives F ∝ 1/r² (confirmed numerically)
- [ ] Validate force direction: opposite charges attract, same charges repel (emergent, not imposed)
- [ ] Plot energy density landscape along axis at various separations
- [ ] Dual-treatment boundary: near-field raw phasor (lock-in) vs far-field smoothed envelope (Coulomb)

> **All linear scalar candidates exhausted (10/10 ruled out, including Phase 1a signed disturbance).** Three remaining paths (1b, 1c, 1d). Paths C and D are deeply connected — non-linear toroidal dynamics naturally produce vector patterns whose directional properties may carry charge information. They may converge into a single solution.

### ✅ [PHASE 1a: Signed Disturbance (forced charge sign)](01a_signed.md#phase-1a-signed-disturbance-forced-charge-sign) — RULED OUT

- ✅ Implemented signed disturbance model in 1D sandbox (equation #6): A₀ + q·δ(r) with `BASE_AMPLITUDE_RATIO`
- ✅ Tested δ(r) = 1/(1+kr) and 1/√(1+(kr)²) — smooth, 1/r far-field decay
- ✅ Same charge repulsion: 9/9 correct direction, near-constant Coulomb ratio
- ❌ Opposite charge attraction: 0/9 — asymmetric energy landscape, Newton's 3rd law violated
- ❌ **Charge sign is NOT emergent** — q = cos(phase) acts as a ±1 label on smooth potential, equivalent to previously ruled-out imposed-sign approach. Not genuine force emergence from wave interference

### ✅ [PHASE 1b: Base Wave + WC Energy Redistribution](01b_base_disturbance.md#phase-1b-base-wave--wc-energy-redistribution)

- Step 1 — Base wave modeling (`wave_engine_1D_v3.py`):
  - ✅ Implement 5 base wave candidate models: uniform, standing, stochastic, quadrature, laplacian
  - ✅ Validate energy uniformity: uniform ✓, quadrature ✓ (flat), standing ✗ (nodes at λ/2), stochastic ✓ (broadband fix)
  - ✅ Laplacian self-stabilizes to standing wave → **standing wave is the physically correct 1D base wave form** (Laplacian retired)
  - ✅ Stochastic monochromatic bug: `Σ cos(kx+φᵢ)` at same k collapses to single standing wave — fixed with broadband k spread
  - ✅ Quadrature: flat energy + traveling wave direction flips with temporal offset sign → possible charge/spin encoding via complex sinusoid channels

- Step 2 — WC disturbance and contender selection:
  - ✅ **Step 2a**: Node-locking charge hypothesis — FALSIFIED. Charge as spatial property (even/odd node position) does not predict force direction. 7/30 match (23%). Actual force has 2λ periodicity, not λ/2. Even separations produce net translation, not repulsion
  - ✅ **Step 2b**: Migrate WC disturbance from v2 to v3 — COMPLETED, additive model ruled out. Base wave + WC additive superposition produces same sinc oscillation as v2. Energy normalization (Option B) conserves ΣE but doesn't change spatial pattern → forces unchanged. WCs must warp the energy field non-additively (reflection, scattering, multiplicative)
  - ✅ **Step 2c**: Non-additive WC disturbance models — passive and elastic tested:
    - Passive (M2 re-validated in 1D):
      - ❌ Option A (multiplicative): 48/48 direction, energy conserved, but charge imposed ±1 (not emergent)
      - ❌ Option B (normalized additive): uniform scaling, no spatial change
      - ❌ Option C (scattering): sinc reintroduced via re-emitted wave, random direction
      - ❌ Option D (absorber): charge-blind, symmetric drain
    - Elastic (new territory, NOT in M2):
      - ❌ Option E (amplitude modulation): charge-blind — symmetric scaling, always repels 24/24
      - ❌ Option F (phase/λ warp): near-zero forces — rotation preserves RMS, no gradient. Requires variable-λ energy equation (Phase 1d) to produce force
      - ❌ Option G (L→T spin): **CHARGE SENSITIVE** — first model to distinguish charges. Opposite: 12/24 oscillates, same: 24/24 unclear. Quadrature proxy limited — needs true two-component displacement (Phase 1c)
  - ✅ **Step 2d**: Dual-channel base wave (π-apart) — dual_uniform: charge-blind (always repels, per-channel energy symmetric). dual_standing: partial (12/24 oscillates). Both: perfect energy conservation + Newton's 3rd. Root cause: `E_ch1 + E_ch2` is symmetric w.r.t. which channel is boosted — needs cross-channel coupling (like L→T) to break symmetry
  - ✅ **Step 2e**: Physics discussion and path decision — COMPLETED. 10 WC models tested: only L→T spin (Option G) distinguishes charges. 1D scalar sandbox has fundamental limitation (can't represent true L/T). Quadrature confirmed as strongest base wave (L/T duality). Path: Phase 1c first (vector displacement for L→T spin), then 1d (variable λ). Both converge into one solution

> **Phase 1b CONCLUSION**: the base wave exists (standing wave, physically validated). WCs must interact with it through **elastic disturbance** (changing wave character, not reflecting). The L→T spin conversion is the only charge-sensitive mechanism found — it needs true vector displacement (Phase 1c) and variable λ in the energy equation (Phase 1d) to fully work. Carry-forward tasks: energy redistribution, far-field drainage, force emergence, Coulomb validation — all require Phase 1d/1c capabilities.

### 🔶 [PHASE 1c: Vector Wave Force](01c_vector_wave.md#phase-1c-vector-wave-force)

> **From Phase 1b**: L→T spin conversion (Option G) is the ONLY charge-sensitive mechanism found (10 models tested). Quadrature phasor proxy showed charge discrimination but is limited — needs true independent L/T displacement.
>
> **Research strategy**: math-only numpy scripts (no visualization). Compute → sweep → analyze → document. Scripts in `scripts_vector_wave/`. When wave equations validated → port to M4 Taichi engine for 3D visualization.
>
> **Force mechanism**: `F = -∇E_total = -∇E_L - ∇E_T`. One force, two directions: ∇E_L → electric (longitudinal/radial), ∇E_T → magnetic (transverse/perpendicular). L/T defined relative to radial from WC. Transverse has 360° freedom → magnetic requires alignment/coherence to not cancel. Gravity = residual total energy deficit.

Step 1 — 3D Vector Base Wave (`step1_base_wave.py`):

- ✅ 3D grid with vector displacement: `ψ(r) = (ψ_x, ψ_y, ψ_z)` — 64³ grid, 8λ extent, 200 Fibonacci sphere sources
- ✅ Isotropic base wave from all directions → mean energy matches theory (ratio 1.003). Energy CV = 0.577 = 1/√3 (chi-squared(6) speckle — fundamental, not artifact)
- ✅ L/T decomposition: `A_L = |ψ · r̂|`, `A_T = |ψ - (ψ · r̂)r̂|` — E_L/E = 1/3, E_T/E = 2/3 (isotropic prediction). Holds at all reference points. E_L + E_T = E exact to machine precision
- ✅ Verify: base wave alone → force is speckle noise only (matches CV·E·k₀ estimate). No large-scale gradient. Null baseline confirmed

Step 2 — WC as L→T Converter (Spin) (`step2_single_wc.py`):

- ✅ L→T conversion at WC: spherical out-wave `sinc(kr)·exp(+i·phase)·[√(1-η)·r̂ + √η·q·(ẑ×r̂)]`. Sweep η from 0 to 1. Physical η = α ≈ 1/137 (fine structure constant)
- ✅ Energy concentration: 1.98x at WC core, returns to 1.0x beyond 1λ. Standing wave formation from in-wave (base) + out-wave (WC) interference
- ✅ L/T ratio shifts at WC core: E_L/E goes from baseline 0.33 up to 0.64 (η=0) or down to 0.17 (η=1). Shift is local (< 1λ)
- ✅ CW/CCW produce identical energy for single WC. Spin sign matters in Step 3 (two WCs)

Step 2a — Key Findings (Spin Scale & Sinc Resolution):

- ✅ **Sinc oscillation = correct K=1 physics**: neutrino is neutral, no spin, no Coulomb. Lock-in IS the strong force / particle formation mechanism. The "main blocker" was never a bug
- ✅ **Spin only at K≥10**: per EWT, single WC doesn't do L→T. Electron (K=10, tetrahedral) has spin from off-node WC repositioning
- ✅ **Annihilation from sinc**: opposite phase wells at r=0 (deepest), barriers at λ/2 (positronium). Same phase wells at λ/2 (lock-in)
- ✅ **No neutrino observational data**: all validation must be at K≥10 (Coulomb, annihilation, magnetic)
- ✅ **M3 electron unstable**: tetrahedral geometry has 15/45 pairs at non-node distances (√3×λ/2, √2×λ/2). Needs variable λ (Phase 1d) and/or vector forces (M4)
- ✅ **Jeff Yee + Dieter Hauger engaged**: Hauger (wavelength shells co-author) may have insights on standing→traveling transition

Step 3 — Two-WC Force Test:

- [ ] Two WCs with spin (η=α) at variable separation, compute `F = -∇(E_L + E_T)`
- [ ] Sweep: does force direction depend on spin sign (CW vs CCW)? Emergent from T component?
- [ ] Test if T component breaks sinc oscillation for far-field Coulomb

Step 4 — Coulomb Validation:

- [ ] Force magnitude vs Coulomb reference at multiple separations
- [ ] 1/r² scaling check
- [ ] Newton's 3rd law (equal and opposite)

Step 5 — Convergence with Phase 1d:

- [ ] Add variable λ(r) to energy: `E = ρV(c·A/λ(r))²`
- [ ] Test combined: vector displacement + variable λ
- [ ] Evaluate "one force, different directions" — electric (L), magnetic (T), gravitational (deficit)
- [ ] Connect to elliptical rotation handedness (CW = electron, CCW = positron)

### 🚧 [PHASE 1d: Non-Linear Wave Equations](01d_non_linear.md#phase-1d-non-linear-wave-equations)

> **From Phase 1b**: elastic phase warp (Option F) produces zero force because `E = ρV(fA)²` uses constant f — can't see λ variation. Phase 1d must implement `E = ρV(c·A/λ(r))²` where `∇λ` creates force from wavelength gradients. Converges with Phase 1c.

- [ ] Implement variable λ(r) in energy equation: `E = ρV(c·A/λ(r))²` — the `∇λ` force term
- [ ] Implement λ(r) profile from Yee & Hauger discrete wavelength shells, WKB phase integral
- [ ] Implement variable ρ(x) in 1D sandbox (density from granule velocity / wave interference)
- [ ] Test Smoliński Ψ³ cubic non-linearity (NLS soliton stabilizer — F(Ψ, ε_G, |ε_M|, N_ν))
- [ ] Test F = -∇E with spatially variable ρ(x), f(x), A(x) — all three gradients contributing
- [ ] Evaluate Smoliński Push-out Operator P̂Φ = -∇·(η_stat/η_soliton)∇Φ as variable-ρ force formalization
- [ ] Evaluate whether non-linear spatial structure breaks the sinc periodicity and resolves force direction
- [ ] Re-test elastic phase warp (Option F) with variable-λ energy equation
- [ ] If successful, validate against Coulomb reference (direction + 1/r² scaling)

## [PHASE 2: NEAR-FIELD Forces (1D Sandbox only)](02_near_field.md)

- [ ] Same-phase lock-in: verify oscillatory force creates stable energy wells
- [ ] Opposite-phase monotonic attraction: verify consistent attraction to annihilation
- [ ] Characterize near-field → far-field transition boundary
- [ ] Validate dual-treatment: raw phasor (near-field) vs smoothed envelope (far-field)
- [ ] Test with Verlet/leapfrog integrator (energy conservation for lock-in stability)
- [ ] Test with f64 precision (check if numerical drift causes escape from wells)

## [PHASE 3: Forces — 3D Validation (Taichi, port from 1D)](03_forces_3D.md)

> **Conditional scheduling**: Ports validated 1D results from  Phase 1c (vector) and/or Phase 1d (non-linear) to 3D engines. Phases 1c and 1d may converge here — non-linear toroidal dynamics naturally produce vector patterns that carry charge information.

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
