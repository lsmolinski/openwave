# ROADMAP

## ✅ [PHASE 0: TOLLING & CONTEXT](#phase-0-tooling--details)

- ✅ Build 1D wave sandbox (matplotlib, interactive controls)
- ✅ Phasor superposition (analytical amplitude, replaces EMA-RMS)
- ✅ Coulomb reference comparison in sandbox
- ✅ Decision: weighted partial standing wave as primary equation
- ✅ Build physics invariant tests (pytest, boundary limits, energy conservation)
- ✅ Validate phasor RMS: single WC sinc envelope, same/opposite charge interference, EMA-RMS equivalence
- ✅ Parameter sweep: force vs separation from 2λ to 10λ (`sweep_force_vs_separation.py`)
- ✅ Validate 1/r² force scaling — interaction energy E_int ∝ |Z₁|·|Z₂| ∝ 1/r gives F ∝ 1/r² (confirmed numerically)

## ✅ [PHASE 1: EXPLORING WAVE EQUATIONS](1_EXPLORING.md)

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

### ✅ [Phase 1a: Signed Disturbance (forced charge sign)](1a_signed.md#phase-1a-signed-disturbance-forced-charge-sign) — RULED OUT

- ✅ Implemented signed disturbance model in 1D sandbox (equation #6): A₀ + q·δ(r) with `BASE_AMPLITUDE_RATIO`
- ✅ Tested δ(r) = 1/(1+kr) and 1/√(1+(kr)²) — smooth, 1/r far-field decay
- ✅ Same charge repulsion: 9/9 correct direction, near-constant Coulomb ratio
- ❌ Opposite charge attraction: 0/9 — asymmetric energy landscape, Newton's 3rd law violated
- ❌ **Charge sign is NOT emergent** — q = cos(phase) acts as a ±1 label on smooth potential, equivalent to previously ruled-out imposed-sign approach. Not genuine force emergence from wave interference

### ✅ [Phase 1b: Base Wave + WC Energy Redistribution](1b_base_disturbance.md#phase-1b-base-wave--wc-energy-redistribution)

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

### ✅ [Phase 1c: Vector Wave Force](1c_vector_wave.md#phase-1c-vector-wave-force)

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

Step 3 — Two-WC Force Test (`step3_two_wc_force.py`):

- ❌ Tested 4 spin configs (CW-CW, CCW-CCW, CW-CCW, CCW-CW) at η=α, 0.1, 0.5 across separations 2-6λ
- ❌ ALL configs show MIXED force directions — no consistent charge-dependent pattern
- ❌ Root cause: on-axis limitation. T component is perpendicular to WC axis → creates magnetic force, not electric. L component still has sinc oscillation → no Coulomb
- ✅ **Conclusion**: spin creates magnetic force (Phase 4), NOT electric force. Electric force needs variable λ(r) (Phase 1d)
- ✅ **LaFreniere phase shift discovery**: electron core (1λ diameter) creates λ/2 phase shift from medium compression (7x). This phase shift — not spin — is what creates charge sign. Connects directly to Phase 1d variable λ(r)

Steps 4-5 — merged into Phase 1d:

- [ ] Phase 1c vector infrastructure (Steps 1-2) ready for Phase 1d integration
- [ ] Variable λ(r) from Yee & Hauger shells + LaFreniere core phase shift → electric force
- [ ] Combined vector displacement + variable λ → full force decomposition (electric from ∇λ, magnetic from spin T, gravitational from deficit)

### ✅ [Phase 1d: Non-Linear Wave Equations](1d_non_linear.md#phase-1d-non-linear-wave-equations)

> **From Phase 1c**: spin → magnetic (not electric). Variable λ(r) → electric force candidate. Converges with Phase 1c vector infrastructure.

1D Variable-λ Test (`step1d_variable_lambda.py`):

- ✅ Implement variable λ(r) from Yee & Hauger shells + WKB phase integral
- ✅ Implement variable-λ energy equation: `E = ρV(c·A/λ(r))²`
- ✅ Test K=1 (neutrino): no λ variation → sinc persists → neutral. Correct.
- ✅ Test K=10 (electron): ∇λ term IS active (reverses force vs const-λ at some separations)
- ❌ Charge-blind: ∇λ creates same force for all phase configs (λ depends on K, not charge)
- ⚠️ **1D on-axis limitation**: sinc oscillation dominates on-axis. But 3D spherical integration may average out the sinc through different oscillation periods at different angles

Sinc Elimination Tests:

- ❌ **#5 Statistical averaging** (step1d_averaged_force.py): sinc perfectly symmetric → averaging gives 50/50, not a direction
- ❌ **#1 Broadband shells** (step1d_broadband.py): multiple cosines still oscillate, no convergence
- ⚠️ **#4 1D Flux** (step1d_flux_force.py): same-charge K=10 = 50/50 REP (CONSISTENT, first ever!), but opp-charge also REP (baseline push dominates in 1D)
- ⚠️ **#4 2D Flux** (step1d_flux_force_2d.py): **100% charge discrimination** — same ≠ opposite at ALL 22 separations. Sinc persists (direction flips every λ/2) but charge sensitivity is perfect. Coulomb correct at half-integer λ separations.

> **Phase 1 CONCLUSION**: the sinc oscillation cos(k·Δr) is intrinsic to monochromatic spherical wave interference. It cannot be removed by variable λ, spin, decomposition, averaging, or broadband. However, the 2D radiation pressure (flux) produces 100% charge discrimination — same and opposite ALWAYS get opposite force directions. The sinc determines WHERE the equilibrium/lock-in points are; the charge difference determines WHICH direction at each point.

Carry-Over to Phase 2:

- [ ] Investigate flux-based force (`S = -c²·ψ·∇ψ`) as the primary Coulomb mechanism instead of `F = -∇E`
- [ ] Test 3D flux integration — does full spherical flux produce consistent Coulomb at any K?
- [ ] Explore whether sinc lock-in + charge discrimination IS the unified force (near-field = strong, far-field = Coulomb envelope)
- [ ] Variable ρ(x), Ψ³ non-linearity, combined vector + variable λ (deferred from Phase 1d)

## 🔶 [PHASE 2: ENERGY LAYERS](2_ENERGY_LAYERS.md)

Follows the energy layers hierarchy from [0_OVERVIEW.md](0_OVERVIEW.md). See [2_ENERGY_LAYERS.md](2_ENERGY_LAYERS.md) for strategy, context sources, and open questions.

### 🔶 LAYERS 1-3: [PARTICLE EMERGENCE](2L3_particle_emergence.md) from Near-Field Standing Waves (M3 Scalar)

M3 scalar method handles longitudinal standing wave physics. Goal: demonstrate stable standalone particle formation. No far-field analysis on M3 — Coulomb needs spin → vector field (M4).

#### Layer 1: Fundamental energy wave

- [ ] Base wave

#### Layer 2: WC disturbance → standing wave (fundamental particle)

- [ ] Validate same-phase lock-in in 3D animation
  - [ ] Two fundamental particles (WCs), same phase, observe oscillation in energy wells
  - [ ] Measure well depth, oscillation period, escape energy
- [ ] Validate opposite-phase annihilation pathway
  - [ ] Two fundamental particles, opposite phase, observe attraction through barriers
  - [ ] Measure barrier heights at λ/2 intervals, compare with positronium lifetime data

#### Layer 3: K=1 lock-in → standalone particle formation

- [ ] Test multi-WC lock-in: K=2, K=3, ..., K=10 progressive build-up
  - [ ] K=2..9 should be unstable (temporary particles that decay)
  - [ ] K=10 should be the first fully stable standalone particle (electron)
  - [ ] Test geometry: self-organize vs prescribed?
  - [ ] Validates EWT prediction: 1-3-6 tetrahedron is simplest 3D geometry where all WCs sit near nodes
- [ ] Stabilize K=10 electron tetrahedron on M3
  - [ ] Test other wave equations
  - [ ] Implement variable λ(r) from Yee & Hauger shells (non-uniform node spacing)
  - [ ] Test whether variable nodes accommodate the 15/45 non-node pairs (√3×λ/2, √2×λ/2)
  - [ ] Leapfrog integrator + damping already in place
- [ ] Characterize near-field → far-field transition boundary (particle radius K²λ)
  - [ ] Where do standing waves end and traveling waves begin?

#### Wrap-Up Layers 1-3

- [ ] Collect proof for particle formation from waves (electron properties)
- ✅ Rename xperiments with near-field results
- [ ] review <https://energywavetheory.com/project/phase1/>
- [ ] review <https://energywavetheory.com/project/phase2/>

### 🚧 LAYER 4: [ELECTROMAGNETISM EMERGENCE](02L4_electromagnetism_emergence.md) from Far-Field Traveling Waves (maybe M3 Scalar enough, or M4 Vector)

M4 vector method may be needed for spin (L→T conversion) → charge, Coulomb, magnetic force. Non-dual standalone particles (electron family) have spin; dual particles (neutrino family) are neutral. See [2_ENERGY_LAYERS.md](2_ENERGY_LAYERS.md) for context sources and references to load when starting this phase.

#### Layer 4a-b: Charge emergence (from non-dual spin)

- [ ] Investigate charge mechanism
  - [ ] Spin
  - [ ] How does spin + 1-3-6 non-dual geometry create the wavefront disturbance?
  - [ ] Role of LaFreniere core phase shift (λ/2 from compression)
  - [ ] Does variable λ(r) inside the standalone particle create the charge sign?

#### Layer 4c: Coulomb / electric force (1st measurable force, longitudinal)

- [ ] Build / extend M4 vector wave engine with L→T spin conversion
  - [ ] Vector displacement `ψ = (ψ_x, ψ_y, ψ_z)` per voxel
  - [ ] L→T conversion at standalone particle WCs (η = α ≈ 1/137)
  - [ ] Elliptical displacement trajectories (6 phasor numbers)
- [ ] Test standalone particle (K=10 electron) as Coulomb source on M4
  - [ ] Place one stabilized electron, measure far-field force on a test particle
  - [ ] Does spinning tetrahedral geometry produce consistent Coulomb direction?
  - [ ] Compare with fundamental particle / neutrino (should be neutral — no Coulomb)
- [ ] Investigate 3D flux-based force
  - [ ] Extend Phase 1 2D flux test to full 3D spherical integration
  - [ ] Does solid-angle averaging smooth the sinc absolute direction?
  - [ ] Test at K=10 standalone particle scale (100λ radius)
- [ ] Validate Coulomb
  - [ ] Direction: same charge repels, opposite attracts at ALL separations
  - [ ] Magnitude: 1/r² scaling
  - [ ] Symmetry: Newton's 3rd law (equal and opposite)
- [ ] review <https://energywavetheory.com/project/phase2/>

#### Layer 4d: Bohr magneton / magnetic force (2nd force, transverse)

- [ ] Model spin as toroidal wave flow (Smoliński Energy Domain, Butto vortex)
- [ ] Separate longitudinal (electric) and transverse (magnetic) force components
  - [ ] Electric: ∇E_L from traveling longitudinal waves beyond particle radius
  - [ ] Magnetic: ∇E_T from transverse wave generated by spin
- [ ] Demonstrate magnetic force from transverse wave interference
- [ ] review <https://energywavetheory.com/project/phase2/>

### 🚧 LAYER 5: [GRAVITY EMERGENCE](02L5_gravity_emergence.md)

- [ ] Gravitational force from spin energy deficit (L→T drainage → amplitude deficit)
  - [ ] Test whether deficit accumulates over K WCs (factor of K×α?)
- [ ] Test wave shading with standalone particle clusters → accumulated deficit → gravity
- [ ] Compare with 10⁻⁴² EM-to-gravitational ratio
- [ ] Validate against Smoliński's Scilab reference values

### 🚧 LAYER 6: [EMERGENT WAVES](2L6_emergent_waves.md)

#### Layer 6a: Electromagnetics Waves

- [ ] Demonstrate photon-like traveling wave packets
- [ ] Validate electromagnetic wave emergence from medium disturbances

#### Layer 6b: Thermal Waves

- [ ] Test thermal energy as standing wave dynamics
- [ ] Validate thermal wave emergence from medium disturbances

### 🚧 LAYER 7: [COMPOSITE PARTICLE](02L7_composite_particle.md) → nuclei → atoms (Orbital Force)

#### Layer 7a: Strong Force (electric at sub-λ between K=10 standalone particles)

- [ ] Test K=10 near-field overlap → combined standing waves, high-energy lock-in
- [ ] Verify ~137× Coulomb strength (= 1/α) at sub-λ distances
- [ ] Investigate gluonic field mechanism (magnetic role?)

#### Layer 7b: Composite Particles

- [ ] Composite particle formation
  - [ ] Proton: 4 electrons + 1 positron at center (tetrahedral)
  - [ ] Neutron: proton + electron at center (neutralizes charge via destructive interference)
  - [ ] Test whether composite particles self-assemble from standalone particles
- [ ] review <https://energywavetheory.com/project/phase3/>
- [ ] review <https://energywavetheory.com/project/phase4/>
- [ ] review <https://energywavetheory.com/project/phase5/>

## 🔶 PHASE 3: LAGRANGIAN FRAMEWORK — see [3_LAGRANGIAN_FRAMEWORK.md](3_LAGRANGIAN_FRAMEWORK.md) · [3a_lagrangian_experiments.md](3a_lagrangian_experiments.md) · [3b_concept_review.md](3b_concept_review.md) · [3c_topological_defect.md](3c_topological_defect.md) · [3d_path_to_m5.md](3d_path_to_m5.md)

A research thread evaluating whether a Lagrangian / topological framework can replace OpenWave's empirical wave-equation search with a first-principles derivation, and produce charge quantization + far-field Coulomb that the M3 scalar method cannot. Sparked by email exchange with Jarek Duda (Jagiellonian) and Robert Close (Clark College) in the "Models of Particles" group.

- ✅ 8 sandbox numerical experiments — validate the core Lagrangian / topological physics before any production engine refactor (COMPLETE 2026-04-17)
  - ✅ Test 1: Sine-Gordon 1D — kink solitons, pair creation/annihilation, Lorentz contraction (build intuition) — passed 2026-04-16
  - ✅ Test 2: Hedgehog energy vs distance — clean 1/d Coulomb, R²=0.993, no sinc — passed 2026-04-16
  - ✅ Test 3: Topological charge quantization — Q=±1 integer, surface-independent, stable under 50% noise — passed 2026-04-16
  - ✅ Test 4: Klein-Gordon from twist dynamics — dispersion ω² = c²k² + m² validated to R² = 0.999982 across 9 modes; slope c² within 0.05%, mass gap within 1.3% — 2026-04-17
  - ⚠️ Test 5: Lagrangian derivation — Smolinski Ψ³ + Noether energy confirmed ✅; M4 sum-form Combined W-L satisfies free-wave EL ✅; doc product form `sin(kr/2)·cos(kr/2−(ωt+φ))/r` does NOT (quadrature term is not a free-wave solution — modeling choice, needs source term to exist) — 2026-04-17
  - ⚠️ Test 6: Three lepton families — E(K) scaling validated (R²=1.0). Three distinct K values reproduce lepton mass² ratios by construction (consistency check, not derivation). Full biaxial Q-tensor derivation deferred to Exp 6.1 — 2026-04-17
  - ⚠️ Test 7: Close's nonlinear vector wave eq — Close's ACTUAL Eq. 19 (linear) and Eq. 21 ("Equation of Everything") correctly implemented after obtaining the paper. Y_l^m seeds disperse (consistent with Close's framework — particles are plane-wave bispinors, not static solitons). Close's equation is a candidate base wave dynamics layer for M5 — 2026-04-17
  - ❌ Test 8: Smolinski's non-linear Ψ³ — K-selectivity hypothesis **FALSIFIED** at Level 1. Ψ³ produces breathing oscillation but no K-dependent geometric stabilization. Nonlinearity alone is insufficient; topology (Exp 2/3) is required for K-selectivity — 2026-04-17
- ✅ Winning combination selected (2026-04-17): topology (hedgehog + winding) + Klein-Gordon wave dynamics + Close's Eq. 19 as base vector wave + M3 near-field standing waves + Skyrme stabilizer. Full recipe in [3a § Winning Approach](3a_lagrangian_experiments.md#winning-approach-for-m5)
- 🔶 M5 / LAGRANGIAN-FIELD METHOD implementation (dir `openwave/xperiments/m5_lagrangian_field/`) — see [3d_path_to_m5.md](3d_path_to_m5.md) for full sub-phase plan
  - ✅ **M5.0 Scaffold** — all 11 sub-phases complete as of 2026-05-08
    - ✅ M5.0a Module rename + alias (`wave_engine.py` → `lagrangian_engine.py`, `ewave` → `lagrange`)
    - ✅ M5.0b Triple buffer (`psi_prev_am` / `psi_am` / `psi_new_am`) + AMR-ready field-storage abstraction
    - ✅ M5.0c Vector Laplacian via 6-point stencil (port + simplify from M2)
    - ✅ M5.0d.1 Leapfrog kernel `evolve_psi` + standing-wave eigenmode test
    - ✅ M5.0d.2 CFL evaluation + plane-wave seed (Gaussian-windowed packet) + tracker EMA + Hamiltonian dashboard + delete legacy M4 `propagate_wave`
    - ✅ M5.0d.3 Drop `scale_factor` / `ewave_res` / EWT-default constants; introduce xperiment-driven `wave_res`
    - ✅ M5.0e Curl, divergence, curl-curl operators (`∇×(∇×ψ) = ∇(∇·ψ) − ∇²ψ` identity form, 2-cell halo) — analytical checks pass on linear / rigid-rotation / Gaussian fields
    - ✅ M5.0f Storage-units decision + natural-units deferral (decision-record sub-phase, no code refactor): storage stays `_am` / `_rs` / `_rHz` for M2/M3/M4 consistency; kernel-internal natural-unit scaling deferred to M5.2 where it lands alongside the nonlinear physics (Klein-Gordon, Close Eq. 23, LdG) that benefits from textbook-readable couplings
    - ✅ M5.0g Per-voxel energy density (Hamiltonian) + force-computation switch (F = −∇E replaces M4's postulated `E = ρV(fA)²` formula); naming convention adopted: quantity is *energy*, `_H` suffix tags the Hamiltonian formula (parallels `_L` Lagrangian, `_K` kinetic in future); V_psi hook in place for M5.2 nonlinear couplings
    - ✅ M5.0h Dispersion gating test — leapfrog reproduces full discrete dispersion (`sin(ω·dt/2)=(c·dt/2)·√K`) within ±0.5% c² recovery across 5 modes (voxels/λ_x ∈ {31,21,16,10,8}). KG mass-term flavor deferred to M5.2 where the mass term lands. Two persisted lessons (Taichi Metal auto-reduce hits atomic contention; dispersion fits must use full space+time relation) captured in `lagrangian_engine.py` AUTO-REDUCE CAVEAT block + auto-memory feedback entries
    - ✅ M5.0i Performance profile baseline (no Tier 2 opts landed; not justified by current budget). 384³ runs at 51 fps (well over 20 fps target). Two persisted lessons during the investigation: (a) rotating-pointer `swap_buffers` is NOT the 30-min change it looks like — Taichi `ti.template()` caches attribute lookups, so attribute rotation is invisible to cached kernels; proper fix requires passing fields explicitly (2-3 hr refactor). (b) **Fusion is the biggest Tier 2 lever** (~50 % win, would close most of the 2.1× gap vs M2); proven prior art is M2's `propagate_wave` in `m2_laplace_propagation/wave_engine.py:603` which fuses leapfrog + tracker EMA + zero-crossing freq + buffer swap into a single ndrange. Deferred to M5.2 re-profile when V(ψ) gets real body
  - ✅ M5.1 Port topology from Exps 2, 3 — **8 of 8 tasks complete** as of 2026-05-11: ✅ `seed_vacuum`, ✅ `seed_hedgehog` (N-defect), ✅ director-glyph visualizer (3-plane line glyphs, palette `1−n_z`, blueprint+ironbow, half-arrowhead barbs, VIZ_STRIDE consolidated), ✅ Frank elastic `compute_energyF_density` + WAVE_MENU=5 (new `FieldObservables` class), ✅ gradient-descent `relax_director_step` (auto-relax on seed), ✅ **1/d Coulomb gating test (M5.2 gate)** — R² = 0.978 attractive + monotone (threshold 0.95; Dirichlet BC limits vs periodic) + **visual confirmation** in [`research_hub/3f_coulomb_visual_geometry.md`](research_hub/3f_coulomb_visual_geometry.md) (dumbbell F-density bridge for opposite charges; pinched perpendicular splay for same charges — matches classical EM field-line geometry, derived from pure topology), ✅ winding-number tracker `compute_winding_number` (Q=±0.996 on seeded ±1 hedgehogs). **Bonus discovery during task 8**: same-direction-pin + Dirichlet-BC setup loses topology after ~100 relax steps — Coulomb test still works (blend-zone elastic dominates F(d)), but defect stability needs M5.5 Skyrme to fully fix. Winding tracker now serves as the diagnostic gate for that issue during M5.2. **M5.1 effectively complete; M5.2 unblocked.**
  - 🔶 M5.2 Wave dynamics — **Close's Eq. 23** as the particle equation (preserves `∇·s = 0`, per Close's 2026-04-18 guidance) + Eq. 19 as V=0 linear limit + Klein-Gordon mass term; validate against Exp 4 dispersion; resonance-hunt amplitude sweep. Kernel-internal natural-unit scaling (`c=1, λ_C=1, ℏ=1`) lands here for the nonlinear couplings (deferred from M5.0f — linear kernels are dimensionally self-balancing and don't need it). **Progress (2026-05-12)**: ✅ Step 1 natural-unit scaffold, ✅ Step 2 KG mass term (`V = ½m²|ψ|²`), ⚠️ Step 3 defect-survival NEGATIVE under plain KG, ⚠️ Step 4a Mexican-hat φ⁴ PARTIAL (`V += ¼λ(|ψ|²−1)²` — damps `|ψ|` excursions to unit sphere but does NOT preserve Q; texture coherence lost in 2-5 propagation steps regardless of pre-relax). External-comms trigger NOT met; 🚧 Step 4b Skyrme term `+ ½κ|∇ψ × ∇ψ|²` next (textbook 1962 fix for 3D hedgehog Derrick collapse — brings forward M5.5 work)
  - [ ] M5.3 Hamiltonian energy density `H = ½ψ̇² + ½c²(∇ψ)² + V(ψ)` replaces postulated `E = ρV(fA)²`
  - [ ] M5.4 **Headline test**: single biaxial hedgehog (electron) is a long-lived resonance; hedgehog + anti-hedgehog pair reproduces dynamic 1/d² Coulomb + annihilation. Single-defect framing per Dr. Duda's lepton-axis hierarchy (electron ≠ K=10 tetrahedron in M5). Additional success criterion: measure the **constant magnetic dipole moment μ** of the relaxed hedgehog — this is the direct answer to Duda's 2026-05-11 graphene-thread question ("how is the electron's constant μ obtained in your model?"), where his own framing is hairy-ball-theorem + clock-propulsion (the latter then probed by M5.7-M5.8)
  - [ ] M5.5 Skyrme stabilizer (conditional on M5.4 showing defect collapse under Derrick's theorem)
  - [ ] M5.6 Biaxial LdG Q-tensor with `(δ, 1, g)` axis hierarchy (long-term, for lepton mass derivation; Exp 6.1 equivalent)
  - [ ] M5.7 Cornell potential / quark confinement (topological vortex string, `V(r) = −α/r + σ·r` with σ ≈ 1 GeV/fm) — per Dr. Duda's 2026-04-17 guidance
  - [ ] M5.8 De Broglie clock / Zitterbewegung test (`ω = 2mc²/ℏ` for electron + neutrino) — per Dr. Duda's 2026-04-17 guidance
- [ ] **Revisit Energy Layers hierarchy** with Lagrangian perspective — insert **Layer 0: Vacuum at rest** before Layer 1; rewrite Layer 1-3 semantics in [2_ENERGY_LAYERS.md](2_ENERGY_LAYERS.md), [0_OVERVIEW.md](0_OVERVIEW.md), and [3z_review_energy_layers.md](3z_review_energy_layers.md) (the content was moved out of the README pending re-evaluation; once revised, decide whether to put it back) — rationale in [3b_concept_review.md](3b_concept_review.md)

## [PHASE 4: TIME DYNAMICS](5_TIME_DYNAMICS.md)

- [ ] Implement variable λ per voxel (local dt)
- [ ] Demonstrate time dilation from energy starvation mechanism
- [ ] Connect λ modulation → granule velocity → pressure → gravity
- [ ] Test force control via frequency/spin manipulation

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
