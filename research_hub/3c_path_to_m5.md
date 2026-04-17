# PATH TO M5 — LAGRANGIAN-WAVE METHOD

Implementation plan for **M5 / LAGRANGIAN-WAVE METHOD** (directory `openwave/xperiments/m5_lagrangian_wave/`), the production wave engine that graduates sandbox-validated Lagrangian / topological physics onto the GPU-accelerated OpenWave platform. This document references concrete code in the current engines and defines what M5 inherits, replaces, and adds.

**Naming**: the method name follows the existing OpenWave convention — *adjective + Wave + Method* (parallels M4 *Vector-Wave Method*). "Lagrangian" refers to the formalism (L = T − V, Euler-Lagrange equations, action principle) used to derive the wave equation being solved; "Wave" because it remains a wave engine; "Method" because this is one specific simulator implementation within the broader Lagrangian framework.

**Spec inputs**:

- [3_LAGRANGIAN_FRAMEWORK.md](3_LAGRANGIAN_FRAMEWORK.md) — Lagrangian framework evaluation, 8 experiments, Duda/Close context
- [3b_lagrangian_experiments.md](3b_lagrangian_experiments.md) — sandbox numerical results (experiment-by-experiment)
- [0_WAVE_EQUATION.md](0_WAVE_EQUATION.md) — the M2/M3 vs Lagrangian comparison and why the equation is the *consequence*, not the goal

**Production code references**:

- `openwave/xperiments/m2_laplace_propagation/wave_engine.py` — existing PDE-stepping engine (free wave, scalar)
- `openwave/xperiments/m4_vector_wave/wave_engine.py` — existing analytical-superposition engine (vector, WPSW)

---

## WHY M5 IS A NEW ENGINE, NOT AN M4 UPGRADE

M4's architecture is **analytical superposition**:

- Each voxel evaluates the closed-form Weighted Partial Standing Wave directly (`m4_vector_wave/wave_engine.py:197–211`)
- Total field = linear sum over active WCs (`wave_engine.py:209`, `+= oscillator * direction`)
- Each WC is a query parameter (position + phase offset) that every voxel reads each frame
- Energy is postulated: `E = ρ·V·(f·A)²` from EMA on RMS amplitude (`wave_engine.py:287–293`)

This pattern has a **fundamental incompatibility** with the Lagrangian framework:

> Linear superposition `ψ = Σ ψ_i` only works when the wave equation is linear. Any nonlinear potential V(ψ) (Smolinski `−k·ψ³`, Sine-Gordon `sin(φ)`, Landau-de Gennes, Close's vector equation) destroys superposition. The field must be evolved *as a whole*, not composed from per-WC contributions.

Further, topology requires the field to carry **winding information**. In M4, the field at each voxel is reconstructed fresh every frame from WC positions — the field has no memory of its own topological state. A topological defect (hedgehog, kink, knot) is a property of the *field configuration*, which M4 doesn't persist.

M5 therefore replaces the core physics kernel while keeping the production infrastructure (Taichi, grid, visualization, force motion, GUI).

---

## WHAT M5 INHERITS FROM M2

M2 (`m2_laplace_propagation/wave_engine.py`) already solves the **free wave equation** via leapfrog time-stepping on a 3D grid. It is, in essence, a linear-Lagrangian PDE engine — the simplest case of what M5 becomes. The following patterns port directly to M5:

| M2 element | Location | Reuse in M5 |
| --- | --- | --- |
| 6-point discrete Laplacian | `compute_laplacianL`, `wave_engine.py:527–562` | Same stencil, computes `∇²ψ` per component (x, y, z in vector case) |
| Triple-buffer pattern (`psi_old`, `psi`, `psi_new`) | Field allocation in `WaveField` + update at `wave_engine.py:649–660` + swap at line 743 | Same buffering needed for leapfrog; extended to vector component by component |
| Leapfrog update formula `ψ_new = 2ψ − ψ_old + (c·dt)²·∇²ψ` | `wave_engine.py:649–653` | Replace RHS with `(c·dt)²·∇²ψ − dt²·dV/dψ` (add nonlinear potential term) |
| Dirichlet BC via skipped boundary voxels | `wave_engine.py:641` (`ndrange((1, nx-1), …)`) | Same approach for clamping ψ=0 at edges, or absorbing layer (commented block at lines 662–699 is a partial PML start) |
| Zero-crossing frequency detection with EMA | `wave_engine.py:724–741` | Transfers directly for per-voxel frequency tracking |
| EMA on ψ² for RMS amplitude | `wave_engine.py:711–722` | Same tracker — one per vector component or per scalar field |
| `charge_*` initial-condition kernels (Gaussian, 1/r falloff, radial sinusoidal, 1λ bubble) | `wave_engine.py:30–257` | Pattern is directly applicable, but **renamed in M5 to `seed_*`** (see note below) — kernels like `seed_hedgehog`, `seed_sine_gordon_kink` configure the initial topological state instead of a spherical energy pulse |
| CFL-aware timestep via `c_amrs * dt_rs` | `wave_engine.py:606–607, 652` | Same dt budgeting |

**Key conceptual point**: M2 is to a linear Lagrangian what M5 is to a full Lagrangian. The only physics difference is the absence of the `V(ψ) ≠ 0` term. M5 is "M2 with teeth" — nonlinear potential, vector field, topological ICs.

---

## WHAT M5 INHERITS FROM M4

M4's vector-field infrastructure and production polish transfer wholesale:

| M4 element | Location | Reuse in M5 |
| --- | --- | --- |
| Vector-valued displacement (3 components per voxel) | `displacement_am[i,j,k]` as `ti.Vector([0,0,0])`, `wave_engine.py:152` | Same — but M5 evolves each component via independent PDE stepping (no longer constrained to `* direction` radial projection at line 211) |
| Flux-mesh visualization (XY/XZ/YZ planes, colormaps) | `update_flux_mesh_values`, `wave_engine.py:544–703` | Unchanged — visualizes any scalar tracker regardless of source |
| Granule position rendering from displacement | `sample_position_to_render`, `wave_engine.py:379–403` | Unchanged |
| 3-plane sampling for global averages | `sample_avg_trackers`, `wave_engine.py:410–536` | Unchanged — computes RMS / freq / energy averages |
| WC position + phase offset + active flag tracking | `WaveCenter` data class (referenced at `wave_engine.py:52–58`) | **Role flips**: WC is no longer a query parameter but a *tracked measurement* of where a topological defect sits in the field |
| Force = −∇E via neighbor sampling | `select_voxels`, `wave_engine.py:31–101` (and downstream in `force_motion.py`) | Keep the gradient computation; the *source* of E changes (now Hamiltonian density, not postulated `ρV(fA)²`) |
| GUI, xperiments scaffold, scale_factor, constants | throughout | Unchanged |

**What M4 gives up**: the analytical WPSW kernel and the selective-voxel shortcut (`select_voxels`, `propagate_wave_neighbors`). The PDE is spatially coupled (neighbors via Laplacian), so voxel skipping is no longer valid — M5 runs full grid, same as M2 does today.

---

## WHAT M5 ADDS (NEW PHYSICS)

Everything below is net-new infrastructure, gated on sandbox experiment results.

### 1. Nonlinear potential term in the EoM

The M5 leapfrog step becomes:

```text
ψ_new = 2·ψ − ψ_old + dt²·[c²·∇²ψ − ∂V/∂ψ]
```

The specific `V(ψ)` is selected from sandbox results:

| Candidate | Potential term | Status |
| --- | --- | --- |
| Sine-Gordon family | `V(φ) = (m²c⁴/ℏ²)(1 − cos φ)` → `∂V/∂φ = (m²c⁴/ℏ²) sin(φ)` | ✅ validated (Exp 1) — useful for 1D analogs, not 3D hedgehogs |
| Smolinski Ψ³ | `V(ψ) = (k/4)·ψ⁴` → `∂V/∂ψ = k·ψ³` | ❌ falsified for K-selectivity (Exp 8); usable only as in-configuration stabilizer |
| Landau-de Gennes | `V(M) = a·Tr(M²) − b·Tr(M³) + c·(Tr M²)²` | ⚠️ mechanism validated (Exp 6); specific parameters deferred to M5.6 |
| Close's elastic solid | Eq. 19 `∂²Q = −c²·∇×∇×Q` + optional Eq. 21 nonlinear terms `−u·∇s + w×s` | ✅ validated (Exp 7 v2) — **selected as M5's base wave dynamics layer** |
| Klein-Gordon mass term | `V(ψ) = ½m²·ψ²` → `∂V/∂ψ = m²·ψ` | ✅ validated (Exp 4) — **selected as M5's perturbation-mass mechanism** |
| Mixed / composite | linear combination (topology seeding + Close dynamics + KG mass + optional Skyrme) | **adopted** — full recipe in [3b § Winning Approach](3b_lagrangian_experiments.md#winning-approach-for-m5) |

### 2. True independent vector components (no radial constraint)

M4 collapses its vector to radial (`oscillator * direction`, line 211). M5 evolves each of `(ψ_x, ψ_y, ψ_z)` independently via its own PDE step. No radial projection; topology emerges from the full vector field.

**Storage**: either three coupled scalar fields, or one `ti.Vector.field(3, dtype=ti.f32, shape=(nx,ny,nz))` with a triple buffer (`psi_old`, `psi`, `psi_new`).

### 3. Director field extraction (for topology)

Optional dedicated field `n(i,j,k) = unit vector` extracted either:

- From the major axis of the local displacement ellipse (the 6-phasor data M4 already implicitly carries), or
- From a separate director field evolved in tandem (LdG style)

This is the field whose winding integer = topological charge.

### 4. Topological initial-condition kernels

New `seed_*` kernels in the M2 style (renamed from M2's `charge_*` to better match IC formalism — these configure the field's topological state, not inject energy). Every scenario is a composition: **first `seed_vacuum` to fill the grid with the Lagrangian's ground state, then `seed_<defect>` to overwrite a localized region with the topological feature**:

- `seed_vacuum(n0)` → fill entire grid with ground state (e.g. `n(x) = ẑ` for a uniaxial vacuum, or `ψ = 0` for a scalar field). **Always the first step** — topology requires a reference state to wind around
- `seed_hedgehog(center, sign)` → `n(x) = ±(x − center)/|x − center|` near center (Experiment 2 ported to production)
- `seed_sine_gordon_kink(position, velocity)` → 1D / axisymmetric kink seeding (Experiment 1 ported to production)
- `seed_biaxial_hedgehog(center, axis)` → three-mass-scale lepton-family seeding (Experiment 6)
- `seed_spherical_harmonic(l, m)` → Close's spherical-harmonic seed for Experiment 7

See the **[Architectural Decision: Background Vacuum + Defects](#architectural-decision-background-vacuum--defects-m2-philosophy)** section below for why the vacuum step is mandatory and how it differs from M2's oscillating base wave.

**Naming note**: M2's kernels are called `charge_*` (e.g. `charge_gaussian`, `charge_falloff`) because they inject an energy pulse into the field — "charging" it up. M5's kernels don't inject energy; they configure the initial *topological configuration* (which way the director points at each voxel, which kink profile fills the domain). The right verb is **seed** (same as random-number generators, initial-condition generators in PDE literature). Keeping `charge_*` would be misleading — "seeding a hedgehog" is accurate; "charging a hedgehog" is not.

### 5. Hamiltonian energy density

Replace `E = ρV(fA)²` with the Hamiltonian density derived from the Lagrangian:

```text
H(i,j,k) = ½|ψ_dot|²  +  ½c²|∇ψ|²  +  V(ψ)
            kinetic     gradient     potential
```

where `ψ_dot = (ψ − ψ_old)/dt`. This energy is conserved by Noether's theorem (barring numerical drift), and its gradient gives force exactly the same way M4 does today (`F = −∇E`). The `rho_qgam * dx³ * (f·A)²` postulate becomes an *approximation* that holds in the linear limit.

### 6. Winding-number tracker

New diagnostic: compute topological charge `Q = (1/4π) ∮_S (∂_u n × ∂_v n) · n du dv` on spheres around each tracked WC location, as per Experiment 3. Returns integer values; becomes a new tracker alongside `amp_local_emarms_am`, `freq_local_cross_rHz`, `energy_local_aJ`.

---

## ARCHITECTURAL DECISION: BACKGROUND VACUUM + DEFECTS (M2-PHILOSOPHY)

A question that has to be settled before any kernel is written: **what is the state of M5's grid at t=0 when nothing is happening?**

OpenWave's existing methods answer this differently:

- **M2** — the field is a self-stabilized standing wave throughout the grid. "Emptiness" is still a wave field in equilibrium. WCs are disturbances *in* this field
- **M3** — the field is literally zero everywhere except where WCs emit. "Emptiness" is nothing; each WC generates its own wave into the void
- **M4** — same as M3 (analytical, zero background), just with a vector rather than a scalar output

**M5 must commit to M2's philosophy** — the field exists as a **background vacuum state** that is the ground state of the Lagrangian's potential V(ψ). A topological defect (hedgehog, kink) is a *deformation of this vacuum*, and the winding number is measured *relative to the vacuum orientation*. This is a hard mathematical requirement, not a stylistic preference:

> **A winding number only makes sense against a reference state.** If there is no vacuum field to wind around, there is no topological charge to measure. M3's "WCs emit into void" approach is structurally incompatible with topology.

### Static Vacuum vs. Oscillating Base Wave

Duda's vacuum is **not** the same as M2's base wave. The distinction:

| Aspect | M2 base wave (current) | M5 vacuum (planned) |
| --- | --- | --- |
| Field state when "empty" | `ψ_base = A₀·cos(kx)·cos(ωt)` — oscillating | `n(x) = ẑ` — static, aligned everywhere |
| What causes oscillation | Assumed intrinsic (all matter in universe, EWT) | Emerges from defect-vacuum interaction (time crystal, phi-4 mechanism) |
| Time dependence | Built in (ω·t) | None — oscillations arise from the dynamics, not the ground state |
| Mathematical role | Initial condition + driving term | Minimum of V(ψ) — the potential well |

So M5's vacuum is **more fundamental** than M2's base wave: it explains *why* fields oscillate (defects oscillate because they have mass and interact with V(ψ)) rather than assuming the universe comes pre-oscillating. See [3_LAGRANGIAN_FRAMEWORK.md § Impact on Base Wave Architecture](3_LAGRANGIAN_FRAMEWORK.md#impact-on-base-wave-architecture--m2-vs-m3-philosophy) for the full comparison and how this resolves Duda's time-crystal insight.

### Concrete Consequence for M5 `seed_*` Kernels

Every M5 initial-condition kernel follows a two-step pattern:

```text
Step 1: seed_vacuum(n0)            # fill the entire grid with the ground state
                                     (e.g. n(x) = ẑ everywhere, or ψ = 0)

Step 2: seed_<defect>(center, ...)  # overwrite a localized region with the defect
                                     (hedgehog, kink, spherical harmonic, etc.)
```

This mirrors how the sandbox experiments seed their fields (a hedgehog is "uniform director + radial override near the center"). The vacuum step is separate and shared — it's what Duda's framework requires to make topology measurable.

### Why This Does NOT Invalidate M3's Results

M3's lock-in, annihilation, and K=10 stability are real near-field wave physics — they happen regardless of whether the far field is modeled as vacuum or void. What M3 cannot do, because of its "no background" philosophy, is measure topological charge or produce far-field Coulomb without sinc barriers. M5 keeps M3's near-field physics as validated phenomena and adds the background-vacuum architecture where it matters (topology, charge quantization, far-field interactions). Cf. [3_LAGRANGIAN_FRAMEWORK.md § Practical Implication](3_LAGRANGIAN_FRAMEWORK.md).

---

## CONCRETE CODE MAPPING

Translating M4's `compute_voxel_wave` → M5's `step_voxel_pde`:

| M4 current | M5 equivalent | Delta |
| --- | --- | --- |
| Loop over WCs, evaluate closed-form per-WC (`wave_engine.py:161–211`) | No WC loop — voxel reads only its own `ψ_prev`, `ψ`, and 6 neighbors for Laplacian | Replace analytical eval with leapfrog step |
| `oscillator * direction` summation | `∇²ψ` via 6-point stencil (port from M2 `compute_laplacianL`, `wave_engine.py:549–560`) | Nonlocal via neighbors |
| N/A (linear wave) | `+ V'(ψ)` nonlinear term | New |
| `displacement_am[i,j,k] = A·oscillator·direction` (analytical result) | `psi_new[i,j,k] = 2·psi − psi_old + dt²·(c²·∇²ψ − V'(ψ))` | PDE step replaces direct evaluation |
| EMA on RMS amplitude (`wave_engine.py:244–255`) | **Same** — keep as tracker, but not as energy source | Unchanged |
| `E = ρV(fA)²` postulate (`wave_engine.py:287–293`) | `H = ½·ψ̇² + ½c²·(∇ψ)² + V(ψ)` Hamiltonian density | New computation, same scalar output |

Buffer layout change:

```text
M4 (current):
    displacement_am[nx,ny,nz]     ← single Vector(3) buffer, overwritten each frame

M5 (target):
    psi_prev[nx,ny,nz]            ← Vector(3) at t−dt
    psi[nx,ny,nz]                 ← Vector(3) at t (current)
    psi_new[nx,ny,nz]             ← Vector(3) at t+dt (computed)
    n_director[nx,ny,nz]          ← Vector(3) optional (for topology)
    # Swap at end of each step, same pattern as M2 wave_engine.py:743
```

---

## SANDBOX VERDICTS (all 8 experiments complete)

Phase 3 sandbox is complete (2026-04-16 / 2026-04-17). The architectural decisions for M5 are now resolved:

| Decision | Resolution | Source |
| --- | --- | --- |
| Does nonlinearity alone give K-selectivity? | ❌ **No** — Smolinski Ψ³ produces breathing oscillation only, no K-dependence | Exp 8 |
| Does topology alone give clean 1/r Coulomb? | ✅ **Yes** — R² = 0.993 post-relax, no sinc | Exp 2 |
| Does topology give integer charge quantization under perturbation? | ✅ **Yes** — Q = ±1 surface-independent, robust to 50% noise | Exp 3 |
| Does the biaxial-hedgehog mass hierarchy match lepton ratios? | ⚠️ Mechanism validated (E ∝ K linear); specific ratios achievable by choice of K, not derived. Full Q-tensor derivation = Exp 6.1 (deferred) | Exp 6 |
| Does Close's nonlinear vector equation produce localized solitons? | ⚠️ No static solitons from harmonic seeds — but Close's framework doesn't actually require them. Close's Eq. 19 IS a valid transverse vector wave equation for M5's base dynamics | Exp 7 v2 |
| Can our Combined W-L be derived from a Lagrangian, or must we replace it? | ⚠️ The M4 code's sum form IS a free-wave solution. **The documented product form `2A·sin(kr/2)·cos(kr/2−(ωt+φ))/r` is NOT** — it has a quadrature residual. M5 keeps the sum form, discards the product form | Exp 5 |
| Does the PDE evolution produce Klein-Gordon dispersion? | ✅ **Yes** — ω² = c²k² + m² validated to R² = 0.999982 across 9 modes | Exp 4 |
| Is leapfrog + nonlinear potential stable with Lorentz-correct kinematics? | ✅ **Yes** — kink v = 0.4997c (input 0.5c), width L/γ matched, energy drift 1.5e-6 | Exp 1 |

**Headline finding**: topology (Exps 1, 2, 3) is load-bearing; pure nonlinearity (Exps 7, 8) is insufficient on its own; Klein-Gordon wave dynamics (Exp 4) are the correct perturbative layer. Full context in [3b_lagrangian_experiments.md § OVERALL CONCLUSIONS](3b_lagrangian_experiments.md#overall-conclusions).

---

## IMPLEMENTATION PHASES (post-sandbox, ready to start)

The winning recipe is now known: **topology + Klein-Gordon dynamics + Close's vector wave equation + M3 near-field physics + Skyrme stabilizer + (eventually) LdG biaxial potential**. Phases below are scoped accordingly.

### Phase M5.0 — Scaffold

- [ ] Create `openwave/xperiments/m5_lagrangian_wave/` directory (mirror m4 structure)
- [ ] Copy M4's `WaveField`, `WaveCenter`, `WaveTrackers` data classes; extend with `psi_prev`, `psi_new` buffers for leapfrog
- [ ] Copy M4's flux-mesh visualization, granule rendering, 3-plane sampling (all unchanged)
- [ ] Port M2's 6-point Laplacian stencil (`compute_laplacianL` at `m2_laplace_propagation/wave_engine.py:527–562`) for the scalar-component case, then generalize to a 3-vector field for Close's `Q`
- [ ] Implement curl, divergence, and `curl(curl())` via `∇(∇·F) − ∇²F` (per Exp 7's implementation)
- [ ] **Physics invariant test**: with `V(ψ) = 0`, M5 must reproduce M2's free-wave behavior AND Exp 4's Klein-Gordon dispersion `ω² = c²k² + m²` for a quadratic potential. Fail this → there's a bug in the core loop

### Phase M5.1 — Port topology (from Exps 2, 3)

- [ ] Implement `seed_vacuum()` — fill grid with ground-state `n = ẑ`
- [ ] Implement `seed_hedgehog(center, sign)` — port Exp 2's weighted superposition + renormalization
- [ ] Implement Frank elastic energy `H = (K/2) · ∫ |∇n|² d³r` on the vector field
- [ ] Implement gradient-descent relaxation (tangent projection + unit-length renormalization + soft core pinning) — directly from Exp 2
- [ ] Validate: reproduce Exp 2's hedgehog-pair 1/d Coulomb on Taichi. Target R² > 0.99 across a separation sweep
- [ ] Implement `winding_number(center, radius)` tracker — port Exp 3's trilinear-sphere-sample + finite-difference surface integral

### Phase M5.2 — Wave dynamics (Close's Eq. 19 + Klein-Gordon)

- [ ] Implement time-stepping leapfrog for `∂²Q/∂t² = −c²·∇×∇×Q + (optional mass term −m²Q)` — Close's Eq. 19, generalized to include a mass term
- [ ] Validate: reproduce Exp 4's Klein-Gordon dispersion on the GPU. FFT-extract ω(k), fit ω² = c²k² + m²
- [ ] Validate: reproduce Exp 7's transverse wave dispersion (dipole/quadrupole seeds disperse as transverse elastic-solid waves)
- [ ] Add Close's nonlinear terms `−u·∇s + w×s` (from Eq. 21) as optional runtime flag

### Phase M5.3 — Hamiltonian energy + force

- [ ] Replace M4's postulated `E = ρV(fA)²` with the Hamiltonian density `H = ½|∂ₜψ|² + ½c²|∇ψ|² + V(ψ)` derived from the Lagrangian (matches Exp 5's Noether result)
- [ ] Verify `F = −∇E` still produces the expected particle motion — the mechanism survives, the energy source changes
- [ ] Cross-check: Exps 2 and 4 measured Hamiltonian energies; M5 must reproduce them

### Phase M5.4 — Multi-defect K=10 test (the headline goal)

- [ ] Seed K=10 hedgehog arrangement on M5 (1-3-6 tetrahedron geometry, per EWT)
- [ ] Measure stability under perturbation — does topology give perturbation-robust K=10 uniqueness?
- [ ] Compare to M3 Combined W-L baseline: does the sinc far-field barrier disappear when topology is active?
- [ ] Measure far-field force on a test particle: does clean 1/d² Coulomb emerge? (Integrate Exp 2's E(d) result with 2-defect M3 motion dynamics)
- [ ] **This is the full Phase 2 → Phase 3 → M5 validation loop.** If K=10 is uniquely stable here, we've done it

### Phase M5.5 — Skyrme stabilizer (if M5.4 reveals defect collapse)

- [ ] Per-Derrick's theorem, a bare topological defect in 3D is unstable to scale change. Add a Skyrme higher-derivative term `(∂_μ s × ∂^μ s)²` to stabilize
- [ ] Scan Skyrme coefficient to find physically meaningful range
- [ ] Rerun M5.4 with Skyrme term; compare stability

### Phase M5.6 — Biaxial LdG (deferred long-term)

- [ ] Full 3×3 Q-tensor dynamics with LdG potential `V(Q) = a·Tr(Q²) − b·Tr(Q³) + c·(Tr Q²)²`
- [ ] Goal: three distinct lepton energy scales emerge from a single set of (a, b, c) LdG parameters
- [ ] This is the "electron, muon, tau from biaxial geometry" experiment — a significant undertaking beyond initial M5

### What M5 does NOT implement

From the sandbox findings, these are ruled out or de-prioritized:

- **The documented Combined W-L product form `2A·sin(kr/2)·cos(kr/2−(ωt+φ))/r`** — Exp 5 showed this isn't a free-wave solution. M5 uses the M4 *sum form* `A·[sin(kr+ωt+φ) + sin(kr−ωt−φ)]/(kr)` if analytical standing waves are needed
- **Smolinski's Ψ³ as primary K-selectivity mechanism** — Exp 8 falsified. Ψ³ may still be a stabilizer *inside* a topologically protected configuration, but M5 does not rely on it for geometric selectivity
- **Seeded-soliton emergence from spherical harmonics as a physics mechanism** — Exp 7 v2 confirmed Close's framework doesn't predict this. Particles in M5 = topological defects (from seed), not emergent solitons from wave seeds

---

## STATUS

- [x] Architecture analysis complete (this document)
- [x] **All 8 sandbox experiments complete** — see [3b_lagrangian_experiments.md](3b_lagrangian_experiments.md):
  - ✅ Exp 1 (Sine-Gordon kinks), ✅ Exp 2 (Hedgehog Coulomb), ✅ Exp 3 (Winding quantization), ✅ Exp 4 (Klein-Gordon dispersion)
  - ⚠️ Exp 5 (Lagrangian derivation — W-L product form falsified), ⚠️ Exp 6 (lepton mechanism; specific ratios deferred), ⚠️ Exp 7 v2 (Close's actual equations implemented)
  - ❌ Exp 8 (Smolinski Ψ³ K-selectivity falsified)
- [x] **Winning recipe identified**: topology + Klein-Gordon + Close's Eq. 19 + M3 near-field + Skyrme stabilizer
- [ ] M5.0 — Scaffold (Taichi structure, triple buffer, Laplacian, curl/div operators)
- [ ] M5.1 — Port topology from Exps 2, 3 (`seed_vacuum`, `seed_hedgehog`, Frank energy, winding tracker)
- [ ] M5.2 — Wave dynamics from Close's Eq. 19 + Klein-Gordon mass term, validate Exp 4 dispersion
- [ ] M5.3 — Hamiltonian energy (replaces postulated `E = ρV(fA)²`)
- [ ] M5.4 — **Headline test**: K=10 hedgehog stability under perturbation + far-field Coulomb recovery
- [ ] M5.5 — Skyrme stabilizer (conditional on M5.4)
- [ ] M5.6 — Biaxial LdG Q-tensor (long-term; for lepton mass derivation)

**Next action**: begin **M5.0 scaffold**. The sandbox is done; the recipe is known. Move to production implementation on the Taichi engine.
