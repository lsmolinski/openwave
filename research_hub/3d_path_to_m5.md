# PATH TO M5 — LAGRANGIAN-FIELD METHOD

Implementation plan for **M5 / LAGRANGIAN-FIELD METHOD** (directory `openwave/xperiments/m5_lagrangian_field/`), the production field engine that graduates sandbox-validated Lagrangian / topological physics onto the GPU-accelerated OpenWave platform. This document references concrete code in the current engines and defines what M5 inherits, replaces, and adds.

**Naming**: **M5 / LAGRANGIAN-FIELD METHOD** (renamed 2026-04-19 from the earlier "Lagrangian-Wave Method" working title). The rename reflects what the method actually does: it is a full **Lagrangian field-theory** simulator, not just a wave engine. "Lagrangian" refers to the variational formalism (L = T − V, Euler-Lagrange equations, action principle) from which the equation of motion is derived; "Field" because the engine integrates a unified PDE `∂²_tψ = c²∇²ψ − ∂V/∂ψ` that simultaneously handles wave propagation *and* preserves topology via the potential `V(ψ)` — two channels, one equation. The module that implements this is `lagrangian_engine.py` (not `wave_engine.py`) for the same reason. See [3b § What wave equation does M5 solve?](3b_concept_review.md#what-wave-equation-does-m5-solve-is-force-still-e) for the full rationale — particularly the 7-item breakdown of what the engine does (only one item is strictly wave propagation). The name distinguishes M5 from M1–M4's "Wave Method" naming (those really are wave engines; M5 is a field-theory engine).

**Spec inputs**:

- [3_LAGRANGIAN_FRAMEWORK.md](3_LAGRANGIAN_FRAMEWORK.md) — Lagrangian framework evaluation, 8 experiments, Duda/Close context
- [3a_lagrangian_experiments.md](3a_lagrangian_experiments.md) — sandbox numerical results (experiment-by-experiment)
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
| Triple-buffer pattern (`psi_prev`, `psi`, `psi_new`) | Field allocation in `WaveField` + update at `wave_engine.py:649–660` + swap at line 743 | Same buffering needed for leapfrog; extended to vector component by component |
| Leapfrog update formula `ψ_new = 2ψ − ψ_prev + (c·dt)²·∇²ψ` | `wave_engine.py:649–653` | Replace RHS with `(c·dt)²·∇²ψ − dt²·dV/dψ` (add nonlinear potential term) |
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
ψ_new = 2·ψ − ψ_prev + dt²·[c²·∇²ψ − ∂V/∂ψ]
```

The specific `V(ψ)` is selected from sandbox results:

| Candidate | Potential term | Status |
| --- | --- | --- |
| Sine-Gordon family | `V(φ) = (m²c⁴/ℏ²)(1 − cos φ)` → `∂V/∂φ = (m²c⁴/ℏ²) sin(φ)` | ✅ validated (Exp 1) — useful for 1D analogs, not 3D hedgehogs |
| Smolinski Ψ³ | `V(ψ) = (k/4)·ψ⁴` → `∂V/∂ψ = k·ψ³` | ❌ falsified for K-selectivity (Exp 8); usable only as in-configuration stabilizer |
| Landau-de Gennes | `V(M) = a·Tr(M²) − b·Tr(M³) + c·(Tr M²)²` | ⚠️ mechanism validated (Exp 6); specific parameters deferred to M5.6 |
| Close's elastic solid | Eq. 19 `∂²Q = −c²·∇×∇×Q` (linear limit) + **Eq. 23** (particle equation, preserves `∇·s = 0`, Close's recommendation 2026-04-18) — optional Eq. 21 nonlinear terms `−u·∇s + w×s` for comparison | ✅ Eq. 19 validated (Exp 7 v2); Eq. 23 selected per Close's explicit guidance — **M5's base wave dynamics layer** |
| Klein-Gordon mass term | `V(ψ) = ½m²·ψ²` → `∂V/∂ψ = m²·ψ` | ✅ validated (Exp 4) — **selected as M5's perturbation-mass mechanism** |
| Mixed / composite | linear combination (topology seeding + Close dynamics + KG mass + optional Skyrme) | **adopted** — full recipe in [3a § Winning Approach](3a_lagrangian_experiments.md#winning-approach-for-m5) |

### 2. True independent vector components (no radial constraint)

M4 collapses its vector to radial (`oscillator * direction`, line 211). M5 evolves each of `(ψ_x, ψ_y, ψ_z)` independently via its own PDE step. No radial projection; topology emerges from the full vector field.

**Storage**: either three coupled scalar fields, or one `ti.Vector.field(3, dtype=ti.f32, shape=(nx,ny,nz))` with a triple buffer (`psi_prev`, `psi`, `psi_new`).

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

where `ψ_dot = (ψ − ψ_prev)/dt`. This energy is conserved by Noether's theorem (barring numerical drift), and its gradient gives force exactly the same way M4 does today (`F = −∇E`). The `rho_qgam * dx³ * (f·A)²` postulate becomes an *approximation* that holds in the linear limit.

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
| `displacement_am[i,j,k] = A·oscillator·direction` (analytical result) | `psi_new[i,j,k] = 2·psi − psi_prev + dt²·(c²·∇²ψ − V'(ψ))` | PDE step replaces direct evaluation |
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

**Headline finding**: topology (Exps 1, 2, 3) is load-bearing; pure nonlinearity (Exps 7, 8) is insufficient on its own; Klein-Gordon wave dynamics (Exp 4) are the correct perturbative layer. Full context in [3a_lagrangian_experiments.md § OVERALL CONCLUSIONS](3a_lagrangian_experiments.md#overall-conclusions).

---

## IMPLEMENTATION PHASES (post-sandbox, ready to start)

The winning recipe is now known: **topology + Klein-Gordon dynamics + Close's vector wave equation + M3 near-field physics + Skyrme stabilizer + (eventually) LdG biaxial potential**. Phases below are scoped accordingly.

---

## RESOLUTION & PERFORMANCE PLAN

A separate, mandatory engineering plan that determines which phases run on the local M4 Max and which need rented GPU time. M2 / M3 hit a hard ceiling at ~350M voxels because they had to resolve the EWT base wave at `λ₀ = 2.85 × 10⁻¹⁷ m`. M5 changes this fundamentally.

### Why M5's resolution budget is much friendlier than M2/M3's

| Constraint | M2 / M3 | M5 |
| --- | --- | --- |
| Vacuum is | An oscillating base wave at `f₀ = 1.05 × 10²⁵ Hz` | A **static** ground state of `V(ψ)` (no oscillation in vacuum) |
| Grid must resolve | `λ₀ = 2.85 × 10⁻¹⁷ m` (the EWT base wave) | `λ_C = ℏ/(mc)` of the **defect** (its Compton wavelength) |
| Min `dx` (Nyquist ~10–12× / wavelength) | `~3 × 10⁻¹⁸ m` | `~3 × 10⁻¹⁴ m` for an electron — **10⁴× larger** |
| Min `dt` (CFL `dt ≤ dx / (c·√3)`) | `~6 × 10⁻²⁷ s` | `~6 × 10⁻²³ s` for an electron — **10⁴× larger** |
| Reachable particle on M4 Max | Below electron (couldn't reach electron radius) | Electron well in scope; muon tight; tau requires AMR |

In one sentence: **M5 doesn't simulate the EWT vacuum's oscillations — there are none; the vacuum is a static ground state. Only the defect's intrinsic length scale needs to be resolved.**

### Per-particle defect scales (the "what does M5 actually need to resolve" table)

| Particle | Mass | Compton λ_C = ℏ/(mc) | Zitterbewegung ω = 2mc²/ℏ | T_Z = 2π/ω |
| --- | --- | --- | --- | --- |
| **Electron** | 0.511 MeV | **3.86 × 10⁻¹³ m** (~386 fm) | 1.55 × 10²¹ rad/s | 4.05 × 10⁻²¹ s |
| **Muon** | 105.7 MeV | 1.87 × 10⁻¹⁵ m | 3.21 × 10²³ rad/s | 1.96 × 10⁻²³ s |
| **Tau** | 1777 MeV | 1.11 × 10⁻¹⁶ m | 5.39 × 10²⁴ rad/s | 1.17 × 10⁻²⁴ s |
| **Up quark** (constituent ~2.2 MeV) | 2.2 MeV | 8.97 × 10⁻¹⁴ m | 6.7 × 10²¹ rad/s | 9.4 × 10⁻²² s |
| **Proton** (composite) | 938 MeV | 2.10 × 10⁻¹⁶ m | 2.85 × 10²⁴ rad/s | 2.20 × 10⁻²⁴ s |
| **Neutrino** (~0.1 eV) | 0.1 eV | ~2 × 10⁻⁶ m | ~3 × 10¹⁴ rad/s | ~2 × 10⁻¹⁴ s |

### Per-phase grid budgets (target on M4 Max 48 GB / 18.4 TFLOPS)

| Phase | Particle | dx (≈λ_C/12) | Domain (≈8·λ_C) | N (per side) | Voxels | Memory* | Steps for 10 T_Z | Wall-clock estimate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| M5.0 / M5.1 / M5.2 (KG dispersion, Coulomb) | n/a (geometric tests) | natural-units `1` | `64` | 256 | 17M | ~1 GB | n/a (relaxation, not Zitterbewegung) | minutes |
| **M5.4** (electron stability, e⁺/e⁻ pair) | electron | 3.2e-14 m | 3.1e-12 m | 256–384 | 17M–57M | 1–4 GB | ~500k steps | **~12–24h** (uniform grid, no AMR) |
| M5.7 (Cornell, light quarks) | u/d quark | ~7.5e-15 m (≈λ_C(u)/12) | ~7e-13 m | 96–128 | 1M–2M | <1 GB | ~50k–500k | hours |
| M5.6 (muon, tau) | muon, tau | ~1.5e-16 m / ~9e-18 m | tight | uniform grid impractical | — | — | — | **needs AMR** |
| M5.8 (Zitterbewegung table for all 7 particles) | all | per-particle `dx` | per-particle | varies | varies | varies | 100×T_Z | **needs AMR + cloud burst** |

*Memory assumes Vector(3) for Q + 3 buffers (prev/curr/new) + winding/director field + Hamiltonian tracker, ~64 bytes/voxel.

### Performance optimizations to bake into M5.0 from Day 1

#### Tier 1 — architectural decisions (fix before any kernel is written)

1. **Native `ti.Vector.field(3, ...)`** with single triple buffer (`psi_prev`, `psi`, `psi_new`) — halves memory traffic vs three independent scalar fields and lets Taichi vectorize 3-component operations
2. **Unit system** — two complementary scalings to keep numerics in float32's sweet spot:
    - **Inherit OpenWave's existing scaled SI** from `openwave/common/constants.py` (lines 10–29): spatial in **attometers** (`1e-18 m`, suffix `_am`) and temporal in **rontoseconds** (`1e-27 s`, suffixes `_rs`, `_amrs`, `_rHz`). Already proven to keep 6–7 significant digits at f32 across M2/M3/M4 — base wavelength `~28.5 am`, wave speed `~0.3 am/rs`, etc. M5 must **continue using these units** for all stored fields and all I/O — don't fork a parallel unit system, don't reinvent precision-protection
    - **Add a separate, internal natural-unit kernel scaling** (`c = 1`, `λ_C(defect-of-interest) = 1`, `ℏ = 1`) for *intra-kernel* arithmetic only. The defect's Compton wavelength is the Lagrangian's natural reference length (electron `λ_C ≈ 386,000 am` is no longer near-1.0 in attometers — leaving f32 sweet spot if used directly in CFL / dispersion / amplitude-sweep math). Conversion happens at kernel entry/exit; no field is ever stored in natural units
    - **Net effect**: data layer stays in proven `_am` / `_rs` units (compatible with all OpenWave tooling, diagnostics, visualization); kernel arithmetic runs in natural units (CFL trivially `dt ≤ dx/√3`, amplitude sweeps `A ∈ {0.25, 0.5, 1, 2}` are O(1), Klein-Gordon mass term is O(1)). Both scalings are needed — neither alone covers M5's precision range across the lepton mass hierarchy `0 < δ ≪ 1 ≪ g`
3. **AMR-ready storage layer** — even if M5.0 ships with a uniform grid, design the field-storage abstraction so a future swap to octree-based AMR (M5.6 / M5.8) doesn't require a rewrite. Don't bake in fixed `(nx, ny, nz)` everywhere
4. **Domain size by defect physics**, not by EWT vacuum — far-field Coulomb cleanly captured at ~8·λ_C, so the box stays small

#### Tier 2 — kernel-level performance hacks

5. **Tile-and-reuse for the 6-point stencil** via Taichi `BlockLocal` / shared-memory loads. The Laplacian is bandwidth-bound, not flop-bound; tiling cuts global-memory traffic ~6×. Free 3–5× speedup
6. **Symplectic / Verlet integrator** instead of vanilla leapfrog. Same cost per step, but ~2–4× larger usable `dt` while preserving energy conservation. **Critical for M5.8** (Zitterbewegung needs energy-clean integration over hundreds of periods)
7. **Operator splitting**: `c²∇²ψ` (cheap, vectorized) and `−∂V/∂ψ` (potentially expensive nonlinear) in separate sub-steps. Each operator individually has better stability → larger overall `dt`
8. **Float32 by default; float64 only for diagnostic snapshots**. 2× speedup + 2× memory headroom; the leapfrog is well-conditioned at f32 in our scale
9. **Skip diagnostic trackers (Hamiltonian + winding) on most frames** — compute every Nth frame (N=10–100). They're observables, not load-bearing on the dynamics
10. **Dirty-tile mask** (the M5 reincarnation of M4's "selective voxel" optimization). Voxels far from any defect are at exactly the vacuum value `ψ_vacuum`; the leapfrog update there is a no-op. Use a chunk-level mask (e.g. 16³ blocks) — skip blocks whose `‖ψ - ψ_vacuum‖` is below epsilon. Realistic 2–10× speedup for sparse single-defect scenarios. Note: unlike M4, this is an *optimization*, not a different physics — a "skipped" block must produce identical output to a computed one (modulo floating-point epsilon)

#### Tier 3 — long-run techniques (deferred to M5.5+ as needed)

11. **Adaptive Mesh Refinement (AMR)** — fine grid near the defect core, coarse grid in the far field. **The single biggest win for M5.6 (muon/tau) and M5.8 (full Zitterbewegung table)**. Cell size scales with local `λ_C`. Realistic 10²–10³× voxel reduction vs uniform fine-grid for a localized defect
12. **Comoving frame for moving defects** — track the defect's center, evolve in its rest frame. Eliminates needing the grid to span the defect's full traveled distance
13. **Multi-scale time stepping** — subcycle small-`dx` (defect core) at fine `dt`, subcycle coarse far-field at large `dt`. Pairs naturally with AMR
14. **Spectral / pseudo-spectral solver for the linear part** (Klein-Gordon) — exact dispersion at any `dx`, FFT-based, ~10× larger usable `dx`. Worth evaluating as an alternative validator for M5.2
15. **Absorbing boundaries / PML** (M2 has a partial start at `wave_engine.py:662–699`). Lets us shrink the domain without reflections — typically ~2× linear / ~8× volume reduction

### Compute strategy — when to leave the M4 Max

The M4 Max is excellent for M5.0–M5.4 (uniform grid, electron-scale runs). The 12–24 hour electron headline run is borderline (long jobs at high GPU utilization risk thermal events; ties up the dev machine for a day). **A laptop is not the right tool for multi-day production runs.**

Per the existing analysis in [`dev_docs/distributed_computing/m4_max_vs_aws.md`](../dev_docs/distributed_computing/m4_max_vs_aws.md) and [`system_upgrade_options.md`](../dev_docs/distributed_computing/system_upgrade_options.md), key constraints:

- **Taichi is single-GPU bound**. Multi-GPU AWS instances (`p3.8xlarge`, `p4d.24xlarge`) only use 1 of N GPUs without an MPI / domain-decomposition rewrite (~2–3 months engineering)
- **Single-GPU AWS instances are barely faster than M4 Max** — V100 is ~1.16× M4 Max FP32 at $3.06/hr; A10G is ~2.3× at $1.01/hr. None are step-change improvements
- **The right cloud play is parameter sweeps**, not single-run speedup. Run 100 amplitude/seed/coupling configurations in parallel on 100 cheap spot instances, finish in the wall-clock time of one
- **The hardware step-change targets are M5 Ultra (2026, ~50–60 TFLOPS, 512 GB) or AMD MI300X (~163 TFLOPS, 192 GB, $15k+)** — both single-GPU, no code rewrite

#### Decision tree per M5 phase

| Phase | Workload | Recommended platform | Rationale |
| --- | --- | --- | --- |
| M5.0 / M5.1 / M5.2 / M5.3 | Small-grid validation, dispersion fits, Coulomb sweep | **M4 Max local** | Minutes-to-hours; tight dev loop is the limiting factor, not compute |
| M5.4 (electron) — single run | 17M–57M voxels, ~500k steps, ~12–24h | **M4 Max overnight, with caveats** | Run when machine isn't needed; monitor for thermal throttling. Pre-validate on small grid first to catch bugs before long run |
| M5.4 — parameter sweep (Y_1^0 amplitude, axis tuning, lifetime characterization) | 4–20 variations × ~6h each | **AWS spot fleet (`g5.xlarge` × N)** | Embarrassingly parallel; total wall-clock = single-run time. Spot pricing (~$0.30/hr discounted) makes a 20-config sweep ~$36 |
| M5.5 (Skyrme) | Conditional; same scale as M5.4 | **Same as M5.4** | — |
| M5.6 (muon, tau, biaxial LdG, Faber regularization) | AMR essential; without it `dx ~ 10⁻¹⁸ m`, voxels astronomical | **AWS A100 (`p4d.24xlarge` 1 GPU)** OR **wait for M5 Ultra (2026)** | Higher-mass particles need AMR + more memory. AMR development on M4 Max, production runs on cloud-or-future-Ultra |
| M5.7 (Cornell potential) | 1D-line topology, lighter than M5.4 | **M4 Max local** | u/d quark scale (~few MeV) is friendlier than muon |
| M5.8 (Zitterbewegung table for 7 particles) | AMR essential; long runs (≥10 T_Z each) for FFT freq extraction | **Cloud or M5 Ultra** | Each particle is a separate run; parallelize across cloud spot instances |

#### Cloud burst recipe (for when local isn't enough)

When a phase requires cloud compute, the playbook is:

1. **Develop and validate on M4 Max** at small scale (256³ or smaller). Confirm correctness end-to-end
2. **Containerize** the run (Taichi + CUDA backend on a Docker image; switch from `ti.metal` to `ti.cuda` is a single-line change)
3. **Use AWS spot instances** for the heavy run. `g5.xlarge` (1× A10G, 24 GB, 31 TFLOPS) at spot pricing ~$0.30/hr is the price/performance sweet spot for OpenWave's single-GPU workloads
4. **Persist results to S3**, pull diagnostics back to M4 Max for analysis. Don't try to render or interactively debug remotely
5. **For parameter sweeps**: AWS Batch or Step Functions launches N spot instances simultaneously, each running one config. Total wall-clock = single-run time

**Estimated cloud spend for the full M5 roadmap**: at spot pricing, the entire M5.4 + M5.7 + M5.8 production run sequence (including parameter sweeps and the 7-particle Zitterbewegung table) is **a few hundred to ~$2,000** total. Not a budget killer, but worth tracking as a separate cost line.

#### Hardware upgrade decision points

| Trigger | Action |
| --- | --- |
| M5.4 production runs become routine (multiple per week) | Evaluate Mac Studio M3 Ultra 80-core / 256–512 GB or wait for M5 Ultra |
| M5.6 / M5.8 require routine production runs and AMR is implemented | Same — single-GPU step change is the right play, not multi-GPU |
| Project moves to fluent multi-particle ensemble simulation (M6.x, M7.x) | Re-evaluate the multi-GPU rewrite; the 2–3 month MPI investment becomes economically viable when single-GPU runs become a regular bottleneck |
| Need >48 GB memory specifically (not compute) | M4 Max 128 GB upgrade (~$4,000) is the cheapest unlock to ~1.5B voxels |

> **Bottom line**: M5.0–M5.4 are M4-Max-feasible with the Tier 1+2 optimizations in place. M5.6 and M5.8 are the phases where AMR is mandatory and cloud burst (or an M5 Ultra) becomes the right tool. Plan for that pivot now so the AMR-ready storage layer is in M5.0's scaffold rather than a retrofit later.

---

### Phase M5.0 — Scaffold

- [ ] Create `openwave/xperiments/m5_lagrangian_field/` directory (mirror m4 structure)
- [ ] **Rename the engine module**: `wave_engine.py` → `lagrangian_engine.py` for M5. Rationale: M5's core loop integrates a Lagrangian-derived PDE (`∂²_tψ = c²∇²ψ − ∂V/∂ψ`) that simultaneously handles wave propagation *and* preserves topology via the potential `V(ψ)`. "Wave" is only one of the two channels the engine produces, so `lagrangian_engine.py` reflects what the module actually is to a new reader. M1–M4 keep `wave_engine.py` (they really are wave engines, no topology layer). See [3b § What wave equation does M5 solve?](3b_concept_review.md#what-wave-equation-does-m5-solve-is-force-still-e) for the reasoning
- [ ] Copy M4's `WaveField`, `WaveCenter`, `WaveTrackers` data classes; extend with `psi_prev`, `psi_new` buffers for leapfrog
- [ ] Copy M4's flux-mesh visualization, granule rendering, 3-plane sampling (all unchanged)
- [ ] Port M2's 6-point Laplacian stencil (`compute_laplacianL` at `m2_laplace_propagation/wave_engine.py:527–562`) for the scalar-component case, then generalize to a 3-vector field for Close's `Q`
- [ ] Implement curl, divergence, and `curl(curl())` via `∇(∇·F) − ∇²F` (per Exp 7's implementation)
- [ ] **Inherit `openwave/common/constants.py` scaled SI** (`_am` / `_rs` / `_rHz`) for ALL stored fields and I/O. This is OpenWave's production-validated f32-precision infrastructure (lines 10–29 of `constants.py`); M5 must not fork its own units
- [ ] **Add a kernel-internal natural-unit scaling** (`c = 1`, `λ_C(electron) = 1`, `ℏ = 1`) for intra-kernel arithmetic only — convert at kernel entry/exit. Required because the electron's `λ_C ≈ 3.86e5 am` is no longer near-1.0 in attometers, which hurts CFL/dispersion math at f32. See [Resolution & Performance Plan § Tier 1](#tier-1--architectural-decisions-fix-before-any-kernel-is-written)
- [ ] **Use native `ti.Vector.field(3, ...)`** with single triple buffer (`psi_prev`, `psi`, `psi_new`) — not three independent scalar fields. See Resolution & Performance Plan § Tier 1
- [ ] **Design field-storage layer to be AMR-ready** (no fixed `(nx, ny, nz)` constants embedded in kernels). M5.0 ships uniform-grid; M5.6 / M5.8 retrofit AMR without a rewrite
- [ ] **Implement chunk-level dirty-tile mask** (e.g. 16³ blocks) to skip vacuum regions. Mandatory for sparse single-defect runs (M5.4 onward)
- [ ] **Symplectic / Verlet integrator** instead of vanilla leapfrog — required for M5.8 long-run energy conservation. See Resolution & Performance Plan § Tier 2
- [ ] **Tile-and-reuse the 6-point Laplacian** via Taichi `BlockLocal` — bandwidth-bound stencil benefits 3–5× from tiling
- [ ] **Physics invariant test**: with `V(ψ) = 0`, M5 must reproduce M2's free-wave behavior AND Exp 4's Klein-Gordon dispersion `ω² = c²k² + m²` for a quadratic potential. Fail this → there's a bug in the core loop

### Phase M5.1 — Port topology (from Exps 2, 3)

- [ ] Implement `seed_vacuum()` — fill grid with ground-state `n = ẑ`
- [ ] Implement `seed_hedgehog(center, sign)` — port Exp 2's weighted superposition + renormalization
- [ ] Implement Frank elastic energy `H = (K/2) · ∫ |∇n|² d³r` on the vector field
- [ ] Implement gradient-descent relaxation (tangent projection + unit-length renormalization + soft core pinning) — directly from Exp 2
- [ ] Validate: reproduce Exp 2's hedgehog-pair 1/d Coulomb on Taichi. Target R² > 0.99 across a separation sweep
- [ ] Implement `winding_number(center, radius)` tracker — port Exp 3's trilinear-sphere-sample + finite-difference surface integral

### Phase M5.2 — Wave dynamics (Close's Eq. 23 + Eq. 19 linear limit + Klein-Gordon)

> **Refinement from Robert Close (2026-04-18)**: use **Eq. 23 as the particle equation**, not Eq. 21. Eq. 23 preserves `∇·s = 0` (zero divergence of spin density), which is a physical invariant Exp 7's Eq. 21 implementation did not enforce. Eq. 19 remains the linear free-wave limit.

- [ ] Implement time-stepping leapfrog for Close's **Eq. 23** as the particle equation, enforcing `∇·s = 0` at each step (divergence-cleaning projection, or vector-potential `s = ∇×A` formulation that makes zero-div automatic)
- [ ] Keep Eq. 19 `∂²Q/∂t² = −c²·∇×∇×Q + (optional mass term −m²Q)` as the V=0 linear limit
- [ ] Validate: reproduce Exp 4's Klein-Gordon dispersion on the GPU. FFT-extract ω(k), fit ω² = c²k² + m²
- [ ] Validate: reproduce Exp 7's transverse wave dispersion (dipole/quadrupole seeds disperse as transverse elastic-solid waves)
- [ ] Add Close's nonlinear terms `−u·∇s + w×s` (from Eq. 21) as optional runtime flag for comparison
- [ ] **Resonance-hunt protocol** (per Close's recommendation): seed `Y_1^0` (electric-dipole harmonic) at amplitudes `A ∈ {0.25λ, 0.5λ, 1.0λ, 2.0λ}` (displacement-comparable-to-wavelength regime); measure localization lifetime. Particles in Close's framework are *unstable resonances* at specific amplitude/wavelength ratios, not static solitons. Success = extended-lifetime localization, not perfect stability

### Phase M5.3 — Hamiltonian energy + force

- [ ] Replace M4's postulated `E = ρV(fA)²` with the Hamiltonian density `H = ½|∂ₜψ|² + ½c²|∇ψ|² + V(ψ)` derived from the Lagrangian (matches Exp 5's Noether result)
- [ ] Verify `F = −∇E` still produces the expected particle motion — the mechanism survives, the energy source changes
- [ ] Cross-check: Exps 2 and 4 measured Hamiltonian energies; M5 must reproduce them

### Phase M5.4 — Electron stability + dynamic Coulomb recovery (the headline goal)

> **Paradigm shift from EWT (2026-04-19)**: the electron is NOT a K=10 tetrahedron of constituents in M5. Per Dr. Duda's framework, the electron is a **single biaxial hedgehog** — one topological defect with a specific axis choice (the δ ~ ℏ axis of the LdG order parameter). Muon and tau are single hedgehogs on the other two axes (unity and g scales respectively, tested in M5.6). EWT's "K" as a combinatorial parameter does not apply at the lepton scale. K re-emerges only as a *descriptive count* for composite-particle configurations (quarks in nucleons, nucleons in nuclei, electrons in atomic orbitals) — best-practice labeling for those deferred until M5.7+.
>
> **Refinement from Robert Close (2026-04-18)**: "Unless you have a good way to model an infinite system, I doubt that you will find completely stable non-radiating solutions." Success criterion is therefore **long-lived resonance with measurable lifetime**, not perfect stability.

**Single-particle stability test (electron as a defect)**:

- [ ] Seed a single biaxial hedgehog with the electron-axis choice (δ ~ ℏ)
- [ ] Evolve the field under the full Lagrangian (Close's Eq. 23 + KG mass term + LdG biaxial potential + Frank elastic + Skyrme if needed)
- [ ] Measure defect stability under perturbation: is the hedgehog a long-lived resonance? What's its measurable lifetime?
- [ ] Compare to variants (different axis choices, topological-vortex-loop alternatives for neutrino): does the δ-axis hedgehog have a measurably longer lifetime than alternatives, as a proxy for "electron is the lightest stable charged lepton"?

**Pair-interaction test (dynamic Coulomb + annihilation)**:

- [ ] Seed a hedgehog + anti-hedgehog pair (electron + positron) at various separations
- [ ] Compare to M3 Combined W-L baseline: does the sinc far-field barrier disappear when topology is active?
- [ ] Measure far-field force *dynamically* (not just statically as in Exp 2): does clean 1/d² Coulomb emerge as the pair moves toward each other?
- [ ] Test annihilation: when the pair closes to near-zero separation, do the winding numbers cancel and the stored energy radiate outward as Klein-Gordon waves at ~511 keV each?

**Mechanism-distinguishing tests (3b Q&A → M5.4 testing points)**:

- [ ] **Topology-vs-wave annihilation test**: seed a *same-winding* pair (Q = +1 and Q = +1) and force them together. They must NOT annihilate (topology forbids `Q_total = +2 → 0` smoothly; only `Q + (−Q) = 0` is topologically allowed). This directly distinguishes M5's "annihilation = topological cancellation" from M3's "annihilation = wave cancellation" — same-winding pairs would annihilate in M3 if waves were forced to overlap, but cannot in M5
- [ ] **Phase-as-derived test**: seed two electrons (same winding Q = −1) with *different* initial oscillation phases. Confirm both behave as electrons (same charge sign) regardless of phase choice. This validates that charge sign comes from winding only, not from emitted-wave phase (which is now a derived output, not a free parameter)
- [ ] **Same-mass Zitterbewegung lock-in test (positronium)**: at sub-λ_C separations, e⁺/e⁻ should form standing-wave bound states from direct Zitterbewegung-frequency interference (matching frequencies → coherent standing waves at λ_Z/2 wells). Lifetime should match positronium's ~125 ps (para) and ~142 ns (ortho) — this validates that *same-mass* pairs DO use direct Zitterbewegung lock-in (in contrast to mass-mismatched pairs like hydrogen, see "Beyond M5.8" cross-mass-class note below)

**Success metrics**:

- [ ] Single biaxial hedgehog (electron-axis) is a long-lived resonance under moderate perturbation (lifetime substantially longer than decoy variants)
- [ ] Hedgehog + anti-hedgehog pair dynamics reproduce 1/d² Coulomb attraction + annihilation; same-winding pairs do NOT annihilate; e⁺/e⁻ standing-wave lock-in reproduces positronium lifetimes
- [ ] **This is the full Phase 2 → Phase 3 → M5 validation loop.** If the single biaxial hedgehog is a long-lived electron-like resonance AND the e⁺/e⁻ pair reproduces Coulomb + annihilation + positronium lock-in dynamically AND same-winding pairs cannot annihilate, the loop is closed

### Phase M5.5 — Skyrme stabilizer (if M5.4 reveals defect collapse)

- [ ] Per-Derrick's theorem, a bare topological defect in 3D is unstable to scale change. Add a Skyrme higher-derivative term `(∂_μ s × ∂^μ s)²` to stabilize
- [ ] Scan Skyrme coefficient to find physically meaningful range
- [ ] Rerun M5.4 with Skyrme term; compare stability

### Phase M5.6 — Biaxial LdG (deferred long-term)

> **Refinement from Jarek Duda (2026-04-17)**: the lepton mass hierarchy is **not ad-hoc tuning** — it comes from the natural scale separation `0 < δ << 1 << g` where `δ ~ ℏ` (QM scale, twist eigenvalue), `1` (unity / matter scale, tilt eigenvalue), `g` (gravity scale, boost eigenvalue). The Exp 6 requirement for ~3477:1 axis ratio is *physically motivated* by this three-scale hierarchy. Additionally, "the main mass contributions of leptons come from LdG-like potential in regularization, it is the most difficult to include in simulation."
>
> **Follow-up refinement from Jarek Duda (2026-04-19)**: `(δ, g)` are **calibration parameters**, not ab-initio derivations: *"while these delta, g parameters describe Lagrangian contributions of QM and gravity, their exact choice seem to require numerical simulations."* No analytical form to pull from — M5.6 iterates them against observed physics. On regularization: the details are still an open research problem, but **Manfried Faber's scheme (slightly different potential, produces running-coupling effect) is the recommended starting point** — Faber et al., *Universe* 11 (2025) 113 (<https://www.mdpi.com/2218-1997/11/4/113>) and arxiv:2604.12021.

- [ ] Full 3×3 Q-tensor dynamics with LdG potential `V(Q) = a·Tr(Q²) − b·Tr(Q³) + c·(Tr Q²)²` using traces of powers from original Landau-de Gennes theory
- [ ] Parameterize eigenvalues as `D = diag(δ, 1, g)` and **calibrate (δ, g) numerically** against observed lepton-mass ratios; they are not derivable ab-initio per Jarek's 2026-04-19 guidance
- [ ] **Port Manfried Faber's regularization scheme** as the baseline (arxiv:2604.12021 + Universe 11/2025/113). Faber's form uses a slightly different potential but demonstrably produces the **running-coupling effect** — adapt to LdG-with-Skyrme rather than reinventing. This derisks the "hardest numerical step" from blank-slate design to port + adapt
- [ ] Implement the core-singularity handling carefully — soft core smoothing + adaptive time step near the defect, on top of Faber's regularization
- [ ] Validate: recover running-coupling effect (charge strength varies with distance / energy scale) as an independent check that the regularization is correctly ported
- [ ] Goal: three distinct lepton energy scales emerge from the calibrated `(δ, 1, g)` hierarchy plus the chosen LdG parameters `(a, b, c)`, with Faber's regularization active
- [ ] This is the "electron, muon, tau from biaxial geometry" experiment — a significant undertaking beyond initial M5

### Phase M5.7 — Cornell potential / quark confinement (new, per Jarek's guidance)

> **Refinement from Jarek Duda (2026-04-17)**: after M5.4 succeeds, the next validation target is **recreation of the [Cornell potential](https://en.wikipedia.org/wiki/Cornell_potential)** for quarks. Quarks are excitations of a **quark string / topological vortex**; fractional charges add, on top of Coulomb, a **linear ~1 GeV/fm confinement energy** from the conflict produced by fractional charges. This is the topological analog of QCD confinement.

- [ ] Seed a topological vortex string (1D defect line, not point hedgehog) connecting two fractional-charge end points
- [ ] Measure the interaction energy `V(r)` as a function of end-point separation
- [ ] Validate the Cornell form: `V(r) = −α/r + σ·r` with `σ ≈ 1 GeV/fm` (string tension)
- [ ] Compare to QCD phenomenology (linear confinement, asymptotic freedom at small r)
- [ ] **Why this matters**: it's the strong force / QCD analog in the Lagrangian-topological framework — demonstrating that the same ingredients that give lepton Coulomb (Exp 2) also give quark confinement when extended from point defects to string defects

### Phase M5.8 — De Broglie clock / Zitterbewegung test (new, per Jarek's guidance)

> **Refinement from Jarek Duda (2026-04-17)**: the 1+1D phi-4 kink toy-model (arxiv 2501.04036) validates the *mechanism*; for actual electron and neutrino, we need to run it in **full LdGS** (Landau-de Gennes + Skyrme) dynamics.
>
> **Conceptual background**: the full mechanism (why defects oscillate, why the oscillation is rotational not linear, how mass + spin + de Broglie wavelength all derive from one rotation) is documented in [3c_topological_defect.md](3c_topological_defect.md). This phase is the empirical-validation step that confirms M5's specific implementation reproduces the experimentally-established Zitterbewegung frequency.

- [ ] Seed a single LdGS defect (electron: non-dual hedgehog; neutrino: dual hedgehog / closed vortex loop)
- [ ] Let it evolve under full M5 dynamics (no external driving)
- [ ] Measure the intrinsic oscillation frequency of the defect core (FFT of director rotation at the core)
- [ ] Validate `ω = 2·m·c²/ℏ` (Zitterbewegung / de Broglie clock) — this is the mass-gap mechanism from Exp 4 extended to the full LdGS field
- [ ] Compare electron (SO(2) ~ U(1), 2D rotation → de Broglie clock) vs neutrino (SO(3) ~ SU(2), 3D rotation → neutrino oscillations)
- [ ] Confirm oscillation is **rotational** (director rotates around an axis) rather than linear position-bouncing (per [3c § What the oscillation looks like](3c_topological_defect.md#what-the-oscillation-looks-like--rotation-not-translation))
- [ ] **Mass → frequency table** (target values to validate; per [3c § Concrete particle table](3c_topological_defect.md#concrete-particle-table--masses-to-zitterbewegung-frequencies)):

| Particle | Defect type | Target ω = 2mc²/ℏ |
| --- | --- | --- |
| Electron | Point hedgehog (δ axis) | 1.55 × 10²¹ rad/s |
| Muon | Point hedgehog (1 axis) | 3.21 × 10²³ rad/s |
| Tau | Point hedgehog (g axis) | 5.39 × 10²⁴ rad/s |
| Neutrino | Closed vortex loop | ~10¹⁵ rad/s (sub-eV mass) |
| Up quark | Vortex string endpoint | ~7 × 10²¹ rad/s |
| Down quark | Vortex string endpoint | ~1.4 × 10²² rad/s |
| Proton (composite) | 3-string Y-config | 2.85 × 10²⁴ rad/s |

- [ ] **Cross-particle test**: seed defects of different masses (electron + muon, or electron + tau) at far separation. Measure each one's intrinsic frequency independently. Confirm each ticks at its own mass-derived `ω` (validates that frequency is set by *each* defect's stored energy, not a coupling between them)
- [ ] **Why this matters**: mass-driven oscillation is the origin of the wave-particle duality; validating it numerically in the full LdGS closes the loop from Exp 4 linear-order validation to full nonlinear particle dynamics. Once the table above is reproduced, M5 has empirically derived particle masses from geometric defect parameters — the Standard Model's free-mass-parameter problem becomes a calibration problem instead

### What M5 does NOT implement

From the sandbox findings, these are ruled out or de-prioritized:

- **The documented Combined W-L product form `2A·sin(kr/2)·cos(kr/2−(ωt+φ))/r`** — Exp 5 showed this isn't a free-wave solution. M5 uses the M4 *sum form* `A·[sin(kr+ωt+φ) + sin(kr−ωt−φ)]/(kr)` if analytical standing waves are needed
- **Smolinski's Ψ³ as primary K-selectivity mechanism** — Exp 8 falsified. Ψ³ may still be a stabilizer *inside* a topologically protected configuration, but M5 does not rely on it for geometric selectivity
- **Seeded-soliton emergence from spherical harmonics as a physics mechanism** — Exp 7 v2 confirmed Close's framework doesn't predict this. Particles in M5 = topological defects (from seed), not emergent solitons from wave seeds

---

## LAYERED VALIDATION ROADMAP — VACUUM TO ATOMS

**OpenWave's distinguishing value**: an **integrated simulator that demonstrates each layer of physics built from previously validated primitives** in a single platform — vacuum → leptons → mesons → nucleons → nuclei → atoms. Each layer depends on the layer below being already validated; cross-layer dependency checks catch inconsistencies between scales that single-layer simulators (QCD-only, atomic-only, liquid-crystal-only) cannot. This integrated cross-scale validation is what makes OpenWave unique.

The full physics hierarchy maps cleanly onto M5 phases:

| Layer | Physics primitive | M5 phase that validates it | Status |
| --- | --- | --- | --- |
| **0 — Vacuum** | LdG ground state, static field configuration | M5.0 (scaffold) + M5.1 (`seed_vacuum`) | 🚧 next |
| **1a — Single point defect (lepton-axis δ)** | Hedgehog as single-defect electron | M5.4 single-particle test (electron stability) | 🚧 |
| **1b — Lepton family hierarchy** | Three biaxial axes `(δ, 1, g)` → e/μ/τ from same Lagrangian | M5.6 (biaxial Q-tensor + Faber regularization) | 🚧 |
| **1c — Closed vortex loop** | Neutrino topology variants (SO(3)~SU(2)) | M5.6 alternative seed + M5.8 (de Broglie clock) | 🚧 |
| **1d — Two-defect interactions** | Dynamic Coulomb 1/d², annihilation | M5.4 pair test (e⁺/e⁻) | 🚧 |
| **1e — Intrinsic oscillation** | Zitterbewegung at `ω = 2mc²/ℏ` | M5.8 | 🚧 |
| **2 — Open vortex string + 2 endpoints** | Quark-antiquark (meson), Cornell potential `V(r) = −α/r + σ·r` with σ ≈ 1 GeV/fm | M5.7 (Cornell potential) | 🚧 |
| **3 — 3 string endpoints (baryon configuration)** | Proton (uud), neutron (udd) — color-neutral 3-quark composites; mass dominated by string energy | **Post-M5.8 (Phase 6 candidate)** | not yet planned |
| **4 — Multi-nucleon nucleus** | Nucleus binding via residual strong force | **Post-M5.8 (Phase 6+ candidate)** | not yet planned |
| **5 — Atom** | Z electrons in standing-wave orbital shells around nucleus, M3-style interference at atomic scale | **Post-M5.8 (Phase 6++ candidate)** | not yet planned |
| **6 — Molecules / bulk matter** | Multi-atom bonding | **Long-term** | not planned |

### Layer-by-layer dependency chain

Each higher layer requires the layer below to be working correctly:

- **Layer 1** (lepton) requires Layer 0 (vacuum that supports topological defects) — checked by M5.0/M5.1 invariant tests
- **Layer 2** (meson / vortex string) requires Layer 1 (defect machinery + force diagnostics) — checked by M5.7 Cornell-potential validation passing
- **Layer 3** (nucleon) requires Layer 2 (validated string-tension physics + 3-endpoint Y-configuration handling) — gated on M5.7 success
- **Layer 4** (nucleus) requires Layer 3 (working nucleon + residual-force model)
- **Layer 5** (atom) requires Layer 4 (working nucleus) + Layer 1 (working electron) + the standing-wave near-field physics (already validated in M3, retained in M5)

**The cross-layer integration is the unique value**: a successful Layer 5 (atom simulation) provides simultaneous validation that the lepton physics, string-string-tension physics, residual-force binding, and orbital-force standing-wave interference all work *together* — something a single-scale simulator cannot test by construction.

### Beyond M5.8 — composite-particle roadmap (sketched, not yet planned)

The current 9-phase plan (M5.0–M5.8) covers the vacuum, lepton, and meson layers (Layers 0–2). Layers 3–5 are major undertakings deferred to a future Phase 6+ scope. Sketched phases:

| Phase candidate | Layer | Headline |
| --- | --- | --- |
| **M6.1 / Nucleon assembly** | Layer 3 | Seed 3-string Y-configuration with color-neutral axis assignment; verify proton/neutron forms a stable bound state with mass dominated by string energy |
| **M6.2 / Color confinement test** | Layer 3 | Attempt to separate one quark from a 3-quark configuration; verify string tension prevents isolation (linear energy growth) |
| **M6.3 / Nuclear binding** | Layer 4 | Two-nucleon and few-nucleon bound states; measure residual strong force; reproduce binding-energy curve |
| **M6.4 / Atomic orbitals** | Layer 5 | Single electron orbiting a nucleus; verify discrete orbital shells emerge from standing-wave interference at atomic scales |
| **M6.5 / Multi-electron atom** | Layer 5 | Z electrons in a single-nucleus configuration; verify shell structure (Pauli-like exclusion via wave interference + topology) |

#### Cross-mass-class machinery required for M6.4+ (per [3c § How different-frequency Zitterbewegung emissions interfere](3c_topological_defect.md#how-different-frequency-zitterbewegung-emissions-interferehydrogen-vs-positronium))

Atom-scale simulations (M6.4 onward) need additional architectural capability beyond the same-mass-class physics validated by M5.4–M5.7. The reason: in hydrogen-like systems, the proton ticks at `ω_p ≈ 1836 × ω_e`, so direct e-p Zitterbewegung interference does NOT form coherent standing waves. Instead, three layered mechanisms must work together:

- [ ] **Topological 1/d Coulomb** between e (Q=−1) and p (effective Q=+1) — provides the static potential well (already validated in M5.4)
- [ ] **Electron self-de-Broglie standing waves** — the electron's Zitterbewegung (ω_e at 10²¹ rad/s) + translational drift produces de Broglie wavelength `λ_dB = h/(m_e v)`. Standing waves of the electron's *own* emission (interfering with reflections off the central proton's potential well) quantize discrete Bohr orbital shells at `n·λ_dB = 2π·r_n`
- [ ] **Quasi-static heavy-center treatment** for the proton — because m_p ≫ m_e, the proton's intrinsic clock (10²⁴ rad/s) is invisible to the electron's slower 10²¹ rad/s response and time-averages out. M6.4's solver must treat the proton as a fixed Coulomb center on the electron's time scale (with proton dynamics on its own slower scale only when needed)

This is structurally different from the same-mass cases: positronium (M5.4) and quark-quark binding (M5.7) DO use direct Zitterbewegung interference; hydrogen-like atoms (M6.4) cannot. The transition criterion is `ω_a ≈ ω_b` (frequency-matched, direct interference works) vs `|ω_a − ω_b| ≫ min(ω_a, ω_b)` (frequency-mismatched, requires the three-mechanism layering above).

This testing-architecture distinction is **gating** for atom-scale simulations: M6.4 cannot succeed by simply scaling up M5.4's mechanism. The cross-mass machinery is a separate engineering effort.

### Beyond M6 — thermal mechanics pathway (Phase 7+)

A working hypothesis for the heat output domain: **per-defect amplitude excess above the topological ground state is the thermal degree of freedom**. Each topological defect (= each particle in M5) has its own intrinsic ground-state oscillation amplitude `A_0 ≈ ℏ/(mc)` (its Compton wavelength); any excess above `A_0` is what aggregates statistically into macroscopic thermodynamic temperature. This extends thermodynamics from ensemble-only statistical mechanics (where heat is a collective property requiring many particles) to a per-particle quantum-mechanical degree of freedom — the same way M5 makes other things "more fundamental, not less" by deriving collective phenomena from defect-level dynamics.

The closest precursor framing in OpenWave is the steepness-conservation / energy-starvation work in [5_TIME_DYNAMICS.md](5_TIME_DYNAMICS.md). The topological-defect deep dive in [3c_topological_defect.md](3c_topological_defect.md) covers the defect's ground-state oscillation; this Phase 7 pathway is what extends it to the excited (thermally-excited) regime.

#### Phase 7 — Thermal mechanics from defect amplitude statistics (proposed)

| Phase candidate | Layer / domain | Headline test |
| --- | --- | --- |
| **M7.1 / Single-defect (A, ω) modulation** | Per-defect thermal | Take a single biaxial hedgehog from M5.4. **Bidirectional perturbation of both amplitude AND frequency, from the defect's *current* state** (not only from the topological ground state — a real defect's instantaneous state is its rest-energy ground `(A₀, ω₀)` *plus* current thermal excess; modulation perturbs from wherever it is now). Five sub-experiments share the same kink + Klein-Gordon background and instrumentation: **(7.1a) AM-up** — bump A above current; measure relaxation pathway, breather modes, emission spectrum. **(7.1b) AM-down** — pull A below current via destructive interference at the defect; measure whether the surrounding field refills the kink (the surroundings cool — the engineering operation). **(7.1c) FM-up** — drive ω above current via tuned external wave; measure A response per wave-steepness conservation `A/λ = const`. **(7.1d) FM-down** — slow ω below current via tuned interference; measure A rise and energy-source direction. **(7.1e) coupling cross-check** — validate `A/λ = const` holds at per-defect level under all four perturbation directions. The cheapest go/no-go check on the joint (A, ω) thermal hypothesis |
| **M7.2 / Soliton-breather comparison** | Per-defect thermal | Compare M7.1's measured amplitude oscillations to known soliton-breather modes (Sine-Gordon, φ⁴). If the framework is correct, these breather modes should match observed thermal excitation patterns. Cross-validation against established field-theory math |
| **M7.3 / Multi-defect amplitude statistics** | Defect ensemble thermal | Seed N defects (10², 10³, 10⁴) with varying initial amplitudes. Run to thermodynamic equilibrium. Extract amplitude distribution. Predict: should match Boltzmann (classical) or Bose-Einstein (quantized) statistics for the appropriate temperature definition |
| **M7.4 / Specific heat reproduction** | Defect ensemble thermal | From M7.3 statistics, derive specific heat C_V(T). Validate against experimental measurements: Dulong-Petit at high T, Einstein-Debye at low T, electronic heat-capacity scaling for free electrons. Stiff prediction — wrong scaling = hypothesis falsified |
| **M7.5 / Blackbody spectrum** | Thermal-EM coupling | Measure emission spectrum from a thermalized defect ensemble. Validate Wien's displacement law and Stefan-Boltzmann scaling. The heat → light channel — physics that connects the matter / forces / EM / heat output domains |
| **M7.6 / Phase-coherence transition** | Quantum thermal | At low T, do defect ensembles transition into phase-coherent ground states (analogous to superconductivity / BEC)? Reproduce critical temperatures for known materials. The cleanest "novel hypothesis" validation |

These Phase 7 phases are deliberately speculative — they assume the joint (A, ω) framework holds and cascade from there. If M7.1 falsifies the hypothesis, the rest don't run as planned. But if M7.1 confirms it, OpenWave has identified a new mechanism for thermal physics: heat as a single-particle quantum-mechanical phenomenon rather than an ensemble-only statistical one.

**M7.1 four-outcome decision matrix**:

| AM (7.1a + 7.1b) | FM (7.1c + 7.1d) + coupling (7.1e) | Interpretation |
| --- | --- | --- |
| ✅ | ✅ + `A/λ = const` validated | Joint (A, ω) hypothesis validated; proceed to M7.2 → M8 |
| ✅ | ❌ | Revise: heat is amplitude-only; engineering pivots to AM-pure pathway; doc updates needed |
| ❌ | ✅ | Revise: heat is frequency-only or coupling inverted; rethink the physics framing |
| ❌ | ❌ | Hypothesis falsified at the per-defect level; redirect Phase 7 effort |

Even if the framework needs refinement (first hypotheses rarely survive intact), the *framework for asking the question* — per-defect amplitude as a thermal degree of freedom alongside topology and mass — is what M5+ matter physics enables. **No other framework currently lets us pose the question this way.** That's the deeper value: even partial validation moves thermal physics from "ensemble-only statistical" to "single-defect quantum-mechanical".

These plans are kept intentionally rough — concrete phase work waits until M5.0–M5.8 establish the foundations.

### Why the hierarchy matters for the project's positioning

This integrated layered scope is what differentiates OpenWave from comparable simulators:

| Simulator class | What it covers | What it misses |
| --- | --- | --- |
| **QCD lattice** (mainstream particle physics) | Quark/gluon dynamics, confinement, hadron masses | Leptons, atoms, vacuum-as-LdG-medium, classical interpretability |
| **Liquid-crystal solvers** (LdG mathematics) | Topological defects, hedgehogs, disclinations | Particle-physics calibration, atom-scale dynamics, lepton families |
| **QED / atomic simulators** | Atoms, molecules, EM | Quark structure, lepton ⊃ defect identity, vacuum mechanics |
| **OpenWave (M5 → M6)** | All of the above, **integrated, classical-field-theoretic** | Calibration to high-precision data still in progress (validation against experimental ratios) |

The integration is the value. Reading the layered table top-to-bottom is also a reading-order suggestion for new contributors: validated primitives at the top, frontier work at the bottom.

For the conceptual companion to this roadmap, see [3b_concept_review.md § Where do quarks, protons, nuclei, and atoms fit?](3b_concept_review.md#where-do-quarks-protons-nuclei-and-atoms-fit) — same hierarchy, framed as a Q&A explanation rather than an implementation checklist.

### Beyond M5.8 — phase, Berry, and entanglement experiments

A separate research direction prompted by the 2026-04 Models of Particles thread on the Orion–Akkermans paper *"Topological sum rule for geometric phases of quantum gates"* (arxiv:2603.29795). The paper's headline corollary is that nontrivial Hamiltonian topology (`ν_H ≠ 0`) is a **necessary condition for quantum entanglement** — and M5's defect framework satisfies it by construction (every defect has nonzero winding). This opens up a class of experiments that test whether M5's twist degree of freedom (see [3c_topological_defect.md § The twist degree of freedom](3c_topological_defect.md#the-twist-degree-of-freedom--quantum-phase-as-a-derived-field-state)) reproduces the geometric-phase / entanglement structure of standard QM.

These are exploratory targets, not committed milestones. They naturally slot in after M5.4 (electron stability) provides the validated single-defect substrate and M5.7 (Cornell potential, vortex strings) provides validated string-defect dynamics.

| Phase candidate | Headline test | Layered on |
| --- | --- | --- |
| **Single-defect Berry phase** | Drive a single defect adiabatically along a closed loop in parameter space (e.g., rotate the LdG axis assignment). Measure the accumulated twist phase around the loop. Should equal the geometric solid angle × topological factor — the Berry phase formula falls out of the field's own dynamics, not as a postulate | M5.4 (single biaxial hedgehog stable) + M5.6 (axis-rotation machinery) |
| **Hydrodynamic Aharonov-Bohm analog** | Seed a topological vortex string (M5.7 setup); send a plane wave past it; measure the phase shift on the far side as a function of impact parameter. Should reproduce the Berry-1980 hydrodynamic AB result and the recent 2026 reproduction (phys.org/news/2026-04-simulation-famous-quantum-effect-reveals.html) | M5.7 (vortex string seeding) |
| **Two-defect concurrence analog** | Seed a positronium-like e⁺/e⁻ pair (M5.4 setup). Vary the initial conditions over a complete basis of "computational states." Measure the redistribution of geometric phases across initial states. Compare to the Orion sum rule prediction `(1/2π) Σ γ_n = 2m·ν_H` and the Wootters-concurrence-vs-phase pattern from the paper | M5.4 (e⁺/e⁻ pair test) |
| **Topological sum rule test** | Implement the same "gate" (e.g. defect-pair SWAP-analog) via two M5 evolution paths from different topological classes. Measure the per-state Berry phases; verify the sum is the same (`2m·ν_H` invariant) but the redistribution differs. Direct analog of the SWAP₁ vs SWAP₂ comparison in the paper's Fig. 3 | M5.4 + post-M5 (gate-analog protocols defined) |
| **CPT symmetry validation** | M5's Lagrangian is Hermitian and Lorentz-invariant by construction → must satisfy CPT (theorem). Verify numerically that pair creation and annihilation operate as CPT mirrors, that same-winding pairs are CPT-forbidden from annihilating (already a planned M5.4 test, here recognized as a CPT consequence), and that absorption / stimulated-emission analogs are CPT-paired (per Duda's 2026-04 deck on negative radiation pressure). Failure mode = a hidden T-symmetry-breaking term in the implementation | M5.4 (pair dynamics) + M7.x (thermal absorption / stimulated emission) |

#### Why this matters for OpenWave

These experiments do three things at once:

- **Validate that M5's twist DoF correctly carries quantum-phase physics.** If the Berry phase falls out of the M5 dynamics with the standard geometric-phase formula, that's strong evidence the framework is structurally complete: charge (topology), mass (amplitude), AND quantum phase (twist) all derived from one Lagrangian
- **Engage productively with the topology-and-entanglement literature.** The Orion paper sits in mainstream condensed-matter physics; M5 reproducing its predictions from a different starting point (classical-field-theoretic, topological defects) is the kind of cross-validation that makes the framework credible to physicists outside the wave-mechanics community
- **Connect to hydrodynamic quantum analogs.** Berry's hydrodynamic Aharonov-Bohm (1980, recreated 2026) and the broader walking-droplet HQA literature have struggled with multi-particle entanglement as the missing piece. M5's defect framework — where multi-defect entanglement-like correlations are automatic from shared topology — may be the structural bridge. If validated, this connects HQA, EWT, LdG-topology, and standard QM into one coherent picture

#### Open questions

- Does the twist field in M5 produce **exactly 2π-quantized** Berry phase, or only approximately? (Quantization should be geometric, but needs numerical confirmation)
- Does the two-defect Wootters concurrence analog quantitatively match the paper's predictions for entangling Hamiltonians, or only the qualitative pattern?
- How does the result depend on the choice of LdG potential, regularization (Faber), and Skyrme coefficient? Are there parameter regimes where the topology-entanglement link breaks?
- Do the standing-wave + topology + twist primitives also reproduce the **double-slit** and **Casimir** patterns Marc Fleury cites — i.e., is the framework consistent with the broader EM-fluid / standing-wave-electron literature beyond just the Orion paper?
- For Bell-inequality experimental analogs: does the twist-field path-dependence reproduce the experimental phase-compensation patterns, supporting Marc's claim that "Bell completely disappears with phase modulation"?

These are deferred research questions — concrete numerical work waits until M5.4 + M5.7 establish the defect substrate. Documented here so the framework is in scope; specific experiments will be planned when the prerequisites land.

### Beyond matter — forces, EM waves, heat

The matter-particle hierarchy above (Layers 0–6) is OpenWave's *foundational layer*. On top of validated matter primitives, the simulator's full scope covers four output domains:

| Output domain | What M5+ must compute | Where it sits in the roadmap |
| --- | --- | --- |
| **MATTER** (foundation) | Particle emergence (leptons → quarks → nucleons → atoms) via topological defects + wave dynamics | Layers 0–5+ above |
| **FORCES** | Electric (topology), strong (string tension + standing waves), magnetic (transverse waves from spin), gravitational (density deficit / 4D boost-axis topology) — all derived from a single Lagrangian | Validated cumulatively across M5.4 (electric), M5.7 (strong), Phase 4 (magnetic), Phase 5 (gravity) |
| **ELECTROMAGNETIC WAVES** | Photon emission/absorption, EM dispersion, polarization, antenna-medium coupling | Phase 4 (alongside magnetic-force validation) |
| **HEAT** | Thermal energy at the wave / spin-coherence level (not bulk kinetic temperature); thermal-EM coupling; Wien's law / blackbody emergence | Phase 6+ (long-term, after the matter layers are stable) |

**Why these four are co-equal scientific goals**:

- The matter hierarchy alone is incomplete physics. You can simulate an electron beautifully and still not understand how it interacts with surrounding fields, photons, or thermal environments
- Force / EM / heat are how matter *interacts* with the rest of the world — these are the observable, measurable outputs against which a wave-mechanics simulator is validated
- A simulator that does only matter is a particle-physics tool. A simulator that does matter + forces + EM + heat from one Lagrangian is a **complete classical-field-theoretic description** — that's the integrated value of OpenWave

### How matter layers feed forces / EM / heat outputs

Each output domain depends on validated matter layers being in place. The dependency chain:

| Output | Depends on matter Layer(s) | M5 phase chain |
| --- | --- | --- |
| Electric force (Coulomb 1/d²) | Layer 1 (single defects) | M5.1 + M5.4 |
| Strong force / confinement | Layer 1 + Layer 2 (vortex strings) | M5.1 + M5.7 |
| Magnetic force | Layer 1 + spin dynamics | Phase 4 (post-M5) |
| Gravitational force | Layer 1 + density-deficit / 4D extension | Phase 5+ (post-M5) |
| Photon / EM wave dynamics | Layer 1 + tilt-axis curl (Maxwell from director field) | Phase 4 (alongside magnetic) |
| Thermal energy / heat | Layer 5 (atoms) + spin-coherence dynamics | Phase 6+ (post-M5, post-atoms) |

Heat is *deepest* in the dependency chain because it requires simulated atoms to model thermal phenomena meaningfully (temperature is a statistical property of many-particle systems).

### Applied-technology counterpart

The applied-technology counterpart of OpenWave's open-source physics work is the SABER project (private, separate repo at `OPENWAVE-LABS/SABER` / `neptunyalabs/SABER`). The split is intentional: OpenWave publishes scientific findings (open-source, eventually peer-reviewable); SABER develops engineering implementations (private, patentable). SABER consumes the physics outputs of M5+ matter / forces / EM / heat work. **OpenWave's scope is the science; SABER's scope is the engineering.** This document is OpenWave's roadmap and stops at the physics.

---

## STATUS

- [x] Architecture analysis complete (this document)
- [x] **All 8 sandbox experiments complete** — see [3a_lagrangian_experiments.md](3a_lagrangian_experiments.md):
  - ✅ Exp 1 (Sine-Gordon kinks), ✅ Exp 2 (Hedgehog Coulomb), ✅ Exp 3 (Winding quantization), ✅ Exp 4 (Klein-Gordon dispersion)
  - ⚠️ Exp 5 (Lagrangian derivation — W-L product form falsified), ⚠️ Exp 6 (lepton mechanism; specific ratios deferred), ⚠️ Exp 7 v2 (Close's actual equations implemented)
  - ❌ Exp 8 (Smolinski Ψ³ K-selectivity falsified)
- [x] **Winning recipe identified**: topology + Klein-Gordon + Close's Eq. 23 + M3 near-field + Skyrme stabilizer
- [x] **Group feedback integrated (2026-04-19)** — Jarek, Jeff, and Robert reviewed the sandbox summary; refinements captured in this document (Eq. 23 over Eq. 21, axis-hierarchy for lepton masses, Cornell potential and de Broglie clock added as M5.7/M5.8 targets, resonance-lifetime success criterion)
- [ ] M5.0 — Scaffold (Taichi structure, triple buffer, Laplacian, curl/div operators)
- [ ] M5.1 — Port topology from Exps 2, 3 (`seed_vacuum`, `seed_hedgehog`, Frank energy, winding tracker)
- [ ] M5.2 — Wave dynamics from **Close's Eq. 23** (with `∇·s = 0` enforced) + Klein-Gordon mass term, validate Exp 4 dispersion, amplitude-sweep resonance hunt
- [ ] M5.3 — Hamiltonian energy (replaces postulated `E = ρV(fA)²`)
- [ ] M5.4 — **Headline test**: single biaxial hedgehog (electron) is a long-lived resonance; hedgehog + anti-hedgehog pair reproduces dynamic 1/d² Coulomb + annihilation
- [ ] M5.5 — Skyrme stabilizer (conditional on M5.4)
- [ ] M5.6 — Biaxial LdG Q-tensor with `(δ, 1, g)` hierarchy (long-term; for lepton mass derivation)
- [ ] M5.7 — Cornell potential / quark confinement (topological vortex string, `V(r) = −α/r + σ·r`)
- [ ] M5.8 — De Broglie clock / Zitterbewegung test (`ω = 2mc²/ℏ`) for electron + neutrino

**Next action**: begin **M5.0 scaffold**. The sandbox is done; the recipe is known; group feedback is integrated. Move to production implementation on the Taichi engine.

---

## GROUP FEEDBACK (2026-04-17/18) — REFINEMENTS TO M5 PLAN

Replies from Jarek Duda, Jeff Yee, and Robert Close to the Apr 17 sandbox-summary email produced four targeted refinements to the M5 plan. Full email thread in [3_LAGRANGIAN_FRAMEWORK.md § EMAIL THREAD](3_LAGRANGIAN_FRAMEWORK.md#email-thread).

### Refinement 1 — Use Close's Eq. 23 as the particle equation (Robert Close)

Exp 7 implemented Close's Eq. 21 ("Equation of Everything"). Close's recommendation for the **particle equation** is actually **Eq. 23** — it preserves `∇·s = 0` (zero divergence of spin density). M5.2 adopts Eq. 23 as the particle equation, with Eq. 19 as the V=0 linear limit, and enforces the divergence constraint per time step.

### Refinement 2 — Particles are resonances, not stable solitons (Robert Close)

Close: *"Unless you have a good way to model an infinite system, I doubt that you will find completely stable non-radiating solutions."* M5's success criterion is reframed: **long-lived resonance with measurable lifetime**, not perfect stability. M5.2 includes an amplitude-sweep protocol (`Y_1^0` seed at `A ∈ {0.25, 0.5, 1.0, 2.0}·λ`) to hunt for these resonances.

### Refinement 3 — Lepton mass hierarchy from `0 < δ << 1 << g` axis lengths (Jarek Duda)

The ~3477:1 biaxial ratio the sandbox required is **physically motivated**, not ad-hoc. Three naturally-separated physical scales map to the three axis eigenvalues: `δ ~ ℏ` (QM / twist), `1` (unity / tilt / matter), `g` (gravity / boost). M5.6 parameterizes the Q-tensor this way from the start. Jarek also flagged that **LdG regularization is the hardest numerical step** — a known challenge that M5.6 budgets accordingly.

### Refinement 4 — New validation phases M5.7 and M5.8 (Jarek Duda)

- **M5.7 — Cornell potential**: after M5.4 succeeds, extend from point hedgehogs to a **topological vortex string** with fractional-charge end points. Measure `V(r) = −α/r + σ·r` and confirm the `~1 GeV/fm` linear confinement coefficient. This is the QCD confinement analog in the Lagrangian-topological framework
- **M5.8 — De Broglie clock / Zitterbewegung**: measure the intrinsic oscillation frequency of a single LdGS defect (electron: non-dual; neutrino: dual) and validate `ω = 2·m·c²/ℏ`. The 1+1D phi-4 toy (arxiv:2501.04036) showed the *mechanism*; M5.8 validates it in full LdGS for real particles

### Refinement 5 — M3 near-field carries three force regimes (Jeff Yee)

Jeff confirmed the M3 + topology coexistence decision and added a crucial scope expansion: M3's standing-wave lock-in is the mechanism for **three** force regimes, not just one — (a) intra-particle binding of K=1 WCs into standalone particles, (b) intra-nucleus strong force between K=10 particles, and (c) **orbital force** (electron-nucleus binding in atoms). This strengthens the rationale for keeping M3 physics intact in M5 — it's load-bearing well beyond the K=10 electron problem.

> **Note on vocabulary (2026-04-19)**: Jeff's regime enumeration uses EWT's K-count vocabulary (K=1 = neutrino, K=10 = electron) as was standard before the Duda framework was adopted. In the M5 paradigm, the mapping shifts: (a) single-particle stability is *not* a wave lock-in problem in M5 (the particle is a single defect, held together by topology, not by K=1-to-K=10 wave interference); (b) intra-nucleus strong force is between *quarks* (single defects of a different type) bound into nucleons, not between electrons; (c) orbital force unchanged — electrons around nuclei. Jeff's underlying point stands: M3's near-field physics is load-bearing for composite-particle and atomic-scale dynamics in M5.

### Summary impact table

| Area | Change | Source |
| --- | --- | --- |
| M5.2 Wave dynamics | Close's **Eq. 23** (not Eq. 21) as particle equation; enforce `∇·s = 0`; amplitude-sweep resonance hunt | Robert |
| M5.4 Headline test | Success criterion = long-lived resonance (lifetime ratio), not perfect stability | Robert |
| M5.6 Biaxial LdG | Axis hierarchy `0 < δ << 1 << g` (physics-motivated); LdG regularization budgeted as hardest step | Jarek |
| **New M5.7** | Cornell potential via topological vortex string; `V(r) = −α/r + σ·r` with `σ ≈ 1 GeV/fm` | Jarek |
| **New M5.8** | De Broglie clock / Zitterbewegung test in full LdGS for electron + neutrino | Jarek |
| M3 retention rationale | M3 near-field covers *three* force regimes: intra-particle, strong, and orbital | Jeff |
| M5.6 `(δ, g)` treatment | Calibrated numerically (no ab-initio form exists); iterate against observed lepton-mass ratios | Jarek (2026-04-19) |
| M5.6 regularization baseline | Port Manfried Faber's scheme (arxiv:2604.12021, Universe 11/2025/113) — slightly different potential but produces running-coupling effect; adapt rather than reinvent | Jarek (2026-04-19) |
