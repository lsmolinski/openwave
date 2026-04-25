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

### Phase M5.0 — Scaffold

- [ ] Create `openwave/xperiments/m5_lagrangian_field/` directory (mirror m4 structure)
- [ ] **Rename the engine module**: `wave_engine.py` → `lagrangian_engine.py` for M5. Rationale: M5's core loop integrates a Lagrangian-derived PDE (`∂²_tψ = c²∇²ψ − ∂V/∂ψ`) that simultaneously handles wave propagation *and* preserves topology via the potential `V(ψ)`. "Wave" is only one of the two channels the engine produces, so `lagrangian_engine.py` reflects what the module actually is to a new reader. M1–M4 keep `wave_engine.py` (they really are wave engines, no topology layer). See [3b § What wave equation does M5 solve?](3b_concept_review.md#what-wave-equation-does-m5-solve-is-force-still-e) for the reasoning
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
- [ ] Test annihilation: when the pair closes to near-zero separation, do the winding numbers cancel and the stored energy radiate outward as Klein-Gordon waves?

**Success metrics**:

- [ ] Single biaxial hedgehog (electron-axis) is a long-lived resonance under moderate perturbation (lifetime substantially longer than decoy variants)
- [ ] Hedgehog + anti-hedgehog pair dynamics reproduce 1/d² Coulomb attraction + annihilation
- [ ] **This is the full Phase 2 → Phase 3 → M5 validation loop.** If the single biaxial hedgehog is a long-lived electron-like resonance AND the e⁺/e⁻ pair reproduces Coulomb + annihilation dynamically, the loop is closed

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

- [ ] Seed a single LdGS defect (electron: non-dual hedgehog; neutrino: dual hedgehog)
- [ ] Let it evolve under full M5 dynamics (no external driving)
- [ ] Measure the intrinsic oscillation frequency of the defect core
- [ ] Validate `ω = 2·m·c²/ℏ` (Zitterbewegung / de Broglie clock) — this is the mass-gap mechanism from Exp 4 extended to the full LdGS field
- [ ] Compare electron (SO(2) ~ U(1), 2D rotation → de Broglie clock) vs neutrino (SO(3) ~ SU(2), 3D rotation → neutrino oscillations)
- [ ] **Why this matters**: mass-driven oscillation is the origin of the wave-particle duality; validating it numerically in the full LdGS closes the loop from Exp 4 linear-order validation to full nonlinear particle dynamics

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

These are kept intentionally rough — concrete phase plans wait until M5.0–M5.8 establish the foundations.

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

### Beyond matter — forces, EM waves, heat (the SABER design-parameter outputs)

The matter-particle hierarchy above (Layers 0–6) is OpenWave's *foundational layer*. On top of validated matter primitives, the simulator must also compute four additional output classes that are the actual deliverables for downstream applied technology:

| Output domain | What M5 must compute | Where it sits in the M5 plan | Why SABER needs it |
| --- | --- | --- | --- |
| **MATTER** (foundation) | Particle emergence (leptons → quarks → nucleons → atoms) via topological defects + wave dynamics | Layers 0–5+ above | Provides the substrate every other output is measured against |
| **FORCES** | Electric (topology), strong (string tension + standing waves), magnetic (transverse waves from spin), gravitational (density deficit / 4D boost-axis topology) — all derived from a single Lagrangian | Validated cumulatively across M5.4 (electric), M5.7 (strong), Phase 4 (magnetic), Phase 5 (gravity) | SABER's energy-conversion designs need quantitative force-law predictions at the device scale |
| **ELECTROMAGNETIC WAVES** | Photon emission/absorption, EM dispersion, polarization, antenna-medium coupling | Phase 4 (alongside magnetic-force validation) | RF/microwave/optical interaction with media is core to multiple SABER mechanisms |
| **HEAT** | Thermal energy at the wave / spin-coherence level (not bulk kinetic temperature); thermal-EM coupling; Wien's law / blackbody emergence | Phase 6 (long-term, after the matter layers are stable) | SABER's flagship goal — ocean-heat → electricity conversion — requires first-principles thermal mechanics |

**Why these four are co-equal goals, not just "nice to have"**:

- The matter hierarchy alone is not actionable for tech. You can simulate an electron beautifully and learn nothing about how to build a device
- Force / EM / heat are what *interfaces* matter with engineering — every device exploits one or more of these to do useful work
- A simulator that does only matter is a particle-physics tool. A simulator that does matter + forces + EM + heat from one Lagrangian is a **design-parameter generator** — and that's the unique value proposition for the integrated SABER programme

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

### SABER as the applied-technology counterpart

OpenWave (this repo) is the open-source scientific simulator. **SABER** (SubAtomic Based Energy Research) is the applied-technology programme that consumes OpenWave's predictions and engineers them into devices. The split:

- **OpenWave**: classical-field theory + topology + waves + heat → first-principles physics outputs (forces, EM, heat tables; design parameters; validated invariants)
- **SABER**: device design + experimental validation + commercialization (most prominently ocean-heat → electricity conversion, but also other tech using subatomic-scale predictions)

The dependency runs OpenWave → SABER. SABER cannot proceed past prototypes without the design parameters M5+ produces. M5's matter physics is therefore the gating step for SABER's longer-term applied work.

This is why M5's roadmap goes beyond "simulate the electron well" — it must reach all the way to producing actionable force-EM-heat outputs at engineering scales.

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
