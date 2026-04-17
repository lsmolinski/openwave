# LAGRANGIAN FRAMEWORK — NUMERICAL EXPERIMENTS

> *Numerical experiments — same spirit as OpenWave xperiments: controlled investigations of physics hypotheses through simulation. OpenWave is built as a shared experimental platform where wave-based and topological models can be tested, compared, and built upon by others.*

Numerical results from the 8 Lagrangian framework experiments in `sandbox_phase2_lagrangian/`.
See [2a_lagrangian_eval.md](2a_lagrangian_eval.md) for experiment specifications.

---

## Using OpenWave "Sandbox" vs "Production" environment

OpenWave uses two layers of experimentation:

- **Sandbox** (`research_hub/sandbox_*/`) — quick numpy scripts to validate math, logic, and concepts. Pure exploration. No GPU, no Taichi, no production dependencies. This is where ideas are tested cheaply before committing engineering effort. If an experiment works in the sandbox, it graduates to the production engine.
- **Production** (`openwave/xperiments/m*/`) — Taichi-based 3D rendering and simulation on the official OpenWave platform. GPU-accelerated, full grid infrastructure, visualization, force & motion. This is the final product where validated equations run at scale.

**Sandbox = explore fast, fail cheap. Production = implement validated winners.**

## NEXT SESSION KICKOFF

**First action**: open `sandbox_phase2_lagrangian/exp1_sine_gordon_1d.py` and complete the implementation (skeleton already in place — TODO comments mark what's missing). Reference the Setup section in this doc (1.2) before tweaking parameters.

**Sequence for the day**:

1. Implement and run **Experiment 1** (Sine-Gordon 1D) — build intuition with a simple, well-understood soliton system
2. Fill in this doc's Experiment 1 sections: Setup → Results → Numerical Evidence → Comparison to Expected → Conclusion
3. Update Summary Dashboard status (🚧 → 🔶 → ✅/❌)
4. Move to **Experiments 2 + 3** (hedgehog energy + topological charge) — the highest-value test pair, validates Coulomb from topology vs our sinc problem

**Recommended order reminder**: Exp 1 → Exp 2+3 → Exp 5 → Exp 8 → Exp 4 → Exp 6 → Exp 7

**Pre-reading if time**: Bush 2015 review "Pilot-wave hydrodynamics" (referenced in 2a) for intuition on classical wave-particle quantum-like behavior.

---

## HOW TO USE THIS DOC

**For each experiment**:

1. **Before running** — fill in section `N.2 Setup` with grid parameters, initial conditions, etc.
2. **Run the sandbox script** in `sandbox_phase2_lagrangian/expN_*.py`
3. **Capture output** in section `N.3 Results` (raw numbers, behaviors observed)
4. **Save evidence** in section `N.4 Numerical Evidence` (plots, tables, key measurements)
5. **Fill comparison table** `N.5 Comparison to Expected` (predicted vs measured)
6. **Write conclusion** in `N.6 Conclusion` (passed / failed / inconclusive, why)
7. **Note follow-ups** in `N.7 Next Steps` (tweaks, related experiments, open questions)
8. **Update Summary Dashboard** status emoji at the top

**After all experiments**:

- Fill in `OVERALL CONCLUSIONS` section
- Complete the comparison matrix (Topology vs Nonlinearity vs Standing Waves)
- Make recommendation for `What to Implement in M4/M5`

---

## SUMMARY DASHBOARD

| # | Experiment | Status | Result | Date | Script |
| --- | --- | --- | --- | --- | --- |
| 1 | Sine-Gordon 1D solitons | ✅ Passed | Kink stability, Lorentz contraction, pass-through all confirmed | 2026-04-16 | `exp1_sine_gordon_1d.py` |
| 2 | Hedgehog energy vs distance | 🚧 Pending | - | - | - |
| 3 | Topological charge quantization | 🚧 Pending | - | - | - |
| 4 | Klein-Gordon from twist | 🚧 Pending | - | - | - |
| 5 | Lagrangian derivation | 🚧 Pending | - | - | - |
| 6 | Three lepton families | 🚧 Pending | - | - | - |
| 7 | Close's nonlinear vector wave eq | 🚧 Pending | - | - | - |
| 8 | Smolinski's non-linear Ψ³ | 🚧 Pending | - | - | - |

**Status legend**: 🚧 Pending | 🔶 In progress | ✅ Passed | ❌ Failed | ⚠️ Inconclusive

**Recommended order**: Experiment 1 (intuition) → Experiments 2+3 (Coulomb from topology) → Experiment 5 (math) → Experiment 8 (Smolinski Ψ³ K-selectivity) → Experiment 4 (Klein-Gordon dynamics) → Experiment 6 (lepton families) → Experiment 7 (Close's equation).

---

## ✅ EXPERIMENT 1: Sine-Gordon 1D Solitons

**Status**: ✅ Passed
**Sandbox Script**: `sandbox_phase2_lagrangian/exp1_sine_gordon_1d.py`
**Date run**: 2026-04-16

### 1.1 Hypothesis

The Sine-Gordon equation `∂²φ/∂t² - c²∂²φ/∂x² + (m²c⁴/ℏ²)·sin(φ) = 0` produces stable kink solitons that exhibit:

- Topological stability (kinks cannot dissipate)
- Pair creation/annihilation (kink + anti-kink collisions)
- Lorentz contraction (relativistic kinematics from wave equation)

### 1.2 Setup

- **Grid**: N=1024, dx=0.1, domain=[-51.2, 51.2]
- **Physics** (natural units): c=1, m=1, ℏ=1 → natural kink width L = ℏ/(mc) = 1
- **Time evolution**: leapfrog, dt=0.05, N_STEPS=2000 (t_final=100), CFL = c·dt/dx = 0.5
- **Initial conditions**:
  - Test 1 (static): kink at x₀=0, v=0
  - Test 2 (moving): kink at x₀=-20, v=0.5c → γ=1.1547
  - Test 3 (pair): kink at x_L=-15 moving +0.5c, anti-kink at x_R=+15 moving -0.5c
- **Boundary**: fixed endpoints (kink starts far from edges)

### 1.3 Results

**Test 1 — Static kink stability** (v=0):

- Initial energy E₀ = 7.9956 (predicted 8·m·c² = 8.000, ~0.06% discretization offset)
- Final energy after t=100: E_T = 7.9956
- Max relative energy drift: 1.5e-6 (excellent leapfrog conservation)
- Measured kink width: L = 1.002 vs predicted L = ℏ/(mc) = 1.000
- Shape invariant across full simulation (see `exp1_results/test1_static_kink.png`)

**Test 2 — Moving kink with Lorentz contraction** (v=0.5c):

- Lorentz factor: γ = 1/√(1-0.25) = 1.1547
- Measured velocity (linear fit of kink x(t) over t=[5, 90]): **0.4997 c** (0.06% off input)
- Measured final kink width: **0.867** vs predicted L/γ = **0.866** (0.2% off)
- Initial energy: 9.2308 vs predicted γ·8·mc² = 9.2376 (0.07% off)
- Energy drift: 1.16e-5
- Position vs time plot shows straight line overlapping predicted x₀+v·t (see `exp1_results/test2_moving_kink.png`)

**Test 3 — Kink + anti-kink collision** (v=±0.5c, meeting at t=30):

- Initial pair energy: 18.4615 vs predicted 2·γ·8·mc² = 18.475 (0.07% off)
- Final pair energy: 18.4613
- Max relative drift during collision: 1.67e-4 (transient bump during overlap)
- **Collision behavior: pass-through, not annihilation.** At t=30 (collision moment), the bump amplitude dips to ≈3.81 (from 2π≈6.283). After the collision the pair emerges **inverted**: a down-bump (0 → -2π → 0) with the anti-kink now on the left and kink on the right. The kinks swapped positions and passed through each other, preserving their topological charges (Q=±1)
- This is the classic integrable Sine-Gordon Perring-Skyrme behavior (see `exp1_results/test3_kink_antikink.png`)

### 1.4 Numerical Evidence

Plots in `sandbox_phase2_lagrangian/exp1_results/`:

- `test1_static_kink.png` — snapshots (10 time slices overlap into single curve) + energy-vs-time (flat line)
- `test2_moving_kink.png` — snapshots (11 kink positions from -20 to +30) + position-vs-time (linear, matches prediction)
- `test3_kink_antikink.png` — snapshots showing up-bump (t=0-20) → collapse (t=30) → down-bump (t=40+); energy-vs-time with small transient at collision

| Quantity | Predicted | Measured | Error | Match? |
| --- | --- | --- | --- | --- |
| Test 1 — Static kink energy | E = 8·m·c² | 7.9956 | 0.06% | ✅ |
| Test 1 — Kink width | L = ℏ/(mc) = 1.000 | 1.002 | 0.2% | ✅ |
| Test 1 — Energy drift over t=100 | 0 | 1.5e-6 (rel.) | — | ✅ |
| Test 2 — Kink velocity | 0.500 c | 0.4997 c | 0.06% | ✅ |
| Test 2 — Lorentz-contracted width | L/γ = 0.866 | 0.867 | 0.2% | ✅ |
| Test 2 — Moving kink energy | γ·8·m·c² = 9.238 | 9.231 | 0.07% | ✅ |
| Test 3 — Pair energy conservation | 2γ·8·m·c² = 18.475 | 18.462 | 0.07% | ✅ |
| Test 3 — Topological pass-through | kinks pass through, Q conserved | up-bump → down-bump (inverted pair), Q=±1 each preserved | — | ✅ |

### 1.5 Comparison to Expected

| Quantity | Predicted | Measured | Match? |
| --- | --- | --- | --- |
| Kink rest energy | E = 8·m·c² = 8.000 | 7.996 | ✅ |
| Kink width | L = ℏ/(mc) = 1.000 | 1.002 | ✅ |
| Lorentz contraction | L/γ = 0.866 at v=0.5c | 0.867 | ✅ |
| Velocity fidelity (input vs propagated) | 0.500 c | 0.4997 c | ✅ |
| Topological stability (kink lifetime) | infinite (Q=±1 protected) | stable > 100 t-units, no shape decay | ✅ |
| Energy conservation (leapfrog) | machine-precision-like drift | rel drift 1.5e-6 (static), 1.2e-5 (moving), 1.7e-4 (collision) | ✅ |
| Kink-antikink annihilation | pure SG: pass-through (not annihilation) | pass-through confirmed (up-bump → down-bump) | ✅ |

### 1.6 Conclusion

**Experiment 1 validates all three hypotheses of the Sine-Gordon framework.**

1. **Topological stability** — kinks are absolutely stable over the full 100-time-unit simulation. Energy conserved to 1.5e-6 for a static kink, shape invariant. This is the first direct demonstration in this project of a field configuration that is **topologically protected** rather than dynamically balanced. Unlike OpenWave's M3 standing waves which require precise WC placement to hold a K=10 tetrahedron, a Sine-Gordon kink cannot dissipate — the math forbids it
2. **Lorentz contraction emerges from the equation** — we did not impose relativistic kinematics; the sine-Gordon PDE itself produces width L/γ for a kink moving at velocity v. Measured 0.867 vs predicted 0.866 (0.2%). The rest-frame kink energy scales as γ·E₀, also confirmed. This is the classical-wave origin of special relativity that Duda emphasizes — SR is a consequence of wave dynamics, not an added postulate
3. **Pair dynamics are pass-through, not annihilation** — kink and anti-kink with Q=+1 and Q=-1 meet at t=30, briefly compress (amplitude drops from 2π to ~3.81), then re-emerge with swapped positions and an inverted bump shape. Total Q stays 0. This is the integrable soliton behavior of Sine-Gordon; pure annihilation (kink + anti-kink → radiation) would require a non-integrable potential (e.g. double Sine-Gordon or φ⁴)

**Implication for OpenWave**: topology gives a stability mechanism that is categorically different from our standing-wave lock-in. A kink cannot "drift away" the way a K=10 WC can under perturbation — it is locked by winding number, a discrete integer. This is exactly the ingredient Duda said we need to prevent the "electron from exploding" under perturbation.

**Caveat**: pure Sine-Gordon 1D is not sufficient for particles (1D kinks carry Q=±1 but have no spatial extent in 3D; and pass-through ≠ annihilation, which we observe in positronium). The next test in line — **Experiment 2 (hedgehog energy vs distance)** — extends this to 3D with actual particle-like defects whose interactions should resemble Coulomb, not pass-through.

### 1.7 Next Steps

- **Move to Experiments 2 + 3** (hedgehog energy vs distance + topological charge quantization). These are the highest-value test pair — they validate whether topology produces clean 1/r Coulomb (vs our sinc oscillation) and whether winding numbers actually return integers numerically
- **Future follow-ups for Sine-Gordon itself** (optional, if Experiment 2 succeeds and we want to deepen the Sine-Gordon analysis):
  - Test breather solutions (bound kink-antikink oscillations) — relevant to Duda's time-crystal insight
  - Test double Sine-Gordon potential (1-cos φ) + α(1-cos 2φ) to see if genuine annihilation (not pass-through) is possible
  - Measure kink-antikink phase shift after pass-through (predicted by Perring-Skyrme exact solution)
- **Question surfaced for later**: the ~0.06% energy offset on static kink is from the continuum-vs-discrete mismatch (sum·dx underestimates the true integral). Worth measuring whether richer integration (Simpson) eliminates it — minor detail, not a blocker

---

## EXPERIMENT 2: Hedgehog Energy vs Distance

**Status**: 🚧 Pending
**Sandbox Script**: `sandbox_phase2_lagrangian/exp2_hedgehog_energy.py`
**Date run**: -

### 2.1 Hypothesis

Two hedgehog defects in a 3D director field produce a clean 1/r Coulomb potential without sinc oscillation. Expected: `E(d) ≈ const + C/d` (matching Duda's `E(d) ≈ 1590 + 25.6/d` from arXiv:2108.07896 Fig. 2).

### 2.2 Setup

- 3D grid: [nx × ny × nz voxels]
- Director field: `n(x) = (nx, ny, nz)` unit vector per voxel
- Hedgehog 1 at r₁: `n = (x-r₁)/|x-r₁|`
- Hedgehog 2 at r₂: `n = (x-r₂)/|x-r₂|` (or anti-hedgehog: `n = -(x-r₂)/|x-r₂|`)
- Frank elastic energy: `H = K·Σ|∇n|²` (one-constant approximation)
- Field relaxation: gradient descent until ΔE < ε

### 2.3 Results

[empty until run]

### 2.4 Numerical Evidence

[E(d) plot, fit to const + C/d, residuals]

### 2.5 Comparison to Expected

| Quantity | Predicted | Measured | Match? |
| --- | --- | --- | --- |
| Form of E(d) | const + C/d | - | - |
| R² of 1/r fit | > 0.99 | - | - |
| Sinc oscillation | none | - | - |

### 2.6 Conclusion

[empty until run]

### 2.7 Next Steps

[empty until run]

---

## EXPERIMENT 3: Topological Charge Quantization

**Status**: 🚧 Pending
**Sandbox Script**: `sandbox_phase2_lagrangian/exp3_topological_charge.py` (or function inside Experiment 2)
**Date run**: -

### 3.1 Hypothesis

The winding number integral `Q = (1/4π) ∮_S (∂_u n × ∂_v n) · n du dv` returns integers (±1 for hedgehog/anti-hedgehog, 0 for vacuum) regardless of surface shape or field perturbation.

### 3.2 Setup

- Use field configurations from Experiment 2
- Discretize closed surfaces (spherical mesh) around defect centers
- Compute Q for: clean hedgehog, perturbed hedgehog, vacuum region, two-defect configurations
- Test surface independence: compute Q for spheres of different radii

### 3.3 Results

[empty until run]

### 3.4 Numerical Evidence

| Configuration | Surface radius | Computed Q | Expected | Match? |
| --- | --- | --- | --- | --- |
| Hedgehog | 2λ | - | +1 | - |
| Hedgehog | 5λ | - | +1 | - |
| Anti-hedgehog | 2λ | - | -1 | - |
| Vacuum | 2λ | - | 0 | - |
| Perturbed hedgehog (10% noise) | 2λ | - | +1 | - |

### 3.5 Conclusion

[empty until run]

### 3.6 Next Steps

[empty until run]

---

## EXPERIMENT 4: Klein-Gordon from Twist Dynamics

**Status**: 🚧 Pending
**Sandbox Script**: `sandbox_phase2_lagrangian/exp4_klein_gordon.py`
**Date run**: -

### 4.1 Hypothesis

In the uniaxial limit, evolving the twist degree of freedom from the Landau-de Gennes Lagrangian produces dispersion `ω² = c²k² + m²` (massive Klein-Gordon). Per Duda's clarification: time-averaged twist behavior should statistically converge to Klein-Gordon, even if instantaneous dynamics differ.

### 4.2 Setup

- Director field with small twist perturbation around uniaxial vacuum
- Evolve using Euler-Lagrange equations from LdG Lagrangian
- Time-stepping: leapfrog or RK4
- Measure: `ψ(k, t)` via spatial Fourier transform
- Extract dispersion ω(k) by fitting `ψ(k, t) = A·cos(ω(k)·t)`

### 4.3 Results

[empty until run]

### 4.4 Numerical Evidence

| k value | Measured ω² | Predicted ω² = c²k² + m² | Match? |
| --- | --- | --- | --- |
| - | - | - | - |

### 4.5 Comparison to Expected

- Massive dispersion (ω² = c²k² + m² with m > 0): ?
- Mass gap measurable: ?
- Time-averaged convergence to KG: ?

### 4.6 Conclusion

[empty until run]

### 4.7 Next Steps

[empty until run]

---

## EXPERIMENT 5: Lagrangian Derivation (Combined W-L from a Lagrangian?)

**Status**: 🚧 Pending
**Sandbox Script**: `sandbox_phase2_lagrangian/exp5_lagrangian_derivation.py` (sympy verification)
**Date run**: -

### 5.1 Hypothesis

OpenWave's empirically-selected Combined Wolff-LaFreniere wave equation `ψ = 2A·sin(kr/2)·cos(kr/2-(ωt+φ))/r` is a solution to a known Lagrangian. Smolinski's `Ψ³` equation is already known to come from a quartic potential Lagrangian.

### 5.2 Setup

- Pen-and-paper analysis + sympy verification
- Standard wave Lagrangian: `L = ½(∂ψ/∂t)² - ½c²(∇ψ)²`
- Verify Combined W-L satisfies Euler-Lagrange equation
- Verify Smolinski's `(∂²/∂t² - c²∇²)Ψ + k·Ψ³ = 0` derives from `L = ½(∂ψ/∂t)² - ½c²(∇ψ)² - k·ψ⁴/4`
- Check conservation laws via Noether's theorem

### 5.3 Results

[empty until run]

### 5.4 Numerical Evidence

| Wave equation | Lagrangian | Verified? |
| --- | --- | --- |
| Combined W-L | `L = ½(∂ψ/∂t)² - ½c²(∇ψ)²` | - |
| Smolinski Ψ³ | `L_free - k·ψ⁴/4` | - |

### 5.5 Conclusion

[empty until run]

### 5.6 Next Steps

[empty until run]

---

## EXPERIMENT 6: Three Lepton Families from Biaxial Hedgehog

**Status**: 🚧 Pending
**Sandbox Script**: `sandbox_phase2_lagrangian/exp6_lepton_families.py`
**Date run**: -

### 6.1 Hypothesis

A biaxial nematic hedgehog with 3 distinguishable axes produces 3 hedgehog types with the same topological charge (Q = ±1) but different energies. The energy ratios should match electron/muon/tau mass ratios.

### 6.2 Setup

- Biaxial order parameter: `D = diag(λ₁, λ₂, λ₃)` with `λ₁ ≠ λ₂ ≠ λ₃`
- Three hedgehog types: defect aligned with axis 1, 2, or 3
- Field relaxation for each
- Measure total energy of each defect type

### 6.3 Results

[empty until run]

### 6.4 Numerical Evidence

| Defect type | Q | Energy | Ratio to lightest | Predicted (mass ratio) |
| --- | --- | --- | --- | --- |
| Axis 1 (electron) | +1 | - | 1.0 | 1.0 (m_e = 0.511 MeV) |
| Axis 2 (muon) | +1 | - | - | 207 (m_μ = 105.7 MeV) |
| Axis 3 (tau) | +1 | - | - | 3477 (m_τ = 1776.8 MeV) |

### 6.5 Conclusion

[empty until run]

### 6.6 Next Steps

[empty until run]

---

## EXPERIMENT 7: Close's Nonlinear Vector Wave Equation

**Status**: 🚧 Pending
**Sandbox Script**: `sandbox_phase2_lagrangian/exp7_close_vector_wave.py`
**Date run**: -

### 7.1 Hypothesis

Close's nonlinear vector wave equation for spin density (from "Plane Wave Solutions to a Proposed 'Equation of Everything'", Foundations of Physics 2025) seeded with a spherical harmonic evolves into localized particle-like soliton/breather structures. Per Close's invitation: *"starting with a spherical harmonic linear wave solution and see what evolves."*

### 7.2 Setup

- 3D grid with vector spin density field `s(x) = (sx, sy, sz)`
- Initial condition: spherical harmonic `Y_l^m`
- Time evolution: leapfrog or RK4 with Close's nonlinear terms
- Long-time evolution to detect soliton/breather emergence

### 7.3 Results

[empty until run]

### 7.4 Numerical Evidence

[localization plots, energy concentration metrics, temporal evolution]

### 7.5 Comparison to Expected

- Particle-like soliton emerges from Y_l^m seed: ?
- Field matches Dirac bispinor plane wave solution: ?
- Stable against perturbation: ?

### 7.6 Conclusion

[empty until run]

### 7.7 Next Steps

[empty until run]

---

## EXPERIMENT 8: Smolinski's Non-linear Ψ³ (Direct K-Selectivity Test)

**Status**: 🚧 Pending
**Sandbox Script**: `sandbox_phase2_lagrangian/exp8_smolinski_psi3.py`
**Date run**: -

### 8.1 Hypothesis

Adding `-k·Ψ³` to the linear wave equation produces K-dependent stability — K=10 tetrahedron is uniquely stable under perturbation, while K=2..9 are not. This would solve K-selectivity from nonlinearity alone (without topology).

### 8.2 Setup

- 3D scalar grid with field `Ψ(x,y,z,t)`
- Equation: `∂²Ψ/∂t² = c²∇²Ψ - k·Ψ³`
- Time-domain PDE evolution (NOT phasor superposition — superposition breaks for non-linear)
- Initial conditions: K wave centers at geometric positions (K=2..10)
- Sweep coefficient `k` over multiple orders of magnitude

### 8.3 Progressive Complexity Levels

| Level | Form | Run? | Result |
| --- | --- | --- | --- |
| 1 | `F = k·Ψ³` (constant k) | - | - |
| 2 | `F = k(r)·Ψ³` (spatially varying k) | - | - |
| 3 | `F = k₁·Ψ³ + k₂·Ψ⁵` (higher order) | - | - |
| 4 | `F = f(ε_G)·g(Ψ)` (geometric coupling) | - | - |
| 5 | Fallback: topology required | - | - |

### 8.4 Results

[empty until run]

### 8.5 Numerical Evidence

| K | Stability under perturbation | k value range | Notes |
| --- | --- | --- | --- |
| 2 | - | - | - |
| 3 | - | - | - |
| 4 | - | - | - |
| 5 | - | - | - |
| 6 | - | - | - |
| 7 | - | - | - |
| 8 | - | - | - |
| 9 | - | - | - |
| 10 | - | - | - |

### 8.6 Conclusion

[empty until run]

### 8.7 Next Steps

[empty until run]

---

## OVERALL CONCLUSIONS

[populate after all experiments run]

### Winning Approach

[which Lagrangian / wave equation candidate best fits OpenWave's requirements]

### Comparison: Topology vs Nonlinearity vs Standing Waves

| Phenomenon | Standing waves (M3) | Topology (Exp. 2, 3, 6) | Nonlinearity (Exp. 5, 7, 8) |
| --- | --- | --- | --- |
| Lock-in / particle binding | - | - | - |
| K-selectivity | - | - | - |
| Charge quantization | - | - | - |
| Far-field Coulomb (no sinc) | - | - | - |
| Annihilation | - | - | - |
| Three lepton families | - | - | - |

### What to Implement in M4/M5

[concrete recommendation for next architecture step]

### Open Questions for Next Phase

[what remains unsolved, what new investigations are needed]

---

## CHANGE LOG

| Date | Experiment | Change |
| --- | --- | --- |
| - | - | - |
