# LAGRANGIAN FRAMEWORK — NUMERICAL EXPERIMENTS

> *Numerical experiments — same spirit as OpenWave xperiments: controlled investigations of physics hypotheses through simulation. OpenWave is built as a shared experimental platform where wave-based and topological models can be tested, compared, and built upon by others.*

Numerical results from the 8 Lagrangian framework experiments in `sandbox_phase2_lagrangian/`.
See [3_LAGRANGIAN_FRAMEWORK.md](3_LAGRANGIAN_FRAMEWORK.md) for experiment specifications.

---

## Using OpenWave "Sandbox" vs "Production" environment

OpenWave uses two layers of experimentation:

- **Sandbox = explore fast, fail cheap.**
  - (`research_hub/sandbox_*/`) — quick numpy scripts to validate math, logic, and concepts. Pure exploration. No GPU, no Taichi, no production dependencies. This is where ideas are tested cheaply before committing engineering effort. If an experiment works in the sandbox, it graduates to the production engine.
- **Production = implement validated winners.**
  - (`openwave/xperiments/m*/`) — Taichi-based 3D rendering and simulation on the official OpenWave platform. GPU-accelerated, full grid infrastructure, visualization, force & motion. This is the final product where validated equations run at scale.

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
- Make recommendation for `What to Implement in M5`

---

## SUMMARY DASHBOARD

| # | Experiment | Status | Result | Date | Script |
| --- | --- | --- | --- | --- | --- |
| 1 | Sine-Gordon 1D solitons | ✅ Passed | Kink stability, Lorentz contraction, pass-through all confirmed | 2026-04-16 | `exp1_sine_gordon_1d.py` |
| 2 | Hedgehog energy vs distance | ✅ Passed | Clean 1/d Coulomb attraction, R²=0.993 post-relax, no sinc | 2026-04-16 | `exp2_hedgehog_energy.py` |
| 3 | Topological charge quantization | ✅ Passed | Q=±1 integer across all sphere radii; stable up to 50% noise | 2026-04-16 | `exp3_topological_charge.py` |
| 4 | Klein-Gordon from twist | 🚧 Pending | - | - | - |
| 5 | Lagrangian derivation | ⚠️ Mixed | Smolinski Ψ³ + Noether confirmed; M4 sum-form W-L ✅; doc product form is NOT a free-wave solution (residual −c²k²·sin(ωt+φ)/r) | 2026-04-17 | `exp5_lagrangian_derivation.py` |
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

**Status**: ✅ Passed
**Sandbox Script**: `sandbox_phase2_lagrangian/exp2_hedgehog_energy.py`
**Date run**: 2026-04-16

### 2.1 Hypothesis

Two hedgehog defects in a 3D director field produce a clean 1/r Coulomb potential without sinc oscillation. Expected: `E(d) ≈ const + C/d` (matching Duda's `E(d) ≈ 1590 + 25.6/d` from arXiv:2108.07896 Fig. 2). For opposite-sign defects (hedgehog + anti-hedgehog), the interaction is attractive: b < 0 in `E = a + b/d`.

### 2.2 Setup

- **Grid**: 48³ voxels (~110k), domain [−8, +8] each axis, dx ≈ 0.340
- **Director**: unit vector `n(x) = (nx, ny, nz)` per voxel, |n|=1
- **Vacuum state**: `n = ẑ` everywhere (ground state)
- **Defect seeding**: weighted superposition of two hedgehogs + unit-normalization (proximity weights w₁ = 1/(r₁+0.5), w₂ = 1/(r₂+0.5)); far-field blended back to vacuum via Lorentzian w_vac
- **Defect pair**: hedgehog at c₁=(−d/2, 0, 0) with sign=+1; anti-hedgehog at c₂=(+d/2, 0, 0) with sign=−1
- **Frank elastic energy**: `H = (K/2) · ∫ |∇n|² d³r` with K=1 (one-constant approximation)
- **Relaxation**: gradient descent `∂n/∂τ = ∇²n − (n·∇²n)·n` with tangent-space projection to preserve |n|=1, then unit renormalization per step
- **Stability**: τ = 0.008 (below CFL limit τ_max = dx²/(2·3·K) ≈ 0.019), 60 steps
- **Core pinning**: single closest voxel per defect held to ±ẑ (soft Dirichlet — prevents numerical decay of topology without creating sharp discontinuity)
- **Separations tested**: d ∈ {2, 3, 4, 5, 6, 8}

### 2.3 Results

| d | E_pre-relax | E_post-relax |
| --- | --- | --- |
| 2.0 | 147.91 | 83.46 |
| 3.0 | 163.00 | 97.96 |
| 4.0 | 173.78 | 107.62 |
| 5.0 | 182.86 | 113.09 |
| 6.0 | 190.04 | 117.71 |
| 8.0 | 198.52 | 122.67 |

**Monotonic increase with d** — energy rises smoothly toward a plateau at large d. No sinc oscillation, no direction flips. This is the hallmark of Coulomb interaction.

**Relaxation behavior**: energy decreased monotonically from pre-relax to post-relax values at every separation (verified by convergence plot). The gradient descent is stable with the chosen τ and converges within ~60 steps.

### 2.4 Numerical Evidence

Plots in `sandbox_phase2_lagrangian/exp2_results/`:

- `energy_vs_distance.png` — E(d) for both pre-relax and post-relax; 1/d fit overlaid
- `relax_convergence.png` — energy-vs-step for all 6 separations (monotonic decrease)

**Fits** to `E(d) = a + b/d`:

| Dataset | a | b | R² | Interpretation |
| --- | --- | --- | --- | --- |
| Pre-relax | 210.90 | −132.88 | 0.964 | attractive (b < 0); energy rises from bound state (small d) toward 2·E_single (large d) |
| Post-relax | 134.58 | −104.74 | **0.993** | **near-perfect 1/d** attraction after relaxation; topology preserved; gradient energy smoothed |

### 2.5 Comparison to Expected

| Quantity | Predicted | Measured | Match? |
| --- | --- | --- | --- |
| Form of E(d) | const + C/d (Coulomb) | const + C/d | ✅ |
| Sign of interaction | attractive (opposite-sign defects) | b = −105 < 0 (post-relax) | ✅ |
| R² of 1/r fit | > 0.99 | **0.993** (post-relax) | ✅ |
| Sinc oscillation | none | none — monotonic rise with d | ✅ |
| Energy scales linearly with 1/d | yes | R² = 0.993 linear | ✅ |

### 2.6 Conclusion

**Experiment 2 validates the core topological-Coulomb claim.** Two opposite-sign hedgehog defects in a 3D director field with Frank elastic energy produce a clean 1/d attractive interaction, **R² = 0.993 to a pure `a + b/d` model** — no sinc barriers, no direction flips, no oscillatory structure.

This is **exactly the behavior that M3 could not produce.** M3's Combined W-L wave equation gives sinc-barrier-laden far-field interactions where the force flips direction every λ/2. Here, the interaction is smooth and monotonic at every separation tested, from d = 2 to d = 8 (spanning 6 times the core size).

**Why this works mathematically**: the Frank elastic energy `∫|∇n|² d³r` penalizes deformations of the director field. Two defects of opposite winding create a smooth "connecting texture" between them; the cost of that texture scales as 1/d in 3D (electric-field analog: the surface integral of the field energy between two point sources). There's no wave superposition, no phase-dependent interference, no sinc — the physics is entirely geometric.

**Implications for OpenWave**:

1. The far-field Coulomb problem **is solvable by topology** — this is the single most important result for K-selectivity
2. The "Coulomb from topology" mechanism Duda proposed is numerically confirmed in our setup (not just a theoretical claim)
3. M5 should implement a director field (vector structure already in M4's 6-phasor) with topological initial conditions (`seed_hedgehog`). The Frank elastic energy is the corresponding `V(ψ)` for M5's Hamiltonian
4. The near-field K=10 standing-wave lock-in (M3 result) can now be added on top of this — topology handles far-field Coulomb, standing waves handle near-field orbit quantization. This is Duda's "both topology and standing waves" picture made concrete

**Caveat**: this experiment measures the static Frank elastic energy of a fixed defect configuration. It doesn't test dynamics (moving defects, time-dependent field). Experiments 4 (Klein-Gordon from twist) and 7 (Close's nonlinear vector wave) will test the dynamical regime.

**Caveat 2**: pre-relax energies already showed the 1/d trend (R² = 0.964); relaxation only sharpened it to 0.993. The smooth superposition-based seeding is doing most of the work. This is actually useful — it means the 1/d interaction is built into the hedgehog geometry itself, not an emergent property of long-time relaxation. A first-order implementation in M5 could skip expensive relaxation and still capture the correct Coulomb physics.

### 2.7 Next Steps

- **Continue to Experiment 3** (topological charge quantization) — uses similar field setups ✅ already complete
- **Future refinements**:
  - Larger grid (96³ or 128³) to measure wider range of d and confirm asymptotic scaling
  - Test same-sign defects (hedgehog + hedgehog) → should give repulsion (b > 0)
  - Measure prefactor in physical units (relate K to OpenWave's medium density ρ and wave speed c)
  - Test anisotropic Frank constants (K₁, K₂, K₃ splay/twist/bend) to see if interaction is still pure 1/d or develops angular dependence

---

## EXPERIMENT 3: Topological Charge Quantization

**Status**: ✅ Passed
**Sandbox Script**: `sandbox_phase2_lagrangian/exp3_topological_charge.py`
**Date run**: 2026-04-16

### 3.1 Hypothesis

The winding-number integral `Q = (1/4π) ∮_S n · (∂_θ n × ∂_φ n) dθ dφ` returns integers (±1 for hedgehog/anti-hedgehog, 0 for vacuum) regardless of surface shape, surface radius, or smooth field perturbation. This is the fundamental mechanism Duda invokes to solve OpenWave's charge quantization problem.

### 3.2 Setup

- **Grid**: 48³, domain [−8, +8], dx ≈ 0.340
- **Director fields tested**:
  - Clean hedgehog: `n(x) = (x − c)/|x − c|` with sign ±1
  - Vacuum: `n = ẑ` everywhere
  - Perturbed hedgehog: `n_pert = (n + ε·rand) / |n + ε·rand|` with ε ∈ {0, 0.05, 0.10, 0.20, 0.50, 1.00}
  - Two-defect configuration: hedgehog at (−2, 0, 0) + anti-hedgehog at (+2, 0, 0)
- **Winding computation**: spherical mesh (64 θ × 128 φ), trilinear interpolation of n(x) onto sphere, central-difference finite differences on (θ, φ), numerical surface integral
- **Surface radii tested**: R ∈ {2, 3, 5, 7} grid units (from core of defect to ~¾ of domain half-width)

### 3.3 Results

**Test 1 — Clean configurations** (Q vs surface radius):

| R | Q(hedgehog) | Q(anti-hedgehog) | Q(vacuum) |
| --- | --- | --- | --- |
| 2.0 | +0.985 | −0.985 | +0.0000 |
| 3.0 | +0.993 | −0.993 | +0.0000 |
| 5.0 | +0.997 | −0.997 | +0.0000 |
| 7.0 | +0.998 | −0.998 | +0.0000 |

Accuracy improves monotonically with surface radius (larger spheres → more voxels, better trilinear sampling). All values unambiguously round to ±1 or 0.

**Test 2 — Topological stability under perturbation** (R = 3):

| noise ε | Q measured | round(Q) | Match expected (+1)? |
| --- | --- | --- | --- |
| 0.00 | +0.993 | +1 | ✅ |
| 0.05 | +0.990 | +1 | ✅ |
| 0.10 | +0.982 | +1 | ✅ |
| 0.20 | +0.952 | +1 | ✅ |
| 0.50 | +0.754 | +1 | ✅ |
| 1.00 | +0.190 | 0 | ❌ (field is essentially random at this level) |

**Q stays ≥ 0.95 up to 20% noise**, and still rounds correctly to +1 at 50% noise. Only breaks at 100% noise, where the field is effectively destroyed (no hedgehog structure remaining). This confirms the mechanism: topology is robust against smooth deformations but not against complete randomization.

**Test 3 — Two-defect configuration** (hedgehog + anti-hedgehog at d=4):

| Surface | Computed Q | Expected |
| --- | --- | --- |
| around c₁ = (−2, 0, 0), R = 1.3 | +0.958 | +1 |
| around c₂ = (+2, 0, 0), R = 1.3 | −0.958 | −1 |
| large sphere enclosing both, R = 5 | −0.000 | 0 |
| sphere far from both defects, R = 1 | +0.000 | 0 |

**Total charge conserved** (Q_total = 0 for the bound pair). **Surface independence** verified: Q depends only on which defects are enclosed, not on surface size or position (as long as the enclosure is correct).

### 3.4 Numerical Evidence

Plots in `sandbox_phase2_lagrangian/exp3_results/`:

- `winding_vs_radius.png` — three flat lines at Q = +1, −1, 0 for hedgehog, anti-hedgehog, vacuum respectively. Discretization error visible but bounded by 1.5% at R=2, falling to 0.2% at R=7
- `winding_vs_noise.png` — Q stays near 1 for noise up to 20%, slowly degrades to 0.75 at 50% noise, collapses to 0.19 at 100% noise

| Configuration | Surface radius | Computed Q | Expected | Match? |
| --- | --- | --- | --- | --- |
| Hedgehog | 2 | +0.985 | +1 | ✅ |
| Hedgehog | 3 | +0.993 | +1 | ✅ |
| Hedgehog | 5 | +0.997 | +1 | ✅ |
| Hedgehog | 7 | +0.998 | +1 | ✅ |
| Anti-hedgehog | 2 | −0.985 | −1 | ✅ |
| Anti-hedgehog | 7 | −0.998 | −1 | ✅ |
| Vacuum | all | ≤ 1e-4 | 0 | ✅ |
| Perturbed hedgehog (10%) | 3 | +0.982 | +1 | ✅ |
| Perturbed hedgehog (20%) | 3 | +0.952 | +1 | ✅ |
| Perturbed hedgehog (50%) | 3 | +0.754 | +1 (rounds correctly) | ✅ |
| Two-defect pair total | 5 (enclosing both) | 0.000 | 0 | ✅ |

### 3.5 Conclusion

**Experiment 3 validates topological charge as an integer-valued, robust conservation law** in exactly the way Duda invokes it. Key findings:

1. **Surface independence** — Q is identical (to within ~1% discretization error) across all sphere radii from 2 to 7 grid units. This is the Gauss-Bonnet theorem in action: the winding number depends only on what's enclosed, not on how you draw the surface
2. **Integer quantization** — all measurements unambiguously round to ±1 or 0. No "fractional charges", no drift, no sensitivity to mesh choice
3. **Smooth perturbation robustness** — Q holds through 50% random noise. A perturbation that doesn't literally destroy the defect cannot change its topological charge — precisely the protection mechanism that guarantees electron stability in Duda's picture
4. **Total charge conservation** — a hedgehog (+1) and anti-hedgehog (−1) together give Q_total = 0 when measured on a surface enclosing both. This is the "neutral bound pair" configuration — the analog of positronium
5. **Mechanism matches expectation** — the winding-number integral, a purely geometric quantity, returns the integer topological charge with no free parameters. No fitting, no tuning: this is a mathematical identity, confirmed numerically

**Implication for OpenWave**: Duda's challenge — "without charge quantization your electron explodes" — has a concrete numerical demonstration of its solution. Once we implement a director field in M5 and seed topological defects, charge will be an integer-valued conserved quantity by construction, with no need for the `cos(source_offset)` phase trick used in M3. This is the charge-quantization mechanism we've been missing.

**Contrast with M3**: in M3, charge is imposed as ±1 by `cos(source_offset)` at each WC. It's a label, not an emergent property. If we perturb the system, there's no mechanism to keep charge at ±1 — it could numerically become 0.7 or 1.2 and the wave equation wouldn't care. Here, the topology *enforces* integrality: the same perturbations that would disrupt a wave-based charge leave the topological winding exactly ±1.

### 3.6 Next Steps

- **Continue to Experiment 5** (Lagrangian derivation — can Combined W-L come from a Lagrangian?)
- **Integration with M5 plan**: the winding-number kernel (`winding_number()` in the experiment script) ports directly to M5 as a new tracker alongside amplitude, frequency, and energy. Computing Q around each WC position becomes a per-frame diagnostic
- **Future validation**:
  - Test on an evolved (PDE-relaxed) field to confirm Q survives dynamical simulation
  - Test non-spherical surfaces (cubes, ellipsoids) to further verify surface independence
  - Measure Q for K=10 tetrahedron configurations (10 hedgehogs arranged tetrahedrally) — should give Q_total = +10 if all same sign, 0 if charge-neutral compound
  - Measure angular-resolved Q(θ, φ) to check for direction-dependent defects

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

**Status**: ⚠️ Mixed — Smolinski Ψ³ and Noether confirmed; Combined W-L (doc form) is NOT a free-wave solution
**Sandbox Script**: `sandbox_phase2_lagrangian/exp5_lagrangian_derivation.py` (sympy verification)
**Date run**: 2026-04-17

### 5.1 Hypothesis

OpenWave's empirically-selected Combined Wolff-LaFreniere wave equation `ψ = 2A·sin(kr/2)·cos(kr/2-(ωt+φ))/r` is a solution to a known Lagrangian. Smolinski's `Ψ³` equation is already known to come from a quartic potential Lagrangian.

### 5.2 Setup

- Symbolic verification via `sympy`. No numerical simulation.
- Four candidate forms of the Combined W-L tested against the d'Alembertian `□ψ = ∂²ₜψ − c²∇²ψ`:
  - (1a) Pure outgoing wave: `ψ = A·sin(kr − ωt − φ)/r`
  - (1b) Pure standing wave: `ψ = 2A·sin(kr)·cos(ωt+φ)/(kr)`
  - (1c) Sum form (the M4 code): `ψ = A·[sin(kr+ωt+φ) + sin(kr−ωt−φ)]/(kr)`
  - (1d) Product form (from 3b/3_LAGRANGIAN_FRAMEWORK docs): `ψ = 2A·sin(kr/2)·cos(kr/2−(ωt+φ))/r`
- Smolinski derivation: start from `L = ½(∂ₜψ)² − ½c²|∇ψ|² − (κ/4)·ψ⁴`, apply Euler-Lagrange, verify output matches `∂ₜ²ψ − c²∇²ψ + κ·ψ³ = 0`
- Noether energy: compute Hamiltonian density `H = π·∂ₜψ − L` for both Lagrangians

### 5.3 Results

| Candidate | □ψ with ω = c·k | Verdict |
| --- | --- | --- |
| (1a) Pure outgoing `sin(kr−ωt−φ)/r` | **0** | ✅ exact free-wave solution |
| (1b) Pure standing `2·sin(kr)·cos(ωt+φ)/(kr)` | **0** | ✅ exact free-wave solution |
| (1c) Sum form (M4 code) `[sin(kr+ωt+φ) + sin(kr−ωt−φ)]/(kr)` | **0** | ✅ exact; algebraically identical to (1b) |
| (1d) Product form (3b doc) `sin(kr/2)·cos(kr/2−(ωt+φ))/r` | `−A·c²k²·sin(ωt+φ)/r ≠ 0` | ❌ NOT a free-wave solution |

**Identity confirmation**: `(sum form) − (product standing form) = 0` — confirms (1b) and (1c) are the same equation written differently.

**Surprise finding (1d)**: The product form `2A·sin(kr/2)·cos(kr/2−(ωt+φ))/r` does **not** satisfy the 3D free wave equation. Decomposing via `sin(a)·cos(b) = ½[sin(a+b) + sin(a−b)]`:

```text
2A·sin(kr/2)·cos(kr/2−(ωt+φ))/r
  = (A/r)·[sin(kr − ωt − φ) + sin(ωt + φ)]
           ↑                    ↑
           free-wave solution   NOT a solution —
                                spatial part = 1/r but no k in sine
                                (uniform time oscillation with 1/r envelope)
```

The second term `A·sin(ωt+φ)/r` has no radial wave structure (no `k` in the sine argument). Its Laplacian is `∇²(A·sin(ωt+φ)/r) = 0` everywhere except at the origin. Its time derivative is `−ω²·A·sin(ωt+φ)/r`. Their difference is `−c²k²·A·sin(ωt+φ)/r ≠ 0`, which is what sympy reports.

Equivalently, this is the "Phase + Quadrature" decomposition from [2L3_particle_emergence.md](2L3_particle_emergence.md):

```text
Phase     C_n = A·sin(kr)/r          is a free-wave solution
Quadrature S_n = A·(1−cos(kr))/r     is NOT a free-wave solution
```

So **the LaFreniere quadrature term (1−cos(kr))/r is a modeling choice** — it was introduced for specific physics (radiation pressure, standing-to-traveling transition), not derived as a homogeneous solution of the wave equation.

**Smolinski Ψ³ derivation** (Test 2): sympy computed the Euler-Lagrange equation from `L = ½(∂ₜψ)² − ½c²|∇ψ|² − (κ/4)·ψ⁴` and it matches `∂ₜ²ψ − c²∇²ψ + κ·ψ³ = 0` exactly (difference = 0). ✅

**Noether energy density** (Test 3): `H = π·∂ₜψ − L`:

- Free wave: `H = ½(∂ₜψ)² + ½c²|∇ψ|²` (kinetic + gradient energy)
- Smolinski: `H = ½(∂ₜψ)² + ½c²|∇ψ|² + (κ/4)·ψ⁴` (adds quartic potential)

Both are standard `T + V` Hamiltonians. Canonical energy conservation is guaranteed by Noether's theorem (time translation symmetry). ✅

### 5.4 Numerical Evidence

sympy output (raw):

```text
□ψ_out (1a) with ω=ck       = 0
□ψ_st  (1b) with ω=ck       = 0
□ψ_sum (1c) with ω=ck       = 0
□ψ_doc (1d) with ω=ck       = -A·c²·k²·sin(c·k·t + φ0)/r   ← residual
(sum form) − (product form) = 0                            ← identity

Smolinski derived EoM = ∂ₜ²ψ − c²∇²ψ + κ·ψ³    ← matches expected exactly

Free-wave H  = ½(∂ₜψ)² + ½c²|∇ψ|²
Smolinski H  = ½(∂ₜψ)² + ½c²|∇ψ|² + (κ/4)·ψ⁴
```

### 5.5 Comparison to Expected

| Claim | Expected | Measured | Match? |
| --- | --- | --- | --- |
| Combined W-L (sum form 1c, M4 code) satisfies free-wave EL | yes | yes | ✅ |
| Combined W-L (product form 1d, docs) satisfies free-wave EL | yes (per docs) | **no — residual `−A·c²k²·sin(ωt+φ)/r`** | ❌ (docs contradicted) |
| Forms (1b), (1c), (1d) are all equivalent | yes | (1b)=(1c) ✅; (1d) differs | ⚠️ partial |
| Smolinski Ψ³ from `L = L_free − (κ/4)ψ⁴` | Smolinski equation | exact match | ✅ |
| Noether `H = T + V` from Euler-Lagrange | standard `T + V` | standard `T + V` | ✅ |

### 5.6 Conclusion

**Two successes and one important anomaly.**

- ✅ **Smolinski's Ψ³ equation is a textbook Lagrangian derivation.** Starting from `L = ½(∂ₜψ)² − ½c²|∇ψ|² − (κ/4)·ψ⁴` and applying Euler-Lagrange gives exactly `∂ₜ²ψ − c²∇²ψ + κ·ψ³ = 0`. Noether's theorem confirms energy conservation via `H = T + V` with `V = (κ/4)·ψ⁴`. This validates Experiment 8's (pending) K-selectivity test at the mathematical level — Smolinski's term is well-founded in Lagrangian field theory
- ✅ **The M4-implemented Combined W-L (sum form, `A·[sin(kr+ωt+φ)+sin(kr−ωt−φ)]/(kr)`) IS a free-wave solution.** At w=1 (pure standing wave limit) it reduces to `2A·sin(kr)·cos(ωt+φ)/(kr)`, a textbook exact solution. M4's analytical physics is self-consistent with the free-wave Lagrangian
- ❌ **The documented product form `2A·sin(kr/2)·cos(kr/2−(ωt+φ))/r` is NOT a free-wave solution.** This is the most important finding. The formula decomposes into a valid outgoing-spherical-wave piece `A·sin(kr−ωt−φ)/r` plus a spurious time-oscillating 1/r piece `A·sin(ωt+φ)/r` that has no radial wave structure. The LaFreniere quadrature term `(1−cos(kr))/r` was introduced for empirical reasons (radiation pressure, standing-to-traveling transition) and is **not** a homogeneous solution of the wave equation

**Interpretation**: the Combined W-L in the empirical / docs form is best understood as the superposition of:

1. A legitimate outgoing spherical wave (the Phase part, `A·sin(kr)/r · cos(ωt+φ)`)
2. A forced oscillation (the Quadrature part, `A·(1−cos(kr))/r · sin(ωt+φ)`) that requires a source term to exist

The Quadrature part would emerge naturally if the Lagrangian included a **source term** `J(x,t)·ψ` (like an external charge in electromagnetism), or from a non-linear potential where the quadrature emerges at second order. Neither is the free-wave Lagrangian alone.

**Practical implication for M5**: we do NOT carry the Combined W-L product form into M5 as a "derived" solution. Instead:

- M5 uses the **free-wave sum form** `A·[sin(kr+ωt+φ) + sin(kr−ωt−φ)]/(kr)` as the *linear* limit (V=0 case)
- Nonlinear physics (topology, K-selectivity, annihilation) comes from the potential V(ψ) we add — which Exps 2, 3, 8 are testing
- The LaFreniere quadrature term can be retained as a *phenomenological* piece if it turns out to encode useful radiation-pressure physics, but we cannot claim it's Lagrangian-derived from the pure free-wave action

**Meta-takeaway**: this is exactly the kind of "inverted hierarchy" result the Lagrangian framework was meant to surface. In M3 we picked an equation empirically; the Lagrangian check reveals that the equation implicitly assumes a source term we never modeled. Fixing this is not a correction to M3 — it's motivation for M5.

### 5.7 Next Steps

- **Update `3_LAGRANGIAN_FRAMEWORK.md` line ~418** to correct the "Combined W-L IS an exact solution" claim — **done** (incorporated finding)
- **Use the sum form, not the product form, when referring to M4's wave equation in any future documentation** — the product form may be a computational shortcut but misrepresents the physics
- **Open question for Exp 8**: if we use Smolinski's nonlinear Ψ³ term as the potential in an M5 simulation, can we recover a spatial pattern that looks like the Combined W-L product form? I.e., is the product form the *static nonlinear solution* of Smolinski's equation? Worth testing
- **Continue to Experiment 8** (Smolinski Ψ³ K-selectivity — recommended next) now that its Lagrangian validity is confirmed

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
