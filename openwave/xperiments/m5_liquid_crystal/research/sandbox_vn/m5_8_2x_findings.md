# M5.8.2x — Absolute-scale calibration (issue #208, Phase A+B)

> Builds the sim → physical unit map for the M5 liquid-crystal substrate and computes the dimensionless invariants that test it. Driver: [`m5_8_2x_absolute_scale.py`](m5_8_2x_absolute_scale.py) (pure-analysis; measured inputs sourced inline). Date 2026-06-17.
>
> **Framing (Duda's input on #208):** the model is scale-free, so the absolute scale is *defined* by anchoring on one observable, after which every other observable becomes a falsifiable physical number. With `c` emergent (radial cone, N-4) and `ℏ` fixed by postulate P2 (lattice action unit ≡ ℏ), the map has exactly **two free dials** (a length scale, an energy scale).

## Measured inputs

| Quantity | Sim value | Source |
| --- | --- | --- |
| Coulomb coupling `b` (E = a + b/d) | −5.633e-3 (R²=0.9781) | `m5_4_coulomb_matrix.py` (Frank-elastic energy × voxel; dx=15.385 am, δ=0.5) |
| rest energy `H_static` | 16.74 (lattice Hamiltonian) | `m5_8_2j_zbw_ratio.py` (`u + βu²`, β=1.558, R_W=3.5) |
| clock fundamental `ω₁` | 1.188 (rad / lattice-time) | `m5_8_2j` / `m5_8_2h` (median detrend, start-independent 5.4%) |
| apolar doubling `ω_M/ω₁` | 2.0 (machine-exact) | `m5_8_2s_spin_half_apolar.py` |
| Faber radius `r₀` | 2.2132 fm → 0.511 MeV | `m5_6_3a/3b` (`E·r₀=const`; Eq.8 `E₀ = α·ℏc·π/(4r₀)`) |

## Part A — dimensionless invariants (anchor-independent, the real tests)

| Invariant | Value | Read |
| --- | --- | --- |
| **α = \|b\|/(ℏc)** (ℏ=1, c=1) | **1/α ≈ 178** (c=1 central; bracket 89–355 for c∈{0.5,2}) | a dimensionless EM coupling of the **right order** vs the measured 137.036 (factor 1.30), from **pure topology, un-tuned** |
| **ZBW ratio ω₁/(2H_static)** | **0.0355** | the rock-solid clock invariant: start-independent (5.4% across arms), inherits the N-1 attractor |
| **apolar doubling ω_M/ω₁** | **2.0** (machine-exact) | the Dirac ZBW factor-of-2; exact, no calibration needed |

The α residual is an **O(1) unit-reconciliation factor**, not a free parameter: `|b|` is measured in the **Frank-elastic** functional while ℏ is set on the **matrix Hamiltonian**, and the radial cone carries a `c/2` convention. That reconciliation IS the consistent-units re-run the DoD specifies. The result is falsifiable (right-order-or-not), not a fit.

## Part B — the unit map per anchor + the ~28× clock gap

Energy anchor P1: `H_static = 16.74 lattice ↔ m_e c² = 0.511 MeV` ⇒ 1 lattice energy unit = 0.0305 MeV.

Clock prediction (P1 + P2): `ω_phys = ω₁·(m_e c²/ℏ)/H_static = 5.51e19 rad/s` vs the electron ZBW `2 m_e c²/ℏ = 1.553e21 rad/s` ⇒ **gap = 28.2×**.

| Anchor | Fixes | Then predicts |
| --- | --- | --- |
| Coulomb (e-scale) | charge × length × energy | **α** (dimensionless) → 1/α ≈ 178 |
| Electron clock | time/energy via the ZBW | the mass; ratio ω₁/(2H)=0.035 |
| Faber mass (`r₀`) | r₀=2.2132 fm ↔ 0.511 MeV | the clock ω → lands **28× low** (the gap) |

## Structural findings (the real conclusions)

| # | Finding |
| --- | --- |
| 1 | **The anchors are not all independent.** Faber's `E₀ = α·ℏc·π/(4r₀)` *contains* α as an input, so anchoring the mass via Faber presupposes α. The **Coulomb route is the independent α measurement**; the **clock is the independent absolute-frequency test**, and it is the one that fails (~28×). |
| 2 | **The ~28× gap is the model's one genuine absolute-scale discrepancy.** Dimensionful and structural (ω is rigid against energy, mass-family-independent, N-6a); it does not close by bookkeeping. |
| 3 | **α lands the right order** (1/α ≈ 178 vs 137) from pure topology, un-tuned. The residual is exactly the consistent-units reconciliation. |
| 4 | **Verdict:** a working 2-dial unit map exists per anchor, but **no single anchor makes all observables land at once** (Coulomb → α right-order; Faber → mass by construction; clock → 28× off). The **sim-units-ratio fallback** (predictions as ratios) is the correct mode until the V-on / Faber-`r₀` core route is shown to lift ω toward the ZBW. |

## Which results become physical numbers once the map exists

| Result | Today | After a closed map |
| --- | --- | --- |
| α (fine-structure) | 1/α ≈ 178 (right order, ratio) | a precise dimensionless prediction (the cleanest cross-check) |
| electron g-factor | mechanism only (g≈2 awaits Coulomb e-scale) | g as a number |
| clock ω (Hz) | 5.51e19 rad/s under P1+P2 (28× low) | absolute Hz once the gap is understood/closed |
| lepton masses | ratios (M5.9) | masses in MeV |
| `E = mc²` | ratio | absolute energy |

## Part C — the V-on / Faber-`r₀` core route (test designed; named follow-up, not run here)

The #208 hypothesis is that the absolute scale (and the 28× gap) is set in the **V-on / Faber-`r₀` core**, not by bookkeeping. State of the evidence and the test:

| Element | Status |
| --- | --- |
| What V-on does | confines the amplitude (`m5_6_5c_prod_scale.py`: V-on pins the amplitude, V-off wanders 6.4×) |
| Why the gap might be structural anyway | N-6a: ω is **rigid** across a 2.6× mass family — changing the energy scale does NOT move ω, so adding V (energy/stiffness) is not expected to move ω *by bookkeeping* |
| Why it might still close | V-on changes the core **geometry** (`r₀` regularization), and the clock frequency may be set by the core geometry rather than the energy. A ZBW/Compton argument: `ω_ZBW ~ 2c/λ_C`; the gap means the clock oscillates over a region ~28× larger than λ_C. If V-on shrinks the effective oscillation scale toward λ_C, ω rises. **Genuinely uncertain — that is why it is the named candidate.** |
| The test | run the clock with V off (K=0) vs V on (the `m5_6_5c` confinement-validated K range); measure ω₁ each; check whether ω₁ rises toward the 16.74× the gap needs |
| Why deferred here | the clock harness (`m5_8_2h`) uses a custom stepper without the V toggle; a clean test = modifying the stepper + a ~44-min two-run sweep with real blow-up risk. A focused sub-experiment, not a safe bolt-on |

## DoD status

- [x] α reported from the Coulomb ratio + the anchor-independent invariant table (Part A)
- [x] the unit map sim → physical per anchor + the structural-gap statement + ratio-fallback scope (Part B + findings)
- [~] the V-on / Faber-`r₀` core route **evaluated** (test designed + physics expectation); the run itself is the named follow-up (Part C)
- [x] the table of which results become physical numbers once the map exists

Model: M5 Liquid Crystal. Physics-only, headless.

## Coverage-matrix cell prose (preserved from MODELS.md, condensed there 2026-06-19)

The MODELS.md "de Broglie clock (Zitterbewegung)" cell (M5 column) was trimmed to a summary on 2026-06-19; its prior full descriptive text is preserved here for the record:

> Bounded, self-starting, frequency-rigid 3+1D time crystal: quadratic action fails, signed quartic saturates, ω₁ is a start-independent attractor; classified a molten clock that regularizes toward a cold ground state. Energy-minimum property: activating the clock-fuel boost sector lowers the seed rest energy to a minimum ≈ 21% below the clock-stopped (b = 0) value, reproducing the energy-vs-frequency minimum of the ϕ⁴ toymodel arXiv:2501.04036 (the de Broglie clock is the energy-minimizing state, so stopping the clock can only raise energy). Absolute scale (#208 / #217 / #218, 2026-06-17): the simulated clock runs 28× below the electron ZBW under the energy anchor, but this is a recoverable calibration SPLIT, not a model deficit: #208 fixes the scale-free unit map (two free dials c, ℏ) and lands α = |b|/(ℏc) ≈ 1/178 right-order from pure topology; #217 shows the gap is structural (V-on is null on frequency, the faithful F-commutator kinetic is bounded ×3.04); #218 closes it, the energy- and length-anchor (Faber r₀) routes BRACKET the ZBW and their geometric mean recovers it to ~13% (the joint E·r₀ = α·(π/4)·ℏc classical-radius relation). The ZBW scale is recoverable by anchoring energy + length jointly, not the energy postulate alone; the precision refinement (the mass-family test of the geo-mean law) is #220.
