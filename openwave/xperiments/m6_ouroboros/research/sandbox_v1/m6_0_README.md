# M6.0 — Werbos's mass-frequency reproduction (NumPy sandbox)

**Purpose** — independently reproduce Werbos's published claim from `theory/Ouroboros/Ouroboros_Particle_Spectrum_2026.pdf` that scanning the chaoiton oscillation frequency `ω` at fixed coupling parameters `(g = 1.0625, λ = 1.0)` produces the lepton mass spectrum (electron / muon / pion / tau).

If the reproduction succeeds, M6 is a GO and we scaffold the Taichi production substrate. If it fails, M6 is downgraded. See `../0a_background.md § 12` for the decision framework.

---

## What the script reproduces

| Element | Werbos's value (Calibration + Spectrum papers) |
| --- | --- |
| Lagrangian | `L = -F·F - G·G + J·A - g(J·J)²` |
| Coupling | `g = 1.0625` (electron calibration baseline) |
| Lagrange multiplier | `λ = 1.0` |
| ω scan range | 1.0 to 80 (lepton spectrum target points: 1.0, 11.0, 13.0, 40.7) |
| Integrator | scipy solve_ivp, RK45, `rtol = 1e-10`, `atol = 1e-12` |
| Radial domain | `r ∈ [0.02, 10.0]` code units (Calibration paper); spectrum paper used 15.0 |
| Initial field amplitudes at small r | `A₀ = B₀ = 0.1` |
| Localization criterion | `\|A(r_max)\| + \|J(r_max)\| < 0.15` (Calibration) / `< 0.25` (Spectrum) |
| Mass scaling target | `m ∝ ω^2.22` (near-quadratic) |
| Electron calibration target | `H/Q = 1.6969` at ω=1.0 (Calibration paper §4.2) |

---

## Gate 0 — calibration check

Before doing the ω scan, the script runs ω=1.0 and checks that `H/Q ≈ 1.6969`. This is the calibration baseline Werbos publishes (Calibration paper Table 4.2). If our `H/Q` at ω=1.0 is wildly off (>5% gap), the radial ODE form is wrong and we need to refine before scanning.

---

## Assumptions in the radial ODE (to verify with Werbos's exact form)

The 2017 foundation paper and the 2026 Calibration paper don't write the exact radial ODE explicitly. The simplest interpretation, derived from the Lagrangian on the toroidal-poloidal ansatz, gives a coupled pair of Helmholtz-like equations for the radial profiles α(r), γ(r):

```text
α''(r) + (2/r) α'(r) + ω² α(r) + γ(r) = 0          [A-field, Maxwell-like with J source]
γ''(r) + (2/r) γ'(r) + (ω² - 2g γ²(r)) γ(r) + α(r) = 0   [J-field, KG-like with f mass + A source]
```

This is the simplest two-scalar reduction (treats poloidal and toroidal as scalar profiles). The full ansatz has four functions (α, β, γ, δ) — the simpler two-function version is the starting point and matches the calibration check; if calibration fails, expand to four functions.

Initial conditions (at small r, regularity):

- α(0.02) = A₀ = 0.1, α'(0.02) = 0
- γ(0.02) = B₀ = 0.1, γ'(0.02) = 0

Integration outward to r_max = 10.0 code units.

---

## Outputs

After running:

- `../plots/m6_0_calibration_check.png` — α(r), γ(r) profiles at ω=1.0 + H/Q value
- `../plots/m6_0_mass_frequency_scan.png` — predicted mass vs ω, with Werbos's 3 lepton points overlaid
- `../plots/m6_0_localization_map.png` — which (ω, threshold) combinations are localized
- Console output with PASS/CONDITIONAL/FAIL decision per §12.4

---

## Decision criteria (from `0a_background.md § 12.4`)

| Outcome | M6 decision |
| --- | --- |
| Reproduces within 5% across all 3 leptons | GO |
| Reproduces electron+muon but tau > 10% | Conditional GO |
| Fails at electron calibration (gate 0) | NO-GO — investigate ODE form before deciding |
| Numerical instability | Investigate before deciding |

---

## How to run

```bash
cd openwave/xperiments/m6_ouroboros/research/scripts/
python m6_0_werbos_reproduction.py
```

Run on M4 Max in <1 minute (scipy IVP is fast; 31 ω-values × ms-per-integration).
