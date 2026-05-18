# M6.1 — Q=0 neutral chaoiton (sandbox v2)

**Purpose** — compute the ground-state mass `m_χ` of the lightest neutral chaoiton (`Q_CS = 0`), per Werbos's 2026-05-17 ask. This is the missing piece for the Ouroboros dark-matter prediction `Ω_χ h² ≈ 0.12` (Planck).

If `m_χ ≥ MeV` → standard freeze-out gives the cold-DM relic density.
If `m_χ ≪ MeV` → dark-radiation candidate.

See `../0a_background.md § 12` for the decision framework and §12.3 for the full ask.

---

## What's new vs sandbox v1

`sandbox_v1/` reproduced electron H/Q with a 2-function, 3D-Laplacian, single-field-cubic ODE. Werbos's 2026-05-18 benchmark document (`theory/Numerical Benchmark for the Ouroboros Lagrangian.docx`) writes the radial ODE explicitly for the first time:

```text
Component          | Benchmark form (v2)
-------------------|--------------------------------------------------
Laplacian          | 2D radial: Δ_r f = f''(r) + (1/r) f'(r)
Function count     | 4 (V, A, Q, J)
V equation         | Δ_r V = Q
A equation         | Δ_r A = J
Q equation         | Δ_r Q = V + m_J² Q + λ Q(Q² − J²)
J equation         | Δ_r J = A − m_J² J − λ J(Q² − J²)
f(s)               | (1/2) m_J² s + (λ/4) s²
Time ansatz        | V·sin(ωt), A·cos(ωt), Q·cos(ωt), J·sin(ωt)
BCs at r=0         | V'(0) = A'(0) = Q'(0) = J'(0) = 0
BCs at r → ∞       | fields → 0, exponential K_0(κr) decay
```

This script implements the full 4-function ODE; the sandbox v1 single-field reduction is structurally **incorrect** at this level.

---

## The locked-ansatz reduction

For `Q_CS = 0`, Werbos suggests the locked ansatz:

```text
J(r) = −A(r)
Q(r) = −V(r)
```

Verifies `Q_CS = (2πR) ∫ r dr [A·J' − J·A']`:
With `J = −A`, `J' = −A'`, this gives `∫ r dr [A·(−A') − (−A)·A'] = ∫ r dr [−AA' + AA'] = 0`. ✅

The locked ansatz is used as the **initial condition** (at `r = r_min`), not as an enforced identity throughout. The full 4-function ODE is integrated outward and we observe whether the constraint stays approximately satisfied.

---

## Open structural questions (carried over from §12.4)

The benchmark document doesn't fully pin down the ODE form. This script exposes them as parameters to sweep:

```text
Question                    | Parameter        | Default | Sweep range
----------------------------|------------------|---------|---------------
How does ω enter?           | OMEGA_INSERTION  | "none"  | none / m_J(ω) / explicit
What formula gives m_J(g)?  | M_J_SQUARED      | 0.0     | 0, g, g², ω², ωg
2D or 3D Laplacian?         | LAPLACIAN_DIM    | 2       | 2 (benchmark) / 3 (sandbox v1)
Source structure            | (locked in code) | V↔Q,A↔J | (benchmark spec)
```

Each open question maps to a CLI parameter. Run with defaults first; sweep when something looks off.

---

## Outputs

- `../plots/m6_1_neutral_chaoiton_profile.png` — V(r), A(r), Q(r), J(r) profiles
- `../plots/m6_1_neutral_chaoiton_energy.png` — energy density vs r
- `../plots/m6_1_neutral_chaoiton_summary.json` — m_χ value + parameters + convergence
- Console summary with m_χ in code units AND MeV (via R^phys = 191 fm)

---

## Decision criteria

```text
Outcome                                | Interpretation
---------------------------------------|--------------------------------
m_χ in 0.5–10 MeV (Werbos's quartic-   | Cold dark matter candidate.
  raised prediction)                   | Compute σv → Ω_χ h² freeze-out.
                                       |   Compare to Planck 0.120.
m_χ in 0.003–0.015 MeV (Werbos's       | Dark radiation candidate.
  linearized estimate)                 | Different cosmological signature.
m_χ < 1 keV or > 100 MeV               | Out of expected range; check
                                       |   parameter sweep + ODE form.
No localized solution found            | Locked ansatz incompatible with
                                       |   the full 4-fn ODE; need to
                                       |   relax constraint or change m_J.
```

---

## How to run

```bash
cd openwave/xperiments/m6_ouroboros/research/sandbox_v2/
python m6_1_neutral_chaoiton.py

# Sweep open questions:
python m6_1_neutral_chaoiton.py --m-j-sq g       # m_J² = g
python m6_1_neutral_chaoiton.py --m-j-sq omega2  # m_J² = ω²
python m6_1_neutral_chaoiton.py --laplacian 3    # sandbox-v1 convention
```

Expected runtime: <1 minute on M4 Max per parameter point.

---

## Reporting back to Werbos

When this script produces a clean m_χ value, the email back to Werbos should include:

- The m_χ number (scaled units AND MeV)
- The convergence-check result (H/Q stable to <0.1% on doubling N_r and r_max)
- The L/Q = ω sanity-check value
- The four open structural questions (§12.4 of background) we had to make choices on
- Plot of V(r), A(r) profiles
- Decision branch (cold DM vs dark radiation)

Per §12.8 collaboration norm: deliver Priority A alone before opening the muon/tau/proton/Sawada threads.
