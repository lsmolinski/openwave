# M6 sandbox v3 — Lean-ODE-corrected chaoiton

**Status:** 2026-05-18. Supersedes sandbox v2's attempts (IVP locked-ansatz, BVP locked-ansatz).

---

## What changed from v2 (sources: new theory files received 2026-05-18 17:20)

| Issue | v2 diagnosis | v3 correction |
| --- | --- | --- |
| ODE form | 2D scalar (1/r) or 3D scalar (2/r) Laplacian | VECTOR Laplacian: f'' + (1/r)f' − f/r² (from chaoiton_theorem.lean.txt) |
| Origin BC | α(r_min) = A₀ (field value) | α(r_min) = A₀·r_min, α'(r_min) = A₀ (slope — field starts at ZERO) |
| Q=0 approach | Q_CS=0 via locked ansatz J=±A | α=0 exactly + ω=0 (A-field off, no angular momentum) |
| Linear damping | not present | −λβ in β equation (Lagrange multiplier) |
| Coupling coefficient | quarter (c/4) | full (β sources α, α sources β with coefficient 1) |

## Canonical ODE (from chaoiton_theorem.lean.txt)

```text
α'' + (1/r)α' − (1/r²)α + ω²α = β
β'' + (1/r)β' − (1/r²)β + ω²β = α − λβ − 4gβ³
```

BCs at origin: α ~ A₀·r, β ~ B₀·r (slopes, not values)
Localization: |α(r)| + |β(r)| < 0.1 for r ≥ 8

## Neutral chaoiton (Q=0 dark matter)

With α=0, ω=0, the β equation becomes:

```text
β'' + (1/r)β' − β/r² + λβ + 4gβ³ = 0
```

Single nonlinear ODE. BC: β ~ B₀·r at origin. Tune B₀ for localization.

From `Claude bosonization May18 4PM.txt`: 1,208 solutions found in this sector.
Mass scales as `m_χ ≈ λ × 0.511 MeV` for the lightest solutions.

## v3 first-run results (2026-05-18)

```text
Neutral chaoiton scan (α=0, ω=0):
  λ=1.0, B0=0.010 → m_χ = 0.508 MeV  ← matches λ × 0.511 MeV to 0.6%
  λ=0.5, B0=0.010 → m_χ = 0.729 MeV
  λ=0.1, B0=0.010 → m_χ = 1.55 MeV

23 localized solutions across (λ, B0) parameter space. ✅

Electron calibration (ω=1.0, g=1.0625, λ=1.0):
  Best H/Q = 1.723 (A0=2, B0=10) — 1.56% gap from target 1.6969
  None yet localized (tail > 0.1). Energy functional needs tuning. ⚠️
```

## Scripts

```text
m6_v3_chaoiton.py   — main: electron calibration + neutral scan
```

## Next steps

1. Fix energy functional for charged chaoiton (calibration localization issue)
2. Add stability test (Gelfand-Fomin conjugate-point per Lean theorem)
3. Electron + neutral confirmed → email Werbos with v3 results + Q8-Q12
