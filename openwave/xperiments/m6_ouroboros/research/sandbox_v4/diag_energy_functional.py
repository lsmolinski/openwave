"""Diagnostic for T1 H-scaling problem.

T1 first run (--quick) found localized solutions but with H decreasing in ω,
opposite to Werbos's bosonization-session trend (m_μ/m_e = 207 at ω=12.78/1.0).

Hypothesis: our energy density is missing the time-derivative term
ω²(α²+β²) from |∂_t Aμ|² in the Ouroboros Lagrangian for an
e^{-iωt} ansatz. Without it, H comes only from the spatial-gradient
+ potential parts, which shrink as the chaoiton becomes more localized
at high ω (decay scale ~1/ω).

Test:
    Run integrate twice per (ω, A0, B0) — current h_dens, and
    h_dens + ω²(α²+β²). Pin (A0, B0)=(0.1, 0.1) (Werbos's calib point)
    and sweep ω. Check whether the "+ω²" variant gives:
      - H_e ≈ 0.494 at ω=1.0
      - H_μ/H_e ≈ 207 at ω=12.78
      - H_τ/H_e ≈ 3477 at ω≈40.7
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

G_COUPLING = 1.0625
LAMBDA = 1.0
R_MIN = 0.01
R_MAX = 20.0
N_GRID = 1000
RTOL = 1e-10
ATOL = 1e-12
LOC_R_CHECK = 8.0
LOC_THRESHOLD = 0.1
M_E_MEV = 0.51100
H_CODE_E_CALIB = 0.494


def ode_lean(r, y, g, lam, omega):
    a, da, b, db = y
    inv_r = 1.0 / r
    inv_r2 = inv_r * inv_r
    w2 = omega * omega
    d2a = b - w2 * a - inv_r * da + a * inv_r2
    d2b = a - lam * b - 4.0 * g * b * b * b - w2 * b - inv_r * db + b * inv_r2
    return [da, d2a, db, d2b]


def integrate_both(A0, B0, omega, g=G_COUPLING, lam=LAMBDA, r_max=R_MAX):
    y0 = [A0 * R_MIN, A0, B0 * R_MIN, B0]
    r_eval = np.linspace(R_MIN, r_max, N_GRID)
    sol = solve_ivp(
        fun=lambda r, y: ode_lean(r, y, g, lam, omega),
        t_span=(R_MIN, r_max), y0=y0, t_eval=r_eval,
        method="RK45", rtol=RTOL, atol=ATOL,
    )
    if not sol.success:
        return None
    r = sol.t
    a, da, b, db = sol.y
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return None

    w2 = omega * omega
    mask = r >= LOC_R_CHECK
    tail = float((np.abs(a[mask]) + np.abs(b[mask])).max())

    h_spatial = (da * da + a * a / (r * r)
                 + db * db + b * b / (r * r)
                 - a * b + lam * b * b / 2.0 + g * b ** 4)
    h_kinetic = w2 * (a * a + b * b)
    H_spatial = float(4.0 * np.pi * np.trapezoid(h_spatial * r * r, r))
    H_full = float(4.0 * np.pi * np.trapezoid((h_spatial + h_kinetic) * r * r, r))

    q_dens = omega * (a * a + b * b)
    Q = float(4.0 * np.pi * np.trapezoid(q_dens * r * r, r))

    return {
        "omega": omega, "A0": A0, "B0": B0,
        "H_spatial": H_spatial, "H_full": H_full, "Q": Q,
        "HQ_spatial": H_spatial / Q if Q > 1e-14 else float("nan"),
        "HQ_full": H_full / Q if Q > 1e-14 else float("nan"),
        "tail": tail,
        "loc_strict": tail < LOC_THRESHOLD,
        "loc_loose": tail < 0.15,   # matches Werbos's "0.137 = localized"
    }


def main():
    print("=" * 90)
    print("DIAGNOSTIC: H scaling vs ω at pinned (A0, B0) = (0.1, 0.1)")
    print("=" * 90)
    print(f"{'ω':>6} {'H_spat':>10} {'H_full':>10} {'HQ_spat':>9} {'HQ_full':>9} "
          f"{'tail':>8} {'strict':>7} {'loose':>6}")
    print("-" * 90)

    omegas = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 12.78, 15.0, 20.0,
              25.0, 30.0, 35.0, 40.7, 45.0, 50.0]
    results = []
    for omega in omegas:
        r = integrate_both(0.1, 0.1, omega)
        if r is None:
            print(f"{omega:>6.2f}  [integration failed]")
            continue
        results.append(r)
        print(f"{omega:>6.2f} {r['H_spatial']:>10.4f} {r['H_full']:>10.4f} "
              f"{r['HQ_spatial']:>9.4f} {r['HQ_full']:>9.4f} "
              f"{r['tail']:>8.4f} {str(r['loc_strict']):>7s} {str(r['loc_loose']):>6s}")

    # Calibration: at ω=1.0 we want H_full or H_spat ≈ 0.494 (electron)
    print("\n" + "=" * 90)
    print("CALIBRATION CHECK at ω=1.0 (electron baseline)")
    print("=" * 90)
    e_row = next((r for r in results if abs(r['omega'] - 1.0) < 0.01), None)
    if e_row:
        print(f"H_spatial = {e_row['H_spatial']:.4f}  (target ≈ 0.494)")
        print(f"H_full    = {e_row['H_full']:.4f}  (target ≈ 0.494 IF ω-term is included)")
        print(f"Bosonization log says H_code=11.4, Q=8.7, H/Q=1.31 at this point")

    # Ratio check
    print("\n" + "=" * 90)
    print("MASS-RATIO CHECK (using H ratios as proxy for m ratios)")
    print("=" * 90)
    print(f"{'particle':>11} {'target m/m_e':>14} {'H_spat ratio':>14} "
          f"{'H_full ratio':>14}")
    print("-" * 65)
    if e_row and e_row['H_spatial'] > 0:
        for omega_t, name, mass_t in [(12.78, "muon", 105.66),
                                       (15.0, "pion+", 139.57),
                                       (40.7, "tau", 1776.86)]:
            row = next((r for r in results if abs(r['omega'] - omega_t) < 0.01), None)
            if row is None:
                continue
            target_ratio = mass_t / M_E_MEV
            spat_ratio = row['H_spatial'] / e_row['H_spatial']
            full_ratio = row['H_full'] / e_row['H_full']
            print(f"{name:>11} {target_ratio:>14.2f} {spat_ratio:>14.2f} "
                  f"{full_ratio:>14.2f}")

    # Also try wider (A0, B0) at ω=1.0 to see if localization is achievable
    print("\n" + "=" * 90)
    print("ω=1.0 LOCALIZATION SEARCH (wider amplitude grid)")
    print("=" * 90)
    print(f"{'A0':>6} {'B0':>6} {'H_spat':>10} {'tail':>8} {'strict':>7} {'loose':>6}")
    print("-" * 60)
    for A0 in [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        for B0 in [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
            r = integrate_both(A0, B0, 1.0)
            if r is None:
                continue
            flag = "Y" if r['loc_strict'] else ("loose" if r['loc_loose'] else "")
            if flag:
                print(f"{A0:>6.2f} {B0:>6.2f} {r['H_spatial']:>10.4f} "
                      f"{r['tail']:>8.4f} {str(r['loc_strict']):>7s} "
                      f"{str(r['loc_loose']):>6s}")


if __name__ == "__main__":
    main()
