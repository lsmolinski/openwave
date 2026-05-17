"""M6.0 iterations — try multiple ODE / H / Q variations to find calibration match.

The first-pass reproduction (m6_0_werbos_reproduction.py) hit a 192% gap on the
electron H/Q calibration. This script systematically tries physically-motivated
variations to find which form reproduces Werbos's published H/Q = 1.6969 at
(g = 1.0625, omega = 1.0).

Each variation is a separate function so they can be mixed. The script runs all
and tabulates the H/Q result for ranking. Best variation becomes the canonical
form for the production reproduction script.

Variations explored:

  ODE forms:
    A1. Original: alpha'' + (2/r)alpha' + omega² alpha + gamma = 0
                  gamma'' + (2/r)gamma' + (omega² - 2g gamma²) gamma + alpha = 0
    A2. Sign flip on coupling: ...− gamma = 0, ... − alpha = 0
    A3. Factor-4 in cubic: 4g gamma² in mass term (matches df/dgamma=4g gamma³)
    A4. Quarter coupling: gamma/4 (matches 4× Lagrangian normalization)
    A5. Static A (drop omega² alpha): treats A as Coulomb-like static

  Charge densities:
    Q1. omega * gamma²              (single-field, original)
    Q2. omega * (alpha² + gamma²)   (both fields share U(1))
    Q3. omega * alpha²              (alpha-only, as if J is non-charged)
    Q4. gamma  (Coulomb J^0 directly, ignoring oscillation factor)

  Hamiltonian densities:
    H1. Original (no 1/2 factor)
    H2. With 1/2 factor on all kinetic + gradient + mass terms
    H3. With 1/4 factor (matches Werbos's L = -F·F not -F·F/4 normalization)
"""

from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

# Werbos's published values
G_COUPLING = 1.0625
A_INIT = 0.1
B_INIT = 0.1
R_MIN = 0.02
R_MAX = 10.0
N_GRID = 800
RTOL = 1.0e-10
ATOL = 1.0e-12

TARGET_HQ = 1.6969


# ---------------------------------------------------------------------------
# ODE variants
# ---------------------------------------------------------------------------

def ode_A1(r, y, omega, g):
    """Original: +gamma coupling on alpha, 2g gamma^2 mass."""
    a, da, c, dc = y
    inv_r = 1.0 / r
    d2a = -2*inv_r*da - omega*omega*a - c
    d2c = -2*inv_r*dc - (omega*omega - 2*g*c*c)*c - a
    return [da, d2a, dc, d2c]


def ode_A2(r, y, omega, g):
    """Sign flip on couplings: -gamma in alpha eq, -alpha in gamma eq."""
    a, da, c, dc = y
    inv_r = 1.0 / r
    d2a = -2*inv_r*da - omega*omega*a + c
    d2c = -2*inv_r*dc - (omega*omega - 2*g*c*c)*c + a
    return [da, d2a, dc, d2c]


def ode_A3(r, y, omega, g):
    """4g gamma^3 mass term (proper f'=2gs, df/dgamma=4g gamma^3)."""
    a, da, c, dc = y
    inv_r = 1.0 / r
    d2a = -2*inv_r*da - omega*omega*a - c
    d2c = -2*inv_r*dc - omega*omega*c + 4*g*c*c*c - a
    return [da, d2a, dc, d2c]


def ode_A4(r, y, omega, g):
    """Quarter coupling (4x Lagrangian normalization)."""
    a, da, c, dc = y
    inv_r = 1.0 / r
    d2a = -2*inv_r*da - omega*omega*a - c/4
    d2c = -2*inv_r*dc - omega*omega*c + 4*g*c*c*c - a/4
    return [da, d2a, dc, d2c]


def ode_A5(r, y, omega, g):
    """Static A (drop omega² alpha — A is Coulomb-like, not oscillating)."""
    a, da, c, dc = y
    inv_r = 1.0 / r
    d2a = -2*inv_r*da - c                    # static Poisson for A
    d2c = -2*inv_r*dc - omega*omega*c + 4*g*c*c*c - a
    return [da, d2a, dc, d2c]


def ode_A6(r, y, omega, g):
    """A5 + sign flip on couplings."""
    a, da, c, dc = y
    inv_r = 1.0 / r
    d2a = -2*inv_r*da + c                    # static Poisson with sign flip
    d2c = -2*inv_r*dc - omega*omega*c + 4*g*c*c*c + a
    return [da, d2a, dc, d2c]


ODE_VARIANTS = {
    "A1": ode_A1,
    "A2": ode_A2,
    "A3": ode_A3,
    "A4": ode_A4,
    "A5": ode_A5,
    "A6": ode_A6,
}


# ---------------------------------------------------------------------------
# Charge density variants
# ---------------------------------------------------------------------------

def q_density_Q1(omega, a, c, g):
    return omega * c * c

def q_density_Q2(omega, a, c, g):
    return omega * (a*a + c*c)

def q_density_Q3(omega, a, c, g):
    return omega * a * a

def q_density_Q4(omega, a, c, g):
    return c

Q_VARIANTS = {
    "Q1": q_density_Q1,
    "Q2": q_density_Q2,
    "Q3": q_density_Q3,
    "Q4": q_density_Q4,
}


# ---------------------------------------------------------------------------
# Hamiltonian density variants
# ---------------------------------------------------------------------------

def h_density_H1(omega, a, da, c, dc, g):
    return da*da + omega*omega*a*a + dc*dc + omega*omega*c*c + g*c*c*c*c - a*c

def h_density_H2(omega, a, da, c, dc, g):
    """With 1/2 factor on kinetic+gradient terms."""
    return 0.5*(da*da + omega*omega*a*a + dc*dc + omega*omega*c*c) + g*c*c*c*c - a*c

def h_density_H3(omega, a, da, c, dc, g):
    """1/4 factor."""
    return 0.25*(da*da + omega*omega*a*a + dc*dc + omega*omega*c*c) + g*c*c*c*c - a*c

def h_density_H4(omega, a, da, c, dc, g):
    """H2 with positive coupling sign in H."""
    return 0.5*(da*da + omega*omega*a*a + dc*dc + omega*omega*c*c) + g*c*c*c*c + a*c

H_VARIANTS = {
    "H1": h_density_H1,
    "H2": h_density_H2,
    "H3": h_density_H3,
    "H4": h_density_H4,
}


# ---------------------------------------------------------------------------
# Integrate + measure
# ---------------------------------------------------------------------------

def integrate_and_measure(ode_fn, h_fn, q_fn, omega, g=G_COUPLING):
    y0 = [A_INIT, 0.0, B_INIT, 0.0]
    r_eval = np.linspace(R_MIN, R_MAX, N_GRID)
    sol = solve_ivp(
        fun=lambda r, y: ode_fn(r, y, omega, g),
        t_span=(R_MIN, R_MAX),
        y0=y0,
        t_eval=r_eval,
        method="RK45",
        rtol=RTOL,
        atol=ATOL,
    )
    if not sol.success:
        return None

    r = sol.t
    a, da, c, dc = sol.y[0], sol.y[1], sol.y[2], sol.y[3]

    # tail check (Werbos's localization criterion)
    tail = abs(a[-1]) + abs(c[-1])

    # Integrals
    h_dens = h_fn(omega, a, da, c, dc, g)
    q_dens = q_fn(omega, a, c, g)
    r2 = r * r
    H = 4.0 * np.pi * np.trapezoid(h_dens * r2, r)
    Q = 4.0 * np.pi * np.trapezoid(q_dens * r2, r)

    hq = H / Q if abs(Q) > 1e-12 else float("nan")

    return {
        "H": float(H),
        "Q": float(Q),
        "HQ": float(hq),
        "tail": float(tail),
        "localized": bool(tail < 0.15),
        "field_max": float(max(np.abs(a).max(), np.abs(c).max())),
    }


# ---------------------------------------------------------------------------
# Full grid scan
# ---------------------------------------------------------------------------

def main():
    rows = []
    for a_key, a_fn in ODE_VARIANTS.items():
        for h_key, h_fn in H_VARIANTS.items():
            for q_key, q_fn in Q_VARIANTS.items():
                res = integrate_and_measure(a_fn, h_fn, q_fn, omega=1.0)
                if res is None:
                    rows.append((a_key, h_key, q_key, float("nan"), float("nan"),
                                 float("nan"), float("nan"), False, "FAIL"))
                    continue
                gap = abs(res["HQ"] - TARGET_HQ) / TARGET_HQ
                rows.append((a_key, h_key, q_key, res["H"], res["Q"], res["HQ"],
                             res["tail"], res["localized"], gap))

    # Sort by gap from target
    valid_rows = [r for r in rows if not (isinstance(r[8], str) or np.isnan(r[8]))]
    valid_rows.sort(key=lambda r: r[8])

    print(f"\n{'='*78}")
    print(f"ITERATION SWEEP — searching for H/Q ≈ {TARGET_HQ} at omega=1.0, g={G_COUPLING}")
    print(f"{'='*78}\n")
    print(f"{'ODE':>4} {'H':>3} {'Q':>3} {'H_int':>10} {'Q_int':>10} "
          f"{'H/Q':>10} {'tail':>8} {'loc':>5} {'gap%':>8}")
    print("-" * 78)

    # Print top 30 (best matches)
    for row in valid_rows[:30]:
        a_key, h_key, q_key, H, Q, hq, tail, loc, gap = row
        gap_str = f"{gap*100:7.2f}" if isinstance(gap, float) else str(gap)
        loc_str = "Y" if loc else "n"
        print(f"{a_key:>4} {h_key:>3} {q_key:>3} {H:10.4f} {Q:10.4f} "
              f"{hq:10.4f} {tail:8.4f} {loc_str:>5} {gap_str}")

    print(f"\n{'='*78}")
    print("Best match (lowest gap, localized):")
    print(f"{'='*78}")
    best_localized = [r for r in valid_rows if r[7]]
    if best_localized:
        a_key, h_key, q_key, H, Q, hq, tail, loc, gap = best_localized[0]
        print(f"\n  Variant:  ODE={a_key}, H={h_key}, Q={q_key}")
        print(f"  H/Q:      {hq:.4f} (target {TARGET_HQ})")
        print(f"  Gap:      {gap*100:.2f}%")
        print(f"  Tail:     {tail:.4f} (threshold 0.15)")
        print(f"  Verdict:  {'PASS' if gap < 0.05 else 'CLOSE' if gap < 0.20 else 'FAIL'}")
    else:
        print("\n  No localized solution found across all variants.")

    return valid_rows


if __name__ == "__main__":
    main()
