"""M6 sandbox v4 T7 — Shooting refinement for 4-fn bound state.

After T6 (IVP blowup) and T6→A (BVP finds 2-fn subset), this script
tries one more approach: find the bound state by treating "blowup
distance" as a continuous figure of merit. Solutions closer to the
true eigenvalue blow up LATER; refinement maximizes that distance.

ALGORITHM:
1. Fine 2D scan over (m_J², λ_bench) at fixed V₀=A₀=Q₀=J₀=0.1
2. For each (m_J², λ_bench), integrate IVP with bounded step + blowup
   event (threshold = 10.0, lower than T6's 100 for finer resolution)
3. Record r_blowup = r-value where |field| first exceeds threshold
   (or r_max if no blowup)
4. Identify the (m_J², λ_bench) region with largest r_blowup —
   that's where the bound state lives
5. Refine via golden-section / 2D Newton in that region

Then alternative: shoot on amplitude V₀=A₀=Q₀=J₀=A_0 with fixed
m_J², λ_bench. Vary A_0 looking for the value that maximizes
r_blowup at the candidate (m_J², λ_bench) point.

USAGE:
    python3 m6_v4_4fn_shoot.py --mode coarse   # 2D parameter scan
    python3 m6_v4_4fn_shoot.py --mode fine     # local refinement
    python3 m6_v4_4fn_shoot.py --mode amp      # amplitude shooting
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


OMEGA_DEF = 1.0
R_MIN = 0.05
R_MAX_DEF = 15.0
RTOL = 1e-8
ATOL = 1e-10
MAX_STEP = 0.02
BLOWUP_THRESHOLD = 10.0   # lower than T6's 100 → finer resolution

TARGET_HQ = 1.6969


def ode_4fn(r, y, m_eff_sq, lam_bench):
    V, dV, A, dA, Q, dQ, J, dJ = y
    inv_r = 1.0 / r
    QQ_JJ = Q * Q - J * J
    d2V = -inv_r * dV + Q
    d2A = -inv_r * dA + J
    d2Q = -inv_r * dQ + V + m_eff_sq * Q + lam_bench * Q * QQ_JJ
    d2J = -inv_r * dJ + A - m_eff_sq * J - lam_bench * J * QQ_JJ
    return [dV, d2V, dA, d2A, dQ, d2Q, dJ, d2J]


def integrate_with_blowup(V0, A0, Q0, J0, m_J_sq, lam_bench,
                           omega=OMEGA_DEF, r_max=R_MAX_DEF):
    """Integrate IVP and report r_blowup distance.
    Returns (r_blowup, r_array, fields_array, peak_amplitudes).
    """
    m_eff_sq = m_J_sq - omega * omega
    y0 = [V0, 0.0, A0, 0.0, Q0, 0.0, J0, 0.0]

    def blowup_event(r, y):
        m = max(abs(y[0]), abs(y[2]), abs(y[4]), abs(y[6]))
        return BLOWUP_THRESHOLD - m
    blowup_event.terminal = True
    blowup_event.direction = -1

    try:
        sol = solve_ivp(
            fun=lambda r, y: ode_4fn(r, y, m_eff_sq, lam_bench),
            t_span=(R_MIN, r_max), y0=y0,
            method="RK45", rtol=RTOL, atol=ATOL,
            max_step=MAX_STEP, events=blowup_event,
        )
    except Exception:
        return R_MIN, None, None, None

    if not sol.success:
        return R_MIN, None, None, None

    r_blowup = sol.t[-1]
    fields = sol.y
    peaks = (max(abs(fields[0])), max(abs(fields[2])),
             max(abs(fields[4])), max(abs(fields[6])))
    return r_blowup, sol.t, fields, peaks


def coarse_scan(omega=OMEGA_DEF, r_max=R_MAX_DEF,
                V0=0.1, A0=0.1, Q0=0.1, J0=0.1):
    """2D scan over (m_J², λ_bench). Find region of max blowup distance."""
    print(f"\n{'='*92}")
    print(f"COARSE SCAN — (m_J², λ_bench) at V₀=A₀=Q₀=J₀={V0}, ω={omega}")
    print(f"{'='*92}")
    print(f"Larger r_blowup = closer to bound state.\n")

    m_J_sq_grid = np.linspace(0.5, 8.0, 16)
    lam_grid = np.linspace(0.1, 10.0, 16)

    # Build matrix of r_blowup
    R_grid = np.zeros((len(m_J_sq_grid), len(lam_grid)))
    print(f"{'m_J²':>6} | " + " ".join(f"{l:5.2f}" for l in lam_grid))
    print("-" * (8 + 6 * len(lam_grid)))
    for i, m_J_sq in enumerate(m_J_sq_grid):
        row = []
        for j, lam in enumerate(lam_grid):
            r_b, _, _, _ = integrate_with_blowup(V0, A0, Q0, J0, m_J_sq, lam,
                                                  omega=omega, r_max=r_max)
            R_grid[i, j] = r_b
            row.append(f"{r_b:5.2f}")
        print(f"{m_J_sq:>6.2f} | " + " ".join(row))

    # Find max r_blowup
    i_max, j_max = np.unravel_index(np.argmax(R_grid), R_grid.shape)
    m_J_sq_best = m_J_sq_grid[i_max]
    lam_best = lam_grid[j_max]
    r_b_max = R_grid[i_max, j_max]
    print(f"\nMax r_blowup = {r_b_max:.3f} at (m_J²={m_J_sq_best:.2f}, "
          f"λ_bench={lam_best:.2f})")
    if r_b_max >= r_max - 0.01:
        print(f"  → REACHED r_max without blowup. Likely bound state!")
        print(f"  Run --mode fine for refinement, OR check field values.")

    return {
        "m_J_sq_grid": m_J_sq_grid.tolist(),
        "lam_grid": lam_grid.tolist(),
        "R_grid": R_grid.tolist(),
        "best": {"m_J_sq": float(m_J_sq_best), "lam_bench": float(lam_best),
                 "r_blowup": float(r_b_max)},
    }


def amplitude_shoot(m_J_sq, lam_bench, omega=OMEGA_DEF, r_max=R_MAX_DEF,
                    helicity="symmetric"):
    """Sweep common amplitude A_0; find A_0 maximizing r_blowup.

    helicity:
      'symmetric'   → V₀ = A₀ = Q₀ = J₀ = +A_0    (T7 original; Q_CS=0)
      'asymmetric'  → V₀ = Q₀ = +A_0, A₀ = J₀ = -A_0  (Werbos 2026-05-19;
                       Q_CS=1 helicity)
      'asymmetric2' → V₀ = Q₀ = -A_0, A₀ = J₀ = +A_0  (alternate sign;
                       Q_CS=-1 helicity, should also work)
    """
    print(f"\n{'='*100}")
    print(f"AMPLITUDE SHOOT — helicity={helicity}, (m_J²={m_J_sq}, "
          f"λ_bench={lam_bench}, ω={omega})")
    print(f"{'='*100}")
    A_grid = np.logspace(-3, 0, 30)  # 0.001 to 1.0
    print(f"{'A_0':>8} {'V₀':>8} {'A₀':>8} {'Q₀':>8} {'J₀':>8} "
          f"{'r_blowup':>10} {'note':>22}")
    print("-" * 84)
    best_r = R_MIN
    best_A = A_grid[0]
    for A in A_grid:
        if helicity == "symmetric":
            V0, A0, Q0, J0 = +A, +A, +A, +A
        elif helicity == "asymmetric":
            V0, A0, Q0, J0 = +A, -A, +A, -A
        elif helicity == "asymmetric2":
            V0, A0, Q0, J0 = -A, +A, -A, +A
        else:
            raise ValueError(f"Unknown helicity {helicity}")
        r_b, _, _, _ = integrate_with_blowup(V0, A0, Q0, J0,
                                              m_J_sq, lam_bench,
                                              omega=omega, r_max=r_max)
        note = ""
        if r_b >= r_max - 0.01:
            note = "REACHED r_max"
        elif r_b > best_r:
            best_r = r_b
            best_A = A
        print(f"{A:>8.4f} {V0:>+8.3f} {A0:>+8.3f} {Q0:>+8.3f} {J0:>+8.3f} "
              f"{r_b:>10.3f} {note:>22}")
    print(f"\nMax r_blowup = {best_r:.3f} at A_0 = {best_A:.4f}")
    return {"best_A": float(best_A), "r_blowup": float(best_r),
            "helicity": helicity, "m_J_sq": float(m_J_sq),
            "lam_bench": float(lam_bench)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mode", choices=["coarse", "amp"], default="coarse")
    parser.add_argument("--m-J-sq", type=float, default=2.0)
    parser.add_argument("--lam-bench", type=float, default=1.0)
    parser.add_argument("--omega", type=float, default=OMEGA_DEF)
    parser.add_argument("--r-max", type=float, default=R_MAX_DEF)
    parser.add_argument("--V0", type=float, default=0.1)
    parser.add_argument("--helicity",
                        choices=["symmetric", "asymmetric", "asymmetric2"],
                        default="asymmetric",
                        help="symmetric (V=A=Q=J=+A), asymmetric (V=Q=+A, "
                             "A=J=-A; Q_CS=+1), asymmetric2 (V=Q=-A, "
                             "A=J=+A; Q_CS=-1)")
    args = parser.parse_args()

    t0 = time.time()
    if args.mode == "coarse":
        result = coarse_scan(omega=args.omega, r_max=args.r_max,
                             V0=args.V0, A0=args.V0, Q0=args.V0, J0=args.V0)
    else:
        result = amplitude_shoot(args.m_J_sq, args.lam_bench,
                                  omega=args.omega, r_max=args.r_max,
                                  helicity=args.helicity)
    print(f"\nElapsed: {time.time() - t0:.1f}s")

    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / f"m6_v4_4fn_shoot_{args.mode}_results.json"
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2,
                  default=lambda o: float(o) if hasattr(o, "item") else str(o))
    print(f"[json] {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
