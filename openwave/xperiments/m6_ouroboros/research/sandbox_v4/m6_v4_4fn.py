"""M6 sandbox v4 — 4-function benchmark ODE for charged sector (2026-05-19).

Implements the canonical 4-function (V, A, Q, J) radial ODE from Werbos's
Numerical Benchmark §3, with the m_eff² = m_J² − ω² substitution Werbos
clarified in his 2026-05-19 1:49 PM reply.

ODE SYSTEM (toroidal Laplacian Δ_r f = f''+f'/r, NO -f/r²):
    Δ_r V  =  Q
    Δ_r A  =  J
    Δ_r Q  =  V  +  m_eff²·Q  +  λ·Q·(Q² − J²)
    Δ_r J  =  A  −  m_eff²·J  −  λ·J·(Q² − J²)

    m_eff² = m_J² − ω²   (Werbos clarification 1)

BCs (value at origin, derivative zero):
    V(r_min) = V₀,  V'(r_min) = 0
    A(r_min) = A₀,  A'(r_min) = 0
    Q(r_min) = Q₀,  Q'(r_min) = 0
    J(r_min) = J₀,  J'(r_min) = 0

CANONICAL CALIBRATION ANCHOR (electron):
    g_v5 = 1.0625,  λ_v5 = 1.0,  ω = 1.0
    V₀ = A₀ = Q₀ = J₀ = 0.1
    Target: H/Q = 1.6969  (Werbos, 0.56% off mass-to-charge target 1.6875)

ENERGY FUNCTIONAL (benchmark §5, R=1 per Werbos clarification 3):
    H = (2π)² ∫ r dr [
        (V')² + (A')² + (Q')² + (J')²
      − V·Q + A·J
      + (m_J²/2)·(Q² − J²)
      + (λ/4)·(Q² − J²)²
    ]

    Optional: include time-derivative kinetic term ω²/2·(V²+A²+Q²+J²)
    via --include-omega-kinetic flag (consistent with Hamiltonian for
    e^{-iωt} time dependence; benchmark doc may have absorbed this).

CHARGE FUNCTIONAL (Chern-Simons linking, R=1):
    Q_CS = 2π ∫ r dr [A·J' − J·A']

    Identity: L/Q = ω  (Noether)
    Calibration target: Q_CS = 1 (electron, integer winding)

Werbos's f(s) parametrization mapping (clarification 2):
    v5 paper:    f(s) = g·s²
    benchmark:   f(s) = ½·m_J²·s + ¼·λ·s²
    Equivalent when m_J² = 0, λ = 4g
    For g_v5 = 1.0625:  λ_bench = 4·g_v5 = 4.25 (when m_J² = 0)
    But m_J² = 0 + ω = 1 → m_eff² = -1 (oscillates, no decay)
    → for charged calibration, m_J² > ω² required for K_0 decay.
    Real calibration value of m_J² is unknown; scan to find it.

USAGE:
    python3 m6_v4_4fn.py --mode calibration  # scan (m_eff², λ) at ω=1.0
    python3 m6_v4_4fn.py --mode single --m-eff-sq 1.0 --lam-bench 4.25
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

G_V5 = 1.0625        # Lean / v5 g parameter
LAMBDA_V5 = 1.0      # Lean / v5 λ parameter (linear coupling)
OMEGA_DEF = 1.0

R_MIN = 0.05            # larger than 0.01 to ease 1/r singularity
R_MAX_DEF = 20.0
N_GRID = 1500
RTOL = 1e-8             # relaxed from 1e-10 for stiffness tolerance
ATOL = 1e-10
SOLVE_METHOD = "RK45"
MAX_STEP = 0.05         # cap step size to prevent overshoot of bound trajectory
BLOWUP_THRESHOLD = 100.0  # |field| > this → terminate integration as unbound

LOC_R_CHECK = 8.0
LOC_THRESHOLD = 0.1

TARGET_HQ = 1.6969          # Werbos electron-calibration anchor
TARGET_HQ_PHYS = 1.6875     # m_e / e_nat (physical electron target)
TARGET_LQ = 1.0             # ω = 1.0 → L/Q = ω
TARGET_QCS = 1.0            # integer winding


# ---------------------------------------------------------------------------
# ODE: 4-function benchmark with m_eff² substitution
# ---------------------------------------------------------------------------

def ode_4fn(r, y, m_eff_sq, lam_bench):
    """State: y = [V, V', A, A', Q, Q', J, J']

    Equations:
        V'' = -V'/r + Q
        A'' = -A'/r + J
        Q'' = -Q'/r + V + m_eff²·Q + λ·Q·(Q² − J²)
        J'' = -J'/r + A − m_eff²·J − λ·J·(Q² − J²)
    """
    V, dV, A, dA, Q, dQ, J, dJ = y
    inv_r = 1.0 / r

    QQ_JJ = Q * Q - J * J

    d2V = -inv_r * dV + Q
    d2A = -inv_r * dA + J
    d2Q = -inv_r * dQ + V + m_eff_sq * Q + lam_bench * Q * QQ_JJ
    d2J = -inv_r * dJ + A - m_eff_sq * J - lam_bench * J * QQ_JJ

    return [dV, d2V, dA, d2A, dQ, d2Q, dJ, d2J]


def integrate(V0, A0, Q0, J0, omega, m_J_sq, lam_bench,
              r_max=R_MAX_DEF, n_grid=N_GRID,
              include_omega_kinetic=True):
    """Integrate 4-fn ODE from r=R_MIN with value BCs.

    Werbos: m_eff² = m_J² − ω²
    Returns dict with H, Q_CS, H/Q, L/Q, localization, etc.
    """
    m_eff_sq = m_J_sq - omega * omega
    y0 = [V0, 0.0, A0, 0.0, Q0, 0.0, J0, 0.0]
    r_eval = np.linspace(R_MIN, r_max, n_grid)

    # Termination event: stop integration if any field exceeds threshold
    # (indicates unbound trajectory — saves time vs running to r_max)
    def blowup_event(r, y):
        m = max(abs(y[0]), abs(y[2]), abs(y[4]), abs(y[6]))
        return BLOWUP_THRESHOLD - m
    blowup_event.terminal = True
    blowup_event.direction = -1

    try:
        sol = solve_ivp(
            fun=lambda r, y: ode_4fn(r, y, m_eff_sq, lam_bench),
            t_span=(R_MIN, r_max), y0=y0, t_eval=r_eval,
            method=SOLVE_METHOD, rtol=RTOL, atol=ATOL,
            max_step=MAX_STEP, events=blowup_event,
        )
    except Exception as e:
        return {"success": False, "error": str(e)[:60],
                "m_J_sq": m_J_sq, "lam_bench": lam_bench,
                "m_eff_sq": m_eff_sq}
    if not sol.success:
        return {"success": False, "error": sol.message[:60],
                "m_J_sq": m_J_sq, "lam_bench": lam_bench,
                "m_eff_sq": m_eff_sq}
    # Detect early termination (blowup event)
    blew_up = sol.t[-1] < r_max - 0.01

    r = sol.t
    V, dV, A, dA, Q, dQ, J, dJ = sol.y

    fields_finite = all(np.all(np.isfinite(f)) for f in (V, A, Q, J))
    if not fields_finite:
        return {"success": False, "error": "non-finite fields"}

    # Localization check (any field exceeding threshold at large r)
    mask = r >= LOC_R_CHECK
    if mask.sum() == 0:
        # integration didn't reach LOC_R_CHECK (early blowup termination)
        max_abs = (np.abs(V) + np.abs(A) + np.abs(Q) + np.abs(J)).max()
        tail = float(max_abs)
        localized = False
    else:
        max_abs = (np.abs(V[mask]) + np.abs(A[mask])
                   + np.abs(Q[mask]) + np.abs(J[mask])).max()
        tail = float(max_abs)
        localized = bool(tail < LOC_THRESHOLD) and not blew_up

    # Peak amplitudes (sanity check: are fields nontrivial?)
    peak = max(np.abs(V).max(), np.abs(A).max(),
               np.abs(Q).max(), np.abs(J).max())

    # Energy functional (benchmark §5, R = 1)
    QQ_JJ = Q * Q - J * J
    h_dens = (dV * dV + dA * dA + dQ * dQ + dJ * dJ
              - V * Q + A * J
              + 0.5 * m_J_sq * QQ_JJ
              + 0.25 * lam_bench * QQ_JJ * QQ_JJ)
    if include_omega_kinetic:
        # Time-derivative kinetic term (Hamiltonian for e^{-iωt}):
        # ½·ω²·(V² + A² + Q² + J²)
        h_dens = h_dens + 0.5 * omega * omega * (V * V + A * A
                                                  + Q * Q + J * J)
    H = float((2.0 * np.pi) ** 2 * np.trapezoid(h_dens * r, r))

    # Chern-Simons charge: Q_CS = 2π ∫ r dr (A·∂_r J − J·∂_r A)
    qcs_dens = A * dJ - J * dA
    Q_CS = float(2.0 * np.pi * np.trapezoid(qcs_dens * r, r))

    # Alternative Q definition: L/Q-style Noether current
    # For e^{-iωt}: Q_noether ~ ω·(V²+A²+Q²+J²)
    q_noether_dens = omega * (V * V + A * A + Q * Q + J * J)
    Q_noether = float((2.0 * np.pi) ** 2 * np.trapezoid(q_noether_dens * r, r))

    HQ_CS = H / Q_CS if abs(Q_CS) > 1e-14 else float("nan")
    HQ_noether = H / Q_noether if abs(Q_noether) > 1e-14 else float("nan")

    return {
        "success": True,
        "V0": V0, "A0": A0, "Q0": Q0, "J0": J0,
        "omega": omega, "m_J_sq": m_J_sq, "m_eff_sq": m_eff_sq,
        "lam_bench": lam_bench,
        "include_omega_kinetic": include_omega_kinetic,
        "H": H, "Q_CS": Q_CS, "Q_noether": Q_noether,
        "HQ_CS": HQ_CS, "HQ_noether": HQ_noether,
        "tail": tail, "localized": localized,
        "blew_up": blew_up,
        "r_reached": float(r[-1]),
        "peak": float(peak),
        "r": r, "V": V, "A": A, "Q": Q, "J": J,
    }


# ---------------------------------------------------------------------------
# Calibration scan
# ---------------------------------------------------------------------------

def calibration_scan(omega=OMEGA_DEF,
                     V0=0.1, A0=0.1, Q0=0.1, J0=0.1,
                     m_J_sq_grid=None, lam_grid=None,
                     include_omega_kinetic=True,
                     r_max=R_MAX_DEF):
    """Find (m_J², λ_bench) giving H/Q ≈ 1.6969 at canonical point."""

    if m_J_sq_grid is None:
        # m_J² > ω² needed for decay; scan a range
        m_J_sq_grid = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
    if lam_grid is None:
        # Werbos's hint: λ_bench = 4·g = 4.25 (simplified); also scan around
        lam_grid = [0.5, 1.0, 2.0, 4.25, 10.0]

    print(f"\n{'='*92}")
    print(f"4-FN BENCHMARK CALIBRATION SCAN")
    print(f"{'='*92}")
    print(f"ω = {omega}, V₀=A₀=Q₀=J₀ = {V0}")
    print(f"Target H/Q (CS) = {TARGET_HQ}  (Werbos electron anchor, 0.56% off "
          f"{TARGET_HQ_PHYS} physical)")
    print(f"Target Q_CS = 1 (integer winding)")
    print(f"Target L/Q = ω = {omega}")
    print(f"include_omega_kinetic = {include_omega_kinetic}")
    print()
    print(f"{'m_J²':>6} {'λ_bench':>9} {'m_eff²':>8} {'H':>12} {'Q_CS':>10} "
          f"{'H/Q_CS':>10} {'tail':>9} {'loc':>5} {'peak':>9}")
    print("-" * 92)

    results = []
    for m_J_sq in m_J_sq_grid:
        for lam_bench in lam_grid:
            res = integrate(V0, A0, Q0, J0, omega, m_J_sq, lam_bench,
                            r_max=r_max,
                            include_omega_kinetic=include_omega_kinetic)
            if not res["success"]:
                m_eff_sq = m_J_sq - omega * omega
                print(f"{m_J_sq:>6.2f} {lam_bench:>9.3f} {m_eff_sq:>8.3f} "
                      f"  [FAIL: {res.get('error', '?')[:35]}]")
                continue
            results.append(res)
            print(f"{m_J_sq:>6.2f} {lam_bench:>9.3f} {res['m_eff_sq']:>8.3f} "
                  f"{res['H']:>12.4f} {res['Q_CS']:>10.4f} "
                  f"{res['HQ_CS']:>10.4f} {res['tail']:>9.4f} "
                  f"{'Y' if res['localized'] else 'n':>5} "
                  f"{res['peak']:>9.4f}")

    # Find closest match to H/Q = 1.6969 among localized + non-trivial
    valid = [r for r in results
             if r["localized"] and r["peak"] > 0.01
             and np.isfinite(r["HQ_CS"])]
    if valid:
        best = min(valid, key=lambda r: abs(r["HQ_CS"] - TARGET_HQ))
        print(f"\nBest match to H/Q = {TARGET_HQ}:")
        print(f"  m_J² = {best['m_J_sq']:.3f}, λ_bench = {best['lam_bench']:.3f}")
        print(f"  m_eff² = {best['m_eff_sq']:.3f}")
        print(f"  H/Q_CS = {best['HQ_CS']:.4f}  "
              f"(gap {abs(best['HQ_CS']-TARGET_HQ)/TARGET_HQ*100:.2f}%)")
        print(f"  Q_CS = {best['Q_CS']:.4f}  (target 1)")
        return best, results
    else:
        print(f"\nNo localized non-trivial solutions in this scan window.")
        return None, results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mode", choices=["calibration", "single"],
                        default="calibration")
    parser.add_argument("--omega", type=float, default=OMEGA_DEF)
    parser.add_argument("--V0", type=float, default=0.1)
    parser.add_argument("--A0", type=float, default=0.1)
    parser.add_argument("--Q0", type=float, default=0.1)
    parser.add_argument("--J0", type=float, default=0.1)
    parser.add_argument("--m-J-sq", type=float, default=1.0,
                        help="Single-mode m_J² value")
    parser.add_argument("--lam-bench", type=float, default=4.25,
                        help="Single-mode λ_bench value")
    parser.add_argument("--no-omega-kinetic", action="store_true",
                        help="Skip ω²/2·(V²+A²+Q²+J²) term in H "
                             "(test benchmark literal interpretation)")
    parser.add_argument("--r-max", type=float, default=R_MAX_DEF)
    args = parser.parse_args()

    include_omega_kinetic = not args.no_omega_kinetic

    t0 = time.time()

    if args.mode == "calibration":
        best, results = calibration_scan(
            omega=args.omega,
            V0=args.V0, A0=args.A0, Q0=args.Q0, J0=args.J0,
            include_omega_kinetic=include_omega_kinetic,
            r_max=args.r_max,
        )
    else:
        print(f"\nSingle-mode integration: m_J² = {args.m_J_sq}, "
              f"λ_bench = {args.lam_bench}, ω = {args.omega}")
        res = integrate(args.V0, args.A0, args.Q0, args.J0,
                        args.omega, args.m_J_sq, args.lam_bench,
                        r_max=args.r_max,
                        include_omega_kinetic=include_omega_kinetic)
        if res["success"]:
            print(f"  m_eff² = {res['m_eff_sq']:.4f}")
            print(f"  H = {res['H']:.4f}")
            print(f"  Q_CS = {res['Q_CS']:.4f}  (target 1)")
            print(f"  Q_noether = {res['Q_noether']:.4f}")
            print(f"  H/Q_CS = {res['HQ_CS']:.4f}  (target {TARGET_HQ})")
            print(f"  H/Q_noether = {res['HQ_noether']:.4f}")
            print(f"  tail = {res['tail']:.4f}  "
                  f"({'localized' if res['localized'] else 'NOT localized'})")
            print(f"  peak field amplitude = {res['peak']:.4f}")
        else:
            print(f"  FAILED: {res.get('error', 'unknown')}")
        results = [res]
        best = res if res.get("success") else None

    print(f"\nElapsed: {time.time() - t0:.1f}s")

    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / "m6_v4_4fn_results.json"
    payload = {
        "mode": args.mode,
        "omega": args.omega,
        "config": vars(args),
        "best": ({k: v for k, v in best.items()
                  if k not in ("r", "V", "A", "Q", "J")}
                 if best and best.get("success") else None),
        "n_scan_points": len(results) if isinstance(results, list) else 1,
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2,
                  default=lambda o: float(o) if hasattr(o, "item") else str(o))
    print(f"[json] {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
