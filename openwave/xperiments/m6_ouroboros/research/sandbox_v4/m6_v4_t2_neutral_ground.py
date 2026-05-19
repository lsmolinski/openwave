"""M6 sandbox v4 Track 2 — Neutral chaoiton ground state via B₀ shooting.

Refines v3's neutral-chaoiton scan to find the actual ground state
(minimum-H localized solution) of the Q=0 sector.

NEUTRAL ODE (α=0, ω=0 reduction of canonical Lean 2-fn):
    β'' + (1/r)β' − β/r² + λβ + 4gβ³ = 0
    BC at origin: β ~ B₀·r  (slope; β(0)=0)
    Localization: |β(r)| < 0.1 for r ≥ 8

v3 result (snapshot end of v3):
    23 localized solutions across (λ, B₀). At λ=1.0, B₀=0.010: H_code=0.491
    (gave m_χ=0.508 MeV with the v3 conversion factor).
    Claim: "smallest-B₀ branch gives the lightest state."

T2 questions:
    Q-T2-a: Is there a TRUE minimum of H(B₀) at some B₀_ground, with H
            growing for both B₀ > B₀_ground (more field, more energy) and
            B₀ < B₀_ground (insufficient nonlinear binding → solution
            de-localizes)?
    Q-T2-b: Or does H decrease monotonically as B₀ → 0 (asymptotic limit
            of trivial / nearly-zero solutions)?

Approach:
    1. Wide log-scan: B₀ ∈ [1e-5, 1.0], 60 points, at λ=1.0, g=1.0625.
    2. Plot H(B₀) and tail(B₀); identify minimum / localization edge.
    3. If a real minimum exists → bisect to refine B₀_ground.
    4. Also scan λ ∈ {0.1, 0.5, 1.0, 2.0, 5.0, 10.0} to map the ground-
       state surface H(λ, B₀_ground).

NOTE on mass conversion:
    v3 used H_CODE_ELECTRON_CALIB = 0.494 → m = H × 0.511 / 0.494. But
    T1 diagnostic showed the actual H at Werbos's electron-calibration
    point (g=1.0625, λ=1.0, ω=1.0, A₀=B₀=0.1) is ~10.7, not 0.494.
    So the mass conversion needs re-derivation once the 4-fn ODE is in
    hand. T2 reports raw H_code; mass conversion deferred.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


G_COUPLING = 1.0625
LAMBDA_DEF = 1.0
R_MIN = 0.01
R_MAX_DEF = 30.0     # wider than T1 — small B₀ → slower decay
N_GRID = 1500
RTOL = 1e-11
ATOL = 1e-13
LOC_R_CHECK = 8.0
LOC_THRESHOLD = 0.1


def ode_neutral(r, y, g, lam):
    """β'' + (1/r)β' − β/r² + λβ + 4gβ³ = 0
    State: y = [β, β']
    """
    b, db = y
    inv_r = 1.0 / r
    d2b = -inv_r * db + b * inv_r * inv_r - lam * b - 4.0 * g * b * b * b
    return [db, d2b]


def integrate(B0, g=G_COUPLING, lam=LAMBDA_DEF, r_max=R_MAX_DEF):
    y0 = [B0 * R_MIN, B0]
    r_eval = np.linspace(R_MIN, r_max, N_GRID)
    try:
        sol = solve_ivp(
            fun=lambda r, y: ode_neutral(r, y, g, lam),
            t_span=(R_MIN, r_max), y0=y0, t_eval=r_eval,
            method="RK45", rtol=RTOL, atol=ATOL,
        )
    except Exception:
        return None
    if not sol.success:
        return None

    r = sol.t
    b, db = sol.y
    if not np.all(np.isfinite(b)):
        return None

    mask = r >= LOC_R_CHECK
    if mask.sum() == 0:
        return None
    tail = float(np.abs(b[mask]).max())
    localized = bool(tail < LOC_THRESHOLD)

    # Energy: H = 4π ∫ r² [β'² + β²/r² + λβ²/2 + g β⁴] dr
    h_dens = db * db + b * b / (r * r) + lam * b * b / 2.0 + g * b ** 4
    H = float(4.0 * np.pi * np.trapezoid(h_dens * r * r, r))

    # Peak amplitude (proxy for "is this a nontrivial solution")
    b_max = float(np.abs(b).max())

    # Find the radial position of peak |β| (helps see if solution shape is sane)
    r_peak = float(r[np.argmax(np.abs(b))])

    # Count radial nodes (sign changes of β) — Lean stability spec: ≤4
    sign_changes = int(np.sum(np.diff(np.sign(b[b != 0])) != 0))

    return {
        "B0": B0, "lam": lam, "g": g,
        "H": H, "tail": tail, "localized": localized,
        "b_max": b_max, "r_peak": r_peak,
        "n_nodes": sign_changes,
        "regular": sign_changes <= 4,
    }


def scan_B0(lam, B0_min=1e-5, B0_max=1.0, n_points=60, r_max=R_MAX_DEF):
    """Log-scan of B₀ at fixed λ."""
    B0_grid = np.logspace(np.log10(B0_min), np.log10(B0_max), n_points)
    results = []
    for B0 in B0_grid:
        res = integrate(float(B0), lam=lam, r_max=r_max)
        if res is not None:
            results.append(res)
    return results


def bisect_min(lam, B0_lo, B0_hi, n_iters=20, r_max=R_MAX_DEF):
    """Golden-section search for minimum H within [B0_lo, B0_hi]
    among localized solutions.
    """
    phi = (np.sqrt(5) - 1) / 2  # ~0.618
    a, b = B0_lo, B0_hi
    # work in log space (since B₀ is log-spaced)
    log_a, log_b = np.log10(a), np.log10(b)

    def eval_H(log_B):
        r = integrate(10 ** log_B, lam=lam, r_max=r_max)
        if r is None or not r["localized"]:
            return float("inf"), r
        return r["H"], r

    c = log_b - phi * (log_b - log_a)
    d = log_a + phi * (log_b - log_a)
    H_c, res_c = eval_H(c)
    H_d, res_d = eval_H(d)

    history = [(10 ** c, H_c), (10 ** d, H_d)]
    for _ in range(n_iters):
        if H_c < H_d:
            log_b, d, H_d, res_d = d, c, H_c, res_c
            c = log_b - phi * (log_b - log_a)
            H_c, res_c = eval_H(c)
            history.append((10 ** c, H_c))
        else:
            log_a, c, H_c, res_c = c, d, H_d, res_d
            d = log_a + phi * (log_b - log_a)
            H_d, res_d = eval_H(d)
            history.append((10 ** d, H_d))

    # best of c, d
    if H_c < H_d:
        return 10 ** c, H_c, res_c, history
    return 10 ** d, H_d, res_d, history


def print_scan_table(results, lam):
    print(f"\nλ = {lam:.3f}")
    print(f"{'B0':>10} {'H_code':>10} {'b_max':>9} {'r_peak':>8} "
          f"{'tail':>8} {'loc':>5} {'nodes':>6}")
    print("-" * 65)
    last_was_loc = False
    for r in results:
        flag = "Y" if r["localized"] else "n"
        # Highlight localization edge
        marker = " ★" if r["localized"] != last_was_loc else "  "
        print(f"{r['B0']:>10.5f} {r['H']:>10.4f} {r['b_max']:>9.4f} "
              f"{r['r_peak']:>8.3f} {r['tail']:>8.4f} {flag:>5} "
              f"{r['n_nodes']:>6d}{marker}")
        last_was_loc = r["localized"]


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--lam", type=float, nargs="+",
                        default=[0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
    parser.add_argument("--B0-min", type=float, default=1e-5)
    parser.add_argument("--B0-max", type=float, default=1.0)
    parser.add_argument("--n-points", type=int, default=60)
    parser.add_argument("--r-max", type=float, default=R_MAX_DEF)
    parser.add_argument("--bisect", action="store_true",
                        help="After scan, golden-section search for H minimum")
    args = parser.parse_args()

    print(f"\n{'='*78}")
    print(f"M6 v4 T2 — Neutral Chaoiton Ground State Scan")
    print(f"{'='*78}")
    print(f"ODE: β'' + (1/r)β' − β/r² + λβ + 4gβ³ = 0  (α=0, ω=0 reduction)")
    print(f"g = {G_COUPLING}, r_max = {args.r_max}")
    print(f"B₀ ∈ [{args.B0_min}, {args.B0_max}], {args.n_points} log points")
    print(f"λ values: {args.lam}")

    all_results = {}
    t0 = time.time()
    for lam in args.lam:
        results = scan_B0(lam, args.B0_min, args.B0_max,
                          args.n_points, args.r_max)
        all_results[lam] = results
        print_scan_table(results, lam)

        # find minimum among localized
        loc = [r for r in results if r["localized"]]
        if loc:
            best = min(loc, key=lambda r: r["H"])
            print(f"  → min H_loc = {best['H']:.4f} at B₀={best['B0']:.5f} "
                  f"(b_max={best['b_max']:.3f}, nodes={best['n_nodes']})")
        else:
            print(f"  → no localized solutions in scan window")

    elapsed = time.time() - t0
    print(f"\nScan time: {elapsed:.1f}s")

    # Bisection refinement at λ=1.0 if requested
    if args.bisect:
        print(f"\n{'='*78}")
        print(f"GOLDEN-SECTION REFINEMENT around H minimum (λ=1.0)")
        print(f"{'='*78}")
        results_1 = all_results.get(1.0, [])
        loc_1 = [r for r in results_1 if r["localized"]]
        if len(loc_1) < 3:
            print("Not enough localized points at λ=1.0 to bisect.")
        else:
            # bracket around min
            Hs = [r["H"] for r in loc_1]
            B0s = [r["B0"] for r in loc_1]
            i_min = int(np.argmin(Hs))
            lo = B0s[max(0, i_min - 2)]
            hi = B0s[min(len(B0s) - 1, i_min + 2)]
            print(f"Bracket: B₀ ∈ [{lo:.5f}, {hi:.5f}]")
            B0_g, H_g, res_g, history = bisect_min(1.0, lo, hi, n_iters=18,
                                                    r_max=args.r_max)
            print(f"\nGround-state estimate (λ=1.0):")
            print(f"  B₀_ground = {B0_g:.6f}")
            print(f"  H_ground  = {H_g:.6f}")
            if res_g is not None:
                print(f"  tail      = {res_g['tail']:.5f}")
                print(f"  b_max     = {res_g['b_max']:.5f}")
                print(f"  nodes     = {res_g['n_nodes']}")
            print(f"\nBisection history (B₀, H):")
            for B0, H in history[-10:]:
                print(f"  {B0:.6f}  {H:.6f}")

    # save
    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / "m6_v4_t2_results.json"
    payload = {
        "config": {
            "g": G_COUPLING, "lam_values": list(args.lam),
            "B0_min": args.B0_min, "B0_max": args.B0_max,
            "n_points": args.n_points, "r_max": args.r_max,
        },
        "scans": {f"{lam:.3f}": all_results[lam] for lam in args.lam},
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2,
                  default=lambda o: float(o) if hasattr(o, "item") else str(o))
    print(f"\n[json] {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
