"""M6 sandbox v4 Track 1 — Lepton mass spectrum scan (2026-05-19).

Tests Werbos's hypothesis (DeepSeek reply 2026-05-18 19:39): the three
leptons (electron, muon, tau) are the lowest 3 stable modes in the Q=1
sector of the canonical 2-function Lean ODE.

CANONICAL ODE (chaoiton_theorem.lean.txt §ode_system):
    α'' + (1/r)α' − (1/r²)α + ω²α = β               (eq A)
    β'' + (1/r)β' − (1/r²)β + ω²β = α − λβ − 4gβ³   (eq B)

Fixed:   g = 1.0625, λ = 1.0  (electron-calibration coupling)
Scan:    ω in [0.5, 50] step 0.5 by default
         At each ω: grid over (A₀, B₀) origin slopes
Output:  per-ω lowest-H localized solution; 3 lowest H → candidate leptons

Expected from Werbos's bosonization session (May18 4PM) at same g, λ:
    ω = 1.00  → m = 0.511 MeV  (electron, calibration)
    ω = 12.78 → m = 106.8 MeV  (muon target 105.66, gap 1.1%)
    ω = 15.00 → m = 145.1 MeV  (pion+ target 139.57, gap 3.9%)
    ω ≈ 40-41 → tau predicted (target 1776.86, no prior numerical)

PASS criterion (G1 + empirical G3): 3 lowest-H ω give masses within
4-6% of electron / muon / tau.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


# ---------------------------------------------------------------------------
# Constants (from Lean theorem + Calibration paper)
# ---------------------------------------------------------------------------

G_COUPLING = 1.0625
LAMBDA_DEF = 1.0
R_MIN = 0.01
R_MAX_DEF = 20.0
N_GRID = 1000
RTOL = 1e-10
ATOL = 1e-12
LOC_R_CHECK = 8.0
LOC_THRESHOLD = 0.1

M_E_MEV = 0.51100
H_CODE_ELECTRON_CALIB = 0.494  # Calibration paper code-units electron H

M_MUON_MEV = 105.66
M_TAU_MEV = 1776.86
M_PION_PLUS_MEV = 139.57


# ---------------------------------------------------------------------------
# Canonical 2-function vector-Laplacian ODE
# ---------------------------------------------------------------------------

def ode_lean(r, y, g, lam, omega):
    a, da, b, db = y
    inv_r = 1.0 / r
    inv_r2 = inv_r * inv_r
    w2 = omega * omega
    d2a = b - w2 * a - inv_r * da + a * inv_r2
    d2b = a - lam * b - 4.0 * g * b * b * b - w2 * b - inv_r * db + b * inv_r2
    return [da, d2a, db, d2b]


def integrate(A0, B0, omega, g=G_COUPLING, lam=LAMBDA_DEF, r_max=R_MAX_DEF):
    y0 = [A0 * R_MIN, A0, B0 * R_MIN, B0]
    r_eval = np.linspace(R_MIN, r_max, N_GRID)
    try:
        sol = solve_ivp(
            fun=lambda r, y: ode_lean(r, y, g, lam, omega),
            t_span=(R_MIN, r_max), y0=y0, t_eval=r_eval,
            method="RK45", rtol=RTOL, atol=ATOL,
        )
    except Exception:
        return None
    if not sol.success:
        return None

    r = sol.t
    a, da, b, db = sol.y

    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return None

    mask = r >= LOC_R_CHECK
    if mask.sum() == 0:
        return None
    tail = float((np.abs(a[mask]) + np.abs(b[mask])).max())
    localized = bool(tail < LOC_THRESHOLD)

    h_dens = (da * da + a * a / (r * r)
              + db * db + b * b / (r * r)
              - a * b + lam * b * b / 2.0 + g * b ** 4)
    H = float(4.0 * np.pi * np.trapezoid(h_dens * r * r, r))

    q_dens = omega * (a * a + b * b)
    Q = float(4.0 * np.pi * np.trapezoid(q_dens * r * r, r))

    return {
        "A0": A0, "B0": B0, "omega": omega,
        "H": H, "Q": Q,
        "HQ": float(H / Q) if Q > 1e-14 else float("nan"),
        "tail": tail, "localized": localized,
        "m_MeV": float(H * M_E_MEV / H_CODE_ELECTRON_CALIB),
    }


# ---------------------------------------------------------------------------
# Scan
# ---------------------------------------------------------------------------

def run_scan(omegas, A0_grid, B0_grid, lam=LAMBDA_DEF, r_max=R_MAX_DEF, verbose=True):
    per_omega_best = {}
    per_omega_all_localized = {}
    t0 = time.time()
    n_total = len(omegas) * len(A0_grid) * len(B0_grid)
    done = 0

    if verbose:
        print(f"\n{'ω':>6} {'best A0':>9} {'best B0':>9} {'H_code':>10} "
              f"{'H/Q':>9} {'m (MeV)':>11} {'tail':>8} {'#loc':>5}")
        print("-" * 78)

    for omega in omegas:
        best = None
        loc_list = []
        for A0 in A0_grid:
            for B0 in B0_grid:
                res = integrate(A0, B0, float(omega), lam=lam, r_max=r_max)
                done += 1
                if res is None:
                    continue
                if res["localized"]:
                    loc_list.append(res)
                    if best is None or res["H"] < best["H"]:
                        best = res
        if best is not None:
            per_omega_best[float(omega)] = best
            per_omega_all_localized[float(omega)] = loc_list
            if verbose:
                print(f"{omega:>6.2f} {best['A0']:>9.3f} {best['B0']:>9.3f} "
                      f"{best['H']:>10.4f} {best['HQ']:>9.4f} "
                      f"{best['m_MeV']:>11.3f} {best['tail']:>8.4f} {len(loc_list):>5d}")
        else:
            if verbose and (omega * 2) % 4 == 0:  # report every 2nd missing ω
                print(f"{omega:>6.2f}    [no localized]")

    elapsed = time.time() - t0
    if verbose:
        print(f"\nScan complete: {done}/{n_total} integrations in {elapsed:.1f}s")
        print(f"ω values with localized solutions: {len(per_omega_best)}/{len(omegas)}")

    return per_omega_best, per_omega_all_localized


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def summarize_lepton_match(per_omega_best, top_k=5):
    if not per_omega_best:
        print("\nNo localized solutions found — cannot test lepton hypothesis.")
        return []

    sorted_best = sorted(per_omega_best.values(), key=lambda r: r["H"])

    print(f"\n{'='*78}")
    print(f"TOP {top_k} LOWEST-H LOCALIZED MODES  (candidate lepton spectrum)")
    print(f"{'='*78}")
    print(f"{'rank':>4} {'ω':>7} {'H_code':>10} {'m (MeV)':>11} {'A0':>6} {'B0':>6}")
    print("-" * 58)
    for i, r in enumerate(sorted_best[:top_k]):
        print(f"{i+1:>4} {r['omega']:>7.2f} {r['H']:>10.4f} "
              f"{r['m_MeV']:>11.3f} {r['A0']:>6.3f} {r['B0']:>6.3f}")

    targets = [
        ("electron", 1.00, M_E_MEV),
        ("muon", 12.78, M_MUON_MEV),
        ("pion+", 15.00, M_PION_PLUS_MEV),
        ("tau (est.)", 40.7, M_TAU_MEV),
    ]
    print(f"\n{'='*78}")
    print("WERBOS BOSONIZATION-SESSION TARGET COMPARISON")
    print(f"{'='*78}")
    print(f"{'particle':>11} {'ω_target':>10} {'ω_found':>10} "
          f"{'m_found':>11} {'m_target':>11} {'gap %':>8}")
    print("-" * 70)
    ws = list(per_omega_best.keys())
    matches = []
    for name, omega_t, mass_t in targets:
        nearest = min(ws, key=lambda w: abs(w - omega_t))
        r = per_omega_best[nearest]
        gap = abs(r["m_MeV"] - mass_t) / mass_t * 100.0
        matches.append({"particle": name, "omega_target": omega_t,
                        "omega_found": nearest, "m_found": r["m_MeV"],
                        "m_target": mass_t, "gap_pct": gap})
        print(f"{name:>11} {omega_t:>10.2f} {nearest:>10.2f} "
              f"{r['m_MeV']:>11.3f} {mass_t:>11.3f} {gap:>8.2f}")
    return matches, sorted_best


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--omega-min", type=float, default=0.5)
    parser.add_argument("--omega-max", type=float, default=50.0)
    parser.add_argument("--omega-step", type=float, default=0.5)
    parser.add_argument("--lam", type=float, default=LAMBDA_DEF)
    parser.add_argument("--r-max", type=float, default=R_MAX_DEF)
    parser.add_argument("--A0", type=float, nargs="+",
                        default=[0.05, 0.1, 0.2, 0.5, 1.0, 2.0])
    parser.add_argument("--B0", type=float, nargs="+",
                        default=[0.05, 0.1, 0.2, 0.5, 1.0, 2.0])
    parser.add_argument("--quick", action="store_true",
                        help="Coarser ω scan + smaller (A0,B0) grid for fast feedback")
    args = parser.parse_args()

    if args.quick:
        omegas = np.array([1.0, 5.0, 10.0, 12.78, 15.0, 20.0, 25.0, 30.0, 35.0, 40.7, 45.0])
        A0_grid = [0.05, 0.1, 0.5]
        B0_grid = [0.05, 0.1, 0.5]
    else:
        omegas = np.arange(args.omega_min,
                           args.omega_max + args.omega_step / 2,
                           args.omega_step)
        A0_grid = args.A0
        B0_grid = args.B0

    print(f"\n{'='*78}")
    print(f"M6 v4 T1 — Lepton Spectrum Scan  (canonical Lean ODE)")
    print(f"{'='*78}")
    print(f"g = {G_COUPLING}, λ = {args.lam}")
    print(f"ω: {len(omegas)} values in [{omegas[0]:.2f}, {omegas[-1]:.2f}]")
    print(f"A0 grid: {A0_grid}")
    print(f"B0 grid: {B0_grid}")
    print(f"Total integrations: {len(omegas) * len(A0_grid) * len(B0_grid)}")

    per_omega_best, per_omega_all = run_scan(
        omegas, A0_grid, B0_grid, lam=args.lam, r_max=args.r_max)

    if per_omega_best:
        matches, sorted_best = summarize_lepton_match(per_omega_best)
    else:
        matches, sorted_best = [], []

    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / "m6_v4_lepton_scan_results.json"
    payload = {
        "config": {
            "g": G_COUPLING, "lam": args.lam,
            "omega_min": float(omegas[0]), "omega_max": float(omegas[-1]),
            "omega_step": args.omega_step,
            "A0_grid": list(A0_grid), "B0_grid": list(B0_grid),
            "r_max": args.r_max,
        },
        "per_omega_best": {f"{k:.2f}": v for k, v in per_omega_best.items()},
        "lowest_5_modes": sorted_best[:5],
        "lepton_targets": matches,
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2, default=lambda o: float(o) if hasattr(o, "item") else str(o))
    print(f"\n[json] {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
