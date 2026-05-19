"""M6 sandbox v4 T6→A — 4-function BVP solver for bound chaoiton.

After T6 IVP attempt showed the 4-fn ODE is an eigenvalue/shooting
problem (all initial conditions blew up across the scan grid), this
script uses scipy.integrate.solve_bvp with explicit decay BCs at r_max
to find the bound state directly.

ODE (same as m6_v4_4fn.py):
    V'' = -V'/r + Q
    A'' = -A'/r + J
    Q'' = -Q'/r + V + m_eff²·Q + λ·Q·(Q² − J²)
    J'' = -J'/r + A − m_eff²·J − λ·J·(Q² − J²)
    m_eff² = m_J² − ω²

BCs:
    At r=R_MIN:  V'=0, A'=0, Q'=0, J'=0    (regularity)
    At r=r_max:  V=0,  A=0,  Q=0,  J=0     (decay)

Two failure modes to guard against:
    1. Trivial solution V=A=Q=J=0 also satisfies BCs. Mitigate with
       non-trivial initial guess + check peak amplitude post-solve.
    2. solve_bvp may fail to converge. Mitigate with multiple initial
       guess shapes + relaxed tolerance.

USAGE:
    python3 m6_v4_4fn_bvp.py --m-J-sq 2.0 --lam-bench 1.0   # single point
    python3 m6_v4_4fn_bvp.py --mode scan                    # parameter scan
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp


G_V5 = 1.0625
LAMBDA_V5 = 1.0
OMEGA_DEF = 1.0

R_MIN = 0.05
R_MAX_DEF = 12.0    # shorter than IVP (faster BVP convergence)
N_MESH = 100        # initial mesh size (increased from 80)
BVP_TOL = 1e-3      # relaxed from 1e-4 (eigenvalue searches are sensitive)
BVP_MAX_NODES = 30000

TARGET_HQ = 1.6969
TARGET_HQ_PHYS = 1.6875
TARGET_QCS = 1.0


def ode_4fn_bvp(r, y, p, lam_bench, omega):
    """Vectorized ODE with ONE eigenvalue parameter: m_J².
    λ_bench is FIXED (per Werbos: 1.0 at calibration).

    State: y shape (8, m) — [V, V', A, A', Q, Q', J, J']
    p[0] = m_J² (unknown)
    """
    V, dV, A, dA, Q, dQ, J, dJ = y
    m_J_sq = p[0]
    m_eff_sq = m_J_sq - omega * omega
    inv_r = 1.0 / r
    QQ_JJ = Q * Q - J * J

    d2V = -inv_r * dV + Q
    d2A = -inv_r * dA + J
    d2Q = -inv_r * dQ + V + m_eff_sq * Q + lam_bench * Q * QQ_JJ
    d2J = -inv_r * dJ + A - m_eff_sq * J - lam_bench * J * QQ_JJ

    return np.vstack([dV, d2V, dA, d2A, dQ, d2Q, dJ, d2J])


def bc_4fn(ya, yb, p, V_norm=0.1):
    """Boundary conditions for 8-dim state + 1 unknown (m_J²):
       9 BCs total.

       At r=R_MIN: 4 regularity (V'=A'=Q'=J'=0)
                  + V(R_MIN) = V_norm (normalization, eats trivial)
       At r=r_max: 4 decay (V=A=Q=J=0)

       A(R_MIN) is NOT pinned — determined by BVP convergence. To bias
       toward a non-zero A solution, use Werbos's asymmetric initial
       guess (helicity='werbos': V=+, A=-, Q=+, J=-).
    """
    return np.array([
        ya[1], ya[3], ya[5], ya[7],     # 4 regularity BCs at origin
        ya[0] - V_norm,                  # V(R_MIN) = V_norm
        yb[0], yb[2], yb[4], yb[6],     # 4 decay BCs at infinity
    ])


def initial_guess(r, profile="exp_decay", amplitude=0.1, decay_rate=0.5,
                  phase_QJ=1.0, helicity="werbos"):
    """Construct initial guess for solve_bvp.

    profile: 'exp_decay' | 'sech' | 'gaussian' | 'oscillating'
    helicity:
      'symmetric'   → all 4 fields same sign (V=+, A=+, Q=+, J=+); Q_CS=0
      'werbos'      → V=+, A=-, Q=+, J=-  (V,Q same; A,J same; pair-pair
                       opposite); Werbos 2026-05-19 4:21 PM specified
                       this gives Q_CS=+1
      'werbos2'     → V=-, A=+, Q=-, J=+ ; Q_CS=-1 (same physics, sign
                       flip)
      'phase_QJ'    → legacy: V,A same sign; Q,J flipped by phase_QJ
                       (NOT Werbos's helicity; kept for back-compat)
    """
    n = len(r)
    y = np.zeros((8, n))
    a = amplitude

    if profile == "exp_decay":
        env = a * np.exp(-decay_rate * r)
        denv_dr = -decay_rate * env
    elif profile == "sech":
        env = a * np.cosh(r * decay_rate) ** -2  # sech²
        denv_dr = -2 * a * decay_rate * np.tanh(r * decay_rate) \
                   * np.cosh(r * decay_rate) ** -2
    elif profile == "gaussian":
        env = a * np.exp(-decay_rate * r * r)
        denv_dr = -2 * decay_rate * r * env
    elif profile == "oscillating":
        # Bound state may have a node — try one-node profile
        env = a * np.exp(-decay_rate * r) * (1 - r / 3.0)
        denv_dr = (-decay_rate * env - a * np.exp(-decay_rate * r) / 3.0)
    else:
        raise ValueError(f"Unknown profile {profile}")

    # Force derivative = 0 at origin (matches BC)
    denv_dr[0] = 0.0

    # Sign assignment per helicity
    if helicity == "symmetric":
        sV, sA, sQ, sJ = +1, +1, +1, +1
    elif helicity == "werbos":   # V=Q=+, A=J=-
        sV, sA, sQ, sJ = +1, -1, +1, -1
    elif helicity == "werbos2":  # V=Q=-, A=J=+
        sV, sA, sQ, sJ = -1, +1, -1, +1
    elif helicity == "phase_QJ":
        sV, sA = +1, +1
        sQ, sJ = phase_QJ, phase_QJ
    else:
        raise ValueError(f"Unknown helicity {helicity}")

    y[0] = sV * env;          y[1] = sV * denv_dr        # V
    y[2] = sA * env;          y[3] = sA * denv_dr        # A
    y[4] = sQ * env;          y[5] = sQ * denv_dr        # Q
    y[6] = sJ * env;          y[7] = sJ * denv_dr        # J
    return y


def try_bvp(m_J_sq_init, lam_bench, omega=OMEGA_DEF,
            r_max=R_MAX_DEF, n_mesh=N_MESH,
            profile="exp_decay", amplitude=0.1, decay_rate=0.5,
            phase_QJ=1.0, helicity="werbos", tol=BVP_TOL, verbose=False,
            V_norm=0.1):
    """BVP solve with ONE eigenvalue (m_J²). λ_bench FIXED.
    V(R_MIN)=V_norm is the only normalization.
    Initial guess uses asymmetric helicity to bias A,J non-zero.
    """
    r = np.linspace(R_MIN, r_max, n_mesh)
    y_init = initial_guess(r, profile=profile, amplitude=amplitude,
                            decay_rate=decay_rate, phase_QJ=phase_QJ,
                            helicity=helicity)
    # Override V(R_MIN) in initial guess to match norm
    y_init[0, 0] = V_norm
    p_init = np.array([m_J_sq_init])

    t0 = time.time()
    try:
        sol = solve_bvp(
            lambda r, y, p: ode_4fn_bvp(r, y, p, lam_bench, omega),
            lambda ya, yb, p: bc_4fn(ya, yb, p, V_norm),
            r, y_init, p=p_init, tol=tol, max_nodes=BVP_MAX_NODES,
            verbose=2 if verbose else 0,
        )
    except Exception as e:
        return {"success": False, "error": str(e)[:80],
                "m_J_sq_init": m_J_sq_init,
                "lam_bench": lam_bench,
                "helicity": helicity,
                "elapsed": time.time() - t0}

    elapsed = time.time() - t0
    m_J_sq_fit = float(sol.p[0]) if sol.p is not None else float("nan")
    lam_bench_fit = lam_bench    # was fixed input
    m_eff_sq = m_J_sq_fit - omega * omega

    if not sol.success:
        return {"success": False,
                "error": (sol.message or "BVP failed")[:80],
                "m_J_sq_init": m_J_sq_init,
                "lam_bench": lam_bench,
                "m_J_sq_fit": m_J_sq_fit,
                "helicity": helicity,
                "elapsed": elapsed}

    r_sol = sol.x
    V, dV, A, dA, Q, dQ, J, dJ = sol.y
    m_J_sq = m_J_sq_fit

    # Check for trivial-solution attractor
    peak = max(np.abs(V).max(), np.abs(A).max(),
               np.abs(Q).max(), np.abs(J).max())
    is_trivial = peak < 1e-4

    # Field values at origin (the "calibration" Werbos sets to 0.1)
    V0 = V[0]; A0 = A[0]; Q0 = Q[0]; J0 = J[0]

    # Energy functional (R=1)
    QQ_JJ = Q * Q - J * J
    h_dens = (dV * dV + dA * dA + dQ * dQ + dJ * dJ
              - V * Q + A * J
              + 0.5 * m_J_sq * QQ_JJ
              + 0.25 * lam_bench * QQ_JJ * QQ_JJ
              + 0.5 * omega * omega * (V * V + A * A + Q * Q + J * J))
    H = float((2.0 * np.pi) ** 2 * np.trapezoid(h_dens * r_sol, r_sol))

    # Charge: Q_CS = 2π ∫ r (A·∂_r J − J·∂_r A) dr
    qcs_dens = A * dJ - J * dA
    Q_CS = float(2.0 * np.pi * np.trapezoid(qcs_dens * r_sol, r_sol))

    HQ_CS = H / Q_CS if abs(Q_CS) > 1e-10 else float("nan")

    return {
        "success": True,
        "m_J_sq_init": m_J_sq_init,
        "m_J_sq_fit": m_J_sq,
        "lam_bench": lam_bench,
        "m_eff_sq": m_eff_sq,
        "omega": omega, "r_max": r_max,
        "profile": profile, "amplitude": amplitude,
        "decay_rate": decay_rate, "phase_QJ": phase_QJ,
        "helicity": helicity,
        "V_norm": V_norm,
        "H": H, "Q_CS": Q_CS, "HQ_CS": HQ_CS,
        "V0": float(V0), "A0": float(A0), "Q0": float(Q0), "J0": float(J0),
        "peak": float(peak),
        "is_trivial": is_trivial,
        "n_nodes": len(r_sol),
        "tol_achieved": sol.rms_residuals.max() if sol.rms_residuals is not None
                        else None,
        "elapsed": elapsed,
        "_r": r_sol, "_V": V, "_A": A, "_Q": Q, "_J": J,
    }


def print_result(res, idx=None):
    prefix = f"#{idx:>2}: " if idx is not None else ""
    if not res["success"]:
        m_init = res.get("m_J_sq_init", 0)
        l_init = res.get("lam_bench_init", 0)
        m_fit = res.get("m_J_sq_fit", float("nan"))
        l_fit = res.get("lam_bench_fit", float("nan"))
        print(f"{prefix}init(m_J²={m_init:.1f}, λ={l_init:.1f}) "
              f"fit(m_J²={m_fit:7.3f}, λ={l_fit:7.3f}) "
              f"[FAIL: {res.get('error', '?')[:30]}] "
              f"({res.get('elapsed', 0):.1f}s)")
        return
    flag = "TRIV" if res["is_trivial"] else "----"
    print(f"{prefix}init(m_J²={res['m_J_sq_init']:.1f}, λ={res['lam_bench_init']:.1f}) "
          f"fit(m_J²={res['m_J_sq_fit']:7.3f}, λ={res['lam_bench_fit']:7.3f}) "
          f"m_eff²={res['m_eff_sq']:6.3f} [{flag}] "
          f"peak={res['peak']:6.3f} H={res['H']:8.2f} "
          f"Q={res['Q_CS']:7.3f} H/Q={res['HQ_CS']:8.3f} "
          f"({res['elapsed']:.1f}s)")


def scan(omega=OMEGA_DEF, r_max=R_MAX_DEF):
    """Scan initial guesses for (m_J², λ_bench) eigenvalue pair."""
    m_J_sq_init_grid = [1.5, 2.5, 5.0]
    lam_init_grid = [0.5, 1.0, 4.25]
    profiles = [("exp_decay", 0.5, +1),
                ("exp_decay", 1.0, +1),
                ("exp_decay", 0.5, -1),
                ("sech", 0.5, +1)]

    print(f"\n{'='*100}")
    print(f"BVP CALIBRATION SCAN  (ω = {omega}, target H/Q ≈ {TARGET_HQ})")
    print(f"{'='*100}")
    print(f"Profiles: exp+, exp_steep+, exp-(phase), sech+, gauss+")
    print()

    results = []
    idx = 0
    for m_J_sq_init in m_J_sq_init_grid:
        for lam_init in lam_init_grid:
            for prof_name, decay, phase in profiles:
                idx += 1
                res = try_bvp(m_J_sq_init, lam_init, omega=omega,
                              r_max=r_max,
                              profile=prof_name, decay_rate=decay,
                              phase_QJ=phase)
                results.append(res)
                print_result(res, idx)

    # Find candidates: success + non-trivial
    valid = [r for r in results
             if r.get("success") and not r.get("is_trivial")]
    print(f"\n{'='*100}")
    print(f"NON-TRIVIAL BOUND STATES FOUND: {len(valid)} / {len(results)}")
    print(f"{'='*100}")

    if valid:
        print(f"\n{'m_J²_fit':>9} {'λ_fit':>8} {'profile':>10} {'peak':>9} "
              f"{'H':>11} {'Q_CS':>9} {'H/Q':>10} {'|H/Q-1.6969|':>14}")
        print("-" * 95)
        valid.sort(key=lambda r: abs(r.get("HQ_CS", 1e9) - TARGET_HQ)
                                 if np.isfinite(r.get("HQ_CS", np.nan))
                                 else 1e9)
        for r in valid[:10]:
            gap = (abs(r["HQ_CS"] - TARGET_HQ) if np.isfinite(r["HQ_CS"])
                   else float("inf"))
            print(f"{r['m_J_sq_fit']:>9.4f} {r['lam_bench_fit']:>8.4f} "
                  f"{r['profile']:>10} {r['peak']:>9.4f} "
                  f"{r['H']:>11.3f} {r['Q_CS']:>9.3f} "
                  f"{r['HQ_CS']:>10.3f} {gap:>14.3f}")

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mode", choices=["single", "scan"], default="single")
    parser.add_argument("--m-J-sq", type=float, default=0.5,
                        help="Initial guess for m_J² (Werbos: ~0.5)")
    parser.add_argument("--lam-bench", type=float, default=1.0,
                        help="Initial guess for λ_bench (Werbos: 1.0)")
    parser.add_argument("--V-norm", type=float, default=0.1,
                        help="V(R_MIN) normalization (Werbos: +0.1)")
    parser.add_argument("--A-norm", type=float, default=-0.1,
                        help="A(R_MIN) normalization (Werbos: -0.1 for "
                             "Q_CS=+1 helicity)")
    parser.add_argument("--helicity",
                        choices=["symmetric", "werbos", "werbos2", "phase_QJ"],
                        default="werbos",
                        help="Initial guess sign pattern (default werbos: "
                             "V=Q=+, A=J=-)")
    parser.add_argument("--omega", type=float, default=OMEGA_DEF)
    parser.add_argument("--r-max", type=float, default=R_MAX_DEF)
    parser.add_argument("--profile",
                        choices=["exp_decay", "sech", "gaussian", "oscillating"],
                        default="exp_decay")
    parser.add_argument("--amplitude", type=float, default=0.1)
    parser.add_argument("--decay-rate", type=float, default=0.5)
    parser.add_argument("--phase-QJ", type=float, default=1.0)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    t0 = time.time()

    if args.mode == "single":
        res = try_bvp(args.m_J_sq, args.lam_bench, omega=args.omega,
                      r_max=args.r_max, profile=args.profile,
                      amplitude=args.amplitude, decay_rate=args.decay_rate,
                      phase_QJ=args.phase_QJ, helicity=args.helicity,
                      verbose=args.verbose,
                      V_norm=args.V_norm)
        print(f"\nSingle BVP solve at:")
        print(f"  init: m_J²={args.m_J_sq}, λ_bench={args.lam_bench}, "
              f"ω={args.omega}, V_norm={args.V_norm}, A_norm={args.A_norm}")
        print(f"  profile={args.profile}, decay={args.decay_rate}, "
              f"phase={args.phase_QJ}")
        print()
        if res["success"]:
            print(f"  SUCCESS in {res['elapsed']:.1f}s "
                  f"({res['n_nodes']} mesh nodes)")
            print(f"  m_J²_fit (eigenvalue) = {res['m_J_sq_fit']:.4f}")
            print(f"  λ_bench (fixed) = {res['lam_bench']:.4f}")
            print(f"  m_eff² = {res['m_eff_sq']:.4f}")
            print(f"  peak field = {res['peak']:.6f}  "
                  f"({'TRIVIAL' if res['is_trivial'] else 'non-trivial'})")
            print(f"  V(0) = {res['V0']:.4f}, A(0) = {res['A0']:.4f}, "
                  f"Q(0) = {res['Q0']:.4f}, J(0) = {res['J0']:.4f}")
            print(f"  H = {res['H']:.4f}")
            print(f"  Q_CS = {res['Q_CS']:.4f}  (target 1)")
            print(f"  H/Q_CS = {res['HQ_CS']:.4f}  (target {TARGET_HQ})")
            if not res["is_trivial"] and np.isfinite(res["HQ_CS"]):
                gap = abs(res["HQ_CS"] - TARGET_HQ) / TARGET_HQ * 100
                print(f"  Gap from H/Q target: {gap:.2f}%")
        else:
            print(f"  FAIL in {res['elapsed']:.1f}s: {res.get('error', '?')}")
            if "m_J_sq_fit" in res and np.isfinite(res.get("m_J_sq_fit", np.nan)):
                print(f"  partial fit: m_J²={res['m_J_sq_fit']:.4f}, "
                      f"λ_bench={res.get('lam_bench', 'NA')}")
        results = [res]
    else:
        results = scan(omega=args.omega, r_max=args.r_max)

    print(f"\nTotal elapsed: {time.time() - t0:.1f}s")

    # Save (without large array data)
    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / "m6_v4_4fn_bvp_results.json"
    payload = {
        "mode": args.mode,
        "config": vars(args),
        "results": [{k: v for k, v in r.items() if not k.startswith("_")}
                    for r in results],
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2,
                  default=lambda o: float(o) if hasattr(o, "item") else str(o))
    print(f"[json] {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
