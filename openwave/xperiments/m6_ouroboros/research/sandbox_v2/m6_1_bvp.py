"""M6.1 BVP — Q=0 neutral chaoiton via boundary-value problem.

Sandbox v2 follow-up. The IVP attempt (m6_1_neutral_chaoiton.py) showed that
the locked-ansatz Q=0 problem cannot be solved as an IVP with regularity at
r=0: that BC selects Bessel J_0 (oscillatory) or I_0 (growing), never K_0
(decay). This script switches to BVP per benchmark §10's own recommendation
("Shooting method or collocation, bvp_solve in SciPy").

Problem statement:
    8 ODEs (4 fields × 2 derivatives) on r ∈ [r_min, r_max]
    8 BCs split between endpoints:
        At r_min:  V'=A'=Q'=J'=0   (regularity — 4 BCs)
        At r_max:  V=A=Q=J=0       (decay      — 4 BCs)

Benchmark ODE (§3 of theory/Numerical Benchmark...):
    Δ_r V  = Q                            with Δ_r f = f'' + (1/r) f'
    Δ_r A  = J
    Δ_r Q  = V + m_J² Q + λ Q (Q² − J²)
    Δ_r J  = A − m_J² J − λ J (Q² − J²)

Trivial solution V=A=Q=J≡0 satisfies all BCs. To find non-trivial solutions
we seed the initial guess with a strong locked-ansatz K_0-like envelope and
let scipy.solve_bvp converge. If only the trivial solution exists, this is
the negative result we report to Werbos.

Two sign conventions sweepable:
    locked_sign = +1   →   Q=+V, J=+A   (gives K_0 modified-Bessel decay
                                          in the linearized limit, MATCHES
                                          benchmark §4's predicted tail)
    locked_sign = −1   →   Q=−V, J=−A   (Werbos's stated locked ansatz in
                                          DarkMatterv1 §5; gives Bessel J_0
                                          oscillation in linearized limit)

Q.D from §"Open structural questions" — sign convention is a load-bearing
ambiguity. This script tests both empirically.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp

# ---------------------------------------------------------------------------
# Benchmark parameters
# ---------------------------------------------------------------------------

G_COUPLING = 1.0625
LAMBDA = 1.0
OMEGA = 1.0
R_MIN = 0.02
R_MAX_DEFAULT = 30.0
N_GRID_DEFAULT = 200       # solve_bvp adapts internally; coarse start
R_PHYS_FM = 191.0
HBAR_C_MEVFM = 197.32698
M_E_MEV = 0.51100
H_CODE_ELECTRON = 0.494    # from Calibration paper §5.1


# ---------------------------------------------------------------------------
# ODE — vectorized for solve_bvp
# ---------------------------------------------------------------------------

def ode_rhs(r, y, m_J_sq, lam, laplacian_dim=2):
    """Vectorized ODE: y is shape (8, N), returns dy/dr shape (8, N).

    State: y = [V, V', A, A', Q, Q', J, J']

    Equations (with Δ_r f = f'' + (D-1)/r · f' for D=2 or 3):
        V'' = Q − (D-1)/r · V'
        A'' = J − (D-1)/r · A'
        Q'' = V + m_J²·Q + λ·Q·(Q² − J²) − (D-1)/r · Q'
        J'' = A − m_J²·J − λ·J·(Q² − J²) − (D-1)/r · J'
    """
    V  = y[0]; Vp = y[1]
    A  = y[2]; Ap = y[3]
    Q  = y[4]; Qp = y[5]
    J  = y[6]; Jp = y[7]

    inv_r_factor = (laplacian_dim - 1) / r
    s = Q * Q - J * J

    Vpp = Q                         - inv_r_factor * Vp
    App = J                         - inv_r_factor * Ap
    Qpp = V + m_J_sq * Q + lam * Q * s - inv_r_factor * Qp
    Jpp = A - m_J_sq * J - lam * J * s - inv_r_factor * Jp

    return np.vstack([Vp, Vpp, Ap, App, Qp, Qpp, Jp, Jpp])


def bc_residuals(ya, yb):
    """Boundary conditions: 4 at r_min (regularity), 4 at r_max (decay)."""
    return np.array([
        ya[1],   # V'(r_min) = 0
        ya[3],   # A'(r_min) = 0
        ya[5],   # Q'(r_min) = 0
        ya[7],   # J'(r_min) = 0
        yb[0],   # V(r_max) = 0
        yb[2],   # A(r_max) = 0
        yb[4],   # Q(r_max) = 0
        yb[6],   # J(r_max) = 0
    ])


# ---------------------------------------------------------------------------
# Initial guess — locked-ansatz K_0-like exponential decay
# ---------------------------------------------------------------------------

def initial_guess(r, V0, A0, decay_length, locked_sign):
    """y_init shape (8, N).

    Locked ansatz at seed:  Q = locked_sign · V, J = locked_sign · A.
    Envelope: V(r) = V0 · exp(-(r-r_min)/L),  same for A.
    Derivatives consistent with the envelope.
    """
    r_shift = r - r[0]
    env = np.exp(-r_shift / decay_length)
    d_env = -env / decay_length

    V   = V0 * env
    Vp  = V0 * d_env
    A   = A0 * env
    Ap  = A0 * d_env
    Q   = locked_sign * V
    Qp  = locked_sign * Vp
    J   = locked_sign * A
    Jp  = locked_sign * Ap

    return np.vstack([V, Vp, A, Ap, Q, Qp, J, Jp])


# ---------------------------------------------------------------------------
# Solve one BVP instance
# ---------------------------------------------------------------------------

def solve_one(m_J_sq, lam=LAMBDA, V0=0.1, A0=0.1,
              decay_length=5.0, locked_sign=+1,
              r_min=R_MIN, r_max=R_MAX_DEFAULT, n_grid=N_GRID_DEFAULT,
              laplacian_dim=2, max_nodes=20000, tol=1e-6, verbose=False):
    """Run scipy.solve_bvp; return result dict."""
    r = np.linspace(r_min, r_max, n_grid)
    y_init = initial_guess(r, V0=V0, A0=A0,
                            decay_length=decay_length, locked_sign=locked_sign)

    try:
        sol = solve_bvp(
            fun=lambda r, y: ode_rhs(r, y, m_J_sq, lam, laplacian_dim),
            bc=bc_residuals,
            x=r, y=y_init,
            tol=tol, max_nodes=max_nodes,
            verbose=2 if verbose else 0,
        )
    except Exception as e:
        return {"success": False, "error": str(e)[:120]}

    return {
        "success": bool(sol.success),
        "status": int(sol.status),
        "message": str(sol.message),
        "n_iter": int(getattr(sol, "niter", -1)),
        "rms_residual": float(sol.rms_residuals.max()) if hasattr(sol, "rms_residuals") and sol.rms_residuals is not None else None,
        "sol": sol,
    }


def post_process(sol_obj, m_J_sq, lam, locked_sign, n_eval=2000):
    """Sample solution, compute energy + Q_CS + diagnostics."""
    if sol_obj is None or not sol_obj.success:
        return None

    sol = sol_obj
    r_lo, r_hi = sol.x[0], sol.x[-1]
    r = np.linspace(r_lo, r_hi, n_eval)
    y = sol.sol(r)

    V  = y[0]; Vp = y[1]
    A  = y[2]; Ap = y[3]
    Q  = y[4]; Qp = y[5]
    J  = y[6]; Jp = y[7]

    # Energy density (benchmark §5, prefactor (2π)²R set to 1 here)
    s = Q * Q - J * J
    density = (Vp * Vp + Ap * Ap + Qp * Qp + Jp * Jp
               - V * Q + A * J
               + 0.5 * m_J_sq * s
               + 0.25 * lam * s * s)
    H_code = float((2.0 * np.pi) ** 2 * np.trapezoid(r * density, r))

    # Chern-Simons charge
    Q_CS = float(2.0 * np.pi * np.trapezoid(r * (A * Jp - J * Ap), r))

    # Locked-ansatz drift
    drift_QV = float(np.max(np.abs(Q - locked_sign * V)))
    drift_JA = float(np.max(np.abs(J - locked_sign * A)))

    # Max amplitudes (used to detect trivial solution)
    max_V = float(np.max(np.abs(V)))
    max_A = float(np.max(np.abs(A)))

    is_trivial = max_V < 1e-6 and max_A < 1e-6

    # m_χ via electron calibration (only meaningful if non-trivial)
    m_chi_mev = H_code * M_E_MEV / H_CODE_ELECTRON

    return {
        "r": r, "V": V, "A": A, "Q": Q, "J": J,
        "Vp": Vp, "Ap": Ap, "Qp": Qp, "Jp": Jp,
        "H_code": H_code,
        "m_chi_MeV": m_chi_mev,
        "Q_CS": Q_CS,
        "max_V": max_V, "max_A": max_A,
        "is_trivial": is_trivial,
        "drift_QV": drift_QV, "drift_JA": drift_JA,
    }


# ---------------------------------------------------------------------------
# Sweep harness
# ---------------------------------------------------------------------------

def run_sweep():
    """Try multiple (sign, m_J², decay_length, amplitude) combinations."""
    rows = []

    # Sweep grid — note: V0 = A0 gives V²-A²=0 → cubic vanishes →
    # trivial solution is unique. We MUST break asymmetry to find
    # non-trivial localized states.
    SIGNS = (+1, -1)
    M_J_SQ_VALUES = (0.0, 0.1, 0.5, 1.0, G_COUPLING)
    DECAY_LENGTHS = (3.0, 5.0, 10.0)
    # (V0, A0) pairs — INCLUDE asymmetric pairs so (V²-A²) ≠ 0
    AMPL_PAIRS = (
        (0.1, 0.1),    # symmetric (trivial-only)
        (0.1, 0.3),    # mild asymmetry, A > V
        (0.3, 0.1),    # mild asymmetry, V > A
        (0.1, 1.0),    # strong asymmetry, A >> V
        (1.0, 0.1),    # strong asymmetry, V >> A
        (0.5, 0.5),    # bigger symmetric amplitude
    )

    print(f"\n{'='*112}")
    print(f"M6.1 BVP sweep — Q=0 neutral chaoiton via scipy.solve_bvp")
    print(f"{'='*112}\n")
    print(f"{'sign':>4} {'m_J²':>7} {'L_init':>7} {'V0':>5} {'A0':>5} "
          f"{'OK':>3} {'triv':>4} {'max|V|':>9} {'H':>10} {'mMeV':>10} {'Q_CS':>10} {'iters':>5}")
    print("-" * 112)

    for sign in SIGNS:
        for m_J_sq in M_J_SQ_VALUES:
            for L in DECAY_LENGTHS:
                for V0, A0 in AMPL_PAIRS:
                    res = solve_one(m_J_sq, locked_sign=sign,
                                     V0=V0, A0=A0,
                                     decay_length=L)
                    if not res["success"]:
                        row = {
                            "sign": sign, "m_J_sq": m_J_sq, "L_init": L,
                            "V0": V0, "A0": A0,
                            "ok": False, "err": res.get("error") or res.get("message", "?")[:40],
                        }
                        rows.append(row)
                        print(f"{sign:>+4} {m_J_sq:7.3f} {L:7.1f} {V0:5.2f} {A0:5.2f} "
                              f"{'NO':>3} (FAIL: {row['err'][:50]})")
                        continue

                    pp = post_process(res["sol"], m_J_sq, LAMBDA, sign)
                    if pp is None:
                        continue

                    row = {
                        "sign": sign, "m_J_sq": m_J_sq, "L_init": L,
                        "V0": V0, "A0": A0,
                        "ok": True,
                        "n_iter": res["n_iter"],
                        "rms_res": res["rms_residual"],
                        "is_trivial": pp["is_trivial"],
                        "max_V": pp["max_V"], "max_A": pp["max_A"],
                        "H_code": pp["H_code"],
                        "m_chi_MeV": pp["m_chi_MeV"],
                        "Q_CS": pp["Q_CS"],
                        "drift_QV": pp["drift_QV"], "drift_JA": pp["drift_JA"],
                    }
                    rows.append(row)
                    triv = "Y" if pp["is_trivial"] else "n"
                    print(f"{sign:>+4} {m_J_sq:7.3f} {L:7.1f} {V0:5.2f} {A0:5.2f} "
                          f"{'Y':>3} {triv:>4} {pp['max_V']:9.4f} "
                          f"{pp['H_code']:10.4f} {pp['m_chi_MeV']:10.4f} "
                          f"{pp['Q_CS']:10.3e} {res['n_iter']:5d}")

    print(f"{'='*108}\n")

    # Surface non-trivial converged solutions
    non_triv = [r for r in rows if r.get("ok") and not r.get("is_trivial", True)]
    print(f"Non-trivial converged solutions: {len(non_triv)} / {len(rows)} rows")
    if non_triv:
        print("\nBest (lowest H_code among non-trivial):")
        best = sorted(non_triv, key=lambda r: abs(r.get("H_code", float("inf"))))
        for r in best[:5]:
            print(f"  sign={r['sign']:+d}, m_J²={r['m_J_sq']:.3f}, "
                  f"L_init={r['L_init']:.1f}, V0={r['V0']:.2f}  →  "
                  f"H={r['H_code']:.4f}, m_χ={r['m_chi_MeV']:.4f} MeV, "
                  f"Q_CS={r['Q_CS']:.2e}")
    else:
        print("\n>>> NO non-trivial localized solution found across "
              f"{len(rows)} parameter points.")
        print("    This is decision-grade: the locked-ansatz Q=0 chaoiton")
        print("    does not exist as a BVP solution of the benchmark ODE.")

    # Save JSON
    out = Path(__file__).resolve().parent.parent / "plots" / "m6_1_bvp_sweep.json"
    out.parent.mkdir(exist_ok=True)
    with open(out, "w") as f:
        json.dump(rows, f, indent=2, default=str)
    print(f"\n[json] {out}")

    return rows, non_triv


# ---------------------------------------------------------------------------
# Single-run mode + plot
# ---------------------------------------------------------------------------

def run_single(args):
    A0 = args.A0 if args.A0 is not None else args.V0
    print(f"\nM6.1 BVP — single run")
    print(f"  sign={args.locked_sign:+d}, m_J²={args.m_j_sq}, λ={LAMBDA}")
    print(f"  V0={args.V0}, A0={A0}, L_init={args.decay_length}")
    print(f"  r ∈ [{R_MIN}, {args.r_max}], n_grid={args.n_grid}")

    res = solve_one(
        m_J_sq=args.m_j_sq,
        lam=LAMBDA, V0=args.V0, A0=A0,
        decay_length=args.decay_length,
        locked_sign=args.locked_sign,
        r_min=R_MIN, r_max=args.r_max, n_grid=args.n_grid,
        verbose=args.verbose,
    )

    if not res["success"]:
        print(f"\n[FAIL] {res.get('error') or res.get('message')}")
        return 1

    pp = post_process(res["sol"], args.m_j_sq, LAMBDA, args.locked_sign)
    print(f"\n  iterations          = {res['n_iter']}")
    print(f"  max RMS residual    = {res['rms_residual']}")
    print(f"  max|V|              = {pp['max_V']:.6f}")
    print(f"  max|A|              = {pp['max_A']:.6f}")
    print(f"  is_trivial          = {pp['is_trivial']}")
    print(f"  H_code              = {pp['H_code']:.6f}")
    print(f"  m_χ (MeV)           = {pp['m_chi_MeV']:.6f}")
    print(f"  Q_CS                = {pp['Q_CS']:.6e}")
    print(f"  locked-drift |Q-sV| = {pp['drift_QV']:.6e}")
    print(f"  locked-drift |J-sA| = {pp['drift_JA']:.6e}")

    # Plot
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        ax[0].plot(pp["r"], pp["V"], label="V")
        ax[0].plot(pp["r"], pp["A"], label="A")
        ax[0].plot(pp["r"], pp["Q"], "--", label="Q")
        ax[0].plot(pp["r"], pp["J"], "--", label="J")
        ax[0].legend(); ax[0].set_xlabel("r"); ax[0].set_title("BVP fields")
        s = pp["Q"]**2 - pp["J"]**2
        density = (pp["Vp"]**2 + pp["Ap"]**2 + pp["Qp"]**2 + pp["Jp"]**2
                   - pp["V"]*pp["Q"] + pp["A"]*pp["J"]
                   + 0.5*args.m_j_sq*s + 0.25*LAMBDA*s*s)
        ax[1].plot(pp["r"], pp["r"]*density)
        ax[1].set_xlabel("r"); ax[1].set_ylabel("r·energy density")
        ax[1].set_title(f"H = {pp['H_code']:.4f}")
        plt.tight_layout()
        out_dir = Path(__file__).resolve().parent.parent / "plots"
        out_dir.mkdir(exist_ok=True)
        plot_path = out_dir / "m6_1_bvp_profile.png"
        plt.savefig(plot_path, dpi=120)
        plt.close()
        print(f"\n  [plot] {plot_path}")
    except ImportError:
        pass

    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="mode", required=False)

    p_single = sub.add_parser("single", help="single BVP run")
    p_single.add_argument("--locked-sign", type=int, default=+1, choices=[+1, -1])
    p_single.add_argument("--m-j-sq", type=float, default=0.0)
    p_single.add_argument("--V0", type=float, default=0.1)
    p_single.add_argument("--A0", type=float, default=None,
                          help="if omitted, defaults to V0 (symmetric)")
    p_single.add_argument("--decay-length", type=float, default=5.0)
    p_single.add_argument("--r-max", type=float, default=R_MAX_DEFAULT)
    p_single.add_argument("--n-grid", type=int, default=N_GRID_DEFAULT)
    p_single.add_argument("--verbose", action="store_true")

    p_sweep = sub.add_parser("sweep", help="parameter sweep")

    args = parser.parse_args()

    if args.mode == "single":
        return run_single(args)
    else:
        run_sweep()
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
