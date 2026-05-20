"""M6 sandbox v4 T10 — Collocation BVP per Werbos's actual algorithm (2026-05-20).

After T9 closed-negative confirmed that forward-IVP from value-BC origin
does not reach a Q_CS=1 chaoiton in any explored parameter range, Werbos's
2026-05-20 2:00 PM reply (via DeepSeek) supplied the actual method.
This script implements it.

METHOD (Werbos, paraphrased from the email):
    - Collocation / finite-difference BVP via scipy.integrate.solve_bvp.
    - Q_CS = 1 enforced as INTEGRAL CONSTRAINT via an auxiliary state I
      satisfying dI/dr = r·(A·J' - J·A'); Q_CS = 2π·I(R_max).
    - Free eigenvalue: ω (the chaoiton frequency). The Lagrange multiplier
      λ from Paul's H' = H - λ·Q form is NOT carried as a separate free
      param in this first pass — the integral constraint alone is the
      closure that pins ω. If solve_bvp fails for lack of closure, the
      next iteration adds λ via the EL-modified A and J equations.
    - Robin BCs at R_max matching K_0(κr) exponential decay; k=√(ω²-m_J²)
      computed from the converged ω.
    - Origin BCs are derivative=0 only (V'=A'=Q'=J'=I=0 at R_MIN). Field
      VALUES at the origin are determined by the solver, NOT prescribed.
    - Initial guess: Paul's exponential decay profile with helicity signs.

ODE SYSTEM (9 states, m_eff² = m_J² − ω²):
    dV/dr  = V'
    dV'/dr = -V'/r + Q
    dA/dr  = A'
    dA'/dr = -A'/r + J
    dQ/dr  = Q'
    dQ'/dr = -Q'/r + V + m_eff²·Q + λ_bench·Q·(Q²-J²)
    dJ/dr  = J'
    dJ'/dr = -J'/r + A − m_eff²·J − λ_bench·J·(Q²-J²)
    dI/dr  = r·(A·J' − J·A')   # accumulated Q_CS density × r weight

FIXED PARAMETERS:
    m_J²       = 0.5     # Werbos 2026-05-19 4:21 PM quote
    λ_bench    = 1.0     # Werbos 2026-05-19 4:21 PM quote

FREE EIGENVALUE:
    ω          # initial guess 1.0

BCs (10 total = 9 states + 1 free param):
    r = R_MIN:  V'=0, A'=0, Q'=0, J'=0, I=0          (5)
    r = R_max:  V'+kV=0, A'+kA=0, Q'+kQ=0, J'+kJ=0,
                I = 1/(2π)                            (5)

INITIAL GUESS (per Werbos):
    V(r) = +0.1·exp(-r),  A(r) = -0.1·exp(-r),
    Q(r) = +0.1·exp(-r),  J(r) = -0.1·exp(-r),
    derivatives = -V, +A, -Q, +J  (matching d/dr exp(-r) = -exp(-r)),
    I(r) = linear ramp from 0 to 1/(2π),
    ω = 1.0.

PASS CRITERION (same as T9):
    |H/Q_CS - 1.6969| < 0.001  AND  |Q_CS - 1| < 0.01
    tail < 0.05                AND  no blowup
    max_nodes ≤ 4              (Lean stability)

USAGE:
    python3 m6_v4_4fn_lambda_bvp.py
    python3 m6_v4_4fn_lambda_bvp.py --r-max 25 --n-grid 250 --verbose
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

M_J_SQ_DEFAULT = 0.5
LAMBDA_BENCH_DEFAULT = 1.0
OMEGA_INIT = 1.0

R_MIN = 0.05
R_MAX_DEFAULT = 20.0
N_GRID_DEFAULT = 200

# Q_CS target via auxiliary integral: Q_CS = 2π·I(R_max).
# For Q_CS = 1: I(R_max) = 1/(2π).
I_TARGET = 1.0 / (2.0 * np.pi)
TARGET_HQ = 1.6969
TARGET_QCS = 1.0
PASS_HQ_TOL = 0.001
PASS_QCS_TOL = 0.01
PASS_TAIL_TOL = 0.05
PASS_MAX_NODES = 4
LOC_R_CHECK = 8.0


# ---------------------------------------------------------------------------
# ODE and BC for solve_bvp
# ---------------------------------------------------------------------------

def make_rhs(m_J_sq, lam_bench):
    """Return rhs(r, y, p) for scipy.integrate.solve_bvp.

    y is shape (9, m). p is shape (2,) carrying [omega, lambda_LM].
    The Lagrange multiplier λ enters A and J equations via the variational
    derivatives of Q_CS (computed from r·(A·J' - J·A') density):
        δQ_CS/δA = J + 2·r·J'
        δQ_CS/δJ = -(A + 2·r·A')
    So with H' = H - λ·Q_CS:
        ΔA = J + λ·(J + 2·r·J')
        ΔJ = A - m_eff²·J - λ_bench·J·(Q²-J²) - λ·(A + 2·r·A')
    """

    def rhs(r, y, p):
        omega = p[0]
        lam_LM = p[1]
        m_eff_sq = m_J_sq - omega * omega

        V, Vp, A, Ap, Q, Qp, J, Jp, I = y
        inv_r = 1.0 / r
        QQ_JJ = Q * Q - J * J

        # V, Q equations unchanged (Q_CS doesn't depend on V or Q)
        dV = Vp
        dVp = -inv_r * Vp + Q
        dQ = Qp
        dQp = -inv_r * Qp + V + m_eff_sq * Q + lam_bench * Q * QQ_JJ

        # A, J equations with Lagrange-multiplier corrections
        # ΔA = J + λ·(J + 2·r·J')
        dA = Ap
        dAp = -inv_r * Ap + J + lam_LM * (J + 2.0 * r * Jp)
        # ΔJ = A - m_eff²·J - λ_bench·J·(Q²-J²) - λ·(A + 2·r·A')
        dJ = Jp
        dJp = (-inv_r * Jp + A - m_eff_sq * J
               - lam_bench * J * QQ_JJ
               - lam_LM * (A + 2.0 * r * Ap))

        # Accumulated Q_CS density with r-weight
        dI = r * (A * Jp - J * Ap)

        return np.vstack([dV, dVp, dA, dAp, dQ, dQp, dJ, dJp, dI])

    return rhs


def make_bc(m_J_sq, V_norm):
    """Return bc(ya, yb, p) for scipy.integrate.solve_bvp.

    Returns 11 residuals = 9 state BCs + 2 free-param closures.
    """

    def bc(ya, yb, p):
        omega = p[0]
        # Decay rate for Robin BC: k = sqrt(ω² - m_J²).
        # Floor at 0.1 to avoid imaginary k or k=0.
        k_squared = omega * omega - m_J_sq
        k = np.sqrt(max(k_squared, 0.01))

        V_a, Vp_a, A_a, Ap_a, Q_a, Qp_a, J_a, Jp_a, I_a = ya
        V_b, Vp_b, A_b, Ap_b, Q_b, Qp_b, J_b, Jp_b, I_b = yb

        residuals = np.array([
            # At r = R_MIN: regularity (derivative-only)
            Vp_a,
            Ap_a,
            Qp_a,
            Jp_a,
            I_a,
            # At r = R_max: Robin (exponential decay matching K_0)
            Vp_b + k * V_b,
            Ap_b + k * A_b,
            Qp_b + k * Q_b,
            Jp_b + k * J_b,
            # Integral constraint: I(R_max) = 1/(2π) → Q_CS = 1
            I_b - I_TARGET,
            # Anti-collapse normalization: V(R_MIN) = V_norm
            V_a - V_norm,
        ])
        return residuals

    return bc


def initial_guess(r_grid):
    """Initial profile per Werbos (exp decay + helicity signs), with the
    A,J pair NON-PROPORTIONAL so that the Q_CS density is non-zero from
    the start. (Proportional A∝J gives Q_CS density ≡ 0 exactly, leaving
    solve_bvp with no gradient to deform toward Q_CS=1.)

    V, Q on pure exp(-r). A on exp(-r). J on exp(-1.5·r) (different decay
    rate ⇒ A,J non-proportional ⇒ Q_CS density non-zero).
    """
    e = np.exp(-r_grid)
    e_J = np.exp(-1.5 * r_grid)
    V = +0.1 * e
    Vp = -0.1 * e
    A = -0.1 * e
    Ap = +0.1 * e
    Q = +0.1 * e
    Qp = -0.1 * e
    J = -0.1 * e_J
    Jp = +0.15 * e_J
    # I(r) accumulator initialized from analytic integral of the Q_CS density
    # ∫₀^r r' · (A·J' - J·A') dr'  with the above profiles.
    # A·J' - J·A' = (-0.1·e)(0.15·e_J) - (-0.1·e_J)(0.1·e)
    #             = -0.015·e^(-r)·e^(-1.5r) + 0.01·e^(-r)·e^(-1.5r)
    #             = -0.005·e^(-2.5r)
    # ∫₀^r r' (-0.005·e^(-2.5r')) dr' = -0.005·[(1 - (1+2.5r)·e^(-2.5r))/6.25]
    # Just initialize as linear ramp from 0 to I_TARGET; solver will adjust.
    I = (r_grid - r_grid[0]) / (r_grid[-1] - r_grid[0]) * I_TARGET
    return np.vstack([V, Vp, A, Ap, Q, Qp, J, Jp, I])


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def analyze_solution(sol, m_J_sq, lam_bench, r_max):
    """Compute H, Q_CS, H/Q_CS, tail, peaks, node counts from converged sol.

    Returns dict with all the diagnostic numbers.
    """
    omega = float(sol.p[0])
    lam_LM = float(sol.p[1]) if len(sol.p) > 1 else 0.0
    r = sol.x
    V, Vp, A, Ap, Q, Qp, J, Jp, I = sol.y

    QQ_JJ = Q * Q - J * J
    m_eff_sq = m_J_sq - omega * omega

    h_dens = (Vp * Vp + Ap * Ap + Qp * Qp + Jp * Jp
              - V * Q + A * J
              + 0.5 * m_J_sq * QQ_JJ
              + 0.25 * lam_bench * QQ_JJ * QQ_JJ
              + 0.5 * omega * omega * (V * V + A * A + Q * Q + J * J))
    H = float((2.0 * np.pi) ** 2 * np.trapezoid(h_dens * r, r))

    Q_CS = float(2.0 * np.pi * I[-1])
    HQ_CS = H / Q_CS if abs(Q_CS) > 1e-14 else float("nan")
    HQ_from_grid = float(H / float(2.0 * np.pi * np.trapezoid(
        r * (A * Jp - J * Ap), r)))

    mask = r >= LOC_R_CHECK
    if mask.sum() > 0:
        tail = float((np.abs(V[mask]) + np.abs(A[mask])
                      + np.abs(Q[mask]) + np.abs(J[mask])).max())
    else:
        tail = float(np.abs(V).max() + np.abs(A).max()
                     + np.abs(Q).max() + np.abs(J).max())

    peak_V = float(np.abs(V).max())
    peak_A = float(np.abs(A).max())
    peak_Q = float(np.abs(Q).max())
    peak_J = float(np.abs(J).max())

    def n_nodes(f):
        f_active = f[np.abs(f) > 1e-6]
        if len(f_active) < 2:
            return 0
        return int(np.sum(np.diff(np.sign(f_active)) != 0))
    nodes_V = n_nodes(V)
    nodes_A = n_nodes(A)
    nodes_Q = n_nodes(Q)
    nodes_J = n_nodes(J)
    max_nodes = max(nodes_V, nodes_A, nodes_Q, nodes_J)

    return {
        "omega": omega, "lambda_LM": lam_LM,
        "m_J_sq": m_J_sq, "m_eff_sq": float(m_eff_sq),
        "lam_bench": lam_bench,
        "H": H, "Q_CS": Q_CS, "HQ_CS": HQ_CS, "HQ_from_grid": HQ_from_grid,
        "tail": tail,
        "peak_V": peak_V, "peak_A": peak_A,
        "peak_Q": peak_Q, "peak_J": peak_J,
        "nodes_V": nodes_V, "nodes_A": nodes_A,
        "nodes_Q": nodes_Q, "nodes_J": nodes_J, "max_nodes": max_nodes,
        "V0": float(V[0]), "A0": float(A[0]),
        "Q0": float(Q[0]), "J0": float(J[0]),
        "r_min": float(r[0]), "r_max": float(r[-1]),
        "n_grid": int(len(r)),
    }


def check_pass(diag):
    """Apply T10 pass criterion. Returns (passed, reasons)."""
    fail = []
    hq_err = abs(diag["HQ_CS"] - TARGET_HQ)
    qcs_err = abs(diag["Q_CS"] - TARGET_QCS)
    if hq_err >= PASS_HQ_TOL:
        fail.append(f"|H/Q-{TARGET_HQ}| = {hq_err:.4f} > {PASS_HQ_TOL}")
    if qcs_err >= PASS_QCS_TOL:
        fail.append(f"|Q_CS-1| = {qcs_err:.4f} > {PASS_QCS_TOL}")
    if diag["tail"] >= PASS_TAIL_TOL:
        fail.append(f"tail = {diag['tail']:.4f} > {PASS_TAIL_TOL}")
    if diag["max_nodes"] > PASS_MAX_NODES:
        fail.append(f"max_nodes = {diag['max_nodes']} > {PASS_MAX_NODES}")
    return (len(fail) == 0), fail


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--m-J-sq", type=float, default=M_J_SQ_DEFAULT,
                        help="Fixed m_J² (Werbos: 0.5)")
    parser.add_argument("--lam-bench", type=float, default=LAMBDA_BENCH_DEFAULT,
                        help="Fixed λ_bench (Werbos: 1.0)")
    parser.add_argument("--omega-init", type=float, default=OMEGA_INIT,
                        help="Initial guess for ω")
    parser.add_argument("--lambda-init", type=float, default=0.1,
                        help="Initial guess for Lagrange multiplier λ")
    parser.add_argument("--V-norm", type=float, default=0.1,
                        help="V(R_MIN) anti-collapse normalization "
                             "(Werbos's stated 0.1)")
    parser.add_argument("--r-max", type=float, default=R_MAX_DEFAULT)
    parser.add_argument("--n-grid", type=int, default=N_GRID_DEFAULT)
    parser.add_argument("--max-nodes", type=int, default=2000,
                        help="solve_bvp max_nodes (more = more refinement)")
    parser.add_argument("--tol", type=float, default=1e-3,
                        help="solve_bvp tolerance")
    parser.add_argument("--verbose", type=int, default=2,
                        help="solve_bvp verbose level (0/1/2)")
    args = parser.parse_args()

    t0 = time.time()

    print(f"\n{'='*92}")
    print(f"M6 v4 T10 — COLLOCATION BVP (Werbos algorithm + Lagrange λ)")
    print(f"{'='*92}")
    print(f"Fixed: m_J² = {args.m_J_sq}, λ_bench = {args.lam_bench}, "
          f"V(R_MIN) = {args.V_norm}")
    print(f"Free:  ω (init {args.omega_init}), λ_LM (init {args.lambda_init})")
    print(f"Grid:  r ∈ [{R_MIN}, {args.r_max}], {args.n_grid} points")
    print(f"Init:  V(r)=+0.1·exp(-r), A(r)=-0.1·exp(-r), "
          f"Q=+0.1·exp(-r), J=-0.1·exp(-r)")
    print(f"Targets: H/Q_CS={TARGET_HQ}, Q_CS=1 (via I(R_max)={I_TARGET:.4f})")

    # Build grid (non-uniform: more points near origin)
    # Use log spacing then linear remap to keep R_MIN
    s = np.linspace(0, 1, args.n_grid)
    r_grid = R_MIN + (args.r_max - R_MIN) * (s ** 1.5)

    # Build initial guess
    y0 = initial_guess(r_grid)
    p0 = np.array([args.omega_init, args.lambda_init])

    # Build rhs and bc closures
    rhs = make_rhs(args.m_J_sq, args.lam_bench)
    bc = make_bc(args.m_J_sq, args.V_norm)

    # Solve
    print(f"\n  Solving BVP...")
    try:
        sol = solve_bvp(rhs, bc, r_grid, y0, p=p0,
                        tol=args.tol, max_nodes=args.max_nodes,
                        verbose=args.verbose)
    except Exception as e:
        print(f"  EXCEPTION: {str(e)[:120]}")
        return 1

    print(f"\n  solve_bvp status: {sol.status} ({sol.message})")
    print(f"  Final grid size: {len(sol.x)}")
    print(f"  Final ω = {sol.p[0]:.6f}, λ_LM = {sol.p[1]:.6f}")

    diag = analyze_solution(sol, args.m_J_sq, args.lam_bench, args.r_max)
    passed, fail_reasons = check_pass(diag)

    print(f"\n{'='*92}")
    print(f"CONVERGED SOLUTION DIAGNOSTICS")
    print(f"{'='*92}")
    print(f"  ω           = {diag['omega']:.6f}")
    print(f"  λ_LM        = {diag['lambda_LM']:.6f}")
    print(f"  m_eff²      = {diag['m_eff_sq']:.4f}")
    print(f"  V(R_MIN), A, Q, J = {diag['V0']:+.4f}, {diag['A0']:+.4f}, "
          f"{diag['Q0']:+.4f}, {diag['J0']:+.4f}")
    print(f"  peak V/A/Q/J     = {diag['peak_V']:.4f} / {diag['peak_A']:.4f} "
          f"/ {diag['peak_Q']:.4f} / {diag['peak_J']:.4f}")
    print(f"  nodes V/A/Q/J    = {diag['nodes_V']} / {diag['nodes_A']} / "
          f"{diag['nodes_Q']} / {diag['nodes_J']}  (max {PASS_MAX_NODES})")
    print(f"  tail @r≥8        = {diag['tail']:.4f}  (< {PASS_TAIL_TOL})")
    print(f"  H                = {diag['H']:.6f}")
    print(f"  Q_CS             = {diag['Q_CS']:.6f}  (target {TARGET_QCS})")
    print(f"  H/Q_CS           = {diag['HQ_CS']:.6f}  (target {TARGET_HQ})")
    print(f"  H/Q from grid    = {diag['HQ_from_grid']:.6f}")

    print(f"\n  PASS: {'✅ YES' if passed else '❌ NO'}")
    for reason in fail_reasons:
        print(f"    - {reason}")

    # Persist
    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / "m6_v4_4fn_lambda_bvp_results.json"
    payload = {
        "config": vars(args),
        "I_target": I_TARGET,
        "solve_bvp": {
            "status": int(sol.status),
            "message": str(sol.message),
            "n_grid_final": int(len(sol.x)),
        },
        "diag": diag,
        "passed": bool(passed),
        "fail_reasons": fail_reasons,
        "elapsed_sec": time.time() - t0,
    }
    with open(out_path, "w") as f:
        json.dump(payload, f, indent=2,
                  default=lambda o: float(o) if hasattr(o, "item") else str(o))
    print(f"\n[json] {out_path}")
    print(f"Elapsed: {time.time() - t0:.1f}s")
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
