"""M6 sandbox v7 — Paul-script variant (2026-05-21 PM).

Tests interpretation (B) of Paul's 2026-05-21 reply: λ_LM is a CONSTANT
fixed at 1.0, NOT a free eigenvalue. This is what his "exact runnable
script" literally encodes (`lambda_lm = 1.0` at module scope), even though
that script doesn't run as-is (r=0 singularity + BC count mismatch).

Differences from v7 main script:

  - λ_LM is a constant (1.0), not a free `solve_bvp` parameter.
  - ONLY ω is the free eigenvalue.
  - NO V_norm anchor BC. BC count balances: 10 BCs = 9 states + 1 free param.
  - exp(-r) initial guess per Paul (V=+0.1*exp(-r), A=-0.1*exp(-r),
    Q=+0.1*exp(-r), J=-0.1*exp(-1.5r)). No warm-start option.
  - r grid starts at R_MIN > 0 to avoid the `J/r` divide-by-zero in
    Paul's literal script.
  - ODE structure preserved verbatim from Paul's script — including the
    differences from our v6/v7 form (Klein-Gordon mass −ω²V on every field;
    λ corrections in all 4 equations; quartic in V equation rather than Q).

If this lands Paul's prescribed ground state (V₀>0, Q₀>0, A₀<0, J₀<0,
≤4 nodes, H_electron/Q ≤ 1.6969), we've cracked the algorithmic puzzle.
If it diverges or lands a wrong basin, Paul's script-as-written is
inconsistent with his stated production behavior, and we email back
empirically rather than speculatively.

USAGE:
    python3 m6_v7_paul_script.py
    python3 m6_v7_paul_script.py --lambda-lm 0.5    # try other fixed λ values
    python3 m6_v7_paul_script.py --r-max 25         # larger domain
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp


# ---------------------------------------------------------------------------
# Constants — Paul's script values
# ---------------------------------------------------------------------------

OMEGA_INIT = 1.0
G_DEFAULT = 1.0625
LAMBDA_LM_FIXED = 1.0   # ★ Paul's script: lambda_lm is a constant

R_MIN = 0.05            # Patch: avoid Paul's r=0 J/r singularity
R_MAX_DEFAULT = 20.0    # Paul's script value
N_GRID_DEFAULT = 500    # Paul's script value

I_TARGET = 1.0

TARGET_HQ = 1.6969
TARGET_QCS = 1.0
PASS_HQ_TOL = 0.003
PASS_QCS_TOL = 0.01
PASS_MAX_NODES = 4
LOC_R_CHECK = 8.0


# ---------------------------------------------------------------------------
# ODE — Paul's script VERBATIM (with R_MIN > 0 patch only)
# ---------------------------------------------------------------------------

def make_rhs(omega_const, g_const, lambda_lm):
    """Paul's ODE from 2026-05-21 reply. Note: omega is fixed in the rhs
    closure; solve_bvp's free parameter only sets the BC k = ω in the
    Robin condition. (This matches Paul's literal script which has omega
    as a module-level constant.)

    Different from our v6/v7 ODE: this has Klein-Gordon −ω²V on every
    field equation, λ-corrections in all four equations (not just A/J),
    and the quartic appears in V's equation (not Q's).
    """

    def rhs(r, y, p):
        omega = p[0]  # ω is the only free param
        V, Vp, A, Ap, Q, Qp, J, Jp, I = y
        inv_r = 1.0 / r

        dV = Vp
        dVp = (-(omega**2)*V + Q
               - lambda_lm*(Jp + J*inv_r)
               - g_const*((V**2+Q**2)*V + (V*A - Q*J)*A))

        dA = Ap
        dAp = (-(omega**2)*A - J
               + lambda_lm*(-(Vp + V*inv_r))
               - g_const*((A**2+J**2)*A + (V*A - Q*J)*V))

        dQ = Qp
        dQp = (-(omega**2)*Q + V
               - lambda_lm*(Jp + J*inv_r)
               - g_const*((V**2+Q**2)*Q - (V*A - Q*J)*J))

        dJ = Jp
        dJp = (-(omega**2)*J - A
               + lambda_lm*(-(Vp + V*inv_r))
               - g_const*((A**2+J**2)*J - (V*A - Q*J)*Q))

        dI = r * (A * Jp - J * Ap)

        return np.vstack([dV, dVp, dA, dAp, dQ, dQp, dJ, dJp, dI])

    return rhs


def make_bc():
    """Paul's BCs verbatim. 10 residuals:
        ya: V'(R_MIN)=0, A'=0, Q'=0, J'=0, I(R_MIN)=0       (5)
        yb: V'+ωV=0, A'+ωA=0, Q'+ωQ=0, J'+ωJ=0, I(R_max)=1  (5)
    NO V_norm anchor.

    Total: 10 BCs = 9 states + 1 free param (ω).
    """

    def bc(ya, yb, p):
        omega = p[0]
        k = omega  # Paul's script uses k = omega directly

        V_a, Vp_a, A_a, Ap_a, Q_a, Qp_a, J_a, Jp_a, I_a = ya
        V_b, Vp_b, A_b, Ap_b, Q_b, Qp_b, J_b, Jp_b, I_b = yb

        residuals = np.array([
            Vp_a, Ap_a, Qp_a, Jp_a, I_a,
            Vp_b + k * V_b,
            Ap_b + k * A_b,
            Qp_b + k * Q_b,
            Jp_b + k * J_b,
            I_b - I_TARGET,
        ])
        return residuals

    return bc


def initial_guess(r_grid):
    """Paul's exp(-r) seed verbatim:
        V=+0.1*exp(-r), A=-0.1*exp(-r), Q=+0.1*exp(-r), J=-0.1*exp(-1.5r)
    plus derivatives consistent with those exponentials.
    """
    e = np.exp(-r_grid)
    e_J = np.exp(-1.5 * r_grid)

    V = +0.1 * e
    Vp = -0.1 * e        # d/dr exp(-r) = -exp(-r)
    A = -0.1 * e
    Ap = +0.1 * e
    Q = +0.1 * e
    Qp = -0.1 * e
    J = -0.1 * e_J
    Jp = +0.15 * e_J     # d/dr (-0.1*exp(-1.5r)) = +0.15*exp(-1.5r)

    I = (r_grid - r_grid[0]) / (r_grid[-1] - r_grid[0]) * I_TARGET
    return np.vstack([V, Vp, A, Ap, Q, Qp, J, Jp, I])


# ---------------------------------------------------------------------------
# Energy functionals — same as v7 for cross-comparison
# ---------------------------------------------------------------------------

def compute_H_full(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega, g):
    """Full H with DeepSeek quartic — canonical per Paul Q28."""
    kinetic = 0.5 * (Vp * Vp + Ap * Ap + Qp * Qp + Jp * Jp)
    omega_kin = 0.5 * omega * omega * (V * V + A * A + Q * Q + J * J)
    cross = -V * Q + A * J
    quartic = 0.25 * g * ((V * V + Q * Q) ** 2
                          + (A * A + J * J) ** 2
                          + 2.0 * (V * A - Q * J) ** 2)
    h_dens = kinetic + omega_kin + cross + quartic
    return float(np.trapezoid(h_dens, r))


def compute_H_electron(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega):
    """Drop quartic for electron H — Paul Q28."""
    kinetic = 0.5 * (Vp * Vp + Ap * Ap + Qp * Qp + Jp * Jp)
    omega_kin = 0.5 * omega * omega * (V * V + A * A + Q * Q + J * J)
    cross = -V * Q + A * J
    h_dens = kinetic + omega_kin + cross
    return float(np.trapezoid(h_dens, r))


def n_nodes_excluding_origin(f, r):
    mask = r > (r[0] + 0.1)
    f_active = f[mask]
    f_active = f_active[np.abs(f_active) > 1e-6]
    if len(f_active) < 2:
        return 0
    return int(np.sum(np.diff(np.sign(f_active)) != 0))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--omega-init", type=float, default=OMEGA_INIT)
    parser.add_argument("--lambda-lm", type=float, default=LAMBDA_LM_FIXED,
                        help="Fixed λ_LM constant (Paul's script: 1.0)")
    parser.add_argument("--g", type=float, default=G_DEFAULT)
    parser.add_argument("--r-min", type=float, default=R_MIN)
    parser.add_argument("--r-max", type=float, default=R_MAX_DEFAULT)
    parser.add_argument("--n-grid", type=int, default=N_GRID_DEFAULT)
    parser.add_argument("--max-nodes", type=int, default=5000,
                        help="Paul's script uses 5000")
    parser.add_argument("--tol", type=float, default=1e-6,
                        help="Paul's script uses 1e-6")
    parser.add_argument("--verbose", type=int, default=2)
    args = parser.parse_args()

    t0 = time.time()

    print(f"\n{'='*92}")
    print(f"M6 v7 — PAUL-SCRIPT VARIANT (fixed λ_LM={args.lambda_lm}, ω-only free)")
    print(f"{'='*92}")
    print(f"Fixed:  λ_LM = {args.lambda_lm} (CONSTANT, not free), g = {args.g}")
    print(f"Free:   ω only (init {args.omega_init})")
    print(f"Grid:   r ∈ [{args.r_min}, {args.r_max}], {args.n_grid} points, "
          f"max_nodes={args.max_nodes}, tol={args.tol}")
    print(f"Init:   Paul's exp(-r) seed (V=+0.1e^-r, A=-0.1e^-r, Q=+0.1e^-r, "
          f"J=-0.1e^-1.5r)")
    print(f"BCs:    10 residuals — derivatives=0 + I(R_min)=0 (5) + "
          f"Robin (4) + I(R_max)=1 (1). NO V_norm anchor.")
    print(f"ODE:    Paul's literal script (Klein-Gordon mass on all 4; "
          f"λ-corr in all 4; quartic in V equation)")
    print(f"Target: H_electron/Q = {TARGET_HQ}; nodes ≤ {PASS_MAX_NODES}; "
          f"Q(0)>0; A,J(0)<0")

    # Linear-spaced grid (Paul uses np.linspace, we shift start to R_MIN)
    r_grid = np.linspace(args.r_min, args.r_max, args.n_grid)

    y0 = initial_guess(r_grid)
    p0 = np.array([args.omega_init])

    rhs = make_rhs(args.omega_init, args.g, args.lambda_lm)
    bc = make_bc()

    print(f"\n  Solving BVP...")
    try:
        sol = solve_bvp(rhs, bc, r_grid, y0, p=p0,
                        tol=args.tol, max_nodes=args.max_nodes,
                        verbose=args.verbose)
    except Exception as e:
        print(f"  EXCEPTION: {str(e)[:200]}")
        return 1

    print(f"\n  solve_bvp status: {sol.status} ({sol.message})")
    print(f"  Final grid size: {len(sol.x)}")
    print(f"  Final ω = {sol.p[0]:.6f}")

    omega = float(sol.p[0])
    r = sol.x
    V, Vp, A, Ap, Q, Qp, J, Jp, I = sol.y

    H_full = compute_H_full(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega, args.g)
    H_electron = compute_H_electron(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega)

    Q_CS = float(I[-1])
    Q_CS_grid = float(np.trapezoid(r * (A * Jp - J * Ap), r))

    if abs(Q_CS) > 1e-14:
        HQ_full = H_full / Q_CS
        HQ_electron = H_electron / Q_CS
    else:
        HQ_full = HQ_electron = float("nan")

    mask = r >= LOC_R_CHECK
    if mask.sum() > 0:
        tail = float((np.abs(V[mask]) + np.abs(A[mask])
                      + np.abs(Q[mask]) + np.abs(J[mask])).max())
    else:
        tail = float(np.abs(V).max() + np.abs(A).max()
                     + np.abs(Q).max() + np.abs(J).max())

    nodes_V = n_nodes_excluding_origin(V, r)
    nodes_A = n_nodes_excluding_origin(A, r)
    nodes_Q = n_nodes_excluding_origin(Q, r)
    nodes_J = n_nodes_excluding_origin(J, r)
    max_nodes = max(nodes_V, nodes_A, nodes_Q, nodes_J)

    V0, A0, Q0, J0 = float(V[0]), float(A[0]), float(Q[0]), float(J[0])

    print(f"\n{'='*92}")
    print(f"CONVERGED SOLUTION DIAGNOSTICS")
    print(f"{'='*92}")
    print(f"  ω           = {omega:.6f}")
    print(f"  λ_LM        = {args.lambda_lm} (FIXED, not a solver output)")
    print(f"  V(R_MIN), A, Q, J = {V0:+.4f}, {A0:+.4f}, {Q0:+.4f}, {J0:+.4f}")
    print(f"  peak V/A/Q/J     = {float(np.abs(V).max()):.4f} / "
          f"{float(np.abs(A).max()):.4f} / "
          f"{float(np.abs(Q).max()):.4f} / {float(np.abs(J).max()):.4f}")
    print(f"  nodes V/A/Q/J (excl r=0) = {nodes_V} / {nodes_A} / "
          f"{nodes_Q} / {nodes_J}  (max {PASS_MAX_NODES})")
    print(f"  tail @r≥8        = {tail:.4f}")
    print(f"  Q_CS             = {Q_CS:.6f}  (target {TARGET_QCS})")
    print(f"  Q_CS from grid   = {Q_CS_grid:.6f}")
    print(f"")
    print(f"  H_full     (w/ quartic)  = {H_full:.6f}")
    print(f"  H_electron (no quartic)  = {H_electron:.6f}")
    print(f"  H_full/Q     = {HQ_full:.6f}")
    print(f"  H_electron/Q = {HQ_electron:.6f}  ★ PRIMARY (target {TARGET_HQ})")

    # Check pass criterion
    fail = []
    if abs(HQ_electron - TARGET_HQ) >= PASS_HQ_TOL:
        fail.append(f"|H_electron/Q - {TARGET_HQ}| = {abs(HQ_electron-TARGET_HQ):.4f} > {PASS_HQ_TOL}")
    if abs(Q_CS - TARGET_QCS) >= PASS_QCS_TOL:
        fail.append(f"|Q_CS - 1| = {abs(Q_CS-TARGET_QCS):.4f} > {PASS_QCS_TOL}")
    if max_nodes > PASS_MAX_NODES:
        fail.append(f"max_nodes = {max_nodes} > {PASS_MAX_NODES}")
    if Q0 <= 0:
        fail.append(f"Q(0) = {Q0:.4f} (should be > 0)")
    if A0 >= 0:
        fail.append(f"A(0) = {A0:.4f} (should be < 0)")
    if J0 >= 0:
        fail.append(f"J(0) = {J0:.4f} (should be < 0)")

    passed = len(fail) == 0
    print(f"\n  PASS: {'✅ YES' if passed else '❌ NO'}")
    for reason in fail:
        print(f"    - {reason}")

    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / "m6_v7_paul_script_results.json"
    payload = {
        "config": vars(args),
        "lambda_lm_fixed": args.lambda_lm,
        "solve_bvp": {
            "status": int(sol.status),
            "message": str(sol.message),
            "n_grid_final": int(len(sol.x)),
        },
        "diag": {
            "omega": omega, "lambda_LM_fixed": args.lambda_lm,
            "V0": V0, "A0": A0, "Q0": Q0, "J0": J0,
            "peak_V": float(np.abs(V).max()),
            "peak_A": float(np.abs(A).max()),
            "peak_Q": float(np.abs(Q).max()),
            "peak_J": float(np.abs(J).max()),
            "nodes_V": nodes_V, "nodes_A": nodes_A,
            "nodes_Q": nodes_Q, "nodes_J": nodes_J,
            "max_nodes": max_nodes,
            "tail": tail,
            "Q_CS": Q_CS, "Q_CS_grid": Q_CS_grid,
            "H_full": H_full, "H_electron": H_electron,
            "HQ_full": HQ_full, "HQ_electron": HQ_electron,
        },
        "passed": bool(passed),
        "fail_reasons": fail,
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
