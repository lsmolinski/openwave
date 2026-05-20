"""M6 sandbox v6 — Calibrated BVP with DeepSeek normalizations (2026-05-20).

After v5 attempt 4 converged at Q_CS=1 with H/Q=52.64 (vs target 1.6969, 31×
off), Paul + DeepSeek replied 2026-05-20 ~4:00 PM with answers to Q22/Q23/Q24:

    Q22 — Q_CS normalization: drop the 2π prefactor.
          Q_CS = ∫_0^∞ r·(A·J' − J·A') dr     (no extra factor)
          Set I_TARGET = 1.0 (not 1/(2π)).

    Q23 — H functional kinetic: use (1/2)(V')² (not (V')²).
          Drop the (2π)²·R toroidal prefactor on H.
          Try DeepSeek's quartic structure:
              (g/4)·((V²+Q²)² + (A²+J²)² + 2·(V·A − Q·J)²)
          Falls back to the v5 Numerical-Benchmark form
          (λ_bench/4)·(Q² − J²)² if --quartic=benchmark is passed.

    Q24 — Reference profile: warm-start from DeepSeek's 8-point table
          interpolated onto our grid. Selects ground-state basin.
          (DeepSeek's profile may be fabricated rather than from an actual
          run; we use it as a warm-start regardless — better than naive
          exp(-r) for selecting the ≤4-node ground state.)

Predicted H/Q post-normalization-fix (back-of-envelope from v5's 52.64):
    52.64 × 0.5 (kinetic) × 1/(2π)² (drop H prefactor) × 2π (Q_CS factor)
    ≈ 4.2  (vs target 1.6969)

If we land near 1.6969 → calibration closed. If still 2-3× off → quartic
discrepancy or λ-coefficient discrepancy is the gap; request DeepSeek's
Python script as definitive cross-check.

ODE SYSTEM (9 states, m_eff² = m_J² − ω²):
    dV/dr  = V'
    dV'/dr = -V'/r + Q
    dA/dr  = A'
    dA'/dr = -A'/r + J + λ_LM·(J + 2·r·J')
    dQ/dr  = Q'
    dQ'/dr = -Q'/r + V + m_eff²·Q + λ_bench·Q·(Q²-J²)
    dJ/dr  = J'
    dJ'/dr = -J'/r + A − m_eff²·J − λ_bench·J·(Q²-J²) − λ_LM·(A + 2·r·A')
    dI/dr  = r·(A·J' − J·A')        # raw Q_CS density (no 2π)

FIXED PARAMETERS:
    m_J²       = 0.5     # Werbos 2026-05-19 4:21 PM
    λ_bench    = 1.0     # Werbos 2026-05-19 4:21 PM

FREE EIGENVALUES:
    ω, λ_LM             # initial guess (1.0, 0.1) — same as v5 attempt 4

BCs (11 total = 9 states + 2 free params):
    r = R_MIN:  V'=0, A'=0, Q'=0, J'=0, I=0          (5)
    r = R_max:  V'+kV=0, A'+kA=0, Q'+kQ=0, J'+kJ=0,
                I = 1.0                                (5)   ★ NEW: I_TARGET=1
    Anti-collapse: V(R_MIN) = V_norm                  (1)

PASS CRITERION:
    |H/Q_CS - 1.6969| < 0.001  AND  |Q_CS - 1| < 0.01
    tail < 0.05                AND  max_nodes ≤ 4

USAGE:
    python3 m6_v6_4fn_calibrated_bvp.py
    python3 m6_v6_4fn_calibrated_bvp.py --quartic deepseek  # try deepseek form
    python3 m6_v6_4fn_calibrated_bvp.py --quartic benchmark # fall back
    python3 m6_v6_4fn_calibrated_bvp.py --no-warm-start    # skip ref profile
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp


# ---------------------------------------------------------------------------
# Constants — v6 normalization conventions per DeepSeek's reply
# ---------------------------------------------------------------------------

M_J_SQ_DEFAULT = 0.5
LAMBDA_BENCH_DEFAULT = 1.0
OMEGA_INIT = 1.0
LAMBDA_LM_INIT = 0.1

R_MIN = 0.05
R_MAX_DEFAULT = 20.0
N_GRID_DEFAULT = 200

# DeepSeek Q22: Q_CS = ∫r·(A·J'-J·A')dr directly, no 2π factor.
# So I_TARGET = 1.0 (not 1/(2π) as in v5).
I_TARGET = 1.0

TARGET_HQ = 1.6969
TARGET_QCS = 1.0
PASS_HQ_TOL = 0.001
PASS_QCS_TOL = 0.01
PASS_TAIL_TOL = 0.05
PASS_MAX_NODES = 4
LOC_R_CHECK = 8.0


# DeepSeek Q24: 8-point reference profile for the electron chaoiton.
# Format: (r, V, A, Q, J).
DEEPSEEK_REF_PROFILE = np.array([
    [0.0,  +0.1000, -0.1000, +0.1000, -0.1000],
    [0.2,  +0.0892, -0.0903, +0.0881, -0.0894],
    [0.5,  +0.0670, -0.0691, +0.0650, -0.0670],
    [1.0,  +0.0381, -0.0402, +0.0355, -0.0377],
    [2.0,  +0.0098, -0.0109, +0.0081, -0.0092],
    [3.0,  +0.0023, -0.0027, +0.0017, -0.0021],
    [5.0,  +0.0001, -0.0002, +0.0000, -0.0001],
    [8.0,  +0.0000,  0.0000,  0.0000,  0.0000],
])


# ---------------------------------------------------------------------------
# ODE and BC for solve_bvp
# ---------------------------------------------------------------------------

def make_rhs(m_J_sq, lam_bench):
    """Return rhs(r, y, p) for scipy.integrate.solve_bvp.
    Identical ODE structure to v5 — DeepSeek didn't change the equations,
    only the integral/energy normalizations.
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
        dA = Ap
        dAp = -inv_r * Ap + J + lam_LM * (J + 2.0 * r * Jp)
        dJ = Jp
        dJp = (-inv_r * Jp + A - m_eff_sq * J
               - lam_bench * J * QQ_JJ
               - lam_LM * (A + 2.0 * r * Ap))

        # Accumulated Q_CS density (raw, no 2π factor per DeepSeek Q22)
        dI = r * (A * Jp - J * Ap)

        return np.vstack([dV, dVp, dA, dAp, dQ, dQp, dJ, dJp, dI])

    return rhs


def make_bc(m_J_sq, V_norm):
    """Return bc(ya, yb, p) for scipy.integrate.solve_bvp.

    Returns 11 residuals = 9 state BCs + 2 free-param closures.
    Same structure as v5; only I_TARGET changes (1.0 instead of 1/(2π)).
    """

    def bc(ya, yb, p):
        omega = p[0]
        k_squared = omega * omega - m_J_sq
        k = np.sqrt(max(k_squared, 0.01))

        V_a, Vp_a, A_a, Ap_a, Q_a, Qp_a, J_a, Jp_a, I_a = ya
        V_b, Vp_b, A_b, Ap_b, Q_b, Qp_b, J_b, Jp_b, I_b = yb

        residuals = np.array([
            Vp_a, Ap_a, Qp_a, Jp_a, I_a,
            Vp_b + k * V_b,
            Ap_b + k * A_b,
            Qp_b + k * Q_b,
            Jp_b + k * J_b,
            I_b - I_TARGET,          # DeepSeek Q22: target = 1.0
            V_a - V_norm,
        ])
        return residuals

    return bc


def initial_guess(r_grid, warm_start=True):
    """Build initial-guess profile.

    If warm_start=True: interpolate DeepSeek's reference profile (Q24) onto
    the grid. Derivatives computed by numerical gradient.

    If warm_start=False: fallback to v5's exp(-r) seed with non-proportional
    A,J. (Useful baseline — confirms v6's behavior is from the normalization
    fixes, not the warm-start.)
    """
    if warm_start:
        ref = DEEPSEEK_REF_PROFILE
        V = np.interp(r_grid, ref[:, 0], ref[:, 1])
        A = np.interp(r_grid, ref[:, 0], ref[:, 2])
        Q = np.interp(r_grid, ref[:, 0], ref[:, 3])
        J = np.interp(r_grid, ref[:, 0], ref[:, 4])
        # Numerical derivatives (forward diff at boundary, central elsewhere)
        Vp = np.gradient(V, r_grid)
        Ap = np.gradient(A, r_grid)
        Qp = np.gradient(Q, r_grid)
        Jp = np.gradient(J, r_grid)
    else:
        # v5 fallback seed: exp(-r) decay with non-proportional A,J
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

    # I(r) accumulator initialized as linear ramp from 0 to I_TARGET
    I = (r_grid - r_grid[0]) / (r_grid[-1] - r_grid[0]) * I_TARGET
    return np.vstack([V, Vp, A, Ap, Q, Qp, J, Jp, I])


# ---------------------------------------------------------------------------
# Energy functionals (Q23 — DeepSeek's H form vs Numerical-Benchmark form)
# ---------------------------------------------------------------------------

def compute_H_deepseek(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega, m_J_sq, g):
    """DeepSeek Q23 form (2026-05-20 ~4 PM):
        H = ∫ [ (1/2)(V'² + A'² + Q'² + J'²)
              + (1/2) ω² (V² + A² + Q² + J²)
              − V·Q + A·J
              + (g/4)( (V²+Q²)² + (A²+J²)² + 2(V·A − Q·J)² ) ] dr
    No (2π)²·R prefactor.
    """
    kinetic = 0.5 * (Vp * Vp + Ap * Ap + Qp * Qp + Jp * Jp)
    omega_kin = 0.5 * omega * omega * (V * V + A * A + Q * Q + J * J)
    cross = -V * Q + A * J
    quartic = 0.25 * g * ((V * V + Q * Q) ** 2
                          + (A * A + J * J) ** 2
                          + 2.0 * (V * A - Q * J) ** 2)
    h_dens = kinetic + omega_kin + cross + quartic
    return float(np.trapezoid(h_dens, r))


def compute_H_benchmark(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega, m_J_sq, lam_bench):
    """Numerical-Benchmark / v5 form (extracted in T5):
        H = (2π)² R · ∫ r dr [ (V'² + A'² + Q'² + J'²)
                              − V·Q + A·J
                              + (m_J²/2) (Q² − J²)
                              + (λ_bench/4) (Q² − J²)²
                              + (1/2) ω² (V² + A² + Q² + J²) ]
    With R=1, prefactor (2π)².
    """
    QQ_JJ = Q * Q - J * J
    h_dens = (Vp * Vp + Ap * Ap + Qp * Qp + Jp * Jp
              - V * Q + A * J
              + 0.5 * m_J_sq * QQ_JJ
              + 0.25 * lam_bench * QQ_JJ * QQ_JJ
              + 0.5 * omega * omega * (V * V + A * A + Q * Q + J * J))
    return float((2.0 * np.pi) ** 2 * np.trapezoid(h_dens * r, r))


# ---------------------------------------------------------------------------
# Post-processing
# ---------------------------------------------------------------------------

def analyze_solution(sol, m_J_sq, lam_bench, g, quartic_form, r_max):
    """Compute H (DeepSeek + benchmark forms for comparison), Q_CS, H/Q_CS,
    tail, peaks, node counts.
    """
    omega = float(sol.p[0])
    lam_LM = float(sol.p[1]) if len(sol.p) > 1 else 0.0
    r = sol.x
    V, Vp, A, Ap, Q, Qp, J, Jp, I = sol.y
    m_eff_sq = m_J_sq - omega * omega

    # Both energy functionals for comparison
    H_deepseek = compute_H_deepseek(r, V, Vp, A, Ap, Q, Qp, J, Jp,
                                     omega, m_J_sq, g)
    H_benchmark = compute_H_benchmark(r, V, Vp, A, Ap, Q, Qp, J, Jp,
                                       omega, m_J_sq, lam_bench)
    H = H_deepseek if quartic_form == "deepseek" else H_benchmark

    # Q_CS per DeepSeek Q22: ∫r·(A·J'-J·A')dr, no 2π factor
    Q_CS = float(I[-1])  # I(R_max) directly = Q_CS in DeepSeek convention
    Q_CS_grid = float(np.trapezoid(r * (A * Jp - J * Ap), r))

    HQ_CS = H / Q_CS if abs(Q_CS) > 1e-14 else float("nan")
    # Also compute H/Q for both forms, for diagnostics
    HQ_deepseek = H_deepseek / Q_CS if abs(Q_CS) > 1e-14 else float("nan")
    HQ_benchmark = H_benchmark / Q_CS if abs(Q_CS) > 1e-14 else float("nan")

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
    nodes_V, nodes_A = n_nodes(V), n_nodes(A)
    nodes_Q, nodes_J = n_nodes(Q), n_nodes(J)
    max_nodes = max(nodes_V, nodes_A, nodes_Q, nodes_J)

    return {
        "quartic_form": quartic_form,
        "omega": omega, "lambda_LM": lam_LM,
        "m_J_sq": m_J_sq, "m_eff_sq": float(m_eff_sq),
        "lam_bench": lam_bench, "g": g,
        "H_used": H, "H_deepseek": H_deepseek, "H_benchmark": H_benchmark,
        "Q_CS": Q_CS, "Q_CS_from_grid": Q_CS_grid,
        "HQ_CS": HQ_CS,
        "HQ_deepseek": HQ_deepseek, "HQ_benchmark": HQ_benchmark,
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
    """Apply v6 pass criterion."""
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
    parser.add_argument("--m-J-sq", type=float, default=M_J_SQ_DEFAULT)
    parser.add_argument("--lam-bench", type=float, default=LAMBDA_BENCH_DEFAULT)
    parser.add_argument("--g", type=float, default=1.0625,
                        help="g coefficient for DeepSeek quartic (default 1.0625)")
    parser.add_argument("--omega-init", type=float, default=OMEGA_INIT)
    parser.add_argument("--lambda-init", type=float, default=LAMBDA_LM_INIT)
    parser.add_argument("--V-norm", type=float, default=0.1)
    parser.add_argument("--r-max", type=float, default=R_MAX_DEFAULT)
    parser.add_argument("--n-grid", type=int, default=N_GRID_DEFAULT)
    parser.add_argument("--max-nodes", type=int, default=50000,
                        help="solve_bvp max_nodes (v5 attempt 4 used 50000)")
    parser.add_argument("--tol", type=float, default=5e-3,
                        help="solve_bvp tolerance")
    parser.add_argument("--verbose", type=int, default=2)
    parser.add_argument("--quartic", choices=["deepseek", "benchmark"],
                        default="deepseek",
                        help="Which quartic form to use for H. deepseek = "
                             "(g/4)·((V²+Q²)²+(A²+J²)²+2(VA-QJ)²). "
                             "benchmark = (λ/4)·(Q²-J²)² (v5 form, T5 extract)")
    parser.add_argument("--no-warm-start", action="store_true",
                        help="Skip DeepSeek reference profile; use exp(-r) seed")
    args = parser.parse_args()

    t0 = time.time()

    print(f"\n{'='*92}")
    print(f"M6 v6 — CALIBRATED BVP (DeepSeek normalizations)")
    print(f"{'='*92}")
    print(f"Fixed: m_J²={args.m_J_sq}, λ_bench={args.lam_bench}, g={args.g}, "
          f"V(R_MIN)={args.V_norm}")
    print(f"Free:  ω (init {args.omega_init}), λ_LM (init {args.lambda_init})")
    print(f"Grid:  r ∈ [{R_MIN}, {args.r_max}], {args.n_grid} points")
    print(f"Init:  {'DeepSeek 8-point reference profile' if not args.no_warm_start else 'v5 exp(-r) seed'}")
    print(f"Q_CS:  I_TARGET = {I_TARGET} (DeepSeek Q22: NO 2π factor)")
    print(f"H:     {args.quartic} quartic form, NO (2π)² prefactor (DeepSeek Q23)")
    print(f"Targets: H/Q_CS={TARGET_HQ}, Q_CS=1")

    # Build grid
    s = np.linspace(0, 1, args.n_grid)
    r_grid = R_MIN + (args.r_max - R_MIN) * (s ** 1.5)

    # Build initial guess
    y0 = initial_guess(r_grid, warm_start=(not args.no_warm_start))
    p0 = np.array([args.omega_init, args.lambda_init])

    # Build rhs and bc closures
    rhs = make_rhs(args.m_J_sq, args.lam_bench)
    bc = make_bc(args.m_J_sq, args.V_norm)

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

    diag = analyze_solution(sol, args.m_J_sq, args.lam_bench, args.g,
                             args.quartic, args.r_max)
    passed, fail_reasons = check_pass(diag)

    print(f"\n{'='*92}")
    print(f"CONVERGED SOLUTION DIAGNOSTICS")
    print(f"{'='*92}")
    print(f"  ω           = {diag['omega']:.6f}")
    print(f"  λ_LM        = {diag['lambda_LM']:.6f}")
    print(f"  m_eff²      = {diag['m_eff_sq']:.4f}")
    print(f"  V(R_MIN), A, Q, J = {diag['V0']:+.4f}, {diag['A0']:+.4f}, "
          f"{diag['Q0']:+.4f}, {diag['J0']:+.4f}")
    print(f"  peak V/A/Q/J     = {diag['peak_V']:.4f} / {diag['peak_A']:.4f} / "
          f"{diag['peak_Q']:.4f} / {diag['peak_J']:.4f}")
    print(f"  nodes V/A/Q/J    = {diag['nodes_V']} / {diag['nodes_A']} / "
          f"{diag['nodes_Q']} / {diag['nodes_J']}  (max {PASS_MAX_NODES})")
    print(f"  tail @r≥8        = {diag['tail']:.4f}")
    print(f"  Q_CS             = {diag['Q_CS']:.6f}  (target {TARGET_QCS})")
    print(f"  Q_CS from grid   = {diag['Q_CS_from_grid']:.6f}")
    print(f"")
    print(f"  H (DeepSeek)     = {diag['H_deepseek']:.6f}")
    print(f"  H (benchmark)    = {diag['H_benchmark']:.6f}")
    print(f"  H/Q (DeepSeek)   = {diag['HQ_deepseek']:.6f}  ← active if --quartic=deepseek")
    print(f"  H/Q (benchmark)  = {diag['HQ_benchmark']:.6f}")
    print(f"  H/Q USED         = {diag['HQ_CS']:.6f}  (target {TARGET_HQ})")

    print(f"\n  PASS: {'✅ YES' if passed else '❌ NO'}")
    for reason in fail_reasons:
        print(f"    - {reason}")

    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / "m6_v6_4fn_calibrated_bvp_results.json"
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
