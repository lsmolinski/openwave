"""M6 sandbox v4 T9 — Newton-residual shooter for Q_CS=1 chaoiton (2026-05-20).

After T8 closed with "helicity necessary but not sufficient," this script
attempts the next logical step: drive the asymmetric-helicity 4-fn IVP toward
the canonical electron chaoiton by minimizing a residual vector against
Werbos's two anchors:

    r_HQ  = H/Q_CS - 1.6969   (electron calibration target)
    r_QCS = Q_CS - 1          (integer winding)

with soft regularizers for localization and to prevent the BVP-style (A,J)→0
collapse seen in T8c.

VARIABLES (6 total):
    theta = (m_J², λ_bench, V₀_abs, A₀_abs, Q₀_abs, J₀_abs)
    All positive magnitudes. Helicity signs applied in evaluation wrapper:
        V₀ = +V₀_abs,  A₀ = -A₀_abs,  Q₀ = +Q₀_abs,  J₀ = -J₀_abs.

RESIDUALS (6):
    r1 = H/Q_CS - 1.6969              # calibration anchor
    r2 = Q_CS - 1                     # integer winding
    r3 = tail                         # localization (∥field∥ at r ≥ 8)
    r4 = max(0, 0.01 - peak_A)        # anti-A-collapse penalty
    r5 = max(0, 0.01 - peak_J)        # anti-J-collapse penalty
    r6 = PENALTY · max(0, 1 - r_reached/r_max)   # blowup distance penalty
    All 6 residuals computed even when IVP terminates early — partial
    trajectory still gives smooth gradient. r6 drives optimizer toward
    parameter regions where IVP reaches r_max (i.e., bound state basin).

BOUNDS:
    m_J²       ∈ [0.01, 10]
    λ_bench    ∈ [0.01, 10]
    amplitudes ∈ [0.001, 1.0]

INITIAL GUESS (Werbos canonical):
    (m_J²=0.5, λ_bench=1.0, V₀=Q₀=0.1, A₀=J₀=0.1)   # signs locked → ±0.1

ALGORITHM:
    Stage 1: differential_evolution (gradient-free, global) to find the
             bound-state basin where IVP reaches r_max. Cost = blow-up
             distance only — pure exploration over (m_J², λ, amps).
    Stage 2: Nelder-Mead from Stage 1 winner to refine calibration:
             cost combines (H/Q − 1.6969)² + (Q_CS − 1)² + tail² + collapse
             penalties. Local polish once we are inside the basin.

    Two stages because the landscape has a discontinuous transition at
    r_reached = r_max: outside the basin only r_reached has gradient;
    inside the basin r_reached saturates and calibration residuals
    become meaningful. One-shot least_squares cannot bridge the two.

PASS CRITERION (per reload prompt T9 design):
    |H/Q - 1.6969| < 0.001  AND  |Q_CS - 1| < 0.01  AND
    tail < 0.05             AND  no blowup           AND  nodes ≤ 4

USAGE:
    python3 m6_v4_4fn_newton.py
    python3 m6_v4_4fn_newton.py --r-max 20 --max-nfev 200 --n-starts 4
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution, minimize


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OMEGA_DEF = 1.0
R_MIN = 0.05
R_MAX_DEF = 15.0
N_GRID_DEF = 1500
RTOL = 1e-8
ATOL = 1e-10
MAX_STEP = 0.05
BLOWUP_THRESHOLD = 100.0

LOC_R_CHECK = 8.0
LOC_THRESHOLD = 0.05

TARGET_HQ = 1.6969
TARGET_QCS = 1.0

PENALTY = 50.0    # blowup penalty magnitude (scaled by 1 - r_reached/r_max)
ANTI_COLLAPSE_FLOOR = 0.01   # below this peak amplitude → soft penalty

# Pass criterion
PASS_HQ_TOL = 0.001
PASS_QCS_TOL = 0.01
PASS_TAIL_TOL = 0.05
PASS_MAX_NODES = 4


# ---------------------------------------------------------------------------
# ODE: 4-function benchmark (toroidal, no -f/r² term)
# ---------------------------------------------------------------------------

def ode_4fn(r, y, m_eff_sq, lam_bench):
    V, dV, A, dA, Q, dQ, J, dJ = y
    inv_r = 1.0 / r
    QQ_JJ = Q * Q - J * J
    d2V = -inv_r * dV + Q
    d2A = -inv_r * dA + J
    d2Q = -inv_r * dQ + V + m_eff_sq * Q + lam_bench * Q * QQ_JJ
    d2J = -inv_r * dJ + A - m_eff_sq * J - lam_bench * J * QQ_JJ
    return [dV, d2V, dA, d2A, dQ, d2Q, dJ, d2J]


def integrate(V0, A0, Q0, J0, m_J_sq, lam_bench,
              omega=OMEGA_DEF, r_max=R_MAX_DEF, n_grid=N_GRID_DEF,
              include_omega_kinetic=True):
    """Integrate 4-fn IVP with blowup-event termination.

    Returns dict with H, Q_CS, HQ_CS, tail, peak_*, blew_up, r_reached.
    """
    m_eff_sq = m_J_sq - omega * omega
    y0 = [V0, 0.0, A0, 0.0, Q0, 0.0, J0, 0.0]
    r_eval = np.linspace(R_MIN, r_max, n_grid)

    def blowup_event(r, y):
        m = max(abs(y[0]), abs(y[2]), abs(y[4]), abs(y[6]))
        return BLOWUP_THRESHOLD - m
    blowup_event.terminal = True
    blowup_event.direction = -1

    try:
        sol = solve_ivp(
            fun=lambda r, y: ode_4fn(r, y, m_eff_sq, lam_bench),
            t_span=(R_MIN, r_max), y0=y0, t_eval=r_eval,
            method="RK45", rtol=RTOL, atol=ATOL,
            max_step=MAX_STEP, events=blowup_event,
        )
    except Exception as e:
        return {"success": False, "blew_up": True, "r_reached": R_MIN,
                "error": str(e)[:80]}

    if not sol.success:
        return {"success": False, "blew_up": True,
                "r_reached": float(sol.t[-1]) if len(sol.t) else R_MIN,
                "error": sol.message[:80]}

    r_reached = float(sol.t[-1])
    blew_up = r_reached < r_max - 0.01

    r = sol.t
    V, dV, A, dA, Q, dQ, J, dJ = sol.y

    if not all(np.all(np.isfinite(f)) for f in (V, A, Q, J)):
        return {"success": False, "blew_up": True, "r_reached": r_reached,
                "error": "non-finite fields"}

    # Localization tail at r >= LOC_R_CHECK (skipped if blew_up too early)
    mask = r >= LOC_R_CHECK
    if mask.sum() == 0:
        max_abs = (np.abs(V) + np.abs(A) + np.abs(Q) + np.abs(J)).max()
        tail = float(max_abs)
    else:
        max_abs = (np.abs(V[mask]) + np.abs(A[mask])
                   + np.abs(Q[mask]) + np.abs(J[mask])).max()
        tail = float(max_abs)

    peak_V = float(np.abs(V).max())
    peak_A = float(np.abs(A).max())
    peak_Q = float(np.abs(Q).max())
    peak_J = float(np.abs(J).max())

    # Node counts (sign changes per field)
    def n_nodes(f):
        return int(np.sum(np.diff(np.sign(f[np.abs(f) > 1e-6])) != 0))
    nodes_V = n_nodes(V)
    nodes_A = n_nodes(A)
    nodes_Q = n_nodes(Q)
    nodes_J = n_nodes(J)
    max_nodes = max(nodes_V, nodes_A, nodes_Q, nodes_J)

    # Energy functional (benchmark §5, R=1)
    QQ_JJ = Q * Q - J * J
    h_dens = (dV * dV + dA * dA + dQ * dQ + dJ * dJ
              - V * Q + A * J
              + 0.5 * m_J_sq * QQ_JJ
              + 0.25 * lam_bench * QQ_JJ * QQ_JJ)
    if include_omega_kinetic:
        h_dens = h_dens + 0.5 * omega * omega * (V * V + A * A + Q * Q + J * J)
    H = float((2.0 * np.pi) ** 2 * np.trapezoid(h_dens * r, r))

    # Chern-Simons charge
    qcs_dens = A * dJ - J * dA
    Q_CS = float(2.0 * np.pi * np.trapezoid(qcs_dens * r, r))

    HQ_CS = H / Q_CS if abs(Q_CS) > 1e-14 else float("nan")

    return {
        "success": True,
        "V0": V0, "A0": A0, "Q0": Q0, "J0": J0,
        "m_J_sq": m_J_sq, "m_eff_sq": m_eff_sq, "lam_bench": lam_bench,
        "H": H, "Q_CS": Q_CS, "HQ_CS": HQ_CS,
        "tail": tail, "blew_up": blew_up, "r_reached": r_reached,
        "peak_V": peak_V, "peak_A": peak_A, "peak_Q": peak_Q, "peak_J": peak_J,
        "nodes_V": nodes_V, "nodes_A": nodes_A,
        "nodes_Q": nodes_Q, "nodes_J": nodes_J, "max_nodes": max_nodes,
        "r": r, "V": V, "A": A, "Q": Q, "J": J,
    }


# ---------------------------------------------------------------------------
# Residual function with helicity sign wrapper
# ---------------------------------------------------------------------------

def evaluate_theta(theta, omega, r_max):
    """Run integrate() at theta with helicity signs applied.
    Returns the full result dict (or failure stub).
    """
    m_J_sq, lam_bench, V_abs, A_abs, Q_abs, J_abs = theta
    V0, A0, Q0, J0 = +V_abs, -A_abs, +Q_abs, -J_abs
    return integrate(V0, A0, Q0, J0, m_J_sq, lam_bench,
                     omega=omega, r_max=r_max)


def cost_stage1(theta, omega, r_max, history):
    """Stage 1 cost: reach r_max with bounded fields.

    cost = PENALTY² · deficit²  +  small log-peak regularizer

    Pure deficit dominates outside basin; log-peak adds gentle pressure
    against blowup-delaying tricks (very small λ + delayed divergence).
    """
    res = evaluate_theta(theta, omega, r_max)
    if not res.get("success"):
        r_reached = res.get("r_reached", R_MIN)
        cost = PENALTY * PENALTY * (1.5 - r_reached / r_max) ** 2
        history.append({"stage": 1, "theta": list(theta), "success": False,
                        "r_reached": r_reached, "cost": float(cost)})
        return float(cost)
    r_reached = res["r_reached"]
    deficit = max(0.0, 1.0 - r_reached / r_max)
    max_peak = max(res["peak_V"], res["peak_A"], res["peak_Q"], res["peak_J"])

    cost_reach = (PENALTY * deficit) ** 2
    # gentle regularizer — capped contribution at log(1+100²) ≈ 9
    cost_peak = 1.0 * np.log(1.0 + max_peak * max_peak)
    cost = float(cost_reach + cost_peak)

    history.append({
        "stage": 1, "theta": list(theta), "success": True,
        "blew_up": res["blew_up"], "r_reached": r_reached,
        "max_peak": max_peak, "tail": res["tail"],
        "cost": cost,
    })
    return cost


def cost_stage2(theta, omega, r_max, history):
    """Stage 2 cost: calibration residuals (assuming we are in the basin).
    Includes a soft blow-up penalty so the optimizer cannot wander out.
    """
    res = evaluate_theta(theta, omega, r_max)
    if not res.get("success"):
        r_reached = res.get("r_reached", R_MIN)
        cost = 1e6 * (1.0 - r_reached / r_max + 0.1)
        history.append({"stage": 2, "theta": list(theta), "success": False,
                        "r_reached": r_reached, "cost": float(cost)})
        return float(cost)

    r_reached = res["r_reached"]
    blew_up = res["blew_up"]

    if blew_up:
        # Punish leaving the basin during Stage 2
        cost = 1e4 + PENALTY * PENALTY * (1.0 - r_reached / r_max) ** 2
        history.append({"stage": 2, "theta": list(theta), "success": True,
                        "blew_up": True, "r_reached": r_reached,
                        "cost": float(cost)})
        return float(cost)

    hq = res.get("HQ_CS", float("nan"))
    qcs = res["Q_CS"]
    tail = res["tail"]
    peak_A = res["peak_A"]
    peak_J = res["peak_J"]

    r1 = (hq - TARGET_HQ) if np.isfinite(hq) else 100.0
    r2 = qcs - TARGET_QCS
    r3 = tail
    r4 = max(0.0, ANTI_COLLAPSE_FLOOR - peak_A) * 100  # weight collapses high
    r5 = max(0.0, ANTI_COLLAPSE_FLOOR - peak_J) * 100
    cost = r1 * r1 + r2 * r2 + r3 * r3 + r4 * r4 + r5 * r5
    history.append({
        "stage": 2, "theta": list(theta), "success": True, "blew_up": False,
        "r_reached": r_reached,
        "HQ_CS": hq if np.isfinite(hq) else None,
        "Q_CS": qcs, "tail": tail,
        "peak_A": peak_A, "peak_J": peak_J,
        "cost": float(cost),
    })
    return float(cost)


# ---------------------------------------------------------------------------
# Newton solve (multi-start)
# ---------------------------------------------------------------------------

def run_stage1_global(omega, r_max, max_iter, popsize, seed, history):
    """Stage 1: differential_evolution search for the bound-state basin."""
    bounds = [
        (0.01, 10.0),    # m_J²
        (0.01, 10.0),    # λ_bench
        (0.001, 1.0),    # |V₀|
        (0.001, 1.0),    # |A₀|
        (0.001, 1.0),    # |Q₀|
        (0.001, 1.0),    # |J₀|
    ]
    print(f"\n  [stage1] differential_evolution: "
          f"maxiter={max_iter}, popsize={popsize}, seed={seed}")
    try:
        result = differential_evolution(
            cost_stage1, bounds=bounds, args=(omega, r_max, history),
            maxiter=max_iter, popsize=popsize, seed=seed,
            tol=0.01, mutation=(0.5, 1.0), recombination=0.7,
            init='sobol', polish=False, workers=1,
        )
    except Exception as e:
        print(f"  [stage1] EXCEPTION: {str(e)[:80]}")
        return None
    print(f"  [stage1] terminated after {result.nfev} evals, "
          f"cost {result.fun:.4f}")
    print(f"  [stage1] θ* = {[f'{x:.4f}' for x in result.x]}")
    # How close did we get to reaching r_max?
    res_eval = evaluate_theta(result.x, omega, r_max)
    if res_eval.get("success"):
        print(f"  [stage1] r_reached = {res_eval['r_reached']:.3f} / "
              f"{r_max} (blew_up = {res_eval['blew_up']})")
    return result


def run_stage2_local(theta0, omega, r_max, max_iter, history):
    """Stage 2: Nelder-Mead refinement of calibration residuals."""
    print(f"\n  [stage2] Nelder-Mead from θ₀ = "
          f"{[f'{x:.4f}' for x in theta0]}")
    # Simple initial-simplex perturbation
    simplex = np.array([theta0] + [theta0 + 0.05 * np.eye(6)[i]
                                    for i in range(6)])
    # Keep simplex inside bounds
    simplex = np.clip(simplex,
                      [0.011, 0.011, 0.0011, 0.0011, 0.0011, 0.0011],
                      [9.99, 9.99, 0.99, 0.99, 0.99, 0.99])
    try:
        result = minimize(
            cost_stage2, theta0, args=(omega, r_max, history),
            method="Nelder-Mead",
            options={"initial_simplex": simplex, "maxiter": max_iter,
                     "xatol": 1e-5, "fatol": 1e-6, "adaptive": True},
        )
    except Exception as e:
        print(f"  [stage2] EXCEPTION: {str(e)[:80]}")
        return None
    print(f"  [stage2] terminated after {result.nfev} evals, "
          f"cost {result.fun:.6f}")
    print(f"  [stage2] θ* = {[f'{x:.4f}' for x in result.x]}")
    return result


def evaluate_at_best(theta_star, omega, r_max):
    """Re-integrate at the best theta to get full field shape for reporting."""
    return evaluate_theta(theta_star, omega, r_max)


def coarse_scan_fixed_amps(omega, r_max, V_abs=0.1, A_abs=0.1,
                            Q_abs=0.1, J_abs=0.1,
                            m_J_sq_grid=None, lam_grid=None):
    """Pre-Stage scan: at Werbos's exact amplitudes, map (m_J², λ_bench).

    Goal: does ANY (m_J², λ) reach r_max with bounded fields?
    If yes → seed Stage 2 from there.
    If no → the bound state requires asymmetric magnitudes too.
    """
    if m_J_sq_grid is None:
        # Wider + log-spaced to catch small λ regime
        m_J_sq_grid = np.geomspace(0.05, 10.0, 18)
    if lam_grid is None:
        lam_grid = np.geomspace(0.01, 10.0, 18)

    print(f"\n  [scan] coarse (m_J², λ_bench) grid at |amps|={V_abs}")
    print(f"  [scan] {len(m_J_sq_grid)}×{len(lam_grid)} = "
          f"{len(m_J_sq_grid)*len(lam_grid)} points, r_max={r_max}")

    R = np.zeros((len(m_J_sq_grid), len(lam_grid)))
    peakmat = np.zeros_like(R)
    best_r = R_MIN
    best_pt = None

    for i, m_J_sq in enumerate(m_J_sq_grid):
        for j, lam in enumerate(lam_grid):
            theta = np.array([m_J_sq, lam, V_abs, A_abs, Q_abs, J_abs])
            res = evaluate_theta(theta, omega, r_max)
            if not res.get("success"):
                R[i, j] = R_MIN
                peakmat[i, j] = float("inf")
                continue
            R[i, j] = res["r_reached"]
            peakmat[i, j] = max(res["peak_V"], res["peak_A"],
                                res["peak_Q"], res["peak_J"])
            if res["r_reached"] > best_r:
                best_r = res["r_reached"]
                best_pt = {"theta": theta.tolist(), "result": res}

    print(f"  [scan] max r_reached = {best_r:.3f} (target r_max = {r_max})")
    if best_pt is not None:
        best_theta = best_pt["theta"]
        print(f"  [scan] at θ = m_J²={best_theta[0]:.3f}, "
              f"λ={best_theta[1]:.3f}, "
              f"peak={max(best_pt['result']['peak_V'], best_pt['result']['peak_A'], best_pt['result']['peak_Q'], best_pt['result']['peak_J']):.3f}")
    print(f"  [scan] count(r_reached ≥ r_max - 0.5): "
          f"{int(np.sum(R >= r_max - 0.5))}")
    print(f"  [scan] count(peak < 5):        {int(np.sum(peakmat < 5.0))}")

    return {
        "m_J_sq_grid": m_J_sq_grid.tolist(),
        "lam_grid": lam_grid.tolist(),
        "r_reached_matrix": R.tolist(),
        "peak_matrix": peakmat.tolist(),
        "best": ({"theta": best_pt["theta"],
                  "r_reached": best_r,
                  "peak_V": best_pt["result"]["peak_V"],
                  "peak_A": best_pt["result"]["peak_A"],
                  "peak_Q": best_pt["result"]["peak_Q"],
                  "peak_J": best_pt["result"]["peak_J"]}
                 if best_pt is not None else None),
    }


def check_pass(res):
    """Apply the T9 pass criterion. Returns (passed, reasons)."""
    if not res or not res["success"]:
        return False, ["integration failed"]
    if res["blew_up"]:
        return False, [f"blew up at r={res['r_reached']:.2f}"]
    fail = []
    hq_err = abs(res["HQ_CS"] - TARGET_HQ)
    qcs_err = abs(res["Q_CS"] - TARGET_QCS)
    if hq_err >= PASS_HQ_TOL:
        fail.append(f"|H/Q-{TARGET_HQ}| = {hq_err:.4f} > {PASS_HQ_TOL}")
    if qcs_err >= PASS_QCS_TOL:
        fail.append(f"|Q_CS-1| = {qcs_err:.4f} > {PASS_QCS_TOL}")
    if res["tail"] >= PASS_TAIL_TOL:
        fail.append(f"tail = {res['tail']:.4f} > {PASS_TAIL_TOL}")
    if res["max_nodes"] > PASS_MAX_NODES:
        fail.append(f"max_nodes = {res['max_nodes']} > {PASS_MAX_NODES}")
    return (len(fail) == 0), fail


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--omega", type=float, default=OMEGA_DEF)
    parser.add_argument("--r-max", type=float, default=R_MAX_DEF)
    parser.add_argument("--de-maxiter", type=int, default=30,
                        help="Stage 1 differential_evolution max iterations")
    parser.add_argument("--de-popsize", type=int, default=15,
                        help="Stage 1 differential_evolution popsize (× n_vars)")
    parser.add_argument("--nm-maxiter", type=int, default=400,
                        help="Stage 2 Nelder-Mead max iterations")
    parser.add_argument("--seed", type=int, default=20260520)
    args = parser.parse_args()

    t0 = time.time()

    print(f"\n{'='*92}")
    print(f"M6 v4 T9 — Q_CS=1 CHAOITON SHOOTER (two-stage + pre-scan)")
    print(f"{'='*92}")
    print(f"ω = {args.omega}, r_max = {args.r_max}")
    print(f"Targets: H/Q_CS = {TARGET_HQ}, Q_CS = {TARGET_QCS}, tail < {LOC_THRESHOLD}")
    print(f"Pass: |H/Q-{TARGET_HQ}|<{PASS_HQ_TOL}, |Q_CS-1|<{PASS_QCS_TOL}, "
          f"tail<{PASS_TAIL_TOL}, nodes≤{PASS_MAX_NODES}, no blowup")

    # ---- Pre-scan: does the basin exist at Werbos's exact amps? ----
    scan = coarse_scan_fixed_amps(args.omega, args.r_max)

    # ---- Stage 1: find bound-state basin via global search ----
    stage1_history = []
    s1 = run_stage1_global(args.omega, args.r_max,
                           args.de_maxiter, args.de_popsize,
                           args.seed, stage1_history)
    if s1 is None:
        print("\n  [stage1] FAILED")
        return 1

    s1_eval = evaluate_theta(s1.x, args.omega, args.r_max)
    in_basin = (s1_eval.get("success")
                and not s1_eval.get("blew_up", True))

    # ---- Stage 2: local polish ----
    stage2_history = []
    if in_basin:
        print(f"\n  [stage1] Inside basin (r_reached = {s1_eval['r_reached']:.3f}). "
              f"Refining with Nelder-Mead.")
        s2 = run_stage2_local(s1.x, args.omega, args.r_max,
                              args.nm_maxiter, stage2_history)
        final_theta = s2.x if s2 is not None else s1.x
    else:
        r_r = s1_eval.get("r_reached", R_MIN)
        print(f"\n  [stage1] DID NOT reach r_max (got {r_r:.3f}). "
              f"Skipping Stage 2 (Nelder-Mead would amplify blow-up noise).")
        final_theta = s1.x
        s2 = None

    final = evaluate_at_best(final_theta, args.omega, args.r_max)
    passed, fail_reasons = check_pass(final)

    print(f"\n{'='*92}")
    print(f"FINAL θ* = {[f'{x:.4f}' for x in final_theta]}")
    print(f"{'='*92}")

    if final.get("success"):
        print(f"  m_J²       = {final['m_J_sq']:.4f}")
        print(f"  m_eff²     = {final['m_eff_sq']:.4f}")
        print(f"  λ_bench    = {final['lam_bench']:.4f}")
        print(f"  V₀ A₀ Q₀ J₀ = {final['V0']:+.4f} {final['A0']:+.4f} "
              f"{final['Q0']:+.4f} {final['J0']:+.4f}")
        print(f"  H          = {final['H']:.4f}")
        print(f"  Q_CS       = {final['Q_CS']:.4f}    (target {TARGET_QCS})")
        print(f"  H/Q_CS     = {final['HQ_CS']:.4f}    (target {TARGET_HQ})")
        print(f"  tail       = {final['tail']:.4f}    (target < {PASS_TAIL_TOL})")
        print(f"  peak V/A/Q/J = {final['peak_V']:.3f}/{final['peak_A']:.3f}/"
              f"{final['peak_Q']:.3f}/{final['peak_J']:.3f}")
        print(f"  nodes V/A/Q/J = {final['nodes_V']}/{final['nodes_A']}/"
              f"{final['nodes_Q']}/{final['nodes_J']}  (max {PASS_MAX_NODES})")
        print(f"  r_reached  = {final['r_reached']:.3f}  (r_max {args.r_max})")
        print(f"  blew_up    = {final['blew_up']}")
    else:
        print(f"  Final eval FAILED: {final.get('error', '?')}")
        print(f"  r_reached = {final.get('r_reached', '?')}")

    print(f"\n  PASS: {'✅ YES' if passed else '❌ NO'}")
    if not passed:
        for reason in fail_reasons:
            print(f"    - {reason}")

    # Persist
    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / "m6_v4_4fn_newton_results.json"
    final_dump = {k: v for k, v in (final.items() if final.get("success") else [])
                  if k not in ("r", "V", "A", "Q", "J")}
    payload = {
        "config": {
            "omega": args.omega, "r_max": args.r_max,
            "de_maxiter": args.de_maxiter, "de_popsize": args.de_popsize,
            "nm_maxiter": args.nm_maxiter, "seed": args.seed,
        },
        "targets": {"H/Q_CS": TARGET_HQ, "Q_CS": TARGET_QCS,
                    "tail_max": PASS_TAIL_TOL, "max_nodes": PASS_MAX_NODES},
        "prescan": scan,
        "stage1": {
            "theta_star": s1.x.tolist(), "cost": float(s1.fun),
            "nfev": int(s1.nfev), "n_evals_recorded": len(stage1_history),
            "in_basin": in_basin,
        },
        "stage2": ({"theta_star": s2.x.tolist(), "cost": float(s2.fun),
                    "nfev": int(s2.nfev),
                    "n_evals_recorded": len(stage2_history)}
                   if s2 is not None else None),
        "final": final_dump,
        "passed": passed,
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
