"""M6 sandbox v7 — Ground-state BVP per Paul's 2026-05-21 reply (Q28/Q29/Q30).

Forks v6.6. ODE structure, normalizations, Lagrange-multiplier formulation —
all unchanged from v6 (validated by the 30× gap closure to H/Q=1.778, and the
0.84% gap from dropping the quartic in step (8) diagnostic).

Three changes from v6.6 driven by Paul's reply:

  Q28 — Quartic interpretation:
    "for the electron, the quartic term is negligible (<0.1% of H). For muon,
     tau the quartic becomes significant. For neutral chaoiton (Q=0), the
     quartic is essential. So the canonical form for H is the full expression
     (with quartic), but for the electron alone, the quartic term is a
     fine-tuning that can be ignored."

    v7 keeps the FULL quartic in the equations of motion (gives correct field
    shapes per Paul). For H computation, we report three values:
      H_full        = with DeepSeek quartic (canonical, for muon/tau/neutral)
      H_electron    = without quartic       (Paul's electron-specific form)
      H_benchmark   = v5 Numerical-Benchmark form (cross-check only)
    Primary acceptance metric: H_electron / Q_CS within 0.2% of 1.6969.

  Q29 — Node counts:
    "Our converged ground state for the electron has V and Q with exactly 4
     zero crossings (excluding r=0). The Lean ≤4-node spec refers to the
     radial functions excluding the origin. If your v6.6 has 5 crossings, you
     are on the first excited state."

    v7 tightens the acceptance criterion to max_nodes ≤ 4 strictly. If the
    solver still lands on the 5-node mode, we fall back to a mode-selector
    (node-penalty or Q(0)>0 hard BC).

  Q30 — Sign of Q at origin:
    "Our ground state has Q(0) = +0.1 (positive), same sign as V(0). The
     asymmetric helicity prescription is V₀ > 0, Q₀ > 0, A₀ < 0, J₀ < 0.
     Your v6.6 giving Q(0) negative suggests you are in a different basin —
     likely because your initial guess or the sign of the Lagrange multiplier
     term flipped the relative sign. To correct, set the initial profiles
     such that V(r) and Q(r) have the same sign (both positive) at r=0, and
     A(r) and J(r) have the same negative sign. The reference profile I gave
     earlier has exactly that."

    v7 ENABLES warm-start by default (Paul says the 8-point profile lands the
    ground state; v6.1 failed with warm-start at tighter r_max — v7 uses
    v6.6's loose config + warm-start, the combination not yet tested).
    Initial Lagrange multiplier: ω=1.0, λ_LM=12.0 (v6.6's converged value,
    biases toward the correct basin per the λ_LM sweep diagnostic).

ODE SYSTEM (9 states, unchanged from v6 — m_eff² = m_J² − ω²):
    dV/dr  = V'
    dV'/dr = -V'/r + Q
    dA/dr  = A'
    dA'/dr = -A'/r + J + λ_LM·(J + 2·r·J')
    dQ/dr  = Q'
    dQ'/dr = -Q'/r + V + m_eff²·Q + λ_bench·Q·(Q²-J²)
    dJ/dr  = J'
    dJ'/dr = -J'/r + A − m_eff²·J − λ_bench·J·(Q²-J²) − λ_LM·(A + 2·r·A')
    dI/dr  = r·(A·J' − J·A')        # Q_CS density per DeepSeek Q22

FIXED PARAMETERS:
    m_J²       = 0.5
    λ_bench    = 1.0
    g          = 1.0625    (only affects H_full diagnostic, not the ODE)

FREE EIGENVALUES:
    ω           init 1.0
    λ_LM        init 12.0    (v6.6 converged value — biases right basin)

BCs (11 total = 9 states + 2 free params):
    r = R_MIN:  V'=0, A'=0, Q'=0, J'=0, I=0          (5)
    r = R_max:  V'+kV=0, A'+kA=0, Q'+kQ=0, J'+kJ=0,
                I = 1.0                                (5)
    Anti-collapse: V(R_MIN) = V_norm                  (1)

PASS CRITERION (tighter than v6):
    |H_electron/Q_CS - 1.6969| < 0.003   (0.2% per Paul's tolerance)
    |Q_CS - 1| < 0.01
    max_nodes ≤ 4   (strict — Paul's Q29)
    Q(0) > 0        (strict — Paul's Q30)
    A(0) < 0, J(0) < 0  (helicity)

USAGE:
    python3 m6_v7_4fn_ground_state_bvp.py
    python3 m6_v7_4fn_ground_state_bvp.py --no-warm-start
    python3 m6_v7_4fn_ground_state_bvp.py --lambda-init 0.1   # v5/v6 default
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp


# ---------------------------------------------------------------------------
# Constants — v7 = v6 normalizations + Paul's Q28/Q29/Q30 guidance
# ---------------------------------------------------------------------------

M_J_SQ_DEFAULT = 0.5
LAMBDA_BENCH_DEFAULT = 1.0
OMEGA_INIT = 1.0
LAMBDA_LM_INIT = 12.0   # v6.6 converged value — biases toward correct basin
G_DEFAULT = 1.0625

R_MIN = 0.05
R_MAX_DEFAULT = 15.0      # v6.6's working r_max
N_GRID_DEFAULT = 500      # v6.6's working n_grid

I_TARGET = 1.0            # DeepSeek Q22

TARGET_HQ = 1.6969
TARGET_QCS = 1.0
PASS_HQ_TOL = 0.003       # 0.2% of 1.6969 per Paul's tolerance
PASS_QCS_TOL = 0.01
PASS_MAX_NODES = 4        # strict per Paul's Q29
LOC_R_CHECK = 8.0


# Paul's 8-point reference profile (DeepSeek Q24 — qualitative shape per Paul's Q30:
# "single peak, monotonic decay"). Signs already correct: V₀=+0.1, Q₀=+0.1,
# A₀=−0.1, J₀=−0.1 (asymmetric helicity).
PAUL_REF_PROFILE = np.array([
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
# ODE and BC — unchanged from v6
# ---------------------------------------------------------------------------

def make_rhs(m_J_sq_fixed, lam_bench_fixed, free_param="none"):
    """Return rhs(r, y, p) for scipy.integrate.solve_bvp.
    Full v6 ODE — quartic kept in EoM per Paul Q28.

    free_param:
        "none"      — 2 free params: ω, λ_LM (v6 default)
        "m_J"       — 3 free params: ω, λ_LM, m_J²
        "lam_bench" — 3 free params: ω, λ_LM, λ_bench
    """

    def rhs(r, y, p):
        omega = p[0]
        lam_LM = p[1]
        if free_param == "m_J":
            m_J_sq = p[2]
            lam_bench = lam_bench_fixed
        elif free_param == "lam_bench":
            m_J_sq = m_J_sq_fixed
            lam_bench = p[2]
        else:
            m_J_sq = m_J_sq_fixed
            lam_bench = lam_bench_fixed
        m_eff_sq = m_J_sq - omega * omega

        V, Vp, A, Ap, Q, Qp, J, Jp, I = y
        inv_r = 1.0 / r
        QQ_JJ = Q * Q - J * J

        dV = Vp
        dVp = -inv_r * Vp + Q
        dQ = Qp
        dQp = -inv_r * Qp + V + m_eff_sq * Q + lam_bench * Q * QQ_JJ

        dA = Ap
        dAp = -inv_r * Ap + J + lam_LM * (J + 2.0 * r * Jp)
        dJ = Jp
        dJp = (-inv_r * Jp + A - m_eff_sq * J
               - lam_bench * J * QQ_JJ
               - lam_LM * (A + 2.0 * r * Ap))

        dI = r * (A * Jp - J * Ap)

        return np.vstack([dV, dVp, dA, dAp, dQ, dQp, dJ, dJp, dI])

    return rhs


def make_bc(m_J_sq_fixed, V_norm, anchor="V", free_param="none"):
    """Return bc(ya, yb, p) for scipy.integrate.solve_bvp.

    anchor:
        "V" — pin V(R_MIN)=+V_norm only        (v6 default — 11 residuals,
              2 free params: ω, λ_LM)
        "Q" — pin Q(R_MIN)=+V_norm only        (v7 single-Q-pin — 11 res,
              2 free params: ω, λ_LM)
        "VQ" — pin BOTH V(R_MIN)=+V_norm AND   (v7 ground-state mode-selector
              Q(R_MIN)=+V_norm                  — 12 residuals, 3 free params:
                                                 ω, λ_LM, m_J²)

    The VQ anchor forces Werbos's same-sign V/Q positive helicity (Q30), and
    treats m_J² as an eigenvalue to absorb the extra constraint. If a clean
    ground state exists at some m_J² near 0.5, the solver finds it.
    """

    def bc(ya, yb, p):
        omega = p[0]
        if free_param == "m_J":
            m_J_sq = p[2]
        else:
            m_J_sq = m_J_sq_fixed
        k_squared = omega * omega - m_J_sq
        k = np.sqrt(max(k_squared, 0.01))

        V_a, Vp_a, A_a, Ap_a, Q_a, Qp_a, J_a, Jp_a, I_a = ya
        V_b, Vp_b, A_b, Ap_b, Q_b, Qp_b, J_b, Jp_b, I_b = yb

        base = [
            Vp_a, Ap_a, Qp_a, Jp_a, I_a,
            Vp_b + k * V_b,
            Ap_b + k * A_b,
            Qp_b + k * Q_b,
            Jp_b + k * J_b,
            I_b - I_TARGET,
        ]

        if anchor == "V":
            base.append(V_a - V_norm)
        elif anchor == "Q":
            base.append(Q_a - V_norm)
        elif anchor == "VQ":
            base.append(V_a - V_norm)
            base.append(Q_a - V_norm)
        else:
            raise ValueError(f"Unknown anchor: {anchor}")

        return np.array(base)

    return bc


def initial_guess(r_grid, warm_start=True):
    """Build initial-guess profile.

    warm_start=True: interpolate Paul's 8-point ref profile onto grid.
    Beyond r=8 (ref profile's last point), tail is already ≈0 so np.interp
    extrapolation as constant zero is fine.

    warm_start=False: fallback to v5/v6 exp(-r) seed.
    """
    if warm_start:
        ref = PAUL_REF_PROFILE
        V = np.interp(r_grid, ref[:, 0], ref[:, 1])
        A = np.interp(r_grid, ref[:, 0], ref[:, 2])
        Q = np.interp(r_grid, ref[:, 0], ref[:, 3])
        J = np.interp(r_grid, ref[:, 0], ref[:, 4])
        Vp = np.gradient(V, r_grid)
        Ap = np.gradient(A, r_grid)
        Qp = np.gradient(Q, r_grid)
        Jp = np.gradient(J, r_grid)
    else:
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

    I = (r_grid - r_grid[0]) / (r_grid[-1] - r_grid[0]) * I_TARGET
    return np.vstack([V, Vp, A, Ap, Q, Qp, J, Jp, I])


# ---------------------------------------------------------------------------
# Energy functionals — three forms per Paul's Q28
# ---------------------------------------------------------------------------

def compute_H_full(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega, g):
    """DeepSeek Q23 H form with the full quartic.
    Canonical for muon/tau/neutral chaoiton per Paul Q28.
        H = ∫ [ (1/2)(V'² + A'² + Q'² + J'²)
              + (1/2) ω² (V² + A² + Q² + J²)
              − V·Q + A·J
              + (g/4)( (V²+Q²)² + (A²+J²)² + 2(V·A − Q·J)² ) ] dr
    """
    kinetic = 0.5 * (Vp * Vp + Ap * Ap + Qp * Qp + Jp * Jp)
    omega_kin = 0.5 * omega * omega * (V * V + A * A + Q * Q + J * J)
    cross = -V * Q + A * J
    quartic = 0.25 * g * ((V * V + Q * Q) ** 2
                          + (A * A + J * J) ** 2
                          + 2.0 * (V * A - Q * J) ** 2)
    h_dens = kinetic + omega_kin + cross + quartic
    return float(np.trapezoid(h_dens, r))


def compute_H_electron(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega):
    """Paul Q28 — electron-specific H: drop the quartic.
    "For the electron alone, the quartic term is a fine-tuning that can be
     ignored."
        H_electron = ∫ [ (1/2)(V'² + A'² + Q'² + J'²)
                       + (1/2) ω² (V² + A² + Q² + J²)
                       − V·Q + A·J ] dr
    """
    kinetic = 0.5 * (Vp * Vp + Ap * Ap + Qp * Qp + Jp * Jp)
    omega_kin = 0.5 * omega * omega * (V * V + A * A + Q * Q + J * J)
    cross = -V * Q + A * J
    h_dens = kinetic + omega_kin + cross
    return float(np.trapezoid(h_dens, r))


def compute_H_benchmark(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega, m_J_sq, lam_bench):
    """Numerical-Benchmark v5 form — cross-check only."""
    QQ_JJ = Q * Q - J * J
    h_dens = (Vp * Vp + Ap * Ap + Qp * Qp + Jp * Jp
              - V * Q + A * J
              + 0.5 * m_J_sq * QQ_JJ
              + 0.25 * lam_bench * QQ_JJ * QQ_JJ
              + 0.5 * omega * omega * (V * V + A * A + Q * Q + J * J))
    return float((2.0 * np.pi) ** 2 * np.trapezoid(h_dens * r, r))


# ---------------------------------------------------------------------------
# Post-processing — node counting now strict (Paul Q29: excluding r=0)
# ---------------------------------------------------------------------------

def n_nodes_excluding_origin(f, r):
    """Count zero crossings of f(r) EXCLUDING the origin (per Paul Q29).
    Skip points near r=R_MIN to avoid counting an origin-region crossing.
    """
    mask = r > (r[0] + 0.1)  # skip first ~0.1 in r above R_MIN
    f_active = f[mask]
    f_active = f_active[np.abs(f_active) > 1e-6]
    if len(f_active) < 2:
        return 0
    return int(np.sum(np.diff(np.sign(f_active)) != 0))


def analyze_solution(sol, m_J_sq, lam_bench, g):
    """Compute all three H forms + diagnostics."""
    omega = float(sol.p[0])
    lam_LM = float(sol.p[1])
    r = sol.x
    V, Vp, A, Ap, Q, Qp, J, Jp, I = sol.y
    m_eff_sq = m_J_sq - omega * omega

    H_full = compute_H_full(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega, g)
    H_electron = compute_H_electron(r, V, Vp, A, Ap, Q, Qp, J, Jp, omega)
    H_benchmark = compute_H_benchmark(r, V, Vp, A, Ap, Q, Qp, J, Jp,
                                       omega, m_J_sq, lam_bench)

    Q_CS = float(I[-1])
    Q_CS_grid = float(np.trapezoid(r * (A * Jp - J * Ap), r))

    if abs(Q_CS) > 1e-14:
        HQ_full = H_full / Q_CS
        HQ_electron = H_electron / Q_CS
        HQ_benchmark = H_benchmark / Q_CS
    else:
        HQ_full = HQ_electron = HQ_benchmark = float("nan")

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

    nodes_V = n_nodes_excluding_origin(V, r)
    nodes_A = n_nodes_excluding_origin(A, r)
    nodes_Q = n_nodes_excluding_origin(Q, r)
    nodes_J = n_nodes_excluding_origin(J, r)
    max_nodes = max(nodes_V, nodes_A, nodes_Q, nodes_J)

    return {
        "omega": omega, "lambda_LM": lam_LM,
        "m_J_sq": m_J_sq, "m_eff_sq": float(m_eff_sq),
        "lam_bench": lam_bench, "g": g,
        "H_full": H_full, "H_electron": H_electron, "H_benchmark": H_benchmark,
        "Q_CS": Q_CS, "Q_CS_from_grid": Q_CS_grid,
        "HQ_full": HQ_full, "HQ_electron": HQ_electron,
        "HQ_benchmark": HQ_benchmark,
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
    """v7 pass criterion (tightened per Paul's reply).
    Primary metric: H_electron / Q_CS within 0.2% of 1.6969.
    Strict node + sign constraints.
    """
    fail = []
    hq_err = abs(diag["HQ_electron"] - TARGET_HQ)
    qcs_err = abs(diag["Q_CS"] - TARGET_QCS)
    if hq_err >= PASS_HQ_TOL:
        fail.append(f"|H_electron/Q-{TARGET_HQ}| = {hq_err:.4f} > {PASS_HQ_TOL}")
    if qcs_err >= PASS_QCS_TOL:
        fail.append(f"|Q_CS-1| = {qcs_err:.4f} > {PASS_QCS_TOL}")
    if diag["max_nodes"] > PASS_MAX_NODES:
        fail.append(f"max_nodes (excl r=0) = {diag['max_nodes']} > {PASS_MAX_NODES}")
    if diag["Q0"] <= 0:
        fail.append(f"Q(0) = {diag['Q0']:.4f} (should be > 0 per Paul Q30)")
    if diag["A0"] >= 0:
        fail.append(f"A(0) = {diag['A0']:.4f} (should be < 0 per helicity)")
    if diag["J0"] >= 0:
        fail.append(f"J(0) = {diag['J0']:.4f} (should be < 0 per helicity)")
    return (len(fail) == 0), fail


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--m-J-sq", type=float, default=M_J_SQ_DEFAULT)
    parser.add_argument("--lam-bench", type=float, default=LAMBDA_BENCH_DEFAULT)
    parser.add_argument("--g", type=float, default=G_DEFAULT)
    parser.add_argument("--omega-init", type=float, default=OMEGA_INIT)
    parser.add_argument("--lambda-init", type=float, default=LAMBDA_LM_INIT,
                        help="Lagrange multiplier init (default 12.0 = v6.6 converged)")
    parser.add_argument("--V-norm", type=float, default=0.1)
    parser.add_argument("--r-max", type=float, default=R_MAX_DEFAULT)
    parser.add_argument("--n-grid", type=int, default=N_GRID_DEFAULT)
    parser.add_argument("--max-nodes", type=int, default=100000,
                        help="solve_bvp max_nodes (default 100k = v6.6 sweet spot)")
    parser.add_argument("--tol", type=float, default=5e-3)
    parser.add_argument("--verbose", type=int, default=2)
    parser.add_argument("--no-warm-start", action="store_true",
                        help="Skip Paul's 8-point profile; use exp(-r) seed")
    parser.add_argument("--anchor", choices=["V", "Q", "VQ"], default="V",
                        help="Anti-collapse anchor: V pins V(R_MIN)=V_norm; "
                             "Q pins Q(R_MIN)=V_norm; VQ pins BOTH (requires "
                             "--free-param to release a 3rd eigenvalue)")
    parser.add_argument("--free-param", choices=["none", "m_J", "lam_bench"],
                        default="none",
                        help="3rd free eigenvalue when --anchor=VQ. "
                             "m_J releases m_J²; lam_bench releases λ_bench.")
    args = parser.parse_args()

    t0 = time.time()

    print(f"\n{'='*92}")
    print(f"M6 v7 — GROUND-STATE BVP (Paul's 2026-05-21 Q28/Q29/Q30 reply)")
    print(f"{'='*92}")
    print(f"Fixed:  m_J²={args.m_J_sq}, λ_bench={args.lam_bench}, "
          f"g={args.g} (for H_full only), V(R_MIN)={args.V_norm}")
    print(f"Free:   ω (init {args.omega_init}), λ_LM (init {args.lambda_init})")
    print(f"Grid:   r ∈ [{R_MIN}, {args.r_max}], {args.n_grid} points")
    print(f"Init:   {'Paul 8-point reference profile' if not args.no_warm_start else 'exp(-r) seed'}")
    print(f"Anchor: {args.anchor}(R_MIN) = {args.V_norm}  "
          f"({'v6 anti-collapse' if args.anchor == 'V' else 'v7 ground-state mode-selector (Q sign-pinned)'})")
    print(f"Q_CS:   I_TARGET = {I_TARGET}")
    print(f"H:      ALL THREE FORMS computed; H_electron/Q is the primary metric")
    print(f"Target: H_electron/Q = {TARGET_HQ} (±{PASS_HQ_TOL/TARGET_HQ*100:.1f}%); "
          f"nodes ≤ {PASS_MAX_NODES} (excl r=0); Q(0)>0; A,J(0)<0")

    s = np.linspace(0, 1, args.n_grid)
    r_grid = R_MIN + (args.r_max - R_MIN) * (s ** 1.5)

    y0 = initial_guess(r_grid, warm_start=(not args.no_warm_start))
    if args.anchor == "VQ":
        if args.free_param == "m_J":
            p0 = np.array([args.omega_init, args.lambda_init, args.m_J_sq])
            print(f"Free³:  m_J² (init {args.m_J_sq})")
        elif args.free_param == "lam_bench":
            p0 = np.array([args.omega_init, args.lambda_init, args.lam_bench])
            print(f"Free³:  λ_bench (init {args.lam_bench})")
        else:
            raise SystemExit("--anchor=VQ requires --free-param m_J or lam_bench")
    else:
        p0 = np.array([args.omega_init, args.lambda_init])

    rhs = make_rhs(args.m_J_sq, args.lam_bench, free_param=args.free_param)
    bc = make_bc(args.m_J_sq, args.V_norm, anchor=args.anchor,
                 free_param=args.free_param)

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
    m_J_sq_used = args.m_J_sq
    lam_bench_used = args.lam_bench
    if len(sol.p) > 2:
        if args.free_param == "m_J":
            m_J_sq_used = float(sol.p[2])
            print(f"  Final ω = {sol.p[0]:.6f}, λ_LM = {sol.p[1]:.6f}, "
                  f"m_J² = {m_J_sq_used:.6f}")
        elif args.free_param == "lam_bench":
            lam_bench_used = float(sol.p[2])
            print(f"  Final ω = {sol.p[0]:.6f}, λ_LM = {sol.p[1]:.6f}, "
                  f"λ_bench = {lam_bench_used:.6f}")
    else:
        print(f"  Final ω = {sol.p[0]:.6f}, λ_LM = {sol.p[1]:.6f}")

    diag = analyze_solution(sol, m_J_sq_used, lam_bench_used, args.g)
    passed, fail_reasons = check_pass(diag)

    print(f"\n{'='*92}")
    print(f"CONVERGED SOLUTION DIAGNOSTICS")
    print(f"{'='*92}")
    print(f"  ω           = {diag['omega']:.6f}")
    print(f"  λ_LM        = {diag['lambda_LM']:.6f}")
    print(f"  m_eff²      = {diag['m_eff_sq']:.4f}")
    print(f"  V(R_MIN), A, Q, J = {diag['V0']:+.4f}, {diag['A0']:+.4f}, "
          f"{diag['Q0']:+.4f}, {diag['J0']:+.4f}")
    print(f"  peak V/A/Q/J      = {diag['peak_V']:.4f} / {diag['peak_A']:.4f} / "
          f"{diag['peak_Q']:.4f} / {diag['peak_J']:.4f}")
    print(f"  nodes V/A/Q/J (excl r=0) = {diag['nodes_V']} / {diag['nodes_A']} / "
          f"{diag['nodes_Q']} / {diag['nodes_J']}  (max {PASS_MAX_NODES})")
    print(f"  tail @r≥8         = {diag['tail']:.4f}")
    print(f"  Q_CS              = {diag['Q_CS']:.6f}  (target {TARGET_QCS})")
    print(f"  Q_CS from grid    = {diag['Q_CS_from_grid']:.6f}")
    print(f"")
    print(f"  H_full     (w/ quartic)  = {diag['H_full']:.6f}")
    print(f"  H_electron (no quartic)  = {diag['H_electron']:.6f}")
    print(f"  H_benchmark              = {diag['H_benchmark']:.6f}")
    print(f"  H_full/Q     = {diag['HQ_full']:.6f}")
    print(f"  H_electron/Q = {diag['HQ_electron']:.6f}  ★ PRIMARY (target {TARGET_HQ})")
    print(f"  H_bench/Q    = {diag['HQ_benchmark']:.6f}")

    print(f"\n  PASS: {'✅ YES' if passed else '❌ NO'}")
    for reason in fail_reasons:
        print(f"    - {reason}")

    out_dir = Path(__file__).resolve().parent
    out_path = out_dir / "m6_v7_4fn_ground_state_bvp_results.json"
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
