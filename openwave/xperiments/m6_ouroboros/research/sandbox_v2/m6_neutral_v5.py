"""M6 neutral chaoiton — v5-paper-corrected (2026-05-18).

Third attempt at the Q=0 chaoiton. Sandbox v2's first two attempts (IVP and
BVP with locked ansatz J=±A, Q=±V) both failed. Werbos's "Law of Everything
v5" PDF (received 2026-05-18 17:20) revealed three structural corrections:

  1. ASYMMETRIC 4-function ansatz (v5 §4 eqs 4.2, 4.3):
        A_μ(x) = (φ(r), 0, 0, α(r)) · e^{-iωt}
        J_μ(x) = (ρ(r), 0, 0, β(r)) · e^{-iωt}
     Four INDEPENDENT radial functions. NOT a locked ansatz.

  2. ω enters as −ω²·field in the radial EOM (v5 §5.2):
        "the e^{-iωt} factor contributes an effective term +ω²|field|²
         to the radial equations of motion"

  3. f(s) = g·s² with g=1.0625, NO m_J² mass term (v5 §2):
        "For the numerical work in [2] we take f(s) = gs²"

     The benchmark doc Werbos sent earlier (with f(s)=½m_J²s + λ/4·s²) is
     superseded by v5.

Q=0 in the v5 framework means Q_A=0 (the EM charge from Lorenz constraint
on A), NOT Q_CS=0. Q_A and Q_J are TWO INDEPENDENT conserved charges
(v5 §8). The neutral chaoiton has Q_A=0 but Q_J≠0 — a WIMP profile.

ODE derivation (linearized far-field, then with cubic):
    □A_0 = J_0  →  Δ_r φ + ω²φ = ρ  →  Δ_r φ = ρ − ω²φ
    □A_3 = J_3  →  Δ_r α + ω²α = β  →  Δ_r α = β − ω²α
    □J_0 = A_0 − f'(J²)·J_0  →  Δ_r ρ = φ − ω²ρ − 2g·s·ρ
    □J_3 = A_3 − f'(J²)·J_3  →  Δ_r β = α − ω²β − 2g·s·β

where s = ρ² − β² (J·J with Minkowski signature, time-component minus
spatial; treated as real for the static radial reduction).

NOTE: this is a best-effort derivation from v5 §2 EL equations + §5.2 ω
prescription. The exact form of s and the cross-coupling deserves
confirmation from Werbos (would-be Q.v5 follow-up question).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_bvp

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

G_COUPLING = 1.0625          # v5 §13b
OMEGA = 1.0                  # electron calibration
R_MIN = 0.02
R_MAX_DEFAULT = 30.0
N_GRID_DEFAULT = 200
HBAR_C_MEVFM = 197.32698
M_E_MEV = 0.51100
H_CODE_ELECTRON = 0.494
R_PHYS_FM = 191.0


# ---------------------------------------------------------------------------
# ODE — v5-corrected 4-function system
# ---------------------------------------------------------------------------

def ode_rhs(r, y, g, omega, laplacian_dim=2):
    """8-state vector [φ, φ', α, α', ρ, ρ', β, β'].

    Equations (best-effort derivation from v5 paper §2, §5.2):
        Δ_r φ = ρ − ω²φ
        Δ_r α = β − ω²α
        Δ_r ρ = φ − ω²ρ − 2g·s·ρ
        Δ_r β = α − ω²β − 2g·s·β
    where s = ρ² − β² (J·J under Minkowski +--- signature)
    and Δ_r f = f''(r) + (D-1)/r · f'(r)   (D=2 thin-torus reduction)
    """
    phi   = y[0]; phi_p = y[1]
    alpha = y[2]; alp_p = y[3]
    rho   = y[4]; rho_p = y[5]
    beta  = y[6]; bet_p = y[7]

    inv_r_f = (laplacian_dim - 1) / r
    s = rho * rho - beta * beta
    w2 = omega * omega

    phi_pp = (rho - w2 * phi)            - inv_r_f * phi_p
    alp_pp = (beta - w2 * alpha)         - inv_r_f * alp_p
    rho_pp = (phi - w2 * rho - 2*g*s*rho) - inv_r_f * rho_p
    bet_pp = (alpha - w2 * beta - 2*g*s*beta) - inv_r_f * bet_p

    return np.vstack([phi_p, phi_pp, alp_p, alp_pp, rho_p, rho_pp, bet_p, bet_pp])


def bc_residuals(ya, yb):
    """Regularity at r_min + decay at r_max."""
    return np.array([
        ya[1],  # φ'(r_min) = 0
        ya[3],  # α'(r_min) = 0
        ya[5],  # ρ'(r_min) = 0
        ya[7],  # β'(r_min) = 0
        yb[0],  # φ(r_max) = 0
        yb[2],  # α(r_max) = 0
        yb[4],  # ρ(r_max) = 0
        yb[6],  # β(r_max) = 0
    ])


# ---------------------------------------------------------------------------
# Initial guess
# ---------------------------------------------------------------------------

def initial_guess(r, phi0, alpha0, rho0, beta0, decay_length):
    """Gaussian envelopes, independent amplitudes (NOT locked)."""
    r_shift = r - r[0]
    env = np.exp(-r_shift / decay_length)
    d_env = -env / decay_length
    return np.vstack([
        phi0 * env,   phi0 * d_env,
        alpha0 * env, alpha0 * d_env,
        rho0 * env,   rho0 * d_env,
        beta0 * env,  beta0 * d_env,
    ])


# ---------------------------------------------------------------------------
# Solve + post-process
# ---------------------------------------------------------------------------

def solve_one(g, omega, phi0, alpha0, rho0, beta0,
              decay_length=5.0, r_min=R_MIN, r_max=R_MAX_DEFAULT,
              n_grid=N_GRID_DEFAULT, tol=1e-6, max_nodes=20000):
    r = np.linspace(r_min, r_max, n_grid)
    y_init = initial_guess(r, phi0, alpha0, rho0, beta0, decay_length)
    try:
        sol = solve_bvp(
            fun=lambda r, y: ode_rhs(r, y, g, omega),
            bc=bc_residuals, x=r, y=y_init,
            tol=tol, max_nodes=max_nodes, verbose=0,
        )
    except Exception as e:
        return {"success": False, "error": str(e)[:120]}
    return {"success": bool(sol.success), "sol": sol,
            "message": sol.message, "n_iter": int(getattr(sol, "niter", -1))}


def post_process(sol_obj, g, omega, n_eval=2000):
    """Extract energy + Q_A + Q_J + diagnostics."""
    if sol_obj is None or not sol_obj.success:
        return None
    sol = sol_obj
    r = np.linspace(sol.x[0], sol.x[-1], n_eval)
    y = sol.sol(r)
    phi, phi_p   = y[0], y[1]
    alpha, alp_p = y[2], y[3]
    rho, rho_p   = y[4], y[5]
    beta, bet_p  = y[6], y[7]

    s = rho * rho - beta * beta
    w2 = omega * omega

    # Energy density — derived from L = -F·F - G·G + J·A - f(J²)
    # for the 4-function ansatz. Treating as time-averaged Hamiltonian.
    # Provisional form; exact normalization deserves Werbos confirmation.
    F_F = phi_p**2 + alp_p**2 + w2*phi**2 + w2*alpha**2     # |dA|² contributions
    G_G = rho_p**2 + bet_p**2 + w2*rho**2  + w2*beta**2     # |dJ|² contributions
    J_A = phi*rho - alpha*beta                              # J·A with metric sign
    f_JJ = g * s * s                                        # f(J²) = g·(J²)²
    density = F_F + G_G - J_A + f_JJ
    H_code = float((2.0 * np.pi)**2 * np.trapezoid(r * density, r))

    # Q_A: Noether charge from U(1) phase on A. For ψ=φ·e^{-iωt}:
    #   j_A ~ 2ω · (φ² − α²)   (time-component minus spatial in Minkowski)
    # Q_A = ∫ j_A d³x  ∝  ω · ∫ r·(φ² − α²)·dr   (2D radial)
    Q_A = float(2.0 * np.pi * omega * np.trapezoid(r * (phi*phi - alpha*alpha), r))
    Q_J = float(2.0 * np.pi * omega * np.trapezoid(r * (rho*rho - beta*beta), r))

    # Field max amplitudes (for trivial-solution detection)
    max_phi = float(np.max(np.abs(phi)))
    max_alp = float(np.max(np.abs(alpha)))
    max_rho = float(np.max(np.abs(rho)))
    max_bet = float(np.max(np.abs(beta)))
    is_trivial = max(max_phi, max_alp, max_rho, max_bet) < 1e-6

    m_chi_mev = H_code * M_E_MEV / H_CODE_ELECTRON

    return {
        "r": r, "phi": phi, "alpha": alpha, "rho": rho, "beta": beta,
        "H_code": H_code, "m_chi_MeV": m_chi_mev,
        "Q_A": Q_A, "Q_J": Q_J,
        "max_phi": max_phi, "max_alpha": max_alp,
        "max_rho": max_rho, "max_beta": max_bet,
        "is_trivial": is_trivial,
    }


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def run_sweep():
    print(f"\n{'='*112}")
    print("M6 neutral chaoiton — v5-paper-corrected sweep")
    print(f"{'='*112}\n")
    print(f"  g = {G_COUPLING},  ω = {OMEGA},  f(s) = g·s²  (no m_J² mass term)")
    print(f"  4-function asymmetric ansatz: (φ, α, ρ, β) all independent")
    print(f"  Target: Q_A ≈ 0 (electrically neutral)  AND  Q_J ≠ 0 (nuclear)\n")
    print(f"{'φ0':>5} {'α0':>5} {'ρ0':>5} {'β0':>5} {'OK':>3} {'triv':>4} "
          f"{'maxφ':>7} {'maxα':>7} {'maxρ':>7} {'maxβ':>7} "
          f"{'Q_A':>9} {'Q_J':>9} {'H':>10} {'mMeV':>9}")
    print("-" * 112)

    rows = []
    AMPL_VALUES = (0.0, 0.1, 0.3, 0.5, 1.0)

    # Seek Q_A ≈ 0 configurations: bias toward |φ| ≈ |α| (so j_A ∝ φ²-α² ≈ 0)
    # while having nonzero ρ, β for nontrivial Q_J.
    for phi0 in AMPL_VALUES:
        for alpha0 in AMPL_VALUES:
            for rho0 in AMPL_VALUES:
                for beta0 in AMPL_VALUES:
                    # Skip all-zero
                    if phi0 == 0 and alpha0 == 0 and rho0 == 0 and beta0 == 0:
                        continue
                    res = solve_one(G_COUPLING, OMEGA,
                                    phi0, alpha0, rho0, beta0,
                                    decay_length=5.0)
                    if not res["success"]:
                        continue
                    pp = post_process(res["sol"], G_COUPLING, OMEGA)
                    if pp is None:
                        continue
                    if pp["is_trivial"]:
                        continue  # only report non-trivial
                    triv = "Y" if pp["is_trivial"] else "n"
                    print(f"{phi0:>5.2f} {alpha0:>5.2f} {rho0:>5.2f} {beta0:>5.2f} "
                          f"{'Y':>3} {triv:>4} "
                          f"{pp['max_phi']:>7.3f} {pp['max_alpha']:>7.3f} "
                          f"{pp['max_rho']:>7.3f} {pp['max_beta']:>7.3f} "
                          f"{pp['Q_A']:>9.3f} {pp['Q_J']:>9.3f} "
                          f"{pp['H_code']:>10.3f} {pp['m_chi_MeV']:>9.3f}")
                    rows.append({
                        "phi0": phi0, "alpha0": alpha0, "rho0": rho0, "beta0": beta0,
                        "H_code": pp["H_code"], "m_chi_MeV": pp["m_chi_MeV"],
                        "Q_A": pp["Q_A"], "Q_J": pp["Q_J"],
                        "max_phi": pp["max_phi"], "max_alpha": pp["max_alpha"],
                        "max_rho": pp["max_rho"], "max_beta": pp["max_beta"],
                    })

    print(f"{'='*112}\n")
    print(f"Non-trivial converged: {len(rows)} rows")

    # Find Q_A ≈ 0 candidates (neutral chaoitons)
    neutral = [r for r in rows if abs(r["Q_A"]) < 0.01]
    print(f"Q_A ≈ 0 candidates (|Q_A| < 0.01): {len(neutral)}")
    if neutral:
        print("\nNEUTRAL CHAOITON CANDIDATES (sorted by lowest H):")
        best = sorted(neutral, key=lambda x: x["H_code"])
        for r in best[:10]:
            print(f"  φ0={r['phi0']:.2f}, α0={r['alpha0']:.2f}, ρ0={r['rho0']:.2f}, β0={r['beta0']:.2f}  "
                  f"→ Q_A={r['Q_A']:.4f}, Q_J={r['Q_J']:.4f}, "
                  f"H={r['H_code']:.4f}, m_χ={r['m_chi_MeV']:.4f} MeV")

    out = Path(__file__).resolve().parent.parent / "plots" / "m6_neutral_v5_sweep.json"
    out.parent.mkdir(exist_ok=True)
    with open(out, "w") as f:
        json.dump(rows, f, indent=2, default=str)
    print(f"\n[json] {out}")

    return rows, neutral


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="mode", required=False)
    p_single = sub.add_parser("single")
    p_single.add_argument("--phi0", type=float, default=0.1)
    p_single.add_argument("--alpha0", type=float, default=0.1)
    p_single.add_argument("--rho0", type=float, default=0.1)
    p_single.add_argument("--beta0", type=float, default=0.1)
    p_single.add_argument("--decay-length", type=float, default=5.0)
    p_single.add_argument("--r-max", type=float, default=R_MAX_DEFAULT)
    p_sweep = sub.add_parser("sweep")
    args = parser.parse_args()

    if args.mode == "single":
        res = solve_one(G_COUPLING, OMEGA,
                        args.phi0, args.alpha0, args.rho0, args.beta0,
                        decay_length=args.decay_length, r_max=args.r_max)
        if not res["success"]:
            print(f"FAIL: {res.get('error') or res.get('message')}")
            return 1
        pp = post_process(res["sol"], G_COUPLING, OMEGA)
        print(f"\nM6 neutral chaoiton — v5-corrected single run")
        print(f"  φ0={args.phi0}, α0={args.alpha0}, ρ0={args.rho0}, β0={args.beta0}")
        print(f"  max|φ|={pp['max_phi']:.4f},  max|α|={pp['max_alpha']:.4f}")
        print(f"  max|ρ|={pp['max_rho']:.4f},  max|β|={pp['max_beta']:.4f}")
        print(f"  Q_A = {pp['Q_A']:.6f}  (target ≈ 0)")
        print(f"  Q_J = {pp['Q_J']:.6f}  (target ≠ 0)")
        print(f"  H_code = {pp['H_code']:.6f}")
        print(f"  m_χ (MeV) = {pp['m_chi_MeV']:.6f}")
        print(f"  is_trivial = {pp['is_trivial']}")
    else:
        run_sweep()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
