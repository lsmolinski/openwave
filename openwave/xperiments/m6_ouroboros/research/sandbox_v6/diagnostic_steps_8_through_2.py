"""Steps (8), (11), (4), (2) diagnostic runner — run AFTER v6.6 to test
candidate fixes against the converged field profile.

(8) H sign + kinetic-factor variants
(11) field profile inspection
(4)  (m_J², λ_bench) sweep
(2)  lepton-scan ratio-invariance trial
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
from scipy.integrate import solve_bvp

from m6_v6_4fn_calibrated_bvp import (
    R_MIN, R_MAX_DEFAULT, I_TARGET, TARGET_HQ,
    make_rhs, make_bc, initial_guess,
)


def solve_v6_at(omega_init=1.0, lambda_init=0.1, m_J_sq=0.5, lam_bench=1.0,
                V_norm=0.1, r_max=15.0, n_grid=500, max_nodes=100000, tol=5e-3,
                warm_start=False):
    """Re-solve the v6 BVP at given config. Returns (sol, ω, λ_LM) or None."""
    s = np.linspace(0, 1, n_grid)
    r_grid = R_MIN + (r_max - R_MIN) * (s ** 1.5)
    y0 = initial_guess(r_grid, warm_start=warm_start)
    p0 = np.array([omega_init, lambda_init])
    rhs = make_rhs(m_J_sq, lam_bench)
    bc = make_bc(m_J_sq, V_norm)
    try:
        sol = solve_bvp(rhs, bc, r_grid, y0, p=p0,
                        tol=tol, max_nodes=max_nodes, verbose=0)
        return sol
    except Exception as e:
        print(f"  EXCEPTION: {str(e)[:100]}")
        return None


def compute_H_variants(sol, omega, g=1.0625, m_J_sq=0.5, lam_bench=1.0):
    """Compute H under multiple sign + factor conventions for the same field
    profile. Returns dict: variant_label -> H value."""
    r = sol.x
    V, Vp, A, Ap, Q, Qp, J, Jp, I = sol.y

    # Core terms
    kin_half = 0.5 * (Vp**2 + Ap**2 + Qp**2 + Jp**2)        # (1/2)(V')² etc.
    kin_full = (Vp**2 + Ap**2 + Qp**2 + Jp**2)              # (V')² etc.
    omega_kin = 0.5 * omega**2 * (V**2 + A**2 + Q**2 + J**2)
    omega_kin_full = omega**2 * (V**2 + A**2 + Q**2 + J**2)
    cross_minus = -V * Q + A * J
    cross_plus = +V * Q - A * J  # SIGN FLIPPED
    quartic_deepseek = 0.25 * g * ((V**2 + Q**2)**2 + (A**2 + J**2)**2
                                   + 2.0 * (V * A - Q * J)**2)
    QQ_JJ = Q * Q - J * J
    quartic_benchmark = (0.5 * m_J_sq * QQ_JJ
                          + 0.25 * lam_bench * QQ_JJ * QQ_JJ)

    variants = {}

    # v6.6 reference (DeepSeek form as we implemented it)
    variants["v6.6_baseline"] = float(np.trapezoid(
        kin_half + omega_kin + cross_minus + quartic_deepseek, x=r))

    # ISSUE A only: flip cross-term sign
    variants["A_flip_cross"] = float(np.trapezoid(
        kin_half + omega_kin + cross_plus + quartic_deepseek, x=r))

    # ISSUE B only: full kinetic (factor 1 not 1/2)
    variants["B_kinetic_full"] = float(np.trapezoid(
        kin_full + omega_kin + cross_minus + quartic_deepseek, x=r))

    # ISSUE A+B: both fixes together
    variants["AB_both_fixes"] = float(np.trapezoid(
        kin_full + omega_kin + cross_plus + quartic_deepseek, x=r))

    # ω-kinetic doubled (ω² not (1/2)ω²)
    variants["omega_kin_full"] = float(np.trapezoid(
        kin_half + omega_kin_full + cross_minus + quartic_deepseek, x=r))

    # Drop quartic entirely (test if quartic is the residual)
    variants["no_quartic"] = float(np.trapezoid(
        kin_half + omega_kin + cross_minus, x=r))

    # Drop ω-kinetic (test if absorbed into m_eff² differently)
    variants["no_omega_kin"] = float(np.trapezoid(
        kin_half + cross_minus + quartic_deepseek, x=r))

    # Toroidal: add r weight on all terms (recovers full toroidal H)
    variants["toroidal_r_weight"] = float(np.trapezoid(
        r * (kin_half + omega_kin + cross_minus + quartic_deepseek), x=r))

    # Toroidal + 2π² prefactor (v5 form)
    variants["v5_full"] = float((2 * np.pi)**2 * np.trapezoid(
        r * (kin_full + omega_kin + cross_minus + quartic_benchmark), x=r))

    return variants


def step_11_profile_inspection(sol, label="v6.6"):
    """Step (11): print field profile at key r values."""
    print(f"\n{'='*80}")
    print(f"STEP (11): FIELD PROFILE INSPECTION — {label}")
    print(f"{'='*80}")
    r = sol.x
    V, Vp, A, Ap, Q, Qp, J, Jp, I = sol.y
    # Sample at standard radii
    sample_r = [0.05, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 15.0]
    print(f"\n{'r':>6} {'V':>10} {'A':>10} {'Q':>10} {'J':>10} {'I (Q_CS acc.)':>16}")
    print("-" * 80)
    for r_target in sample_r:
        if r_target > r[-1]:
            continue
        idx = np.argmin(np.abs(r - r_target))
        print(f"{r[idx]:>6.2f} {V[idx]:>+10.4f} {A[idx]:>+10.4f} "
              f"{Q[idx]:>+10.4f} {J[idx]:>+10.4f} {I[idx]:>+16.6f}")
    print(f"\nPeak amplitudes: V={np.abs(V).max():.4f}, A={np.abs(A).max():.4f}, "
          f"Q={np.abs(Q).max():.4f}, J={np.abs(J).max():.4f}")
    print(f"Zero crossings: V={int(np.sum(np.diff(np.sign(V[np.abs(V)>1e-6])) != 0))}, "
          f"A={int(np.sum(np.diff(np.sign(A[np.abs(A)>1e-6])) != 0))}, "
          f"Q={int(np.sum(np.diff(np.sign(Q[np.abs(Q)>1e-6])) != 0))}, "
          f"J={int(np.sum(np.diff(np.sign(J[np.abs(J)>1e-6])) != 0))}")


def step_8_h_variants(sol, omega, label="v6.6"):
    """Step (8): test H sign + kinetic-factor variants."""
    print(f"\n{'='*80}")
    print(f"STEP (8): H FUNCTIONAL VARIANTS — {label}")
    print(f"{'='*80}")
    Q_CS = float(sol.y[8, -1])  # I(R_max) directly = Q_CS in v6 convention
    variants = compute_H_variants(sol, omega)
    print(f"\nQ_CS (from I-state) = {Q_CS:.6f}")
    print(f"Target H/Q_CS       = {TARGET_HQ}")
    print(f"\n{'Variant':<22} {'H':>14} {'H/Q_CS':>10} {'% over 1.6969':>15}")
    print("-" * 70)
    for label_v, H in variants.items():
        ratio = H / Q_CS if abs(Q_CS) > 1e-14 else float('nan')
        pct = (ratio - TARGET_HQ) / TARGET_HQ * 100
        marker = " ★" if abs(pct) < 5 else ""
        print(f"{label_v:<22} {H:>14.4f} {ratio:>10.4f} {pct:>+14.2f}%{marker}")


def step_4_m_lambda_sweep(omega_init=1.0, lambda_init=0.1):
    """Step (4): (m_J², λ_bench) parameter sweep at v6 config."""
    print(f"\n{'='*80}")
    print(f"STEP (4): (m_J², λ_bench) SWEEP at v6 config")
    print(f"{'='*80}")
    m_J_sq_grid = [0.2, 0.5, 1.0, 2.0]
    lam_bench_grid = [0.5, 1.0, 2.0, 4.0]
    print(f"\nFor each (m_J², λ_bench), solve v6 BVP. Report ω, λ_LM, "
          f"H/Q_CS (DeepSeek form).\n")
    print(f"{'m_J²':>6} {'λ_bench':>9} {'status':>7} {'ω':>9} {'λ_LM':>9} "
          f"{'H/Q':>10} {'% off target':>13}")
    print("-" * 75)
    best = None
    for m_J_sq in m_J_sq_grid:
        for lam_bench in lam_bench_grid:
            sol = solve_v6_at(omega_init=omega_init, lambda_init=lambda_init,
                              m_J_sq=m_J_sq, lam_bench=lam_bench,
                              max_nodes=50000, tol=5e-3)
            if sol is None:
                print(f"{m_J_sq:>6.2f} {lam_bench:>9.2f}    FAIL")
                continue
            omega = float(sol.p[0])
            lam_LM = float(sol.p[1])
            r = sol.x
            V, Vp, A, Ap, Q, Qp, J, Jp, I = sol.y
            Q_CS = float(I[-1])
            variants = compute_H_variants(sol, omega, m_J_sq=m_J_sq,
                                           lam_bench=lam_bench)
            HQ = variants["v6.6_baseline"] / Q_CS if abs(Q_CS) > 1e-14 else float('nan')
            pct = (HQ - TARGET_HQ) / TARGET_HQ * 100
            print(f"{m_J_sq:>6.2f} {lam_bench:>9.2f} {sol.status:>7} "
                  f"{omega:>9.4f} {lam_LM:>9.3f} {HQ:>10.4f} {pct:>+12.2f}%")
            if abs(HQ - TARGET_HQ) / TARGET_HQ < (
                abs(best[2] - TARGET_HQ) / TARGET_HQ if best else 1.0
            ):
                best = (m_J_sq, lam_bench, HQ, omega, lam_LM, sol.status)
    if best:
        print(f"\nBest: m_J²={best[0]}, λ_bench={best[1]}, H/Q={best[2]:.4f} "
              f"(ω={best[3]:.4f}, λ_LM={best[4]:.3f}, status={best[5]})")


def step_2_lepton_scan(omega_targets=(1.0, 5.0, 12.0, 20.0, 40.7),
                       m_J_sq=0.5, lam_bench=1.0):
    """Step (2): lepton scan — fix omega at target values and look for the
    Q_CS=1 chaoiton. If the muon mode (ω≈12) gives H ratio ≈ 207, the
    absolute calibration gap doesn't affect mass ratios.

    Implementation: ω is a free eigenvalue in v6, so we vary the INITIAL
    GUESS to nudge the solver toward different ω. For each target ω_init,
    the solver may snap back to ω≈1 or settle near the initial value.
    """
    print(f"\n{'='*80}")
    print(f"STEP (2): LEPTON SCAN — vary ω init to find different stable modes")
    print(f"{'='*80}")
    print(f"\nm_J²={m_J_sq}, λ_bench={lam_bench}, varying ω initial guess:\n")
    print(f"{'ω_init':>7} {'final ω':>9} {'λ_LM':>9} {'status':>7} "
          f"{'H':>12} {'Q_CS':>9} {'H/Q':>10} {'rel. to ω≈1':>12}")
    print("-" * 90)
    base_H = None
    results = []
    for omega_init in omega_targets:
        sol = solve_v6_at(omega_init=omega_init, lambda_init=0.1,
                          m_J_sq=m_J_sq, lam_bench=lam_bench,
                          max_nodes=100000, tol=5e-3)
        if sol is None:
            print(f"{omega_init:>7.2f}    FAIL")
            continue
        omega = float(sol.p[0])
        lam_LM = float(sol.p[1])
        Q_CS = float(sol.y[8, -1])
        variants = compute_H_variants(sol, omega, m_J_sq=m_J_sq,
                                       lam_bench=lam_bench)
        H = variants["v6.6_baseline"]
        HQ = H / Q_CS if abs(Q_CS) > 1e-14 else float('nan')
        if base_H is None and abs(omega - 1.0) < 0.5:
            base_H = H
        rel = H / base_H if base_H else float('nan')
        results.append({"omega_init": omega_init, "omega": omega, "lam_LM": lam_LM,
                        "status": sol.status, "H": H, "Q_CS": Q_CS, "HQ": HQ,
                        "rel": rel})
        print(f"{omega_init:>7.2f} {omega:>9.4f} {lam_LM:>9.3f} "
              f"{sol.status:>7} {H:>12.4f} {Q_CS:>9.4f} {HQ:>10.4f} "
              f"{rel:>12.4f}")
    print(f"\nTarget mass ratios:")
    print(f"  m_μ/m_e = 207 (muon/electron)")
    print(f"  m_τ/m_e ≈ 3477 (tau/electron)")
    print(f"\nIf H ratios above match these, the 4.8% absolute gap is")
    print(f"irrelevant for ApJ deliverables (only ratios matter).")
    return results


def main():
    # Reproduce v6.6 best run
    print("=" * 80)
    print("REPRODUCING v6.6 BEST RUN")
    print("=" * 80)
    sol_66 = solve_v6_at(omega_init=1.0, lambda_init=0.1, r_max=15.0,
                          n_grid=500, max_nodes=100000, tol=5e-3)
    if sol_66 is None:
        print("v6.6 reproduction FAILED")
        return
    omega_66 = float(sol_66.p[0])
    lam_LM_66 = float(sol_66.p[1])
    Q_CS_66 = float(sol_66.y[8, -1])
    print(f"  Reproduced: ω={omega_66:.4f}, λ_LM={lam_LM_66:.3f}, "
          f"Q_CS={Q_CS_66:.4f}, status={sol_66.status}")

    # Step (8) — H functional variants
    step_8_h_variants(sol_66, omega_66, label="v6.6 reproduction")

    # Step (11) — field profile inspection
    step_11_profile_inspection(sol_66, label="v6.6 reproduction")

    # Step (4) — (m_J², λ_bench) sweep
    step_4_m_lambda_sweep()

    # Step (2) — lepton scan
    step_2_lepton_scan()


if __name__ == "__main__":
    main()
