#!/usr/bin/env python3
"""
neutral_bvp_solver_fixed_BC.py

Modification of DeepSeek's 2026-05-22 template `neutral_bvp_solver.py`.

DeepSeek's verbatim script has β'(R_MIN)=0 as the origin BC, but for
3D spherical l=1 (p-wave) the regular solution behaves as β(r) ~ B0·r
near origin, so β(0) = 0 (not β'(0) = 0). With the verbatim BCs, every
B0_init converges to the trivial β ≡ 0 solution.

This modification:
  - Origin BC: β(R_MIN) ≈ 0 + slope B0 free parameter
  - Far-field BC: Robin (unchanged)
  - Free parameter: B0 (slope at origin)

3D spherical l=1:
  β'' + (2/r)β' - (2/r²)β - m_J²β + 4gβ³ = 0

Near origin (r small): the linear regular solution is β ~ B0·r.
  β = B0·r → β' = B0 → β'' = 0
  → 0 + 2·B0/r - 2·B0·r/r² - m_J²·B0·r + 4g·B0³·r³
  = 2·B0/r - 2·B0/r - m_J²·B0·r + ... = -m_J²·B0·r + ...
  → consistent with β ~ B0·r at leading order

Far field (large r): β ~ A·exp(-m_J·r)/r  (Yukawa/k_1 form for l=1, 3D)
  Robin BC: β'(R_max) + m_J·β(R_max) = 0 captures this.
"""

import numpy as np
from scipy.integrate import solve_bvp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Parameters
g = 1.0
m_J = 1.0
Rmin = 0.02      # small but not zero (avoid 1/r² singularity)
Rmax = 20.0
N = 500


def ode(r, y, p):
    """y = [β, β']; p[0] = B0 (slope at origin, free parameter)."""
    beta, betap = y
    betapp = -(2.0/r)*betap + (2.0/(r*r))*beta + (m_J*m_J)*beta - 4.0*g*beta**3
    return np.vstack((betap, betapp))


def bc(ya, yb, p):
    """
    Origin (R_MIN): β(R_MIN) = B0 · R_MIN (l=1 regular behavior β ~ B0·r)
    Far-field (R_max): Robin β' + m_J·β = 0 (k_1 decay)
    Free param: B0 (slope at origin)
    """
    B0 = p[0]
    return np.array([
        ya[0] - B0 * Rmin,        # β(R_MIN) = B0 · R_MIN
        ya[1] - B0,               # β'(R_MIN) = B0
        yb[1] + m_J * yb[0],      # Robin decay
    ])


def compute_observables(r, b):
    """Sonnet-style cylindrical r·dr or 3D-spherical r²·dr?"""
    dr = np.diff(r)
    rm = 0.5*(r[:-1] + r[1:])
    bm = 0.5*(b[:-1] + b[1:])
    dbm = np.diff(b)/dr

    # 3D spherical integration (per Q44 fix)
    Q_J_3D = np.sum(bm**2 * rm**2 * dr)
    H_3D = np.sum(
        (dbm**2 + 2*bm**2/rm**2 + m_J**2*bm**2 + 4*g*bm**4) * rm**2 * dr
    )
    HQ = H_3D / Q_J_3D if Q_J_3D > 1e-15 else float('nan')

    tail = abs(b[-1])
    peak_idx = np.argmax(abs(b))
    peak = abs(b[peak_idx])
    peak_r = r[peak_idx]
    sign_changes = np.sum(np.diff(np.sign(b)) != 0)

    return {
        'H': float(H_3D),
        'Q_J': float(Q_J_3D),
        'HQ': float(HQ),
        'tail': float(tail),
        'peak': float(peak),
        'peak_r': float(peak_r),
        'sign_changes': int(sign_changes),
        'tail_over_peak': float(tail/peak) if peak > 1e-12 else float('inf'),
    }


def run_bvp(B0_init):
    r_mesh = np.linspace(Rmin, Rmax, N)
    beta_guess = B0_init * r_mesh * np.exp(-m_J * r_mesh)
    betap_guess = B0_init * (1.0 - m_J * r_mesh) * np.exp(-m_J * r_mesh)
    y_guess = np.vstack((beta_guess, betap_guess))

    sol = solve_bvp(
        ode, bc, r_mesh, y_guess,
        p=[B0_init],
        tol=1e-6, max_nodes=20000,
    )
    return sol


def main():
    print("=" * 70)
    print("3D spherical l=1 BVP — DeepSeek's ODE + fixed l=1 BC + B0 free param")
    print(f"  g={g}, m_J={m_J}, R ∈ [{Rmin}, {Rmax}]")
    print("=" * 70)

    B0_list = [0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0]
    for B0_init in B0_list:
        sol = run_bvp(B0_init)
        if not sol.success:
            print(f"\nB0_init={B0_init}: FAILED — {sol.message}")
            continue

        B0_conv = sol.p[0]
        r = sol.x
        b = sol.y[0]
        obs = compute_observables(r, b)

        print(f"\nB0_init={B0_init}")
        print(f"  Converged B0 = {B0_conv:+.6e}")
        print(f"  β(R_MIN) = {b[0]:+.6e}    β'(R_MIN) ≈ {(b[1]-b[0])/(r[1]-r[0]):+.6e}")
        print(f"  Peak β   = {obs['peak']:+.6e} at r = {obs['peak_r']:.3f}")
        print(f"  β(R_MAX) = {b[-1]:+.6e}")
        print(f"  tail/peak = {obs['tail_over_peak']:.4e}  sign changes = {obs['sign_changes']}")
        print(f"  H/Q = {obs['HQ']:.6f}  Q_J = {obs['Q_J']:.6e}  H = {obs['H']:.6e}")

        # β profile sample
        print(f"  β sample:")
        for rs in [0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 19.0]:
            idx = np.argmin(np.abs(r - rs))
            print(f"    r={r[idx]:6.2f}  β={b[idx]:+.6e}")


if __name__ == '__main__':
    main()
