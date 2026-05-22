#!/usr/bin/env python3
"""
neutral_bvp_solver_mJ_free.py

Third variant of DeepSeek's 2026-05-22 template.

DeepSeek noted: "If you prefer to keep m_J as a free eigenvalue, we can
modify the script to treat it as an unknown parameter — just say the word."

Modification:
  - B0 (amplitude at origin) FIXED to a specified value
  - m_J FREE (eigenvalue — solver finds the value that supports the soliton)
  - l=1 BC: β(R_MIN) = B0·R_MIN + β'(R_MIN) = B0
  - Robin decay at R_max: β'(R_max) + m_J·β(R_max) = 0

For a 3D spherical l=1 cubic NLS (subcritical), there should be a discrete
set of m_J eigenvalues for each fixed B0 (or equivalently, a discrete set
of B0 for each fixed m_J). The lowest m_J gives the ground state.

3D spherical l=1 ODE (DeepSeek's confirmed form per Q43+Q44):
  β'' + (2/r)β' - (2/r²)β - m_J²β + 4gβ³ = 0
"""

import numpy as np
from scipy.integrate import solve_bvp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Parameters
g = 1.0
Rmin = 0.02
Rmax = 20.0
N = 500


def ode(r, y, p):
    """y = [β, β']; p[0] = m_J (free eigenvalue)."""
    beta, betap = y
    m_J = p[0]
    betapp = -(2.0/r)*betap + (2.0/(r*r))*beta + (m_J*m_J)*beta - 4.0*g*beta**3
    return np.vstack((betap, betapp))


def make_bc(B0):
    """BCs at R_MIN (l=1 regular) + Robin at R_max. m_J = p[0] free."""
    def bc(ya, yb, p):
        m_J = p[0]
        return np.array([
            ya[0] - B0 * Rmin,         # β(R_MIN) = B0 · R_MIN
            ya[1] - B0,                # β'(R_MIN) = B0  (l=1 slope)
            yb[1] + m_J * yb[0],       # Robin decay at R_max
        ])
    return bc


def compute_observables(r, b, m_J):
    """3D spherical r²·dr integration (per Q44 fix)."""
    dr = np.diff(r)
    rm = 0.5*(r[:-1] + r[1:])
    bm = 0.5*(b[:-1] + b[1:])
    dbm = np.diff(b)/dr

    Q_J = np.sum(bm**2 * rm**2 * dr)
    H = np.sum(
        (dbm**2 + 2*bm**2/rm**2 + m_J**2*bm**2 + 4*g*bm**4) * rm**2 * dr
    )
    HQ = H / Q_J if Q_J > 1e-15 else float('nan')

    tail = abs(b[-1])
    peak_idx = np.argmax(abs(b))
    peak = abs(b[peak_idx])
    peak_r = r[peak_idx]
    sign_changes = np.sum(np.diff(np.sign(b)) != 0)

    return {
        'H': float(H), 'Q_J': float(Q_J), 'HQ': float(HQ),
        'tail': float(tail), 'peak': float(peak), 'peak_r': float(peak_r),
        'sign_changes': int(sign_changes),
        'tail_over_peak': float(tail/peak) if peak > 1e-12 else float('inf'),
    }


def run_bvp(B0, m_J_init):
    r_mesh = np.linspace(Rmin, Rmax, N)
    beta_guess = B0 * r_mesh * np.exp(-m_J_init * r_mesh)
    betap_guess = B0 * (1.0 - m_J_init * r_mesh) * np.exp(-m_J_init * r_mesh)
    y_guess = np.vstack((beta_guess, betap_guess))

    sol = solve_bvp(
        ode, make_bc(B0), r_mesh, y_guess,
        p=[m_J_init],
        tol=1e-6, max_nodes=20000,
    )
    return sol


def main():
    print("=" * 70)
    print("3D spherical l=1 BVP — B0 fixed, m_J free eigenvalue")
    print(f"  g={g}, R ∈ [{Rmin}, {Rmax}]")
    print("=" * 70)

    HBAR_C = 197.3269804  # MeV·fm

    # Scan over (B0 fixed, m_J init) pairs — tighter near the ground-state regime
    test_cases = [
        # (B0_fixed, m_J_init, label)
        (0.40, 1.0, "B0=0.40, m_J_init=1.0"),
        (0.45, 1.0, "B0=0.45, m_J_init=1.0"),
        (0.50, 1.0, "B0=0.50, m_J_init=1.0  ← previous ground state"),
        (0.55, 1.0, "B0=0.55, m_J_init=1.0"),
        (0.60, 1.0, "B0=0.60, m_J_init=1.0"),
        (0.70, 1.0, "B0=0.70, m_J_init=1.0"),
        (0.80, 1.0, "B0=0.80, m_J_init=1.0"),
    ]
    for B0, m_J_init, label in test_cases:
        sol = run_bvp(B0, m_J_init)
        print(f"\n--- {label} ---")
        if not sol.success:
            print(f"  FAILED: {sol.message}")
            continue

        m_J_conv = sol.p[0]
        r = sol.x
        b = sol.y[0]
        obs = compute_observables(r, b, m_J_conv)

        print(f"  Converged m_J = {m_J_conv:+.6f}  (physical: m_J·ℏc/R_phys would scale)")
        print(f"  β(R_MIN) = {b[0]:+.6e}    β'(R_MIN) ≈ {(b[1]-b[0])/(r[1]-r[0]):+.6e}")
        print(f"  Peak β   = {obs['peak']:+.6e} at r = {obs['peak_r']:.3f}")
        print(f"  β(R_MAX) = {b[-1]:+.6e}")
        print(f"  tail/peak = {obs['tail_over_peak']:.4e}  sign changes = {obs['sign_changes']}")
        print(f"  H/Q = {obs['HQ']:.6f}  Q_J = {obs['Q_J']:.6e}  H = {obs['H']:.6e}")
        print(f"  β sample:")
        for rs in [0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 19.0]:
            idx = np.argmin(np.abs(r - rs))
            print(f"    r={r[idx]:6.2f}  β={b[idx]:+.6e}")


if __name__ == '__main__':
    main()
