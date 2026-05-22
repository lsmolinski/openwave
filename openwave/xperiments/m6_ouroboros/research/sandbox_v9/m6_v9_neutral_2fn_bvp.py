"""
v9 — Neutral chaoiton BVP, Interpretation B (2-function canonical, consistent with Q34).

Per Paul Werbos's 2026-05-21 PM reply (via DeepSeek):
  "apply the same `solve_bvp` method (Robin BCs at large r, integral constraint
   Q_CS=0) to the neutral sector. The ODEs are the same as for the charged
   case but with the topological charge fixed to zero..."

Interpretation B takes "same ODEs as the charged case" to mean the canonical
2-function (α, β) reduction in Sonnet's `ouroboros_benchmark.py`. For the
neutral chaoiton, α=0 (A-field absent), ω=0 (no time periodicity), reducing
to the single β-equation:

  β'' + β'/r - β/r² + λβ + 4gβ³ = 0   (same as Sonnet's neutral_ode)

We replace Sonnet's forward `solve_ivp` (slope BC β'(0)=B0, integrates outward
with growing oscillating tail at λ>0) with `solve_bvp` enforcing β(R_max)=0
(Dirichlet decay BC) or Robin β'(R_max) + k·β(R_max) = 0.

The free parameter is B0 (the linear slope of β at the origin). For each
(g, λ), the BVP solver finds the B0 that makes β(R_max)=0 — i.e., the
amplitude that produces a localized state with β tail vanishing at R_max.

If a true ground state exists, this gives:
  - a definite m_χ = H/Q × m_e (no integration-window dependence)
  - a clean far-field structure
  - reliable Yukawa coupling C from far-field amplitude
"""

import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
import warnings
warnings.filterwarnings('ignore')

# Physical constants (from Sonnet's benchmark)
HBAR_C   = 197.3269804   # MeV·fm
M_E_MEV  = 0.51100       # MeV
R_PHYS   = 191.0         # fm
M_SCALE  = HBAR_C / R_PHYS  # 1.033 MeV per natural unit

# Calibration point: electron-calibrated parameters
G_VAL   = 1.0000
LAM_VAL = 1.0

# Grid
R_MIN = 0.02
R_MAX = 12.0
N_GRID = 200


def neutral_ode_bvp(r, y, p, g, lam):
    """
    State y = [β, β']; parameter p[0] = B0 (slope at origin).
    ODE: β'' + β'/r - β/r² + λβ + 4gβ³ = 0
    Returns dy/dr as (2, n) array.
    """
    b, db = y
    d2b = -db/r + b/r**2 - lam*b - 4*g*b**3
    return np.vstack((db, d2b))


def neutral_bc_dirichlet(ya, yb, p):
    """
    Dirichlet decay BC: β(R_max) = 0.
    3 BCs total (2 states + 1 free parameter).
    """
    B0 = p[0]
    return np.array([
        ya[0] - B0 * R_MIN,   # β(R_MIN) = B0 · R_MIN  (linear-from-origin)
        ya[1] - B0,           # β'(R_MIN) = B0
        yb[0]                  # β(R_max) = 0  (Dirichlet decay)
    ])


def neutral_bc_robin(ya, yb, p, k_decay):
    """
    Robin decay BC: β'(R_max) + k·β(R_max) = 0.
    For exponential decay β(r) ~ exp(-k·r): β'/β = -k, so β' + k·β = 0.
    """
    B0 = p[0]
    return np.array([
        ya[0] - B0 * R_MIN,
        ya[1] - B0,
        yb[1] + k_decay * yb[0]   # β'(R_max) + k·β(R_max) = 0
    ])


def compute_observables(r, b, g, lam):
    """H/Q × m_e on Sonnet's cylindrical r·dr integration."""
    dr = np.diff(r)
    rm = 0.5*(r[:-1] + r[1:])
    bm = 0.5*(b[:-1] + b[1:])
    dbm = np.diff(b)/dr

    Q_J = np.sum(bm**2 * rm * dr)
    H = np.sum((dbm**2 + bm**2/rm**2 + lam*bm**2 + 4*g*bm**4) * rm * dr)
    HQ = H / Q_J if Q_J > 1e-15 else float('nan')

    r_rms = np.sqrt(np.sum(bm**2 * rm**3 * dr) / Q_J) * R_PHYS if Q_J > 1e-15 else float('nan')

    tail = abs(b[-1])
    peak = max(abs(b))
    peak_idx = np.argmax(abs(b))

    return {
        'H': float(H),
        'Q_J': float(Q_J),
        'HQ': float(HQ),
        'm_chi_MeV': float(HQ * M_E_MEV),
        'r_rms_fm': float(r_rms),
        'tail': float(tail),
        'peak': float(peak),
        'peak_r': float(r[peak_idx]),
        'tail_over_peak': float(tail/peak) if peak > 1e-12 else float('inf'),
    }


def run_dirichlet_bvp(B0_init, g=G_VAL, lam=LAM_VAL, r_max=R_MAX, n_grid=N_GRID):
    """Run BVP with Dirichlet decay BC β(R_max)=0."""
    r_mesh = np.linspace(R_MIN, r_max, n_grid)

    # Initial guess: localized profile shape
    b_init = B0_init * r_mesh * np.exp(-r_mesh)
    db_init = np.gradient(b_init, r_mesh)
    y_init = np.vstack((b_init, db_init))

    sol = solve_bvp(
        lambda r, y, p: neutral_ode_bvp(r, y, p, g, lam),
        neutral_bc_dirichlet,
        r_mesh,
        y_init,
        p=[B0_init],
        max_nodes=50000,
        tol=1e-5,
    )

    return sol


def run_robin_bvp(B0_init, k_decay, g=G_VAL, lam=LAM_VAL, r_max=R_MAX, n_grid=N_GRID):
    """Run BVP with Robin decay BC β'(R_max) + k·β(R_max) = 0."""
    r_mesh = np.linspace(R_MIN, r_max, n_grid)

    b_init = B0_init * r_mesh * np.exp(-r_mesh)
    db_init = np.gradient(b_init, r_mesh)
    y_init = np.vstack((b_init, db_init))

    sol = solve_bvp(
        lambda r, y, p: neutral_ode_bvp(r, y, p, g, lam),
        lambda ya, yb, p: neutral_bc_robin(ya, yb, p, k_decay),
        r_mesh,
        y_init,
        p=[B0_init],
        max_nodes=50000,
        tol=1e-5,
    )

    return sol


def report_solution(sol, label, g, lam):
    """Print converged solution diagnostics."""
    print(f"\n--- {label} ---")
    print(f"solve_bvp.status = {sol.status} ({sol.message})")
    print(f"solve_bvp.success = {sol.success}")
    if not sol.success:
        print(f"  CONVERGENCE FAILED — skipping observables")
        return None
    r = sol.x
    b = sol.y[0]
    B0_converged = sol.p[0]
    obs = compute_observables(r, b, g, lam)

    print(f"Final mesh: {len(r)} points on [{r[0]:.3f}, {r[-1]:.3f}]")
    print(f"Converged B0 = {B0_converged:+.6e}")
    print(f"β(R_MIN) = {b[0]:+.6e}    β'(R_MIN) ≈ {(b[1]-b[0])/(r[1]-r[0]):+.6e}")
    print(f"β(R_MAX) = {b[-1]:+.6e}   (should be ~0 for Dirichlet)")
    print(f"Peak:    β = {obs['peak']:+.6e}  at r = {obs['peak_r']:.3f}")
    print(f"tail/peak ratio: {obs['tail_over_peak']:.4e}  (lower = better localization)")
    print(f"H/Q = {obs['HQ']:.6f}")
    print(f"m_chi = {obs['m_chi_MeV']:.6f} MeV")
    print(f"r_rms = {obs['r_rms_fm']:.2f} fm")

    # β profile sample
    print(f"\n  β profile sample:")
    print(f"  {'r':>6}  {'β(r)':>14}")
    for rs in [0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 10.0, R_MAX*0.99]:
        if rs <= R_MAX:
            idx = np.argmin(np.abs(r - rs))
            print(f"  {r[idx]:>6.3f}  {b[idx]:>+14.6e}")

    # Sign-change count
    sign_changes = np.sum(np.diff(np.sign(b)) != 0)
    print(f"\nSign changes: {sign_changes}  (ground state = 0 nodes)")

    return obs


def main():
    print("=" * 70)
    print(f"v9 Interpretation B — Neutral 2-function BVP")
    print(f"(g={G_VAL}, λ={LAM_VAL}, R_MIN={R_MIN}, R_MAX={R_MAX})")
    print("=" * 70)
    print(f"ODE: β'' + β'/r - β/r² + λβ + 4gβ³ = 0")
    print(f"     (same as Sonnet's neutral_ode in `ouroboros_benchmark.py`)")
    print()

    # === Dirichlet variant ===
    print("\n" + "=" * 70)
    print("DIRICHLET DECAY BC — β(R_max) = 0")
    print("=" * 70)
    for B0_init in [0.001, 0.01, 0.05, 0.1, 0.3]:
        sol = run_dirichlet_bvp(B0_init)
        report_solution(sol, f"Dirichlet, B0_init={B0_init}", G_VAL, LAM_VAL)

    # === Robin variant ===
    print("\n" + "=" * 70)
    print("ROBIN DECAY BC — β'(R_max) + k·β(R_max) = 0")
    print("=" * 70)
    for k_decay in [0.5, 1.0, 2.0]:
        print(f"\n>>> k_decay = {k_decay}")
        for B0_init in [0.01, 0.05]:
            sol = run_robin_bvp(B0_init, k_decay)
            report_solution(sol, f"Robin k={k_decay}, B0_init={B0_init}", G_VAL, LAM_VAL)

    # === Negative λ scan — does TRUE exponential decay produce a soliton? ===
    print("\n" + "=" * 70)
    print("NEGATIVE λ SCAN — try λ < 0 where linear far-field has K_1 decay")
    print("=" * 70)
    print(f"  Linear far-field: β'' + λβ ≈ 0 → if λ < 0, β ~ exp(-√(-λ)·r)")
    print(f"  For λ = -1: decay rate √1 = 1 in natural units (m_J = 1.033 MeV)")
    print()

    for lam_neg in [-0.1, -0.5, -1.0, -2.0]:
        print(f"\n>>> λ = {lam_neg}  (linear K_1 decay rate = √(-λ) = {np.sqrt(-lam_neg):.3f})")
        for B0_init in [0.05, 0.1, 0.3, 0.5]:
            sol = run_dirichlet_bvp(B0_init, g=G_VAL, lam=lam_neg)
            obs = report_solution(sol, f"λ={lam_neg}, B0_init={B0_init} (Dirichlet)", G_VAL, lam_neg)

    # === Sanity: compare to Sonnet's IVP at same (g, λ, B0) ===
    print("\n" + "=" * 70)
    print("SANITY: Sonnet's IVP at (g={}, λ={}, B0=0.001) for comparison".format(G_VAL, LAM_VAL))
    print("=" * 70)

    from scipy.integrate import solve_ivp
    def ode_ivp(r, y, g, lam):
        b, db = y
        d2b = -db/r + b/r**2 - lam*b - 4*g*b**3
        return [db, d2b]

    B0_ivp = 0.001
    r_eval = np.linspace(R_MIN, R_MAX, 800)
    y0 = [B0_ivp * R_MIN, B0_ivp]
    sol_ivp = solve_ivp(
        ode_ivp, [R_MIN, R_MAX], y0,
        args=(G_VAL, LAM_VAL),
        t_eval=r_eval,
        method='RK45', rtol=1e-9, atol=1e-11, max_step=0.02
    )
    if sol_ivp.success:
        obs_ivp = compute_observables(sol_ivp.t, sol_ivp.y[0], G_VAL, LAM_VAL)
        print(f"IVP: m_chi = {obs_ivp['m_chi_MeV']:.6f} MeV  (Sonnet's published 0.998 at λ=1)")
        print(f"     tail/peak = {obs_ivp['tail_over_peak']:.4e}  (compare to BVP localization)")


if __name__ == '__main__':
    main()
