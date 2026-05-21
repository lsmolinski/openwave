"""
Ouroboros Lagrangian — Canonical Numerical Benchmark
=====================================================
Authors: Paul J. Werbos, Claude (Anthropic AI, Claude Sonnet 4.6)
Zenodo companion code to DOIs: 10.5281/zenodo.20030162, 20218067, 20242421, 20298669

CANONICAL ODE SYSTEM (vector Laplacian form):
  α'' + α'/r - α/r² + ω²α = β(r)           [A-field]
  β'' + β'/r - β/r² + ω²β = α(r) - λβ - 4gβ³  [J-field]

BOUNDARY CONDITIONS:
  α(0) = 0, α'(0) = A0   (slope BC, NOT value BC)
  β(0) = 0, β'(0) = B0   (slope BC, NOT value BC)
  α, β → 0 as r → ∞     (localization)

NOTE: The correct Laplacian for the φ-component of a vector field
in cylindrical/toroidal coordinates is:
  Δ_vec f = f'' + f'/r - f/r²    (NOT f'' + 2f'/r which is scalar spherical)

NEUTRAL CHAOITON (Q_A = 0, spin-0):
  Set A0 = 0, ω = 0. The α-equation decouples (α=0 exact solution).
  Remaining single equation:
  β'' + β'/r - β/r² + λβ + 4gβ³ = 0
  with β(0)=0, β'(0)=B0.
"""

import numpy as np
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings('ignore')

# ── Physical constants ───────────────────────────────────────────
hbar_c   = 197.3269804   # MeV·fm
m_e_MeV  = 0.51100       # MeV
R_phys   = 191.0         # fm  (from electron calibration)
alpha_fs = 1/137.036
e_nat    = np.sqrt(4*np.pi*alpha_fs)  # = 0.30282 (dimensionless charge)
HQ_target = m_e_MeV / e_nat           # = 1.6875 (electron target)

# ── Integration parameters ───────────────────────────────────────
R_MAX   = 12.0
R_INNER = 0.02    # start away from r=0 singularity
N_GRID  = 800

r_grid = np.linspace(R_INNER, R_MAX, N_GRID)

# ════════════════════════════════════════════════════════════════
# CANONICAL ODE (vector Laplacian, slope BCs)
# ════════════════════════════════════════════════════════════════

def charged_odes(r, y, g, omega, lam):
    """
    Full 2-function ODE for charged chaoiton.
    State: y = [α, α', β, β']
    Vector Laplacian: f'' + f'/r - f/r²
    """
    a, da, b, db = y
    # α-equation: α'' + α'/r - α/r² + ω²α = β
    d2a = -da/r + a/r**2 - omega**2 * a + b
    # β-equation: β'' + β'/r - β/r² + ω²β = α - λβ - 4gβ³
    d2b = -db/r + b/r**2 - omega**2 * b + a - lam*b - 4*g*b**3
    return [da, d2a, db, d2b]

def neutral_ode(r, y, g, lam):
    """
    Single-function ODE for neutral chaoiton (α=0, ω=0).
    β'' + β'/r - β/r² + λβ + 4gβ³ = 0
    State: y = [β, β']
    """
    b, db = y
    d2b = -db/r + b/r**2 - lam*b - 4*g*b**3
    return [db, d2b]

# ════════════════════════════════════════════════════════════════
# SOLVERS
# ════════════════════════════════════════════════════════════════

def solve_charged(g, omega, lam, A0, B0,
                  r_max=R_MAX, r_inner=R_INNER, n_grid=N_GRID,
                  rtol=1e-9, atol=1e-11):
    """
    Solve charged chaoiton ODE.
    SLOPE boundary conditions: α'(r_inner)=A0, β'(r_inner)=B0
    Returns (r, alpha, beta) or None if integration fails.
    """
    r_eval = np.linspace(r_inner, r_max, n_grid)
    # Slope BCs: field ~ A0*r near origin → value = A0*r_inner, slope = A0
    y0 = [A0 * r_inner, A0, B0 * r_inner, B0]
    sol = solve_ivp(
        charged_odes, [r_inner, r_max], y0,
        args=(g, omega, lam),
        t_eval=r_eval,
        method='RK45', rtol=rtol, atol=atol, max_step=0.02
    )
    if not sol.success:
        return None
    return sol.t, sol.y[0], sol.y[2]

def solve_neutral(g, lam, B0,
                  r_max=R_MAX, r_inner=R_INNER, n_grid=N_GRID,
                  rtol=1e-9, atol=1e-11):
    """
    Solve neutral chaoiton ODE (α=0, ω=0).
    Returns (r, beta) or None if integration fails.
    """
    r_eval = np.linspace(r_inner, r_max, n_grid)
    y0 = [B0 * r_inner, B0]
    sol = solve_ivp(
        neutral_ode, [r_inner, r_max], y0,
        args=(g, lam),
        t_eval=r_eval,
        method='RK45', rtol=rtol, atol=atol, max_step=0.02
    )
    if not sol.success:
        return None
    return sol.t, sol.y[0]

# ════════════════════════════════════════════════════════════════
# OBSERVABLES
# ════════════════════════════════════════════════════════════════

def compute_charged_observables(r, a, b, g, omega):
    """Compute H, Q_A, Q_J, L, H/Q, L/Q for charged chaoiton."""
    dr = np.diff(r)
    rm = 0.5*(r[:-1] + r[1:])
    am = 0.5*(a[:-1] + a[1:])
    bm = 0.5*(b[:-1] + b[1:])
    dam = np.diff(a)/dr
    dbm = np.diff(b)/dr

    Q_A = np.sum(am**2 * rm * dr)
    Q_J = np.sum(bm**2 * rm * dr)
    if Q_J < 1e-12:
        return None

    H = np.sum((dam**2 + dbm**2
                + am**2/rm**2 + bm**2/rm**2
                + 4*g*bm**4) * rm * dr)
    L = omega * Q_J

    # RMS radius of J-field
    r_rms = np.sqrt(np.sum(bm**2 * rm**3 * dr) / Q_J) * R_phys

    return {
        'H': float(H), 'Q_A': float(Q_A), 'Q_J': float(Q_J),
        'L': float(L), 'HQ': float(H/Q_J), 'LQ': float(omega),
        'r_rms_fm': r_rms,
        'tail': float(abs(a[-1]) + abs(b[-1])),
        'mass_pred_MeV': float(H/Q_J * m_e_MeV),
    }

def compute_neutral_observables(r, b, g, lam):
    """Compute H, Q_J, H/Q for neutral chaoiton."""
    dr = np.diff(r)
    rm = 0.5*(r[:-1] + r[1:])
    bm = 0.5*(b[:-1] + b[1:])
    dbm = np.diff(b)/dr

    Q_J = np.sum(bm**2 * rm * dr)
    if Q_J < 1e-12:
        return None

    H = np.sum((dbm**2 + bm**2/rm**2
                + lam*bm**2 + 4*g*bm**4) * rm * dr)

    r_rms = np.sqrt(np.sum(bm**2 * rm**3 * dr) / Q_J) * R_phys

    return {
        'H': float(H), 'Q_J': float(Q_J), 'Q_A': 0.0,
        'HQ': float(H/Q_J),
        'r_rms_fm': r_rms,
        'tail': float(abs(b[-1])),
        'mass_pred_MeV': float(H/Q_J * m_e_MeV),
    }

def is_localized(obs, tail_thresh=0.15):
    """Check localization criterion."""
    return obs is not None and obs['tail'] < tail_thresh

# ════════════════════════════════════════════════════════════════
# ELECTRON CALIBRATION SCAN
# ════════════════════════════════════════════════════════════════

def electron_calibration_scan(
        g_vals=None, lam=1.0, omega=1.0,
        A0_vals=None, B0_vals=None):
    """
    Scan for electron calibration: target H/Q = 1.6875
    Best result: g=1.0625, lam=1.0, omega=1.0, A0=0.1, B0=0.1
    gives H/Q = 1.6969 (0.56% gap).
    """
    if g_vals is None:
        g_vals = np.linspace(0.8, 1.4, 25)
    if A0_vals is None:
        A0_vals = [0.05, 0.1, 0.2]
    if B0_vals is None:
        B0_vals = [0.05, 0.1, 0.2]

    print(f"Electron calibration scan: target H/Q = {HQ_target:.4f}")
    print(f"{'g':>8} {'A0':>6} {'B0':>6} {'H/Q':>9} {'gap%':>7} {'tail':>7}")
    print("─" * 50)

    best = None; best_gap = 1e10
    for g in g_vals:
        for A0 in A0_vals:
            for B0 in B0_vals:
                sol = solve_charged(g, omega, lam, A0, B0)
                if sol is None: continue
                r, a, b = sol
                obs = compute_charged_observables(r, a, b, g, omega)
                if not is_localized(obs): continue
                gap = abs(obs['HQ'] - HQ_target) / HQ_target * 100
                print(f"{g:>8.4f} {A0:>6.3f} {B0:>6.3f} "
                      f"{obs['HQ']:>9.4f} {gap:>7.3f}% {obs['tail']:>7.4f}")
                if gap < best_gap:
                    best_gap = gap
                    best = {'g': g, 'A0': A0, 'B0': B0, **obs}

    if best:
        print(f"\nBest: g={best['g']:.4f}, gap={best_gap:.3f}%")
        print(f"  H/Q = {best['HQ']:.6f} (target {HQ_target:.6f})")
        print(f"  L/Q = {best['LQ']:.3f} (2L/Q = {2*best['LQ']:.3f}, "
              f"electron g_e = 2.00232)")
    return best

# ════════════════════════════════════════════════════════════════
# NEUTRAL CHAOITON SCAN
# ════════════════════════════════════════════════════════════════

def neutral_chaoiton_scan(
        g_vals=None, lam_vals=None, B0_vals=None):
    """
    Scan for neutral chaoiton ground state.
    Mass prediction: m ≈ λ × 0.511 MeV
    At λ=1.0, lightest found: 0.508 MeV (Griesi 2026, independent).
    """
    if g_vals is None:
        g_vals = [0.5, 1.0, 1.0625, 2.0, 3.0, 5.0]
    if lam_vals is None:
        lam_vals = [0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0]
    if B0_vals is None:
        # Fine B0 scan including very small values for ground state search
        B0_vals = [0.001, 0.003, 0.005, 0.008, 0.01, 0.02,
                   0.05, 0.1, 0.2, 0.4, 0.6, 0.8]

    print("Neutral chaoiton scan (α=0, ω=0)")
    print(f"Mass prediction: m ≈ λ × {m_e_MeV} MeV")
    print(f"{'g':>6} {'λ':>6} {'B0':>8} {'H/Q':>9} {'mass(MeV)':>11} "
          f"{'m/λ×me':>9} {'r_rms fm':>10} {'tail':>7}")
    print("─" * 72)

    solutions = []
    for g in g_vals:
        for lam in lam_vals:
            for B0 in B0_vals:
                sol = solve_neutral(g, lam, B0)
                if sol is None: continue
                r, b = sol
                obs = compute_neutral_observables(r, b, g, lam)
                if not is_localized(obs, tail_thresh=0.12): continue
                # Mass ratio to λ×m_e
                m_ratio = obs['mass_pred_MeV'] / (lam * m_e_MeV)
                solutions.append({'g': g, 'lam': lam, 'B0': B0, **obs})
                print(f"{g:>6.3f} {lam:>6.2f} {B0:>8.4f} "
                      f"{obs['HQ']:>9.4f} {obs['mass_pred_MeV']:>11.4f} "
                      f"{m_ratio:>9.4f} {obs['r_rms_fm']:>10.1f} "
                      f"{obs['tail']:>7.4f}")

    print(f"\nTotal neutral solutions: {len(solutions)}")
    if solutions:
        lightest = min(solutions, key=lambda x: x['mass_pred_MeV'])
        print(f"Lightest at λ=1.0: "
              f"{[s for s in solutions if abs(s['lam']-1.0)<0.01]}")
        print(f"\nOverall lightest: g={lightest['g']:.3f}, "
              f"λ={lightest['lam']:.2f}, B0={lightest['B0']:.4f}")
        print(f"  mass = {lightest['mass_pred_MeV']:.4f} MeV")
        print(f"  = {lightest['mass_pred_MeV']/m_e_MeV:.3f} × m_e")
    return solutions

# ════════════════════════════════════════════════════════════════
# LEPTON SPECTRUM SCAN
# ════════════════════════════════════════════════════════════════

def lepton_spectrum_scan(g=1.0625, lam=1.0, A0=0.1, B0=0.1,
                          omega_vals=None):
    """
    Scan ω at fixed (g, λ) to map lepton mass spectrum.
    Results: electron (ω=1), muon (ω≈11, gap 4.3%), tau (ω≈41, gap 4.9%)
    Mass scaling: H ∝ ω^2.22
    """
    if omega_vals is None:
        omega_vals = list(np.linspace(1.0, 3.0, 5)) + \
                     list(np.linspace(3.0, 15.0, 12)) + \
                     list(np.linspace(15.0, 50.0, 15)) + \
                     list(np.linspace(50.0, 80.0, 6))
        omega_vals = sorted(set([round(o, 2) for o in omega_vals]))

    targets = {
        'electron': 0.511, 'muon': 105.658,
        'pion+': 139.570, 'tau': 1776.86
    }

    # Get electron reference
    e_sol = solve_charged(g, 1.0, lam, A0, B0)
    if e_sol is None:
        print("ERROR: electron reference solution not found")
        return
    r_e, a_e, b_e = e_sol
    obs_e = compute_charged_observables(r_e, a_e, b_e, g, 1.0)
    HQ_e = obs_e['HQ']
    print(f"Electron reference: H/Q = {HQ_e:.4f}, "
          f"R_phys = {R_phys} fm")
    print()
    print(f"{'ω':>7} {'H/Q':>9} {'pred MeV':>10} {'best match':>12} "
          f"{'gap%':>7} {'r_rms fm':>10}")
    print("─" * 62)

    results = []
    for omega in omega_vals:
        sol = solve_charged(g, omega, lam, A0, B0)
        if sol is None: continue
        r, a, b = sol
        obs = compute_charged_observables(r, a, b, g, omega)
        if not is_localized(obs): continue
        mass_pred = obs['HQ'] / HQ_e * m_e_MeV
        best = min(targets.items(), key=lambda x: abs(x[1]-mass_pred)/x[1])
        gap = abs(best[1]-mass_pred)/best[1]*100
        marker = " ★" if gap < 5 else (" ◆" if gap < 10 else "")
        print(f"{omega:>7.2f} {obs['HQ']:>9.4f} {mass_pred:>10.3f} "
              f"{best[0]:>12} {gap:>7.2f}%{marker} "
              f"{obs['r_rms_fm']:>8.1f}")
        results.append((omega, mass_pred, obs, best[0], gap))

    good = [(o, m, obs, n, g2) for o, m, obs, n, g2 in results if g2 < 8]
    print(f"\nMatches within 8%: {len(good)}")
    for o, m, obs, name, gap2 in good:
        print(f"  {name}: ω={o:.2f}, pred={m:.3f} MeV, gap={gap2:.2f}%")
    return results

# ════════════════════════════════════════════════════════════════
# QUICK DEMO
# ════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("=" * 62)
    print("OUROBOROS BENCHMARK — Quick demo")
    print("=" * 62)
    print()

    # 1. Electron reference point
    print("1. Electron calibration point:")
    sol = solve_charged(1.0625, 1.0, 1.0, 0.1, 0.1)
    if sol:
        r, a, b = sol
        obs = compute_charged_observables(r, a, b, 1.0625, 1.0)
        if obs:
            print(f"   g=1.0625, ω=1.0, λ=1.0, A0=0.1, B0=0.1")
            print(f"   H/Q = {obs['HQ']:.4f}  (target {HQ_target:.4f}, "
                  f"gap {abs(obs['HQ']-HQ_target)/HQ_target*100:.2f}%)")
            print(f"   2L/Q = {2*obs['LQ']:.4f}  (g_e = 2.00232, "
                  f"gap {abs(2*obs['LQ']-2.00232)/2.00232*100:.3f}%)")
            print(f"   r_rms = {obs['r_rms_fm']:.1f} fm")
    print()

    # 2. Neutral chaoiton at λ=1.0
    print("2. Neutral chaoiton at λ=1.0, g=0.5:")
    sol_n = solve_neutral(0.5, 1.0, 0.01)
    if sol_n:
        r_n, b_n = sol_n
        obs_n = compute_neutral_observables(r_n, b_n, 0.5, 1.0)
        if obs_n and obs_n['tail'] < 0.15:
            print(f"   B0=0.01, g=0.5, λ=1.0")
            print(f"   mass = {obs_n['mass_pred_MeV']:.4f} MeV  "
                  f"(prediction λ×m_e = {1.0*m_e_MeV:.4f} MeV)")
            print(f"   r_rms = {obs_n['r_rms_fm']:.1f} fm")
            print(f"   tail = {obs_n['tail']:.5f}")
    print()
    print("For full scans call:")
    print("  electron_calibration_scan()")
    print("  neutral_chaoiton_scan()")
    print("  lepton_spectrum_scan()")
