#!/usr/bin/env python3
"""
m6_v11_neutral_gl_scan.py

Answers Paul/DeepSeek's 2026-06-05 request ("an urgent problem after all"):
run a parameter scan over (g, lambda) in the neutral sector to find the
values that yield m_chi = 0.460 MeV and m_J = 0.618 MeV, deliver the
definitive soliton profile beta(r) for the J-A mixing integral, and
resolve the "BVP inconsistency" flagged in the DM paper v5a Technical Note.

KEY ANALYTIC RESULT (verified numerically below): the neutral-sector ODE

    beta'' + (2/r) beta' - (2/r^2) beta + lambda*beta + 4g*beta^3 = 0,
    lambda = -m_J^2 < 0                                  (Q43+Q44 form)

has an EXACT scaling symmetry. With  x = m_J * r,  a = m_J / (2 sqrt(g)):

    beta(r) = a * u(x),   u'' + (2/x)u' - (2/x^2)u - u + u^3 = 0

where u(x) is a single UNIVERSAL profile with no parameters. Consequences:

  1. The entire (g, lambda<0) quarter-plane is one solution family, all
     exact rescalings of u(x).
  2. eta = m_J * eta_u  and  H/Q = m_J^2 * K_u, where
     eta_u = (int u^2 x dx)/(int u^2 x^2 dx) and K_u are universal numbers.
  3. m_J_phys = (m_J/eta) * m_e = m_e / eta_u  -- INDEPENDENT of (g, lambda).
     This is why v10 found m_J/eta = 1.21024 family-invariant: it is the
     scaling identity, and m_J = 0.6184 MeV is a parameter-free prediction.
  4. m_chi = eta*(H/Q)*m_e = (-lambda)^(3/2) * eta_u * K_u * m_e
     -- depends ONLY on lambda, not on g. m_chi = 0.460 MeV selects
     lambda* = -0.2647 (i.e. m_J_raw = 0.5145) for ANY g.
  5. g is a pure amplitude scale for the masses (beta ~ 1/sqrt(g)); it does
     NOT move m_chi or m_J. It DOES set the absolute amplitude of beta(r),
     so it matters for the J-A mixing integral; electron universality
     (Q45: same Lagrangian parameters) fixes g = 1.0.
  6. At fixed B0 (v10's chart), m_J = (2 sqrt(g) B0 / u'(0))^(1/2), i.e.
     m_J ~ g^(1/4) -- reproduces the v10 g-scan table exactly.

Also reproduces the two failure modes in the v5a Technical Note:
  - "Griesi-1 finds near-trivial solutions (peak ~1e-11) at all B0
    initializations": at fixed (g, lambda) the trivial solution beta=0
    satisfies the BCs, and Newton collapses to it from a cold-start
    amplitude guess. The scaling law gives the correct seed amplitude
    a = m_J/(2 sqrt(g)); with it, the same code converges to the soliton.
    (Our actual v9/v10 production code never had this failure: it pins the
    amplitude via B0 as a hard BC and frees m_J as the eigenvalue.)
  - Wrong l=1 regularity class (template BC beta'(RMIN)=0): collapses to
    trivial from cold seeds (the v9 phase 2 observation). With a
    scaling-quality seed the BC error at RMIN is only O(RMIN^2) and the
    solver still lands on the soliton — empirically the trivial Newton
    basin (cold-start amplitude) is the DOMINANT failure mechanism.

Outputs: scan tables + canonical profile CSV (natural + physical units)
for the alpha_JN / J-A mixing-integral calculation.
"""

import csv
import os

import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import interp1d
from scipy.optimize import brentq, curve_fit

# --- physics constants ---
M_E_MEV = 0.5109989461           # electron mass in MeV (PDG)
HBAR_C_MEV_FM = 197.3269804      # hbar*c in MeV*fm
R_PHYS_FM = HBAR_C_MEV_FM / M_E_MEV   # ~386.16 fm (electron Compton length)

RMIN = 0.02
N_GRID = 800

OUTDIR = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------
# v9/v10 production formulation: B0 fixed (hard BC), m_J free eigenvalue
# ----------------------------------------------------------------------

def run_bvp_mJ_free(B0, g=1.0, m_J_init=1.0, rmax=20.0):
    """The v9 phase 2 / v10 canonical solver (regression anchor)."""

    def ode(r, y, p):
        beta, betap = y
        m_J = p[0]
        betapp = (-(2.0/r)*betap + (2.0/(r*r))*beta
                  + (m_J*m_J)*beta - 4.0*g*beta**3)
        return np.vstack((betap, betapp))

    def bc(ya, yb, p):
        m_J = p[0]
        return np.array([
            ya[0] - B0 * RMIN,        # beta(RMIN) = B0*RMIN   (l=1: beta ~ B0*r)
            ya[1] - B0,               # beta'(RMIN) = B0
            yb[1] + m_J * yb[0],      # Robin K_1 decay at RMAX
        ])

    r_mesh = np.linspace(RMIN, rmax, N_GRID)
    beta_guess = B0 * r_mesh * np.exp(-m_J_init * r_mesh)
    betap_guess = B0 * (1.0 - m_J_init * r_mesh) * np.exp(-m_J_init * r_mesh)
    return solve_bvp(ode, bc, r_mesh, np.vstack((beta_guess, betap_guess)),
                     p=[m_J_init], tol=1e-8, max_nodes=40000)


# ----------------------------------------------------------------------
# Direct (g, lambda) formulation: the scan DeepSeek asked for.
# lambda fixed => m_J = sqrt(-lambda) fixed, NO free parameter.
# ----------------------------------------------------------------------

def run_bvp_fixed_lambda(g, lam, seed='scaling', B0_init=None, bc_variant='l1',
                         u_interp=None):
    """Solve the neutral BVP at FIXED (g, lambda).

    seed='scaling'  -> amplitude-correct guess a*u(m_J r) built from the
                       universal profile interpolant u_interp (converges to
                       the soliton in a few Newton steps)
    seed='cold'     -> naive B0_init * r * exp(-m_J r) guess (Technical-Note
                       failure mode: Newton collapses to trivial beta=0)
    bc_variant='l1' -> correct l=1 regularity (beta - r*beta' = 0 at RMIN)
    bc_variant='template' -> DeepSeek's verbatim template beta'(RMIN)=0
                       (wrong regularity class; collapses to trivial)
    """
    assert lam < 0.0, "neutral sector requires lambda < 0 (K_1 decay, Q43)"
    m_J = np.sqrt(-lam)
    rmax = max(20.0, 14.0 / m_J)

    def ode(r, y):
        beta, betap = y
        betapp = (-(2.0/r)*betap + (2.0/(r*r))*beta
                  - lam*beta - 4.0*g*beta**3)
        return np.vstack((betap, betapp))

    def bc_l1(ya, yb):
        return np.array([
            ya[0] - RMIN * ya[1],     # l=1 regularity: beta ~ B0*r (B0 free)
            yb[1] + m_J * yb[0],      # Robin K_1 decay at RMAX
        ])

    def bc_template(ya, yb):
        return np.array([
            ya[1],                    # beta'(RMIN) = 0  (WRONG class: beta~r^2)
            yb[1] + m_J * yb[0],
        ])

    bc = bc_l1 if bc_variant == 'l1' else bc_template
    r_mesh = np.linspace(RMIN, rmax, N_GRID)

    if seed == 'scaling':
        # scaling law: beta(r) = a * u(m_J r), amplitude a = m_J/(2 sqrt(g))
        assert u_interp is not None, "scaling seed needs the universal profile"
        a = m_J / (2.0*np.sqrt(g))
        beta_guess = a * u_interp(m_J * r_mesh)
        betap_guess = np.gradient(beta_guess, r_mesh)
    else:
        B0 = 0.5 if B0_init is None else B0_init
        beta_guess = B0 * r_mesh * np.exp(-m_J * r_mesh)
        betap_guess = B0 * (1.0 - m_J * r_mesh) * np.exp(-m_J * r_mesh)

    return solve_bvp(ode, bc, r_mesh, np.vstack((beta_guess, betap_guess)),
                     tol=1e-8, max_nodes=40000)


# ----------------------------------------------------------------------
# Observables (same definitions as v10)
# ----------------------------------------------------------------------

def observables(r, b, m_J, g):
    dr = np.diff(r)
    rm = 0.5*(r[:-1] + r[1:])
    bm = 0.5*(b[:-1] + b[1:])
    dbm = np.diff(b)/dr

    Q = np.sum(bm**2 * rm**2 * dr)
    H = np.sum((dbm**2 + 2*bm**2/rm**2 + m_J**2*bm**2 + 4*g*bm**4) * rm**2 * dr)
    int_b2_r = np.sum(bm**2 * rm * dr)
    eta = int_b2_r / Q if Q > 1e-15 else float('nan')
    HQ = H / Q if Q > 1e-15 else float('nan')

    peak_idx = np.argmax(np.abs(b))
    peak = abs(b[peak_idx])
    out = {
        'H': H, 'Q': Q, 'HQ': HQ, 'eta': eta,
        'peak': peak, 'peak_r': r[peak_idx],
        'tail_over_peak': abs(b[-1])/peak if peak > 1e-12 else float('inf'),
        'sign_changes': int(np.sum(np.diff(np.sign(b[np.abs(b) > 1e-12])) != 0)),
        'B0_out': b[0]/r[0],                       # recovered slope at origin
        'mJ_over_eta': m_J/eta,                    # should be 1/eta_u = 1.21024
        'HQ_over_mJ2': HQ/m_J**2,                  # should be K_u (universal)
        'm_chi_MeV': eta*HQ*M_E_MEV,               # Q46 recipe
        'm_J_MeV': (m_J/eta)*M_E_MEV,              # Q46 recipe
    }
    return out


def extract_C(r, b, m_J, r_min_fit_factor=0.4):
    """Far-field K_1 amplitude fit (same as v10): beta ~ C e^{-mr}(1+1/(mr))/r."""
    r_min_fit = r_min_fit_factor * r[-1]
    mask = (r >= r_min_fit) & (np.abs(b) > 1e-12)
    if mask.sum() < 5:
        return float('nan'), float('nan')
    rr, bb = r[mask], b[mask]

    def k1_asym(r, C):
        return C * np.exp(-m_J*r) * (1.0 + 1.0/(m_J*r)) / r

    popt, _ = curve_fit(k1_asym, rr, bb, p0=[1.0], maxfev=4000)
    rel = np.max(np.abs(bb - k1_asym(rr, popt[0])) / np.maximum(np.abs(bb), 1e-15))
    return float(popt[0]), float(rel)


# ----------------------------------------------------------------------
# Main work program
# ----------------------------------------------------------------------

def main():
    np.set_printoptions(precision=6)
    print("=" * 78)
    print("M6 v11 — neutral sector (g, lambda) scan  [Paul/DeepSeek 2026-06-05 request]")
    print("=" * 78)

    # ---- STEP 1: regression anchor — reproduce v9/v10 canonical point ----
    print("\nSTEP 1 — Regression anchor: v9/v10 solver at (g=1.0, B0=0.5)")
    print("-" * 78)
    sol = run_bvp_mJ_free(B0=0.5, g=1.0)
    assert sol.success, sol.message
    m_J_c = float(sol.p[0])
    obs_c = observables(sol.x, sol.y[0], m_J_c, 1.0)
    print(f"  m_J_raw = {m_J_c:+.6f}   (v10: +0.514460)")
    print(f"  eta = {obs_c['eta']:.6f}   (v10: 0.4251)")
    print(f"  (H/Q)_spherical = {obs_c['HQ']:.6f}   (v10: 2.1173)")
    print(f"  m_chi = {obs_c['m_chi_MeV']:.4f} MeV   m_J = {obs_c['m_J_MeV']:.4f} MeV")
    print(f"  peak beta = {obs_c['peak']:.4f} @ r = {obs_c['peak_r']:.3f}   "
          f"nodes = {obs_c['sign_changes']}   tail/peak = {obs_c['tail_over_peak']:.2e}")
    print(f"  >> NOTE vs v5a Technical Note: this IS direct BVP integration of the")
    print(f"     ODE; peak amplitude is O(1), not 1e-11. Near-trivial collapse only")
    print(f"     occurs in the fixed-lambda formulation with a cold seed (Step 5).")

    # universal profile u(x) = beta(x/m_J)/a from the anchor (g=1 -> a = m_J/2)
    a_c = m_J_c / 2.0
    r_dense = np.linspace(RMIN, sol.x[-1], 4000)
    u_interp = interp1d(m_J_c * r_dense, sol.sol(r_dense)[0] / a_c,
                        bounds_error=False, fill_value=0.0)

    # universal constants from the anchor
    eta_u = obs_c['eta'] / m_J_c
    K_u = obs_c['HQ'] / m_J_c**2
    print(f"\n  Universal constants of u(x):  eta_u = {eta_u:.6f}   K_u = {K_u:.6f}")
    print(f"  Scaling-law predictions:")
    print(f"    m_J_phys = m_e/eta_u = {M_E_MEV/eta_u:.6f} MeV   for ALL (g, lambda)")
    print(f"    m_chi(lambda) = (-lambda)^(3/2) * eta_u*K_u*m_e "
          f"= {eta_u*K_u*M_E_MEV:.6f} * (-lambda)^(3/2) MeV   (g-independent)")

    # ---- STEP 2: the requested (g, lambda) grid scan ----
    print("\nSTEP 2 — Direct (g, lambda) grid scan (fixed-lambda BVP, scaling seed)")
    print("-" * 78)
    lam_star_guess = -m_J_c**2
    g_values = [0.5, 1.0, 2.0]
    lam_values = [-0.09, -0.16, lam_star_guess, -0.36, -0.64, -1.0]
    print(f"  {'g':>5} {'lambda':>10} {'m_J_raw':>9} {'B0_out':>9} {'peak':>8} "
          f"{'nodes':>5} {'m_chi(MeV)':>11} {'m_J(MeV)':>9} {'H/Q/mJ^2':>9} {'mJ/eta':>8}")
    grid = []
    for g in g_values:
        for lam in lam_values:
            s = run_bvp_fixed_lambda(g, lam, seed='scaling', u_interp=u_interp)
            if not s.success:
                print(f"  {g:>5.2f} {lam:>10.6f}   FAILED: {s.message[:45]}")
                continue
            m_J = np.sqrt(-lam)
            o = observables(s.x, s.y[0], m_J, g)
            grid.append((g, lam, o))
            print(f"  {g:>5.2f} {lam:>10.6f} {m_J:>9.5f} {o['B0_out']:>9.5f} "
                  f"{o['peak']:>8.4f} {o['sign_changes']:>5d} {o['m_chi_MeV']:>11.4f} "
                  f"{o['m_J_MeV']:>9.4f} {o['HQ_over_mJ2']:>9.5f} {o['mJ_over_eta']:>8.5f}")
    print(f"\n  >> m_J(MeV) column is constant and m_chi depends only on lambda:")
    print(f"     the masses CANNOT independently fix (g, lambda). lambda fixes m_chi;")
    print(f"     g only rescales the profile amplitude (beta ~ 1/sqrt(g)).")

    # ---- STEP 3: locus where m_chi = 0.460 MeV ----
    print("\nSTEP 3 — Locus m_chi = 0.460 MeV (brentq on lambda, per g)")
    print("-" * 78)
    M_CHI_TARGET = 0.460

    def m_chi_of_lambda(lam, g):
        s = run_bvp_fixed_lambda(g, lam, seed='scaling', u_interp=u_interp)
        if not s.success:
            return float('nan')
        return observables(s.x, s.y[0], np.sqrt(-lam), g)['m_chi_MeV'] - M_CHI_TARGET

    lam_analytic = -(M_CHI_TARGET / (eta_u*K_u*M_E_MEV))**(2.0/3.0)
    print(f"  Analytic (scaling law): lambda* = {lam_analytic:.6f} "
          f"(m_J_raw = {np.sqrt(-lam_analytic):.6f})")
    for g in g_values:
        try:
            lam_s = brentq(m_chi_of_lambda, -0.36, -0.16, args=(g,), xtol=1e-7)
            print(f"  g = {g:>4.2f}:  lambda* = {lam_s:.6f}   "
                  f"(m_J_raw = {np.sqrt(-lam_s):.6f})")
        except ValueError as e:
            print(f"  g = {g:>4.2f}:  bracket failed ({e})")

    # ---- STEP 4: canonical solution at (g=1, lambda*) + profile export ----
    print("\nSTEP 4 — Definitive profile at canonical (g=1, lambda*) for alpha_JN")
    print("-" * 78)
    g_canon = 1.0
    lam_canon = brentq(m_chi_of_lambda, -0.36, -0.16, args=(g_canon,), xtol=1e-9)
    s = run_bvp_fixed_lambda(g_canon, lam_canon, seed='scaling', u_interp=u_interp)
    assert s.success
    m_J = np.sqrt(-lam_canon)
    r_hi = np.linspace(RMIN, s.x[-1], 4000)
    b_hi = s.sol(r_hi)[0]
    bp_hi = s.sol(r_hi)[1]
    o = observables(r_hi, b_hi, m_J, g_canon)
    C_nat, C_resid = extract_C(r_hi, b_hi, m_J)
    C_MeV_fm = C_nat * R_PHYS_FM * M_E_MEV
    print(f"  (g, lambda) = ({g_canon}, {lam_canon:.6f})   m_J_raw = {m_J:.6f}")
    print(f"  B0_out = {o['B0_out']:.6f}  (v9/v10 canonical B0 = 0.5 — same point)")
    print(f"  m_chi = {o['m_chi_MeV']:.4f} MeV   m_J = {o['m_J_MeV']:.4f} MeV")
    print(f"  C = {C_MeV_fm:.1f} MeV*fm  (K_1 fit residual {C_resid:.2e})")
    print(f"  peak beta = {o['peak']:.4f} @ r = {o['peak_r']:.3f} "
          f"({o['peak_r']*R_PHYS_FM:.1f} fm)   nodes = {o['sign_changes']}   "
          f"tail/peak = {o['tail_over_peak']:.2e}")

    csv_path = os.path.join(OUTDIR, 'v11_canonical_beta_profile.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow([f'# M6 v11 canonical neutral chaoiton profile beta(r)'])
        w.writerow([f'# ODE: beta\'\' + (2/r)beta\' - (2/r^2)beta + lambda*beta'
                    f' + 4g*beta^3 = 0'])
        w.writerow([f'# g = {g_canon}, lambda = {lam_canon:.8f},'
                    f' m_J_raw = sqrt(-lambda) = {m_J:.8f}'])
        w.writerow([f'# eta = {o["eta"]:.8f}, H/Q_spherical = {o["HQ"]:.8f}'])
        w.writerow([f'# m_chi = {o["m_chi_MeV"]:.6f} MeV, m_J = {o["m_J_MeV"]:.6f} MeV,'
                    f' C = {C_MeV_fm:.2f} MeV*fm'])
        w.writerow([f'# natural length unit R_phys = hbar*c/m_e = {R_PHYS_FM:.4f} fm'])
        w.writerow(['r_natural', 'r_fm', 'beta_natural', 'dbeta_dr_natural'])
        for ri, bi, bpi in zip(r_hi, b_hi, bp_hi):
            w.writerow([f'{ri:.8f}', f'{ri*R_PHYS_FM:.4f}',
                        f'{bi:.10e}', f'{bpi:.10e}'])
    print(f"  Profile exported: {csv_path}  ({len(r_hi)} points)")

    # ---- STEP 5: reproduce the Technical-Note failure modes ----
    print("\nSTEP 5 — Reproducing the v5a Technical-Note failure modes")
    print("-" * 78)
    print("  5a. Fixed-lambda BVP, COLD seed (the 'Griesi-1 near-trivial' mode):")
    for B0i in [0.1, 0.5, 1.0]:
        s5 = run_bvp_fixed_lambda(1.0, -1.0, seed='cold', B0_init=B0i)
        pk = np.max(np.abs(s5.y[0])) if s5.success else float('nan')
        print(f"      B0_init = {B0i:>4.1f}:  success={s5.success}  "
              f"peak = {pk:.3e}  ({'TRIVIAL collapse' if pk < 1e-6 else 'soliton'})")
    s5 = run_bvp_fixed_lambda(1.0, -1.0, seed='scaling', u_interp=u_interp)
    pk = np.max(np.abs(s5.y[0]))
    o5 = observables(s5.x, s5.y[0], 1.0, 1.0)
    print(f"      scaling seed:   success={s5.success}  peak = {pk:.4f}  "
          f"nodes = {o5['sign_changes']}  (SOLITON — same code, correct seed)")
    print(f"      >> required amplitude at lambda=-1, g=1 is a = m_J/(2 sqrt(g)) = 0.5;")
    print(f"         cold r*exp(-r) seeds sit in the trivial Newton basin.")
    print(f"  5b. Wrong l=1 regularity (template BC beta'(RMIN)=0):")
    for sd, kw in [('cold', dict(seed='cold', B0_init=0.5)),
                   ('scaling', dict(seed='scaling', u_interp=u_interp))]:
        s5b = run_bvp_fixed_lambda(1.0, -1.0, bc_variant='template', **kw)
        pk_b = np.max(np.abs(s5b.y[0])) if s5b.success else float('nan')
        print(f"      {sd:>7} seed:  success={s5b.success}  peak = {pk_b:.3e}  "
              f"({'TRIVIAL collapse' if pk_b < 1e-6 else 'soliton'})")
    print(f"      >> wrong BC collapses from cold seeds (v9 phase 2 observation);")
    print(f"         with a scaling-quality seed the RMIN BC error is O(RMIN^2) and")
    print(f"         the soliton still wins — the trivial Newton basin (cold-start")
    print(f"         amplitude) is the dominant failure mechanism.")
    print(f"  5c. The v9/v10 production formulation avoids both: B0 is a hard BC")
    print(f"      (amplitude pinned, trivial solution excluded) and m_J is the free")
    print(f"      eigenvalue. Shooting/IVP failure is separate — documented Q42.")

    # ---- STEP 6: scaling-law cross-checks ----
    print("\nSTEP 6 — Scaling-law cross-checks")
    print("-" * 78)
    mJoEta = [o['mJ_over_eta'] for _, _, o in grid]
    HQm2 = [o['HQ_over_mJ2'] for _, _, o in grid]
    print(f"  m_J/eta over grid:  mean = {np.mean(mJoEta):.6f}  "
          f"spread = {np.max(mJoEta)-np.min(mJoEta):.2e}   (v10 invariant: 1.21024)")
    print(f"  (H/Q)/m_J^2 over grid:  mean = {np.mean(HQm2):.6f}  "
          f"spread = {np.max(HQm2)-np.min(HQm2):.2e}")
    print(f"  >> (H/Q)/m_J^2 = {np.mean(HQm2):.4f} ~ 8: candidate exact virial")
    print(f"     (Pohozaev) identity H = 8 m_J^2 Q for the l=1 cubic NLS — would")
    print(f"     explain v10's Q47 invariance analytically. For Paul/DeepSeek.")
    # g^(1/4) law at fixed B0 (v10 chart)
    print(f"\n  g^(1/4) law at fixed B0=0.5 (v10 g-scan reproduction):")
    for g in [0.5, 1.0, 1.6]:
        sg = run_bvp_mJ_free(B0=0.5, g=g)
        pred = m_J_c * (g/1.0)**0.25
        print(f"    g = {g:>4.2f}:  m_J = {float(sg.p[0]):+.5f}   "
              f"scaling prediction m_J(1)*g^0.25 = {pred:+.5f}")

    print("\n" + "=" * 78)
    print("SUMMARY — answers to DeepSeek's three asks")
    print("=" * 78)
    print(f"""
  1. Neutral-sector parameters from first principles:
       The (g, lambda) scan is analytically degenerate (exact scaling family).
       m_J = {M_E_MEV/eta_u:.4f} MeV is parameter-free (invariant over the whole plane).
       m_chi = 0.460 MeV fixes lambda* = {lam_canon:.6f} (any g).
       g fixed at 1.0 by electron universality (Q45) -> amplitude determined.
       Canonical: (g, lambda) = (1.0, {lam_canon:.4f})  [= v9/v10 point (g=1, B0=0.5)]
  2. Definitive beta(r) profile: v11_canonical_beta_profile.csv (natural + fm).
  3. BVP inconsistency resolved:
       - 'near-trivial at all B0': cold-seed Newton basin at fixed lambda
         (5a) or wrong regularity BC (5b). Not present in v9/v10 formulation.
       - 'shooting code wrong masses': IVP windowing artifact, documented Q42.
       - v5a Technical Note correction needed: our 0.460/0.618 WERE direct BVP
         integration (+ Q46 eta recipe), not an energy-functional fit.
""")


if __name__ == '__main__':
    main()
