"""M6.0 full 4-function ansatz — toroidal/poloidal mutual confinement.

Werbos's 2017 paper specifies the chaoiton ansatz as 4 radial profile
functions on the toroidal/poloidal coordinate system:

    A_0(x) = alpha(r)         (poloidal scalar)
    A_phi(x) = beta(r)        (poloidal vector — wraps around central circle)
    J_0(x) = gamma(r)         (toroidal scalar)
    J_phi(x) = delta(r)       (toroidal vector — current loop)

The vector phi-components carry the angular momentum and live with the
mutual confinement mechanism Werbos called "Ouroboros" (yin-yang). The
2-function scalar reduction (alpha, gamma only) gave H/Q scaling ~ omega^2.04,
8% off from Werbos's claimed omega^2.22. This script tests whether the
4-function ansatz closes the remaining gap.

Lagrangian: L = -F_munu F^munu - G_munu G^munu + J^mu A_mu - g(J^mu J_mu)^2

In Minkowski metric (-,+,+,+): J^mu J_mu = -gamma^2 + delta^2
So f(J·J) = g*(-gamma^2 + delta^2)^2 = g*(gamma^2 - delta^2)^2

Variational derivatives in 4-function ansatz:
  delta L / delta alpha = gamma   (J^0 source for A^0)
  delta L / delta beta  = delta   (J_phi source for A_phi)
  delta L / delta gamma = alpha   (A^0 source for J^0) - 4g(gamma^2-delta^2)*gamma
  delta L / delta delta = beta    (A_phi source for J_phi) + 4g(gamma^2-delta^2)*delta

Vector Laplacian for phi-component in spherical (axisymmetric) coordinates:
  nabla^2 psi_phi = psi_phi'' + (2/r) psi_phi' - psi_phi / (r sin theta)^2
  averaged radially: nabla^2 ~ psi'' + (2/r) psi' - 2 psi / r^2

Hamiltonian (Werbos 2017 eq 18 approximation):
  h(r) = alpha'^2 + beta'^2 + gamma'^2 + delta'^2
       + omega^2 (alpha^2 + beta^2 + gamma^2 + delta^2)
       + g (gamma^2 - delta^2)^2
       - (alpha gamma - beta delta)   (J·A = J^0 A_0 - J^phi A_phi)
"""

from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


G_COUPLING = 1.0625
R_MIN = 0.02
R_MAX = 10.0
N_GRID = 800
RTOL = 1.0e-10
ATOL = 1.0e-12

INIT_AMP = 0.1
LOCALIZATION_THRESHOLD = 0.15

LEPTON_TARGETS = [
    ("electron", 1.0, 0.511),
    ("muon", 11.0, 105.7),
    ("pion+", 13.0, 139.6),
    ("tau", 40.7, 1777.0),
]
E_NAT = 0.30282

TARGET_HQ_AT_OMEGA_1 = 1.6969
TARGET_SCALING = 2.22


# ---------------------------------------------------------------------------
# 4-function ODE
# ---------------------------------------------------------------------------

def radial_ode_4fn(r, y, omega, g=G_COUPLING):
    """y = [alpha, alpha', beta, beta', gamma, gamma', delta, delta']

    Scalar (alpha, gamma): standard radial Laplacian (2/r) term
    Vector phi (beta, delta): vector radial Laplacian (2/r) term plus -2/r^2 curvature

    Couplings:
      alpha equation: source from gamma (J^0 sources A^0)
      beta  equation: source from delta (J_phi sources A_phi)
      gamma equation: source from alpha + mass from f
      delta equation: source from beta  + mass from f (opposite sign)

    f mass term: from f(s) = g s^2 with s = gamma^2 - delta^2
      df/dgamma = 4g (gamma^2 - delta^2) gamma
      df/ddelta = -4g (gamma^2 - delta^2) delta

    Using quarter-coupling per A4 calibration.
    """
    alpha, da, beta, db, gamma, dc, delta, dd = y
    inv_r = 1.0 / r
    s_minkowski = gamma * gamma - delta * delta
    mass_g = 4.0 * g * s_minkowski * gamma
    mass_d = -4.0 * g * s_minkowski * delta

    # alpha: scalar Laplacian + omega² + source
    d2a = -2*inv_r*da - omega*omega*alpha - gamma/4

    # beta: vector Laplacian (with -2/r^2 curvature) + omega² + source
    d2b = -2*inv_r*db + 2*beta*inv_r*inv_r - omega*omega*beta - delta/4

    # gamma: scalar Laplacian + omega² - mass + source
    d2c = -2*inv_r*dc - omega*omega*gamma + mass_g - alpha/4

    # delta: vector Laplacian - mass term + source
    d2d = -2*inv_r*dd + 2*delta*inv_r*inv_r - omega*omega*delta + mass_d - beta/4

    return [da, d2a, db, d2b, dc, d2c, dd, d2d]


def integrate_4fn(omega, amp=INIT_AMP, g=G_COUPLING):
    y0 = [amp, 0.0, amp, 0.0, amp, 0.0, amp, 0.0]
    r_eval = np.linspace(R_MIN, R_MAX, N_GRID)
    sol = solve_ivp(
        fun=lambda r, y: radial_ode_4fn(r, y, omega, g),
        t_span=(R_MIN, R_MAX), y0=y0, t_eval=r_eval,
        method="RK45", rtol=RTOL, atol=ATOL,
    )
    return sol


def compute_obs_4fn(sol, omega, g=G_COUPLING):
    r = sol.t
    a, da, b, db, c, dc, d, dd = sol.y
    s_mink = c*c - d*d

    # Hamiltonian density (Werbos 2017 eq 18)
    h_dens = (
        da*da + db*db + dc*dc + dd*dd
        + omega*omega*(a*a + b*b + c*c + d*d)
        + g * s_mink * s_mink
        - (a*c - b*d)
    )

    # Q: both fields contribute, no omega factor (per dual-criteria search winner)
    q_dens = a*a + b*b + c*c + d*d

    r2 = r * r
    H = 4.0 * np.pi * np.trapezoid(h_dens * r2, r)
    Q = 4.0 * np.pi * np.trapezoid(q_dens * r2, r)
    hq = H / Q if abs(Q) > 1e-12 else float("nan")
    tail = abs(a[-1]) + abs(b[-1]) + abs(c[-1]) + abs(d[-1])
    return {
        "H": float(H), "Q": float(Q), "HQ": float(hq), "tail": float(tail),
        "localized": bool(tail < LOCALIZATION_THRESHOLD),
        "field_max": float(np.max(np.abs(np.stack([a, b, c, d])))),
    }


# ---------------------------------------------------------------------------
# Run scan
# ---------------------------------------------------------------------------

def main():
    here = Path(__file__).parent
    plots_dir = here.parent / "plots"
    plots_dir.mkdir(exist_ok=True)

    print("=" * 78)
    print("FULL 4-FUNCTION ANSATZ — toroidal/poloidal mutual confinement")
    print("=" * 78)
    print(f"  g = {G_COUPLING}, initial amplitudes = {INIT_AMP} (all 4 fields)\n")

    print(f"  {'omega':>8}  {'H':>10}  {'Q':>10}  {'H/Q':>10}  "
          f"{'tail':>8}  {'loc':>5}  {'field_max':>10}")
    print("  " + "-" * 78)

    scan_omegas = [1.0, 2.0, 5.0, 8.0, 10.0, 11.0, 13.0, 15.0, 20.0, 30.0, 40.7, 50.0]
    results = []
    for omega in scan_omegas:
        sol = integrate_4fn(omega)
        if not sol.success:
            print(f"  {omega:8.2f}  integration failed")
            continue
        obs = compute_obs_4fn(sol, omega)
        results.append({"omega": omega, **obs, "sol": sol})
        loc_str = "Y" if obs["localized"] else "n"
        print(f"  {omega:8.2f}  {obs['H']:10.4f}  {obs['Q']:10.4f}  {obs['HQ']:10.4f}  "
              f"{obs['tail']:8.4f}  {loc_str:>5}  {obs['field_max']:10.4f}")

    # Scaling exponent
    omegas = np.array([r["omega"] for r in results])
    hqs = np.array([r["HQ"] for r in results])
    valid = (hqs > 0) & np.isfinite(hqs)
    slope, intercept = np.polyfit(np.log(omegas[valid]), np.log(hqs[valid]), 1)

    print(f"\n  Calibration: H/Q at ω=1 = {hqs[0]:.4f} (target {TARGET_HQ_AT_OMEGA_1})")
    print(f"  Scaling: H/Q ~ ω^{slope:.3f}  (target ω^{TARGET_SCALING}, "
          f"gap {abs(slope-TARGET_SCALING)/TARGET_SCALING*100:.2f}%)")

    # Lepton predictions
    print(f"\n  Lepton mass predictions (m = H/Q * e_nat):")
    print(f"  {'lepton':>10}  {'ω':>6}  {'H/Q':>10}  {'m_pred (MeV)':>14}  "
          f"{'m_obs (MeV)':>13}  {'gap%':>7}")
    print("  " + "-" * 75)
    for name, om_t, mass_obs in LEPTON_TARGETS:
        idx = next((i for i, o in enumerate(omegas) if abs(o - om_t) < 0.05), None)
        if idx is None:
            continue
        m_pred = hqs[idx] * E_NAT
        gap = abs(m_pred - mass_obs) / mass_obs * 100
        print(f"  {name:>10}  {om_t:6.2f}  {hqs[idx]:10.4f}  {m_pred:14.4f}  "
              f"{mass_obs:13.4f}  {gap:7.2f}")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Profiles at omega=1
    sol_e = next(r["sol"] for r in results if abs(r["omega"] - 1.0) < 0.01)
    r = sol_e.t
    a, b, c, d = sol_e.y[0], sol_e.y[2], sol_e.y[4], sol_e.y[6]
    ax1.plot(r, a, label="α (A poloidal scalar)", color="C0")
    ax1.plot(r, b, label="β (A poloidal vector)", color="C2")
    ax1.plot(r, c, label="γ (J toroidal scalar)", color="C1")
    ax1.plot(r, d, label="δ (J toroidal vector)", color="C3")
    ax1.set_xlabel("r (code units)")
    ax1.set_ylabel("amplitude")
    ax1.set_title(f"4-function chaoiton, ω=1.0\nH/Q = {hqs[0]:.4f} "
                  f"(target {TARGET_HQ_AT_OMEGA_1:.4f})")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Mass-frequency
    localized_mask = np.array([r["localized"] for r in results])
    ax2.scatter(omegas[~localized_mask], hqs[~localized_mask], color="lightgray",
                alpha=0.5, s=40, label="not localized")
    ax2.scatter(omegas[localized_mask], hqs[localized_mask], color="C0", s=80,
                label="localized")
    for name, om_t, _ in LEPTON_TARGETS:
        ax2.axvline(om_t, color="red", linestyle=":", alpha=0.5)
    x_fit = np.linspace(min(omegas), max(omegas), 100)
    y_fit = np.exp(intercept) * x_fit ** slope
    ax2.plot(x_fit, y_fit, color="C2", linestyle="--",
             label=f"fit: H/Q ~ ω^{slope:.3f}  (Werbos: ω^{TARGET_SCALING})")
    ax2.set_xlabel("omega")
    ax2.set_ylabel("H/Q")
    ax2.set_title(f"4-function H/Q scan, g={G_COUPLING}")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.legend(loc="lower right")
    ax2.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    out_path = plots_dir / "m6_0_full_ansatz.png"
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"\n  Plot saved: {out_path}")

    # Verdict
    print(f"\n{'=' * 78}")
    print("VERDICT")
    print(f"{'=' * 78}")
    cal_gap = abs(hqs[0] - TARGET_HQ_AT_OMEGA_1) / TARGET_HQ_AT_OMEGA_1
    scale_gap = abs(slope - TARGET_SCALING) / TARGET_SCALING
    print(f"\n  Calibration gap:  {cal_gap*100:.2f}%")
    print(f"  Scaling gap:      {scale_gap*100:.2f}%")
    if cal_gap < 0.05 and scale_gap < 0.05:
        print("  → 4-function ansatz REPRODUCES Werbos. M6 GO.")
    elif cal_gap < 0.05 and scale_gap < 0.15:
        print("  → Improvement from scalar reduction; partial reproduction.")
    elif cal_gap > 0.20:
        print("  → 4-function ansatz BREAKS calibration. Wrong direction.")
    else:
        print("  → 4-function ansatz doesn't close the gap. Framework limitations confirmed.")


if __name__ == "__main__":
    main()
