"""M6.0 — Werbos's Ouroboros mass-frequency reproduction.

Independently reproduces the published lepton mass-spectrum claim from
Werbos's Particle Spectrum paper (Zenodo 20242421, May 2026): at fixed
(g = 1.0625, lambda = 1.0), scanning the chaoiton oscillation frequency
omega produces electron (omega=1.0), muon (omega~11.0), pion (omega~13.0),
tau (omega~40.7) with mass scaling m ~ omega^2.22.

Gate 0 is the electron calibration check: at omega=1.0, the published
ratio H/Q = 1.6969 must reproduce within 5% before the scan is trusted.

See ../0a_background.md section 12 for the decision framework around
this reproduction's pass/fail.
"""

from pathlib import Path
import json

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Werbos's published parameters (from Calibration + Spectrum papers)
# ---------------------------------------------------------------------------

G_COUPLING = 1.0625      # f(s) = G_COUPLING * s^2, with s = J . J
LAMBDA = 1.0             # Lagrange multiplier (chemical potential in H' = H - lambda*Q)
A_INIT = 0.1             # initial A-field amplitude at r = R_MIN
B_INIT = 0.1             # initial J-field amplitude at r = R_MIN

R_MIN = 0.02             # inner radial cutoff (regularity)
R_MAX = 10.0             # outer integration limit (Calibration paper used 10.0; Spectrum 15.0)
N_GRID = 800             # grid points for output

RTOL = 1.0e-10
ATOL = 1.0e-12

LOCALIZATION_THRESHOLD = 0.15   # |A(r_max)| + |J(r_max)| < threshold for "stable chaoiton"

# Calibration target (Calibration paper section 4.2)
CALIBRATION_TARGET_HQ = 1.6969
CALIBRATION_TOLERANCE = 0.05    # 5% gap allowed at gate 0

# Lepton spectrum targets (Spectrum paper)
LEPTON_TARGETS = [
    # (name, omega, mass_MeV_observed)
    ("electron", 1.0, 0.511),
    ("muon", 11.0, 105.7),
    ("pion+", 13.0, 139.6),
    ("tau", 40.7, 1777.0),
]


# ---------------------------------------------------------------------------
# Radial ODE — simplest two-scalar reduction of the toroidal-poloidal ansatz
# ---------------------------------------------------------------------------
# Lagrangian: L = -F.F - G.G + J.A - g(J.J)^2
#
# With oscillatory ansatz psi(t,r) = psi(r) * exp(-i omega t), the wave
# equation becomes a Helmholtz-like equation for the spatial profile.
# In 3D spherical radial form, nabla^2 -> psi'' + (2/r) psi'.
#
# The simplest scalar reduction (treats poloidal A_0 and toroidal J_0 as
# scalar profiles alpha(r), gamma(r); ignores the d_phi vector components
# in this first-pass scaffold):
#
#   alpha'' + (2/r) alpha' + omega^2 alpha + gamma = 0       (Maxwell-like)
#   gamma'' + (2/r) gamma' + (omega^2 - 2g gamma^2) gamma + alpha = 0   (KG-like)
#
# State vector: y = [alpha, alpha', gamma, gamma']

def radial_ode(r, y, omega, g):
    alpha, dalpha, gamma, dgamma = y

    # Avoid r=0 singularity via R_MIN > 0 in the IVP setup
    inv_r = 1.0 / r

    # Mass-like term from f(s) = g s^2 with s = gamma^2 (scalar reduction):
    # df/ds = 2g s, so f'(gamma^2) * d(gamma^2)/d(gamma) = 2g gamma^2 * (effective coupling)
    mass_term = 2.0 * g * gamma * gamma

    d2alpha = -2.0 * inv_r * dalpha - omega * omega * alpha - gamma
    d2gamma = -2.0 * inv_r * dgamma - (omega * omega - mass_term) * gamma - alpha

    return [dalpha, d2alpha, dgamma, d2gamma]


# ---------------------------------------------------------------------------
# Integration + observables
# ---------------------------------------------------------------------------

def integrate_chaoiton(omega, g=G_COUPLING, a0=A_INIT, b0=B_INIT,
                      r_min=R_MIN, r_max=R_MAX, n_grid=N_GRID):
    """Integrate the radial ODE outward from r_min to r_max.

    Returns (r_grid, alpha, dalpha, gamma, dgamma) on an evenly spaced grid.
    """
    y0 = [a0, 0.0, b0, 0.0]
    r_eval = np.linspace(r_min, r_max, n_grid)

    sol = solve_ivp(
        fun=lambda r, y: radial_ode(r, y, omega, g),
        t_span=(r_min, r_max),
        y0=y0,
        t_eval=r_eval,
        method="RK45",
        rtol=RTOL,
        atol=ATOL,
    )
    if not sol.success:
        return None

    return sol.t, sol.y[0], sol.y[1], sol.y[2], sol.y[3]


def is_localized(alpha, gamma, threshold=LOCALIZATION_THRESHOLD):
    """Werbos's localization criterion: |A(r_max)| + |J(r_max)| < threshold."""
    return abs(alpha[-1]) + abs(gamma[-1]) < threshold


def compute_observables(r, alpha, dalpha, gamma, dgamma, omega, g=G_COUPLING):
    """Compute the chaoiton observables: H (energy), Q (charge), L (angular momentum).

    All in code units. Volume element for 3D radial integral: 4 pi r^2 dr.

    Hamiltonian density (radial, in scalar reduction):
        h(r) = alpha'^2 + omega^2 alpha^2     (A-field kinetic + oscillation)
             + gamma'^2 + omega^2 gamma^2     (J-field kinetic + oscillation)
             + g gamma^4                     (f(J.J) potential)
             - alpha gamma                   (coupling, sign by convention)

    Charge density (radial, U(1) Noether on the J-field oscillation):
        q(r) = omega gamma^2

    Angular momentum density:
        l(r) = omega^2 gamma^2 = omega * q(r)
        (consistent with L/Q = omega identity from Calibration paper section 3.1)
    """
    r2 = r * r
    h_density = (
        dalpha * dalpha + omega * omega * alpha * alpha
        + dgamma * dgamma + omega * omega * gamma * gamma
        + g * gamma * gamma * gamma * gamma
        - alpha * gamma
    )
    q_density = omega * gamma * gamma
    l_density = omega * q_density

    # 4 pi r^2 weighted integral over [r_min, r_max]
    H = 4.0 * np.pi * np.trapezoid(h_density * r2, r)
    Q = 4.0 * np.pi * np.trapezoid(q_density * r2, r)
    L = 4.0 * np.pi * np.trapezoid(l_density * r2, r)

    return H, Q, L


def run_chaoiton(omega, g=G_COUPLING, verbose=False):
    """Integrate + compute observables + localization for one omega."""
    result = integrate_chaoiton(omega, g)
    if result is None:
        return None

    r, alpha, dalpha, gamma, dgamma = result
    H, Q, L = compute_observables(r, alpha, dalpha, gamma, dgamma, omega, g)
    localized = is_localized(alpha, gamma)
    tail = abs(alpha[-1]) + abs(gamma[-1])

    if verbose:
        print(f"  omega = {omega:7.3f}: H = {H:9.4f}, Q = {Q:8.4f}, "
              f"H/Q = {H/Q if abs(Q) > 1e-12 else float('nan'):8.4f}, "
              f"L/Q = {L/Q if abs(Q) > 1e-12 else float('nan'):8.4f}, "
              f"tail = {tail:7.4f}, localized = {localized}")

    return {
        "omega": omega,
        "r": r, "alpha": alpha, "dalpha": dalpha,
        "gamma": gamma, "dgamma": dgamma,
        "H": float(H), "Q": float(Q), "L": float(L),
        "HQ": float(H / Q) if abs(Q) > 1e-12 else float("nan"),
        "LQ": float(L / Q) if abs(Q) > 1e-12 else float("nan"),
        "tail": float(tail),
        "localized": bool(localized),
    }


# ---------------------------------------------------------------------------
# Gate 0 — electron calibration check
# ---------------------------------------------------------------------------

def gate_0_calibration():
    """Verify electron calibration: H/Q at omega=1.0 must match 1.6969 within 5%."""
    print("\n=== GATE 0 — electron calibration check ===")
    print(f"  Target: H/Q = {CALIBRATION_TARGET_HQ} at (omega=1.0, g={G_COUPLING})")
    print(f"  Tolerance: {CALIBRATION_TOLERANCE * 100:.1f}%\n")

    res = run_chaoiton(omega=1.0, verbose=True)
    if res is None:
        print("  RESULT: integration failed — investigate before scanning.")
        return False, None

    hq = res["HQ"]
    gap = abs(hq - CALIBRATION_TARGET_HQ) / CALIBRATION_TARGET_HQ
    passed = gap < CALIBRATION_TOLERANCE

    print(f"\n  H/Q (computed): {hq:.4f}")
    print(f"  Target:          {CALIBRATION_TARGET_HQ:.4f}")
    print(f"  Gap:             {gap * 100:.2f}%")
    print(f"  Verdict:         {'PASS' if passed else 'FAIL'}")

    if not passed:
        print("\n  Gate 0 failed. The radial ODE form needs refinement before the scan.")
        print("  Likely refinements:")
        print("    (1) Use the full four-function ansatz (alpha, beta, gamma, delta)")
        print("    (2) Adjust coupling sign / convention in radial_ode()")
        print("    (3) Adjust H, Q definitions in compute_observables()")
        print("    (4) Reach out to Werbos for his actual scipy code")

    return passed, res


# ---------------------------------------------------------------------------
# Omega scan — mass spectrum
# ---------------------------------------------------------------------------

def omega_scan(omegas, g=G_COUPLING):
    """Scan over omega values, return list of result dicts."""
    print(f"\n=== OMEGA SCAN ({len(omegas)} values, g={g}) ===\n")
    results = []
    for omega in omegas:
        res = run_chaoiton(omega, g, verbose=True)
        if res is not None:
            results.append(res)
    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_calibration(res, out_path):
    """Plot alpha(r) and gamma(r) profiles at omega=1.0."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(res["r"], res["alpha"], label="alpha (A-field)", color="C0")
    ax1.plot(res["r"], res["gamma"], label="gamma (J-field)", color="C1")
    ax1.set_xlabel("r (code units)")
    ax1.set_ylabel("amplitude")
    ax1.set_title(f"Chaoiton radial profiles at omega=1.0\n"
                  f"H/Q = {res['HQ']:.4f} (target {CALIBRATION_TARGET_HQ:.4f})")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Log scale on second axis to show tail decay
    ax2.semilogy(res["r"], np.abs(res["alpha"]) + 1e-16, label="|alpha|", color="C0")
    ax2.semilogy(res["r"], np.abs(res["gamma"]) + 1e-16, label="|gamma|", color="C1")
    ax2.axhline(LOCALIZATION_THRESHOLD, color="red", linestyle="--",
                label=f"localization threshold ({LOCALIZATION_THRESHOLD})")
    ax2.set_xlabel("r (code units)")
    ax2.set_ylabel("|field| (log)")
    ax2.set_title(f"Tail decay (localized: {res['localized']})")
    ax2.legend()
    ax2.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out_path}")


def plot_mass_frequency_scan(results, out_path):
    """Plot computed mass-equivalent vs omega, with Werbos's lepton targets overlaid."""
    omegas_localized = [r["omega"] for r in results if r["localized"]]
    h_localized = [r["H"] for r in results if r["localized"]]
    omegas_all = [r["omega"] for r in results]
    h_all = [r["H"] for r in results]

    fig, ax = plt.subplots(figsize=(10, 7))

    # All points (faded if not localized)
    ax.scatter(omegas_all, h_all, color="lightgray", alpha=0.5, s=30,
               label="unstable (delocalized)")
    ax.scatter(omegas_localized, h_localized, color="C0", s=60,
               label="localized (chaoiton)")

    # Overlay Werbos's lepton targets, scaled to our calibration
    # Werbos: m_observed_MeV / m_electron_MeV = H_predicted / H_at_omega_1
    res_electron = next((r for r in results if abs(r["omega"] - 1.0) < 0.01), None)
    if res_electron and res_electron["localized"]:
        H_electron = res_electron["H"]
        mass_to_H_scale = 0.511 / H_electron   # MeV per code H unit
        for name, omega_target, mass_obs in LEPTON_TARGETS:
            ax.axvline(omega_target, color="red", linestyle=":", alpha=0.5)
            ax.text(omega_target, ax.get_ylim()[1] * 0.95, f"{name}\nm={mass_obs} MeV",
                    rotation=0, ha="center", va="top", fontsize=9, color="red")

        # Fit m ~ omega^p on the localized points
        if len(omegas_localized) >= 3:
            log_omega = np.log(np.array(omegas_localized))
            log_h = np.log(np.array(h_localized))
            slope, intercept = np.polyfit(log_omega, log_h, 1)
            print(f"\n  Mass scaling fit: m ~ omega^{slope:.3f}")
            print(f"  Werbos's published exponent: 2.22")
            x_fit = np.linspace(min(omegas_localized), max(omegas_localized), 100)
            y_fit = np.exp(intercept) * x_fit ** slope
            ax.plot(x_fit, y_fit, color="C2", linestyle="--",
                    label=f"fit: m ~ omega^{slope:.2f} (Werbos: omega^2.22)")

    ax.set_xlabel("omega (code units)")
    ax.set_ylabel("H (code units, proportional to mass)")
    ax.set_title(f"Ouroboros mass-frequency scan, g={G_COUPLING}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(loc="lower right")
    ax.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"  saved: {out_path}")


# ---------------------------------------------------------------------------
# Decision-grade summary
# ---------------------------------------------------------------------------

def evaluate_decision(scan_results, gate_0_passed):
    """Compare scan results to Werbos's lepton targets and emit M6 decision."""
    print("\n=== M6 DECISION ===\n")

    if not gate_0_passed:
        print("  Gate 0 (electron calibration) FAILED.")
        print("  Decision: NO-GO at current ODE form — refine before continuing.")
        return "NO-GO (gate 0)"

    res_electron = next((r for r in scan_results if abs(r["omega"] - 1.0) < 0.01), None)
    if res_electron is None or not res_electron["localized"]:
        print("  Electron baseline not localized.")
        return "NO-GO (no electron baseline)"

    H_electron = res_electron["H"]
    mass_scale = 0.511 / H_electron     # MeV per code unit

    print(f"  Mass scale: {mass_scale:.6f} MeV per code H unit")
    print(f"  (calibrated by setting H(omega=1.0) = m_electron = 0.511 MeV)\n")

    gaps = []
    for name, omega_target, mass_obs in LEPTON_TARGETS:
        res = next((r for r in scan_results if abs(r["omega"] - omega_target) < 0.1), None)
        if res is None:
            print(f"  {name} (omega={omega_target}): not in scan range")
            continue
        if not res["localized"]:
            print(f"  {name} (omega={omega_target}): NOT localized")
            continue

        m_predicted = res["H"] * mass_scale
        gap = abs(m_predicted - mass_obs) / mass_obs
        gaps.append(gap)
        print(f"  {name:8} (omega={omega_target:5.1f}): predicted = {m_predicted:8.2f} MeV, "
              f"observed = {mass_obs:8.2f} MeV, gap = {gap * 100:5.2f}%")

    if not gaps[1:]:    # exclude electron (calibration anchor)
        print("\n  Decision: INSUFFICIENT (need at least muon to be localized)")
        return "INSUFFICIENT"

    max_lepton_gap = max(gaps[1:])     # exclude electron
    print(f"\n  Max non-electron lepton gap: {max_lepton_gap * 100:.2f}%")

    if max_lepton_gap < 0.05:
        print("  Decision: GO (within 5% on all leptons)")
        return "GO"
    elif max_lepton_gap < 0.10:
        print("  Decision: CONDITIONAL GO (electron+muon within 5-10%; investigate tau)")
        return "CONDITIONAL GO"
    else:
        print("  Decision: NO-GO (gaps exceed 10% — framework's numerics don't reproduce)")
        return "NO-GO (lepton gaps)"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    here = Path(__file__).parent
    plots_dir = here.parent / "plots"
    plots_dir.mkdir(exist_ok=True)

    # Gate 0 — calibration check at omega=1.0
    gate_0_passed, res_electron = gate_0_calibration()

    if res_electron is not None:
        plot_calibration(res_electron, plots_dir / "m6_0_calibration_check.png")

    # Omega scan — even if Gate 0 fails, run the scan to see qualitative behavior
    # (Werbos scanned 1.0 to 80; we use targets + nearby values)
    omegas = np.array([1.0, 2.0, 5.0, 8.0, 10.0, 11.0, 12.0, 13.0, 15.0,
                       20.0, 25.0, 30.0, 35.0, 40.7, 50.0, 60.0, 70.0, 80.0])
    scan_results = omega_scan(omegas)

    plot_mass_frequency_scan(scan_results, plots_dir / "m6_0_mass_frequency_scan.png")

    # Decision
    decision = evaluate_decision(scan_results, gate_0_passed)

    # Persist results for downstream analysis
    summary = {
        "g_coupling": G_COUPLING,
        "lambda": LAMBDA,
        "gate_0_passed": gate_0_passed,
        "calibration_HQ_target": CALIBRATION_TARGET_HQ,
        "calibration_HQ_computed": res_electron["HQ"] if res_electron else None,
        "decision": decision,
        "scan": [{k: v for k, v in r.items()
                  if k not in ("r", "alpha", "dalpha", "gamma", "dgamma")}
                 for r in scan_results],
    }
    summary_path = plots_dir / "m6_0_summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n  summary written: {summary_path}")

    return decision


if __name__ == "__main__":
    main()
