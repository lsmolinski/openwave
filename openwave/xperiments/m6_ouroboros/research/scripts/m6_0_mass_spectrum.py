"""M6.0 mass spectrum — full omega scan with calibrated ODE form.

Uses the variant identified by m6_0_iterations.py: ODE A4 + H1 + Q2.
This combination reproduces Werbos's electron calibration H/Q = 1.6969
to 0.30% at omega=1.0, g=1.0625.

Variant rationale:
  ODE A4 — quarter coupling. Werbos's Lagrangian L = -F·F (no 1/4 factor)
           gives source-on-RHS divided by 4 in the field equation.
  H1     — Hamiltonian density without 1/2 factor. The 4× Lagrangian
           normalization gives 4× the standard Maxwell energy density,
           equivalent to dropping the 1/2 from the standard form.
  Q2     — Both A and J contribute to U(1) Noether charge. Both fields
           share the e^(-i omega t) oscillation, so both transform under
           the same U(1) and both contribute to the Noether current.

Result is the chaoiton mass-frequency relationship — checks the lepton
spectrum claims from Werbos's Particle Spectrum paper.
"""

from pathlib import Path
import json

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# Calibrated parameters
G_COUPLING = 1.0625
A_INIT = 0.1
B_INIT = 0.1
R_MIN = 0.02
R_MAX = 10.0
N_GRID = 800
RTOL = 1.0e-10
ATOL = 1.0e-12
LOCALIZATION_THRESHOLD = 0.15
CALIBRATION_TARGET_HQ = 1.6969
CALIBRATION_TOLERANCE = 0.05

LEPTON_TARGETS = [
    ("electron", 1.0, 0.511),
    ("muon", 11.0, 105.7),
    ("pion+", 13.0, 139.6),
    ("tau", 40.7, 1777.0),
]


# ---------------------------------------------------------------------------
# Calibrated ODE — A4 (quarter coupling)
# ---------------------------------------------------------------------------

def radial_ode(r, y, omega, g=G_COUPLING):
    """Calibrated radial ODE with quarter coupling.

    State: y = [alpha, alpha', gamma, gamma']

    alpha'' + (2/r) alpha' + omega^2 alpha + gamma/4 = 0        (Maxwell-like, J/4 source)
    gamma'' + (2/r) gamma' + omega^2 gamma - 4g gamma^3 + alpha/4 = 0  (Klein-Gordon-like)

    Quarter factor comes from Werbos's Lagrangian L = -F·F + J·A - f(J·J)
    where the standard EM normalization would give -F·F / 4.  His 4x
    normalization sends the source-of-A term to J/4 in the field equation.
    """
    a, da, c, dc = y
    inv_r = 1.0 / r
    d2a = -2*inv_r*da - omega*omega*a - c/4
    d2c = -2*inv_r*dc - omega*omega*c + 4*g*c*c*c - a/4
    return [da, d2a, dc, d2c]


def hamiltonian_density(omega, a, da, c, dc, g=G_COUPLING):
    """H1 — no 1/2 factor.

    h(r) = alpha'^2 + omega^2 alpha^2     [A-field kinetic + oscillation]
         + gamma'^2 + omega^2 gamma^2     [J-field kinetic + oscillation]
         + g gamma^4                      [f(J·J) potential]
         - alpha gamma                    [J·A coupling, sign from Lagrangian]
    """
    return (da*da + omega*omega*a*a
            + dc*dc + omega*omega*c*c
            + g*c*c*c*c
            - a*c)


def charge_density(omega, a, c):
    """Q2 — both fields contribute to U(1) Noether charge."""
    return omega * (a*a + c*c)


# ---------------------------------------------------------------------------
# Integration + observables
# ---------------------------------------------------------------------------

def run_chaoiton(omega, g=G_COUPLING, verbose=False):
    y0 = [A_INIT, 0.0, B_INIT, 0.0]
    r_eval = np.linspace(R_MIN, R_MAX, N_GRID)

    sol = solve_ivp(
        fun=lambda r, y: radial_ode(r, y, omega, g),
        t_span=(R_MIN, R_MAX),
        y0=y0,
        t_eval=r_eval,
        method="RK45",
        rtol=RTOL,
        atol=ATOL,
    )
    if not sol.success:
        return None

    r = sol.t
    a, da, c, dc = sol.y[0], sol.y[1], sol.y[2], sol.y[3]
    tail = abs(a[-1]) + abs(c[-1])
    localized = tail < LOCALIZATION_THRESHOLD

    h_dens = hamiltonian_density(omega, a, da, c, dc, g)
    q_dens = charge_density(omega, a, c)
    r2 = r * r
    H = 4.0 * np.pi * np.trapezoid(h_dens * r2, r)
    Q = 4.0 * np.pi * np.trapezoid(q_dens * r2, r)
    hq = H / Q if abs(Q) > 1e-12 else float("nan")

    if verbose:
        loc_str = "Y" if localized else "n"
        print(f"  omega = {omega:7.3f}: H = {H:9.4f}, Q = {Q:9.4f}, "
              f"H/Q = {hq:8.4f}, tail = {tail:7.4f}, localized = {loc_str}")

    return {
        "omega": omega,
        "r": r, "alpha": a, "gamma": c,
        "H": float(H), "Q": float(Q), "HQ": float(hq),
        "tail": float(tail), "localized": bool(localized),
    }


# ---------------------------------------------------------------------------
# Mass spectrum
# ---------------------------------------------------------------------------

def main():
    here = Path(__file__).parent
    plots_dir = here.parent / "plots"
    plots_dir.mkdir(exist_ok=True)

    # Gate 0
    print("=" * 78)
    print("GATE 0 — electron calibration check")
    print("=" * 78)
    print(f"  Target: H/Q = {CALIBRATION_TARGET_HQ} at (omega=1.0, g={G_COUPLING})")
    print(f"  Tolerance: {CALIBRATION_TOLERANCE*100:.1f}%\n")
    res_electron = run_chaoiton(1.0, verbose=True)
    gap_0 = abs(res_electron["HQ"] - CALIBRATION_TARGET_HQ) / CALIBRATION_TARGET_HQ
    print(f"\n  Gap: {gap_0*100:.2f}%   Verdict: {'PASS' if gap_0 < CALIBRATION_TOLERANCE else 'FAIL'}")

    # Mass-frequency scan
    print(f"\n{'=' * 78}")
    print("MASS-FREQUENCY SCAN")
    print(f"{'=' * 78}\n")
    omegas = sorted(set([1.0, 2.0, 5.0, 8.0, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5,
                         13.0, 13.5, 14.0, 15.0, 20.0, 25.0, 30.0, 35.0,
                         40.0, 40.7, 41.0, 45.0, 50.0, 60.0]))
    results = [run_chaoiton(omega, verbose=True) for omega in omegas]
    results = [r for r in results if r is not None]

    # Mass-scaling fit
    localized = [r for r in results if r["localized"]]
    omegas_loc = np.array([r["omega"] for r in localized])
    H_loc = np.array([r["H"] for r in localized])

    log_omega = np.log(omegas_loc)
    log_H = np.log(H_loc)
    slope, intercept = np.polyfit(log_omega, log_H, 1)

    print(f"\n{'=' * 78}")
    print(f"MASS-FREQUENCY SCALING")
    print(f"{'=' * 78}")
    print(f"\n  Fit: H ~ omega^{slope:.3f}")
    print(f"  Werbos's published: omega^2.22")
    print(f"  Scaling gap: {abs(slope - 2.22) / 2.22 * 100:.2f}%")

    # Lepton mass comparison
    print(f"\n{'=' * 78}")
    print("LEPTON MASS COMPARISON")
    print(f"{'=' * 78}")
    H_electron = next((r["H"] for r in results if abs(r["omega"] - 1.0) < 0.01), None)
    if H_electron is None:
        print("  Electron baseline missing — cannot compute lepton masses")
        return

    mass_scale = 0.511 / H_electron  # MeV per code H unit
    print(f"\n  Mass scale: {mass_scale:.6f} MeV per code H unit\n")
    print(f"  {'particle':>10}  {'omega':>8}  {'predicted (MeV)':>16}  "
          f"{'observed':>10}  {'gap':>8}")
    print("  " + "-" * 60)

    gaps = []
    for name, omega_target, mass_obs in LEPTON_TARGETS:
        res = next((r for r in results if abs(r["omega"] - omega_target) < 0.05), None)
        if res is None:
            print(f"  {name:>10}  {omega_target:8.2f}  not in scan")
            continue
        if not res["localized"]:
            print(f"  {name:>10}  {omega_target:8.2f}  NOT LOCALIZED")
            continue
        m_pred = res["H"] * mass_scale
        gap = abs(m_pred - mass_obs) / mass_obs
        gaps.append(gap)
        print(f"  {name:>10}  {omega_target:8.2f}  {m_pred:16.2f}  "
              f"{mass_obs:10.2f}  {gap*100:6.2f}%")

    # Decision
    print(f"\n{'=' * 78}")
    print("M6 DECISION")
    print(f"{'=' * 78}\n")
    if gap_0 > 0.05:
        decision = "NO-GO (Gate 0)"
    elif not gaps[1:]:
        decision = "INSUFFICIENT"
    else:
        max_lepton_gap = max(gaps[1:])
        if max_lepton_gap < 0.05:
            decision = "GO"
        elif max_lepton_gap < 0.10:
            decision = "CONDITIONAL GO"
        else:
            decision = "NO-GO (lepton gaps)"
    print(f"  Decision: {decision}\n")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Calibration profiles
    ax1.plot(res_electron["r"], res_electron["alpha"], label="alpha (A-field)", color="C0")
    ax1.plot(res_electron["r"], res_electron["gamma"], label="gamma (J-field)", color="C1")
    ax1.set_xlabel("r (code units)")
    ax1.set_ylabel("amplitude")
    ax1.set_title(f"Calibrated chaoiton, omega=1.0\nH/Q = {res_electron['HQ']:.4f} "
                  f"(target {CALIBRATION_TARGET_HQ:.4f}, gap {gap_0*100:.2f}%)")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Mass-frequency
    omegas_all = [r["omega"] for r in results]
    H_all = [r["H"] for r in results]
    ax2.scatter(omegas_all, H_all, color="lightgray", alpha=0.5, s=40,
                label="all points")
    ax2.scatter(omegas_loc, H_loc, color="C0", s=80, label="localized")
    for name, omega_t, mass_obs in LEPTON_TARGETS:
        ax2.axvline(omega_t, color="red", linestyle=":", alpha=0.5)
        ax2.text(omega_t, ax2.get_ylim()[1]*0.92 if ax2.get_ylim()[1] > 0 else 1,
                 name, ha="center", fontsize=9, color="red")
    x_fit = np.linspace(min(omegas_loc), max(omegas_loc), 100)
    y_fit = np.exp(intercept) * x_fit ** slope
    ax2.plot(x_fit, y_fit, color="C2", linestyle="--",
             label=f"fit: H ~ omega^{slope:.3f}  (Werbos: omega^2.22)")
    ax2.set_xlabel("omega (code units)")
    ax2.set_ylabel("H (code units, proportional to mass)")
    ax2.set_title(f"Mass-frequency scan, g={G_COUPLING}\nDecision: {decision}")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.legend(loc="lower right")
    ax2.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    out_path = plots_dir / "m6_0_calibrated_scan.png"
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"  Plot saved: {out_path}")

    # JSON summary
    summary = {
        "ode_variant": "A4 + H1 + Q2",
        "g_coupling": G_COUPLING,
        "gate_0_HQ_target": CALIBRATION_TARGET_HQ,
        "gate_0_HQ_computed": res_electron["HQ"],
        "gate_0_gap_pct": gap_0 * 100,
        "gate_0_passed": gap_0 < CALIBRATION_TOLERANCE,
        "mass_scaling_exponent_fit": float(slope),
        "mass_scaling_exponent_werbos": 2.22,
        "mass_scaling_exponent_gap_pct": abs(slope - 2.22) / 2.22 * 100,
        "decision": decision,
        "scan": [{k: v for k, v in r.items()
                  if k not in ("r", "alpha", "gamma")}
                 for r in results],
    }
    summary_path = plots_dir / "m6_0_calibrated_summary.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)
    print(f"  JSON summary saved: {summary_path}")


if __name__ == "__main__":
    main()
