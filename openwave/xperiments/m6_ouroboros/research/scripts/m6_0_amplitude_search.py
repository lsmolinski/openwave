"""M6.0 amplitude search — does varying initial amplitude per omega
reproduce Werbos's mass spectrum?

The fixed-amplitude scan (m6_0_mass_spectrum.py) showed H/Q ~ omega^1.2
instead of Werbos's claimed omega^2.22. The likely missing ingredient:
the f(J·J) = g·gamma^4 term is 4th-order in gamma, while kinetic and
oscillation terms are 2nd-order. Scaling amplitudes up at higher omega
makes the 4th-order term dominate, boosting H/Q faster.

This script tests that hypothesis: at each lepton's target omega, can we
find an initial amplitude that gives the lepton's m/q ratio = H/Q?

If yes → Werbos's spectrum holds, but he has a hidden 4th parameter
        (amplitude per omega), so the "3 parameters total" claim is misleading.
If no  → The mass-spectrum mechanism doesn't reproduce; M6 is genuinely
        a NO-GO on its headline claim.
"""

from pathlib import Path
import json

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib.pyplot as plt


# Calibrated parameters (from m6_0_mass_spectrum.py)
G_COUPLING = 1.0625
R_MIN = 0.02
R_MAX = 10.0
N_GRID = 800
RTOL = 1.0e-10
ATOL = 1.0e-12
LOCALIZATION_THRESHOLD = 0.15

# Calibration target = electron m/q ratio = 0.511 / 0.30282 = 1.6875
# Werbos's solution at (g=1.0625, omega=1.0, amp=0.1) gives H/Q = 1.6969
ELECTRON_HQ_TARGET = 1.6969
E_NAT = 0.30282           # elementary charge in natural units (sqrt(4*pi*alpha))

# Lepton spectrum targets: (name, omega, mass_MeV, target_HQ)
LEPTON_TARGETS = [
    ("electron", 1.0, 0.511, 0.511 / E_NAT),     # = 1.6875, ~= 1.6969 our cal
    ("muon", 11.0, 105.7, 105.7 / E_NAT),        # = 349.1
    ("pion+", 13.0, 139.6, 139.6 / E_NAT),       # = 461.0
    ("tau", 40.7, 1777.0, 1777.0 / E_NAT),       # = 5867.7
]


# ---------------------------------------------------------------------------
# Calibrated ODE — A4 (quarter coupling, see m6_0_iterations.py)
# ---------------------------------------------------------------------------

def radial_ode(r, y, omega, g=G_COUPLING):
    a, da, c, dc = y
    inv_r = 1.0 / r
    d2a = -2*inv_r*da - omega*omega*a - c/4
    d2c = -2*inv_r*dc - omega*omega*c + 4*g*c*c*c - a/4
    return [da, d2a, dc, d2c]


def integrate(omega, amplitude, g=G_COUPLING):
    """Integrate at given omega and initial amplitude (A_0 = B_0 = amplitude)."""
    y0 = [amplitude, 0.0, amplitude, 0.0]
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
    return sol


def compute_HQ(sol, omega, g=G_COUPLING):
    """Compute H, Q, H/Q using calibrated definitions (H1 + Q2)."""
    r = sol.t
    a, da, c, dc = sol.y[0], sol.y[1], sol.y[2], sol.y[3]
    tail = abs(a[-1]) + abs(c[-1])
    h_dens = da*da + omega*omega*a*a + dc*dc + omega*omega*c*c + g*c*c*c*c - a*c
    q_dens = omega * (a*a + c*c)
    r2 = r * r
    H = 4.0 * np.pi * np.trapezoid(h_dens * r2, r)
    Q = 4.0 * np.pi * np.trapezoid(q_dens * r2, r)
    hq = H / Q if abs(Q) > 1e-12 else float("nan")
    field_max = float(max(np.abs(a).max(), np.abs(c).max()))
    return H, Q, hq, tail, field_max


# ---------------------------------------------------------------------------
# Amplitude search: at given omega, find amp such that H/Q ≈ target
# ---------------------------------------------------------------------------

def hq_at_amplitude(amplitude, omega, g=G_COUPLING):
    """H/Q as a function of amplitude (for the root finder)."""
    sol = integrate(omega, amplitude, g)
    if sol is None:
        return float("nan")
    _, _, hq, _, _ = compute_HQ(sol, omega, g)
    return hq


def find_amplitude_for_HQ(omega, hq_target, g=G_COUPLING,
                          amp_min=0.001, amp_max=20.0, n_probe=80):
    """Find the smallest amplitude such that H/Q ≈ hq_target.

    First probes the (amp_min, amp_max) range to find a bracket where H/Q
    crosses hq_target, then uses brentq for precise root.
    """
    # Probe
    amps_probe = np.geomspace(amp_min, amp_max, n_probe)
    hqs = []
    for amp in amps_probe:
        hq = hq_at_amplitude(amp, omega, g)
        hqs.append(hq)
    hqs = np.array(hqs)

    # Find first bracket where H/Q crosses target
    sign = np.sign(hqs - hq_target)
    brackets = []
    for i in range(len(sign) - 1):
        if not np.isnan(sign[i]) and not np.isnan(sign[i+1]) and sign[i] * sign[i+1] < 0:
            brackets.append((amps_probe[i], amps_probe[i+1]))

    if not brackets:
        return None, None, hqs, amps_probe

    # Bisect on the first bracket
    amp_lo, amp_hi = brackets[0]
    try:
        amp_found = brentq(
            lambda a: hq_at_amplitude(a, omega, g) - hq_target,
            amp_lo, amp_hi, xtol=1e-6, maxiter=100
        )
    except (ValueError, RuntimeError):
        return None, None, hqs, amps_probe

    sol = integrate(omega, amp_found, g)
    if sol is None:
        return None, None, hqs, amps_probe
    H, Q, hq, tail, field_max = compute_HQ(sol, omega, g)

    return amp_found, {
        "amp": float(amp_found),
        "H": float(H), "Q": float(Q), "HQ": float(hq),
        "tail": float(tail), "field_max": float(field_max),
        "localized": bool(tail < LOCALIZATION_THRESHOLD),
    }, hqs, amps_probe


# ---------------------------------------------------------------------------
# Main: search amplitude for each lepton
# ---------------------------------------------------------------------------

def main():
    here = Path(__file__).parent
    plots_dir = here.parent / "plots"
    plots_dir.mkdir(exist_ok=True)

    print("=" * 78)
    print("AMPLITUDE SEARCH — can we hit Werbos's lepton H/Q at each omega?")
    print("=" * 78)
    print(f"  g = {G_COUPLING}, ODE = A4 (quarter coupling)")
    print(f"  Searching amp ∈ [0.001, 20.0] for each lepton's target H/Q\n")

    print(f"  {'lepton':>10}  {'omega':>6}  {'target H/Q':>10}  "
          f"{'found amp':>10}  {'achieved H/Q':>13}  "
          f"{'field max':>10}  {'tail':>8}  {'loc':>5}")
    print("  " + "-" * 95)

    results = []
    probes = []
    for name, omega, mass_mev, hq_target in LEPTON_TARGETS:
        amp, res, hqs, amps_probe = find_amplitude_for_HQ(omega, hq_target)
        probes.append({"name": name, "omega": omega, "amps": amps_probe.tolist(),
                       "hqs": hqs.tolist(), "target": hq_target})
        if res is None:
            print(f"  {name:>10}  {omega:6.2f}  {hq_target:10.2f}  "
                  f"{'NO BRACKET':>10}  {'-':>13}  {'-':>10}  {'-':>8}  {'-':>5}")
            results.append({"name": name, "omega": omega, "target_HQ": hq_target,
                            "amp_found": None, "result": None})
            continue
        loc_str = "Y" if res["localized"] else "n"
        print(f"  {name:>10}  {omega:6.2f}  {hq_target:10.2f}  "
              f"{res['amp']:10.4f}  {res['HQ']:13.4f}  "
              f"{res['field_max']:10.4f}  {res['tail']:8.4f}  {loc_str:>5}")
        results.append({"name": name, "omega": omega, "target_HQ": hq_target,
                        "mass_mev": mass_mev,
                        "amp_found": res["amp"], "result": res})

    # Amplitude-vs-omega scaling
    valid = [r for r in results if r["amp_found"] is not None]
    if len(valid) >= 3:
        omegas = np.array([r["omega"] for r in valid])
        amps = np.array([r["amp_found"] for r in valid])
        log_omega = np.log(omegas[omegas > 0])
        log_amp = np.log(amps[amps > 0])
        if len(log_omega) >= 2:
            slope, intercept = np.polyfit(log_omega, log_amp, 1)
            print(f"\n  Amplitude scaling: amp ~ omega^{slope:.3f}")
            print(f"  Intercept: amp(omega=1) ≈ {np.exp(intercept):.4f}")

    # Localization audit
    print(f"\n  Localization audit:")
    all_localized = all(r["result"]["localized"] for r in results if r["result"])
    print(f"    All lepton chaoitons localized? {'YES' if all_localized else 'NO'}")
    for r in results:
        if r["result"]:
            print(f"    {r['name']:>10}: tail = {r['result']['tail']:.4f}  "
                  f"field_max = {r['result']['field_max']:.4f}")

    # Plot amplitude vs omega and H/Q vs amplitude curves
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # H/Q vs amplitude curves (for each omega)
    for p in probes:
        amps_arr = np.array(p["amps"])
        hqs_arr = np.array(p["hqs"])
        mask = ~np.isnan(hqs_arr)
        ax1.loglog(amps_arr[mask], hqs_arr[mask], label=f"{p['name']} (omega={p['omega']})")
        ax1.axhline(p["target"], color="gray", linestyle=":", alpha=0.5)
    ax1.set_xlabel("initial amplitude (A_0 = B_0)")
    ax1.set_ylabel("H/Q")
    ax1.set_title("H/Q vs amplitude at each lepton omega")
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3, which="both")

    # Amplitude vs omega scaling
    if len(valid) >= 2:
        omegas = np.array([r["omega"] for r in valid])
        amps = np.array([r["amp_found"] for r in valid])
        ax2.loglog(omegas, amps, "o-", color="C0", markersize=10,
                   label="found amplitudes")
        # Show Werbos's nominal A_0=0.1 as a horizontal line
        ax2.axhline(0.1, color="red", linestyle="--", alpha=0.7,
                    label="Werbos nominal A_0=0.1")
        for r in valid:
            ax2.annotate(r["name"], (r["omega"], r["amp_found"]),
                         textcoords="offset points", xytext=(10, 5), fontsize=9)
        ax2.set_xlabel("omega")
        ax2.set_ylabel("required initial amplitude")
        ax2.set_title("Hidden 4th parameter — amplitude per omega?")
        ax2.legend()
        ax2.grid(True, alpha=0.3, which="both")

    plt.tight_layout()
    out_path = plots_dir / "m6_0_amplitude_search.png"
    plt.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"\n  Plot saved: {out_path}")

    # JSON summary
    summary = {
        "g_coupling": G_COUPLING,
        "leptons": [
            {"name": r["name"], "omega": r["omega"],
             "target_HQ": r["target_HQ"],
             "amp_found": r.get("amp_found"),
             "result": r.get("result")}
            for r in results
        ],
    }
    if len(valid) >= 2:
        summary["amp_scaling"] = {
            "exponent": float(slope) if len(log_omega) >= 2 else None,
            "intercept_amp": float(np.exp(intercept)) if len(log_omega) >= 2 else None,
        }
    summary_path = plots_dir / "m6_0_amplitude_search.json"
    with summary_path.open("w") as f:
        json.dump(summary, f, indent=2)
    print(f"  JSON saved: {summary_path}")

    # Verdict
    print(f"\n{'=' * 78}")
    print("VERDICT")
    print(f"{'=' * 78}\n")
    if all_localized and all(r["amp_found"] is not None for r in results):
        print("  Hypothesis VALIDATED: varying initial amplitude per omega reproduces")
        print("  the lepton spectrum. This is a HIDDEN 4th PARAMETER in Werbos's")
        print("  framework. His '3 parameters total' claim is misleading by 1.")
        print("\n  M6 still has technical merit (electron calibration is real), but the")
        print("  mass-spectrum claim is parameter-fitting in disguise.")
    elif all(r["amp_found"] is not None for r in results):
        print("  PARTIAL: amplitudes exist that match target H/Q, but localization")
        print("  may be questionable. Inspect tail / field_max columns.")
    else:
        print("  NOT REPRODUCIBLE: even with amplitude search, the lepton H/Q values")
        print("  are out of reach in the chaoiton localization regime.")
        print("\n  This is a stronger M6 NO-GO — the framework's headline claim does")
        print("  not reproduce under either fixed-amplitude or amplitude-search scans.")


if __name__ == "__main__":
    main()
