"""M6.0 dual-criteria search — find (H, Q) definitions matching BOTH:
  1. H/Q = 1.6969 at omega=1.0, g=1.0625 (electron calibration)
  2. H/Q scaling ~ omega^2.22 across the lepton spectrum (Werbos's claim)

The fixed-amplitude scan showed H/Q ~ omega^1.2 with the original Q
definition. Dropping omega from Q gave H/Q ~ omega^2.04 — close to
Werbos's 2.22 but still off. This script systematically scans more
H and Q variants looking for the combination that matches both
criteria.

Each variant of H and Q is tried in combination. Output is ranked by
"composite score" = 0.5*(calibration_gap%) + 0.5*(scaling_gap%) plus
the lepton mass predictions per variant.
"""

from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
import json


G_COUPLING = 1.0625
R_MIN = 0.02
R_MAX = 10.0
N_GRID = 800
RTOL = 1.0e-10
ATOL = 1.0e-12

TARGET_HQ_AT_OMEGA_1 = 1.6969
TARGET_SCALING_EXPONENT = 2.22
SCAN_OMEGAS = [1.0, 2.0, 5.0, 11.0, 13.0, 20.0, 40.7]

LEPTON_TARGETS = [
    ("electron", 1.0, 0.511),
    ("muon", 11.0, 105.7),
    ("pion+", 13.0, 139.6),
    ("tau", 40.7, 1777.0),
]


# ---------------------------------------------------------------------------
# Calibrated ODE (A4 - quarter coupling, kept fixed)
# ---------------------------------------------------------------------------

def radial_ode(r, y, omega, g=G_COUPLING):
    a, da, c, dc = y
    inv_r = 1.0 / r
    d2a = -2*inv_r*da - omega*omega*a - c/4
    d2c = -2*inv_r*dc - omega*omega*c + 4*g*c*c*c - a/4
    return [da, d2a, dc, d2c]


def integrate(omega, amplitude=0.1, g=G_COUPLING):
    y0 = [amplitude, 0.0, amplitude, 0.0]
    r_eval = np.linspace(R_MIN, R_MAX, N_GRID)
    sol = solve_ivp(
        fun=lambda r, y: radial_ode(r, y, omega, g),
        t_span=(R_MIN, R_MAX), y0=y0, t_eval=r_eval,
        method="RK45", rtol=RTOL, atol=ATOL,
    )
    return sol


# ---------------------------------------------------------------------------
# H density variants
# ---------------------------------------------------------------------------

def h_v1(omega, a, da, c, dc, g):
    """Original: gradient + ω² field + g γ^4 - coupling"""
    return da*da + omega*omega*a*a + dc*dc + omega*omega*c*c + g*c*c*c*c - a*c

def h_v2(omega, a, da, c, dc, g):
    """v1 with +coupling instead of -coupling"""
    return da*da + omega*omega*a*a + dc*dc + omega*omega*c*c + g*c*c*c*c + a*c

def h_v3(omega, a, da, c, dc, g):
    """No coupling term"""
    return da*da + omega*omega*a*a + dc*dc + omega*omega*c*c + g*c*c*c*c

def h_v4(omega, a, da, c, dc, g):
    """v1 with ω·coupling (time-mediated J·A)"""
    return da*da + omega*omega*a*a + dc*dc + omega*omega*c*c + g*c*c*c*c - omega*a*c

def h_v5(omega, a, da, c, dc, g):
    """Only oscillation (drop gradients) - assumes localization tightens fast enough"""
    return omega*omega*a*a + omega*omega*c*c + g*c*c*c*c

def h_v6(omega, a, da, c, dc, g):
    """v1 with f(s) where s = γ² (assume scalar dominates) but with ω² inside f"""
    return da*da + omega*omega*a*a + dc*dc + omega*omega*c*c + g*omega*omega*c*c*c*c - a*c

def h_v7(omega, a, da, c, dc, g):
    """Just gamma kinetic + mass + coupling (alpha is Coulomb-static, no kinetic)"""
    return dc*dc + omega*omega*c*c + g*c*c*c*c - a*c

def h_v8(omega, a, da, c, dc, g):
    """Squared coupling: -(a·c)² as a nonlinear term"""
    return da*da + omega*omega*a*a + dc*dc + omega*omega*c*c + g*c*c*c*c - a*c*a*c


H_VARIANTS = {
    "v1": h_v1, "v2": h_v2, "v3": h_v3, "v4": h_v4,
    "v5": h_v5, "v6": h_v6, "v7": h_v7, "v8": h_v8,
}


# ---------------------------------------------------------------------------
# Q density variants — focus on omega-scaling structure
# ---------------------------------------------------------------------------

def q_v1(omega, a, c, g):
    """ω·(α²+γ²) — original (Werbos-claim)"""
    return omega * (a*a + c*c)

def q_v2(omega, a, c, g):
    """α²+γ² — no omega factor"""
    return a*a + c*c

def q_v3(omega, a, c, g):
    """α² alone (only A contributes)"""
    return a*a

def q_v4(omega, a, c, g):
    """γ² alone (only J contributes)"""
    return c*c

def q_v5(omega, a, c, g):
    """ω·γ² (only J, with omega factor)"""
    return omega * c * c

def q_v6(omega, a, c, g):
    """(α+γ)² — sum of fields squared (cross-coupled)"""
    return (a + c)**2

def q_v7(omega, a, c, g):
    """α·γ (product, not sum of squares)"""
    return a * c

def q_v8(omega, a, c, g):
    """ω/(ω²+1) · (α²+γ²) — saturating omega dependence"""
    return omega / (omega*omega + 1) * (a*a + c*c)


Q_VARIANTS = {
    "qv1": q_v1, "qv2": q_v2, "qv3": q_v3, "qv4": q_v4,
    "qv5": q_v5, "qv6": q_v6, "qv7": q_v7, "qv8": q_v8,
}


# ---------------------------------------------------------------------------
# Evaluate a (H, Q) pair against both criteria
# ---------------------------------------------------------------------------

def evaluate_variant(h_fn, q_fn, name=""):
    hq_values = []
    H_values = []
    Q_values = []

    for omega in SCAN_OMEGAS:
        sol = integrate(omega)
        if not sol.success:
            return None
        r = sol.t
        a, da, c, dc = sol.y[0], sol.y[1], sol.y[2], sol.y[3]
        h_dens = h_fn(omega, a, da, c, dc, G_COUPLING)
        q_dens = q_fn(omega, a, c, G_COUPLING)
        r2 = r * r
        H = 4.0 * np.pi * np.trapezoid(h_dens * r2, r)
        Q = 4.0 * np.pi * np.trapezoid(q_dens * r2, r)
        if abs(Q) < 1e-12:
            return None
        hq = H / Q
        if not np.isfinite(hq):
            return None
        hq_values.append(hq)
        H_values.append(H)
        Q_values.append(Q)

    hq_values = np.array(hq_values)
    H_values = np.array(H_values)
    omegas = np.array(SCAN_OMEGAS)

    # Criterion 1: H/Q at omega=1
    hq_at_1 = hq_values[0]
    cal_gap = abs(hq_at_1 - TARGET_HQ_AT_OMEGA_1) / TARGET_HQ_AT_OMEGA_1

    # Criterion 2: log-log fit slope on H/Q vs omega
    pos_mask = (hq_values > 0) & (omegas > 0)
    if pos_mask.sum() < 3:
        return None
    slope, _ = np.polyfit(np.log(omegas[pos_mask]), np.log(hq_values[pos_mask]), 1)
    scaling_gap = abs(slope - TARGET_SCALING_EXPONENT) / TARGET_SCALING_EXPONENT

    # Lepton mass gaps using H/Q itself as proxy for m/q (Werbos's scale-free ratio)
    # m_predicted = H/Q * e_nat (since H/Q in code = m/q in natural units)
    e_nat = 0.30282
    mass_predictions = {}
    for name_t, om_t, mass_obs in LEPTON_TARGETS:
        idx = list(SCAN_OMEGAS).index(om_t) if om_t in SCAN_OMEGAS else None
        if idx is None:
            continue
        m_pred = hq_values[idx] * e_nat
        mass_gap = abs(m_pred - mass_obs) / mass_obs
        mass_predictions[name_t] = {"omega": om_t, "m_pred": float(m_pred),
                                    "m_obs": mass_obs, "gap_pct": mass_gap * 100}

    # Composite score: lower is better
    composite = 0.4 * cal_gap + 0.4 * scaling_gap + 0.2 * np.mean(
        [m["gap_pct"]/100 for m in mass_predictions.values()] or [1.0]
    )

    return {
        "name": name,
        "hq_at_1": float(hq_at_1),
        "cal_gap_pct": cal_gap * 100,
        "slope": float(slope),
        "scaling_gap_pct": scaling_gap * 100,
        "composite_score": float(composite),
        "leptons": mass_predictions,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    results = []
    for h_key, h_fn in H_VARIANTS.items():
        for q_key, q_fn in Q_VARIANTS.items():
            res = evaluate_variant(h_fn, q_fn, name=f"{h_key}+{q_key}")
            if res is not None:
                results.append(res)

    results.sort(key=lambda r: r["composite_score"])

    print("=" * 90)
    print("DUAL-CRITERIA SEARCH")
    print(f"  Target: H/Q={TARGET_HQ_AT_OMEGA_1} at ω=1, scaling ~ ω^{TARGET_SCALING_EXPONENT}")
    print("=" * 90)
    print(f"  {'variant':>10}  {'H/Q@ω=1':>9}  {'cal gap%':>9}  "
          f"{'slope':>7}  {'sc gap%':>8}  {'composite':>10}")
    print("  " + "-" * 78)
    for r in results[:30]:
        print(f"  {r['name']:>10}  {r['hq_at_1']:9.4f}  {r['cal_gap_pct']:9.2f}  "
              f"{r['slope']:7.3f}  {r['scaling_gap_pct']:8.2f}  "
              f"{r['composite_score']*100:10.2f}")

    print("\n" + "=" * 90)
    print("TOP 5 — lepton predictions")
    print("=" * 90)
    for r in results[:5]:
        print(f"\n  Variant: {r['name']}")
        print(f"    H/Q at ω=1: {r['hq_at_1']:.4f}  (target {TARGET_HQ_AT_OMEGA_1}, "
              f"gap {r['cal_gap_pct']:.2f}%)")
        print(f"    Scaling: ω^{r['slope']:.3f}  (target ω^{TARGET_SCALING_EXPONENT}, "
              f"gap {r['scaling_gap_pct']:.2f}%)")
        for name, lp in r["leptons"].items():
            print(f"    {name:>10} (ω={lp['omega']:.2f}): predicted {lp['m_pred']:.2f} MeV  "
                  f"obs {lp['m_obs']:.2f}  gap {lp['gap_pct']:.2f}%")

    # JSON
    here = Path(__file__).parent
    out_path = here.parent / "plots" / "m6_0_dual_criteria_search.json"
    with out_path.open("w") as f:
        json.dump({"results": results}, f, indent=2)
    print(f"\n  JSON written: {out_path}")


if __name__ == "__main__":
    main()
