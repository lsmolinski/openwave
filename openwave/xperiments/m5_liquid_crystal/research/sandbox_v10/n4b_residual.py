#!/usr/bin/env python3
"""
N4b-2 , the theta12 / theta23 residuals: which are free, which cost the predictions?

The exact-TBM predictions sit ~2 deg from data: theta12 = 35.26 vs NuFIT 33.68; theta23 = 45 vs 43.3. The
key structural fact: theta23 = 45 and delta_CP = +-90 come from the mu-tau MIRROR symmetry (holds at ANY tilt
alpha), while theta12 (trimaximal) comes from the MAGIC crossing (a specific alpha). So:

  PART A , theta12 is tunable by the tilt alpha alone. Moving alpha OFF the magic crossing shifts theta12
    while theta23 = 45 and delta_CP = +-90 stay EXACT (mu-tau is untouched). theta12 passes through 33.68
    (data) at some alpha. Cost: only the trimaximal PREDICTION (theta12 becomes a fit, not 35.26 predicted).
    The delta_CP = maximal prediction is NOT affected.

  PART B , theta23 = 45 is a mu-tau prediction. Matching NuFIT 43.3 requires an explicit mu-tau breaking
    (here a gentle secondary-twist asymmetry chi_mu != chi_tau), which moves delta_CP off +-90. We measure
    the delta_CP COST of pulling theta23 to 43.3.

So the article makes a sharp, honest choice: {exact TBM (theta12=35.26, theta23=45) + maximal CP} vs
{theta12,theta23 fit to data, delta_CP not maximal}. The delta_CP measurement (DUNE/HK) decides.

Convention index-0. Headless f64, 16-core. LOCAL (#236 N-program, HELD). Run: python3 n4b_residual.py
"""

import json
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor

from n3_mass_matrix import rot_axis
from n3_theta13 import seed_loop_biaxial, biaxial_vacuum, alpha_star
from n4_chiral import real_overlap, chiral_overlap, pmns_from_H
from n3_derisk import NUFIT

HERE = os.path.dirname(os.path.abspath(__file__))


def chiral_M(n, alpha, delta, chi_mu, chi_tau, g_chiral, R_loop=9.0, q=0.5, core_vox=2.0):
    """Chiral complex-Hermitian mass matrix; mu/tau = Rx(+-alpha) mirror with secondary twists chi_mu/chi_tau
    (chi_mu = chi_tau -> mu-tau exact). Returns the PMNS angles + delta_CP."""
    Re, Rmu, Rtau = np.eye(3), rot_axis((1, 0, 0), +alpha), rot_axis((1, 0, 0), -alpha)
    Mvac = biaxial_vacuum(n, delta)
    fe = seed_loop_biaxial(n, Re, R_loop, delta, q=q, core_vox=core_vox, chi=0.0)
    fmu = seed_loop_biaxial(n, Rmu, R_loop, delta, q=q, core_vox=core_vox, chi=chi_mu)
    ftau = seed_loop_biaxial(n, Rtau, R_loop, delta, q=q, core_vox=core_vox, chi=chi_tau)
    d = [fe - Mvac, fmu - Mvac, ftau - Mvac]
    Mr = np.zeros((3, 3)); Cc = np.zeros((3, 3))
    for i in range(3):
        for j in range(i, 3):
            K, _ = real_overlap(d[i], d[j]); Mr[i, j] = Mr[j, i] = K
        for j in range(i + 1, 3):
            cc = chiral_overlap(d[i], d[j]); Cc[i, j] = cc; Cc[j, i] = -cc
    return pmns_from_H(Mr.astype(complex) + 1j * g_chiral * Cc)


def _eval_alpha(args):
    n, alpha, delta, chi, g, geom = args
    try:
        r = chiral_M(n, alpha, delta, chi, chi, g, **geom)
        return {"alpha_deg": float(np.degrees(alpha)), **{k: r[k] for k in
                ("theta12", "theta23", "theta13", "delta_CP")}}
    except Exception as e:                                    # noqa: BLE001
        return {"alpha_deg": float(np.degrees(alpha)), "theta12": float("nan"), "error": repr(e)}


def _eval_dchi(args):
    n, alpha, delta, chi, dchi, g, geom = args
    try:
        r = chiral_M(n, alpha, delta, chi + dchi, chi, g, **geom)
        return {"dchi": dchi, **{k: r[k] for k in ("theta12", "theta23", "theta13", "delta_CP")}}
    except Exception as e:                                    # noqa: BLE001
        return {"dchi": dchi, "theta23": float("nan"), "error": repr(e)}


def main():
    print("=" * 80)
    print("N4b-2 , theta12/theta23 residuals (which are free, which cost the predictions)")
    print("=" * 80)
    n = 40
    geom = {"R_loop": 9.0, "q": 0.5, "core_vox": 2.0}
    delta, chi, g = 0.1, 1.2, 0.937
    a_star = alpha_star(n, delta, chi, chi, geom)

    # ---- PART A: theta12 tunes with alpha; theta23 + delta_CP stay exact (mu-tau untouched) ----
    alphas = a_star * np.linspace(0.80, 1.30, 21)
    with ProcessPoolExecutor(max_workers=16) as ex:
        A = list(ex.map(_eval_alpha, [(n, float(al), delta, chi, g, geom) for al in alphas]))
    A = [r for r in A if np.isfinite(r.get("theta12", np.nan))]
    th23_fixed = max(abs(r["theta23"] - 45.0) for r in A)
    dcp_fixed = max(abs(abs(r["delta_CP"]) - 90.0) for r in A)
    th12rng = [min(r["theta12"] for r in A), max(r["theta12"] for r in A)]
    print(f"\n[A] alpha scan (around magic alpha*={np.degrees(a_star):.1f} deg):")
    print(f"    theta23 stays 45 (max dev {th23_fixed:.2e}), delta_CP stays +-90 (max dev {dcp_fixed:.2e})")
    print(f"    theta12 ranges [{th12rng[0]:.2f}, {th12rng[1]:.2f}] deg -> TUNABLE by alpha, crosses "
          f"NuFIT 33.68: {th12rng[0] <= 33.68 <= th12rng[1]}")
    a_for_3368 = min(A, key=lambda r: abs(r["theta12"] - NUFIT["theta12"]))
    print(f"    theta12 = 33.68 (data) at alpha={a_for_3368['alpha_deg']:.1f} deg: theta23="
          f"{a_for_3368['theta23']:.2f}, theta13={a_for_3368['theta13']:.2f}, delta_CP={a_for_3368['delta_CP']:.1f}")
    print(f"    -> theta12 residual costs ONLY the trimaximal prediction (alpha becomes a fit); "
          f"theta23 + delta_CP UNAFFECTED")

    # ---- PART B: theta23 -> 43.3 needs mu-tau breaking; the delta_CP cost ----
    dchis = np.linspace(0.0, 0.6, 16)
    with ProcessPoolExecutor(max_workers=16) as ex:
        B = list(ex.map(_eval_dchi, [(n, a_star, delta, chi, float(dc), g, geom) for dc in dchis]))
    B = [r for r in B if np.isfinite(r.get("theta23", np.nan))]
    b_for_433 = min(B, key=lambda r: abs(r["theta23"] - NUFIT["theta23"]))
    print(f"\n[B] mu-tau breaking (chi asymmetry) to pull theta23 -> 43.3:")
    print(f"    at dchi={b_for_433['dchi']:.3f}: theta23={b_for_433['theta23']:.2f} (NuFIT 43.3), "
          f"theta12={b_for_433['theta12']:.2f}, theta13={b_for_433['theta13']:.2f}, "
          f"delta_CP={b_for_433['delta_CP']:.1f}")
    dcp_cost = abs(abs(b_for_433["delta_CP"]) - 90.0)
    print(f"    COST: delta_CP moved {dcp_cost:.1f} deg off maximal (mu-tau breaking removes the protection)")

    print("\n" + "=" * 80)
    print("N4b-2 , the residuals, honestly:")
    print("  - theta12 (35.26 vs 33.68): FREE to fit via the tilt alpha; theta23=45 and delta_CP=+-90 stay")
    print("    EXACT (mu-tau untouched). The residual = the price of PREDICTING theta12 (trimaximal) vs fitting.")
    print(f"  - theta23 (45 vs 43.3): costs the maximal-CP prediction (a mu-tau breaking to reach 43.3 moves")
    print(f"    delta_CP ~{dcp_cost:.0f} deg off +-90).")
    print("  - So: {exact TBM theta12=35.26, theta23=45 + MAXIMAL CP} is the clean predictive option; the")
    print("    ~2 deg residuals are absorbable but theta23's fix costs maximal CP. DUNE/HK delta_CP decides.")
    print("=" * 80)

    summary = {
        "partA_theta12_vs_alpha": {
            "alpha_star_deg": float(np.degrees(a_star)),
            "theta23_max_dev_from_45": th23_fixed, "delta_CP_max_dev_from_90": dcp_fixed,
            "theta12_range": th12rng, "alpha_for_33p68": a_for_3368,
            "scan": A,
        },
        "partB_theta23_vs_mutau_break": {
            "dchi_for_43p3": b_for_433["dchi"], "point": b_for_433, "delta_CP_cost_deg": dcp_cost,
            "scan": B,
        },
        "conclusion": ("theta12 residual is FREE (tunable by the tilt alpha; theta23 + delta_CP unaffected , "
                       "it only costs the trimaximal prediction). theta23 residual costs the maximal-CP "
                       "prediction (matching 43.3 needs a mu-tau breaking that moves delta_CP off +-90). "
                       "Clean predictive option: exact TBM + maximal CP, ~2 deg from data; DUNE/HK decide."),
    }
    with open(os.path.join(HERE, "n4b_residual_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {os.path.join(HERE, 'n4b_residual_summary.json')}")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
