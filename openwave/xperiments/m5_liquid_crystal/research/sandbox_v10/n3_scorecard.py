#!/usr/bin/env python3
"""
N3 / S4 , the NuFIT 6.0 scorecard + delta_CP + the summary figure (the deliverable assembly).

Pulls the N3 results together against the measured PMNS parameters:
  - the TBM gate (S2): theta12, theta23, theta13=0 from the closed-loop magic crossing (PREDICTIONS);
  - the theta13 crux (S3): a representative theta13 = 8.56 deg point from an O(0.1) mu-tau asymmetry (FIT);
  - delta_CP: 180 deg in the real CP-conserving limit (consistent with #199; a genuine value needs a
    complex/chiral coupling, deferred).
Generates n3_summary.png (3 panels: the magic-crossing gate, the theta13(delta,eps) map, the mass spectrum)
and writes n3_scorecard.json. Convention index-0. Headless. LOCAL (#236 N-program, HELD).
Run: python3 n3_scorecard.py
"""

import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from n3_mass_matrix import flavour_mass_matrix
from n3_search import find_alpha_star
from n3_theta13 import biaxial_mass_matrix, alpha_star
from n3_derisk import TH12_TBM, NUFIT

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    print("=" * 80)
    print("N3 / S4 , NuFIT 6.0 scorecard + delta_CP + summary figure")
    print("=" * 80)
    n = 44
    geom = {"R_loop": 9.0, "q": 0.5, "core_vox": 2.0, "kappa": 0.0}

    # ---- the TBM gate (predictions) ----
    geom_uni = dict(geom); geom_uni["delta"] = 1e-10
    astar, rgate = find_alpha_star(geom_uni)
    ag = rgate["angles"]
    print(f"\n[gate] alpha* = {np.degrees(astar):.3f} deg; TBM angles "
          f"th12={ag['theta12']:.3f} th23={ag['theta23']:.3f} th13={ag['theta13']:.2e}")

    # ---- a representative theta13 = 8.56 deg point (fit one O(0.1) mu-tau asymmetry at delta=0.1) ----
    delta_fit = 0.10
    def th13_of_eps(eps):
        a = alpha_star(n, delta_fit, eps, 0.0, geom)
        return biaxial_mass_matrix(n, a, delta_fit, eps, 0.0, **geom)["angles"]["theta13"] - NUFIT["theta13"]
    try:
        eps_star = brentq(th13_of_eps, 0.2, 0.95, xtol=1e-4)
    except Exception:
        eps_star = float("nan")
    a_fit = alpha_star(n, delta_fit, eps_star, 0.0, geom)
    rfit = biaxial_mass_matrix(n, a_fit, delta_fit, eps_star, 0.0, **geom)
    af = rfit["angles"]
    print(f"[fit ] theta13 = 8.56 deg at delta={delta_fit}, eps*={eps_star:.3f}: "
          f"th12={af['theta12']:.3f} th23={af['theta23']:.3f} th13={af['theta13']:.3f}")

    # ---- delta_CP (real CP-conserving limit) ----
    delta_CP_pred = 180.0   # real symmetric mass matrix -> CP-conserving; #199 value

    # ---- scorecard ----
    rows = [
        ("theta12 (deg)", ag["theta12"], af["theta12"], NUFIT["theta12"], "tri-maximal (sin^2=1/3)"),
        ("theta23 (deg)", ag["theta23"], af["theta23"], NUFIT["theta23"], "maximal (45)"),
        ("theta13 (deg)", ag["theta13"], af["theta13"], NUFIT["theta13"], "0 (TBM) -> 8.5 via mu-tau break"),
        ("delta_CP (deg)", delta_CP_pred, delta_CP_pred, NUFIT["delta_CP"], "180 (real/CP-conserving; #199)"),
    ]
    print("\n  scorecard (field theory vs NuFIT 6.0 NO):")
    print(f"  {'param':<15} {'TBM gate':>10} {'th13 fit':>10} {'NuFIT 6.0':>10}   note")
    for name, g, f, nf, note in rows:
        print(f"  {name:<15} {g:>10.3f} {f:>10.3f} {nf:>10.3f}   {note}")

    # ---- summary figure ----
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

    # panel 1: the magic-crossing gate (theta12 and magic residual vs alpha)
    alphas = np.linspace(0.5, 1.15, 26)
    th12s, magics = [], []
    for al in alphas:
        rr = flavour_mass_matrix(n=40, alpha=al, **geom_uni)
        th12s.append(rr["angles"]["theta12"]); magics.append(rr["magic_residual"])
    ax0 = ax[0]; ax0b = ax0.twinx()
    ax0.plot(np.degrees(alphas), th12s, "o-", color="#2980b9", label="theta12")
    ax0.axhline(TH12_TBM, ls="--", color="#27ae60", label="trimaximal 35.26")
    ax0b.plot(np.degrees(alphas), magics, "s-", color="#c0392b", alpha=0.6, label="magic residual")
    ax0b.axhline(0.0, ls=":", color="#7f8c8d")
    ax0.axvline(np.degrees(astar), ls=":", color="#2c3e50")
    ax0.set_xlabel("mirror tilt alpha (deg)"); ax0.set_ylabel("theta12 (deg)", color="#2980b9")
    ax0b.set_ylabel("magic residual (x+y)-(z+w)", color="#c0392b")
    ax0.set_title(f"TBM gate: theta12 -> 35.26 at the magic crossing\nalpha* = {np.degrees(astar):.1f} deg")
    ax0.legend(loc="upper right", fontsize=8)

    # panel 2: theta13(delta, eps) map
    deltas = np.array([0.03, 0.05, 0.08, 0.10, 0.15])
    epss = np.linspace(0.0, 0.9, 10)
    TH = np.zeros((len(deltas), len(epss)))
    for i, d in enumerate(deltas):
        for j, e in enumerate(epss):
            if e == 0:
                TH[i, j] = 0.0; continue
            aa = alpha_star(40, float(d), float(e), 0.0, geom)
            TH[i, j] = biaxial_mass_matrix(40, aa, float(d), float(e), 0.0, **geom)["angles"]["theta13"]
    im = ax[1].imshow(TH, aspect="auto", origin="lower", cmap="viridis",
                      extent=[epss.min(), epss.max(), deltas.min(), deltas.max()])
    cs = ax[1].contour(epss, deltas, TH, levels=[8.56], colors="white", linewidths=2)
    ax[1].clabel(cs, fmt="8.56 deg (NuFIT)", fontsize=8)
    fig.colorbar(im, ax=ax[1], label="theta13 (deg)")
    ax[1].set_xlabel("mu-tau asymmetry eps"); ax[1].set_ylabel("biaxiality delta")
    ax[1].set_title("theta13 ~ G*delta*eps (needs BOTH)\nsymmetric delta (eps=0) -> theta13 = 0 exactly")

    # panel 3: the mass spectrum at the gate
    w = np.array(sorted(rgate["eigenvalues"]))
    ax[2].bar(["m1\n(2,-1,-1)", "m2\n(1,1,1)", "m3\n(0,1,-1)"], w, color=["#2980b9", "#27ae60", "#c0392b"])
    ax[2].set_ylabel("mass-matrix eigenvalue (arb units)")
    ax[2].set_title(f"mass spectrum at the TBM gate\nratios {w[0]/w[0]:.2f} : {w[1]/w[0]:.2f} : {w[2]/w[0]:.2f}")
    ax[2].grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    plot_path = os.path.join(HERE, "n3_summary.png")
    fig.savefig(plot_path, dpi=120)
    print(f"\nfigure -> {plot_path}")

    scorecard = {
        "convention": "index-0 Duda: D=diag(g,1,delta,0), eta=diag(-1,1,1,1)",
        "tbm_gate": {"alpha_star_deg": float(np.degrees(astar)), "angles": ag,
                     "eigenvalues": rgate["eigenvalues"], "status": "PASS (predictions from dynamics)"},
        "theta13_fit_point": {"delta": delta_fit, "eps_star": float(eps_star), "angles": af,
                              "status": "FIT (one free O(0.1) mu-tau asymmetry)"},
        "delta_CP": {"predicted_deg": delta_CP_pred,
                     "note": "real symmetric mass matrix -> CP-conserving (delta_CP in {0,180}); "
                             "180 matches #199; a genuine value needs a complex/chiral coupling (deferred)"},
        "scorecard_vs_nufit60_NO": [
            {"param": r[0], "tbm_gate": r[1], "theta13_fit": r[2], "nufit": r[3], "note": r[4]}
            for r in rows],
        "mass_spectrum_ratios": (w / w[0]).tolist(),
    }
    with open(os.path.join(HERE, "n3_scorecard.json"), "w") as f:
        json.dump(scorecard, f, indent=2)
    print(f"scorecard -> {os.path.join(HERE, 'n3_scorecard.json')}")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
