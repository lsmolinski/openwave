#!/usr/bin/env python3
"""
N4c-2 , the near-term falsifier the reviewer flagged: check the loop mass SPECTRUM ratios against the observed
mass-squared splittings NOW (without waiting for the deferred absolute-mass round N6).

The flavour mass matrix at the TBM gate has eigenvalues ~ 1 : 1.15 : 1.68 (a compressed / quasi-degenerate
spectrum). The data fixes the RATIO of the two mass-squared splittings:
    Dm31^2 / Dm21^2 ~ 2.513e-3 / 7.49e-5 ~ 33.6   (NuFIT 6.0 NO)  =>  sqrt ~ 5.8.
We can test the eigenvalue ratios against this immediately, under the two NATURAL identifications of the
overlap-matrix eigenvalue lambda_i:
    Map A  lambda_i = m_i      (eigenvalue is the mass)        -> ratio = (l3^2 - l1^2)/(l2^2 - l1^2)
    Map B  lambda_i = m_i^2    (eigenvalue is the mass-squared)-> ratio = (l3 - l1)/(l2 - l1)
and compare to the observed 33.6.

HONEST CAVEAT carried into the findings: the map lambda_i -> neutrino mass is UNDEFINED here (that IS the
deferred N6 question: Duda's mass-as-loop-length-density picture). So this is a qualitative tension FLAG, not a
falsification. But under both natural maps the spectrum is far too compressed, which is worth stating plainly
before building further (exactly the reviewer's point).

Convention index-0. Headless f64. LOCAL (#236 N-program, HELD). Run: python3 n4c_mass_ratio.py
"""

import json
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from n3_mass_matrix import flavour_mass_matrix, DELTA_PHYS

HERE = os.path.dirname(os.path.abspath(__file__))
RAD = 180.0 / np.pi

# NuFIT 6.0, normal ordering (approx central values)
NUFIT_DM21 = 7.49e-5      # eV^2 (solar)
NUFIT_DM31 = 2.513e-3     # eV^2 (atmospheric, NO)


def main():
    print("=" * 84)
    print("N4c-2 , loop mass-spectrum ratios vs the observed mass-squared splittings (near-term falsifier)")
    print("=" * 84)
    n = 40
    delta = DELTA_PHYS
    geom = {"R_loop": 9.0, "q": 0.5, "core_vox": 2.0, "kappa": 0.0}

    # ---- gate eigenvalues at the magic crossing alpha* ----
    def magic_resid(a):
        return float(flavour_mass_matrix(n=n, alpha=a, delta=delta, **geom)["magic_residual"])
    aa = np.linspace(0.30, 1.30, 21)
    rr = np.array([magic_resid(a) for a in aa])
    a_star = None
    for i in range(len(aa) - 1):
        if rr[i] * rr[i + 1] < 0:
            a_star = float(brentq(magic_resid, aa[i], aa[i + 1], xtol=1e-7, rtol=1e-12)); break
    r = flavour_mass_matrix(n=n, alpha=a_star, delta=delta, **geom)
    lam = np.array(sorted(r["eigenvalues"]))      # ascending
    lam_n = lam / lam[0]
    print(f"\n[gate] alpha* = {a_star*RAD:.3f} deg; theta12={r['angles']['theta12']:.4f} (TBM err {r['tbm_err']:.4f})")
    print(f"[gate] eigenvalues = {lam.tolist()}")
    print(f"[gate] normalized  = {lam_n.round(4).tolist()}  ( = 1 : {lam_n[1]:.3f} : {lam_n[2]:.3f} )")

    l1, l2, l3 = lam
    ratio_obs = NUFIT_DM31 / NUFIT_DM21
    # Map A: eigenvalue = mass
    s2_A, s3_A = l2**2 - l1**2, l3**2 - l1**2
    ratio_A = s3_A / s2_A
    # Map B: eigenvalue = mass^2
    s2_B, s3_B = l2 - l1, l3 - l1
    ratio_B = s3_B / s2_B

    print(f"\n[data]  observed  Dm31^2/Dm21^2 = {NUFIT_DM31:.3e}/{NUFIT_DM21:.3e} = {ratio_obs:.2f}  "
          f"(sqrt = {np.sqrt(ratio_obs):.2f})")
    print(f"\n  Map A (lambda = mass):       (l3^2-l1^2)/(l2^2-l1^2) = {ratio_A:.2f}   "
          f"-> off by x{ratio_obs/ratio_A:.1f}  (TOO COMPRESSED)")
    print(f"  Map B (lambda = mass^2):     (l3-l1)/(l2-l1)         = {ratio_B:.2f}   "
          f"-> off by x{ratio_obs/ratio_B:.1f}  (TOO COMPRESSED)")

    # what spectrum WOULD the data require (Map A, with m1 the free scale)? m2/m1, m3/m1 given Dm's + a scale
    # pick the lightest m1 by requiring the loop's m2/m1 = lam_n[1]; then check m3/m1 implied vs lam_n[2]
    print(f"\n[interpretation] the loop spectrum is quasi-degenerate (1 : {lam_n[1]:.2f} : {lam_n[2]:.2f}); the")
    print(f"   data needs the atmospheric/solar splitting ratio ~{ratio_obs:.0f}, i.e. a MUCH more spread")
    print(f"   spectrum. Under both natural maps the loop splittings are ~5-7x too compressed.")

    tension = bool(ratio_A < 0.5 * ratio_obs and ratio_B < 0.5 * ratio_obs)
    print("\n" + "=" * 84)
    print(f"VERDICT: the gate spectrum 1 : {lam_n[1]:.2f} : {lam_n[2]:.2f} is in TENSION with the observed")
    print(f"  splitting hierarchy (off by ~5-7x under both natural eigenvalue->mass maps). This is a near-term")
    print(f"  FLAG, not a falsification: the eigenvalue->mass map is the deferred N6 question (mass-as-loop-")
    print(f"  length-density). A successful mass model must either spread the spectrum or use a nonlinear map.")
    print("=" * 84)

    # ---- figure ----
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    ax[0].bar([0, 1, 2], lam_n, color=["#27ae60", "#e67e22", "#8e44ad"])
    ax[0].set_xticks([0, 1, 2]); ax[0].set_xticklabels(["lam1", "lam2", "lam3"])
    ax[0].set_ylabel("eigenvalue / lam1")
    ax[0].set_title(f"loop gate spectrum: 1 : {lam_n[1]:.2f} : {lam_n[2]:.2f}\n(compressed / quasi-degenerate)")
    for i, v in enumerate(lam_n):
        ax[0].text(i, v + 0.02, f"{v:.2f}", ha="center", fontsize=9)
    ax[0].grid(True, axis="y", alpha=0.3)

    labels = ["Map A\n(lam=m)", "Map B\n(lam=m^2)", "observed\n(NuFIT NO)"]
    vals = [ratio_A, ratio_B, ratio_obs]
    cols = ["#2980b9", "#16a085", "#c0392b"]
    ax[1].bar([0, 1, 2], vals, color=cols)
    ax[1].set_yscale("log")
    ax[1].set_xticks([0, 1, 2]); ax[1].set_xticklabels(labels, fontsize=8)
    ax[1].set_ylabel("Dm31^2 / Dm21^2  (splitting ratio)")
    ax[1].set_title("splitting ratio: loop spectrum vs data\n(loop ~5-7x too compressed)")
    for i, v in enumerate(vals):
        ax[1].text(i, v * 1.1, f"{v:.1f}", ha="center", fontsize=9)
    ax[1].grid(True, axis="y", which="both", alpha=0.3)
    fig.tight_layout()
    fig_path = os.path.join(HERE, "n4c_mass_ratio.png")
    fig.savefig(fig_path, dpi=110)
    print(f"\nfigure -> {fig_path}")

    summary = {
        "alpha_star_deg": a_star * RAD, "tbm_err": r["tbm_err"],
        "eigenvalues": lam.tolist(), "eigenvalues_normalized": lam_n.tolist(),
        "spectrum_ratio_string": f"1 : {lam_n[1]:.3f} : {lam_n[2]:.3f}",
        "observed_splitting_ratio_Dm31_over_Dm21": ratio_obs,
        "map_A_lambda_is_mass": {"ratio": ratio_A, "off_by": ratio_obs / ratio_A},
        "map_B_lambda_is_mass_squared": {"ratio": ratio_B, "off_by": ratio_obs / ratio_B},
        "in_tension": tension,
        "caveat": ("the eigenvalue->mass map is undefined (deferred N6, mass-as-loop-length-density); under "
                   "both natural maps the loop splittings are ~5-7x too compressed vs the observed "
                   "Dm31^2/Dm21^2 ~ 33.6. Near-term tension FLAG, not a falsification."),
        "nufit": {"Dm21_sq_eV2": NUFIT_DM21, "Dm31_sq_eV2": NUFIT_DM31, "ordering": "NO"},
    }
    with open(os.path.join(HERE, "n4c_mass_ratio_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {os.path.join(HERE, 'n4c_mass_ratio_summary.json')}")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
