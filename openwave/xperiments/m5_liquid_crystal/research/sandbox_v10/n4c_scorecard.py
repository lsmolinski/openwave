#!/usr/bin/env python3
"""
N4c-3 , the HONEST scorecard figure (replaces the over-claimed "3 predictions" headline).

After the decisive alpha-energy test (n4c_alpha_energy.py) and the mass-ratio falsifier (n4c_mass_ratio.py),
the honest provenance of the four PMNS parameters is:
  - theta_23 = 45, theta_13 = 0 (baseline), |delta_CP| = 90  : CONSEQUENCES of ONE imposed mu-tau symmetry
                                                                (Harrison-Scott), not independent predictions.
  - theta_12 = 35.26                                          : geometrically PINNED by the magic condition
                                                                (a derived locus), but NOT energetically
                                                                selected (E_self flat to 0.09%). Between fit
                                                                and prediction.
  - theta_13 = 8.56                                           : one-parameter FIT (g_chiral ~ 0.94).
  - delta_CP SIGN                                             : UNDETERMINED (handedness-degenerate).

Panel 1: the provenance, colour-coded. Panel 2: a pull plot (deviation in sigma) vs NuFIT 6.0 NO with error
bars -> theta_12 is ~2.3 sigma off (NOT "close"), theta_23 marginal, theta_13 set (pull 0 by construction),
delta_CP consistent (the range is enormous).

Headless f64. LOCAL (#236 N-program, HELD). Run: python3 n4c_scorecard.py
"""

import json
import os
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# this work
PRED = {"theta12": 35.26, "theta23": 45.00, "theta13": 8.56, "delta_CP": 270.0}
# NuFIT 6.0 NO (approx central +- 1sigma; theta23 + delta_CP carry octant/large ranges)
NUFIT = {
    "theta12": (33.68, 0.70),
    "theta23": (43.3, 1.0),
    "theta13": (8.56, 0.11),
    "delta_CP": (212.0, 30.0),
}
# honest provenance category per parameter
PROV = [
    ("theta_23 = 45", "imposed mu-tau symmetry", "#7f8c8d"),
    ("theta_13 = 0 (baseline)", "imposed mu-tau symmetry", "#7f8c8d"),
    ("|delta_CP| = 90", "imposed mu-tau reflection", "#95a5a6"),
    ("theta_12 = 35.26", "geometrically PINNED (magic); not energy-selected", "#2980b9"),
    ("theta_13 = 8.56", "FREE coupling (g_chiral ~ 0.94)", "#e67e22"),
    ("delta_CP SIGN", "UNDETERMINED (handedness-degenerate)", "#c0392b"),
]


def main():
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.0))

    # ---- panel 1: honest provenance ----
    y = np.arange(len(PROV))[::-1]
    for yi, (label, prov, col) in zip(y, PROV):
        ax[0].barh(yi, 1.0, color=col, alpha=0.85)
        ax[0].text(0.02, yi, f"  {label}", va="center", ha="left", fontsize=10, color="white", weight="bold")
        ax[0].text(0.98, yi, f"{prov}  ", va="center", ha="right", fontsize=8.5, color="white")
    ax[0].set_xlim(0, 1); ax[0].set_ylim(-0.6, len(PROV) - 0.4)
    ax[0].set_xticks([]); ax[0].set_yticks([])
    ax[0].set_title("Honest provenance: 1 imposed symmetry (3 consequences)\n+ 1 pinned angle + 1 free coupling "
                    "+ 1 undetermined sign", fontsize=10)

    # ---- panel 2: pull plot (deviation in sigma) ----
    names = ["theta12", "theta23", "theta13", "delta_CP"]
    disp = ["theta_12", "theta_23", "theta_13", "delta_CP"]
    pulls, notes = [], []
    for k in names:
        c, s = NUFIT[k]
        pull = (PRED[k] - c) / s
        pulls.append(pull)
    colors = ["#2980b9", "#7f8c8d", "#e67e22", "#95a5a6"]
    yb = np.arange(len(names))[::-1]
    ax[1].barh(yb, pulls, color=colors, alpha=0.85)
    for thr in (1, 2, 3):
        ax[1].axvline(thr, ls=":", color="#bbb"); ax[1].axvline(-thr, ls=":", color="#bbb")
    ax[1].axvline(0, color="#2c3e50", lw=1)
    note = {"theta12": "+2.3 sigma (NOT close)", "theta23": "~1.7 sigma (octant-ambiguous)",
            "theta13": "SET (fit -> pull 0)", "delta_CP": "consistent (range enormous)"}
    for yi, k, p in zip(yb, names, pulls):
        ax[1].text(p + (0.1 if p >= 0 else -0.1), yi, note[k], va="center",
                   ha="left" if p >= 0 else "right", fontsize=8.5)
    ax[1].set_yticks(yb); ax[1].set_yticklabels(disp, fontsize=10)
    ax[1].set_xlabel("(this work - NuFIT 6.0 NO) / 1 sigma   [pull]")
    ax[1].set_xlim(-4, 5)
    ax[1].set_title("Pull vs data with error bars\n(theta_13 set; theta_12 ~2.3 sigma; delta_CP consistent)",
                    fontsize=10)
    ax[1].grid(True, axis="x", alpha=0.2)

    fig.tight_layout()
    fig_path = os.path.join(HERE, "n4c_scorecard.png")
    fig.savefig(fig_path, dpi=120)
    print(f"figure -> {fig_path}")

    summary = {
        "predicted": PRED, "nufit_central_1sigma": NUFIT,
        "pulls_sigma": {k: (PRED[k] - NUFIT[k][0]) / NUFIT[k][1] for k in names},
        "honest_provenance": [{"param": p[0], "provenance": p[1]} for p in PROV],
        "note": ("theta_23/theta_13=0/|delta_CP|=90 are consequences of ONE imposed mu-tau symmetry; "
                 "theta_12 geometrically pinned (magic) but not energy-selected; theta_13 a free coupling; "
                 "delta_CP sign undetermined. theta_12 is ~2.3 sigma from data; delta_CP is consistent."),
    }
    with open(os.path.join(HERE, "n4c_scorecard.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(f"summary -> {os.path.join(HERE, 'n4c_scorecard.json')}")
    return True


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
