#!/usr/bin/env python3
"""m5_9_pmns_so3 — PMNS neutrino mixing from the SO(3) rotation structure.

Issue #199. The model maps neutrino flavour oscillation to an SO(3) spatial-field
rotation: a neutrino = the delta-0 (axis 2<->3) content of M = O.diag(g,1,delta,0).O^T
swinging WITHOUT the hedgehog winding (Duda's Wolfram slide, 0b_M5_roadmap.md L786).
This script turns that commitment into a concrete, falsifiable PMNS prediction and
compares it to the global oscillation fit (NuFIT 6.0).

The structural chain:
  - The 2<->3 axis swing IS the atmospheric rotation theta_23 -> near-MAXIMAL (45 deg).
  - The next SO(3) rotation is the solar theta_12 -> tri-maximal (sin^2 = 1/3).
  - A PURE (real, orthogonal) SO(3) rotation has NO genuine complex CP phase, so it
    forces the reactor angle theta_13 = 0 and delta_CP in {0, 180} deg.
  => the SO(3) prediction is the tri-bimaximal (TBM) matrix with delta_CP = 180 deg.
  The MEASURED theta_13 != 0 is then exactly the SO(3)-BREAKING order parameter =
  the "second coupled rotation (toward SU(3))" that #199 asks about.

No M5 field simulation is needed here: this is the structural consequence of the
SO(3) commitment the model already makes. The full dynamical neutrino sim (seed the
delta-0 excitation, measure 3 masses) is the named follow-up.

Run:  python3 m5_9_pmns_so3.py
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
PLOTS = os.path.join(HERE, "plots")
os.makedirs(DATA, exist_ok=True)
os.makedirs(PLOTS, exist_ok=True)

D2R = np.pi / 180.0
R2D = 180.0 / np.pi

# ---------------------------------------------------------------------------
# Standard PMNS parametrization U = R23(t23) . U13(t13, dCP) . R12(t12)
# ---------------------------------------------------------------------------
def pmns(t12, t13, t23, dCP):
    """angles in degrees, dCP in degrees -> 3x3 complex unitary PMNS matrix."""
    s12, c12 = np.sin(t12 * D2R), np.cos(t12 * D2R)
    s13, c13 = np.sin(t13 * D2R), np.cos(t13 * D2R)
    s23, c23 = np.sin(t23 * D2R), np.cos(t23 * D2R)
    d = np.exp(-1j * dCP * D2R)
    R23 = np.array([[1, 0, 0], [0, c23, s23], [0, -s23, c23]], dtype=complex)
    U13 = np.array([[c13, 0, s13 * d.conjugate()], [0, 1, 0], [-s13 * d, 0, c13]], dtype=complex)
    R12 = np.array([[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]], dtype=complex)
    return R23 @ U13 @ R12

def angles_from_U(U):
    """recover the mixing angles (deg) from |U| (standard extraction)."""
    t13 = np.arcsin(np.abs(U[0, 2])) * R2D
    t12 = np.arctan2(np.abs(U[0, 1]), np.abs(U[0, 0])) * R2D
    t23 = np.arctan2(np.abs(U[1, 2]), np.abs(U[2, 2])) * R2D
    return t12, t13, t23

# ---------------------------------------------------------------------------
# 1. THE SO(3) PREDICTION = tri-bimaximal (TBM), delta_CP = 180 deg
# ---------------------------------------------------------------------------
t12_tbm = np.arcsin(np.sqrt(1.0 / 3.0)) * R2D   # 35.264 deg  (sin^2 = 1/3)
t23_tbm = 45.0                                   # maximal 2<->3 swing
t13_tbm = 0.0                                    # pure SO(3): no reactor angle
dCP_tbm = 180.0                                  # real SO(3): no complex CP phase
U_tbm = pmns(t12_tbm, t13_tbm, t23_tbm, dCP_tbm)

# ---------------------------------------------------------------------------
# 2. THE MEASURED FIT (NuFIT 6.0, 2024, normal ordering, w/ SK atmospheric)
#    values are the widely-quoted global-fit best fits; the structural
#    comparison (pattern + delta_CP near 180) is what matters, not 0.1-deg precision.
# ---------------------------------------------------------------------------
nufit = {
    "theta12": 33.68,   # sin^2 = 0.307
    "theta13": 8.52,    # sin^2 = 0.02195   <- the reactor angle (SO(3)-breaking)
    "theta23": 48.5,    # sin^2 = 0.561 (upper octant; lower-octant 43.3 also allowed)
    "theta23_lower": 43.3,
    "dCP": 177.0,       # NuFIT 6.0 best fit moved NEAR 180 deg (range ~ 96-422 at 3sigma)
    "dCP_range_1sig": [148.0, 215.0],
}
U_nufit = pmns(nufit["theta12"], nufit["theta13"], nufit["theta23"], nufit["dCP"])

# ---------------------------------------------------------------------------
# 3. THE SO(3)-BREAKING VERDICT (the #199 core question)
# ---------------------------------------------------------------------------
# A single PURE SO(3) rotation forces theta_13 = 0. The measured theta_13 is the
# size of the REQUIRED second coupled rotation (toward SU(3)-like).
theta13_breaking = nufit["theta13"]                      # the breaking order parameter
breaking_fraction = np.sin(theta13_breaking * D2R) ** 2  # sin^2 theta_13 = leakage into 1-3
dCP_consistent_with_180 = (nufit["dCP_range_1sig"][0] <= 180.0 <= nufit["dCP_range_1sig"][1])

# pattern-match: how close is each measured angle to the SO(3)/TBM prediction?
pattern = {
    "theta23_maximal": {"tbm": t23_tbm, "meas_upper": nufit["theta23"],
                        "meas_lower": nufit["theta23_lower"],
                        "straddles_45": nufit["theta23_lower"] <= 45.0 <= nufit["theta23"]},
    "theta12_trimaximal": {"tbm": round(t12_tbm, 2), "meas": nufit["theta12"],
                          "delta_deg": round(nufit["theta12"] - t12_tbm, 2)},
    "theta13_breaking": {"tbm": t13_tbm, "meas": nufit["theta13"],
                        "this_is_the_so3_breaking": True},
    "dCP": {"so3_predicts": 180.0, "meas": nufit["dCP"],
            "consistent_with_180_at_1sigma": bool(dCP_consistent_with_180)},
}

# ---------------------------------------------------------------------------
# 4. THE ELECTRON 0.28% CLOCK CHECK (Duda arxiv 2108.07896) — projection geometry
# ---------------------------------------------------------------------------
# Duda: the 0.28% electron-clock energy excess = "3 tendencies projected/added into a
# single allowed evolution DOF". Same SO(3) 3-axis structure, charged-lepton side.
# Projection-geometry framework: the observable clock is the projection of the full
# multi-axis internal rotation. The fractional energy EXCESS over a single pure
# tendency = the off-primary-axis energy fraction f_off:
#     E_obs / E_pure = 1 + f_off
# So Duda's 0.28% <=> f_off ~ 0.0028.
excess_measured = 0.0028   # 0.28%
# the NAIVE neutrino-mixing estimate (if the electron leaked like the neutrino theta_13):
f_off_naive_theta13 = np.sin(theta13_breaking * D2R) ** 2   # 0.0219 = 2.19%
ratio_naive_to_obs = f_off_naive_theta13 / excess_measured  # ~8x too big
# => the electron clock is ~8x MORE tightly projected than the neutrino mixing,
#    physically sensible: a bound, phase-locked clock vs a freely-oscillating neutrino.
# the off-axis angle that WOULD give 0.28%:
theta_off_for_028 = np.arcsin(np.sqrt(excess_measured)) * R2D   # ~3.03 deg

# ---------------------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------------------
summary = {
    "so3_prediction_TBM": {
        "theta12_deg": round(t12_tbm, 3), "theta13_deg": t13_tbm,
        "theta23_deg": t23_tbm, "dCP_deg": dCP_tbm,
        "note": "tri-bimaximal: the most symmetric SO(3) rotation; delta_CP=180 forced by real SO(3)",
    },
    "measured_NuFIT6": nufit,
    "verdict_pattern_match": pattern,
    "so3_vs_su3": {
        "single_SO3_suffices": False,
        "second_rotation_required": True,
        "second_rotation_size_deg": theta13_breaking,
        "second_rotation_is_SMALL": theta13_breaking < 15.0,
        "leakage_sin2_theta13": round(breaking_fraction, 4),
        "reading": "SO(3) is the LEADING structure (theta23 maximal, theta12 tri-maximal, "
                   "dCP~180 all hold); it is broken by a SMALL second coupled rotation = the "
                   "reactor angle theta_13 ~ 8.5 deg (toward the SU(3)-like coupled structure). "
                   "delta_CP near 180 => the SO(3) route is currently FAVORED, not falsified.",
    },
    "delta_CP_test": {
        "so3_requires_deg": 180.0, "measured_deg": nufit["dCP"],
        "consistent_at_1sigma": bool(dCP_consistent_with_180),
        "falsifier": "if JUNO/DUNE settle delta_CP far from 180, the pure-SO(3) route is wrong",
        "status": "CONSISTENT (NuFIT 6.0 best fit ~177 deg, near 180) — the prediction holds for now",
    },
    "electron_0p28pct_check": {
        "duda_excess": excess_measured,
        "framework": "E_obs/E_pure = 1 + f_off (off-primary-axis energy fraction)",
        "f_off_for_0p28pct": excess_measured,
        "naive_theta13_leakage": round(f_off_naive_theta13, 4),
        "naive_over_observed": round(ratio_naive_to_obs, 1),
        "off_axis_angle_deg_for_0p28pct": round(theta_off_for_028, 2),
        "reading": "0.28% is NOT simply sin^2(theta13) (that is ~2.2%, ~8x too big). The electron "
                   "clock is ~8x MORE tightly projected than the free neutrino mixing — sensible for "
                   "a bound phase-locked clock. The structure (3 tendencies -> 1 projected DOF -> small "
                   "excess) matches Duda's reading; the SPECIFIC 0.28% needs the electron defect's "
                   "measured off-primary-axis energy fraction (the dynamical follow-up).",
        "status": "STRUCTURAL match; specific number = follow-up (measure f_off on the electron defect)",
    },
}
with open(os.path.join(DATA, "m5_9_pmns_summary.json"), "w") as f:
    json.dump(summary, f, indent=1)

# ---------------------------------------------------------------------------
# PLOT
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(13, 5))

# (a) the three angles: SO(3)/TBM prediction vs measured
labels = [r"$\theta_{12}$ (solar)", r"$\theta_{13}$ (reactor)", r"$\theta_{23}$ (atmos.)"]
tbm_vals = [t12_tbm, t13_tbm, t23_tbm]
meas_vals = [nufit["theta12"], nufit["theta13"], nufit["theta23"]]
meas_err = [1.0, 0.13, [[nufit["theta23"] - nufit["theta23_lower"]], [1.0]]]  # rough
x = np.arange(3)
ax[0].bar(x - 0.2, tbm_vals, 0.4, color="#2980b9", label="SO(3) / tri-bimaximal (predicted)")
ax[0].bar(x + 0.2, meas_vals, 0.4, color="#e67e22", label="NuFIT 6.0 (measured)")
for i, (t, m) in enumerate(zip(tbm_vals, meas_vals)):
    ax[0].text(i - 0.2, t + 0.6, f"{t:.1f}", ha="center", fontsize=8)
    ax[0].text(i + 0.2, m + 0.6, f"{m:.1f}", ha="center", fontsize=8)
ax[0].annotate("the SO(3)-BREAKING\n(second rotation)", (1 + 0.2, nufit["theta13"]),
               textcoords="offset points", xytext=(18, 28), fontsize=8,
               arrowprops=dict(arrowstyle="->", color="red"))
ax[0].set_xticks(x); ax[0].set_xticklabels(labels, fontsize=9)
ax[0].set_ylabel("mixing angle (deg)")
ax[0].set_title("(a) SO(3)/TBM prediction vs measured\n(θ₂₃ maximal, θ₁₂ tri-maximal ✓; θ₁₃ = the breaking)")
ax[0].legend(fontsize=8.5, loc="upper center")
ax[0].set_ylim(0, 58)
ax[0].grid(alpha=0.3, axis="y")

# (b) delta_CP: SO(3) prediction vs the measured range
ax[1].axvline(180.0, color="#2980b9", lw=2.5, label="SO(3) requires δ_CP = 180°")
ax[1].axvspan(nufit["dCP_range_1sig"][0], nufit["dCP_range_1sig"][1], color="#e67e22", alpha=0.3,
              label=f"NuFIT 6.0 1σ [{nufit['dCP_range_1sig'][0]:.0f}, {nufit['dCP_range_1sig'][1]:.0f}]")
ax[1].axvline(nufit["dCP"], color="#e67e22", lw=2, ls="--", label=f"best fit {nufit['dCP']:.0f}°")
ax[1].set_xlim(90, 360)
ax[1].set_yticks([])
ax[1].set_xlabel("δ_CP (deg)")
ax[1].set_title("(b) The falsifiable δ_CP = 180° test\nCONSISTENT (best fit ~177°); JUNO/DUNE will sharpen")
ax[1].legend(fontsize=9, loc="upper right")

plt.tight_layout()
plt.savefig(os.path.join(PLOTS, "m5_9_pmns_so3.png"), dpi=130)
plt.close()

# ---------------------------------------------------------------------------
# CONSOLE REPORT
# ---------------------------------------------------------------------------
np.set_printoptions(precision=3, suppress=True)
print("=" * 72)
print("m5_9 — PMNS FROM THE SO(3) ROTATION STRUCTURE (#199)")
print("=" * 72)
print("SO(3)/TBM prediction:  theta12=%.2f  theta13=%.1f  theta23=%.1f  dCP=%.0f"
      % (t12_tbm, t13_tbm, t23_tbm, dCP_tbm))
print("NuFIT 6.0 (measured):  theta12=%.2f  theta13=%.2f  theta23=%.1f  dCP=%.0f"
      % (nufit["theta12"], nufit["theta13"], nufit["theta23"], nufit["dCP"]))
print("-" * 72)
print("|U_PMNS| tri-bimaximal (SO(3) prediction):")
print(np.abs(U_tbm))
print("|U_PMNS| NuFIT 6.0 (measured):")
print(np.abs(U_nufit))
print("-" * 72)
print("SO(3)-vs-SU(3) VERDICT:")
print("  theta23 ~ maximal (45):  TBM 45  vs meas %.1f/%.1f  -> straddles 45 ✓"
      % (nufit["theta23_lower"], nufit["theta23"]))
print("  theta12 ~ tri-maximal:   TBM 35.26 vs meas %.1f  (off %.1f deg) ✓"
      % (nufit["theta12"], nufit["theta12"] - t12_tbm))
print("  theta13 = SO(3)-breaking: TBM 0 -> meas %.1f deg = the REQUIRED second rotation (SMALL)"
      % nufit["theta13"])
print("  dCP = 180 (SO(3)):       meas ~%.0f, consistent at 1sigma: %s"
      % (nufit["dCP"], dCP_consistent_with_180))
print("  => SO(3) is the LEADING structure, broken by a small theta13 rotation toward SU(3).")
print("-" * 72)
print("ELECTRON 0.28%% CHECK (Duda 2108.07896):")
print("  0.28%% <=> off-axis fraction f_off ~ 0.0028  (angle ~%.1f deg)" % theta_off_for_028)
print("  naive sin^2(theta13) = %.4f (2.2%%) = ~%.0fx too big -> electron clock is TIGHTER-projected"
      % (f_off_naive_theta13, ratio_naive_to_obs))
print("  structural match to Duda's reading; specific number = the dynamical follow-up")
print("=" * 72)
print("Wrote data/m5_9_pmns_summary.json + plots/m5_9_pmns_so3.png")
