#!/usr/bin/env python3
"""m5_9b_theta13_breaking — where does theta_13 come from? (the SO(3)-breaking)

Issue #199 extension. Duda's follow-up on the PMNS result: "theta13 is clearly
incorrect - could it be a matter of approximations, choice of potential, etc.?"

He is right that theta13 = 0 is wrong: it is the PURE-SO(3) (tri-bimaximal) limit.
The measured theta13 ~ 8.5deg is the SO(3)-BREAKING. This script computes that
breaking at the EFFECTIVE-MODEL level: the SO(3) tri-bimaximal neutrino rotation is
corrected by a charged-lepton-sector rotation (PMNS = U_e^dagger . U_nu), the standard
quark-lepton-complementarity mechanism, and we read off theta13 + the correlations.

The question this answers concretely:
  (1) which correction generates theta13, and in which plane?
  (2) is its size CABIBBO-SCALE (not 0, not maximal)?
  (3) what does it do to the other angles + delta_CP?

What this is NOT: the full M5 field derivation of the correction SIZE (that needs the
charged-lepton mass matrix #200 + V-on on the actual defect = the field follow-up).
Here the correction angle is a parameter; the result is the MECHANISM + the scale.

Run:  python3 m5_9b_theta13_breaking.py
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

# Cabibbo angle (quark sector), the natural correction scale
THETA_C = np.degrees(np.arcsin(0.2250))   # ~13.0 deg  (sin theta_C ~ 0.225)

# NuFIT 6.0 (NO) targets
NUFIT = {"t12": 33.68, "t13": 8.52, "t23": 48.5, "dCP": 177.0}

# ---------------------------------------------------------------------------
# tri-bimaximal neutrino rotation (the SO(3) prediction, theta13 = 0)
# ---------------------------------------------------------------------------
U_TBM = np.array([
    [ np.sqrt(2/3),  np.sqrt(1/3),  0.0          ],
    [-np.sqrt(1/6),  np.sqrt(1/3),  np.sqrt(1/2) ],
    [-np.sqrt(1/6),  np.sqrt(1/3), -np.sqrt(1/2) ],
], dtype=complex)

def rot(plane, theta_deg, phase_deg=0.0):
    """a complex rotation in the given (i,j) plane by theta, with a phase."""
    i, j = plane
    c, s = np.cos(theta_deg * D2R), np.sin(theta_deg * D2R)
    ph = np.exp(1j * phase_deg * D2R)
    U = np.eye(3, dtype=complex)
    U[i, i] = c; U[j, j] = c
    U[i, j] = s * ph.conjugate(); U[j, i] = -s * ph
    return U

def angles_from_U(U):
    """extract (theta12, theta13, theta23) in deg from a unitary U (standard)."""
    t13 = np.degrees(np.arcsin(min(1.0, abs(U[0, 2]))))
    t12 = np.degrees(np.arctan2(abs(U[0, 1]), abs(U[0, 0])))
    t23 = np.degrees(np.arctan2(abs(U[1, 2]), abs(U[2, 2])))
    return t12, t13, t23

def jarlskog(U):
    """rephasing-invariant CP measure J = Im(U_e1 U_mu2 U_e2* U_mu1*)."""
    return np.imag(U[0, 0] * U[1, 1] * np.conj(U[0, 1]) * np.conj(U[1, 0]))

def dCP_from_U(U):
    """delta_CP in [0,360) deg, full quadrant via sin (Jarlskog) + cos (|U_mu1|^2).
    nan if theta13 or theta12 ~ 0 (delta undefined there)."""
    t12, t13, t23 = (a * D2R for a in angles_from_U(U))
    c12, s12, c13, s13, c23, s23 = (np.cos(t12), np.sin(t12), np.cos(t13),
                                    np.sin(t13), np.cos(t23), np.sin(t23))
    denom_s = c12 * s12 * c23 * s23 * c13 ** 2 * s13      # Jarlskog denom -> sin d
    denom_c = 2.0 * s12 * c23 * c12 * s23 * s13           # |U_mu1|^2 relation -> cos d
    if abs(denom_s) < 1e-9 or abs(denom_c) < 1e-9:
        return float("nan")
    sin_d = np.clip(jarlskog(U) / denom_s, -1.0, 1.0)
    cos_d = np.clip((abs(U[1, 0]) ** 2 - s12 ** 2 * c23 ** 2 - c12 ** 2 * s23 ** 2 * s13 ** 2)
                    / denom_c, -1.0, 1.0)
    return np.degrees(np.arctan2(sin_d, cos_d)) % 360.0

# ---------------------------------------------------------------------------
# 1. WHICH PLANE? apply a Cabibbo-size charged-lepton rotation in each plane
#    PMNS = U_e^dagger . U_TBM
# ---------------------------------------------------------------------------
planes = {"1-2 (Cabibbo-like)": (0, 1), "1-3": (0, 2), "2-3": (1, 2)}
plane_results = {}
for name, pl in planes.items():
    U_e = rot(pl, THETA_C, phase_deg=0.0)          # real correction first
    U = U_e.conj().T @ U_TBM
    t12, t13, t23 = angles_from_U(U)
    plane_results[name] = {"theta12": round(t12, 2), "theta13": round(t13, 2),
                           "theta23": round(t23, 2)}

# ---------------------------------------------------------------------------
# 2. THE 1-2 (Cabibbo-like) CORRECTION SWEEP: theta13 vs the correction angle
#    -> theta13 ~ theta_e / sqrt(2)  (the QLC relation)
# ---------------------------------------------------------------------------
theta_e_scan = np.linspace(0.0, 18.0, 60)
t13_scan, t12_scan, t23_scan = [], [], []
for te in theta_e_scan:
    U = rot((0, 1), te, 0.0).conj().T @ U_TBM
    a12, a13, a23 = angles_from_U(U)
    t13_scan.append(a13); t12_scan.append(a12); t23_scan.append(a23)
t13_scan = np.array(t13_scan); t12_scan = np.array(t12_scan); t23_scan = np.array(t23_scan)

# at the Cabibbo point:
U_cab = rot((0, 1), THETA_C, 0.0).conj().T @ U_TBM
t12_c, t13_c, t23_c = angles_from_U(U_cab)
# the analytic QLC relations for reference:
theta13_qlc = THETA_C / np.sqrt(2.0)            # ~9.2 deg
solar_sum = t12_c + THETA_C                       # the "theta12 + theta_C ~ 45" sum rule

# ---------------------------------------------------------------------------
# 2b. THE 2D BEST FIT: a real 1-2 correction over-rotates theta12; a PHASE fixes
#     it (and sets delta_CP). Scan (theta_e, delta_e) -> best match to NuFIT.
# ---------------------------------------------------------------------------
sig = {"t12": 1.0, "t13": 0.5, "t23": 2.0, "dCP": 25.0}   # rough 1-sigma weights
best = {"chi2": 1e9}
te_grid = np.linspace(0.0, 20.0, 81)
de_grid = np.linspace(0.0, 360.0, 145)
for te in te_grid:
    for de in de_grid:
        U = rot((0, 1), te, de).conj().T @ U_TBM
        a12, a13, a23 = angles_from_U(U)
        adCP = dCP_from_U(U)
        if np.isnan(adCP):
            continue
        # delta_CP is circular; fold to nearest to target
        ddCP = (adCP - NUFIT["dCP"] + 180) % 360 - 180
        chi2 = ((a12 - NUFIT["t12"]) / sig["t12"]) ** 2 + ((a13 - NUFIT["t13"]) / sig["t13"]) ** 2 \
            + ((a23 - NUFIT["t23"]) / sig["t23"]) ** 2 + (ddCP / sig["dCP"]) ** 2
        if chi2 < best["chi2"]:
            best = {"chi2": chi2, "theta_e": te, "delta_e": de,
                    "t12": a12, "t13": a13, "t23": a23, "dCP": adCP}

# ---------------------------------------------------------------------------
# 3. delta_CP: with a complex charged-lepton phase, delta_CP moves off 180
# ---------------------------------------------------------------------------
phase_scan = np.linspace(0.0, 360.0, 73)
dCP_scan, t13_phase = [], []
for ph in phase_scan:
    U = rot((0, 1), THETA_C, ph).conj().T @ U_TBM
    dCP_scan.append(dCP_from_U(U))
    t13_phase.append(angles_from_U(U)[1])
dCP_scan = np.array(dCP_scan)
# the phase that best matches NuFIT delta_CP ~ 177:
best_ph = phase_scan[int(np.nanargmin(np.abs(dCP_scan - NUFIT["dCP"])))]

# ---------------------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------------------
summary = {
    "question": "where does theta13 come from? (Duda: approximations / choice of potential?)",
    "answer_short": "the SO(3)-breaking = a charged-lepton (second-rotation) correction to TBM. ROBUST: "
                    "a CABIBBO-SIZE 1-2 correction gives theta13 ~ 9 deg ~ the measured 8.5 deg (the answer "
                    "to 'where does theta13 come from'). TENSION (honest): a SINGLE 1-2 correction cannot fit "
                    "all of NuFIT at once -- a real correction keeps delta_CP=180 but over-rotates theta12; "
                    "fitting theta12 needs a phase that pushes delta_CP toward ~90-100; and theta23 stays "
                    "~45 (the 1-2 rotation cannot reach the upper-octant 48.5). The full fit needs more "
                    "structure (a 2-3 piece, or a non-TBM neutrino sector) -> the field model's job.",
    "cabibbo_angle_deg": round(THETA_C, 2),
    "best_fit_2D": {
        "theta_e_deg": round(best["theta_e"], 2), "delta_e_deg": round(best["delta_e"], 1),
        "chi2": round(best["chi2"], 2),
        "resulting": {"t12": round(best["t12"], 2), "t13": round(best["t13"], 2),
                      "t23": round(best["t23"], 2), "dCP": round(best["dCP"], 1)},
        "nufit": NUFIT,
        "reading": "best single-1-2-correction fit: theta_e ~ Cabibbo-scale, theta12 + theta13 land near "
                   "NuFIT, but delta_CP comes out ~100 (not 180) and theta23 stays ~45 (misses 48.5). "
                   "The (theta12, theta13, delta_CP) trade-off is REAL: you can have any two, not all three, "
                   "from a single 1-2 correction. delta_CP=180 (the SO(3) signature) needs a near-REAL "
                   "correction, which over-rotates theta12. So delta_CP becomes a SHARP discriminator: "
                   "NuFIT's delta_CP ~177 (near 180) currently FAVORS the near-real / SO(3)-preserving "
                   "breaking, NOT a large complex correction.",
    },
    "which_plane": plane_results,
    "plane_verdict": "the 1-2 (Cabibbo-like) charged-lepton rotation reproduces the observed pattern "
                     "(theta13 ~ Cabibbo/sqrt2, theta23 stays ~45, theta12 shifts down by the sum rule); "
                     "a 1-3 rotation overshoots theta13 and misses the solar sum rule; 2-3 mostly moves theta23.",
    "cabibbo_point": {
        "theta_e_deg": round(THETA_C, 2),
        "theta13_deg": round(t13_c, 2), "theta13_qlc_analytic": round(theta13_qlc, 2),
        "theta12_deg": round(t12_c, 2), "theta23_deg": round(t23_c, 2),
        "solar_sum_theta12_plus_thetaC": round(solar_sum, 2),
        "vs_nufit": NUFIT,
    },
    "theta13_is_cabibbo_scale": bool(3.0 < t13_c < 12.0),
    "exact_thetaC_over_sqrt2_excluded": "yes (gives ~9.2 deg; NuFIT ~8.5 deg now excludes the EXACT "
                                        "relation) -> target is CABIBBO-SCALE, not the exact relation",
    "delta_CP": {
        "real_correction_gives_deg": round(dCP_from_U(U_cab), 2) if not np.isnan(dCP_from_U(U_cab)) else 180.0,
        "phase_needed_for_nufit_177": round(best_ph, 1),
        "reading": "a REAL charged-lepton correction keeps delta_CP at 0/180 (still the SO(3) signature); "
                   "a COMPLEX correction phase moves delta_CP continuously -> delta_CP measures the phase "
                   "of the second rotation. NuFIT ~177 needs only a small phase.",
    },
    "v_on_argument": "the calibration thread (#217) found V rotation-invariant on the SINGLE clock "
                     "(dV/dclock=0), so V-on does NOT obviously break the rotation structure -> the "
                     "charged-lepton / second-rotation correction is the more likely theta13 source than "
                     "the potential. The inter-flavour V-on case needs the field sim (the follow-up).",
    "scope": "EFFECTIVE-MODEL: the correction ANGLE is a parameter; the result is the MECHANISM + the "
             "Cabibbo scale. The field derivation of the correction SIZE (charged-lepton mass matrix #200 "
             "+ V-on on the defect) is the named follow-up.",
}
with open(os.path.join(DATA, "m5_9b_theta13_summary.json"), "w") as f:
    json.dump(summary, f, indent=1)

# ---------------------------------------------------------------------------
# PLOT
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

# (a) theta13 vs the 1-2 correction angle -> QLC line
ax[0].plot(theta_e_scan, t13_scan, "b-", lw=2.2, label=r"$\theta_{13}$ (computed)")
ax[0].plot(theta_e_scan, theta_e_scan / np.sqrt(2), "k--", lw=1.2, label=r"$\theta_e/\sqrt{2}$ (QLC)")
ax[0].axvline(THETA_C, color="gray", ls=":", lw=1)
ax[0].axhline(NUFIT["t13"], color="#27ae60", ls=":", lw=1.5, label=f"NuFIT θ₁₃={NUFIT['t13']}°")
ax[0].plot([THETA_C], [t13_c], "ro", ms=8)
ax[0].annotate(f"Cabibbo point\nθ₁₃={t13_c:.1f}°", (THETA_C, t13_c),
               textcoords="offset points", xytext=(-80, 10), fontsize=8)
ax[0].set_xlabel(r"charged-lepton correction $\theta_e$ (deg)")
ax[0].set_ylabel(r"$\theta_{13}$ (deg)")
ax[0].set_title("(a) θ₁₃ from the charged-lepton correction\n(θ₁₃ ≈ θ_e/√2; Cabibbo → ~9° ≈ measured)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

# (b) which plane: the 3 corrections vs NuFIT pattern
labels = [r"$\theta_{12}$", r"$\theta_{13}$", r"$\theta_{23}$"]
x = np.arange(3)
colors = {"1-2 (Cabibbo-like)": "#c0392b", "1-3": "#8e44ad", "2-3": "#16a085"}
w = 0.22
for k, (name, r) in enumerate(plane_results.items()):
    vals = [r["theta12"], r["theta13"], r["theta23"]]
    ax[1].bar(x + (k - 1) * w, vals, w, label=name, color=colors[name])
ax[1].plot(x, [NUFIT["t12"], NUFIT["t13"], NUFIT["t23"]], "k*", ms=14, label="NuFIT 6.0")
ax[1].set_xticks(x); ax[1].set_xticklabels(labels)
ax[1].set_ylabel("angle (deg)")
ax[1].set_title("(b) which correction plane?\n1-2 (Cabibbo) matches the pattern")
ax[1].legend(fontsize=7.5); ax[1].grid(alpha=0.3, axis="y")

# (c) delta_CP vs the correction phase
ax[2].plot(phase_scan, dCP_scan, "b-", lw=2)
ax[2].axhline(180, color="#2980b9", ls=":", lw=1.5, label="SO(3) real: δ_CP=180°")
ax[2].axhline(NUFIT["dCP"], color="#e67e22", ls="--", lw=1.5, label=f"NuFIT ~{NUFIT['dCP']:.0f}°")
ax[2].set_xlabel(r"charged-lepton correction phase $\delta_e$ (deg)")
ax[2].set_ylabel(r"$\delta_{CP}$ (deg)")
ax[2].set_title("(c) δ_CP from the correction phase\n(real → 0/180; complex → moves it)")
ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(PLOTS, "m5_9b_theta13_breaking.png"), dpi=130)
plt.close()

# ---------------------------------------------------------------------------
# CONSOLE
# ---------------------------------------------------------------------------
print("=" * 72)
print("m5_9b — WHERE DOES theta13 COME FROM? (the SO(3)-breaking, #199 ext)")
print("=" * 72)
print(f"Cabibbo angle theta_C = {THETA_C:.2f} deg")
print("-" * 72)
print("WHICH PLANE (Cabibbo-size correction, PMNS = U_e^dag . U_TBM):")
for name, r in plane_results.items():
    print(f"  {name:20s}: theta12={r['theta12']:5.1f}  theta13={r['theta13']:5.1f}  theta23={r['theta23']:5.1f}")
print(f"  NuFIT 6.0           : theta12={NUFIT['t12']:5.1f}  theta13={NUFIT['t13']:5.1f}  theta23={NUFIT['t23']:5.1f}")
print("-" * 72)
print(f"1-2 CABIBBO POINT: theta13 = {t13_c:.2f} deg  (analytic theta_C/sqrt2 = {theta13_qlc:.2f})")
print(f"  -> CABIBBO-SCALE, ~ the measured {NUFIT['t13']} deg (exact theta_C/sqrt2 EXCLUDED, scale matches)")
print(f"  solar sum rule: theta12 + theta_C = {solar_sum:.1f} deg (~45 = bimaximal)")
print(f"  theta23 stays {t23_c:.1f} (~maximal); theta12 shifts to {t12_c:.1f} (NuFIT {NUFIT['t12']})")
print("-" * 72)
print(f"delta_CP: real correction -> {dCP_from_U(U_cab):.0f}; phase for NuFIT 177 -> delta_e ~ {best_ph:.0f} deg")
print("-" * 72)
print("2D BEST FIT (1-2 correction, angle + phase) vs NuFIT:")
print(f"  theta_e = {best['theta_e']:.1f} deg (Cabibbo-scale), delta_e = {best['delta_e']:.0f} deg, chi2 = {best['chi2']:.2f}")
print(f"  -> theta12={best['t12']:.1f} (NuFIT {NUFIT['t12']}), theta13={best['t13']:.1f} ({NUFIT['t13']}), "
      f"theta23={best['t23']:.1f} ({NUFIT['t23']}), dCP={best['dCP']:.0f} ({NUFIT['dCP']})")
print("  => theta12+theta13 fit at Cabibbo scale, but dCP~100 (not 180) + theta23~45 (misses 48.5):")
print("     a SINGLE 1-2 correction gives any TWO of (theta12, theta13, dCP=180), not all three.")
print("     dCP=180 needs a NEAR-REAL correction -> NuFIT dCP~177 FAVORS the SO(3)-preserving breaking.")
print("-" * 72)
print("V-on: #217 found V rotation-invariant on the clock -> charged-lepton correction is the likelier source")
print("=" * 72)
print("Wrote data/m5_9b_theta13_summary.json + plots/m5_9b_theta13_breaking.png")
