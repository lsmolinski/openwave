"""
M5.8.2j — N-3: THE FIRST ZBW RATIO — `ω·ℏ / (2·H_rest)` on the saturated
breather, in lattice-natural units.

THE TARGET: Duda's Zitterbewegung clock runs at ω = 2mc²/ℏ — the ratio
ω·ℏ/(2H_rest) = 1 for a calibrated ground-state particle clock.

PRE-REGISTERED HONEST FRAMING (roadmap N-3):
  · this is *the first measurement of the ratio*, NOT a pass/fail vs 1 —
    our state is an EXCITED breather, not a ground clock, and
  · the sandbox units are UNCALIBRATED (the classical action never fixed
    ℏ; the lattice (h, dt) are arbitrary): the number below is the ratio
    in lattice-natural units (ℏ = 1 by interpretation). The ABSOLUTE
    ratio needs the Faber `r₀` mass calibration — the N-6 ZBW program (promoted 2026-06-07).
    What IS unit-free here: the ratio's spread across arms (start-
    independence — inherits N-1's attractor) and across H_rest choices.

INPUTS (all SAVED data — CPU-only, no Taichi, runs alongside N-2's GPU):
  ω₁  — the N-1 attractor fundamental per arm (the m5_8_2h detrend-stable
        p2 peak; N-1 step 1 established the saturated state does NOT carry
        the linear 5.86 fast clock — ω₁ is the only resolved carrier)
  H_rest, two candidates (transparency — both reported):
    A. STATIC: H of the dressed seed = Σ_act (u + βu²)·h³ with T = 0 exact
       (the defect+dressing rest energy — closest to "mc²")
    B. DYNAMIC: median H over the breather's dt-converged early window
       (the state the clock actually rides — includes excitation energy)

RESULTS (2026-06-07 — N-3 COMPLETE, the first measurement on record):
  H_rest A (static dressed seed, T = 0): 16.74 quartic (9.21 quadratic part).
  | arm | ω₁ (median over detrend dial) | ω₁/(2H_static) | ω₁/(2H_dyn) |
  | R0-kicked | 1.092 | 0.03261 | 0.01147 |
  | R2-seed   | 1.148 | 0.03430 | 0.01263 |
  | R3-jit    | 1.089 | 0.03255 | 0.01186 |
  · cross-arm spread of ω₁/(2H_static) = 5.4% — the UNIT-FREE statement:
    the ratio is start-independent (inherits the N-1 attractor).
  · the breather is strongly EXCITED: H_dyn/H_static ≈ 2.7–2.8 (the state
    carries ~1.8× the rest energy as excitation — quantifying the "excited
    breather, not ground clock" caveat).
  · extractor note: the MEDIAN over the detrend dial is the robust ω₁
    readout — at single dial settings the 2ω₁ harmonic can outrank the
    fundamental (R0 max 2.10, R3 max 2.64 = rank flips, not drift; R2 is
    razor-tight 1.148–1.154).
  · the M5.8.3-spec ℏ↔δ convention (δ = 0.3): δ·ω₁/(2H_static) ≈ 0.0098–0.0103.
  · G7 factor-2 bookkeeping respected: the M-field probes (p2/dM) read
    ω_M = 2ω_clock = the ZBW frequency directly — no double-correction.
  · WHY "VALIDATED" IS NOT CLAIMABLE WITHOUT CALIBRATION (the sharp form):
    rescaling the classical action L → λL rescales H → λH but leaves the
    EOM (hence ω) unchanged ⇒ ω·ℏ/(2H) = 1 is the DEFINITION of the
    calibrated ℏ (equivalently δ), not a test, until ℏ/δ is fixed
    INDEPENDENTLY (the Faber r₀ → 0.511 MeV chain — N-6b). The unit-free
    ZBW tests available pre-calibration: (a) start-independence of the
    ratio (measured here, 5.4%); (b) the ZBW LAW ω ∝ H_rest across a mass
    family (vary r₀/core scale, ratio constant) — N-6a (was NG-1/NG-4).
  · the carry-forward number for N-6b (the absolute calibration): in lattice-
    natural units ω₁/(2H_static) ≈ 0.033 on the 24³ β = 1.558 stack — the
    Faber r₀ calibration must supply a units factor of ~30× to reach the
    ZBW ω = 2mc²/ℏ; whether that is physical or exposes the model is
    exactly N-6b's question. NO pass/fail claimed here.

N-6b ADDENDUM (2026-06-07 — the first ABSOLUTE ω, analysis):
  Under two EXPLICIT postulates — (P1) H_static(R_W=3.5)=16.74 ↔ 0.511 MeV,
  (P2) lattice action unit ↔ ℏ — the measured ω₁ = 1.188 maps to
  **ω = 5.51 × 10¹⁹ rad/s vs the electron ZBW 1.553 × 10²¹ — a factor
  ≈28 short**. With N-6a's rigidity (ω insensitive to dressing energy) the
  gap is STRUCTURAL, not closable by energy bookkeeping: candidate physics
  = the V-on/Faber-r₀ core sector, the faithful kinetic (frequency-level,
  5d), or P2 itself (no independent ℏ exists in a classical sim; the
  length-anchor alternative needs the V-on core — the V=0 frame is
  scale-free). First measurement recorded, anchors explicit, no pass/fail.

USAGE:  python m5_8_2j_zbw_ratio.py
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
V8 = HERE.parent / "sandbox_v8"
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

# ti-free chain only (N-2 owns the GPU right now)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v6.m5_6_2a_biaxial_hedgehog import (  # noqa: E402
    DELTA,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_v8.m5_8_2c1_full_evolution import (  # noqa: E402
    central, tw,
)
from openwave.xperiments.m5_liquid_crystal.research.sandbox_vn.m5_8_2h_omega_attractor import (  # noqa: E402
    np_commf, spec, peaks, healthy_end,
)

BETA = 1.558


def u_density(M, h):
    u = 0.0
    Mi = [central(M, ax, h) for ax in range(3)]
    for i in range(3):
        for j in range(i + 1, 3):
            F = np_commf(Mi[i], Mi[j])
            u = u + 2.0 * np.einsum("...ab,...ab->...", F, tw(F))
    return u


def omega1(z, tag):
    """The N-1 pipeline verbatim: healthy window, detrend dial, top p2 peak."""
    dsamp = float(z["dense"]) * float(z["dt"])
    ne = healthy_end(z[f"{tag}_umin"])
    y = z[f"{tag}_p2"][:ne]
    tops = []
    for T in (2.5, 5.0, 10.0):
        n = max(int(T / dsamp) | 1, 3)
        pf = peaks(*spec(y, dsamp, detrend_n=n), 1.0, np.inf, n=1)
        if pf:
            tops.append(pf[0][0])
    return (float(np.median(tops)), float(np.min(tops)), float(np.max(tops)),
            ne)


def main():
    d = np.load(V8 / "_m5_8_2cb_ref.npz")
    M0, act, h = d["M_seed"], d["act"] > 0.5, float(d["h"])
    u0 = u_density(M0, h)
    H_static = float((u0 + BETA * u0 * u0)[act].sum()) * h ** 3
    H_quad = float(u0[act].sum()) * h ** 3
    print("=" * 78)
    print("M5.8.2j — the first ZBW ratio  ω₁/(2·H_rest)  [lattice-natural"
          " units, ℏ ≡ 1]")
    print("=" * 78)
    print(f"  H_rest A (static dressed seed, T = 0 exact):  quartic"
          f" {H_static:.4f}   (quadratic part {H_quad:.4f})")
    z = dict(np.load(HERE / "_m5_8_2h_dense.npz"))
    print("\n  | arm | ω₁ median [min–max over detrend dial] | H_dyn"
          " (median, early dt-converged window) | ω₁/(2H_static) |"
          " ω₁/(2H_dyn) |")
    print("  | --- | --- | --- | --- | --- |")
    ratios = []
    for tag in ("R0-kicked", "R2-seed", "R3-jit"):
        om, om_lo, om_hi, ne = omega1(z, tag)
        # the dynamic H over the matching healthy window (H probe is sparser)
        tH = z[f"{tag}_ht"]
        t_end = z[f"{tag}_t"][ne - 1]
        Hwin = z[f"{tag}_hH"][tH <= t_end]
        # early third of the healthy window = the dt-converged regime
        H_dyn = float(np.median(Hwin[:max(len(Hwin) // 3, 1)]))
        rs = om / (2 * H_static)
        rd = om / (2 * abs(H_dyn))
        ratios.append(rs)
        print(f"  | {tag} | {om:.3f} [{om_lo:.3f}–{om_hi:.3f}] |"
              f" {H_dyn:+.3f} | {rs:.5f} | {rd:.5f} |")
    print(f"\n  cross-arm spread of ω₁/(2H_static): "
          f"{min(ratios):.5f} – {max(ratios):.5f}"
          f"  ({100 * (max(ratios) / min(ratios) - 1):.1f}% — the unit-free"
          " start-independence statement)")
    print(f"  M5.8.3-spec convention ℏ↔δ (δ = {DELTA}):  δ·ω₁/(2H_static) ="
          f" {DELTA * min(ratios):.5f} – {DELTA * max(ratios):.5f}")
    print("  G7 factor-2 bookkeeping (apolar doubling, ω_M = 2ω_clock"
          " machine-exact): the p2/dM")
    print("  probes read the M-FIELD oscillation, which IS the ZBW frequency"
          " 2mc²/ℏ directly —")
    print("  per the M5.8.3 rule, NO double-correction is applied (and none"
          " would be valid).")
    print("\n  READ HONESTLY: the ratio's ABSOLUTE value is meaningless until"
          " the Faber r₀")
    print("  calibration (N-6b) fixes the lattice units — ω = 2mc²/ℏ is a"
          " statement about")
    print("  CALIBRATED units. What this measurement establishes: (i) the"
          " ratio exists and is")
    print("  start-independent (the N-1 attractor); (ii) the saturated"
          " breather's clock is the")
    print("  breathing fundamental, NOT the linear 5.86 mode (N-1 step 1);"
          " (iii) the number to")
    print("  carry into N-6b: ω₁/(2H_static) above, on the 24³ β = 1.558"
          " stack.")
    # ── N-6b: THE FIRST ABSOLUTE ω (explicit-postulate calibration) ────────
    MEV_OVER_HBAR = 0.511e6 * 1.602176634e-19 / 1.054571817e-34  # rad/s per 0.511 MeV
    om_lat = 1.188                       # the validated member (2m RW3.5)
    om_phys = om_lat * MEV_OVER_HBAR / H_static
    zbw = 2.0 * MEV_OVER_HBAR
    print("\n[N-6b] the FIRST ABSOLUTE ω — under TWO EXPLICIT POSTULATES:")
    print("  (P1) energy anchor: H_static(R_W=3.5) = 16.74 lattice ↔ m_e c²"
          " = 0.511 MeV")
    print("  (P2) action anchor: the lattice action unit ↔ ℏ (no independent"
          " ℏ exists in a classical sim)")
    print(f"  ⇒ ω₁ = {om_phys:.3e} rad/s   vs the electron ZBW"
          f" 2m_ec²/ℏ = {zbw:.3e} rad/s")
    print(f"  ⇒ a factor {zbw / om_phys:.1f} SHORT — and with the N-6a"
          " RIGIDITY finding this cannot be")
    print("  closed by adding dressing energy: the gap is STRUCTURAL —"
          " candidate physics: the V-on/")
    print("  Faber-r₀ core sector (the true mass family), the faithful"
          " kinetic (5d: frequency-level),")
    print("  or postulate P2 itself (the weakest link; the length-anchor"
          " alternative needs the V-on core).")
    print("  Recorded as the first measurement, anchor-dependence explicit"
          " — NO pass/fail claimed.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
