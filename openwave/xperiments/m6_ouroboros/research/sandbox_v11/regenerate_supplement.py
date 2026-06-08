"""Regenerate the DM-paper supplement from the canonical neutral-chaoiton profile.

Triggered by Werbos's 2026-06-08 normalization question (relayed from DeepSeek):
the raw mixing integral M_JA(0) is 271x the draft value, which -- plugged into
g_JA_eff = M_JA / V_Yukawa -- gives an unphysical mixing fraction (>1) and would
push sigma_p above the LZ limit. This script:

  1. recomputes the amplitude-INVARIANT quantities that actually carry the
     spatial physics (F(q), S_vert) -- both are ratios in which the beta scale
     cancels, so the 8x/271x amplitude/extent difference does NOT propagate;
  2. demonstrates the blow-up of the raw (non-invariant) g_JA_eff so the
     resolution is documented, not hidden;
  3. writes a normalization_resolution.md the paper can cite for an
     end-to-end-consistent Section 5.

All numbers reproduce from v11_canonical_beta_profile.csv. NumPy 2.x (trapezoid).
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

trap = np.trapezoid

# --- constants (canonical anchors, 0d_canonical.md sec 5) ---
HBARC = 197.3269804  # MeV*fm
M_E = 0.5109989461  # MeV
R_PHYS = HBARC / M_E  # 386.16 fm
M_J = 0.618427  # MeV
R_PROTON = 0.84  # fm (proton charge radius)
ALPHA_J = 1.21
ALPHA_EM = 1.0 / 137.036

HERE = __file__.rsplit("/", 1)[0]
CSV = f"{HERE}/v11_canonical_beta_profile.csv"
OUT = f"{HERE}/dm_paper_supplement"


def load_profile():
    rows = []
    with open(CSV) as f:
        for line in f:
            s = line.strip()
            if not s or s[0] in "#\"" or s.startswith("r_natural"):
                continue
            p = s.split(",")
            if len(p) == 4:
                rows.append([float(x) for x in p])
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]  # r_nat, r_fm, beta, dbeta_dr


def j0(x):
    return np.where(np.abs(x) < 1e-12, 1.0, np.sin(x) / np.where(x == 0, 1, x))


def form_factor(beta, r_fm, q_MeV):
    x = (q_MeV / HBARC) * r_fm
    return trap(beta * j0(x) * r_fm**2, r_fm) / trap(beta * r_fm**2, r_fm)


def main():
    r_nat, r_fm, beta, dbeta = load_profile()
    ipk = int(np.argmax(beta))
    beta_peak = beta[ipk]
    B0_nat = dbeta[0]
    B0_fm = B0_nat / R_PHYS

    # --- amplitude-INVARIANT small-r vertex suppression ---
    beta_084 = B0_fm * R_PROTON
    ratio = beta_084 / beta_peak
    S_vert = ratio**2

    # --- raw (non-invariant) mixing integral, for the record ---
    I_tab = trap(beta * r_fm**2, r_fm)
    I_smallr = B0_fm * r_fm[0] ** 4 / 4.0  # analytic [0, r_fm[0]] with beta~B0 r
    I_total = I_tab + I_smallr
    M_JA = 4 * np.pi * I_total
    V_Yukawa = (4 * np.pi / 3) * (HBARC / M_J) ** 3
    g_raw = M_JA / V_Yukawa  # >>1 : the bug

    # --- amplitude-invariant alternative (J-norm normalized) ---
    N_J = trap(beta**2 * r_fm**2, r_fm)  # scales as amp^2

    # --- form factor curve ---
    q_grid = [0, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.2, 0.3, 0.4,
              0.47, 0.5, 0.7, 1.0, 3.0, 10.0, 30.0, 100.0]
    Fq = {q: form_factor(beta, r_fm, q) for q in q_grid}

    # ---------------- write Fq_curve.csv ----------------
    with open(f"{OUT}/Fq_curve.csv", "w") as f:
        f.write("q_MeV,F_q\n")
        for q in q_grid:
            f.write(f"{q:g},{Fq[q]:.6e}\n")

    # ---------------- write beta_small_r_comparison.csv ----------------
    mask = (r_fm >= 1) & (r_fm <= 7800)  # tabulated small-r region
    with open(f"{OUT}/beta_small_r_comparison.csv", "w") as f:
        f.write("r_fm,beta_numerical,beta_linear_fit_B0r\n")
        for rr, bb in zip(r_fm[mask][:200], beta[mask][:200]):
            f.write(f"{rr:.4f},{bb:.6e},{B0_fm*rr:.6e}\n")

    # ---------------- plots ----------------
    plt.figure(figsize=(7, 4.5))
    qs = np.array([q for q in q_grid if q > 0])
    fs = np.array([Fq[q] for q in q_grid if q > 0])
    plt.semilogx(qs, fs, "o-", lw=1.4, ms=4)
    plt.axhline(0, color="k", lw=0.6)
    plt.axvline(3e-4, color="g", ls="--", lw=1, label="q_gal = 3e-4 MeV  (F=1.0000)")
    plt.axvline(M_J, color="r", ls="--", lw=1, label=f"m_J = {M_J:.3f} MeV")
    plt.xlabel("q (MeV)")
    plt.ylabel("F(q)")
    plt.title("Neutral chaoiton form factor F(q) (canonical beta profile)")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{OUT}/Fq_curve.png", dpi=130)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.plot(r_fm[mask][:200], beta[mask][:200], "o", ms=3, label="numerical beta")
    plt.plot(r_fm[mask][:200], B0_fm * r_fm[mask][:200], "-", lw=1.2,
             label=f"linear fit beta = B0 r, B0 = {B0_fm:.3e}/fm")
    plt.xlabel("r (fm)")
    plt.ylabel("beta(r)")
    plt.title("Small-r linearity (l=1 p-wave): beta ~ B0 r")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{OUT}/beta_small_r_comparison.png", dpi=130)
    plt.close()

    # ---------------- normalization_resolution.md ----------------
    md = f"""# Mixing-integral normalization resolution (2026-06-08)

Response to Werbos/DeepSeek's 2026-06-08 question: the raw mixing integral
M_JA(0) computed on the canonical profile is 271x the DM-paper draft value, which
plugged into `g_JA_eff = M_JA / V_Yukawa` gives a mixing fraction > 1 (unphysical)
and a cross section above LZ. This file records why that is a formula artifact, not
a physics exclusion. All values reproduce from `v11_canonical_beta_profile.csv`.

## Canonical anchors

| Quantity | Value |
| --- | --- |
| hbar*c | {HBARC} MeV*fm |
| m_e | {M_E} MeV |
| R_phys = hbar*c/m_e | {R_PHYS:.4f} fm |
| m_J | {M_J} MeV |
| proton charge radius | {R_PROTON} fm |

## Amplitude-INVARIANT spatial quantities (these carry the physics)

Both are ratios in which the overall beta scale cancels exactly, so the 8x
amplitude difference (peak 0.697 vs draft 0.085) and the 271x integral difference
do NOT propagate into sigma_p.

| Quantity | Value | Note |
| --- | --- | --- |
| beta_peak | {beta_peak:.5f} @ r = {r_fm[ipk]:.1f} fm | physical amplitude, fixed by a = m_J/(2 sqrt(g)) |
| origin slope B0 | {B0_fm:.4e} /fm ({B0_nat:.5f} natural) | l=1 regularity beta ~ B0 r |
| beta(0.84 fm) | {beta_084:.4e} | = B0 * 0.84 |
| ratio beta(0.84)/beta_peak | {ratio:.4e} | amplitude cancels |
| **S_vert = ratio^2** | **{S_vert:.4e}** | draft used 1.54e-4 (wrong peak) -> 63x too large |
| F(q_gal = 3e-4 MeV) | {Fq[3e-4]:.7f} | normalized F(0)=1; ~1 confirms spatial-only suppression |
| F(q = 3 MeV) | {Fq[3.0]:.3e} | strongly suppressed at LZ scales |

## Raw (NON-invariant) g_JA_eff -- the bug, for the record

`g_JA_eff = M_JA / V_Yukawa` is linear in the beta amplitude (M_JA is proportional
to beta; V_Yukawa is a fixed geometric volume), so it is NOT a valid mixing
fraction. With the physical amplitude it explodes:

| Quantity | Value |
| --- | --- |
| int beta r^2 dr | {I_total:.4e} fm^3  (nat {trap(beta*r_nat**2, r_nat):.4f}) |
| M_JA(0) = 4*pi*int | {M_JA:.4e} fm^3  (nat {4*np.pi*trap(beta*r_nat**2, r_nat):.3f}) |
| V_Yukawa = (4pi/3)(hbar c/m_J)^3 | {V_Yukawa:.4e} fm^3 |
| g_JA_eff RAW | {g_raw:.2f}  (>>1, unphysical) |
| ratio vs draft M_JA (4.23e7) | {M_JA/4.23e7:.0f}x  =  {beta_peak/0.08487:.2f}x amplitude  *  {(M_JA/4.23e7)/(beta_peak/0.08487):.1f}x extent |

## Fix (Werbos side): make g_JA_eff amplitude-invariant

Option A -- Lagrangian-level: the J-A mixing is set by the fixed J^mu A_mu coupling
and m_J^2, independent of the soliton amplitude; beta enters only via F(q) and S_vert.

Option B -- normalize M_JA by the soliton's own J-norm sqrt(N_J), N_J = int beta^2 r^2 dr
= {N_J:.4e} fm^3 (scales as amp^2), instead of the external V_Yukawa, so the scale cancels.

## Open structural point (l=1)

The chaoiton is l=1 (p-wave). M_JA(0) = 4*pi*int beta r^2 dr implicitly projects
onto l=0; for a pure l=1 field the angular integral of Y_1 over the sphere vanishes,
so the true zero-momentum monopole coupling is zero -- arguably the meaning of
"neutral." If so, M_JA(0) is the wrong object and the surviving coupling is
dipole/higher (automatically suppressed), with S_vert as the spatial stand-in.

Regenerate: `python regenerate_supplement.py`
"""
    with open(f"{OUT}/normalization_resolution.md", "w") as f:
        f.write(md)

    # ---------------- console summary ----------------
    print("REGENERATED supplement in", OUT)
    print(f"  beta_peak           = {beta_peak:.5f} @ {r_fm[ipk]:.1f} fm")
    print(f"  B0                  = {B0_fm:.4e} /fm")
    print(f"  S_vert (corrected)  = {S_vert:.4e}   [draft 1.54e-4, 63x too large]")
    print(f"  F(q_gal)            = {Fq[3e-4]:.7f}")
    print(f"  M_JA(0)             = {M_JA:.4e} fm^3  ({M_JA/4.23e7:.0f}x draft)")
    print(f"  g_JA_eff RAW        = {g_raw:.2f}  (non-invariant -> the bug)")
    print("  files: Fq_curve.csv/.png, beta_small_r_comparison.csv/.png,")
    print("         normalization_resolution.md")


if __name__ == "__main__":
    main()
