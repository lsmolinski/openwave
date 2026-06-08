# M6 — DM Paper Section-5 Verification (Werbos, June 2026)

**Purpose:** record of the numerical verification of how the OpenWave/M6 results
are reported in Paul Werbos's dark-matter paper — *"The Neutral Chaoiton as Dark
Matter"* (v8, June 2026, to be submitted to ApJ). Triggered by Werbos's 2026-06-08
request (DeepSeek flagged that Section 5 must report our results + references to
Griesi accurately). Every β-derived number in §5 was re-computed directly from the
canonical profile [`sandbox_v11/v11_canonical_beta_profile.csv`](sandbox_v11/v11_canonical_beta_profile.csv).

**Last updated:** 2026-06-08.

**Bottom line:** the core (parameters, profile shape, form factor) is accurate; one
solid correction lands in the paper's favor (S_vert 63× too large → σ_p further
below LZ); one item needs Werbos's normalization convention before it can be
finalized (the absolute mixing integral); two attribution items need fixing.

---

## 1. Paper under review

| Field | Value |
| --- | --- |
| Title | The Neutral Chaoiton as Dark Matter: A Falsifiable Candidate from the Ouroboros Lagrangian |
| Author | Paul J. Werbos (with contributions: Griesi / Claude Sonnet+Code / DeepSeek) |
| Version | v8, June 2026; Zenodo preprint → ApJ submission |
| Our contribution | the canonical neutral-chaoiton β(r) profile + the BVP solution + the scaling symmetry (sandbox_v9 → v11); cited as Reference [7] = the GitHub repo |
| §5 = where it lives | "J-Field–Nucleon Coupling: Numerical Computation" — the β-profile-derived coupling chain |

---

## 2. What is ACCURATE (confirm as-is)

| Item | Paper | Recompute from canonical profile | |
| --- | --- | --- | --- |
| m_χ / m_J / C / λ / g | 0.460 MeV / 0.618 MeV / 770 MeV·fm / −0.2647 / 1.0 | identical (CSV header: m_chi 0.460004, m_J 0.618427, C 770.16) | ✅ |
| β-profile peak location | r ≈ 688 fm | 687.8 fm | ✅ |
| small-r slope B₀ | 1.254×10⁻³ /fm | 1.295×10⁻³ /fm (R² = 1.0000; origin slope dβ/dr = 0.5 natural) | ✅ |
| β(0.84 fm) (proton radius) | 1.053×10⁻³ | 1.088×10⁻³ | ✅ |
| Form factor F(q_gal = 3×10⁻⁴ MeV) | 0.9999995 | 0.999998 | ✅ — the F≈1 assumption in §5.5 holds |
| Technical Note (BVP "inconsistency") | numerical artifact: trivial solution satisfies BCs, Newton collapses from low-amplitude seed; correct amplitude seeding β(0)∝m_J/(2√g) yields the soliton | matches the v11 "cold-seed artifact" finding (`0d_canonical.md §3.7`, `0c_sandbox_v11.md`) | ✅ |
| Repo citation [7] + URL | github.com/openwave-labs/openwave | correct, reproducible | ✅ |

---

## 3. What is WRONG and CORRECTED (solid — send)

The §5.4 nuclear-scale vertex suppression is wrong because of a normalization
mix-up; the correction is **amplitude-independent**, so it is certain.

| Quantity | Paper §5 | Corrected (canonical profile) | Factor |
| --- | --- | --- | --- |
| β_peak | 0.08487 | **0.69707** | draft 8.21× too small |
| ratio β(0.84 fm)/β_peak | 0.01241 | **1.561×10⁻³** | |
| **S_vert = ratio²** | **1.54×10⁻⁴** | **2.44×10⁻⁶** | **draft 63× too large** |

**Why the correction is certain:** S_vert is a *ratio* of two β-values from the
same profile, so it is independent of the overall β normalization. The draft's
small-r value (β(0.84) ≈ 1.05×10⁻³) is correct and matches ours; it was divided by
a peak (0.085) that is 8× smaller than the canonical peak (0.697). Under *any*
consistent normalization the canonical profile pins S_vert at 2.44×10⁻⁶. (Tell in
the draft itself: §5.4 notes the new ratio is "8× larger than the earlier
hand-extrapolation (1.56×10⁻³)" — and 1.56×10⁻³ is exactly what the canonical
profile gives.)

**Consequence:** σ_p = σ_p⁰ × S_vert drops by ~63× → roughly **1×10⁻³⁴ cm²** (if
σ_p⁰ holds) instead of 6.5×10⁻³³ — i.e. *further* below the LZ limit. The correction
**strengthens** the paper's "not excluded by direct detection" conclusion.

---

## 4. What needs RECONCILIATION (do not overwrite — await Werbos)

The absolute mixing integral does not reproduce from the raw profile, and a raw
recompute is unphysical — so Werbos's side is applying a Yukawa-volume normalization
not visible in the draft.

| Quantity | Paper §5 | Raw recompute (canonical profile) | Status |
| --- | --- | --- | --- |
| ∫β r² dr | 3.37×10⁶ fm³ | 9.13×10⁸ fm³ | ~270× larger |
| M_JA(0) = 4π∫β r² dr | 4.23×10⁷ fm³ | 1.15×10¹⁰ fm³ | ~270× larger |
| g_JA_eff = M_JA/V_Yukawa | 0.311 | raw → **> 1 (unphysical)** | needs V_Yukawa convention |

⚠️ Because α_eff (§5.3) and σ_p⁰ (§5.5) ride on M_JA via V_Yukawa, they **cannot be
finalized** until Werbos supplies the exact `V_Yukawa` definition (or the script
that produced M_JA = 4.23×10⁷). Requested in the 2026-06-08 reply.

**Form-factor fine detail (minor, grid-sensitive):** F(q_gal) ≈ 1 is confirmed (the
load-bearing claim). The high-q tail differs from the draft (zero crossing ours ~0.67
MeV vs draft ~0.47 MeV; F(q_LZ) ours −6.1×10⁻⁵ vs draft −1.15×10⁻³) — both confirm F
is strongly suppressed at LZ scales; the exact values are sensitive to the radial
grid/cutoff. Use the attached `Fq_curve.csv` as the reference if consistency matters.

---

## 5. Attribution issues (factual — flagged in the reply)

| Where | Issue | Fix |
| --- | --- | --- |
| §5.9 + §4 | "Rodrigo's exact Table 2" / "When Rodrigo supplies the precise halo model streams" | The six Gaia Dark Shard streams are Werbos's astrophysics input (O'Hare 2019), **not** OpenWave's. Remove Griesi's name. (This is why `Phi_time_series.csv` / `Phi_annual_modulation.png` were not in our repo.) |
| Acknowledgments | "discovering the scaling symmetry … computing the key cross-sections" over-credits Griesi personally | OpenWave provided the BVP solution + β(r) profile + scaling symmetry, computed on the platform with the AI agents (Claude Code/Sonnet, DeepSeek) under Griesi's engineering direction; the §5 cross-sections were Werbos-side. Don't cast Griesi as personally deriving the physics (consistent with the `feedback_role_split_platform_engineer` + `feedback_writing_role_declined` stance). |
| Abstract vs Table 1 | "three free parameters" (abstract) vs "two parameters (g, λ)" (Table 1 caption) | reconcile the count |

---

## 6. Supplementary files generated (from the canonical profile)

Written to [`sandbox_v11/dm_paper_supplement/`](sandbox_v11/dm_paper_supplement/);
these cover 4 of the 6 supplementary files Werbos listed (the other two, `Phi_*`,
are his Gaia-stream files):

| File | Contents |
| --- | --- |
| `beta_small_r_comparison.csv` / `.png` | r_fm, numerical β, linear fit β = B₀r (B₀ = 1.295×10⁻³/fm); the l=1 small-r linearity |
| `Fq_curve.csv` / `.png` | q (MeV) vs F(q) form factor, computed F(q) = ∫β j₀(qr) r² dr / ∫β r² dr |
| (source) `v11_canonical_beta_profile.csv` | the canonical β(r) profile (r_natural, r_fm, β, dβ/dr) |

**Recompute method (reproducible):** parse the canonical CSV; ℏc = 197.327 MeV·fm,
R_phys = ℏc/m_e = 386.16 fm; peak/slope/ratio read directly; S_vert = (β(0.84 fm)/β_peak)²;
F(q) via Simpson/trapezoid with argument x = q·r/ℏc. The one-off script lived inline
(2026-06-08); regenerate from the CSV if needed.

---

## 7. Open items / next steps

| Item | Status |
| --- | --- |
| Werbos to supply `V_Yukawa` definition / M_JA script | 🚧 requested 2026-06-08 — needed to finalize α_eff, σ_p⁰, σ_p |
| Apply corrected S_vert (2.44×10⁻⁶) + β_peak (0.697) in §5.4 | 🔶 sent; Werbos to integrate |
| Fix §5.9/§4 stream attribution + Acknowledgments wording | 🔶 sent |
| Open canonical questions Q48 (H = 8 m_J² Q identity), Q49 (α_JN definition) | see [`0d_canonical.md §8`](0d_canonical.md) |

---

## 8. Cross-references

- [`0d_canonical.md`](0d_canonical.md) — canonical numerical spec (the source of the verified numbers)
- [`0c_sandbox_v11.md`](0c_sandbox_v11.md) — (g, λ) scan + scaling symmetry + β-profile export
- [`sandbox_v11/v11_canonical_beta_profile.csv`](sandbox_v11/v11_canonical_beta_profile.csv) — the canonical profile
- Werbos Paper I (Ouroboros Lagrangian v11) — Zenodo 20357670
