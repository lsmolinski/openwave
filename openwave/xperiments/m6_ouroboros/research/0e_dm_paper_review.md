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
below LZ). As of Werbos's 2026-06-08 reply, the S_vert correction and all attribution
fixes are **accepted**, and the exact mixing integral + β convention have been **sent**.

**2026-06-08 (2nd round) — normalization resolved.** Werbos/DeepSeek flagged that the
raw mixing integral is 271× the draft, so `g_JA_eff = M_JA/V_Yukawa` → 84.3 (>1,
unphysical) → σ_p ~10⁻¹⁹ (excluded). **Diagnosis: the formula is the bug, not the
profile.** `g_JA_eff = M_JA/V_Yukawa` is *linear in the β amplitude* (M_JA ∝ β;
V_Yukawa is a fixed geometric volume) so it is not a valid rescaling-invariant mixing
fraction. The β amplitude is **physical, not a free normalization** — the chaoiton is
neutral (Q_CS=0, no charge to normalize to); the scale is pinned by the field equation
via a = m_J/(2√g) (the same solution giving m_χ = 0.460). A correctly-built σ_p depends
on β **only through amplitude-invariant ratios** — F(q) (F(0)=1 normalization,
F(q_gal)=0.99999985) and S_vert (ratio of two β-values, corrected 2.44×10⁻⁶) — so the
8×/271× has **zero** effect on σ_p; corrected S_vert pushes σ_p *down* ~63× to
~1×10⁻³⁴ cm² (further below LZ). Fix on Werbos's side: define g_JA_eff Lagrangian-level
(fixed JᵘAᵘ + m_J², amplitude-independent) or as an amplitude-cancelling ratio (normalize
M_JA by the soliton's own J-norm √∫β²r²dr, not external V_Yukawa). **Structural flag:**
the chaoiton is l=1 — M_JA(0)=4π∫βr²dr projects onto l=0, but ∫Y₁dΩ=0, so the true
zero-momentum monopole coupling vanishes (= what "neutral" should mean); M_JA(0) is
likely the wrong object and the surviving coupling is dipole/higher, with S_vert as the
spatial stand-in. Full record + reproducible numbers in
[`sandbox_v11/dm_paper_supplement/normalization_resolution.md`](sandbox_v11/dm_paper_supplement/normalization_resolution.md).

**2026-06-09 — adopted + closed.** Werbos published v11 to Zenodo (DOI 20612600) and
sent a June-9-AM cleanup that adopted our real profile, S_vert = 2.44×10⁻⁶, and the l=1
monopole-vanishing selection rule (cross-section now dipole). We caught one new error (the
r³ dipole moment was tail-truncated: full ∫βr³dr = 1.58×10¹² fm⁴ ⇒ R_dipole = 1.73×10³ fm,
σ_p ≈ 9×10⁻⁴¹ cm²); DeepSeek accepted it. We flagged a possible S_vert × |F_dipole|²
double-count; DeepSeek disagrees (asserts independent factors, will clarify §5.5). σ_p is
far below LZ under either reading. **Griesi is not replying further** — our verification
role is complete; remaining residuals are Werbos-side. Details in §8.

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

## 4. The 271× mixing-integral gap — RESOLVED (normalization is a formula artifact)

The absolute mixing integral is 271× the draft. **Root cause: `g_JA_eff = M_JA/V_Yukawa`
is not a valid mixing fraction** — it is linear in the β amplitude (M_JA ∝ β; V_Yukawa is
a fixed geometric volume independent of β), so it is not invariant under rescaling of β,
which is precisely the property DeepSeek (correctly) says a physical mixing fraction must
have. Plugging the physical amplitude in gives g_JA_eff = 84.3 and σ_p⁰ ≈ 2.3×10⁻¹⁹ cm²
(excluded) — a signal the *estimator* is wrong, not that the chaoiton is excluded.

| Quantity | Paper §5 (draft) | Canonical profile — fm | natural | |
| --- | --- | --- | --- | --- |
| ∫β r² dr | 3.37×10⁶ fm³ | 9.134×10⁸ fm³ | 15.863 | |
| M_JA(0) = 4π∫β r² dr | 4.23×10⁷ fm³ | 1.148×10¹⁰ fm³ | 199.34 | 271× draft |
| g_JA_eff = M_JA/V_Yukawa | 0.311 | **84.35 (>1, unphysical)** | — | non-invariant ⇒ the bug |
| σ_p⁰ (∝ g_JA_eff⁴) | 4.25×10⁻²⁹ cm² | **2.3×10⁻¹⁹ cm² (excluded)** | — | if raw value used |

**Decomposition of the 271×:** 8.21× is amplitude (peak 0.697 vs draft 0.085), the
remaining **33× is shape/extent** — so the draft's M_JA did *not* come from the canonical
profile even amplitude-matched. The draft's 0.311 / 6.5×10⁻³³ are convention artifacts
that happened to land below LZ, not robust numbers.

**β normalization convention:** β(r) is the raw, dimensionless BVP soliton amplitude
(peak 0.697), **not** normalized to unit norm or to a charge — the chaoiton is neutral
(Q_CS = 0, nothing to normalize to), so the amplitude is fixed by the field equation via
the v11 scaling symmetry a = m_J/(2√g) (the same solution giving m_χ = 0.460), not by a
choice. It **cannot** be rescaled down to recover 0.311. Near the origin it is the l=1
p-wave β = B₀r, B₀ = 0.5 natural (= 1.295×10⁻³/fm). R_phys = ℏc/m_e = 386.16 fm.

**The resolution (amplitude-invariant σ_p):** a correctly-built cross-section depends on
β only through ratios where the scale cancels — F(q) (F(0)=1; F(q_gal)=0.99999985) and
S_vert = (β(0.84)/β_peak)² (corrected 2.44×10⁻⁶). The 8×/271× has **zero** effect; the
corrected S_vert *lowers* σ_p ~63× to ~1×10⁻³⁴ cm² (further below LZ, conclusion
strengthened).

**Fix on Werbos's side (sent):** define g_JA_eff so the amplitude drops out — either
Lagrangian-level (J–A mixing set by the fixed JᵘAᵘ coupling and m_J², amplitude-independent;
β enters only via F(q) and S_vert), or as an explicit amplitude-cancelling ratio (normalize
M_JA by the soliton's own J-norm √∫β²r²dr = 4.4×10⁵ fm³ scale, not the external V_Yukawa).

**Structural flag (likely the cleanest answer):** the chaoiton is l=1 (p-wave). M_JA(0) =
4π∫β r² dr implicitly projects onto l=0, but ∫Y₁ dΩ = 0 over the sphere, so the true
zero-momentum **monopole coupling vanishes** — arguably the meaning of "neutral." If so
M_JA(0) is the wrong object entirely; the surviving coupling is dipole/higher (auto-suppressed),
with S_vert the spatial stand-in. Posed to Werbos + DeepSeek as a question (his framework
to confirm).

**Form-factor fine detail (minor, grid-sensitive):** F(q_gal) ≈ 1 is confirmed (the
load-bearing claim). The high-q tail differs from the draft (zero crossing ours ~0.67
MeV vs draft ~0.47 MeV; F(q_LZ) ours −6.1×10⁻⁵ vs draft −1.15×10⁻³) — both confirm F
is strongly suppressed at LZ scales; the exact values are sensitive to the radial
grid/cutoff. Use the attached `Fq_curve.csv` as the reference if consistency matters.

---

## 5. Attribution issues (factual — flagged, Werbos ACCEPTED 2026-06-08)

Werbos's 2026-06-08 reply confirmed all three: he will revise §5.9/§4 to credit the
six-peak modulation to astrophysical halo models (Gaia streams) not OpenWave data, update
the Acknowledgments to the correct division of labour (OpenWave: BVP solution + β(r)
profile + scaling symmetry; Werbos: cross-sections), and add a plain-language summary
linking each §5 number to its repo script.

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
| `normalization_resolution.md` | the 2026-06-08 mixing-integral resolution + 2026-06-09 dipole section: corrected S_vert, raw-g_JA_eff blow-up, amplitude-invariant fix, l=1 monopole flag, full-tail ∫βr³dr / R_dipole / σ_p, S_vert×\|F_dipole\|² double-count question |
| `beta_small_r_comparison.csv` / `.png` | r_fm, numerical β, linear fit β = B₀r (B₀ = 1.295×10⁻³/fm); the l=1 small-r linearity |
| `Fq_curve.csv` / `.png` | q (MeV) vs F(q) monopole form factor (F(0)=1), computed F(q) = ∫β j₀(qr) r² dr / ∫β r² dr |
| `Fq_dipole.csv` | q (MeV) vs F_dipole(q) l=1 form factor (F(0)=0), ∫β j₁(qr) r² dr / ∫β r² dr — the 2026-06-09 dipole channel |
| (regen) `regenerate_supplement.py` | reproducible regenerator for all of the above from the canonical CSV (monopole + dipole blocks) |
| (source) `v11_canonical_beta_profile.csv` | the canonical β(r) profile (r_natural, r_fm, β, dβ/dr) |

**Recompute method (reproducible):** parse the canonical CSV; ℏc = 197.327 MeV·fm,
R_phys = ℏc/m_e = 386.16 fm; peak/slope/ratio read directly; S_vert = (β(0.84 fm)/β_peak)²;
F(q) via Simpson/trapezoid with argument x = q·r/ℏc. The one-off script lived inline
(2026-06-08); regenerate from the CSV if needed.

---

## 7. Open items / next steps

| Item | Status |
| --- | --- |
| Corrected S_vert (2.44×10⁻⁶) + β_peak (0.697) in §5.4/5.5 | ✅ Werbos accepted 2026-06-08; σ_p ~1×10⁻³⁴ cm² (further below LZ), conclusion unchanged |
| Mixing integral ∫β r² dr + 4π form + β convention | ✅ sent 2026-06-08 (9.134×10⁸ fm³ / 1.148×10¹⁰ fm³; nat 15.863 / 199.34) |
| **Normalization (271× gap) — g_JA_eff non-invariance** | ✅ resolved + reply sent 2026-06-08 (2nd round): formula is linear in β ⇒ not a valid mixing fraction; σ_p rides on amplitude-invariant F(q) + S_vert only; raw plug-in (g=84) is the bug, not the profile. See §4 + `normalization_resolution.md` |
| §5.9/§4 stream attribution + Acknowledgments + plain-language summary | ✅ Werbos accepted 2026-06-08; he will revise. ⚠️ Desktop docx still credits Griesi with "computing the key cross-sections" — re-flag (cross-sections are Werbos-side) |
| Werbos to redefine g_JA_eff amplitude-invariantly (Lagrangian-level or √∫β²r²dr norm) → finalize α_eff, σ_p⁰ | ⚠️ partial: June-9-AM version drops the g_JA_eff line but σ_p⁰ still rides on α_eff = 8.56×10⁻⁴ (= old g_JA_eff = 0.311); the absolute scale still carries the normalization. See §8 |
| l=1 monopole vanishing (∫Y₁dΩ=0) ⇒ M_JA(0) is wrong object? | ✅ Werbos ADOPTED 2026-06-09 (§5.2 derives ∫dΩ Y₁ₘ=0, F(0)=0); cross-section now dipole. See §8 |
| **2026-06-09: ∫βr³dr truncated → R_dipole / F_dipole / σ_p off ~14×** | ✅ DeepSeek ACCEPTED 2026-06-09 12:02 ("the last numerical loose end"); will update to full-tail ∫βr³dr = 1.58×10¹² fm⁴, R_dipole = 1.73×10³ fm, σ_p ≈ 9×10⁻⁴¹ cm². See §8 |
| **2026-06-09: S_vert × \|F_dipole\|² possible double-count** | ⚠️ DeepSeek DISAGREES (2026-06-09 12:02): claims independent (F_dipole = global/large-r soliton property; S_vert = local amplitude at nucleon), will add a clarifying sentence in §5.5. Our read: soft (S_vert pins nucleon at r=0.84 fm from soliton center = its own charge radius, not its position; the form factor already integrates over r). Below LZ either way (9×10⁻⁴¹ vs 3.8×10⁻³⁵). **Griesi not replying further; closed unless Werbos asks.** See §8 |
| Open canonical questions Q48 (H = 8 m_J² Q identity), Q49 (α_JN definition) | see [`0d_canonical.md §8`](0d_canonical.md) |

---

## 8. 2026-06-09 round — Werbos adopted the corrections + dipole rewrite (new r³ error)

Werbos published v11 to Zenodo (DOI 10.5281/zenodo.20612600, dated 2026-06-09) then
emailed a cleaned-up "DarkJune9AM.docx" ("DeepSeek says it cleaned up everything …
your views certainly welcome"). Sequence + verification:

**The published Zenodo v11 was still inconsistent** (it predates the cleanup): §5.4
used an *analytic piecewise stand-in* (peak 0.1411, S_vert 5.59×10⁻⁵) per its bundled
`v8_numeric_summary.txt` ("analytic piecewise model"), **not** our canonical CSV (peak
0.697) — even though the CSV ships with the record. Abstract / §5.2 / §5.6 / §5.7 still
carried the old v8 numbers (S_vert 1.54×10⁻⁴, M_JA 4.23×10⁷, σ_p 6.5×10⁻³³), so three
S_vert values coexisted. Flagged to Werbos.

**The June-9-AM cleanup fixed the big things** (verified against the canonical CSV):

| Item | Status |
| --- | --- |
| Real canonical profile (peak 0.697, ∫βr²dr = 9.134×10⁸ fm³) | ✅ matches CSV exactly |
| S_vert = 2.44×10⁻⁶ used consistently | ✅ |
| l=1 monopole-vanishing selection rule (§5.2: ∫dΩ Y₁ₘ = 0, F(0)=0) | ✅ adopted our flag; cross-section now dipole |
| Three-S_vert mess removed (old 1.54e-4 / 5.59e-5 / 6.5e-33 / 0.311 / 0.1411 gone from body) | ✅ |
| Abstract → σ_p = 6.7×10⁻⁴² | ✅ |

**New arithmetic error — the r³ dipole moment is tail-truncated.** R_dipole = ∫βr³dr/∫βr²dr
needs the full exponential tail (β decays with ~750 fm scale out to ~10⁴ fm). The paper's
∫βr³dr = 4.28×10¹¹ fm⁴ matches truncating near r ≈ 1450 fm.

| Quantity | June-9-AM paper | Correct (full tail, CSV) | Off by |
| --- | --- | --- | --- |
| ∫βr³ dr | 4.28×10¹¹ fm⁴ | **1.58×10¹² fm⁴** | 3.7× |
| R_dipole | 4.69×10² fm | **1.73×10³ fm** | 3.7× |
| \|F_dipole(q_gal)\|² | 6.5×10⁻⁸ | **8.90×10⁻⁷** (exact j₁ = small-q approx) | ~14× |
| σ_p (their factorization) | 6.7×10⁻⁴² cm² | **9.2×10⁻⁴¹ cm²** | ~14× |

∫βr²dr matches theirs exactly ⇒ the discrepancy is purely the r³ tail. Still ~11 orders
below LZ — conclusion unchanged. q_gal = m_χv/ℏc = 1.63×10⁻⁶ fm⁻¹.

**Physics question — S_vert × \|F_dipole\|² may double-count.** Both trace to β ∼ B₀r at
small r (the l=1 BC). The dipole form factor already integrates β over all radii, so the
extra S_vert factor looks like the l=1 suppression applied twice. Clean readings:
σ_p ≈ σ_p⁰ × \|F_dipole(q)\|² ≈ 4×10⁻³⁵ cm² (form-factor picture) **or** a local-vertex
picture, not the product. Posed to Werbos + DeepSeek; below LZ either way.

**Residual normalization (lower priority):** σ_p⁰ still rides on α_eff = 8.56×10⁻⁴
(Table 1, Appendix B) = old monopole g_JA_eff = 0.311. The dipole rewrite fixed the
angular structure; the dimensionful prefactor still carries the §4 normalization.

Reply sent 2026-06-09 (positive lead + r³ fix + double-count question; no future-work
offers per Griesi's instruction). Numbers reproducible via
[`sandbox_v11/regenerate_supplement.py`](sandbox_v11/regenerate_supplement.py) (dipole block).

**DeepSeek reply (2026-06-09 12:02) + close.** DeepSeek accepted the r³ correction
("we'll update with your full-tail value … the last numerical loose end") and will
propagate R_dipole → 1.73×10³ fm, |F_dipole|² → 8.9×10⁻⁷, σ_p → ~9×10⁻⁴¹ cm². On the
double-count it **disagreed**: it reads F_dipole as a global soliton property (large-r
dominated) and S_vert as a local amplitude at the nucleon position, "independent
factors," and will add a clarifying sentence in §5.5.

Our read (recorded, not sent): the rebuttal is physically soft. S_vert = (β(0.84 fm)/β_peak)²
evaluates β at r = 0.84 fm **from the soliton center**, but 0.84 fm is the proton's own
charge radius, not its position inside the soliton; in standard DM-nucleon scattering the
nucleon's relative position is integrated over, and that integral **is** F_dipole(q)
(which already runs over all r, large-r dominated only because the small-r integrand is
r⁴-suppressed). So the spatial structure is plausibly counted once by F_dipole, and S_vert
at a single radius is not obviously a second independent factor. Not certain (a genuine
two-scale EFT factorization could differ), and it does **not** change the conclusion:
σ_p = 9×10⁻⁴¹ cm² (both factors) or 3.8×10⁻³⁵ cm² (form-factor only), both far below LZ.

**Griesi decision (2026-06-09): not replying further.** It is Werbos's model, he/DeepSeek
engaged the question and will make §5.5 explicit, and nothing observable hinges on it. The
thread is closed on our side unless Werbos asks for more. Our verification role is complete:
profile, S_vert, l=1 selection rule, and the r³ moment are all verified against the canonical
CSV; the only residuals are Werbos-side (the §5.5 factorization wording + the α_eff/g_JA_eff
absolute-scale normalization from §4).

---

## 9. Cross-references

- [`0d_canonical.md`](0d_canonical.md) — canonical numerical spec (the source of the verified numbers)
- [`0c_sandbox_v11.md`](0c_sandbox_v11.md) — (g, λ) scan + scaling symmetry + β-profile export
- [`sandbox_v11/v11_canonical_beta_profile.csv`](sandbox_v11/v11_canonical_beta_profile.csv) — the canonical profile
- [`sandbox_v11/dm_paper_supplement/normalization_resolution.md`](sandbox_v11/dm_paper_supplement/normalization_resolution.md) — mixing-integral + dipole resolution
- Werbos DM paper — Zenodo 20612600 (v11, 2026-06-09); Werbos Paper I (Ouroboros Lagrangian v11) — Zenodo 20357670
