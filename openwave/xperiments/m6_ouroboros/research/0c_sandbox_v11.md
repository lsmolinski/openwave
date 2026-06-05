# 2026-06-05 — Sandbox v11: (g, λ) scan answers Paul/DeepSeek's "urgent problem"

Triggered by Paul Werbos's 2026-06-05 email ("an urgent problem after all"),
relaying DeepSeek's request: run a parameter scan over (g, λ) in the neutral
sector to find the values yielding m_χ = 0.460 MeV and m_J = 0.618 MeV,
deliver a definitive soliton profile β(r) for the J-A mixing integral
(α_JN calculation), and resolve the "BVP inconsistency" documented in the
DM paper v5a Technical Note (Zenodo 20549492).

**Net result:** the scan is analytically solved — the neutral ODE has an
**exact scaling symmetry** that collapses the entire (g, λ < 0)
quarter-plane into rescalings of ONE universal profile. m_J = 0.6184 MeV is
a parameter-free invariant; m_χ = 0.460 MeV fixes λ* = −0.2647 at any g;
g = 1.0 comes from electron universality (Q45). The canonical (g, λ) point
is identical to the v9/v10 point (g=1, B0=0.5). The "BVP inconsistency" is
a cold-seed Newton-basin artifact reproduced exactly (peak ~10⁻¹¹), not a
physics discrepancy. Bonus: v10's Q47 invariance (m_J/η = 1.21024) is now
ANALYTICALLY EXPLAINED, plus a new candidate exact identity H = 8 m_J² Q.

---

## The request (DeepSeek via Paul, verbatim asks)

| # | Ask | v11 answer |
| --- | --- | --- |
| 1 | Determine correct neutral-sector (g, λ) from first principles | λ* = −0.264697 (any g); g = 1.0 by electron universality; masses CANNOT fix g (exact degeneracy) |
| 2 | Definitive β(r) profile for the J-A mixing integral (α_JN) | `sandbox_v11/v11_canonical_beta_profile.csv` — 4000 points, natural + fm units, full metadata header |
| 3 | Resolve the BVP inconsistency (near-trivial "Griesi-1" vs wrong-mass shooting) | Cold-seed trivial-basin collapse reproduced at peak ~10⁻¹¹ (their exact number); same code + scaling seed → soliton. Shooting failure = Q42 windowing artifact (already documented) |

---

## The scaling symmetry (the central result)

The neutral-sector ODE (Q43+Q44 canonical form, λ ≡ −m_J²):

```text
β'' + (2/r) β' - (2/r²) β + λβ + 4g β³ = 0,   λ < 0
```

Substituting β(r) = a·u(x) with x = m_J·r, a = m_J/(2√g), m_J = √(−λ):

```text
u'' + (2/x) u' - (2/x²) u - u + u³ = 0
```

— a **universal, parameter-free equation**. Every (g, λ < 0) ground state is
an exact rescaling of the single profile u(x). Consequences, each verified
numerically on an 18-point grid (g ∈ {0.5, 1.0, 2.0} × λ ∈ [−1.0, −0.09]):

| Identity | Verified |
| --- | --- |
| η = m_J·η_u with η_u = 0.826287 universal | ✅ spread 1.8×10⁻⁵ across grid |
| m_J_phys = (m_J/η)·m_e = m_e/η_u = **0.6184 MeV** for ALL (g, λ) | ✅ constant column in scan table |
| m_χ = (−λ)^{3/2}·η_u·K_u·m_e = 3.3778·(−λ)^{3/2} MeV — depends ONLY on λ | ✅ g-rows identical at fixed λ |
| H/Q = K_u·m_J² with K_u = **7.99999 ≈ 8** | ✅ spread 1.9×10⁻⁴ — candidate EXACT virial/Pohozaev identity H = 8 m_J² Q |
| Amplitude β ~ 1/√g (g = pure amplitude scale; no mass observable moves) | ✅ peak column scales as g^{−1/2} |
| At fixed B0: m_J ∝ g^{1/4} (v10 chart conversion) | ✅ reproduces v10 g-scan to 5 digits |

**This explains Q47 analytically**: v10's empirically-discovered invariance
m_J_corrected = m_J/η = 1.21024 is just 1/η_u — a property of the universal
profile, not a coincidence and not a tunable. DeepSeek's "should scale to
1.0" heuristic could never land; 1.21024 = 1/η_u is exact.

---

## Scan results (Step 2 of `m6_v11_neutral_gl_scan.py`)

Excerpt at g = 1.0 (full 18-point table in `v11_output.txt`):

| λ | m_J_raw | B0_out | peak β | nodes | m_χ (MeV) | m_J (MeV) |
| --- | --- | --- | --- | --- | --- | --- |
| −0.090 | 0.300 | 0.170 | 0.407 | 0 | 0.0912 | 0.6184 |
| −0.160 | 0.400 | 0.302 | 0.542 | 0 | 0.2162 | 0.6184 |
| **−0.2647** | **0.5145** | **0.5000** | **0.697** | **0** | **0.4599** | **0.6184** |
| −0.360 | 0.600 | 0.680 | 0.813 | 0 | 0.7296 | 0.6184 |
| −0.640 | 0.800 | 1.209 | 1.084 | 0 | 1.7294 | 0.6184 |
| −1.000 | 1.000 | 1.889 | 1.355 | 0 | 3.3778 | 0.6184 |

Locus m_χ = 0.460 MeV (brentq on λ, per g):

| g | λ* | m_J_raw | Note |
| --- | --- | --- | --- |
| 0.5 | −0.264696 | 0.514487 | identical to analytic prediction −0.264695 |
| 1.0 | −0.264697 | 0.514487 | **canonical** — B0_out = 0.50005 = the v9/v10 point |
| 2.0 | −0.264698 | 0.514488 | g-independence confirmed |

**Key honest framing:** the two masses cannot independently determine
(g, λ) — m_J is invariant (carries no information) and m_χ only fixes λ.
The scan reveals the degeneracy rather than breaking it. g is fixed by
physics input (Q45: same Lagrangian as electron ⇒ g = 1.0), and g DOES
matter for α_JN because the mixing integral is linear in β (amplitude
∝ 1/√g).

---

## BVP inconsistency — resolved (Step 5)

The v5a Technical Note claims: "The Griesi-1 code (3D spherical l=1, Robin
BC) finds only near-trivial solutions with peak amplitude ~10^-11 at all B0
initializations." Reproduction at fixed (g=1, λ=−1):

| Run | Seed | Peak β | Outcome |
| --- | --- | --- | --- |
| Fixed-λ, correct l=1 BC | cold `B0·r·e^{−r}`, B0=0.1 | 6.1×10⁻¹² | TRIVIAL collapse |
| Fixed-λ, correct l=1 BC | cold, B0=0.5 | 2.9×10⁻²⁴ | TRIVIAL collapse |
| Fixed-λ, correct l=1 BC | cold, B0=1.0 | 5.1×10⁻¹¹ | TRIVIAL collapse — **the Note's exact ~10⁻¹¹** |
| Fixed-λ, correct l=1 BC | scaling seed a·u(m_J r) | 1.3549 | ✅ SOLITON, 0 nodes |
| Fixed-λ, template BC β'(R_MIN)=0 | cold | 1.7×10⁻²⁴ | TRIVIAL collapse (v9 phase 2 observation) |
| Fixed-λ, template BC β'(R_MIN)=0 | scaling seed | 1.355 | soliton — BC error at R_MIN is O(R_MIN²); the trivial Newton basin is the DOMINANT failure mechanism, not the BC |

Diagnosis: at fixed λ the trivial solution β ≡ 0 satisfies both BCs, and
`solve_bvp`'s Newton iteration falls into it unless seeded at the nonlinear
amplitude a = m_J/(2√g). The required amplitude at λ=−1, g=1 is a = 0.5
with B0_req = 1.89 — cold seeds at B0 ≤ 1 are simply in the wrong basin.
**Whoever ran "Griesi-1" ran a fixed-λ formulation with cold seeds.** Our
actual v9/v10 production code (`neutral_bvp_solver_mJ_free.py`) never had
this failure: it pins the amplitude via B0 as a hard BC and frees m_J as
the eigenvalue, which structurally excludes the trivial solution.

The Note's second leg ("v9 shooting code … masses of 1-4 MeV") is the IVP
windowing artifact already closed as Q42 (`0c_sandbox_v8.md`).

**Correction needed in the paper's Technical Note:** it states our
m_χ = 0.460 / m_J = 0.618 "were obtained by a separate fitting procedure to
the soliton energy functional, not by direct BVP integration." Incorrect —
they came from direct `solve_bvp` integration of the l=1 ODE (v9 phase 2)
plus DeepSeek's own Q46 η-recipe (v10). v11 closes the loop: direct
integration AT (g=1, λ*) lands on B0 = 0.500, m_χ = 0.4600, m_J = 0.6184,
C = 770.2 — direct BVP and "energy-functional calibration" are the SAME
calculation and were never inconsistent.

---

## Deliverables at the canonical point (g = 1.0, λ* = −0.264697)

| Quantity | Value | Cross-check vs v10 |
| --- | --- | --- |
| m_χ | 0.4600 MeV | ✅ 0.4599 |
| m_J | 0.6184 MeV | ✅ 0.6184 |
| C | 770.2 MeV·fm | ✅ 770 |
| B0 (recovered, was an output here) | 0.500050 | ✅ 0.5 (v9/v10 input) |
| peak β @ r | 0.6971 @ 1.781 (687.8 fm) | ✅ 0.6970 @ 1.779 |
| nodes / tail-to-peak | 0 / 1.9×10⁻⁷ | ✅ clean K₁ ground state |
| β(r) profile CSV | `sandbox_v11/v11_canonical_beta_profile.csv` | 4000 points; r in natural + fm; β, β′; metadata header |

---

## New questions raised (for Paul/DeepSeek — email v17)

| Q | Content |
| --- | --- |
| Q48 | H = 8 m_J² Q holds to 5 digits across the entire (g, λ) grid — is this an exact Pohozaev/virial identity for the 3D l=1 cubic NLS? (A proof would make m_χ(λ) = 8·η_u·(−λ)^{3/2}·m_e fully analytic except for the single universal number η_u = 0.826287.) |
| Q49 | α_JN: we can compute the J-A mixing integral M_JA(0) from the delivered profile as soon as its exact definition is specified (integrand + measure + normalization). The Note's g_JA_eff ~ 10-50 "strongly dependent on the correct lambda" should collapse to a single number at λ* = −0.2647, g = 1.0. |
| (flag) | v5a Table 1 internal errors: Yukawa range "1/m_J = 0.319 fm" should be ħc/m_J = **319 fm** (10³ off); α_J = 239 in v5a vs 1.21 in the v5 draft — needs reconciliation. |

---

## Files

| File | Content |
| --- | --- |
| `sandbox_v11/m6_v11_neutral_gl_scan.py` | The scan script: anchor regression, 18-point (g, λ) grid, λ* locus, failure-mode reproduction, scaling-law checks, profile export |
| `sandbox_v11/v11_output.txt` | Full run output |
| `sandbox_v11/v11_canonical_beta_profile.csv` | Definitive β(r) for the α_JN / M_JA(0) calculation |

---

## Cross-references

- `0c_sandbox_v10.md` — the v10 canonical point + Q47 empirical invariance this doc explains
- `0c_sandbox_v9.md` — phase 2 ground-state formulation (B0 hard BC + m_J free)
- `0c_sandbox_v8.md` — Q42 shooting/IVP windowing artifact
- `0d_canonical.md` — canonical spec, updated with the scaling symmetry
- `0b_question_tracker.md` — Q47 ANALYTICALLY RESOLVED; Q48/Q49 NEW
- DM paper v5a — Zenodo 20549492 (Technical Note this doc responds to)
