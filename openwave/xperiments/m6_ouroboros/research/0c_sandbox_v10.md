# 2026-05-22 evening — Sandbox v10: Canonical DM paper inputs per Q45+Q46 reply

Triggered by Paul Werbos's 2026-05-22 7:05 PM email (via DeepSeek), answering
Q45 + Q46 from email v13/v14. DeepSeek's recipe:

1. **Q45 — Canonical point:** "same Lagrangian parameters as electron — g=1.0,
   m_J=1.0 (after geometry normalization). Your scan gives m_J ≈ 0.46-0.56,
   likely factor ~2 short due to cylindrical vs spherical mismatch. Once η
   applied, m_J should scale to 1.0."
2. **Q46 — Geometry normalization:**
   `η = (∫ β² r dr) / (∫ β² r² dr)`
   "Multiply your neutral H/Q by η before applying the m_e scaling. For β
   peaking at r ~ 1.5, η ≈ 0.5-0.7."
3. **m_χ = η × (H/Q) × m_e**
4. **C from far-field K_1 amplitude using spherical measure.**

**Net result:** final DM paper inputs delivered at canonical (g=1.0, B0=0.5):
**m_χ = 0.460 MeV, m_J = 0.618 MeV, C = 770 MeV·fm.** Plus one empirical
finding for Paul/DeepSeek to interpret: m_J_corrected = m_J/η is family-
invariant at 1.21024 (NOT 1.0), suggesting a virial-type identity. The
"factor ~2" heuristic was approximate; actual factor 1/η = 2.35.

---

## Workflow — Paul/DeepSeek's recipe applied step-by-step

| Step | Action | Result |
| --- | --- | --- |
| 1 | Re-run v9 phase 2 BVP at canonical (g=1.0, B0=0.5) | m_J_raw = +0.5145, peak β = 0.697 @ r = 1.782, tail/peak = 1.1×10⁻⁵, sign-changes = 0, H = 9.017, Q_J = 4.259, (H/Q)_spherical = 2.117 |
| 2 | Compute η = (∫β²·r·dr) / (∫β²·r²·dr) for the ground state | η = 0.4251 (DeepSeek predicted 0.5-0.7; actual slightly lower due to β peaking at r = 1.78) |
| 3 | m_J_corrected = m_J_raw / η | m_J_corrected = 0.5145 / 0.4251 = **1.2102** (vs DeepSeek target 1.0; 21% deviation) |
| 4a | B0 scan to find canonical where m_J_corrected = 1.0 | **No bracket found in [0.10, 0.65]**; m_J_corrected is empirically INVARIANT at 1.21024 across the entire ground-state family |
| 4b | g-scan at B0=0.5 to find canonical g | **Same invariance** at all g ∈ [0.5, 1.6]; m_J_corrected stays at ±1.21024 (sign-flip indicates excited branch at higher g) |
| 4c | Lock canonical at (g=1.0, B0=0.5) | v9 phase 2 reference point — clean ground state with sign-changes = 0 |
| 5 | m_χ = η × (H/Q) × m_e = 0.4251 × 2.1173 × 0.511 MeV | **m_χ = 0.4599 MeV** |
| 6 | m_J in physical units: m_J (natural) × m_e = 1.2102 × 0.511 | **m_J = 0.6184 MeV** |
| 7 | Fit far-field tail to β(r) ~ C·exp(-m_J·r)·(1 + 1/(m_J·r))/r | **C = 770 MeV·fm** (fit residual 5.2% — caveat) |

---

## Family-invariance of m_J_corrected — empirical finding

Across the v9 phase 2 ground-state family at g = 1.0:

| B0 | m_J_raw | η | m_J_corrected | H/Q_eff | sign_changes |
| --- | --- | --- | --- | --- | --- |
| 0.35 | +0.43043 | 0.3557 | +1.21023 | 0.52714 | 0 |
| 0.40 | +0.46015 | 0.3802 | +1.21023 | 0.64404 | 0 |
| 0.45 | +0.48806 | 0.4033 | +1.21024 | 0.76848 | 0 |
| **0.50** | **+0.51446** | **0.4251** | **+1.21024** | **0.90005** | **0** |
| 0.55 | +0.53957 | 0.4458 | +1.21024 | 1.03837 | 0 |
| 0.60 | +0.56356 | 0.4657 | +1.21024 | 1.18313 | 0 |

And at B0 = 0.5 across g:

| g | m_J_raw | η | m_J_corrected | H/Q_eff | m_χ (MeV) | sign_changes |
| --- | --- | --- | --- | --- | --- | --- |
| 0.5 | +0.43261 | 0.3575 | +1.21024 | 0.5352 | 0.2735 | 0 |
| 0.7 | +0.47057 | 0.3888 | +1.21023 | 0.6888 | 0.3520 | 0 |
| **1.0** | **+0.51446** | **0.4251** | **+1.21024** | **0.9000** | **0.4599** | **0** |
| 1.2 | +0.53845 | 0.4449 | +1.21024 | 1.0319 | 0.5273 | 0 |
| 1.4 | +0.55961 | 0.4624 | +1.21023 | 1.1584 | 0.5920 | 0 |
| 1.6 | +0.57860 | 0.4781 | +1.21023 | 1.2805 | 0.6543 | 0 |

**m_J_corrected = m_J/η = 1.21024 to 6 decimal places, invariant under both
B0 and g.** Looks like a Pohozaev-type virial identity for the 3D spherical
l=1 cubic NLS — not a tunable parameter.

| Interpretation | Reading |
| --- | --- |
| DeepSeek heuristic | "m_J should scale to 1.0 after geometry correction" — was rough; actual factor 1/η ≈ 2.35 gives m_J_corr = 1.21024 |
| Virial reading | 1.21024 IS the canonical value; the "1.0 target" was an approximation |
| Alternative correction formula | m_J_corrected = m_J/η might not be exactly right; some other dimensionless combo could give 1.0 exactly |

**Reported to Paul/DeepSeek in email v15 for their interpretation.** The
deliverable numbers (m_χ, m_J, C) follow either reading; only the
identification of "the canonical match" interpretation depends on it.

---

## Final DM paper inputs (canonical point: g=1.0, B0=0.5)

| Quantity | Value | How computed |
| --- | --- | --- |
| **m_χ** | **0.4599 MeV** | η × (H/Q)_spherical × m_e = 0.4251 × 2.1173 × 0.511 |
| **m_J** | **0.6184 MeV** | m_J_corrected × m_e = 1.2102 × 0.511 (post-Q46 geometry) |
| **C** | **770 MeV·fm** | Fit far-field β(r) ~ C·exp(-m_J·r)·(1+1/(m_J·r))/r → C_natural = 3.903, scale by R_phys × m_e = 197.3 / 0.511 fm × 0.511 MeV |
| **η** | **0.4251** | (∫β²r dr) / (∫β²r² dr) at ground-state β profile |
| (H/Q)_spherical | 2.1173 | Raw H/Q in spherical-measure integration |
| H/Q_eff | 0.9000 | η × (H/Q)_spherical |
| Peak β @ r | 0.6970 @ r = 1.782 | (natural units; r_phys = 1.78 × R_phys ≈ 687 fm) |
| Tail/peak | 1.1×10⁻⁵ | Clean K_1 exponential decay |
| Sign changes | 0 | True ground state (no nodes) |

### Caveats

| Caveat | Detail |
| --- | --- |
| C fit residual | 5.2% max relative residual in far-field fit. Tightens with larger r_min_fit; could improve if needed. |
| m_J_corrected = 1.21 | Family-invariant; doesn't match DeepSeek's "1.0 target" heuristic exactly. Deliverables stand under either interpretation. |
| Geometry correction formula | DeepSeek's recipe `η × (H/Q)` produces m_χ = 0.46 MeV; if the correction is actually different (e.g., m_J_corr should be m_J/√η or m_J × η_other), then C and m_J shift but m_χ likely close. |

---

## Workflow summary — v9 phase 2 → v10 canonical

```text
v9 phase 1   →  BVP exhausted (Townes pathology + wrong helicity for Q_CS=0)
                Email v12 sent w/ Q43 (sign) + Q44 (geometry) → DeepSeek confirms both
v9 phase 2   →  Modified BVP with proper l=1 BC + B0 fixed + m_J free eigenvalue
                neutral_bvp_solver_mJ_free.py finds true nonlinear ground state
                Email v13 sent w/ Q45 (canonical point) + Q46 (normalization)
v14          →  v13 content resent + Acks-update ask bundled (after Paul's PM blip)
Q45/Q46 ←    →  Paul/DeepSeek 7:05 PM reply: g=1.0, m_J=1.0 after eta normalization;
                "factor ~2 from cyl vs spher mismatch"; eta = (∫β²r dr)/(∫β²r² dr)
v10          →  Apply Q46 recipe: compute eta, m_J_corrected, m_chi = eta·(H/Q)·m_e,
                extract C from K_1 far-field amplitude
                FINAL NUMBERS: m_chi = 0.46 MeV, m_J = 0.62 MeV, C = 770 MeV·fm
                PLUS empirical finding: m_J_corr = 1.21 is family-invariant (virial?)
                Email v15 sent w/ numbers + family-invariance observation +
                Acks-update reinforcement (item missed in DeepSeek's v14 reply)
```

---

## Cross-references

- `0b_M6_roadmap.md` — sandbox sequence + v10 update
- `0b_model_gates.md` — G2 status: GROUND-STATE-FOUND-PENDING-CALIBRATION-POINT → DM PAPER INPUTS DELIVERED (pending Acks-update + virial-identity interpretation)
- `0b_question_tracker.md` — Q45/Q46 RESOLVED via DeepSeek 7:05 PM reply; Q47 NEW (virial identity / canonical interpretation)
- `0c_sandbox_v9.md` — v9 phase 2 ground state that feeds v10 numerical results
- `sandbox_v10/m6_v10_canonical_neutral_chaoiton.py` — the working script
- `sandbox_v10/v10_output.txt` — full script output
