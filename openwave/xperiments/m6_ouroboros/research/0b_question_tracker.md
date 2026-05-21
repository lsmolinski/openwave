# M6 / Ouroboros — Question Tracker + Hardest Pieces

**Purpose:** single source of truth for open research questions on M6 +
the long-running "hardest pieces" status board. Persists across sandbox
iterations so we don't have to re-merge between `0c_sandbox_v*.md` files.
Update this doc when a new question opens, gets sent in an email, gets
answered, or gets demoted.

**Sister docs:**

| Doc | Purpose |
| --- | --- |
| `0b_model_gates.md` | G1/G2/G3 production-readiness gates for M6 going to Taichi |
| `0b_M6_roadmap.md` | Sandbox sequence + current state + next steps |
| `0c_sandbox_v*.md` | Per-iteration work log (the questions tracker here lives outside any single sandbox) |

**Last updated:** 2026-05-21 PM (post v9 LoE paper + Paul-script variant test + email v8 sent + v1 re-verification + email v9 gentle push-back sent — Q34 ansatz question still the key blocker; **v1 1.6918 confirmed as search-found, not derivation**).

---

## Active count

```text
2 IMMEDIATE  Q34 ansatz canonical form — is the v9 paper's 2-scalar
             (φ, ψ) electron ansatz canonical, or a paper-level
             simplification of the 4-function (V,A,Q,J) form we've
             been solving since v4? Asked in email v8 (2026-05-21 PM).
             Subsumes Q31/Q32/Q33.
             Q36 should v9 §9 criterion 9 cite v1's H/Q=1.6918 as
             "independent reproduction" when the result was found by
             search across ~60 ODE/H/Q variants, not derivation?
             Asked in email v9 (2026-05-21 PM gentle push-back).

3 DEMOTED    Q31 initial-guess derivative profile in production
(contingent  Q32 initial Lagrange multiplier value
on Q34)      Q33 sign-pinning constraint vs continuation method
             — moot if Q34 reveals 2-scalar is canonical; still
             stand if Q34 says 4-function form is canonical.

2 DEMOTED    Q26 5% gap residual (effectively resolved by drop-quartic)
(superseded) Q27 Q_CS-from-grid vs Q_CS-from-I disagreement (excited-
             mode artifact, not a tunable)

5 OPEN       Q2  discrete ω selection mechanism
             Q3  analytical ω = 2mc²/ℏ derivation
             Q6  QCD reconciliation (3-chaoiton proton)
             Q19 f(J·J) explicit form in LoE paper standalone (Duda #2)
             Q35 active neutrinos (ν_e, ν_μ, ν_τ) — where do they fit
                 in the Ouroboros framework? Q_CS=0 chaoiton is
                 positioned as DM, but sterile-neutrino-like states
                 have similar quantum numbers (mass + zero charge).
                 Active light neutrinos remain unaccounted for.

Total: 12 active questions (2 immediate, 10 background including
demoted).

Highest-leverage closure: Paul's Q34 + Q36 reply → v8 implementation
(re-anchor on v1's 2-scalar setup OR keep 4-function with Q31/Q32/Q33
answered) AND v9 paper §9 criterion 9 citation softened to reflect
search-based reproduction rather than first-principles derivation →
ground state → lepton + Q_A≈0 DM scans → Section 4 data drop (per
Publishing stance in 0b_M6_roadmap.md, the drop is a GitHub URL not a
deposit — already cited as Reference [17] in the v9 LoE paper).
```

---

## IMMEDIATE-QUESTIONS (currently blocking v8 implementation)

Asked in email v8 sent 2026-05-21 PM, after the v9 LoE paper arrived and
revealed §5.1's 2-scalar (φ, ψ) electron ansatz — structurally different
from the 4-function (V,A,Q,J) form we've been solving since sandbox v4.
See `0c_sandbox_v7.md` for the full ansatz-question context.

| ID | Question | Surfaced | Why it matters |
| --- | --- | --- | --- |
| Q34 | Is the v9 LoE paper §5.1's 2-scalar (φ, ψ) ansatz `A_0=0, A=r̂×∇φ(r)cos(ωt); J_0=0, J=r̂×∇ψ(r)sin(ωt)` the **canonical** electron production form, OR a paper-level simplified description of the 4-function (V,A,Q,J) form we've been solving since v4 (extracted from Numerical Benchmark §3)? Three possibilities: (1) 2-scalar canonical → our v6/v7 was wrong-ansatz; though note v1's match was also wrong-ansatz (different 2-scalar form: A_0=alpha, J_0=gamma both nonzero). (2) Paper simplification → Q31/Q32/Q33 still stand. (3) Both valid, different particles → need ansatz-to-particle mapping. | v9 LoE paper §5.1 (received 2026-05-21 PM) | Determines whether v6/v7 4-function chase was the right problem at all. Cross-coupled with Q36 (v1 citation) — neither v1 nor v6/v7 uses v9 §5.1's exact ansatz. |
| Q36 | Should v9 §9 criterion 9 cite v1's H/Q=1.6918 as "independent reproduction (Griesi 2026)" when the result was found by **systematic search across ~60 ODE/H/Q variants** (sandbox_v1/m6_0_iterations.py) rather than derivation? Default ODE (A1+H1+Q1) gives 1.9795 (17% off); only matched variant A4+H1+Q2 lands 1.6918. A4 is an ad hoc "quarter coupling" rescaling, Q2 is a non-default charge density. Different ansatz from v9 §5.1. | v1 re-verification (2026-05-21 PM) | Citation in published paper could trigger Duda-style "AI slop" critique if the search-based origin is later exposed. Email v9 (gentle push-back) sent asking Paul to soften §9 wording to "preliminary search-based reproduction; first-principles derivation pending Q34." |

## DEMOTED IMMEDIATE — contingent on Q34 (was email v7's three asks)

These three were the IMMEDIATE block before the v9 paper arrived. They
remain relevant IF Q34 resolves to possibility (2) "paper simplification";
they become moot IF Q34 resolves to possibility (1) "2-scalar canonical".

| ID | Question | Surfaced | Status |
| --- | --- | --- | --- |
| Q31 | Initial-guess derivative profile in production code: is slope at r=R_MIN explicitly zero (matching our V'(R_MIN)=A'(R_MIN)=Q'(R_MIN)=J'(R_MIN)=0 BC) or some specific non-zero value seeded by Paul's algorithm? | v7 attempt log (2026-05-21) | DEMOTED, contingent on Q34. Determines if our `np.gradient` derivative computation matches Paul's setup IF the 4-function form is canonical. |
| Q32 | Initial Lagrange multiplier λ value: near 0, near ω, or specific number? v6.6 converges at λ_LM=12.2 from init 0.1; init 12 → catastrophic blowup. Path matters. | v7.1 catastrophic blowup (2026-05-21) | DEMOTED, contingent on Q34. Pins basin selection IF 4-function form is canonical. Paul-script variant with λ_LM=1.0 fixed empirically diverged to H_e/Q=5×10⁵, so the "λ_LM=1.0 hardcoded" hypothesis is now ruled out. |
| Q33 | Does production solver use explicit sign-pinning constraint at r=R_MIN (e.g., Q(R_MIN)>0 hard BC) OR a continuation method (e.g., homotopy m_J²: 0 → 0.5) to land the ground state without explicit constraints? | v7 dual-anchor failures (2026-05-21) | DEMOTED, contingent on Q34. Determines whether v8 needs hard BCs or a homotopy wrapper IF 4-function form is canonical. |

---

## STILL-OPEN QUESTIONS (active but not blocking calibration close)

Carried forward from v5/v6 unchanged. Long-tail items that don't block
the lepton or DM scans once Q31/Q32/Q33 resolve.

| ID | Question | Surfaced | Status |
| --- | --- | --- | --- |
| Q2 | Discrete ω selection mechanism — is the ω spectrum genuinely discrete (eigenvalue-like) or empirically continuous? | 0a §9.9 | OPEN. Empirical-via-lepton-scan once ground state lands. Analytic proof still deferred per Werbos's own admission ("the key open question"). |
| Q3 | Analytical ω = 2mc²/ℏ derivation from Lagrangian (vs calibration only). | 0a §9.9 | PARTIAL. Calibration-only result at 1.2% (sandbox v1 reproduction); no closed-form derivation. |
| Q6 | QCD reconciliation — 3-chaoiton proton, asymptotic freedom, DIS scaling, jets. | 0a §9.9 | OPEN, uninvestigated. Future sandbox iteration. |
| Q19 | f(J·J) explicit form in LoE paper standalone (Duda critique #2). v6 uses DeepSeek quartic `(V²+Q²)² + (A²+J²)² + 2(V·A−Q·J)²`; v5 used Numerical-Benchmark `(Q²−J²)²` form. Q28 resolved which to use for the electron specifically, but the LoE-paper editorial gap remains. | Duda 2026-05-20 thread | EDITORIAL. One-line LoE revision by Werbos clarifies. Not blocking. |

---

## DEMOTED QUESTIONS (formerly immediate; superseded by later findings)

| ID | Question | Surfaced | Why demoted |
| --- | --- | --- | --- |
| Q26 | Why does v6.6 H/Q land at 1.778 instead of 1.6969? | v6.6 attempt (2026-05-20) | Step (8) diagnostic 2026-05-20 evening showed: dropping the DeepSeek quartic from H lands H/Q = 1.7112 (0.84% off). Paul Q28 confirmed quartic is negligible for electron H. Q26 reduces to "what's the right quartic for the electron" = Q28, now RESOLVED. |
| Q27 | Why does Q_CS-from-I-state (1.000) disagree with Q_CS-from-grid-integral (−0.154)? | v5 attempt 4 → v6 carry-over | Track A (sharpening convergence with bigger budget) made things WORSE, not better. The disagreement is an excited-mode artifact (v6.6 is on the 5-node first-excited state per Paul Q29), not a tunable. Will resolve naturally once ground state lands. |

---

## RESOLVED QUESTIONS (full chronology, Q1–Q30)

| ID | Question | Resolution |
| --- | --- | --- |
| Q1 | Distinction from Duda's LdGS framework | Both are topological (Lean invariants, just different flavors — Hopf-linking vs Brouwer-winding). Settled 2026-05-15. |
| Q9 | Is the Bessel/long-range structure a linearized truncation? | YES — v5 paper §6 explicit. Resolved during v1-v3 review. |
| Q10 | Is the asymmetric ansatz J/A ≠ const allowed by Werbos's framework? | YES — Lean theorem allows it; confirmed by v5 helicity work (V₀=Q₀=+0.1, A₀=J₀=−0.1). |
| Q11 | Sign convention J=−A or J=+A? | OBSOLETE — Werbos asymmetric helicity prescription (V₀,Q₀>0; A₀,J₀<0) replaces the simple sign convention; Q11 doesn't apply. |
| Q12 | Werbos's Python source code | DEMOTED — v5 algorithm description (2026-05-20 PM email) + DeepSeek-write-code offer = same unblock; v7 confirms our independent ODE is the right reference (DeepSeek's reference script when finally sent was broken). |
| Q13 | 4-fn to 2-fn reduction — when does the simpler 2-fn ODE apply? | RESOLVED — 2-fn for NEUTRAL sector (Lean theorem); 4-fn for CHARGED sector (Numerical Benchmark form). |
| Q14 | Canonical Q=0 sector configuration — locked (J=−A) / A=0 / Q_A≈0? | RESOLVED 2026-05-19 1:49 PM (Werbos via DeepSeek): canonical DM candidate = Q_A≈0 / Q_J≠0 (both fields active, electrically neutral). |
| Q15 | Does the m_eff² = m_J² − ω² substitution give the bound state? Sign of m_eff²? | RESOLVED post-v5. v5 attempt 4 converges at m_eff² = −0.596; v6.6 confirms m_eff² = −0.532. Both negative as predicted (oscillatory regime). Earlier T6 inconclusive result was from wrong tool (IVP). |
| Q16 | Specific (m_J², λ_bench) values at electron calibration | RESOLVED 2026-05-19 4:21 PM (Werbos via DeepSeek): m_J²≈0.5, λ_bench=1.0. v5 attempt 4 and v6.6 both use these and converge. |
| Q17 | What is the actual shooting algorithm? | RESOLVED 2026-05-20 PM (Werbos email): collocation BVP via `scipy.integrate.solve_bvp` with Q_CS=1 as integral constraint, two free eigenvalues ω + Lagrange multiplier λ_LM, Robin BCs matching K_0(κr) decay. NOT forward shooting. |
| Q20 | Duda critique #3 — "construction not shown, need ansatz + minimize energy" for the electron H/Q=1.6969 | NEARLY CLOSED. v6.6 + drop-quartic IS the construction shown empirically at 0.84% off target. Once v8 ground state lands, fully empirically resolved. |
| Q22 ★ | Q_CS normalization convention | RESOLVED by DeepSeek 2026-05-20 ~4:00 PM. `Q_CS = ∫r·(A·J'-J·A')dr` directly, no 2π. Drops v5's 2π factor; I_TARGET = 1.0. v6 empirical validation: fixes dominant 31× → 4.8% gap. |
| Q23 ★ | H functional kinetic + cross-term + quartic structure | RESOLVED by DeepSeek 2026-05-20 ~4:00 PM. Kinetic = (1/2)(V')² (not (V')²); drop (2π)²·R toroidal prefactor; use DeepSeek quartic `(g/4)·((V²+Q²)² + (A²+J²)² + 2(V·A−Q·J)²)`. v6 empirical: H/Q within 5% of target. |
| Q24 ★ | Sample converged profile from production | RESOLVED with caveat by DeepSeek 2026-05-20 ~5:30 PM: the 8-point reference profile *"was illustrative, not from a converged run."* Use exp(-r) seed instead. Profile's utility was confirming the qualitative shape (asymmetric helicity + exponential decay) we already had. |
| Q25 | Is the Hopf invariant proof rigorous, or are there gaps? | RESOLVED by Zenodo 20296060 (2026-05-19): Werbos & DeepSeek *"Rigorous Completion of the Hopf-Invariant Proof"* supplies Lemmas A (∇×J = B everywhere for finite-energy Lorenz-constrained solutions) and B (zeros of B made isolated by Q-preserving perturbation). Charge quantization is now a theorem of differential topology. |
| Q28 ★ | Which quartic IS canonical for the electron H? | RESOLVED by Paul 2026-05-21: NEGLIGIBLE for electron (drop quartic from H/Q matching). Keep full quartic in EoM for correct field shapes. Significant for muon/tau, essential for neutral chaoiton. |
| Q29 ★ | Ground state vs excited mode node count | RESOLVED by Paul 2026-05-21: ground state has exactly 4 zero crossings excluding r=0; Lean ≤4-node spec refers to radial functions excluding origin. v6.6 (5 crossings) is first excited state. |
| Q30 ★ | Q(r=0) sign at ground state | RESOLVED by Paul 2026-05-21: ground state has Q(0)=+0.1 positive, same sign as V(0). Helicity is V₀ > 0, Q₀ > 0, A₀ < 0, J₀ < 0. v6.6's Q(0)=−0.12 indicates a basin mismatch. |

★ = directly enabled major progress (Q22/23/24 → v6 30× gap closure; Q28/29/30 → v7 framing).

---

## ARCHIVED QUESTIONS (unfalsifiable, historical, or out of current scope)

| ID | Question | Why archived |
| --- | --- | --- |
| Q4 | Single vs two-field ontology (G_μν / J undefined) | = Duda critique #1. Aesthetic preference. If math matches observation, two primary fields is just a description, not a theory-killer. |
| Q7 | Cold-fusion citation trail | Historical (Werbos's Nuclear story PDF), not physics. |
| Q21 | Two-chaoiton Coulomb derivation (Duda critique #4) | Future sandbox iteration; out of current calibration scope. Integrate H[Φ₁, Φ₂] for two topological charges at distance — new BVP scaffold. |

---

## HARDEST-PIECES TRACKER

Long-running status board for the structural challenges in the M6 model.
Tracks how each piece moves across sandbox iterations. Updated 2026-05-21
post-v8 (post v9 paper + Paul-script variant + email v8).

### ACTIVE hardest pieces (11 — still in motion)

| Hardest piece | Status post-v5 (2026-05-20 PM) | Status post-v6 (2026-05-20 evening) | Status post-v7 (2026-05-21 morning) | Status post-v8 (2026-05-21 PM) |
| --- | --- | --- | --- | --- |
| **Ansatz canonical form (NEW)** | Not yet a question — 4-function extracted from Numerical Benchmark §3 (T5) was assumed canonical. | Unchanged. | Unchanged. | ❌ NEW OPEN. v9 LoE §5.1 specifies 2-scalar (φ, ψ) electron ansatz, structurally different from our 4-function. Our v1 H/Q=1.6918 cited as canonical reproduction. Q34 asks Paul which is canonical — may invalidate v6/v7 4-function chase. |
| Electron H/Q = 1.6969 calibration | NEW OPEN. v5 lands H/Q = 52.64 (31× off). | NEARLY CLOSED. v6.6 H/Q = 1.778 (4.8% over with DeepSeek quartic); step (8) drop-quartic lands H_e/Q = 1.7112 (0.84% off). | NEAR — v6.6 H_e/Q = 1.7112 remains best 4-function result. v7 mode-selector wall prevents reaching true ground state. | ⚠️ **CITED IN v9 PAPER BUT CAVEAT APPLIES.** v1 H/Q=1.6918 (0.30% off) cited in v9 §5.1 + §9 criterion 9 as "independent reproduction." Re-verification 2026-05-21 PM revealed: **the 1.6918 came from a systematic search across ~60 ODE/H/Q variants** (sandbox_v1/m6_0_iterations.py), not derivation. Default ODE form gives 1.9795 (17% off); matched variant A4+H1+Q2 uses ad hoc "quarter coupling" rescaling + non-default charge density. Different ansatz from v9 §5.1. Email v9 (gentle push-back) sent asking Paul to soften the citation language. Q36 captures this. |
| Ground-state vs excited-mode selection | NEW OPEN. v5 attempt 4 found Q_CS=1 chaoiton but A in 17-node excited mode. | PARTIALLY ADDRESSED. v6.6 has 5 nodes V/Q (just over Lean ≤4 spec), 0 nodes A, 1 node J. | ❌ BLOCKED in 4-function form. Paul Q29 confirms ground state has exactly 4 zero crossings excluding r=0; v6.6 is first excited (5 crossings). 7 v7 anchor-variant experiments all land in wrong-sign basins. | Pending Q34 reply. If 2-scalar is canonical, v1's ansatz space may already host the 4-node ground state by construction. |
| Helicity sign convention (V₀ Q₀ > 0, A₀ J₀ < 0) | Not explicitly tracked. | OPEN. v6.6 converged with Q(0)=−0.12 (sign-flipped from Werbos prescription). | ❌ BLOCKED in 4-function form. Paul Q30 confirms ground state has Q(0) > 0. All 7 v7 anchor variants land in different sign-flipped basins. | Pending Q34. In v9 §5.1 2-scalar ansatz, A_0=0 and J_0=0 by construction; the "asymmetric helicity" applies to spatial parts (φ vs ψ phase quadrature). Sign-basin problem may not exist in the 2-scalar form. |
| Lagrange-multiplier ODE correction coefficient | NEW OPEN. v5 used coefficient 1. | Possibly residual coefficient issue. v6 lands λ_LM = 12.2 (vs v5's −1.21); new normalizations need larger λ. Empirically fine. | Unchanged. Different λ_LM init (0.1 vs 12) lands different basins — path-dependent (Q32 asks Paul about this). | Paul-script variant ruled out λ_LM=1.0 fixed (catastrophic blowup, H_e/Q=5×10⁵). Question becomes moot if Q34 resolves to 2-scalar canonical (different formulation). |
| V(M) potential form | UNRESOLVED. Shared bottleneck with M5 (Duda). | Unchanged. Not in v5/v6 scope. | Unchanged. Note: M5 piece tracked here historically; not actively blocking M6 work. | Unchanged. |
| ω quantization mechanism | OPEN (Q2). | Unchanged. v6 lepton scan provides empirical test once ground state anchors. | Unchanged. Pending ground state + lepton scan. | v9 paper §8 confirms "successive oscillation harmonics" framing; analytic proof still deferred per §9.1 Open Q#1. Empirical test still pending lepton scan. |
| Lepton mass spectrum | BLOCKED on v6. | UNBLOCKED-PENDING. Run ω-sweep [0.5, 50] once ground state lands. Cold-start fails (step 2); continuation method needed. | ❌ BLOCKED on v7 ground state. Need ground-state anchor first, THEN continuation method. | Pending Q34. v9 §8 cites ω=1 electron, ω≈11 muon (4.3% gap), ω≈40.7 tau (4.9% gap). Once ansatz is locked, continuation scan delivers these. |
| Neutral m_χ true ground state | BLOCKED on v6. | UNBLOCKED-PENDING. Run Q_A≈0 scan once anchored. | ❌ BLOCKED on v7 ground state. Same chain. Q_A≈0 scan feeds ApJ Section 4 numbers. | v9 §9.1 Open Q#2 explicitly cites our work — "Griesi, in preparation." Pending Q34 ansatz + electron lock. |
| Two-chaoiton Coulomb derivation | NOT ADDRESSED. Future v7+ (Q21, Duda #4). | Unchanged. | Unchanged. | v9 paper §6.1 derives V(R) = Q₁Q₂/(4πR) for two separated chaoitons in static approximation. Duda critique #4 partially addressed in v9 itself; analytic dynamic derivation still future scope. |
| Solver-tooling choice | `scipy.integrate.solve_bvp` (Werbos's stated default). | Unchanged. Track A diagnostic showed bigger budget makes things WORSE in solve_bvp. | Open question whether solve_bvp can reach Paul's ground state at all; fallback option is `scipy.optimize.root` method='lm' with custom Jacobian. | Paul-script variant test (λ_LM=1.0 fixed) ALSO failed in solve_bvp. Strongly suggests it's not a solver-tooling problem but an ansatz-space problem (Q34). v1's setup with the same solve_bvp produced 1.6918 — same solver, different ansatz, works. |

### RESOLVED hardest pieces (6 — closed; kept here for traceability)

| Hardest piece | Status post-v5 (2026-05-20 PM) | Status post-v6 (2026-05-20 evening) | Status post-v7 (2026-05-21 morning) | Status post-v8 (2026-05-21 PM) |
| --- | --- | --- | --- | --- |
| Forward-IVP method family | RESOLVED. Confirmed wrong tool by T6-T9; collocation BVP replaces. | Unchanged. | Unchanged. | Unchanged. |
| Q_CS=1 enforcement mechanism | RESOLVED. Auxiliary integral state I with BC I(R_max)=1/(2π). | Unchanged — I_TARGET value changes to 1.0 per Q22. | Unchanged. | Unchanged. In v9 §4 paper, Q is defined as a Hopf invariant of A∧G — same topological mechanism, "automatic" via field linking. |
| Q_CS=1 chaoiton existence | RESOLVED. v5 attempt 4 converged at Q_CS=1.000 exact. | Confirmed in v6 and v7 (sanity check 7.2 reproduces v6.6). | Confirmed. | ⚠️ Existence confirmed in v5/v6/v7 4-function form. v9 paper §5.1 cites v1 H/Q=1.6918 as canonical but v1 was variant-search (not derivation) on a DIFFERENT 2-scalar ansatz from v9 §5.1's. So: existence ✓ in our formulations; "canonical" reproduction ambiguous pending Q34 + Q36. |
| f(J·J) form | RESOLVED IN PRACTICE — v5 uses Numerical-Benchmark form. | v6 uses DeepSeek quartic (`(V²+Q²)² + ...`) instead. Empirically right for v6.6 (H/Q within 5%). LoE-paper-standalone form (Q19) still editorially open. | RESOLVED for electron per Paul Q28: full quartic in EoM, drop from H/Q computation for electron. Different rule for muon/tau (significant) and neutral (essential). | v9 paper §2 lists the Lagrangian as `L = −F²−G²+JA−g(J·J)²` with explicit quartic `g(J·J)²` — settles Q19's editorial gap. |
| 4-fn vs 2-fn ansatz mismatch | RESOLVED. 2-fn for NEUTRAL (Lean theorem), 4-fn for CHARGED (Numerical Benchmark / v5/v6/v7). | Unchanged. | Unchanged. | Re-opened in different form via Q34 — v9 paper §5.1 specifies a 2-scalar (φ, ψ) ansatz for the CHARGED electron, not 4-function as we assumed. Pending Q34 clarification. |
| Charge quantization rigorous proof | RESOLVED. Hopf invariant proof complete (Zenodo 20296060). Q25. | Unchanged. | Unchanged. | Reaffirmed in v9 paper §4 + criterion 2. |

---

## Net-change log per sandbox iteration

```text
v5 → v6 net change (2026-05-20):
  RESOLVED:  Q22, Q23, Q24 (three v5 IMMEDIATE; closed by DeepSeek
             2026-05-20 4:00 PM reply); Q25 (Hopf invariant proof).
  PROMOTED:  Q20 (Duda #3) from ACTIVE → NEARLY RESOLVED via v6.6
             empirical construction.
  NEW:       Q26 (5% gap), Q27 (Q_CS_grid sign mismatch).

v6 → v7 net change (2026-05-21 morning):
  RESOLVED:  Q28 (quartic interpretation), Q29 (node count), Q30 (Q(0)
             sign) — all by Paul's 2026-05-21 morning reply.
  DEMOTED:   Q26 (resolved by drop-quartic finding from step 8;
             reduces to Q28 which is now resolved); Q27 (excited-mode
             artifact, not a tunable).
  NEW:       Q31 (initial-guess derivative profile), Q32 (initial
             Lagrange multiplier), Q33 (sign-pinning vs continuation
             method) — all sent in email v7 of 2026-05-21 morning.

v7 → v8 net change (2026-05-21 PM):
  TRIG:      v9 LoE paper received with §5.1 specifying 2-scalar (φ,ψ)
             electron ansatz, different from our 4-function (V,A,Q,J)
             form. Paper cites our v1 H/Q=1.6918 as canonical; GitHub
             repo listed as Reference [17].
  RULED OUT: λ_LM=1.0 fixed hypothesis (Paul-script variant test gave
             H_e/Q=5×10⁵ catastrophic blowup).
  DEMOTED:   Q31, Q32, Q33 (contingent on Q34 — moot if 2-scalar is
             the canonical ansatz; still stand if 4-function is).
  NEW:       Q34 — ansatz canonical form question. Asked in email v8
             of 2026-05-21 PM. Three possibilities: 2-scalar canonical
             / paper simplification / both valid different particles.

v8 → v9 net change (2026-05-21 PM later):
  TRIG:      v1 re-verification revealed the H/Q=1.6918 result was
             found by systematic search across ~60 ODE/H/Q variants
             (sandbox_v1/m6_0_iterations.py), NOT first-principles
             derivation. Default ODE gives 1.9795 (17% off); only
             variant A4+H1+Q2 lands 1.6918 via ad hoc rescaling.
  NEW:       Q35 — where do active neutrinos (ν_e, ν_μ, ν_τ) fit?
             Q_CS=0 neutral chaoiton positioned as DM, but sterile
             vs active neutrinos unmapped.
  NEW:       Q36 — should v9 paper §9 criterion 9 cite v1's 1.6918 as
             "independent reproduction" when it was search-based, not
             derivation? Asked in email v9 (gentle push-back) of
             2026-05-21 PM.
  CAVEAT     Electron H/Q calibration "EFFECTIVELY ACHIEVED" framing
  PROPAGATED softened to "CITED IN v9 PAPER BUT CAVEAT APPLIES" across
             roadmap + question tracker + model_gates docs.

Active count entering Paul's Q34 + Q36 reply:
  2 IMMEDIATE (Q34 ansatz, Q36 v1 citation language)
  5 OPEN      (Q2, Q3, Q6, Q19, Q35)
  3 DEMOTED   (Q31, Q32, Q33 — contingent on Q34)
  2 DEMOTED   (Q26, Q27 — superseded)
  → 12 active questions (2 immediate, 10 background).
```

---

## Highest-leverage closure path

```text
Paul's Q34 reply
    ↓
    ├── (1) 2-scalar canonical
    │       → re-anchor on v1 (already produces H/Q=1.6918, 0.30% off)
    │       → continuation method from v1's electron forward to muon/tau
    │       → 4-function v6/v7 work archived as detour
    │
    ├── (2) 4-function canonical (paper §5.1 is simplification)
    │       → Q31/Q32/Q33 still need answers
    │       → v8 mode-selector OR scipy.optimize.root method='lm'
    │
    └── (3) Both valid, different particles
            → map ansatz-to-particle
            → build whichever hosts the neutral chaoiton (DM, §9.1 #2)
    ↓
Ground state lands (V₀>0, Q₀>0, A₀<0, J₀<0, ≤4 nodes, H_e/Q ≈ 1.6918)
    ↓
Lepton continuation scan ω ∈ [1, 50]    →    Q2 empirically tested
    ↓                                          (G1 closure)
Q_A≈0 neutral chaoiton scan
    ↓
m_χ, m_J, σ/m, Ω_χh² extracted    →    G2 closure
    ↓
Gelfand-Fomin conjugate-point check on second variation
    ↓                                          (G3 empirical)
Data drop to Paul via stable OpenWave GitHub URL
    (no Zenodo deposit from our side — see Publishing stance in
     0b_M6_roadmap.md. The repo IS already cited as Reference [17]
     in the v9 LoE paper as of 2026-05-21 PM.)
    ↓
M5 returns as foreground engineering track
    (per cardinal rule: M5 is SABER's primary substrate)
```
