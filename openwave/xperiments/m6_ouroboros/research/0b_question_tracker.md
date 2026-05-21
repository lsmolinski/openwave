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

**Last updated:** 2026-05-21 morning (post v7 Phase 1 + email v7 sent).

---

## Active count

```text
3 IMMEDIATE  Q31 initial-guess derivative profile in production
             Q32 initial Lagrange multiplier value
             Q33 sign-pinning constraint vs continuation method
             — all asked in email v7 (2026-05-21 morning); block v7
             ground-state convergence.

2 DEMOTED    Q26 5% gap residual (effectively resolved by drop-quartic)
             Q27 Q_CS-from-grid vs Q_CS-from-I disagreement (excited-
             mode artifact, not a tunable)

4 OPEN       Q2 discrete ω selection mechanism
             Q3 analytical ω = 2mc²/ℏ derivation
             Q6 QCD reconciliation (3-chaoiton proton)
             Q19 f(J·J) explicit form in LoE paper standalone (Duda #2)
             — long-tail; not blocking lepton/DM scans

Total: 9 active questions (3 immediate, 6 background including demoted).

Highest-leverage closure: Paul's Q31/Q32/Q33 reply → v8 mode-selector
(or continuation method) → ground state → lepton + Q_A≈0 DM scans →
Section 4 data drop (per Publishing stance in 0b_M6_roadmap.md, the
drop is a GitHub URL not a deposit).
```

---

## IMMEDIATE-QUESTIONS (currently blocking v7 ground-state convergence)

All three asked in email v7 sent 2026-05-21 morning, after v7's 7-variant
mode-selector wall (see `0c_sandbox_v7.md` Attempt log).

| ID | Question | Surfaced | Why it matters |
| --- | --- | --- | --- |
| Q31 | Initial-guess derivative profile in production code: is slope at r=R_MIN explicitly zero (matching our V'(R_MIN)=A'(R_MIN)=Q'(R_MIN)=J'(R_MIN)=0 BC) or some specific non-zero value seeded by Paul's algorithm? | v7 attempt log (2026-05-21) | Determines if our `np.gradient` derivative computation matches Paul's setup. If non-zero, biases basin selection. |
| Q32 | Initial Lagrange multiplier λ value: near 0, near ω, or specific number? v6.6 converges at λ_LM=12.2 from init 0.1; init 12 → catastrophic blowup. Path matters. | v7.1 catastrophic blowup (2026-05-21) | Pins basin selection. Different inits land different basins. |
| Q33 | Does production solver use explicit sign-pinning constraint at r=R_MIN (e.g., Q(R_MIN)>0 hard BC) OR a continuation method (e.g., homotopy m_J²: 0 → 0.5) to land the ground state without explicit constraints? | v7 dual-anchor failures (2026-05-21) | Determines whether v8 needs hard BCs or a homotopy wrapper. |

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
post-v7.

### ACTIVE hardest pieces (10 — still in motion)

| Hardest piece | Status post-v5 (2026-05-20 PM) | Status post-v6 (2026-05-20 evening) | Status post-v7 (2026-05-21 morning) |
| --- | --- | --- | --- |
| Electron H/Q = 1.6969 calibration | NEW OPEN. v5 lands H/Q = 52.64 (31× off). | NEARLY CLOSED. v6.6 H/Q = 1.778 (4.8% over with DeepSeek quartic); step (8) drop-quartic lands H_e/Q = 1.7112 (0.84% off). | NEAR — v6.6 H_e/Q = 1.7112 remains best result. v7 mode-selector wall prevents reaching true ground state. Paul says ground state should give H/Q ≤ 1.6969 (lower than current). |
| Ground-state vs excited-mode selection | NEW OPEN. v5 attempt 4 found Q_CS=1 chaoiton but A in 17-node excited mode. | PARTIALLY ADDRESSED. v6.6 has 5 nodes V/Q (just over Lean ≤4 spec), 0 nodes A, 1 node J — much closer to ground state. | ❌ BLOCKED. Paul Q29 confirms ground state has exactly 4 zero crossings excluding r=0; v6.6 is first excited (5 crossings). 7 v7 anchor-variant experiments all land in wrong-sign basins. Awaiting Paul Q31/Q32/Q33 reply for solver-detail unblock. |
| Helicity sign convention (V₀ Q₀ > 0, A₀ J₀ < 0) | Not explicitly tracked. | OPEN. v6.6 converged with Q(0)=−0.12 (sign-flipped from Werbos prescription). | ❌ BLOCKED. Paul Q30 confirms ground state has Q(0) > 0. All 7 v7 anchor variants land in different sign-flipped basins; each anchor fixes one sign but flips another. |
| Lagrange-multiplier ODE correction coefficient | NEW OPEN. v5 used coefficient 1. | Possibly residual coefficient issue. v6 lands λ_LM = 12.2 (vs v5's −1.21); new normalizations need larger λ. Empirically fine. | Unchanged. Different λ_LM init (0.1 vs 12) lands different basins — path-dependent (Q32 asks Paul about this). |
| V(M) potential form | UNRESOLVED. Shared bottleneck with M5 (Duda). | Unchanged. Not in v5/v6 scope. | Unchanged. Note: M5 piece tracked here historically; not actively blocking M6 work. |
| ω quantization mechanism | OPEN (Q2). | Unchanged. v6 lepton scan provides empirical test once ground state anchors. | Unchanged. Pending ground state + lepton scan. |
| Lepton mass spectrum | BLOCKED on v6. | UNBLOCKED-PENDING. Run ω-sweep [0.5, 50] once ground state lands. Cold-start fails (step 2); continuation method needed. | ❌ BLOCKED on v7 ground state. Need ground-state anchor first, THEN continuation method to nudge ω upward. |
| Neutral m_χ true ground state | BLOCKED on v6. | UNBLOCKED-PENDING. Run Q_A≈0 scan once anchored. | ❌ BLOCKED on v7 ground state. Same chain. Q_A≈0 scan feeds ApJ Section 4 numbers. |
| Two-chaoiton Coulomb derivation | NOT ADDRESSED. Future v7+ (Q21, Duda #4). | Unchanged. | Unchanged. |
| Solver-tooling choice | `scipy.integrate.solve_bvp` (Werbos's stated default). | Unchanged. Track A diagnostic showed bigger budget makes things WORSE in solve_bvp. | Open question whether solve_bvp can reach Paul's ground state at all; fallback option is `scipy.optimize.root` method='lm' with custom Jacobian (Paul mentioned as alternative). |

### RESOLVED hardest pieces (6 — closed; kept here for traceability)

| Hardest piece | Status post-v5 (2026-05-20 PM) | Status post-v6 (2026-05-20 evening) | Status post-v7 (2026-05-21 morning) |
| --- | --- | --- | --- |
| Forward-IVP method family | RESOLVED. Confirmed wrong tool by T6-T9; collocation BVP replaces. | Unchanged. | Unchanged. |
| Q_CS=1 enforcement mechanism | RESOLVED. Auxiliary integral state I with BC I(R_max)=1/(2π). | Unchanged — I_TARGET value changes to 1.0 per Q22. | Unchanged. |
| Q_CS=1 chaoiton existence | RESOLVED. v5 attempt 4 converged at Q_CS=1.000 exact. | Confirmed in v6 and v7 (sanity check 7.2 reproduces v6.6). | Confirmed. |
| f(J·J) form | RESOLVED IN PRACTICE — v5 uses Numerical-Benchmark form. | v6 uses DeepSeek quartic (`(V²+Q²)² + ...`) instead. Empirically right for v6.6 (H/Q within 5%). LoE-paper-standalone form (Q19) still editorially open. | RESOLVED for electron per Paul Q28: full quartic in EoM, drop from H/Q computation for electron. Different rule for muon/tau (significant) and neutral (essential). |
| 4-fn vs 2-fn ansatz mismatch | RESOLVED. 2-fn for NEUTRAL (Lean theorem), 4-fn for CHARGED (Numerical Benchmark / v5/v6/v7). | Unchanged. | Unchanged. |
| Charge quantization rigorous proof | RESOLVED. Hopf invariant proof complete (Zenodo 20296060). Q25. | Unchanged. | Unchanged. |

---

## Net-change log per sandbox iteration

```text
v5 → v6 net change (2026-05-20):
  RESOLVED:  Q22, Q23, Q24 (three v5 IMMEDIATE; closed by DeepSeek
             2026-05-20 4:00 PM reply); Q25 (Hopf invariant proof).
  PROMOTED:  Q20 (Duda #3) from ACTIVE → NEARLY RESOLVED via v6.6
             empirical construction.
  NEW:       Q26 (5% gap), Q27 (Q_CS_grid sign mismatch).

v6 → v7 net change (2026-05-21):
  RESOLVED:  Q28 (quartic interpretation), Q29 (node count), Q30 (Q(0)
             sign) — all by Paul's 2026-05-21 morning reply.
  DEMOTED:   Q26 (resolved by drop-quartic finding from step 8;
             reduces to Q28 which is now resolved); Q27 (excited-mode
             artifact, not a tunable).
  NEW:       Q31 (initial-guess derivative profile), Q32 (initial
             Lagrange multiplier), Q33 (sign-pinning vs continuation
             method) — all sent in email v7 evening of 2026-05-21.

Active count entering v8 / Paul's Q31/Q32/Q33 reply:
  3 IMMEDIATE (Q31, Q32, Q33)
  4 OPEN      (Q2, Q3, Q6, Q19)
  2 DEMOTED   (Q26, Q27)
  → 9 active questions.
```

---

## Highest-leverage closure path

```text
Paul's Q31/Q32/Q33 reply
    ↓
v8 mode-selector OR continuation wrapper
    ↓
Ground state lands (V₀>0, Q₀>0, A₀<0, J₀<0, ≤4 nodes, H_e/Q ≤ 1.6969)
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
     0b_M6_roadmap.md)
    ↓
M5 returns as foreground engineering track
    (per cardinal rule: M5 is SABER's primary substrate)
```
