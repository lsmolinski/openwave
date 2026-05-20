# M6 / Ouroboros — Roadmap

**Status:** 🔶 Active sandbox evaluation. Gate paths revised after v4 first runs invalidated key v3 claim. Production decision delayed pending 4-fn ODE anchoring.

Last updated: 2026-05-19.

See `0c_model_gates.md` for the G1/G2/G3 production criteria.

---

## Overview

M6 evaluates Paul Werbos's Ouroboros Lagrangian as a second production model
to run in OpenWave in parallel with M5 (Liquid Crystal / Duda). The model is
a classical 2-vector-field theory (A_μ, J_μ on Minkowski spacetime) where
particles emerge as *chaoitons* — time-periodic, localized field solutions
that evade Derrick's theorem via oscillation.

```text
Source framework: Werbos (2017, 2026), Ouroboros system
Key papers: theory/ directory (Calibration, Spectrum, System v2, Law v5, Lean theorem)
Primary contact: Paul Werbos <paul.werbos@gmail.com> — active collaborator since 2026-05-15
Groups cc'd: Models-of-Particles, Jeff Yee, Robert Close, Jarek Duda
```

---

## What has been done

### Phase 0 — Initial evaluation (2026-05-15 to 2026-05-16)

| Task | Outcome |
| --- | --- |
| Read all available Ouroboros papers | ✅ Full corpus in `theory/` |
| Evaluated vs Duda's LdGS challenges | ✅ `0a_background.md` §1-8 |
| Side-by-side M5 vs M6 comparison | ✅ `0a_background.md` §3 |
| Verdict: viable M6 candidate? | ✅ YES — with caveats (`0c_model_gates.md`) |

### Phase 1 — Sandbox v1 (2026-05-17)

| Task | Outcome |
| --- | --- |
| Reproduce electron calibration at (g=1.0625, ω=1.0) | ✅ H/Q = 1.6918 (0.30% gap) → ODE form A4+H1+Q2 calibrates |
| Reproduce lepton mass scan ω^2.22 | ❌ Got ω^2.04, lepton gaps 31-44% |
| Sent findings + THE QUESTION to Werbos | ✅ Email 2026-05-17 (in `0b_sandbox_v1.md`) |
| Key finding: ω^2.22 = log(207)/log(11) → ω values appear fitted, not predicted | ✅ Post-hoc fit signature detected |

### Phase 2 — Sandbox v2 (2026-05-18)

| Task | Outcome |
| --- | --- |
| Werbos replied — active collaboration | ✅ 2026-05-17 "most serious feedback" |
| Read benchmark doc (new 4-fn ODE) | ✅ But benchmark doc had wrong form |
| Attempt Q=0 neutral chaoiton via IVP locked-ansatz (J=-A, Q=-V) | ❌ Bessel J_0 oscillation — no decay |
| Attempt Q=0 neutral chaoiton via BVP (`scipy.solve_bvp`) | ❌ Only trivial + r_max artifacts. No physical chaoiton found. |
| Key finding: locked-ansatz approach is fundamentally incorrect | ✅ Wrong approach identified. Both IVP and BVP confirmed this. |

### Phase 3 — Sandbox v3 (2026-05-18)

| Task | Outcome |
| --- | --- |
| Read 11 new files from Werbos (17:20) | ✅ Key finds: Lean ODE + neutral method |
| Found canonical ODE in Lean theorem (`chaoiton_theorem.lean.txt`) | ✅ Vector Laplacian + slope BCs (fields start at zero, not A₀) |
| Found neutral chaoiton method (`Claude bosonization May18 4PM.txt`) | ✅ α=0, ω=0, single β ODE (A-field off, spin-0) |
| Rebuilt ODE with Lean canonical form (vector Laplacian, slope BCs) | ✅ `m6_v3_chaoiton.py` |
| Neutral chaoiton scan (α=0, ω=0) | ✅ 23 localized solutions. m_χ = 0.508 MeV at λ=1.0 (vs Werbos est. 0.003-0.015 MeV) |
| Electron calibration with Lean ODE | ⚠️ H/Q = 1.723 (1.56% gap). Not yet localized — energy functional needs geometry factor |
| Sent v3 results + Q12/Q13 to Werbos | ✅ Email 2026-05-18 (paste-ready draft) |
| Code search on Zenodo (3 records) | ❌ No Python code found anywhere. DOCX files only — code is "available on request from QAGI LLC" |
| Werbos reply (via DeepSeek, 19:39) | ✅ Q12 PROMISED (code in 1-2 days). Q13 partial answer (locked V→Q, A→J). Q2 hypothesis: leptons = lowest 3 stable Q=1 modes. NEW conflict surfaced: locked vs A=0 descriptions of neutral chaoiton → new Q14. |

### Phase 4 plan — Sandbox v4 (planned 2026-05-19)

| Track | Action | Blocked by |
| --- | --- | --- |
| T1 | Lepton scan with Lean 2-fn ODE → test "lowest 3 stable = leptons" | Not blocked |
| T2 | Neutral ground state via B0 shoot | Not blocked |
| T3 | Werbos's Python source (1-2 days) | Werbos reply |
| T4 | Reconcile Q14 (locked vs A=0) | T3 (code) |

### Phase 5 — Sandbox v4 first runs (2026-05-19)

| Track | Outcome | Status |
| --- | --- | --- |
| T1 | Built `m6_v4_lepton_scan.py`, ran ω∈{1,5,..,45} quick scan + `diag_energy_functional.py` for H-scaling test. Found: 2-fn Lean ODE is WRONG TOOL for charged sector. H decreases with ω (opposite Werbos's bosonization); ω=1.0 marginally fails localization; m_μ/m_e ratio target 207 found 0.04. v3's H_CODE_ELECTRON_CALIB=0.494 also confirmed wrong. Per Q13: 2-fn for NEUTRAL, 4-fn (φ,α,ρ,β) for CHARGED. Bosonization muon@1.1% used 4-fn code, not Lean reduction. | ❌ BLOCKED |
| T2 | Built `m6_v4_t2_neutral_ground.py`. Wide log-scan B₀ ∈ [1e-5, 1.0] × 60 points × 6 λ values; golden-section refinement. NEGATIVE RESULT: H is monotonic in B₀; no minimum exists. All "localized" scan points violate Lean ≤4-node regularity (found 13-29 nodes). Linear part β''+β'/r-β/r²+λβ=0 is J_1(√λr) Bessel — oscillatory not bound. v3's 23 solutions + m_χ=0.508 MeV are ARTIFACTS of permissive absolute-tail check. Q=0 ground-state existence now uncertain. Q14 expanded to 3-way: locked / A=0 / Q_A≈0. | ❌ CLOSED-NEGATIVE |
| T5 (new) | Manual 4-fn ODE extraction via general-purpose agent across v5 paper + Numerical Benchmark §3 + bosonization logs. Found canonical 4-fn ODE: Δ_r V=Q, Δ_r A=J, Δ_r Q=V+m_J²Q+λQ(Q²−J²), Δ_r J=A−m_J²J−λJ(Q²−J²) with toroidal Δ_r=f''+f'/r (NO -f/r²), value BC f'(0)=0, toroidal volume (2π)²R·r dr. Calibration anchor: H/Q=1.6969 at (g=1.0625, λ=1.0, ω=1.0, V₀=A₀=Q₀=J₀=0.1). Confidence: medium. Three gaps: ω-in-static, f(s) two-form, R unspecified. | ✅ DONE |
| NEW | Werbos 2026-05-19 10:51 AM sent DM paper draft for ApJ review citing Griesi v3 m_χ=0.508 MeV as load-bearing numerical input. Rodrigo replied 11:30ish with T2 invalidation + 4-fn extraction + honest code ask. Werbos reply 2026-05-19 1:49 PM (via DeepSeek): resolved all 3 implementation gaps (m_eff² formula, f(s) two-form mapping, R cancels in H/Q). Endorsed our Q14 → canonical Q=0 is Q_A≈0, Q_J≠0. Code "as soon as I can" (no committed date). | ✅ REPLY RECEIVED |
| T6 (new) | Built `m6_v4_4fn.py` with Werbos's m_eff² substitution. Tried RK45, LSODA, RK45+max_step+blowup_event. Scanned (m_J², λ_bench) at canonical V₀=A₀=Q₀=J₀=0.1 across 40 grid points. RESULT: ALL 40 BLOW UP. Confirms 4-fn ODE is eigenvalue/shooting problem; generic IVP from canonical initial conditions doesn't find the bound state — that requires the SPECIFIC calibrated (m_J², λ_bench) which we don't yet have. | ❌ NEGATIVE |
| T6→A (new) | Built `m6_v4_4fn_bvp.py` with `scipy.solve_bvp` + eigenvalue parameter(s). 1-eigenvalue mode converges to non-trivial (V,Q) bound state, but (A,J) collapses to zero spontaneously. Eigenvalue r_max-dependent (Bessel-zero artifact). 2-eigenvalue mode fails universally with singular Jacobian. Q_CS=1 is a topological/integral constraint not a BVP-compatible boundary condition. | ⚠️ PARTIAL |
| T7 (new) | Built `m6_v4_4fn_shoot.py` using r_blowup distance as continuous figure of merit. Coarse 2D scan (m_J², λ_bench) ∈ [0.5,8] × [0.1,10] at V₀=A₀=Q₀=J₀=0.1: max r_blowup=6.44, well below r_max=15. m_J²=1.0 row flat (m_eff²=0 degenerate). Amplitude shoot: r_blowup monotone decreasing in A_0; no resonance. Conclusion: NO bound state in symmetric V₀=A₀=Q₀=J₀ regime anywhere we scanned. Werbos's calibration must use asymmetric initial values OR different parameters. | ❌ NEGATIVE (definitive) |
| EMAIL (out) | Email sent ~5 PM 2026-05-19: T6/T6→A/T7 findings + ask for specific m_J², λ_bench, V₀, A₀, Q₀, J₀ values | ✅ SENT |
| EMAIL (in) | Werbos reply ~4:21 PM (via DeepSeek): m_J²≈0.5, λ_bench=1.0, V₀=Q₀=+0.1, A₀=J₀=-0.1 (asymmetric helicity gives Q_CS=1). m_eff²=-0.5 correct. Shooting strategy: *"decay rate at infinity as target"* (algorithm not specified). + v2 of DM paper sent for review. | ✅ RECEIVED |
| T8 (new) | Re-ran 4fn IVP, shoot, BVP scripts with Werbos's asymmetric helicity. FINDING: helicity is necessary (Q_CS goes from 0 → non-zero) but NOT SUFFICIENT. IVP still blows up (Q_CS=4483, not 1). Amplitude shoot: ~12% better than symmetric, but still monotone. BVP with 1-eigenvalue + asymmetric init guess: still collapses to (V,Q) sub bound state at m_J²=5.51 regardless of init. A,J drift to 0 without a BC forcing them. 2-eigenvalue BVP: singular Jacobian. | ⚠️ PARTIAL |
| T9 (new) | 2026-05-20: Built `sandbox_v4/m6_v4_4fn_newton.py`. Two-stage optimizer over 6 vars (m_J², λ, \|V₀\|, \|A₀\|, \|Q₀\|, \|J₀\|) with helicity signs locked. Stage 1 DE global + Stage 2 Nelder-Mead local. Pre-scan: 2D grid at Werbos's exact \|amps\|=0.1 (324 pts, m_J²∈[0.05,10], λ∈[0.01,10], r_max=15): 0/324 reach r_max with peak<5. Max r_reached=7.24, peak=95. Wider scan (576 pts up to m_J²=50, λ down to 0.001, r_max=30): 0/576, max r_reached=8.06. Stage 1 6-var DE (6528 evals): r_reached=8.4, Q_CS=58662 — delayed-blowup regime, not a bound state. SIX independent attempts (T6 / T6→A 1-eig / T6→A 2-eig / T7 / T8 / T9) all confirm: Q_CS=1 chaoiton UNREACHABLE via forward-IVP from value-BC origin in any parameter region we have explored. | ❌ NEGATIVE DEFINITIVE |
| TRIG (new) | 2026-05-20 PM: Paul email — new compact ApJ paper *"The Neutral Chaoiton: A Dark Matter Candidate from the Ouroboros Lagrangian"* arrived. Distinct from v4 DM (Zenodo 20298669); cites both v8 LoE (Zenodo 20313063) and v4 DM. Acknowledges Griesi for helicity-structure dialogue. Reference [14] = *"Griesi & AI in prep."* Paul holds "TRULY new" Zenodo upload pending OUR T9 / ground-state numbers. With T9 negative, the m_χ + m_J + σ/m numbers in his Section 4 depend on first resolving the shooting-algorithm question. | ✅ NOTED |
| EMAIL (out) | Email Paul ~1:00 PM with definitive T9 negative + sharp algorithmic ask listing 4 candidate algorithms. | ✅ SENT |
| DUDA (in) | Public reply 2026-05-20 ~1:15 PM on models-of-particles list to Paul's v8 LoE Zenodo announce. 4 substantive technical critiques: (1) G_μν undefined; (2) f(J·J) unspecified; (3) construction not shown for electron; (4) two-charge Coulomb derivation absent. Closing: *"LLM-generated word salad"*. (3) is what T9 just demonstrated empirically — Werbos has the algorithm offline but never published it. Now-resolved internally (see below); still unresolved in any public document. | ⚠️ NOTED |
| EMAIL (in) | Werbos algorithm reply ~2:00 PM via DeepSeek. Confirms forward-IVP doesn't work. Specifies algorithm:<br>• Collocation BVP via `scipy.integrate.solve_bvp`<br>• Q_CS=1 as INTEGRAL CONSTRAINT (not BC)<br>• Two free eigenvalues: ω and Lagrange multiplier λ<br>• Robin BCs at R_max matching K_0(κr) decay<br>• Origin BCs: derivative=0 only; values free<br>• Initial profile: exp(-r) decay shapes with helicity signs (V,Q=+0.1·exp(-r); A,J=-0.1·exp(-r))<br>• ω=λ=1 initial guess; Gelfand-Fomin post-convergence<br>Method is fundamentally different from any of T6-T9. Should be implementable in scipy. | ✅ RECEIVED |
| SANDBOX v5 | `sandbox_v5/m6_v5_4fn_lambda_bvp.py` built. Werbos's actual algorithm: collocation BVP via `scipy.integrate.solve_bvp` with Q_CS=1 as integral constraint, free eigenvalues ω + λ_LM (Lagrange multiplier from H' = H - λ·Q), Robin BCs at R_max, derivative-only origin BCs + V(R_MIN)=0.1 anti-collapse, non-proportional A,J initial profile. Attempt 4 (final): `solve_bvp.status` = 0 CONVERGED. ω = 1.047 (vs Werbos 1.0, 4.7% over), m_eff² = -0.596 (vs Werbos -0.5, 19% over), Q_CS = 1.000 exact, λ_LM = -1.212 (first time pinned). FIRST EVER clean convergence to a Q_CS=1 chaoiton via Werbos's actual method. v4's "Q_CS=1 unreachable" negative INVERTED. OPEN: H/Q_CS = 52.64 vs Werbos's 1.6969 (31× off, consistent across V_norm sweep) — most likely a Q_CS or H normalization convention mismatch. Also: A-mode has 17 nodes (excited state, not ground state). See `0b_sandbox_v5.md` for full v5 recipe + diagnostics. | ⚠️ PARTIAL |
| NEXT | Email Paul v4: thank-you for the algorithm + three specific normalization questions (Q_CS convention, H functional coefficients, sample converged profile). Hold ApJ paper Zenodo upload pending these answers. After Paul replies → build sandbox_v6 with corrected conventions; if normalization is the only gap, the electron calibration H/Q=1.6969 lands within hours. | 🔶 DRAFTING |

---

## Current state (2026-05-19, afternoon)

| Item | Status |
| --- | --- |
| Werbos collaboration | 🔶 Active — TWO emails today. AM: ApJ paper draft for review. PM: 3-gap clarifications + Q14 resolution + ApJ-revision intent. |
| Email out (v4 findings + code ask) | ✅ Sent 2026-05-19 AM with T2 invalidation + 4-fn surface area |
| Werbos reply (3 clarifications) | ✅ Received 2026-05-19 1:49 PM. m_eff²=m_J²-ω², f(s) mapping, R cancels. Q14 → Q_A≈0/Q_J≠0. |
| Werbos's promised code | ⚠️ "as soon as I can" — no committed date in PM reply |
| T6 — 4-fn IVP attempt | ❌ NEGATIVE — all 40 (m_J², λ) scan points blow up. Confirms eigenvalue structure; can't anchor with IVP from canonical initial conditions. |
| Gate G1 (lepton scan) | ❌ BLOCKED — 4-fn implementation pending T6→A (BVP) and T7 (shoot) |
| Gate G2 (neutral m_χ ground state) | ❌ v3 INVALIDATED. Werbos endorsed Q_A≈0 canonical; still need 4-fn anchor before running it. |
| Gate G3 (discrete ω selection) | ⚠️ Empirical-via-G1 path still blocked (waiting for G1 path). |
| M6 production in Taichi | 🚧 Production decision delayed past 2026-05-27 if T6→A and T7 land this week; or 2026-05-30+ if Werbos's code is the unblock |
| M5 / Liquid Crystal | 🔶 Active — M5.4 substrate migration queued |

---

## Next steps

### Immediate (waiting on Werbos's reply to 2026-05-19 ask)

| Step | Action | Gate | ETA |
| --- | --- | --- | --- |
| 1 | Receive Werbos response to T2 invalidation + 4-fn ODE ask. If code arrives: implement, run lepton scan + Q=0 ground state. If ODE-form only: implement benchmark 4-fn, calibrate to H/Q=1.6969, then run scans. If no response in ~24-48h: solo 4-fn implementation. | G1+G2 | 2026-05-20/21<br>• code → unblocks all<br>• ODE-form only → still need full scan code (~2 day)<br>• no response → fall back to T5 code (Path A) |

### Can run now (solo path, if Werbos slips)

| Step | Action | Gate | ETA |
| --- | --- | --- | --- |
| 2 | Solo T5→code: implement 4-fn benchmark ODE in `sandbox_v4/`. Calibrate against H/Q = 1.6969 first; only proceed to lepton scan + Q=0 scan once anchored. Three gaps need empirical iter: ω in static ODE (try m_J² ← m_J² − ω² substitution); f(s) v5 vs benchmark form; toroidal R value (try R=1). | G1+G2 | 1-2 days |
| 3 | After 4-fn anchored: run charged lepton scan over ω∈[0.5, 50] | G1 | +1 day post-anchor |
| 4 | After 4-fn anchored: run Q=0 scan with Q_A≈0 constraint (the untested 3rd description from Q14) | G2 | +1 day post-anchor |

### If G1 + G2 PASS (GO decision)

| Step | Action | ETA |
| --- | --- | --- |
| 4 | Scaffold M6 in Taichi: Vector(4) × 2 substrate (A, J); Lorenz constraint enforcer; Chern-Simons charge kernel; mirror M5's rendering pipeline. | post-M5.4 |
| 5 | Gate 1 Taichi: Maxwell limit → A-field alone = standard EM | Week 1 of M6 build |
| 6 | Gate 2 Taichi: charge quantization → Q[A,J] integer on seeded configs | Week 2 |
| 7 | Gate 3 Taichi: chaoiton existence → localized time-periodic solution at same (g,λ,ω) as sandbox | Weeks 3-6 |

### Parked (post-G1+G2)

| Step | Action | Notes |
| --- | --- | --- |
| C | 3-body proton bound state | V(R) ~ -C/R⁶ classical 3-body problem; deferred |

---

## Open questions (summary — see 0b_sandbox_v5.md for full v5 tracker)

**v5 IMMEDIATE** (block v6, asked in Werbos email v4 2026-05-20 PM):

| ID | Question | Why it matters |
| --- | --- | --- |
| Q22 | Q_CS normalization convention. | Pin the factor that accounts for the 31× H/Q gap between v5's 52.64 and Werbos's stated 1.6969. |
| Q23 | Exact H functional (kinetic coefficient, ω-kinetic placement, cross-term signs). | Pins the Lagrange-multiplier coefficient. |
| Q24 | Sample converged profile from one of Werbos's runs (5-10 grid points for V, A, Q, J). | Selects ground-state basin vs excited 17-node mode. |

**v5 ACTIVE** (close-able via v6):

| ID | Question | Why it matters |
| --- | --- | --- |
| Q20 | Duda critique #3 — "construction not shown, need ansatz + minimize energy". | HALF ADDRESSED by v5 (method works in script). v6 closes by hitting H/Q=1.6969 exactly. |

**v5 STILL OPEN** (not blocking v6):

| ID | Question | Why it matters |
| --- | --- | --- |
| Q2 | Discrete ω selection mechanism | Empirical-via-lepton-scan once v6 anchors. Analytic proof deferred per Werbos. |
| Q3 | Analytical ω = 2mc²/ℏ derivation | Calibration only (1.2%). |
| Q6 | QCD reconciliation (3-chaoiton proton) | Future v7+. |
| Q19 | f(J·J) explicit form in LoE paper standalone (Duda critique #2). | Editorial — one-line LoE revision by Werbos fixes. |
| Q21 | Two-chaoiton Coulomb derivation (Duda critique #4). | Future sandbox v7+. ARCHIVED for current scope. |

**v5 RESOLVED** (closed by v5 or earlier):

| ID | Resolution |
| --- | --- |
| Q1, Q9, Q10, Q11, Q13 | Earlier rounds (v1-v4). |
| Q12 (Werbos Python code) | DEMOTED — v5 algorithm description + DeepSeek-write offer = same unblock. |
| Q14 (canonical Q=0) | Q_A≈0 / Q_J≠0 (Werbos 1:49 PM 2026-05-19). |
| Q15 (m_eff² substitution) | RESOLVED post-v5. m_eff² = -0.596 in v5 attempt 4 confirms negative is correct. |
| Q16 (m_J², λ_bench) | m_J²≈0.5, λ_bench=1.0 (Werbos 4:21 PM 2026-05-19); v5 uses these. |
| Q17 (shooting algorithm) | RESOLVED — collocation BVP per Werbos 2026-05-20 PM email. |
| Q25 (Hopf invariant proof rigorous?) | RESOLVED — Zenodo 20296060 supplies the two missing lemmas. Charge quantization is now a theorem of differential topology. |

**v5 ARCHIVED** (unfalsifiable / future scope):

| ID | Why archived |
| --- | --- |
| Q4 (Single vs two-field ontology, G_μν / J undefined) | = Duda critique #1. Aesthetic preference. If math matches observation, two primary fields is just a description, not theory-killer. |
| Q7 (Cold-fusion citation trail) | Historical, not physics. |
| Q21 (Two-chaoiton Coulomb derivation, Duda #4) | (Also listed in STILL OPEN as future-v7 milestone.) |

Active count entering v6: **3 IMMEDIATE + 1 ACTIVE + 5 OPEN = 9 questions**.
Highest-leverage closure: Q22-Q24 (Werbos reply) → v6 → Q20 (Duda #3).

---

## Hardest pieces (summary — see 0b_sandbox_v5.md for live tracker)

| Hardest piece | Status post-v5 (2026-05-20) |
| --- | --- |
| Forward-IVP method family | RESOLVED. Confirmed wrong tool by both us (T6-T9) and Werbos. Replaced with collocation BVP in v5. |
| Q_CS=1 enforcement mechanism | RESOLVED. Auxiliary integral state I with BC I(R_max)=1/(2π). Exact Q_CS=1.000 achieved. |
| Q_CS=1 chaoiton existence | RESOLVED. v5 attempt 4 converges (`solve_bvp` status=0) at ω=1.047, m_eff²=-0.596. First empirical chaoiton in OpenWave. |
| Electron H/Q = 1.6969 calibration | NEW OPEN. v5 lands H/Q=52.64 (31× off) across V_norm sweep — structural normalization gap. Q22 (Q_CS convention) most likely closes this. |
| Ground-state vs excited mode selection | NEW OPEN. v5 attempt 4 found Q_CS=1 chaoiton but with A in 17-node excited mode (Lean ≤4 node spec). Q24 (sample profile from Werbos) resolves immediately. |
| Lagrange-multiplier ODE correction coefficient | NEW OPEN. v5 derived coefficient = 1. If Werbos's H kinetic factor differs, this may need adjustment. Q23 pins it. |
| V(M) potential form | UNRESOLVED — shared bottleneck with M5 (Duda). Not in v5/v6 scope. |
| f(J·J) form | RESOLVED IN PRACTICE — v5 uses Numerical Benchmark form f(s) = ½m_J²s + ¼λs². Q19 (LoE paper editorial gap) remains for Werbos to fix in next LoE revision. |
| 4-fn vs 2-fn ansatz mismatch | RESOLVED. 2-fn for NEUTRAL (Lean theorem), 4-fn for CHARGED (Numerical Benchmark / v5). Both are canonical, not reductions. |
| ω quantization mechanism | OPEN (Q2). v6 lepton scan provides empirical test once H/Q calibration anchors. Analytic proof deferred per Werbos's own admission. |
| Lepton mass spectrum | BLOCKED on v6. Once H/Q=1.6969 lands, run ω-sweep [0.5, 50] to test lowest-3-stable hypothesis. |
| Neutral m_χ true ground state | BLOCKED on v6. Once anchored, run Q_A≈0 scan — feeds Section 4 of ApJ Neutral Chaoiton paper. |
| Two-chaoiton Coulomb derivation | NOT ADDRESSED. Future sandbox v7+ (Q21, Duda critique #4). Out of current scope. |
| Charge quantization rigorous proof | RESOLVED. Hopf invariant proof complete (Zenodo 20296060). Q25. |

---

## Resources

**Research docs:**

| Doc | Contents |
| --- | --- |
| `0a_background.md` | Technical evaluation (§1-10) + updates (§11) |
| `0b_sandbox_v1.md` | v1 plan, results, emails, Q-list |
| `0b_sandbox_v2.md` | v2 plan, IVP+BVP attempts, Q-list |
| `0b_sandbox_v3.md` | v3 plan, Lean ODE, neutral scan, Q-list |
| `0b_sandbox_v4.md` | v4 plan, T1-T9 negative results, forward pointer to v5 |
| `0b_sandbox_v5.md` | v5 plan (current), Werbos algorithm reply, T10 partial success, question tracker |
| `0c_model_gates.md` | G1/G2/G3 production decision criteria |
| `0c_roadmap.md` | This file |

**Sandbox scripts:**

| Folder | Contents |
| --- | --- |
| `sandbox_v1/` | 96-variant sweep, calibration, mass scan (wrong ODE) |
| `sandbox_v2/` | IVP + BVP locked-ansatz attempts (superseded) |
| `sandbox_v3/` | Lean ODE, neutral β, calibration (last completed) |
| `sandbox_v4/` | T1 lepton scan (blocked, 2-fn wrong tool); T2 neutral ground state (negative); `diag_energy_functional`; T5 4-fn extracted; T6-T9 negative-result reference |
| `sandbox_v5/` | `m6_v5_4fn_lambda_bvp.py` — Werbos algorithm: collocation BVP + Lagrange multiplier, status=0 converged on Q_CS=1 chaoiton |

**Theory papers:**

| Folder | Contents |
| --- | --- |
| `theory/` | All Werbos corpus (v1-v8 LoE, Lean theorem, calibration, spectrum, Numerical Benchmark, bosonization logs, Hopf invariant proof completion, Neutral Chaoiton ApJ paper, Dark Matter paper) |
