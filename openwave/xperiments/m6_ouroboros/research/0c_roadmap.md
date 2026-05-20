# M6 / Ouroboros — Roadmap

**Status:** 🔶 v5 PARTIAL SUCCESS — first-ever Q_CS=1 chaoiton converged via Werbos's collocation BVP algorithm at ω=1.047, m_eff²=-0.596, Q_CS=1.000 exact (`solve_bvp.status=0`). Open: H/Q=52.64 vs Werbos's stated 1.6969 (31× off, normalization-convention mismatch). Email v4 sent 2026-05-20 PM with three specific normalization questions (Q22 Q_CS convention, Q23 H functional, Q24 sample profile); awaiting reply to build sandbox_v6 and close calibration.

Last updated: 2026-05-20 evening.

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

### Phase 5 — Sandbox v4 attempts (2026-05-19 to 2026-05-20 PM)

Six independent forward-IVP-family attempts (T1, T2, T6, T6→A, T7, T8, T9) plus
two key Werbos clarification emails. **Net result: forward-IVP from value-BC
origin DOES NOT reach a Q_CS=1 chaoiton in any parameter region we explored.**
This negative result was load-bearing — it told Werbos his canonical method had
never been written down in published form, and prompted the v5 algorithm reply.

| Track | Outcome | Status |
| --- | --- | --- |
| T1 | Built `m6_v4_lepton_scan.py`, ran ω∈{1,5,..,45} quick scan + `diag_energy_functional.py`. Found 2-fn Lean ODE is WRONG TOOL for charged sector: H decreases with ω (opposite Werbos's bosonization); ω=1.0 marginally fails localization; m_μ/m_e ratio target 207 found 0.04. v3's H_CODE_ELECTRON_CALIB=0.494 also confirmed wrong. Per Q13: 2-fn for NEUTRAL, 4-fn (φ,α,ρ,β) for CHARGED. | ❌ BLOCKED |
| T2 | Built `m6_v4_t2_neutral_ground.py`. Wide log-scan B₀ ∈ [1e-5, 1.0] × 60 points × 6 λ values; golden-section refinement. NEGATIVE: H is monotonic in B₀; no minimum exists. All "localized" scan points violate Lean ≤4-node regularity (found 13-29 nodes). v3's 23 solutions + m_χ=0.508 MeV are ARTIFACTS of permissive tail check. Q14 expanded to 3-way: locked / A=0 / Q_A≈0. | ❌ CLOSED-NEGATIVE |
| T5 (new) | Manual 4-fn ODE extraction via general-purpose agent across v5 paper + Numerical Benchmark §3 + bosonization logs. Found canonical 4-fn ODE: Δ_r V=Q, Δ_r A=J, Δ_r Q=V+m_J²Q+λQ(Q²−J²), Δ_r J=A−m_J²J−λJ(Q²−J²) with toroidal Δ_r=f''+f'/r. Calibration anchor: H/Q=1.6969 at (g=1.0625, λ=1.0, ω=1.0). Three gaps: ω-in-static, f(s) two-form, R unspecified. | ✅ DONE |
| EMAIL (in) 2026-05-19 1:49 PM | Werbos reply (via DeepSeek) resolved all 3 T5 gaps: m_eff² = m_J²−ω²; f(s) v5 ≡ benchmark with m_J²=0, λ=4g; toroidal R cancels in H/Q ratio. Endorsed Q14 → canonical Q=0 is Q_A≈0, Q_J≠0. Code "as soon as I can" (no committed date). | ✅ RECEIVED |
| T6 (new) | Built `m6_v4_4fn.py` with Werbos's m_eff² substitution. Scanned (m_J², λ_bench) at canonical V₀=A₀=Q₀=J₀=0.1 across 40 grid points. RESULT: ALL 40 BLOW UP. Confirms 4-fn ODE is eigenvalue/shooting problem; generic IVP from canonical initial conditions doesn't find the bound state. | ❌ NEGATIVE |
| T6→A (new) | Built `m6_v4_4fn_bvp.py` with `scipy.solve_bvp` + eigenvalue parameter(s). 1-eigenvalue mode converges to non-trivial (V,Q) bound state, but (A,J) collapses to zero spontaneously. Eigenvalue r_max-dependent (Bessel-zero artifact). 2-eigenvalue mode fails universally with singular Jacobian. Q_CS=1 is integral constraint, not BC. | ⚠️ PARTIAL |
| T7 (new) | Built `m6_v4_4fn_shoot.py` using r_blowup as continuous figure of merit. Coarse 2D scan (m_J², λ_bench) ∈ [0.5,8] × [0.1,10] at V₀=A₀=Q₀=J₀=0.1: max r_blowup=6.44, well below r_max=15. Amplitude shoot monotone decreasing in A_0; no resonance. NO bound state in symmetric ansatz anywhere scanned. | ❌ NEGATIVE (definitive) |
| EMAIL (in) 2026-05-19 4:21 PM | Werbos reply (via DeepSeek): asymmetric helicity gives Q_CS=1. V₀=Q₀=+0.1, A₀=J₀=-0.1; m_J²≈0.5, λ_bench=1.0, m_eff²=-0.5 (negative is correct, time-periodic). Shooting strategy: *"decay rate at infinity as target"* (algorithm not detailed). + v2 of DM paper sent for review. | ✅ RECEIVED |
| T8 (new) | Re-ran 4fn IVP, shoot, BVP scripts with asymmetric helicity. Helicity is NECESSARY (Q_CS goes 0 → 4483) but NOT SUFFICIENT. IVP still blows up, no Q_CS=1. Amplitude shoot ~12% better than symmetric but still monotone. BVP with 1-eigenvalue + asymmetric init: still collapses to (V,Q) sub bound state at m_J²=5.51. A,J drift to 0. | ⚠️ PARTIAL |
| T9 (new) 2026-05-20 PM | Built `m6_v4_4fn_newton.py`. Two-stage optimizer over 6 vars with helicity locked. Pre-scan 324 pts at \|amps\|=0.1: 0/324 reach r_max with peak<5. Wider scan 576 pts: 0/576. Stage 1 6-var DE (6528 evals): r_reached=8.4, Q_CS=58662 — delayed-blowup, not bound state. **DEFINITIVE NEGATIVE:** across 6 attempts, Q_CS=1 chaoiton UNREACHABLE via forward-IVP. Werbos's shooting must be a different method family. | ❌ DEFINITIVE NEGATIVE |

### Phase 6 — Sandbox v5 (2026-05-20 PM)

Triggered by Paul's algorithm-clarification email after T9 negative result. Built
`sandbox_v5/m6_v5_4fn_lambda_bvp.py` implementing Werbos's actual method:
collocation BVP via `scipy.integrate.solve_bvp` with Q_CS=1 as INTEGRAL CONSTRAINT.
**Result: first-ever clean convergence to a Q_CS=1 chaoiton.** Open: 31× H/Q
normalization gap blocking electron calibration.

| Event | Outcome | Status |
| --- | --- | --- |
| TRIG 2026-05-20 PM | Paul email — new compact ApJ paper *"The Neutral Chaoiton: A Dark Matter Candidate from the Ouroboros Lagrangian"* arrived. Distinct from v4 DM (Zenodo 20298669); cites both v8 LoE (Zenodo 20313063) and v4 DM. Acknowledges Griesi for helicity-structure dialogue. Reference [14] = *"Griesi & AI in prep."* Paul holds "TRULY new" Zenodo upload pending OUR ground-state numbers. | ✅ NOTED |
| EMAIL (out) ~1:00 PM | Email Paul with definitive T9 negative + sharp algorithmic ask listing 4 candidate algorithms (backward integration; BVP with Q_CS Lagrange multiplier; tabulated init profile; multi-shot connection). | ✅ SENT |
| DUDA (in) ~1:15 PM | Public reply on models-of-particles list to Paul's v8 LoE announce. 4 substantive technical critiques: (1) G_μν undefined; (2) f(J·J) unspecified; (3) construction not shown for electron; (4) two-charge Coulomb derivation absent. Closing: *"LLM-generated word salad"*. Our internal read: #1 archived (unfalsifiable preference), #2 editorial, #3 v6 closes empirically, #4 future v7+. | ⚠️ NOTED |
| EMAIL (in) ~2:00 PM | Werbos algorithm clarification via DeepSeek. Confirms forward-IVP doesn't work ("consistent with our experience"). Specifies algorithm:<br>• Collocation BVP via `scipy.integrate.solve_bvp` (or `scipy.optimize.root` method='lm')<br>• Q_CS=1 as INTEGRAL CONSTRAINT, not BC<br>• Two free eigenvalues: ω and Lagrange multiplier λ (from H' = H − λ·Q)<br>• Robin BCs at R_max: V'+k·V=0, k=√(ω²-m_J²)<br>• Origin BCs: derivative=0 only; values FREE<br>• Initial profile: V=+0.1·exp(-r), A=-0.1·exp(-r), Q=+0.1·exp(-r), J=-0.1·exp(-r); ω=λ=1 init guess<br>• Post-convergence: Gelfand-Fomin conjugate-point test | ✅ RECEIVED |
| v5 attempt 1 | No λ in ODEs, Q_CS=1 via integral. solve_bvp converged but (V,Q) collapsed to noise; A,J huge; H=12317. Wrong basin. | ❌ wrong basin |
| v5 attempt 2 | Added λ-corrections to A,J via δQ_CS variation + V(R_MIN)=0.1 anti-collapse. Converged with V,A,Q,J non-trivial. H=458, 16 nodes V/Q. Mid-basin progress. | ⚠️ excited mode |
| v5 attempt 3 | Non-proportional A,J initial guess (J on exp(-1.5r) to break Q_CS≡0 degeneracy) + 20k node budget. Converged. ω=0.89, m_eff²=-0.28. H=32.8, 3 nodes V/Q, 0 nodes A/J. Closer to ground state. | ⚠️ partial |
| v5 attempt 4 ★ | r_max=12, n_grid=400, max_nodes=50000, tol=5e-3. **`solve_bvp.status = 0` CONVERGED** (28k final nodes, max residual 1.3e-4). ω=1.047 (vs Werbos 1.0, 4.7% over). m_eff²=-0.596 (vs Werbos -0.5, 19% over). λ_LM=-1.212 (first time pinned). Q_CS=1.000 exact. H/Q-from-I and H/Q-from-grid agree to 0.01% — solver self-consistent. **WIN:** first-ever Q_CS=1 chaoiton in OpenWave; v4's "unreachable" inverted. **OPEN:** H/Q=52.64 vs target 1.6969 (31× off); A has 17 nodes (excited mode, not ground state); tail=0.17. | ⚠️ PARTIAL ★ |
| V_norm sweep | Tested V_norm ∈ {0.5, 0.3, 0.2, 0.15, 0.1, 0.05, 0.03, 0.01}: H/Q_CS bottoms at 52.64 at V_norm=0.1; never approaches 1.6969 anywhere. 31× ratio is STRUCTURAL across all modes, not a tuning issue. Most likely Q_CS or H normalization-convention mismatch (Paul: (1/4π²)·∫ε A ∂J d³x; us: 2π·∫r·(A·J'-J·A')dr radial-toroidal form). | ⚠️ confirms structural gap |
| Q_CS = 1.000 exact achieved | First-ever Q_CS=1 chaoiton in OpenWave. v4's "unreachable" definitively inverted. Chaoiton EXISTS as a `solve_bvp` solution — only the calibration mapping remains open. | ✅ WIN |
| EMAIL (out) v4 ~PM | Sent Paul a thank-you + three specific normalization questions: Q22 (exact Q_CS normalization convention); Q23 (exact H functional kinetic coefficients, ω-kinetic placement, cross-term signs); Q24 (sample converged profile from one of his runs, even 5-10 grid points). Stated explicit commitment to hold ApJ Zenodo upload pending answers. | ✅ SENT |

### Phase 7 plan — Sandbox v6 (planned, awaiting Paul reply to Q22/Q23/Q24)

| Step | Action | Triggers | Estimate |
| --- | --- | --- | --- |
| 1 | Receive Paul reply on Q22/Q23/Q24 normalization questions | required — without this, v6 cannot start. DeepSeek-write-Python fallback if Paul slips ~48h. | 24-48h |
| 2 | Fork v5 → `sandbox_v6/m6_v6_4fn_calibrated_bvp.py` | once #1 lands; mechanical | 5 min |
| 3 | Apply Q22 fix (Q_CS normalization factor) | scalar coefficient adjustment on auxiliary integral state | 15 min |
| 4 | Apply Q23 fix (H functional coefficients + Lagrange-multiplier coefficient in A,J equations) | edit to ODE rhs | 15 min |
| 5 | Apply Q24 warm-start (Paul's sample profile as initial guess, interpolated onto our grid) | replaces exp(-r) seed; resolves ground-state basin selection | 15 min |
| 6 | Re-converge attempt-4-style run | should give H/Q ≈ 1.6969 if normalization is the only gap; verify ≤4 nodes per Lean | 5 min solve + verify |
| 7 | Lepton scan ω ∈ [0.5, 50] at calibrated (m_J², λ_bench) fixed | tests "lowest 3 stable modes = leptons" hypothesis; expect muon ω≈12.78 (1.1% gap), tau ω≈40.7 | ~1 hour |
| 8 | Q_A ≈ 0 neutral chaoiton scan for m_χ | lands DM mass + m_J mediator mass + σ/m self-interaction; feeds ApJ Section 4 | ~1 hour |
| 9 | Gelfand-Fomin stability check (conjugate-point test on converged second variation) | confirms ground state, not excited mode | ~30 min |
| 10 | Handoff (m_χ, m_J, σ/m, relic abundance) to Paul for ApJ Section 4; lift Zenodo upload hold | "Griesi & AI in prep" reference becomes real citation | done |

---

## Current state (2026-05-20 evening)

| Item | Status |
| --- | --- |
| Werbos collaboration | 🔶 Active — three emails 2026-05-20: AM new compact ApJ Neutral Chaoiton paper + DM v4 deposit; PM algorithm reply (~2:00 PM via DeepSeek); evening "DeepSeek agrees we should upload nothing new until we hear from you". |
| Email v4 (out) — Q22/Q23/Q24 | ✅ Sent 2026-05-20 PM. Three specific normalization questions: Q_CS convention, H functional coefficients + ω-kinetic + cross-term signs, sample converged profile from one of his runs. |
| Werbos reply on v4 | 🔶 INCOMING — DeepSeek-mediated reply being reviewed next. Expected to unblock v6. |
| sandbox_v5 Q_CS=1 chaoiton | ✅ DEMONSTRATED. `solve_bvp.status=0` at ω=1.047, m_eff²=-0.596, Q_CS=1.000 exact, λ_LM=-1.212. First-ever clean convergence in OpenWave. |
| Electron H/Q = 1.6969 calibration | ❌ OPEN. v5 lands H/Q=52.64 (31× off). V_norm sweep confirms structural; most likely Q_CS or H normalization convention mismatch. v6 closes if Q22 lands. |
| Ground-state vs excited mode | ❌ OPEN. v5 attempt 4 found Q_CS=1 chaoiton but with A in 17-node excited mode (Lean spec ≤4 nodes). Q24 sample profile resolves. |
| Hopf invariant proof (charge quantization) | ✅ RESOLVED. Zenodo 20296060 supplies the two missing lemmas. Now a theorem of differential topology. |
| Duda critiques (2026-05-20 thread) | #1 G_μν / single-field ontology — ARCHIVED (unfalsifiable). #2 f(J·J) unspecified — editorial OPEN. #3 construction not shown — HALF-ADDRESSED (v5 method works in script); v6 closes if H/Q=1.6969 lands. #4 two-charge Coulomb — FUTURE v7+. |
| Gate G1 (lepton scan) | ❌ BLOCKED on v6 calibration anchor. v6 ω-sweep [0.5, 50] will test "lowest 3 stable = leptons" empirically. |
| Gate G2 (neutral m_χ ground state) | ❌ BLOCKED on v6 calibration anchor. v6 Q_A≈0 scan lands DM candidate m_χ, m_J, σ/m for ApJ Section 4. |
| Gate G3 (discrete ω selection) | ⚠️ BLOCKED on G1 (empirical-via-lepton-scan). Analytic proof still deferred per Werbos's admission. |
| ApJ Neutral Chaoiton paper | 🔶 HOLD. Reference [14] *"Griesi & AI in prep"*. Paul + DeepSeek explicitly agreed: "upload nothing new until we hear from you". Numbers depend on v6 close. |
| M6 production in Taichi | 🚧 Decision deferred until v6 calibration lands. Earliest realistic: 1-2 days post-Paul-reply if normalization is the only gap. |
| M5 / Liquid Crystal | 🔶 M5.4 substrate migration queued post-v6 close. Cardinal rule: M5 is SABER's primary engineering track; M6 stays parallel-research. |

---

## Next steps

### Immediate (waiting on Paul's reply to Q22/Q23/Q24)

| Step | Action | Gate | ETA |
| --- | --- | --- | --- |
| 1 | Receive Werbos reply on Q22/Q23/Q24. Most likely Q22 alone (Q_CS normalization) closes the 31× gap, but all three answers strengthen v6. | G1+G2 anchor | 24-48h post-email |
| 2 | If Paul slips past ~48h: take DeepSeek up on its offer to write Python code directly — has Werbos's conventions in context. Cross-check v6 against its output. | fallback | 48-72h post-email |
| 3 | If still stuck past ~96h: park M6, return to M5.4. M5 is SABER's primary track. | fallback | 96h+ |

### v6 build sequence (triggered by Paul's reply)

| Step | Action | Gate | Estimate |
| --- | --- | --- | --- |
| 1 | Fork v5 → `sandbox_v6/m6_v6_4fn_calibrated_bvp.py`. | — | 5 min |
| 2 | Apply Q22 normalization fix (scalar factor on auxiliary integral state I, or on H prefactor). | — | 15 min |
| 3 | Apply Q23 H-functional fix (kinetic coefficient + Lagrange-multiplier coefficient in A,J ODE corrections). | — | 15 min |
| 4 | Apply Q24 warm-start (Paul's sample profile interpolated onto grid; replaces exp(-r) seed). | — | 15 min |
| 5 | Re-converge target run. Expect H/Q = 1.6969 ± 0.001 at Q_CS=1, all fields ≤4 nodes. | G1+G2 anchor | 5 min solve |
| 6 | Lepton scan ω ∈ [0.5, 50] with calibrated (m_J², λ_bench) fixed. Expected muon ω≈12.78 (1.1% gap target), tau ω≈40.7. | G1 | ~1 hour |
| 7 | Q_A≈0 neutral chaoiton scan. Lands m_χ, m_J mediator mass, σ/m self-interaction for ApJ Section 4. | G2 | ~1 hour |
| 8 | Gelfand-Fomin conjugate-point stability check. Confirms ground state vs excited. | G3 (empirical) | ~30 min |
| 9 | Handoff (m_χ, m_J, σ/m, Ω_χh²) to Paul. Lift ApJ Zenodo upload hold. "Griesi & AI in prep" becomes a real citation. | — | done |

**Total estimated time post-Paul-reply: 2-3 hours of solo work + solve runs.**

### If G1 + G2 PASS (Taichi production GO)

| Step | Action | ETA |
| --- | --- | --- |
| 1 | Scaffold M6 in Taichi: Vector(4) × 2 substrate (A, J); Lorenz constraint enforcer; Chern-Simons charge kernel; mirror M5's rendering pipeline. | post-M5.4 |
| 2 | Gate 1 Taichi: Maxwell limit → A-field alone = standard EM | Week 1 of M6 build |
| 3 | Gate 2 Taichi: charge quantization → Q[A,J] integer on seeded configs | Week 2 |
| 4 | Gate 3 Taichi: chaoiton existence → localized time-periodic solution at same (g,λ,ω) as sandbox | Weeks 3-6 |

### Parked (post-G1+G2; future v7+ work)

| Step | Action | Notes |
| --- | --- | --- |
| C | Two-chaoiton Coulomb derivation (Duda critique #4) | Integrate H[Φ₁, Φ₂] for two topological charges at distance. New BVP scaffold. |
| D | 3-body proton bound state | V(R) ~ -C/R⁶ classical 3-body problem; deferred |
| E | Hopfion candidates for excited neutrino oscillation states | Liu et al. *Nature Physics* 2026 lab anchor; topology-as-particles frontier |

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
| `sandbox_v3/` | Lean ODE, neutral β, calibration (last completed pre-v4) |
| `sandbox_v4/` | T1 lepton scan (blocked); T2 neutral ground state (negative); `diag_energy_functional`; T5 4-fn extracted; T6-T9 forward-IVP attempts — all negative; reference for definitive "forward-IVP wrong tool" result |
| `sandbox_v5/` | `m6_v5_4fn_lambda_bvp.py` — Werbos's collocation BVP algorithm with Lagrange multiplier λ_LM. `solve_bvp.status=0` converged on Q_CS=1 chaoiton at ω=1.047, m_eff²=-0.596. H/Q=52.64 vs target 1.6969 (open: normalization gap). |
| `sandbox_v6/` | PLANNED — `m6_v6_4fn_calibrated_bvp.py` (forks v5 + Q22/Q23/Q24 fixes from Paul's reply). Target: H/Q=1.6969, then lepton scan + Q_A≈0 DM scan. |

**Theory papers:**

| Folder | Contents |
| --- | --- |
| `theory/` | All Werbos corpus (v1-v8 LoE, Lean theorem, calibration, spectrum, Numerical Benchmark, bosonization logs, Hopf invariant proof completion, Neutral Chaoiton ApJ paper, Dark Matter paper) |
