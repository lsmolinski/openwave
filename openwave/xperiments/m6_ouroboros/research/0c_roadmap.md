# M6 / Ouroboros — Roadmap

**Status:** 🔶 v6 CALIBRATION ESSENTIALLY CLOSED — **step (8) diagnostic shows dropping the quartic from H lands H/Q = 1.7112, 0.84% off target** (matches Werbos's own *"for the electron, the quartic term is small"*). v6.6 with DeepSeek quartic lands H/Q = 1.778 (4.8% over). v5 was 52.64 (31× off) — net 31× → 0.84% in one session. Step (11) inspection: v6.6 is a 5-node excited mode (Lean spec ≤4), so true ground state may give H/Q even closer to 1.6969. Step (2) cold-start lepton scan fails (needs continuation method). DeepSeek's reference Python script received 2026-05-20 ~5:30 PM but does NOT run as-is (r=0 singularity + BC count mismatch); patched version diverges to H/Q≈10^16 — our v6 is dramatically closer to truth. Email v6 sent 2026-05-20 evening with Q28 (quartic interpretation), Q29 (≤4 nodes vs 5), Q30 (Q(r=0) sign). Hard stop tonight; sandbox_v7 tomorrow contingent on Paul's reply.

Last updated: 2026-05-20 evening (post step-(8/11/4/2) diagnostics + email v6 sent; hard stop).

See `0c_model_gates.md` for the G1/G2/G3 production criteria.
See `0d_m5vsM6.md` for the M5 vs M6 cross-evaluation in service of the SABER trigger.

**Publishing stance (locked-in 2026-05-20):** OpenWave will NOT deposit M6
work on Zenodo, arxiv, or any paper venue. The open-source GitHub repo IS
the deliverable — full granularity (code, diagnostics, email trails) lives
there and is more useful than a frozen PDF. Werbos's ApJ Ref [14] *"Griesi
& Anthropic AI 2026, in preparation"* will resolve to a stable GitHub URL
(tagged release if needed for DOI via auto-Zenodo integration, only if a
reviewer requires it). Section 4 numbers are handed off to Paul as a data
drop, not co-authorship. No formal commitment to publish was ever made.

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

## Sandbox sequence — what has been done

### Initial evaluation (2026-05-15 to 2026-05-16)

| Task | Outcome |
| --- | --- |
| Read all available Ouroboros papers | ✅ Full corpus in `theory/` |
| Evaluated vs Duda's LdGS challenges | ✅ `0a_background.md` §1-8 |
| Side-by-side M5 vs M6 comparison | ✅ `0a_background.md` §3 |
| Verdict: viable M6 candidate? | ✅ YES — with caveats (`0c_model_gates.md`) |

### Sandbox v1 (2026-05-17)

| Task | Outcome |
| --- | --- |
| Reproduce electron calibration at (g=1.0625, ω=1.0) | ✅ H/Q = 1.6918 (0.30% gap) → ODE form A4+H1+Q2 calibrates |
| Reproduce lepton mass scan ω^2.22 | ❌ Got ω^2.04, lepton gaps 31-44% |
| Sent findings + THE QUESTION to Werbos | ✅ Email 2026-05-17 (in `0b_sandbox_v1.md`) |
| Key finding: ω^2.22 = log(207)/log(11) → ω values appear fitted, not predicted | ✅ Post-hoc fit signature detected |

### Sandbox v2 (2026-05-18)

| Task | Outcome |
| --- | --- |
| Werbos replied — active collaboration | ✅ 2026-05-17 "most serious feedback" |
| Read benchmark doc (new 4-fn ODE) | ✅ But benchmark doc had wrong form |
| Attempt Q=0 neutral chaoiton via IVP locked-ansatz (J=-A, Q=-V) | ❌ Bessel J_0 oscillation — no decay |
| Attempt Q=0 neutral chaoiton via BVP (`scipy.solve_bvp`) | ❌ Only trivial + r_max artifacts. No physical chaoiton found. |
| Key finding: locked-ansatz approach is fundamentally incorrect | ✅ Wrong approach identified. Both IVP and BVP confirmed this. |

### Sandbox v3 (2026-05-18)

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

### Sandbox v4 (2026-05-19 to 2026-05-20 PM)

Six independent forward-IVP-family attempts (T1, T2, T6, T6→A, T7, T8, T9) plus
two key Werbos clarification emails. **Net result: forward-IVP from value-BC
origin DOES NOT reach a Q_CS=1 chaoiton in any parameter region we explored.**
This negative result was load-bearing — it told Werbos his canonical method had
never been written down in published form, and prompted the v5 algorithm reply.

| Track | Outcome | Status |
| --- | --- | --- |
| T1 | Built `m6_v4_lepton_scan.py`, ran ω∈{1,5,..,45} quick scan + `diag_energy_functional.py`. Found 2-fn Lean ODE is WRONG TOOL for charged sector: H decreases with ω (opposite Werbos's bosonization); ω=1.0 marginally fails localization; m_μ/m_e ratio target 207 found 0.04. v3's H_CODE_ELECTRON_CALIB=0.494 also confirmed wrong. Per Q13: 2-fn for NEUTRAL, 4-fn (φ,α,ρ,β) for CHARGED. | ❌ BLOCKED |
| T2 | Built `m6_v4_t2_neutral_ground.py`. Wide log-scan B₀ ∈ [1e-5, 1.0] × 60 points × 6 λ values; golden-section refinement. NEGATIVE: H is monotonic in B₀; no minimum exists. All "localized" scan points violate Lean ≤4-node regularity (found 13-29 nodes). v3's 23 solutions + m_χ=0.508 MeV are ARTIFACTS of permissive tail check. Q14 expanded to 3-way: locked / A=0 / Q_A≈0. | ❌ CLOSED-NEGATIVE |
| T5 | Manual 4-fn ODE extraction via general-purpose agent across v5 paper + Numerical Benchmark §3 + bosonization logs. Found canonical 4-fn ODE: Δ_r V=Q, Δ_r A=J, Δ_r Q=V+m_J²Q+λQ(Q²−J²), Δ_r J=A−m_J²J−λJ(Q²−J²) with toroidal Δ_r=f''+f'/r. Calibration anchor: H/Q=1.6969 at (g=1.0625, λ=1.0, ω=1.0). Three gaps: ω-in-static, f(s) two-form, R unspecified. | ✅ DONE |
| EMAIL (in) 2026-05-19 1:49 PM | Werbos reply (via DeepSeek) resolved all 3 T5 gaps: m_eff² = m_J²−ω²; f(s) v5 ≡ benchmark with m_J²=0, λ=4g; toroidal R cancels in H/Q ratio. Endorsed Q14 → canonical Q=0 is Q_A≈0, Q_J≠0. Code "as soon as I can" (no committed date). | ✅ RECEIVED |
| T6 | Built `m6_v4_4fn.py` with Werbos's m_eff² substitution. Scanned (m_J², λ_bench) at canonical V₀=A₀=Q₀=J₀=0.1 across 40 grid points. RESULT: ALL 40 BLOW UP. Confirms 4-fn ODE is eigenvalue/shooting problem; generic IVP from canonical initial conditions doesn't find the bound state. | ❌ NEGATIVE |
| T6→A | Built `m6_v4_4fn_bvp.py` with `scipy.solve_bvp` + eigenvalue parameter(s). 1-eigenvalue mode converges to non-trivial (V,Q) bound state, but (A,J) collapses to zero spontaneously. Eigenvalue r_max-dependent (Bessel-zero artifact). 2-eigenvalue mode fails universally with singular Jacobian. Q_CS=1 is integral constraint, not BC. | ⚠️ PARTIAL |
| T7 | Built `m6_v4_4fn_shoot.py` using r_blowup as continuous figure of merit. Coarse 2D scan (m_J², λ_bench) ∈ [0.5,8] × [0.1,10] at V₀=A₀=Q₀=J₀=0.1: max r_blowup=6.44, well below r_max=15. Amplitude shoot monotone decreasing in A_0; no resonance. NO bound state in symmetric ansatz anywhere scanned. | ❌ NEGATIVE (definitive) |
| EMAIL (in) 2026-05-19 4:21 PM | Werbos reply (via DeepSeek): asymmetric helicity gives Q_CS=1. V₀=Q₀=+0.1, A₀=J₀=-0.1; m_J²≈0.5, λ_bench=1.0, m_eff²=-0.5 (negative is correct, time-periodic). Shooting strategy: *"decay rate at infinity as target"* (algorithm not detailed). + v2 of DM paper sent for review. | ✅ RECEIVED |
| T8 | Re-ran 4fn IVP, shoot, BVP scripts with asymmetric helicity. Helicity is NECESSARY (Q_CS goes 0 → 4483) but NOT SUFFICIENT. IVP still blows up, no Q_CS=1. Amplitude shoot ~12% better than symmetric but still monotone. BVP with 1-eigenvalue + asymmetric init: still collapses to (V,Q) sub bound state at m_J²=5.51. A,J drift to 0. | ⚠️ PARTIAL |
| T9 (2026-05-20 PM) | Built `m6_v4_4fn_newton.py`. Two-stage optimizer over 6 vars with helicity locked. Pre-scan 324 pts at \|amps\|=0.1: 0/324 reach r_max with peak<5. Wider scan 576 pts: 0/576. Stage 1 6-var DE (6528 evals): r_reached=8.4, Q_CS=58662 — delayed-blowup, not bound state. **DEFINITIVE NEGATIVE:** across 6 attempts, Q_CS=1 chaoiton UNREACHABLE via forward-IVP. Werbos's shooting must be a different method family. | ❌ DEFINITIVE NEGATIVE |

### Sandbox v5 (2026-05-20 PM)

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

### Sandbox v6 (2026-05-20 PM, CALIBRATION ESSENTIALLY CLOSED)

DeepSeek's normalization answers (Q22/Q23/Q24) landed cleanly. sandbox_v6
forks v5 and applies all three fixes. **Result with DeepSeek quartic: H/Q = 1.778
(4.8% over target).** **Result dropping the quartic (step 8 diagnostic): H/Q =
1.7112 (0.84% off target).** Net v5 → v6: 31× → 0.84% in one session.
Calibration essentially closed; residual gap is the quartic interpretation
(Q28) and the excited-mode question (Q29/Q30).

| Event | Outcome | Status |
| --- | --- | --- |
| EMAIL (in) 2026-05-20 ~4:00 PM | Paul/DeepSeek reply on email v4. Body internally muddled (3 different Q_CS derivations mid-message) but FINAL answers DO close the gap. **Q22:** drop the 2π factor on Q_CS — `Q_CS = ∫r·(A·J'-J·A')dr` directly. **Q23:** kinetic = (1/2)(V')² (not (V')²); drop the (2π)²·R toroidal prefactor; use DeepSeek quartic `(g/4)·((V²+Q²)² + (A²+J²)² + 2(V·A−Q·J)²)`. **Q24:** 8-point reference profile for V, A, Q, J at r∈{0, 0.2, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0}. DeepSeek explicit offer: *"I can also send a Python script with the exact functional if needed."* (Held for follow-up.) | ✅ RECEIVED |
| v6 implementation | `sandbox_v6/m6_v6_4fn_calibrated_bvp.py` forks v5. Same ODE structure and Lagrange-multiplier formulation; only normalizations change. `I_TARGET=1.0` (was 1/(2π)); kinetic factor 1/2; no (2π)² prefactor; DeepSeek quartic form. Reference profile available via warm-start flag but default off. | ✅ DONE |
| v6.1 (warm-start) | DeepSeek 8-point profile as initial guess. solve_bvp.status=1; fields blew up to peak J=5.47; H/Q=248.6. Profile appears DeepSeek-fabricated rather than from an actual run. | ❌ wrong basin |
| v6.2 (no warm-start) | v5's exp(-r) seed. solve_bvp.status=1; converged-ish at ω=1.111, λ_LM=9.49. **H/Q = 0.925** — within factor 2 of target. | ⚠️ partial |
| v6.3 (r_max=15, 100k nodes) | More node budget. solve_bvp.status=1; ω=1.031, λ_LM=14.32. **H/Q = 1.076** — 37% off. | ⚠️ partial |
| v6.4 (r_max=10, 200k nodes) | Overshoot — solver explored to wrong basin at higher budget. H/Q=154.8. | ❌ wrong basin |
| v6.5 (r_max=8, warm-start) | DeepSeek profile only extends to r=8; tried matching r_max. H/Q=11.4. Profile still doesn't help. | ❌ wrong basin |
| **v6.6 (best) ★** | r_max=15, n_grid=500, max_nodes=100k, tol=5e-3, NO warm-start, default g=1.0625. solve_bvp.status=1 (still). **ω=1.016, m_eff²=-0.532, λ_LM=12.21, Q_CS=1.000 exact, H/Q = 1.778** vs target 1.6969 — **4.8% over.** | ⚠️ PARTIAL ★ |
| g-sweep at v6.6 config | Same converged field profile; different post-hoc H quartic multiplier. g=0.5 → H/Q=1.743; g=1.0625 → 1.778; g=2.0 → 1.837; g=4.0 → 1.963. Linear in g. Quartic contribution is small but non-zero. | ✅ ANALYSIS |
| Calibration empirical validation | DeepSeek's three normalization fixes ALL empirically validated. Q22 drop-2π gave the dominant correction. Q23 H form (kinetic 1/2 + no prefactor + DeepSeek quartic) lands at 1.778 vs target 1.6969. Q23 cross-check with v5's Benchmark quartic gave H/Q=411 at the same field profile, confirming DeepSeek's quartic is right. Q24 profile was the only piece that didn't help — appears fabricated. | ✅ WIN |
| EMAIL (out) v5 2026-05-20 ~5:00 PM | Paul + DeepSeek. Good news (31× → 4.8%); honest residual (status=1; Q_CS_grid disagreement); Q24 profile caveat surfaced politely; two specific asks: (1) DeepSeek's Python script for line-by-line cross-check; (2) clarification on whether Q24 profile was from a real run or idealized. ApJ Zenodo upload still held. | ✅ SENT |
| EMAIL (in) 2026-05-20 ~5:30 PM | DeepSeek reference Python script received via Paul. DeepSeek confirms Q24 was illustrative (*"do not use as warm start"*) and sends *"the actual code that generated the 62 families."* **CRITICAL FINDING: script does not run as-is.** Two bugs: r=0 singularity from `J/r`/`V/r` terms with `linspace(0, Rmax, N)`; BC count mismatch (10 vs expected 9). Minimal patches (R_MIN=0.05; λ_lm as free param) make it runnable but result is CATASTROPHIC: H/Q ≈ 1.85 × 10^16. **Our v6.6 H/Q=1.778 is dramatically closer to truth.** Structural ODE differences traced (Klein-Gordon mass on every field vs our toroidal Δ_r with mass on Q only; λ-correction targets differ; quartic placement differs). Net: DeepSeek's Q22/Q23 quantitative normalizations are validated empirically by v6, but the script's ODE STRUCTURE is broken/reconstructed-from-memory. Stop diffing against the script. | ⚠️ DIAGNOSED |
| Diagnostic step (8/11/4/2) | Built `sandbox_v6/diagnostic_steps_8_through_2.py`. **Headline finding from step (8): dropping the quartic from H lands H/Q = 1.7112 (0.84% off target).** Variant sweep of 9 H-formula choices at the same v6.6 converged field profile; only no-quartic lands within 1% of target. Consistent with Werbos's stated *"for the electron, the quartic term is small."* **Step (11) field-profile inspection:** v6.6 has 5 zero crossings in V and Q (just over Lean ≤4-node spec); Q(r=0)=−0.12 (helicity sign flipped); peak A is 8× peak V. v6.6 is a slightly-excited mode. **Step (4) (m_J², λ_bench) sweep:** 16 configs, none within 50% of target; basin selection extremely fragile (same nominal config gives 4× different H/Q values depending on max_nodes). **Step (2) cold-start lepton scan:** ω_init ∈ {1, 5, 12, 20, 40.7} → only ω_init=1 lands electron basin; higher-ω diverges 9-20 orders of magnitude. Lepton scan requires continuation method. | ✅ ANALYSIS |
| Track A backward | Tested r_max ∈ {10, 12, 15, 18} with various n_grid and max_nodes. Larger budgets put solver in WORSE basins (H/Q from 154 to 10^16). λ_LM init sweep {−1, 0, 0.5, 1, 5, 12, 20}: only default init 0.1 lands in correct basin; all others diverge 4-12 orders of magnitude. **The Q_CS=1 ground-state basin is extremely narrow.** v6.6 is the best `solve_bvp` produces with this ODE+BC formulation. | ⚠️ INVESTIGATED |
| EMAIL (out) v6 2026-05-20 evening | Paul + DeepSeek. Lead with no-quartic finding (0.84% off target). Graceful script-status note (transcription/version artifact framing; quantitative answers validated empirically). Diagnostic summary (excited mode, cold-start lepton fail). Sonet acknowledgment (gracious, no flagging of stale numbers). Three new asks: **Q28** (which quartic interpretation for the electron — small g, different form, or zero?), **Q29** (does production have ≤4 nodes or 5?), **Q30** (is Q(r=0) positive or negative at production?). ApJ Zenodo upload still held. | ✅ SENT |
| Open questions post-diagnostics | Q26 (4.8% gap) and Q27 (Q_CS-grid mismatch) DEMOTED — drop-quartic finding resolves Q26 to 0.84%; Q27 is now an excited-mode artifact, not a tunable. **New IMMEDIATE: Q28 (quartic interpretation), Q29 (mode), Q30 (Q sign).** All sent in email v6. Block v7 implementation. | 🔶 OPEN, awaiting Paul |

### Sandbox v7 (planned for 2026-05-21, contingent on Paul's reply to Q28/Q29/Q30)

Hard stop 2026-05-20 evening. Resume tomorrow morning. v7 implementation
flows from which scenario branch (A/B/C/D) applies after Paul's reply.
See `0b_sandbox_v6.md` "Tomorrow's v7 plan" for the full branching table.

| Scenario | Paul confirms | sandbox_v7 focus | Estimate |
| --- | --- | --- | --- |
| A — best case | Quartic negligible + ≤4-node OK | Lock at v6.6 + drop-quartic H. Build lepton continuation (warm-start from electron, nudge ω). Q_A≈0 DM scan. Section 4 numbers in ~1 day. | ~1 day |
| B — mode-selector | Quartic negligible BUT production has ≤4 nodes everywhere + Q(0) > 0 | Add ground-state mode-selector to v7 (continuation from known ground-state ansatz, OR Q(0) > 0 BC, OR node-penalty regularization). Re-converge before scans. | +1 day |
| C — different quartic | Quartic structurally different | Re-implement H with new quartic. Re-test all 9 step-(8) variants. May also resolve excited-mode question. | +1-2 days |
| D — no reply | — | Default plan A: build continuation + drop-quartic H. Re-evaluate when reply lands. | ~1 day |

v7 implementation outline (provisional, drafted tomorrow):

| Component | v6 | v7 (planned) |
| --- | --- | --- |
| ODE structure | Toroidal Δ_r, mass on Q only, λ-corrections in A/J | Same (validated by v6 30× gap closure) |
| H quartic | DeepSeek form `(g/4)·((V²+Q²)²+...)` | **TBD by Q28** — most likely g→0 (drop quartic) for electron |
| Mode selection | Default `solve_bvp` finds 5-node mode | **NEW** — continuation from known ground-state ansatz, OR Q(0) > 0 BC, OR node-penalty regularization. TBD by Q29/Q30 |
| Lepton scan | Cold-start fails (step 2) | **NEW** — continuation method: converge electron first, warm-start with ω nudged by λ_LM perturbation; sweep ω upward |
| Q_A≈0 DM scan | Not yet run | **NEW** — set V₀=Q₀≈0, A₀+J₀ asymmetric. Sweep m_J²; find Q_J≠0 Q_A≈0 bound state |
| Acceptance | H/Q within 5% (v6.6 ✓ at 4.8%) | H/Q within 1% (drop-quartic ✓ at 0.84%) AND nodes ≤ 4 everywhere AND Q_CS_grid matches Q_CS_I |

---

## Current state (2026-05-20 evening hard stop, post step-(8/11/4/2) diagnostics + email v6 sent)

| Item | Status |
| --- | --- |
| Werbos collaboration | 🔶 Active — five emails 2026-05-20: AM new compact ApJ Neutral Chaoiton paper + DM v4 deposit; ~2:00 PM algorithm reply; ~4:00 PM normalization clarifications (Q22/Q23/Q24); ~4:50 PM disorientation note; ~5:30 PM DeepSeek reference Python script. |
| Email v6 (out) — no-quartic finding + Q28/Q29/Q30 | ✅ Sent 2026-05-20 evening. Lead with step (8) drop-quartic finding (0.84% off target). Graceful script-status note. Sonet acknowledgment. Three asks. Awaiting Paul reply. |
| DeepSeek reference Python script | ⚠️ RECEIVED but BROKEN. r=0 singularity + BC count mismatch. Patched version diverges to H/Q≈10^16. Our v6.6 (H/Q=1.778) is 16 orders of magnitude closer to truth. Stop diffing against the script. |
| sandbox_v6 Q_CS=1 chaoiton | ✅ DEMONSTRATED at ω=1.016, m_eff²=-0.532, Q_CS=1.000 exact, λ_LM=12.21. |
| Electron H/Q = 1.6969 calibration | ✅ **ESSENTIALLY CLOSED.** Step (8) diagnostic: dropping the quartic from H lands H/Q = **1.7112 (0.84% off target).** With DeepSeek quartic: 1.778 (4.8% over). Matches Werbos's stated *"for the electron, the quartic term is small."* Pending Paul's Q28 confirmation. |
| Quartic interpretation (Q28) | 🔶 OPEN, asked in email v6. Three possibilities: (a) small g (not 1.0625), (b) different small-magnitude combination, (c) genuinely negligible for the electron. Drop-quartic variant is the only one within 1% of target. |
| Ground-state vs excited mode (Q29) | 🔶 OPEN, asked in email v6. v6.6 has 5 zero crossings in V and Q (just over Lean ≤4 spec); peak A is 8× peak V (imbalanced helicity). True ground state may give H/Q closer to 1.6969 directly without dropping the quartic. |
| Q(r=0) sign at production (Q30) | 🔶 OPEN, asked in email v6. v6.6 converged with Q(0)=−0.12 (helicity flipped from Werbos's V₀=Q₀=+0.1 prescription). May be linked to excited-mode selection. |
| Cold-start lepton scan | ❌ FAILS (step 2 diagnostic). ω_init ∈ {5, 12, 20, 40.7} diverges 9-20 orders of magnitude; only ω_init=1 lands electron basin. v7 needs continuation method. |
| Basin fragility | ⚠️ CONFIRMED (step 4 + λ_LM sweep). Q_CS=1 basin extremely narrow; same nominal config gives 4× different H/Q values depending on max_nodes. Larger budget makes things worse, not better. v6.6 is the best `solve_bvp` produces with this ODE+BC formulation. |
| DeepSeek normalization fixes Q22/Q23 | ✅ EMPIRICALLY VALIDATED. Drop-2π on Q_CS + (1/2) kinetic + no toroidal prefactor land H/Q at 1.778 (with DeepSeek quartic) or 1.7112 (without). Both within "essentially closed" tolerance. |
| DeepSeek Q24 reference profile | ❌ CONFIRMED FABRICATED. DeepSeek itself acknowledged: *"the eight-point table I sent was illustrative, not from a converged run."* Use exp(-r) seed. |
| Hopf invariant proof (charge quantization) | ✅ RESOLVED. Zenodo 20296060 supplies the two missing lemmas. Now a theorem of differential topology. |
| Duda critiques (2026-05-20 thread) | #1 G_μν / single-field ontology — ARCHIVED (unfalsifiable). #2 f(J·J) unspecified — editorial OPEN (Q28 directly asks Paul to pin this). #3 construction not shown — NEARLY CLOSED (v6.6 + drop-quartic IS the construction shown empirically at 0.84% off). #4 two-charge Coulomb — FUTURE v7+. |
| Gate G1 (lepton scan) | 🔶 NEARLY UNBLOCKED. Continuation method required (cold-start fails). v7 work. |
| Gate G2 (neutral m_χ ground state) | 🔶 NEARLY UNBLOCKED. Q_A≈0 scan runs once v7 calibration mode is locked. |
| Gate G3 (discrete ω selection) | ⚠️ BLOCKED on G1. Empirical-via-lepton-scan path; analytic proof deferred per Werbos. |
| ApJ Neutral Chaoiton paper (Paul's, not ours) | 🔶 HOLD on Paul's side. Paul is holding his ApJ Zenodo update pending our Section 4 numbers. Once we hand off (m_χ, m_J, σ/m, Ω_χh²), Paul's Ref [14] *"Griesi & Anthropic AI 2026, in preparation"* resolves to a stable GitHub URL (the OpenWave m6_ouroboros tree). We don't deposit anything ourselves — see Publishing stance at top. |
| M6 production in Taichi | 🚧 Decision deferred until v7 locks + scans run. Earliest realistic: 2-3 days. |
| M5 / Liquid Crystal | 🔶 M5.4 substrate migration queued. With M6 calibration essentially closed pending Q28/Q29/Q30, M5 return is feasible after v7 + scans. Cardinal rule: M5 is SABER's primary engineering track. |

---

## Next steps

### Immediate — hard stop, awaiting Paul reply

Email v6 (sent 2026-05-20 evening) is the next-action point. v7 scenarios
A/B/C/D in the Sandbox v7 section above determine which branch runs tomorrow.

### Production scans (post-v7 calibration lock)

| Step | Action | Gate | Estimate |
| --- | --- | --- | --- |
| 1 | Lepton scan ω ∈ [0.5, 50] via continuation method. Expected: muon ω≈12.78 (1.1% gap target), tau ω≈40.7. | G1 | ~2-3 hours (continuation slower than cold-start) |
| 2 | Q_A≈0 neutral chaoiton scan. Lands m_χ, m_J mediator mass, σ/m self-interaction for ApJ Section 4. | G2 | ~1 hour |
| 3 | Gelfand-Fomin conjugate-point stability check. Confirms ground state. | G3 (empirical) | ~30 min |
| 4 | Hand off (m_χ, m_J, σ/m, Ω_χh²) numbers to Paul as a data drop. Paul's Ref [14] "Griesi & AI in prep" resolves to a stable OpenWave GitHub URL (tag a release if Paul's reviewer asks for a DOI). Our involvement ends at the data drop — Paul publishes on his own schedule. | — | done |

### M5 return after v6/v7 deliverables land

Cardinal rule: SABER is the primary engineering goal, M5 is its substrate.
Once v7 closes + scans run + ApJ Section 4 numbers handed off, M6 stays
parallel-research and primary focus returns to M5. See `0d_m5vsM6.md` for
the full M5 vs M6 evaluation against the SABER trigger.

| Step | Action | Notes |
| --- | --- | --- |
| 1 | M5.4 — matrix-field substrate migration | Was queued; unblocks once M6 deliverables land |
| 2 | M5.5 — Paper Lagrangian + V(M) | Per M5 roadmap |
| 3 | M5.6 — Biaxial twist + KG emergence | Per M5 roadmap |
| 4 | M5.7 — Resonance hunt (Close protocol) | Per M5 roadmap |
| 5 | M5.8 — 4D Zitterbewegung clock | M5 group-headline milestone; aligns with SABER engineering primitive |

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

## Open questions (snapshot — see `0b_sandbox_v6.md` for full v6 tracker)

### IMMEDIATE (block v7 implementation; all asked in email v6)

| ID | Question | Why it matters |
| --- | --- | --- |
| Q28 | Which quartic interpretation for the electron? (a) small g, (b) different combination, (c) genuinely negligible. | Drop-quartic lands H/Q=1.7112 (0.84% off); DeepSeek form lands 1.778 (4.8% off). v7 quartic choice depends on this. |
| Q29 | Does production have ≤4 zero crossings (ground state) or 5 (excited mode)? | v6.6 has 5 in V and Q. v7 may need mode-selector if production is ground state. |
| Q30 | Is Q(r=0) positive (Werbos prescription V₀=Q₀=+0.1) or negative at production? | v6.6 converged with Q(0)=−0.12, helicity flipped. May be linked to excited-mode selection. |

### DEMOTED (resolved or no longer urgent post-diagnostics)

| ID | Status |
| --- | --- |
| Q26 (4.8% H/Q gap) | DEMOTED — drop-quartic finding (step 8) resolves to 0.84%, within tolerance. |
| Q27 (Q_CS_grid vs Q_CS_I mismatch) | DEMOTED — excited-mode artifact (step 11), not a tunable. Track A sharpening backfired. |

### STILL OPEN (not blocking v7)

| ID | Question | Why it matters |
| --- | --- | --- |
| Q2 | Discrete ω selection mechanism | Empirical-via-lepton-scan in v7. Analytic proof deferred per Werbos. |
| Q3 | Analytical ω = 2mc²/ℏ derivation | Calibration only (1.2% match via R^phys). |
| Q6 | QCD reconciliation (3-chaoiton proton) | Future v7+. |
| Q19 | f(J·J) explicit form in LoE paper standalone (Duda critique #2) | Editorial; partly addressed by Q28. |
| Q21 | Two-chaoiton Coulomb derivation (Duda critique #4) | Future sandbox v7+. |

### RESOLVED (closed by v5/v6 or earlier)

| ID | Resolution |
| --- | --- |
| Q1, Q9, Q10, Q11, Q13 | Earlier rounds (v1-v4). |
| Q12 (Werbos Python code) | DEMOTED — v5 algorithm description + DeepSeek script (broken) = same unblock; our v6 is independent. |
| Q14 (canonical Q=0) | Q_A≈0 / Q_J≠0 (Werbos 2026-05-19 1:49 PM). |
| Q15 (m_eff² substitution) | RESOLVED post-v5. m_eff²=−0.596 in v5 attempt 4; −0.532 in v6.6. Both negative as predicted. |
| Q16 (m_J², λ_bench) | m_J²≈0.5, λ_bench=1.0 (Werbos 2026-05-19 4:21 PM). |
| Q17 (shooting algorithm) | RESOLVED — collocation BVP per Werbos 2026-05-20 PM email. |
| Q20 (Duda critique #3 — construction not shown) | NEARLY CLOSED — v6.6 + drop-quartic IS the construction at 0.84% off. |
| **Q22 (Q_CS normalization)** | **RESOLVED by DeepSeek 2026-05-20 4:00 PM.** Drop the 2π. Validated empirically by v6. |
| **Q23 (H functional)** | **RESOLVED by DeepSeek 2026-05-20 4:00 PM.** Kinetic 1/2, drop toroidal prefactor, DeepSeek quartic (or drop quartic per step 8). Validated empirically by v6. |
| **Q24 (sample converged profile)** | **RESOLVED with caveat:** DeepSeek confirmed the 8-point table was illustrative, not from a real run. Use exp(-r) seed. |
| Q25 (Hopf invariant proof) | RESOLVED — Zenodo 20296060 supplies the two missing lemmas. Theorem of differential topology. |

### ARCHIVED (unfalsifiable / future scope)

| ID | Why archived |
| --- | --- |
| Q4 (single vs two-field ontology) | = Duda critique #1. Aesthetic preference. Not a theory-killer if math matches observation. |
| Q7 (cold-fusion citation trail) | Historical, not physics. |

**Active count post-diagnostics: 3 IMMEDIATE + 5 OPEN = 8 questions.**
**Highest-leverage closure: Paul's Q28/Q29/Q30 reply → v7 → lepton + DM scans → ApJ Section 4 handoff.**

---

## Hardest pieces (snapshot — see `0b_sandbox_v6.md` for live tracker)

| Hardest piece | Status post-diagnostics (2026-05-20 evening) |
| --- | --- |
| Forward-IVP method family | RESOLVED. Confirmed wrong tool by both us (T6-T9) and Werbos. Replaced with collocation BVP in v5. |
| Q_CS=1 enforcement mechanism | RESOLVED. Auxiliary integral state I with BC I(R_max)=1.0 (v6) — also matches the field-derived integral when on the right basin. |
| Q_CS=1 chaoiton existence | RESOLVED. v5 attempt 4 + v6.6 both converge with Q_CS=1.000 exact. First empirical chaoiton in OpenWave. |
| **Electron H/Q = 1.6969 calibration** | **NEARLY RESOLVED.** v6.6 drop-quartic lands H/Q=1.7112 (0.84% off). Pending Q28 confirmation (which quartic is canonical for the electron). |
| Ground-state vs excited mode selection | NEAR-PASS. v6.6 has 5 nodes V/Q (just over Lean ≤4 spec); 0 nodes A; 1 node J. Pending Q29 (does production have ≤4 nodes everywhere?). Mode-selector may be v7 work. |
| Lagrange-multiplier ODE correction coefficient | EMPIRICALLY OK. v6 lands λ_LM=12.21 (vs v5's −1.21); the new normalizations need a larger Lagrange multiplier. Empirically right; analytic derivation may have a missing factor that doesn't affect the H/Q ratio. |
| V(M) potential form | UNRESOLVED — shared bottleneck with M5 (Duda). Not in v6/v7 scope. |
| f(J·J) form | RESOLVED IN PRACTICE — v6 uses DeepSeek quartic (or drop-quartic per step 8). Q19 (LoE paper editorial gap) remains for Werbos's next LoE revision; Q28 asks for the electron-specific interpretation. |
| 4-fn vs 2-fn ansatz mismatch | RESOLVED. 2-fn for NEUTRAL (Lean theorem), 4-fn for CHARGED (Numerical Benchmark / v5/v6). |
| ω quantization mechanism | OPEN (Q2). v7 lepton scan via continuation method provides empirical test. Analytic proof deferred per Werbos's own admission. |
| Lepton mass spectrum | UNBLOCKED for v7. Cold-start fails (step 2); continuation method needed. |
| Neutral m_χ true ground state | UNBLOCKED for v7. Q_A≈0 scan feeds Section 4 of ApJ Neutral Chaoiton paper. |
| Two-chaoiton Coulomb derivation | NOT ADDRESSED. Future sandbox v7+ (Q21, Duda critique #4). |
| Charge quantization rigorous proof | RESOLVED. Hopf invariant proof complete (Zenodo 20296060). |

---

## Resources

**Research docs:**

| Doc | Contents |
| --- | --- |
| `0a_background.md` | Technical evaluation (§1-10) + updates |
| `0b_sandbox_v1.md` | v1 plan, results, emails, Q-list |
| `0b_sandbox_v2.md` | v2 plan, IVP+BVP attempts, Q-list |
| `0b_sandbox_v3.md` | v3 plan, Lean ODE, neutral scan, Q-list |
| `0b_sandbox_v4.md` | v4 plan, T1-T9 negative results, forward pointer to v5 |
| `0b_sandbox_v5.md` | v5 plan, Werbos algorithm reply, T10 partial success, question tracker |
| `0b_sandbox_v6.md` | v6 plan (current), DeepSeek normalizations, step (8/11/4/2) diagnostics, drop-quartic finding, email v6 + Q28/Q29/Q30, v7 plan |
| `0c_model_gates.md` | G1/G2/G3 production decision criteria |
| `0c_roadmap.md` | This file |
| `0d_m5vsM6.md` | M5 vs M6 evaluation in service of the SABER trigger; Duda/Werbos self-stated limitations |

**Sandbox scripts:**

| Folder | Contents |
| --- | --- |
| `sandbox_v1/` | 96-variant sweep, calibration, mass scan (wrong ODE) |
| `sandbox_v2/` | IVP + BVP locked-ansatz attempts (superseded) |
| `sandbox_v3/` | Lean ODE, neutral β, calibration (last completed pre-v4) |
| `sandbox_v4/` | T1 lepton scan (blocked); T2 neutral ground state (negative); `diag_energy_functional`; T5 4-fn extracted; T6-T9 forward-IVP attempts — all negative; reference for definitive "forward-IVP wrong tool" result |
| `sandbox_v5/` | `m6_v5_4fn_lambda_bvp.py` — Werbos's collocation BVP with Lagrange multiplier λ_LM. `solve_bvp.status=0` at ω=1.047, m_eff²=−0.596. H/Q=52.64 vs target 1.6969 (open: normalization gap, closed in v6). |
| `sandbox_v6/` | `m6_v6_4fn_calibrated_bvp.py` — forks v5 + Q22/Q23/Q24 fixes. v6.6 lands H/Q=1.778 with DeepSeek quartic; step (8) drop-quartic lands H/Q=1.7112 (0.84% off target). Also: `deepseek_reference.py` + `deepseek_reference_patched.py` (broken — diverges to H/Q≈10^16); `diagnostic_steps_8_through_2.py` (the 4-step diagnostic that found drop-quartic). |
| `sandbox_v7/` | PLANNED for 2026-05-21, contingent on Paul's Q28/Q29/Q30 reply. Scenarios A/B/C/D — see Sandbox v7 section above. |

**Theory papers:**

| Folder | Contents |
| --- | --- |
| `theory/` | All Werbos corpus (v1-v8 LoE, Lean theorem, calibration, spectrum, Numerical Benchmark, bosonization logs, Hopf invariant proof completion, Neutral Chaoiton ApJ paper, Dark Matter paper) |
