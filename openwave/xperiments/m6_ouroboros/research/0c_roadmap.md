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

```text
Task                                    | Outcome
----------------------------------------|----------------------------------------
Read all available Ouroboros papers     | ✅ Full corpus in theory/
Evaluated vs Duda's LdGS challenges     | ✅ 0a_background.md §1-8
Side-by-side M5 vs M6 comparison        | ✅ 0a_background.md §3
Verdict: viable M6 candidate?           | ✅ YES — with caveats (0c_model_gates.md)
```

### Phase 1 — Sandbox v1 (2026-05-17)

```text
Task                                    | Outcome
----------------------------------------|----------------------------------------
Reproduce electron calibration          | ✅ H/Q = 1.6918 (0.30% gap)
  at (g=1.0625, ω=1.0)                  |    → ODE form A4+H1+Q2 calibrates
Reproduce lepton mass scan ω^2.22       | ❌ Got ω^2.04, lepton gaps 31-44%
Sent findings + THE QUESTION to Werbos  | ✅ Email 2026-05-17 (in 0b_sandbox_v1.md)
Key finding: ω^2.22 = log(207)/log(11)  | ✅ Post-hoc fit signature detected
  → ω values appear fitted, not predicted|
```

### Phase 2 — Sandbox v2 (2026-05-18)

```text
Task                                    | Outcome
----------------------------------------|----------------------------------------
Werbos replied — active collaboration   | ✅ 2026-05-17 "most serious feedback"
Read benchmark doc (new 4-fn ODE)       | ✅ But benchmark doc had wrong form
Attempt Q=0 neutral chaoiton via IVP    | ❌ Bessel J_0 oscillation — no decay
  locked-ansatz (J=-A, Q=-V)            |
Attempt Q=0 neutral chaoiton via BVP    | ❌ Only trivial + r_max artifacts
  (scipy solve_bvp)                     |    No physical chaoiton found
Key finding: locked-ansatz approach     | ✅ Wrong approach identified. Both IVP
  is fundamentally incorrect            |    and BVP confirmed this.
```

### Phase 3 — Sandbox v3 (2026-05-18)

```text
Task                                    | Outcome
----------------------------------------|----------------------------------------
Read 11 new files from Werbos (17:20)   | ✅ Key finds: Lean ODE + neutral method
Found canonical ODE in Lean theorem     | ✅ Vector Laplacian + slope BCs
  (chaoiton_theorem.lean.txt)           |    (fields start at zero, not A₀)
Found neutral chaoiton method           | ✅ α=0, ω=0, single β ODE
  (Claude bosonization May18 4PM.txt)   |    (A-field off, spin-0)
Rebuilt ODE with Lean canonical form    | ✅ m6_v3_chaoiton.py
  (vector Laplacian, slope BCs)         |
Neutral chaoiton scan (α=0, ω=0)        | ✅ 23 localized solutions
                                        |    m_χ = 0.508 MeV at λ=1.0
                                        |    (vs Werbos est. 0.003-0.015 MeV)
Electron calibration with Lean ODE      | ⚠️  H/Q = 1.723 (1.56% gap)
                                        |    Not yet localized — energy
                                        |    functional needs geometry factor
Sent v3 results + Q12/Q13 to Werbos     | ✅ Email 2026-05-18 (paste-ready draft)
Code search on Zenodo (3 records)       | ❌ No Python code found anywhere
                                        |    DOCX files only — code is "available
                                        |    on request from QAGI LLC"
Werbos reply (via DeepSeek, 19:39)      | ✅ Q12 PROMISED (code in 1-2 days)
                                        |    Q13 partial answer (locked V→Q, A→J)
                                        |    Q2 hypothesis: leptons = lowest 3
                                        |       stable Q=1 modes
                                        |    NEW conflict surfaced: locked vs A=0
                                        |       descriptions of neutral chaoiton
                                        |       → new Q14
```

### Phase 4 plan — Sandbox v4 (planned 2026-05-19)

```text
Track | Action                              | Blocked by
------|-------------------------------------|------------------
T1    | Lepton scan with Lean 2-fn ODE      | Not blocked
      | → test "lowest 3 stable = leptons"  |
T2    | Neutral ground state via B0 shoot   | Not blocked
T3    | Werbos's Python source (1-2 days)   | Werbos reply
T4    | Reconcile Q14 (locked vs A=0)       | T3 (code)
```

### Phase 5 — Sandbox v4 first runs (2026-05-19)

```text
Track | Outcome                                                | Status
------|--------------------------------------------------------|---------
T1    | Built m6_v4_lepton_scan.py, ran ω∈{1,5,..,45} quick    | ❌ BLOCKED
      | scan + diag_energy_functional.py for H-scaling test.   |
      | Found: 2-fn Lean ODE is WRONG TOOL for charged sector. |
      | H decreases with ω (opposite Werbos's bosonization);   |
      | ω=1.0 marginally fails localization; m_μ/m_e ratio     |
      | target 207 found 0.04. v3's H_CODE_ELECTRON_CALIB=0.494|
      | also confirmed wrong. Per Q13: 2-fn for NEUTRAL, 4-fn  |
      | (φ,α,ρ,β) for CHARGED. Bosonization muon@1.1% used     |
      | 4-fn code, not Lean reduction.                         |
T2    | Built m6_v4_t2_neutral_ground.py. Wide log-scan B₀ ∈   | ❌ CLOSED -
      | [1e-5, 1.0] × 60 points × 6 λ values; golden-section   |   NEGATIVE
      | refinement. NEGATIVE RESULT: H is monotonic in B₀;     |
      | no minimum exists. All "localized" scan points violate |
      | Lean ≤4-node regularity (found 13-29 nodes). Linear    |
      | part β''+β'/r-β/r²+λβ=0 is J_1(√λr) Bessel —           |
      | oscillatory not bound. v3's 23 solutions + m_χ=0.508   |
      | MeV are ARTIFACTS of permissive absolute-tail check.   |
      | Q=0 ground-state existence now uncertain. Q14 expanded |
      | to 3-way: locked / A=0 / Q_A≈0.                        |
T5    | Manual 4-fn ODE extraction via general-purpose agent   | ✅ DONE
(new) | across v5 paper + Numerical Benchmark §3 + bosonization|
      | logs. Found canonical 4-fn ODE: Δ_r V=Q, Δ_r A=J,      |
      | Δ_r Q=V+m_J²Q+λQ(Q²−J²), Δ_r J=A−m_J²J−λJ(Q²−J²) with |
      | toroidal Δ_r=f''+f'/r (NO -f/r²), value BC f'(0)=0,    |
      | toroidal volume (2π)²R·r dr. Calibration anchor:       |
      | H/Q=1.6969 at (g=1.0625, λ=1.0, ω=1.0, V₀=A₀=Q₀=J₀=0.1)|
      | Confidence: medium. Three gaps: ω-in-static, f(s)      |
      | two-form, R unspecified.                               |
NEW   | Werbos 2026-05-19 10:51 AM sent DM paper draft for ApJ | ✅ REPLY
      | review citing Griesi v3 m_χ=0.508 MeV as load-bearing  | RECEIVED
      | numerical input. Rodrigo replied 11:30ish with T2      |
      | invalidation + 4-fn extraction + honest code ask.      |
      | Werbos reply 2026-05-19 1:49 PM (via DeepSeek):        |
      | resolved all 3 implementation gaps (m_eff² formula,    |
      | f(s) two-form mapping, R cancels in H/Q). Endorsed     |
      | our Q14 → canonical Q=0 is Q_A≈0, Q_J≠0. Code "as soon |
      | as I can" (no committed date).                         |
T6    | Built m6_v4_4fn.py with Werbos's m_eff² substitution.  | ❌ NEGATIVE
(new) | Tried RK45, LSODA, RK45+max_step+blowup_event. Scanned |
      | (m_J², λ_bench) at canonical V₀=A₀=Q₀=J₀=0.1 across    |
      | 40 grid points. RESULT: ALL 40 BLOW UP. Confirms       |
      | 4-fn ODE is eigenvalue/shooting problem; generic IVP   |
      | from canonical initial conditions doesn't find the     |
      | bound state — that requires the SPECIFIC calibrated    |
      | (m_J², λ_bench) which we don't yet have.               |
T6→A  | Built m6_v4_4fn_bvp.py with scipy.solve_bvp +          | ⚠️ PARTIAL
(new) | eigenvalue parameter(s). 1-eigenvalue mode converges   |
      | to non-trivial (V,Q) bound state, but (A,J) collapses  |
      | to zero spontaneously. Eigenvalue r_max-dependent      |
      | (Bessel-zero artifact). 2-eigenvalue mode fails        |
      | universally with singular Jacobian. Q_CS=1 is a        |
      | topological/integral constraint not a BVP-compatible   |
      | boundary condition.                                    |
T7    | Built m6_v4_4fn_shoot.py using r_blowup distance as    | ❌ NEGATIVE
(new) | continuous figure of merit. Coarse 2D scan (m_J²,     |    (definitive)
      | λ_bench) ∈ [0.5,8] × [0.1,10] at V₀=A₀=Q₀=J₀=0.1:     |
      | max r_blowup=6.44, well below r_max=15. m_J²=1.0 row   |
      | flat (m_eff²=0 degenerate). Amplitude shoot: r_blowup  |
      | monotone decreasing in A_0; no resonance. Conclusion:  |
      | NO bound state in symmetric V₀=A₀=Q₀=J₀ regime         |
      | anywhere we scanned. Werbos's calibration must use     |
      | asymmetric initial values OR different parameters.     |
EMAIL | Email sent ~5 PM 2026-05-19: T6/T6→A/T7 findings + ask | ✅ SENT
(out) | for specific m_J², λ_bench, V₀, A₀, Q₀, J₀ values      |
EMAIL | Werbos reply ~4:21 PM (via DeepSeek): m_J²≈0.5,        | ✅ RECEIVED
(in)  | λ_bench=1.0, V₀=Q₀=+0.1, A₀=J₀=-0.1 (asymmetric        |
      | helicity gives Q_CS=1). m_eff²=-0.5 correct. Shooting  |
      | strategy: "decay rate at infinity as target" (algorithm|
      | not specified). + v2 of DM paper sent for review.      |
T8    | Re-ran 4fn IVP, shoot, BVP scripts with Werbos's       | ⚠️ PARTIAL
(new) | asymmetric helicity. FINDING: helicity is necessary    |
      | (Q_CS goes from 0 → non-zero) but NOT SUFFICIENT.      |
      | IVP still blows up (Q_CS=4483, not 1). Amplitude       |
      | shoot: ~12% better than symmetric, but still monotone. |
      | BVP with 1-eigenvalue + asymmetric init guess: still   |
      | collapses to (V,Q) sub bound state at m_J²=5.51        |
      | regardless of init. A,J drift to 0 without a BC        |
      | forcing them. 2-eigenvalue BVP: singular Jacobian.     |
T9    | 2026-05-20: Built sandbox_v4/m6_v4_4fn_newton.py.      | ❌ NEGATIVE
(new) | Two-stage optimizer over 6 vars (m_J², λ, |V₀|, |A₀|,  |   DEFINITIVE
      | |Q₀|, |J₀|) with helicity signs locked. Stage 1 DE     |
      | global + Stage 2 Nelder-Mead local. Pre-scan: 2D grid  |
      | at Werbos's exact |amps|=0.1 (324 pts, m_J²∈[0.05,10], |
      | λ∈[0.01,10], r_max=15): 0/324 reach r_max with peak<5. |
      | Max r_reached=7.24, peak=95. Wider scan (576 pts up to |
      | m_J²=50, λ down to 0.001, r_max=30): 0/576, max        |
      | r_reached=8.06. Stage 1 6-var DE (6528 evals):         |
      | r_reached=8.4, Q_CS=58662 — delayed-blowup regime, not |
      | a bound state. SIX independent attempts (T6/T6→A 1-eig |
      | /T6→A 2-eig/T7/T8/T9) all confirm: Q_CS=1 chaoiton     |
      | UNREACHABLE via forward-IVP from value-BC origin in    |
      | any parameter region we have explored.                 |
TRIG  | 2026-05-20 PM: Paul email — new compact ApJ paper      | ✅ NOTED
(new) | "The Neutral Chaoiton: A Dark Matter Candidate from    |
      | the Ouroboros Lagrangian" arrived. Distinct from v4 DM |
      | (Zenodo 20298669); cites both v8 LoE (Zenodo 20313063) |
      | and v4 DM. Acknowledges Griesi for helicity-structure  |
      | dialogue. Reference [14] = "Griesi & AI in prep."      |
      | Paul holds "TRULY new" Zenodo upload pending OUR T9 /  |
      | ground-state numbers. With T9 negative, the m_χ + m_J  |
      | + σ/m numbers in his Section 4 depend on first         |
      | resolving the shooting-algorithm question.             |
NEXT  | Email Paul v3: definitive 6-attempt negative + sharp   | 🔶 DRAFTING
      | algorithmic question. Likely candidates to ask about:  |
      | (a) backward integration from K_0 asymptotic at r→∞    |
      | (b) BVP with Q_CS Lagrange multiplier                  |
      | (c) tabulated init profile shape, not just BC values   |
      | (d) multi-shot connection at intermediate r            |
      | Hold any Zenodo upload of either paper until algorithm |
      | clarified.                                             |
```

---

## Current state (2026-05-19, afternoon)

```text
Item                                | Status
------------------------------------|------------------------------------------
Werbos collaboration                | 🔶 Active — TWO emails today.
                                    |    AM: ApJ paper draft for review.
                                    |    PM: 3-gap clarifications + Q14
                                    |    resolution + ApJ-revision intent.
Email out (v4 findings + code ask)  | ✅ Sent 2026-05-19 AM with T2
                                    |    invalidation + 4-fn surface area
Werbos reply (3 clarifications)     | ✅ Received 2026-05-19 1:49 PM.
                                    |    m_eff²=m_J²-ω², f(s) mapping,
                                    |    R cancels. Q14 → Q_A≈0/Q_J≠0.
Werbos's promised code              | ⚠️  "as soon as I can" — no
                                    |    committed date in PM reply
T6 — 4-fn IVP attempt               | ❌ NEGATIVE — all 40 (m_J², λ) scan
                                    |    points blow up. Confirms eigenvalue
                                    |    structure; can't anchor with IVP
                                    |    from canonical initial conditions.
Gate G1 (lepton scan)               | ❌ BLOCKED — 4-fn implementation
                                    |    pending T6→A (BVP) and T7 (shoot)
Gate G2 (neutral m_χ ground state)  | ❌ v3 INVALIDATED. Werbos endorsed
                                    |    Q_A≈0 canonical; still need 4-fn
                                    |    anchor before running it.
Gate G3 (discrete ω selection)      | ⚠️  Empirical-via-G1 path still
                                    |    blocked (waiting for G1 path).
M6 production in Taichi             | 🚧 Production decision delayed past
                                    |    2026-05-27 if T6→A and T7 land
                                    |    this week; or 2026-05-30+ if
                                    |    Werbos's code is the unblock
M5 / Liquid Crystal                 | 🔶 Active — M5.4 substrate migration queued
```

---

## Next steps

### Immediate (waiting on Werbos's reply to 2026-05-19 ask)

```text
Step | Action                             | Gate | ETA
-----|------------------------------------|------|------------------------
1    | Receive Werbos response to T2      | G1+  | 2026-05-20/21
     | invalidation + 4-fn ODE ask        | G2   |   - code → unblocks all
     | → If code arrives: implement, run  |      |   - ODE-form only →
     |   lepton scan + Q=0 ground state   |      |     still need full
     | → If ODE-form only: implement      |      |     scan code (~2 day)
     |   benchmark 4-fn, calibrate to     |      |   - no response →
     |   H/Q=1.6969, then run scans       |      |     fall back to T5→
     | → If no response in ~24-48h: solo  |      |     code (Path A)
     |   4-fn implementation              |      |
```

### Can run now (solo path, if Werbos slips)

```text
Step | Action                              | Gate | ETA
-----|-------------------------------------|------|------------------------
2    | Solo T5→code: implement 4-fn        | G1+  | 1-2 days
     | benchmark ODE in sandbox_v4/.       | G2   |
     | Calibrate against H/Q = 1.6969      |      |
     | first; only proceed to lepton       |      |
     | scan + Q=0 scan once anchored.      |      |
     | Three gaps need empirical iter:     |      |
     |   - ω in static ODE (try m_J² ←    |      |
     |     m_J² − ω² substitution)         |      |
     |   - f(s) v5 vs benchmark form       |      |
     |   - toroidal R value (try R=1)      |      |
3    | After 4-fn anchored: run charged    | G1   | +1 day post-anchor
     | lepton scan over ω∈[0.5, 50]        |      |
4    | After 4-fn anchored: run Q=0        | G2   | +1 day post-anchor
     | scan with Q_A≈0 constraint (the     |      |
     | untested 3rd description from Q14)  |      |
```

### If G1 + G2 PASS (GO decision)

```text
Step | Action                              | ETA
-----|-------------------------------------|------------------------
4    | Scaffold M6 in Taichi               | post-M5.4
     | → Vector(4) × 2 substrate (A, J)    |
     | → Lorenz constraint enforcer        |
     | → Chern-Simons charge kernel        |
     | → Mirror M5's rendering pipeline    |
5    | Gate 1 Taichi: Maxwell limit        | Week 1 of M6 build
     | → A-field alone = standard EM       |
6    | Gate 2 Taichi: charge quantization  | Week 2
     | → Q[A,J] integer on seeded configs  |
7    | Gate 3 Taichi: chaoiton existence   | Weeks 3-6
     | → localized time-periodic solution  |
     | → same (g,λ,ω) as sandbox           |
```

### Parked (post-G1+G2)

```text
Step | Action                              | Notes
-----|-------------------------------------|------------------------
C    | 3-body proton bound state           | V(R) ~ -C/R⁶ classical
     |                                     | 3-body problem; deferred
```

---

## Open questions (summary — see 0b_sandbox_v4.md for full list)

```text
Q2   | Discrete ω selection mechanism       | BLOCKED — was empirical-via-G1;
                                            | G1 now blocked too
Q3   | Analytical ω = 2mc²/ℏ derivation     | PARTIAL — calibration only
Q6   | QCD reconciliation for proton        | OPEN — uninvestigated
Q8   | Neutral m_χ true ground state        | CLOSED-NEGATIVE — T2 showed
                                            | v3's 0.508 MeV was artifact.
                                            | Existence in 2-fn form
                                            | UNCERTAIN.
Q12  | Werbos's Python source code          | STILL PROMISED — ETA was
                                            | 2026-05-20; window open
Q13  | 4-function vs 2-function reduction   | RESOLVED in principle — 2-fn
                                            | for NEUTRAL only, 4-fn for
                                            | CHARGED. Confirmed via T5
                                            | extraction. Full reduction
                                            | mapping still unclear.
Q14  | Canonical Q=0: locked / A=0 /        | RESOLVED (Werbos 1:49 PM):
       Q_A≈0?                               | canonical DM candidate is
                                            | Q_A≈0/Q_J≠0 (both fields
                                            | active, EM-neutral). Locked
                                            | and A=0 are limiting
                                            | approximations.
Q15  | Does m_eff² substitution close       | TESTED via T6 → INCONCLUSIVE
       electron H/Q gap?                    | (eigenvalue structure means
                                            | IVP can't test from generic
                                            | initial conditions)
Q16  | What specific m_J², λ_bench values   | NEW — load-bearing for T6→A
       does Werbos use at calibration?      | / T7. Resolves if Werbos
                                            | answers OR if BVP finds it
                                            | independently.
```

---

## Hardest pieces (summary — see 0b_sandbox_v4.md for the live tracker)

```text
Hardest piece       | Status after v4 first runs (2026-05-19)
--------------------|--------------------------------------------------
V(M) potential form | UNRESOLVED — shared bottleneck with M5 (Duda)
f(J·J) form         | AMBIGUOUS — v5 says g·s²; benchmark says
                    | (m_J²/2)·s + (λ/4)·s². Reconciliation pending.
ω quantization      | BLOCKED — empirical test via G1 also blocked
mechanism           | now until 4-fn ODE implementation lands.
Q=0 chaoiton        | UNCERTAIN — v3's 23 "solutions" invalidated by
existence           | T2 (artifacts). Q14 3-way candidates: locked
                    | (v2 failed), A=0 (v3/T2 disproved), Q_A≈0
                    | (untested). Code OR canonical form resolves.
Lepton mass         | UNTESTED — 2-fn was wrong tool. Need 4-fn
spectrum            | benchmark ODE (T5 extracted; implementation
                    | open).
Electron calibration| OPEN — anchor target H/Q = 1.6969 (corrected
localization        | from v3's wrong H_CODE_ELECTRON_CALIB=0.494).
                    | Match against this once 4-fn implemented.
4-fn vs 2-fn ansatz | NEW — Numerical Benchmark and Lean theorem use
mismatch            | STRUCTURALLY different ODEs (toroidal l=0 vs
                    | vector l=1). Not reductions of each other.
                    | Werbos uses 4-fn for charged; reduction route
                    | for neutral is unclear.
```

---

## Resources

```text
Research docs:
  0a_background.md     Technical evaluation (§1-10) + updates (§11)
  0b_sandbox_v1.md     v1 plan, results, emails, Q-list
  0b_sandbox_v2.md     v2 plan, IVP+BVP attempts, Q-list
  0b_sandbox_v3.md     v3 plan, Lean ODE, neutral scan, Q-list
  0b_sandbox_v4.md     v4 plan (current), Werbos reply, tracks T1-T4
  0c_model_gates.md    G1/G2/G3 production decision criteria
  0c_roadmap.md        This file

Sandbox scripts:
  sandbox_v1/          96-variant sweep, calibration, mass scan (wrong ODE)
  sandbox_v2/          IVP + BVP locked-ansatz attempts (superseded)
  sandbox_v3/          Lean ODE, neutral β, calibration (last completed)
  sandbox_v4/          T1 lepton scan (blocked, 2-fn wrong tool),
                       T2 neutral ground state (negative result),
                       diag_energy_functional, T5 4-fn extracted

Theory papers:
  theory/              All Werbos corpus (v1-v5, Lean, calibration, spectrum,
                       bosonization logs May 18 4PM, 11 new files 17:20)
```
