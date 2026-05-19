# M6 / Ouroboros — Roadmap

**Status:** 🔶 Active sandbox evaluation. Production decision pending gate G1 + G2.

Last updated: 2026-05-18.

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

### Phase 4 plan — Sandbox v4 (starting 2026-05-19)

```text
Track | Action                              | Blocked by
------|-------------------------------------|------------------
T1    | Lepton scan with Lean 2-fn ODE      | Not blocked
      | → test "lowest 3 stable = leptons"  |
T2    | Neutral ground state via B0 shoot   | Not blocked
T3    | Werbos's Python source (1-2 days)   | Werbos reply
T4    | Reconcile Q14 (locked vs A=0)       | T3 (code)
```

---

## Current state (2026-05-18, evening)

```text
Item                                | Status
------------------------------------|------------------------------------------
Werbos collaboration                | 🔶 Active — reply received 19:39
                                    |    Code coming in 1-2 days
Email sent (sandbox v3 results)     | ✅ Sent 2026-05-18 (edited by Rodrigo)
Werbos reply received               | ✅ 19:39 via DeepSeek summary
Gate G1 (lepton scan, Lean ODE)     | 🔶 READY TO RUN — v4 Track 1
                                    |    Tests Werbos's "lowest 3 stable=leptons"
                                    |    hypothesis directly
Gate G2 (neutral m_χ ground state)  | 🔶 CODE PROMISED — ETA 2026-05-20
Gate G3 (discrete ω selection)      | ⚠️  PATH FORWARD — analytic proof
                                    |    deferred; empirical G1 result is the
                                    |    effective gate
M6 production in Taichi             | 🚧 On hold pending G1 + G2
                                    |    Production decision plausible by 2026-05-25
M5 / Liquid Crystal                 | 🔶 Active — M5.4 substrate migration queued
```

---

## Next steps

### Immediate (waiting on Werbos's code)

```text
Step | Action                             | Gate | ETA
-----|------------------------------------|------|------------------------
1    | Receive Werbos's Python source     | G2   | 2026-05-20 (1-2 days
     | (Q=0 scan + calibration code)      |      | per his reply)
     | → close calibration localization   |      |
     | → verify neutral m_χ ground state  |      |
     | → resolve Q14 (locked vs A=0)      |      |
```

### Can run now (not waiting on Werbos)

```text
Step | Action                              | Gate | ETA
-----|-------------------------------------|------|------------------------
2    | Run v4 Track 1: lepton scan with    | G1   | 1-2 days
     | Lean 2-fn ODE                       |      |
     | → test "lowest 3 stable modes =     |      |
     |   leptons" hypothesis (Werbos)      |      |
     | → if ω={1, ~12.78, ~40.7} match     |      |
     |   leptons within 4-6%, G3 is        |      |
     |   empirically closed alongside G1   |      |
3    | Run v4 Track 2: neutral ground      |      | 1 day
     | state via B0 shooting               |      |
     | → tune B0 to minimum-energy β       |      |
     | → confirm m_χ value (cross-check    |      |
     |   with Werbos's code when arrives)  |      |
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
Q2   | Discrete ω selection mechanism       | PATH FORWARD — Werbos: leptons =
                                            | lowest 3 stable modes; v4 T1 tests
Q3   | Analytical ω = 2mc²/ℏ derivation     | PARTIAL — calibration only
Q6   | QCD reconciliation for proton        | OPEN — uninvestigated
Q8   | Neutral m_χ true ground state        | PARTIAL — v4 T2 will refine
Q12  | Werbos's Python source code          | PROMISED — ETA 2026-05-20
Q13  | 4-function vs 2-function reduction   | PARTIAL — DeepSeek answer received;
                                            | code will fully resolve
Q14  | Canonical Q=0: locked or A=0?        | NEW — surfaced 2026-05-18 19:39;
                                            | code (T3) will resolve
```

---

## Hardest pieces (summary — see 0b_sandbox_v4.md for the live tracker)

```text
Hardest piece       | Status at start of v4
--------------------|--------------------------------------------------
V(M) potential form | UNRESOLVED — shared bottleneck with M5 (Duda)
f(J·J) form         | MOSTLY SETTLED — f=gs² confirmed (v5 + Lean);
                    | factor-of-2 vs factor-of-4 check still needed
ω quantization      | OPEN — Werbos hypothesis: leptons = 3 lowest
mechanism           | stable Q=1 modes. v4 Track 1 tests directly.
Q=0 chaoiton        | PARTIAL — 23 solutions; ground state pending
existence           | (Track 2); Q14 locked-vs-A=0 conflict pending
                    | (Track 3 / code)
Lepton mass         | UNTESTED with Lean ODE — Track 1 will settle
spectrum            | in v4
Electron calibration| OPEN — H/Q = 1.723 (1.56% gap) but not yet
localization        | localized in v3. Track 3 (Werbos's energy
                    | functional) is the primary fix; T1 baseline
                    | depends on it.
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
  sandbox_v4/          Lepton scan (T1) + neutral ground state (T2) — planned

Theory papers:
  theory/              All Werbos corpus (v1-v5, Lean, calibration, spectrum,
                       bosonization logs May 18 4PM, 11 new files 17:20)
```
