# M6 / Ouroboros — Production Gates

Three gates must pass before committing M6 to full Taichi production. The
model is a credible scientific candidate; these gates determine whether it is
a credible *engineering* candidate alongside M5.

Last updated: 2026-05-19.

---

## Gate definitions

```text
Gate | Decision question                       | Unblocked by
-----|-----------------------------------------|-------------------------------
G1   | Does the lepton scan reproduce          | Run v4 ourselves (Priority A,
     | ω^2.22 scaling at <5% muon/tau          | ~1 day once Lean ODE confirmed).
     | gaps with the correct Lean ODE?         | Werbos's session shows 1.1%
     |                                         | muon gap at ω=12.78 — we need
     |                                         | to reproduce that ourselves.
G2   | What is the true neutral chaoiton       | Werbos's code (Q12 in email,
     | ground-state m_χ — our 0.508 MeV        | sent 2026-05-18). Our v3 gets
     | or Werbos's 0.003-0.015 MeV?            | 0.508 MeV; his prelim estimate
     |                                         | is 30× lighter. One of us has
     |                                         | a wrong energy functional.
G3   | Does a discrete ω selection             | Werbos reply on Q2 (email) OR
     | mechanism exist that picks              | failed: lepton masses are 3
     | ω = {1, 11, 40.7} from first            | fitted inputs, not predictions.
     | principles?                             | This is Werbos's own "key open
     |                                         | question" (Spectrum §6.1).
```

---

## Pass/fail criteria

```text
Gate | PASS criterion                          | FAIL consequence
-----|-----------------------------------------|-------------------------------
G1   | Muon gap <5%, tau gap <10%, reproduced  | M6 stays in sandbox phase
     | with our own code, not just Werbos's    | indefinitely. Cross-validation
     | bosonization session result             | value is not real if we can't
     |                                         | independently verify the claim.
G2   | m_χ confirmed to within 2× of           | DM prediction in DarkMatterv1
     | Werbos's prelim estimate with           | is at risk. Still worth building
     | same code base / energy functional      | M6 but without the DM headline.
G3   | A rigorous selection mechanism          | Lepton masses are fitted (like
     | (quantization condition) identified     | Standard Model, not a revolution).
     | OR its absence confirmed rigorously     | M6 still worth building as an
     |                                         | alternative substrate; weaker
     |                                         | scientific claim.
```

---

## Gate dependencies

```text
G2 → depends on Werbos's code (waiting reply as of 2026-05-18)
G1 → independent; can run sandbox v4 now with Lean ODE
G3 → parallel; may be resolved by Werbos reply or remain open

Recommended sequence:
  Run G1 first (days) → run G2 once code arrives (days) → G3 (weeks/open)
  Don't wait for G3 before production decision — it may never close analytically.
```

---

## GO / NO-GO decision

```text
Scenario                         | Decision
---------------------------------|------------------------------------------
G1 PASS + G2 PASS                | GO — scaffold M6 in Taichi, parallel to M5.
                                 | Start from Vector(4) × 2 substrate.
                                 | Timeline: post-M5.4 substrate migration.
G1 PASS + G2 FAIL                | CONDITIONAL GO — build M6 but without DM
                                 | headline. Lepton spectrum + nuclear force
                                 | are enough for scientific value.
G1 FAIL                          | HOLD — don't commit engineering time until
                                 | we can reproduce lepton gaps <5% ourselves.
G3 PASS (discrete mechanism)     | Upgrade to STRONG GO at any stage.
G3 FAIL (mechanism absent)       | Lepton masses are fitted inputs.
                                 | M6 is valid as alternative substrate;
                                 | weaker claim than Werbos markets.
```

---

## What M6 offers regardless of gate outcomes

Even in a G1 FAIL scenario, the Ouroboros model provides:

- Maxwell EM recovered exactly (structural fact — no gate needed)
- Charge quantization via Chern-Simons linking (Lean-proved from axioms)
- An independent test of the heat-as-oscillation hypothesis on a different substrate
- A simpler Taichi substrate than M5 (2-function ODE vs matrix M = ODO^T)
- A collaboration anchor with Paul Werbos (NSF X-Labs positioning)

These justify maintaining M6 as a sandbox/research track even if production is deferred.

---

## Gate status log

```text
Date       | Gate | Update
-----------|------|----------------------------------------------------------
2026-05-18 | G1   | PENDING — have correct Lean ODE; haven't run lepton scan
  morning  |      | with it yet. Werbos's bosonization session shows 1.1% muon
           |      | gap at ω=12.78 — promising. Running this is v4 Priority A.
2026-05-18 | G2   | PENDING — code request sent in email to Werbos. Our v3
  morning  |      | gets 0.508 MeV neutral m_χ; his estimate is 0.003-0.015 MeV.
           |      | 30× discrepancy. Likely energy functional difference.
2026-05-18 | G3   | OPEN — three sandbox rounds (v1, v2, v3) and 3 email exchanges
  morning  |      | have not produced a discrete ω selection mechanism. Werbos
           |      | himself flagged this as "the key open question."
-----------|------|----------------------------------------------------------
2026-05-18 | G1   | UPGRADED PRIORITY — Werbos's reply hypothesizes that
  19:39    |      | leptons are the 3 lowest stable modes in Q=1 sector.
           |      | Running v4 Track 1 lepton scan tests this DIRECTLY.
           |      | If our 3 lowest stable ω match leptons within 4-6%,
           |      | G3 is empirically resolved alongside G1.
2026-05-18 | G2   | CODE PROMISED — Werbos: "I'll share the Python source
  19:39    |      | for the Q=0 scan and the calibration run... Give me a
           |      | day or two to get the files to you." ETA: 2026-05-20.
           |      | Will resolve calibration localization + neutral ground
           |      | state + new Q14 (locked vs A=0 description conflict).
2026-05-18 | G3   | RECEIVING PATH FORWARD — Werbos acknowledges "empirical
  19:39    |      | for now"; analytic proof = "major mathematical
           |      | undertaking" deferred. But the v4 Track 1 lepton scan
           |      | provides an EMPIRICAL test: if 3 lowest stable modes =
           |      | leptons, the framework's discrete-spectrum claim has
           |      | numerical backing even without analytic proof. This
           |      | could effectively close G3 for production purposes.
-----------|------|----------------------------------------------------------
2026-05-19 | G1   | ATTEMPTED, BLOCKED. Built sandbox_v4/m6_v4_lepton_scan.py
  PM       |      | with Lean 2-fn ODE; quick ω-sweep found H decreases
           |      | with ω (opposite of Werbos's bosonization). ω=1.0
           |      | marginally fails localization. m_μ/m_e ratio target
           |      | 207, found 0.04. Diagnostic confirmed v3's
           |      | H_CODE_ELECTRON_CALIB=0.494 is wrong (real H ≈ 10.7
           |      | at calibration point). Root cause per Werbos's Q13:
           |      | 2-fn ODE is for NEUTRAL sector only; CHARGED scan
           |      | needs 4-fn (φ,α,ρ,β) toroidal ansatz. G1 cannot
           |      | progress until 4-fn ODE is implemented.
2026-05-19 | G2   | v3 RESULT INVALIDATED. Built sandbox_v4/m6_v4_t2_
  PM       |      | neutral_ground.py with wider B₀ ∈ [1e-5, 1.0] log
           |      | scan + golden-section refinement. NEGATIVE: H is
           |      | monotonic in B₀; no minimum exists. All "localized"
           |      | scan points violate Lean ≤4-node regularity (found
           |      | 13-29 nodes). Linear part is β''+β'/r-β/r²+λβ=0 =
           |      | J_1(√λr) Bessel — oscillatory, not bound. v3's 23
           |      | solutions + m_χ=0.508 MeV are artifacts of permissive
           |      | absolute-tail check. The cited DM result is not
           |      | scientifically defensible until canonical Q=0 ODE is
           |      | settled. Q14 now 3-way: locked (v2 failed) / A=0
           |      | (v3/T2 disproved) / Q_A≈0 (untested).
2026-05-19 | G3   | BLOCKED. Empirical-via-G1 path now blocked alongside
  PM       |      | G1 itself. Analytic-mechanism path still deferred
           |      | per Werbos's own admission. No movement until G1
           |      | unblocks (Werbos's code OR solo 4-fn).
-----------|------|----------------------------------------------------------
2026-05-19 | TRIG | Werbos sent DM paper draft (Dark_Matter_in_Universe
  10:51 AM |      | v0.txt) for ApJ review citing Griesi v3 m_χ=0.508 MeV
           |      | as numerical input to relic-abundance Ω_χh²≈0.1-0.2.
           |      | The cited result is the one T2 invalidated 30min
           |      | earlier. Rodrigo replied same morning with honest
           |      | feedback: artifacts finding, 3-way Q14 conflict,
           |      | 4-fn extraction work-in-progress, code request as
           |      | path to firm up m_χ before submission. Awaiting Paul's
           |      | response.
```

---

## Werbos reply update (2026-05-18 19:39)

The first round of gate-relevant clarification arrived. Key impacts:

- **G1**: Sharper test designed — v4 lepton scan can falsify or confirm
  Werbos's "lowest 3 stable modes = leptons" hypothesis directly.
- **G2**: Code coming in 1-2 days. Will close calibration + verify neutral
  ground state in one round.
- **G3**: Analytic proof deferred indefinitely (Werbos's own admission).
  Empirical G1 result becomes the effective gate.

**Refined GO/NO-GO**: with G1 + G2 both within reach this week, the
production decision can plausibly land by 2026-05-25.

---

## Status revision (2026-05-19, post-v4 first runs)

The 2026-05-18 "production decision by 2026-05-25" timeline is now too
optimistic. Both G1 and G2 hit unexpected blockers in v4 first runs:

```text
Gate | Was (2026-05-18)               | Is (2026-05-19)
-----|--------------------------------|--------------------------------
G1   | READY — run lepton scan w/     | BLOCKED — 2-fn ODE wrong tool
     | Lean 2-fn ODE; expect muon@    | for charged sector. Need 4-fn
     | 1.1% gap (Werbos bosonization) | benchmark ODE first.
G2   | CODE COMING — Werbos's source  | v3 RESULT INVALIDATED. T2 showed
     | will close calibration + m_χ   | 23 solutions are artifacts.
     | ground state in one round      | Existence in 2-fn form uncertain.
G3   | EMPIRICAL VIA G1               | BLOCKED with G1.
```

**Why the optimism was wrong**: we assumed the Lean 2-fn ODE was the
canonical form for both sectors. Per Werbos's Q13 reply, it is the
reduced form for the NEUTRAL sector only; the charged calibration +
lepton scan use the 4-function (φ, α, ρ, β) toroidal ansatz from his
pre-Lean code. The two ODEs are structurally different (toroidal l=0
vs vector l=1), not reductions of each other.

**Revised timeline**: production decision likely 2026-05-23 to
2026-05-27 if Werbos's code arrives this week; 2026-05-25 to
2026-05-30 if we have to implement the 4-fn ODE solo.

**Two unblock paths**:

```text
Path | Trigger                          | Resolves
-----|----------------------------------|---------------------------
A    | Werbos sends his Python source   | G1 + G2 + Q12 + Q14
     | (in response to our 2026-05-19   | in one shot
     | ask)                             |
B    | Solo: implement 4-fn benchmark   | G1 + G2 paths (with risk
     | ODE in sandbox_v4/, calibrate    | from 3 open gaps:
     | to H/Q=1.6969, then run scans    | ω-in-static, f(s) form, R)
```

The 2026-05-19 reply to Werbos was framed around Path A but acknowledges
Path B is the fallback if his code doesn't materialize.
