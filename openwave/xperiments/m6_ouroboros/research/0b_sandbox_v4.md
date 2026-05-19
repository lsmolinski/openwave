# 2026-05-18 evening — Sandbox v4: Lepton scan + awaiting Werbos's code

Triggered by Werbos's reply at 19:39 on 2026-05-18 (via DeepSeek; Paul is 78
and copy-pastes DeepSeek's drafts when tired). The reply addresses three of
the four immediate questions from v3, promises the Python source within 1-2
days, and gives an operational hypothesis for the discrete-ω selection
mechanism that v4 can test directly.

**Source:** Email thread continuation, "Intro & Nice to eMeet You" / reply
chain. Rodrigo responded at 20:10 acknowledging the late hour, will resume
2026-05-19.

---

## What Werbos's reply answered (vs v3 end-of-file items)

```text
v3 ID | v3 status            | Werbos reply (2026-05-18 19:39)
------|----------------------|--------------------------------------------------
Q12   | NOT ON ZENODO        | PROMISED — "I'll share the Python source for the
      | (asked directly)     | Q=0 scan and the calibration run... Give me a day
      |                      | or two to get the files to you." Under QAGI LLC.
Q13   | PENDING — 4→2 fn     | EXPLAINED — Lorenz constraint locks (V,Q) and
      | reduction?           | (A,J) into pairs:
      |                      |   Δ_r V = Q  → Q determined by V
      |                      |   Δ_r A = J  → J determined by A
      |                      | Operational rule: use 4-fn for CHARGED calibration,
      |                      | 2-fn for NEUTRAL scan. ⚠️ But conflicts with the
      |                      | bosonization log's "A-field = 0" approach. → Q14.
Q2    | OPEN (3 sandboxes)   | ACKNOWLEDGED OPEN — "At this stage, the ω values
      | discrete ω           | are empirical." HYPOTHESIS provided: leptons are
      | selection?           | the lowest 3 stable modes in Q=1 sector. Analytic
      |                      | proof = "major mathematical undertaking", deferred.
      |                      | v4 lepton scan tests this hypothesis directly.
Q8    | PARTIAL (v3)         | (not directly addressed; will be settled by code)
Q9    | RESOLVED (v3)        | (no change)
```

### New question surfaced (Q14)

Werbos's DeepSeek reply says for the NEUTRAL chaoiton: `J(r) = -A(r)` AND
`Q(r) = -V(r)` (locked ansatz → single β). But his earlier bosonization
session log (May 18 4PM) said the neutral chaoiton has `A-field = 0, pure
J-field` (no locked relation, just A turned off).

```text
Description         | Source                | Our experiment
--------------------|-----------------------|-------------------------------
Locked: J=-A, Q=-V  | DeepSeek reply 19:39  | Sandbox v2 — FAILED (no soln)
A-field = 0         | Bosonization 4PM log  | Sandbox v3 — 23 SOLUTIONS ✅
```

These are different setups giving different empirical outcomes. v3 (A=0)
worked; v2 (locked) didn't. The DeepSeek summary may be approximating; the
actual code in his bosonization session probably used A=0. The Python source
(Q12) will resolve this definitively.

→ **Q14: Which is the canonical Q=0 description — locked ansatz or A=0?**

---

## v4 implementation plan

```text
Track         | Action                                | Gate | Blocked by
--------------|---------------------------------------|------|----------------
TRACK 1       | Run lepton scan with 2-function       | G1   | Not blocked —
(can run now) | Lean ODE on full (α, β) at varying ω  |      | can run today
              | g=1.0625, λ=1.0, ω scan [0.5, 50]     |      |
              | Test: are ω ≈ 1.0, ~12.78, ~40.7      |      |
              |   the lowest 3 stable modes?          |      |
TRACK 2       | Find neutral chaoiton ground state    |      | Not blocked —
(can run now) | Tune B0 to minimize H_code            |      | refines v3
              | Currently 0.508 MeV at B0=0.010;      |      |
              | smaller B0 → ground state             |      |
TRACK 3       | Receive Werbos's Python source        | G2   | Werbos 1-2 days
(blocked)     | Compare our v3 energy functional to   |      |
              | his; close calibration localization;  |      |
              | resolve Q14 (locked vs A=0)           |      |
TRACK 4       | If Q14 = locked: re-run v3 with the   | -    | Q14 answered
(branched)    | true locked ansatz; if Q14 = A=0:     |      | by code (T3)
              | proceed with v3 as canonical          |      |
```

### Track 1 — Lepton scan (highest-value next test)

Setup:

```text
ODE:       α'' + (1/r)α' − α/r² + ω²α = β
           β'' + (1/r)β' − β/r² + ω²β = α − λβ − 4gβ³
Params:    g = 1.0625, λ = 1.0  (electron-calibration values)
ω scan:    [0.5, 50] in steps of 0.5 (~100 points)
BCs:       α ~ A₀·r, β ~ B₀·r at origin (slopes; field=0 at r=0)
           α, β → 0 at r_max = 30
Stability: at each ω, check conjugate-point (Gelfand-Fomin) per
           Lean theorem §stability_condition
Outputs:   - localized solutions found at each ω
           - H/Q at each localized ω
           - mass predictions via m = H × m_e / H_electron
           - flag the 3 lowest-energy ω values → candidate leptons
Expected:  per Werbos's bosonization (May 18 4PM):
             ω = 1.0   → electron 0.511 MeV (calibration)
             ω = 12.78 → muon 106.8 MeV (1.1% gap)
             ω = 15.0  → pion+ 145.1 MeV (3.9% gap)
             ω ≈ 40-41 → tau (predicted)
Test:      if our 3 lowest stable ω values match leptons →
           Werbos's hypothesis (lepton = stable mode) confirmed
```

Effort: 1-2 days. New script `sandbox_v4/m6_v4_lepton_scan.py`.

**Calibration dependency (carried over from v3 ⚠️):** T1 anchors on ω=1.0
producing the electron baseline. v3 ended with H/Q = 1.723 (1.56% gap)
NOT localized. If T1 still can't localize at ω=1.0 with current energy
functional, we have no baseline to convert H_code → MeV for muon/tau
predictions. Two fallback paths:

```text
Path | Action                              | Pre-condition
-----|-------------------------------------|-----------------------
A    | Wait for T3 (Werbos's code) and    | Best if code arrives
     | use his energy functional with     | soon (1-2 days)
     | toroidal geometry factor (2π)²R    |
B    | Try alternative energy formulas    | If T3 delayed; lower
     | (add (2π)²R prefactor explicitly;  | confidence but
     | rescale H by 1/r at r=8; etc.)     | unblocks T1
```

Recommended: try Path B for a few hours first (cheap experiment);
if localization is still elusive, wait for T3.

### Track 2 — Neutral chaoiton ground state

Refine v3 sandbox to find the minimum-energy localized β solution by
shooting on B0 with bisection. Current: 0.508 MeV at B0=0.010. Expected: a
slightly lighter ground state below B0=0.010 — possibly approaching
Werbos's 0.003-0.015 MeV estimate, possibly not.

Effort: ~1 day. Either new script or augment `sandbox_v3/m6_v3_chaoiton.py`.

### Track 3 — Werbos's code (waiting)

Will give us:

- The exact energy functional with toroidal geometry factor → closes
  the calibration 1.56% gap (currently not localized)
- The canonical Q=0 description (Q14) — locked or A=0
- The full lepton scan code (independent cross-check on Track 1)
- The conjugate-point stability test implementation

ETA: 1-2 days per Werbos's reply.

---

## SANDBOX v4 RESULTS, ROADBLOCKS, NEXT STEPS

(Placeholder — to be filled in after v4 runs)

### IMMEDIATE-QUESTIONS (gates v5)

```text
ID  | Question                                | First | Status at start of v4
    |                                         | surf  |
----|-----------------------------------------|-------|----------------------
Q12 | Werbos's Python source                  | v2    | PROMISED — 1-2 days
Q13 | 4-fn to 2-fn reduction                  | v3    | PARTIAL — DeepSeek
    |                                         |       | answer given, but
    |                                         |       | conflicts with
    |                                         |       | bosonization log
    |                                         |       | (A=0 vs locked).
    |                                         |       | Resolved by code.
Q14 | Canonical Q=0: locked J=-A,Q=-V  OR    | v4    | NEW — surfaced by
    | A-field = 0?                            |       | DeepSeek reply
    |                                         |       | conflict with
    |                                         |       | bosonization log.
    |                                         |       | Code will resolve.
```

### STILL OPEN-QUESTIONS

```text
ID  | Question                                | First | Status at start of v4
    |                                         | surf  |
----|-----------------------------------------|-------|----------------------
Q2  | Discrete ω selection mechanism          | 0a    | OPEN — Werbos says
    |                                         | §9.9  | empirical for now;
    |                                         |       | hypothesis = leptons
    |                                         |       | are lowest 3 stable
    |                                         |       | modes. v4 Track 1
    |                                         |       | will test directly.
Q3  | Analytical ω = 2mc²/ℏ derivation?       | 0a    | PARTIAL — only via
    |                                         | §9.9  | calibration (1.2%)
Q6  | QCD reconciliation (3-chaoiton proton)  | 0a    | OPEN — uninvestigated
    |                                         | §9.9  |
Q8  | Neutral m_χ true ground state           | v2    | PARTIAL — 0.508 MeV;
    |                                         |       | Track 2 will refine
```

### RESOLVED-QUESTIONS (carries forward from v3)

```text
ID  | Question                                | Resolution
----|-----------------------------------------|----------------------
Q1  | Distinction from Duda's LdGS            | Both are topological
Q9  | Linearized Bessel truncation?           | YES (v5 §6)
Q10 | Asymmetric ansatz J/A ≠ const?          | YES (Lean theorem)
Q11 | Sign convention J=−A or J=+A?           | OBSOLETE
```

### ARCHIVED-QUESTIONS

```text
ID  | Question                                | Why archived
----|-----------------------------------------|----------------------
Q4  | Single vs two-field ontology            | Aesthetic preference
Q7  | Cold-fusion citation trail              | Historical, not physics
```

### HARDEST-PIECES TRACKER

```text
Hardest piece     | Source       | v1 | v2 | v3 | v4 | Status start v4
------------------|--------------|----|----|----|----|--------------------------
V(M) potential    | Duda §III    | ✓  | ✓  | ✓  | ✓  | UNRESOLVED (shared M5)
form              |              |    |    |    |    |
f(J·J) form       | Werbos 2017+ | ✓  | ✓  | ~  | ~  | MOSTLY SETTLED — f=gs²;
                  | v5 confirms  |    |    |    |    | factor 2/4 needs check
ω quantization    | Werbos Spec  | ✓  | ✓  | ✓  | ✓  | OPEN — Werbos hypoth.:
mechanism         | §6.1         |    |    |    |    | leptons = 3 lowest
                  |              |    |    |    |    | stable modes. v4 tests.
Q=0 chaoiton      | v2 BVP sweep | —  | ✓  | ~  | ~  | PARTIAL — 23 solutions;
existence         |              |    |    |    |    | ground state pending;
                  |              |    |    |    |    | Q14 conflict to resolve
Lepton mass       | v1 vs v5     | ✓  | ✓  | ✓  | ?  | UNTESTED with Lean
spectrum          |              |    |    |    |    | ODE — Track 1 will
                  |              |    |    |    |    | settle this in v4
Electron          | v3 result    | —  | —  | ✓  | ✓  | OPEN — H/Q=1.723 (1.56%
calibration       |              |    |    |    |    | gap) but not localized.
localization      |              |    |    |    |    | T3 primary fix; T1
                  |              |    |    |    |    | baseline depends on it.
```

Active question count at start of v4: 3 IMMEDIATE (Q12, Q13, Q14) + 4 STILL OPEN (Q2, Q3, Q6, Q8) = **7 load-bearing**.

The single highest-value action: run Track 1 (lepton scan) tomorrow. If our 3 lowest stable ω match leptons within Werbos's 4-6% gap claim, Q2 is effectively resolved empirically — the framework's discrete spectrum claim has direct numerical backing.
