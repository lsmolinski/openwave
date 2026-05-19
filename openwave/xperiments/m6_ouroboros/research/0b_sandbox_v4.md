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

### T1 first run — 2026-05-19 (BLOCKED on wrong ODE form)

Built `sandbox_v4/m6_v4_lepton_scan.py` adapting the canonical Lean 2-fn
ODE from v3 into an ω-sweep. Quick scan run on ω ∈ {1, 5, 10, 12.78, 15,
20, 25, 30, 35, 40.7, 45} at g=1.0625, λ=1.0, (A0, B0) grid {0.05, 0.1,
0.5}². Total 99 integrations in 20.5s.

```text
Outcome              | Status
---------------------|--------------------------------------------------
Script works         | ✅ Localized solutions found at 10/11 ω values
                     |    (ω≥5; ω=1.0 marginally fails localization)
H scaling with ω     | ❌ DECREASING — opposite to Werbos's bosonization
                     |    H_full(12.78)/H_full(1.0) ≈ 0.04 (need 207 for muon)
Calibration H/Q      | ❌ Our H/Q at (ω=1.0, A0=B0=0.1) ≈ 1.04
                     |    Werbos: H/Q = 1.6969 (electron). 38% off.
v3 constant wrong    | ⚠️  H_CODE_ELECTRON_CALIB = 0.494 in v3 is wrong;
                     |    real H_spat at calib point ≈ 10.7 (close to
                     |    bosonization log's 11.4)
```

### Root cause: T1 is using the wrong ODE

Re-reading Q13 of Werbos's reply: the Lorenz constraint locks pairs
(V, Q) and (A, J), giving an operational rule:

```text
Sector     | ODE form used     | Status
-----------|-------------------|--------------------------------
Neutral    | 2-fn (Lean)       | ✅ v3 found 23 solutions
(α=0, ω=0) |                   |
Charged    | 4-fn (v5 paper)   | ❌ T1 used 2-fn — WRONG TOOL
(ω≠0)      | (φ, α, ρ, β)      |
```

The bosonization-session lepton scan that got ω=12.78 → muon 1.1% gap
ran Werbos's **earlier 4-function code**, *before* the May 18 Lean
reduction. The Lean 2-fn ODE is the structure-existence proof; it's
the reduced form for the Q=0 sector. Charged calibration + lepton scan
need the full 4-fn ansatz.

### What T1 needs to actually run

```text
Need                                  | Source / unblock
--------------------------------------|--------------------------------
Explicit 4-fn ODE (φ, α, ρ, β)        | (a) Werbos's Python source (Q12;
  with toroidal volume element and    |     ETA 2026-05-20, his promise)
  electron-calibration energy formula | (b) Manual extraction from
                                      |     v5 paper §6 + Numerical
                                      |     Benchmark docx (medium effort)
Correct H_CODE_ELECTRON_CALIB         | Implied by the above; v3's 0.494
                                      | is structurally wrong
Stability test (Gelfand-Fomin         | Needed even after ODE fix to rank
  conjugate-point)                    | "stable modes" per Werbos hypothesis
```

### Diagnostic kept

`sandbox_v4/diag_energy_functional.py` records the 2-fn H-scaling failure
and the H_full = H_spat + ω·Q identity. Keep for cross-check once 4-fn
code lands — the 2-fn form still applies to the Q=0 sandbox track.

### T1 decision

**HOLD T1.** Pivoting to T2 (neutral ground state via B₀ shooting) is
the right next move — it runs on the 2-fn Lean ODE we already have,
refines a v3 result rather than testing a hypothesis with the wrong
tool, and doesn't burn cycles guessing at Werbos's energy functional.

Resume T1 when:

```text
Trigger                              | Action
-------------------------------------|------------------------------
Werbos's code arrives (T3 unblocks)  | Adopt his 4-fn calibration
                                     | + lepton scan code as ground truth
OR ω=2026-05-20 +1d passes w/o code  | Manual 4-fn extraction from
                                     | v5 paper §6 + Numerical Benchmark
                                     | (slower, more uncertain)
```

---

### T2 first run — 2026-05-19 (NEGATIVE RESULT, corrects v3 overclaim)

Built `sandbox_v4/m6_v4_t2_neutral_ground.py`. Wide log-scan B₀ ∈
[1e-5, 1.0], 60 points, six λ values {0.1, 0.5, 1.0, 2.0, 5.0, 10.0},
r_max=30, plus golden-section bisection at λ=1.0.

```text
What we found              | What it means
---------------------------|----------------------------------------------
H is monotonic in B₀       | H grows like B₀^~3 (cubic-dominant scaling)
  across [1e-5, 1.0]       | for all λ. No minimum exists.
H → 0 as B₀ → 0            | "Ground state" of v3 was just trivial near-
                           | zero solution. Not a physical chaoiton.
Solution de-localizes at   | At λ=1.0: ~B₀=0.35. At larger λ: later.
moderate B₀ (~0.3-0.7)     | This is where the cubic term breaks the
                           | small-amplitude Bessel oscillation pattern.
Radial-node count 13-29    | Lean stability spec: ≤4 nodes. EVERY scan
                           | point violates this. The linear part is
                           | β''+β'/r-β/r²+λβ = J₁(√λ r) Bessel —
                           | inherently oscillatory.
```

**v3's "23 neutral chaoiton solutions" are artifacts.** Localization
was checked with `|β(r)| < 0.1 at r ≥ 8`, but at small B₀ the entire
oscillatory Bessel solution stays below 0.1 in absolute terms — not
because it decays, but because the amplitude itself is small. A
proper localization check needs EITHER node-count regularity (≤4 per
Lean) OR a decay-rate test (β² fall-off at large r relative to peak).

### What this implies for Q=0 / DM scan

The canonical Q=0 / neutral-chaoiton description is now in **triple
conflict**:

```text
Description              | Source                  | Empirical status
-------------------------|-------------------------|------------------
Locked: J=-A, Q=-V       | DeepSeek 19:39 reply    | v2 — FAILED
                         | (single β via locked    |   (no localized
                         |  ansatz)                |    solution found)
A-field=0, ω=0           | Bosonization 4PM log    | v3 claimed ✅;
  (α=0 reduction)        |                         | T2 DISPROVES
                         |                         |   (23 v3 solutions
                         |                         |    are artifacts)
Q_A≈0, Q_J≠0             | Bosonization log later  | UNTESTED — both
  (both fields active,   |   ("Q_A ≈ 0, Q_J ≠ 0")  | A and J fields
   EM charge zero but    |                         | active under
   nuclear charge ≠0)    |                         | constraint
```

→ This sharpens **Q14**: the canonical Q=0 description is one of THREE
candidates, not two. Werbos's code (Q12) is now the only clean path to
resolution. Without it, we cannot run a credible neutral chaoiton scan.

### T2 decision

**T2 CLOSED as negative.** Neutral-chaoiton existence in the canonical
Lean 2-fn form is not confirmed. v3's m_χ=0.508 MeV figure should be
treated as obsolete in roadmap and gates. Re-run only after Q14 is
settled.

### Sandbox v4 status (end of 2026-05-19 first-run)

```text
Track | Status     | Notes
------|------------|----------------------------------------------------
T1    | BLOCKED    | 2-fn ODE wrong tool for charged. Awaiting Werbos
      |            | code (Q12) OR manual 4-fn extraction.
T2    | CLOSED -   | v3 overclaim corrected. Q=0 chaoiton existence
      |            | now uncertain pending Q14 resolution.
T3    | WAITING    | Werbos's Python source, ETA 2026-05-20 (his
      | on Werbos  | promise). If slips → manual 4-fn extraction.
T4    | DEFERRED   | Q14 reconciliation rolls into T3 or 4-fn work.
T5    | DONE -     | Extraction complete (2026-05-19). 4-fn ODE
(new) | EXTRACTED  | found in Numerical Benchmark §3; calibration
      |            | anchor H/Q = 1.6969 confirmed. See section below.
```

### T5 — 4-fn ODE extracted (2026-05-19)

The Numerical Benchmark for the Ouroboros Lagrangian docx was **not
"wrong ODE form" as v3 notes claimed** — it's the canonical 4-fn ODE
source. The "wrong ODE" call was based on a misreading. Findings:

**4-fn benchmark ODE** (with mapping to our 2-fn notation):

```text
Mapping:  φ = V (A_0)   α = A (A_φ)   ρ = Q (J_0)   β = J (J_φ)
Convention: Δ_r f = f'' + (1/r)·f'    (NO -f/r² — toroidal, l=0)

Δ_r V  =  Q                                  (eq 1, Coulomb-type)
Δ_r A  =  J                                  (eq 2, Ampere-type)
Δ_r Q  =  V  +  m_J²·Q  +  λ·Q·(Q² − J²)     (eq 3)
Δ_r J  =  A  −  m_J²·J  −  λ·J·(Q² − J²)     (eq 4)
```

**Boundary conditions:**

```text
At r → 0:  V'(0) = A'(0) = Q'(0) = J'(0) = 0
           (VALUE BC, fields nonzero at origin — NOT slope BC)
           Calibration: V(0) = A(0) = Q(0) = J(0) = 0.1
At r → ∞:  all fields → 0, decay as modified Bessel K_0(κr)
```

**Energy functional** (toroidal volume element, NOT spherical 4πr²):

```text
H = (2π)²·R · ∫_0^∞ r dr [ (V')² + (A')² + (Q')² + (J')²
                          − V·Q + A·J
                          + (m_J²/2)·(Q² − J²)
                          + (λ/4)·(Q² − J²)² ]
```

**Charge functional** (Chern-Simons in toroidal form):

```text
Q_CS = (2π·R) ∫_0^∞ r dr [ A·(∂_r J) − J·(∂_r A) ]
Identity: L/Q = ω    →    2·L/Q = g-factor
```

**Calibration anchor (corrected):**

```text
Param                       | Value
----------------------------|------------------
g                           | 1.0625
λ                           | 1.0
ω                           | 1.0
V₀ = A₀ = Q₀ = J₀           | 0.1
H/Q (target, m_e/e_nat)     | 1.6875
H/Q (Werbos predicted)      | 1.6969 (0.56% gap)
Q_CS                        | 1 (integer, electron)
L/Q                         | ω = 1.000 → g_e = 2.000 (0.116%)
```

**v3 H_CODE_ELECTRON_CALIB = 0.494 is wrong**; the bosonization log's
H=11.4, Q=8.7 figures are from an earlier abandoned Compton calibration
iteration. The single authoritative anchor for charged sector is
**H/Q = 1.6969 at the canonical point.**

### T5 gaps (medium confidence — implementation risks)

```text
Gap                                  | Resolution path
-------------------------------------|------------------------------
ω doesn't appear in the static ODE   | Werbos's code (preferred) OR
  but calibration is at ω=1.0        | trial: replace m_J² with
                                     | (m_J² − ω²) in eq 3,4
f(s) two-form ambiguity              | Werbos's code OR trial:
  v5: f(s) = g·s²                    | benchmark form parametrizes
  benchmark: f = (m_J²/2)·s          | g via (m_J², λ) — need to
              + (λ/4)·s²             | derive the mapping
Toroidal radius R in (2π)²R          | Code-units choice; expect
                                     | R≈1 in normalized units
Numerical Benchmark vs               | Two ODE systems are STRUCTURALLY
  Lean theorem ODE                   | different. Lean: l=1 angular,
                                     | vector Laplacian, slope BC, ω².
                                     | Benchmark: l=0 toroidal, value
                                     | BC, static. Different ansätze.
```

The structural difference between Lean ODE and Benchmark ODE is the
*real* unknown — it's not that one is a "reduction" of the other.
They appear to be different ansatz choices for the chaoiton geometry.
Werbos's pre-Lean charged-scan code uses the Benchmark form; the Lean
theorem proves existence for a different (l=1 angular) ansatz.

### Decision needed (HUMAN)

Two paths forward:

```text
Path | Action                              | Risk
-----|-------------------------------------|----------------------------
A    | Implement 4-fn benchmark ODE now    | Medium — three open gaps
     | (T1 attempt 2 with correct ODE).    | (ω-in-static, f(s) form,
     | Calibrate against H/Q=1.6969 first. | R) need empirical tuning
     | If matches → run lepton scan with   | OR Werbos's code arriving
     | (m_J², λ) varied for ω-equivalents. | to settle them
B    | Wait for Werbos's code (ETA tomorrow| Low — settles all gaps
     | 2026-05-20). Use it as ground truth | cleanly. But delays
     | for both calibration + spectrum.    | progress by 1-2 days,
     |                                     | and code may not arrive
     |                                     | (Werbos is 78, may be
     |                                     | cautious about sharing
     |                                     | under QAGI LLC)
```

---

### Werbos clarifications — 2026-05-19 1:49 PM (via DeepSeek)

Werbos replied to the morning T2/T5 message with explicit answers to
all three implementation gaps and an endorsement of our Q14 finding.
Verbatim quote (DeepSeek-authored, Werbos-forwarded):

> Rodrigo — thank you for this. Finding that the 23 solutions don't
> survive a finer scan is exactly the kind of honest checking that
> separates real science from wishful thinking. I'm grateful.
>
> Three quick technical answers to your implementation gaps:
>
> **1. ω in the static benchmark ODE.** The benchmark document gives
> the static radial equations after separating time. The oscillation
> frequency ω is absorbed into the effective mass term. For the
> electron calibration at ω = 1.0, the mass parameter in the ODE is
> m_eff² = m_J² - ω². With the calibrated parameters, this effective
> mass is what drives the exponential decay. The static ODE in the
> benchmark is correct as written; you just need to use m_eff rather
> than m_J when you set up the solver.
>
> **2. Two parametrizations of f(s).** The v5 paper uses f(s) = g·s²
> as a simplified notation for the pure quartic case. The benchmark
> document uses the full form f(s) = ½ m_J² s + ¼ λ s². They are
> consistent when m_J² = 0 and λ = 4g. For the electron calibration,
> both forms should give the same result if you match the parameters
> correctly. Use the benchmark form — it's more general.
>
> **3. Toroidal major radius R.** The factor (2π)²R in the energy
> functional is a geometric factor from integrating around the torus.
> It cancels out when you compute the ratio H/Q because both the
> energy and the charge contain the same factor. For the mass-to-
> charge ratio, you can set R = 1 in scaled units and the ratio will
> be correct. The absolute mass scale is fixed by matching the
> electron mass, which absorbs any overall scale factor.
>
> **On the neutral chaoiton:** Your discovery that the three candidate
> definitions give different masses is important. This is not a bug —
> it means the neutral sector has a richer structure than we initially
> thought. The correct definition for the dark matter candidate is
> most likely the one where both fields are active but the net EM
> charge is zero (Q_A ≈ 0, Q_J ≠ 0). This is the configuration that
> would actually be produced in the early universe. The pure J-field
> case and the locked ansatz are limiting approximations.
>
> I'll get you the calibration code as soon as I can. In the meantime,
> try the m_eff fix for the electron and see if that closes your 1.56%
> gap.
>
> The structural argument in the ApJ paper — axion subsumption,
> no-free-parameter relic match, sensor architecture — holds
> regardless of the exact numerical value of m_χ. But I agree we
> should get m_χ on solid footing before submission. Once the
> canonical ODE is locked, your offer to run the full ground-state
> scan is exactly what's needed.

### Resolution summary

```text
Gap                  | Pre-Werbos state    | Post-Werbos resolution
---------------------|---------------------|------------------------
ω in static ODE      | Unknown — possibly  | EXPLICIT: replace m_J²
                     | absorbed into m_J²  | with m_eff² where
                     |                     | m_eff² = m_J² − ω²
                     |                     | Direct test: try at the
                     |                     | electron calibration
                     |                     | point, expect 1.56% gap
                     |                     | to close
f(s) two-form        | v5: f = g·s²        | EXPLICIT mapping:
                     | benchmark:          |   v5(g)  ≡  benchmark
                     |   f = (m_J²/2)·s    |   form with m_J²=0,
                     |       + (λ/4)·s²    |   λ = 4·g
                     |                     | Use benchmark form as
                     |                     | the canonical general
                     |                     | parametrization. For
                     |                     | g=1.0625: λ_bench=4.25
                     |                     | when m_J²=0.
Toroidal R           | Unspecified         | CANCELS in H/Q ratio.
                     |                     | Set R = 1 in scaled
                     |                     | units. Electron mass
                     |                     | anchor absorbs the
                     |                     | absolute scale factor.
Q14 (Q=0 canonical)  | 3-way conflict      | RESOLVED in principle:
                     | (locked, A=0,       | canonical DM candidate
                     |  Q_A≈0/Q_J≠0)       | is Q_A ≈ 0, Q_J ≠ 0
                     | from DeepSeek vs    | (both fields active).
                     | bosonization log    | Locked and A=0 are
                     |                     | limiting approximations.
                     |                     | "produced in the early
                     |                     | universe" — physical
                     |                     | early-cosmology framing
Code arrival         | "1-2 days"          | "as soon as I can" —
                     | (committed)         | no specific date;
                     |                     | implies slip from
                     |                     | original ETA
```

### Concrete next test (Werbos proposed)

> "Try the m_eff fix for the electron and see if that closes your
> 1.56% gap."

This is the single highest-leverage validation: implement the 4-fn
benchmark ODE with the m_eff² substitution and run at the canonical
electron calibration point. If H/Q lands at 1.6969 (within Werbos's
0.56% target gap), we have the correct charged ODE form.

### Implementation plan — v4 Path A → CONFIRMED GO

```text
Step | Action                                  | Outcome
-----|-----------------------------------------|----------------------
1    | Document this Werbos reply (DONE)       | ✅ verbatim above
2    | Build m6_v4_4fn.py — 4-fn benchmark     | new script
     | ODE with m_eff² = m_J² − ω², value      |
     | BCs f'(0)=0, R=1 in volume element     |
3    | KEY MILESTONE: calibrate against H/Q =  | success → unblocks
     | 1.6969 at (g, λ, ω, V₀=A₀=Q₀=J₀)        | everything; failure
     | = (1.0625, 1.0, 1.0, 0.1)               | → ask one targeted
     |                                         | Werbos clarification
4    | Charged lepton scan ω ∈ [0.5, 50]       | test "lowest 3 stable
     | with calibrated m_J² fixed              | = leptons" hypothesis
5    | Q_A ≈ 0 neutral scan (canonical DM      | true m_χ for ApJ
     | candidate per Werbos)                   | relic-abundance calc
```

The Q12 (code) ask is still LIVE but no longer blocking. If Werbos's
code arrives mid-implementation, we cross-check against our results.

### Update to question statuses

Several IMMEDIATE-QUESTIONS resolve with this reply (see updated
tables below). Q14 in particular moves from "3-way conflict" to
"resolved-in-principle: Q_A≈0/Q_J≠0 is canonical".

---

### IMMEDIATE-QUESTIONS (gates v5)

```text
ID  | Question                                | First | Status post-Werbos
    |                                         | surf  | 2026-05-19 1:49 PM
----|-----------------------------------------|-------|----------------------
Q12 | Werbos's Python source                  | v2    | LIVE but DEMOTED —
    |                                         |       | "as soon as I can"
    |                                         |       | (no committed date).
    |                                         |       | No longer blocking
    |                                         |       | progress; useful
    |                                         |       | for cross-check.
Q13 | 4-fn to 2-fn reduction                  | v3    | RESOLVED in
    |                                         |       | principle. v5
    |                                         |       | g·s² ≡ benchmark
    |                                         |       | with m_J²=0, λ=4g.
    |                                         |       | Use benchmark form.
Q14 | Canonical Q=0: locked / A=0 / Q_A≈0?    | v4    | RESOLVED. Canonical
    |                                         |       | DM candidate is
    |                                         |       | Q_A≈0 / Q_J≠0
    |                                         |       | (both fields active,
    |                                         |       | EM-neutral). Locked
    |                                         |       | and A=0 are limiting
    |                                         |       | approximations.
Q15 | Does m_eff² = m_J² − ω² substitution    | v4    | NEW — surfaced by
    | close electron calibration H/Q to       |       | Werbos's proposed
    | 1.6969?                                 |       | direct test. Step 3
    |                                         |       | of v4 4-fn build
    |                                         |       | tests this.
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
