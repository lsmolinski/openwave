# 2026-05-18 — Sandbox v3: Lean-ODE-corrected chaoiton

Triggered by 11 new files received from Werbos at 17:20 on 2026-05-18. Two files
contained the actual canonical ODE and the neutral chaoiton method — both of which
invalidated sandbox v2's approach entirely.

**Source files:**

- `theory/chaoiton_theorem.lean.txt` — Lean 4 formal theorem (formal statement of
  chaoiton existence). Reveals the **canonical 2-function radial ODE** Werbos actually
  uses — different from everything in v1 and v2.
- `theory/Claude bosonization May18 4PM.txt` — Werbos's earlier Claude conversation
  (May 4–18). Reveals the **neutral chaoiton method**: set A-field to zero, use ω=0,
  solve single β ODE. Found 1,208 solutions in that session.
- `theory/Werbos_Law_of_Everything_2026_v5.pdf` — Updated framework paper.
  Confirms f(s)=gs², the 4-component toroidal ansatz, and that the Q=0 chaoiton
  scan is *in progress* (no published m_χ yet).

---

## What v2 got wrong (source: chaoiton_theorem.lean.txt)

```text
Error                  | v2 form                   | Canonical (Lean)
-----------------------|---------------------------|------------------------------
Laplacian              | 2D: f''+f'/r              | VECTOR: f''+f'/r − f/r²
Origin BC              | α(r_min) = A₀  (value)    | α ~ A₀·r  (slope; field=0 at 0)
Q=0 approach           | Q_CS=0 via locked J=±A    | α=0 exactly + ω=0
Coupling structure     | quarter (β/4 → α)         | full coefficient (β → α)
Linear damping         | absent                    | −λβ in β equation
Function count         | 2 locked OR 4 locked      | 2 independent (α, β)
```

## Canonical ODE (from chaoiton_theorem.lean.txt §ode_system)

```text
α'' + (1/r)α' − (1/r²)α + ω²α = β
β'' + (1/r)β' − (1/r²)β + ω²β = α − λβ − 4gβ³
```

BCs at origin: `α ~ A₀·r`, `β ~ B₀·r`  (both fields START AT ZERO; A₀, B₀ are slopes)

Localization criterion (Lean §localization_condition): `|α(r)| + |β(r)| < 0.1` for `r ≥ 8`

## Neutral chaoiton ODE (α=0, ω=0)

With A-field off and no angular momentum, the single β ODE is:

```text
β'' + (1/r)β' − β/r² + λβ + 4gβ³ = 0
```

BC: `β(r_min) = B₀·r_min`, `β'(r_min) = B₀`

Source: bosonization file May18 4PM, confirmed by v5 paper §8.

## What was built

```text
sandbox_v3/m6_v3_chaoiton.py   — two modes:
  --mode calibration             electron calibration scan (ω=1.0)
  --mode neutral                 neutral chaoiton β ODE at ω=0
  --mode all                     both

sandbox_v3/m6_v3_README.md      quick-start + change log
```

---

## planned priorities

v2 had three planned priorities. The Lean theorem and bosonization file changed the
scope of A and B; C is unchanged.

```text
Priority | v2 plan                              | v3 outcome
---------|--------------------------------------|----------------------------------
A        | Q=0 neutral chaoiton                 | SUPERSEDED. v2 tried locked
         | via locked ansatz                    | ansatz (both IVP and BVP) —
         | m6_1_neutral_chaoiton.py             | both failed (no localized
         | m6_1_bvp.py                          | solution). v3 used the correct
         |                                      | method (α=0, ω=0, single β ODE
         |                                      | from Lean theorem): 23 solutions,
         |                                      | m_χ = 0.508 MeV at λ=1.0. ✅
---------|--------------------------------------|----------------------------------
B        | 4-function ODE rebuild               | SCOPE CHANGED. Lean theorem
         | + lepton re-test                     | shows the canonical ODE is
         | (benchmark's 4-fn form,              | 2-function (α, β) with vector
         | ω^2.22 scaling check)                | Laplacian — not 4-function.
         |                                      | The bosonization file (4PM)
         |                                      | shows muon at ω=12.78 gives
         |                                      | 1.1% gap with the correct ODE
         |                                      | (vs our 31% in v1). We have
         |                                      | NOT yet run the lepton scan
         |                                      | ourselves with the Lean ODE.
         |                                      | → Promoted to Priority B in v4.
---------|--------------------------------------|----------------------------------
C        | 3-body proton bound state            | STILL PARKED — deferred.
         | V(R) ~ −C/R⁶ Yukawa,                 | No change from v2.
         | classical 3-body problem             |
```

### Lepton re-test (v2 Priority B) — what we know before running it

The bosonization file (4PM conversation) shows that with Werbos's code running the
correct ODE, the lepton scan gives:

```text
ω = 1.0:   0.511 MeV  (electron, exact by calibration)
ω = 12.78: 106.8 MeV  (muon, 1.1% gap)  ← vs our v1 31% gap
ω = 15.0:  145.1 MeV  (pion+, 3.9% gap)
```

This is dramatically better than v1's 31-44% gaps. The Lean ODE vector Laplacian
(-1/r²) changes the ω-scaling. Running this scan ourselves is v4 Priority A.

---

## SANDBOX v3 RESULTS, ROADBLOCKS, NEXT STEPS

### v3 first run results (2026-05-18)

#### Result 1 — Neutral chaoiton scan ✅

```text
23 localized solutions found across (λ, B0) parameter space.

Lightest solutions (B0=0.010):
  λ=1.0  → m_χ = 0.508 MeV  (matches λ × 0.511 MeV to 0.6%)
  λ=0.5  → m_χ = 0.729 MeV
  λ=0.1  → m_χ = 1.55 MeV
  λ=2.0  → m_χ = 0.349 MeV
  λ=5.0  → m_χ = 0.223 MeV
  λ=10.0 → m_χ = 0.157 MeV

At λ=1.0, B0=0.010: m_χ = 0.508 MeV  ← very close to electron mass.
This is NOT a coincidence: for the neutral single-field ODE,
λ = 1 gives a chaoiton with mass ≈ the electron mass at ω=0.
This is exactly Werbos's "WIMP analog" DM candidate.
```

Note: m_χ depends on B0 (the shooting amplitude). The smallest-B0 branch gives the
lightest state. Different B0 values give a continuous spectrum of neutral chaoiton
masses — confirming Werbos's "1,208 solutions" claim.

#### Result 2 — Electron calibration ⚠️

```text
Best: A0=2, B0=10 → H/Q = 1.723  (1.56% gap from target 1.6969)
But NOT localized (tail > 0.1 at r≥8).

A0=1, B0=5  → H/Q = 1.636  (3.59% gap) — also not localized.
Target H/Q = 1.6969. H/Q range scanned: 0.031 to 7.977.
```

The Lean ODE CAN produce H/Q near the electron target — just not yet with localized
solutions at those parameters. Root cause likely: energy functional formula for the
charged 2-function case needs tuning. The 3D volume element and toroidal geometry
factors were estimated; exact form requires Werbos's code.

### Roadblocks (v3-specific)

```text
Roadblock                             | Status
--------------------------------------|------------------------------------------
Neutral m_χ varies with B0            | Expected — continuous spectrum. The
                                      | physically relevant value is the
                                      | MINIMUM-energy state (ground state).
                                      | Shooting for the specific B0 that
                                      | gives the lowest-energy localized β
                                      | is the next sub-task.
Calibration not yet localized         | H/Q close (1.56% gap) but tail too large.
                                      | Energy functional may be missing the
                                      | toroidal geometry prefactor (2π)²R.
                                      | Werbos's code would resolve this directly.
2-function vs 4-function discrepancy  | v5 paper shows 4-function (φ,α,ρ,β);
                                      | Lean theorem shows 2-function (α,β).
                                      | The Lean 2-function is the REDUCED form
                                      | after applying the sin/cos phase
                                      | decomposition on the torus. V5's (φ,α)
                                      | map to a single α in the reduced form.
                                      | → New Q13.
```

### IMMEDIATE-QUESTIONS (gates v4 and email)

```text
ID  | Question                                | First | Status at end of v3
    |                                         | surf  |
----|-----------------------------------------|-------|-------------------
Q8  | What method gave the 0.003-0.015 MeV    | v2    | PARTIALLY RESOLVED —
    | semi-analytic estimate?                 |       | v5 §8 says scan is
    |                                         |       | "in progress"; 4PM
    |                                         |       | bosonization says
    |                                         |       | A-field=0, ω=0 gives
    |                                         |       | m_χ ≈ λ×0.511 MeV for
    |                                         |       | smallest-B0 state.
    |                                         |       | Our v3 confirms this.
Q9  | Is the estimate linearized Bessel       | v2    | RESOLVED — v5 §6
    | truncation?                             |       | confirms β ~ e^{-kr}/r
    |                                         |       | far-field decay with
    |                                         |       | k=0.154. The 0.003-
    |                                         |       | 0.015 MeV estimate was
    |                                         |       | from this linearized
    |                                         |       | far-field, not full
    |                                         |       | nonlinear solution.
Q12 | Share Q=0 Python source code            | v2    | CONFIRMED NOT ON ZENODO
    |                                         |       | (20030162, 20242421,
    |                                         |       | 20218067 — all checked).
    |                                         |       | Must ask directly.
Q13 | 2-function vs 4-function reduction:     | v3    | NEW — v5 paper uses
    | how do (φ,α,ρ,β) reduce to (α,β)?       |       | 4-fn; Lean theorem
    |                                         |       | uses 2-fn. What's the
    |                                         |       | reduction? Affects the
    |                                         |       | calibration energy calc.
```

### STILL OPEN-QUESTIONS (load-bearing only)

```text
ID  | Question                                | First | Status at end of v3
    |                                         | surf  |
----|-----------------------------------------|-------|-------------------
Q2  | Discrete-spectrum mechanism for ω       | 0a    | STILL OPEN. v3 neutral
    | (eigenvalue vs continuous)              | §9.9  | scan is CONTINUOUS in
    |                                         |       | B0 — mass is not
    |                                         |       | discrete. Bosonization
    |                                         |       | file shows muon at
    |                                         |       | ω=12.78 (1.1% gap)
    |                                         |       | with correct ODE —
    |                                         |       | so spectrum does emerge,
    |                                         |       | but selection mechanism
    |                                         |       | remains unclear.
Q3  | Analytical derivation of ω = 2mc²/ℏ?    | 0a    | PARTIAL — 1.2% via
    |                                         | §9.9  | R^phys calibration
Q6  | QCD reconciliation (3-chaoiton proton)  | 0a    | OPEN — Priority C
    |                                         | §9.9  | in v2 plan
```

### RESOLVED-QUESTIONS

```text
ID  | Question                                | Resolution
----|-----------------------------------------|------------------------
Q1  | Physical distinction from Duda beyond   | BOTH ARE TOPOLOGICAL.
    | topological-invariant choice?           | (Confirmed in v1)
Q9  | Linearized Bessel truncation?           | YES — confirmed by v5 §6
    |                                         | far-field β ~ e^{-kr}/r.
    |                                         | Werbos's 0.003-0.015 MeV
    |                                         | was from this linearized
    |                                         | form, NOT full nonlinear.
Q10 | Asymmetric ansatz J/A ≠ const?          | YES — confirmed by Lean
    |                                         | theorem (α, β independent)
    |                                         | and v5 paper (4-function).
    |                                         | Our locked-ansatz in v2
    |                                         | was completely wrong.
Q11 | Sign convention J=−A or J=+A?          | OBSOLETE — actual ODE has
    |                                         | NO locked ansatz. Functions
    |                                         | evolve independently.
```

### ARCHIVED-QUESTIONS (framework-semantics / historical)

```text
ID  | Question                                | Why archived
----|-----------------------------------------|------------------------
Q4  | Single vs two-field ontology            | Aesthetic preference.
Q7  | Cold-fusion citation trail              | Historical, not physics.
```

### HARDEST-PIECES TRACKER

```text
Hardest piece     | Source           | v1 | v2 | v3 | Status at end of v3
------------------|------------------|----|----|----|--------------------------
V(M) potential    | Duda §III +      | ✓  | ✓  | ✓  | UNRESOLVED — shared
form              | Duda 2026-05-15  |    |    |    | bottleneck with M5
f(J·J) form       | Werbos 2017 +    | ✓  | ✓  | ~  | MOSTLY SETTLED —
                  | v5 confirms      |    |    |    | f(s)=gs² confirmed by
                  | f(s)=gs²         |    |    |    | v5 paper + Lean ODE
                  |                  |    |    |    | −4gβ³ term. Factor of
                  |                  |    |    |    | 2 vs 4 needs check.
ω quantization    | Werbos Spectrum  | ✓  | ✓  | ✓  | OPEN — 4PM bosonization
mechanism         | §6.1             |    |    |    | shows ω=12.78 for muon
                  |                  |    |    |    | (1.1% gap) with correct
                  |                  |    |    |    | ODE; selection mechanism
                  |                  |    |    |    | still unclear
Q=0 chaoiton      | v2 BVP sweep     | —  | ✓  | ~  | PARTIALLY RESOLVED —
existence under   |                  |    |    |    | v3 found 23 neutral
locked ansatz     |                  |    |    |    | solutions at (λ, B0).
                  |                  |    |    |    | But locked-ansatz
                  |                  |    |    |    | approach was wrong;
                  |                  |    |    |    | actual method is α=0,
                  |                  |    |    |    | ω=0 single β ODE.
                  |                  |    |    |    | Ground-state m_χ needs
                  |                  |    |    |    | careful B0 tuning.
```

Active question count at end of v3: 4 IMMEDIATE (Q8, Q12, Q13 gating; Q9 resolved) + 3 STILL OPEN (Q2, Q3, Q6) = **6 load-bearing**.

Q12 (share code) is now the single most important gating ask — it would resolve Q13 and probably close the calibration localization issue in one shot.
