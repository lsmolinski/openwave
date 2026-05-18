# 2026-05-18 update — Werbos collaboration + benchmark spec + sandbox v2 plan

Triggered by Werbos's reply on 2026-05-17 to the sandbox-v1 reproduction email. Multiple emails received over the same evening; key contents:

- **Collaboration commitment.** Werbos: *"your feedback was by far the most serious feedback so far on the Ouroboros model"* — accepts collaboration, offers Ouroboros2 as a parallel branch if v1 needs amendments.
- **The benchmark document** (`theory/Numerical Benchmark for the Ouroboros Lagrangian.docx`) — DeepSeek-mediated, but the math is from Werbos. **First time the radial ODE is written explicitly in any Ouroboros document.** Resolves the structural ambiguity that forced sandbox-v1 to derive the equations from the Lagrangian.
- **Dark matter ask** (separate cover note) — Werbos's #1 priority for the next sandbox round: compute the Q=0 (neutral) chaoiton ground-state mass m_χ. This is the missing piece for the dark-matter prediction `Ω_χ h² ≈ 0.12` (Planck).
- **Two new theory papers added** to `theory/`:
  - `DarkMatterv1.pdf` — Werbos + DeepSeek, May 2026. Frames the neutral chaoiton as a unification of axion + WIMP via the Brennan 2024 JHEP "dual axion mass" result.
  - `Werbos_Chaoitons_Ouroboros_2025_with_far_field.pdf` — same content as the .docx already in `theory/`, updated PDF version.
- **Sawada outreach** — Werbos sent a direct note to Tetsuo Sawada on 2026-05-13 claiming vindication via Ouroboros. The contact address bounced; Werbos forwarded the exchange to me on 2026-05-18 as one of "many dialogues."
- **Bosonization conversation log** (`theory/TotalClaudeBosnoizationMay18.txt`) — full Werbos-Claude exchange from 2026-05-04 → 05-07 covering: the two-shoelace / Chern-Simons-linking framing of charge quantization, the mountain pass theorem as proof tool for chaoiton existence, far-field analysis (J-field range 6.5 units, A-field range 4.35 units), and the corrected `A` (not `A²`) coherent enhancement for Sawada's neutron-Pb anomaly.

## The benchmark document — explicit 4-function radial ODE

The benchmark spells out the static radial equations Werbos uses (their absence from the Calibration and Spectrum papers was the source of our sandbox-v1 ODE-form uncertainty):

```text
Component          | Benchmark form
-------------------|--------------------------------------------------
Laplacian          | 2D radial: Δ_r f = f''(r) + (1/r) f'(r)
Function count     | 4 (V, A, Q, J)
V equation         | Δ_r V = Q
A equation         | Δ_r A = J
Q equation         | Δ_r Q = V + m_J² Q + λ Q(Q² − J²)
J equation         | Δ_r J = A − m_J² J − λ J(Q² − J²)
f(s)               | (1/2) m_J² s + (λ/4) s²    [mass + quartic]
Time-dep ansatz    | V·sin(ωt), A·cos(ωt), Q·cos(ωt), J·sin(ωt)
BCs at r=0         | V'(0) = A'(0) = Q'(0) = J'(0) = 0  (regularity)
BCs at r → ∞       | fields → 0, exponential decay K_0(κr) Bessel
```

Calibration values (locked from electron):

| Parameter | Value |
| --- | --- |
| g | 1.0625 |
| λ | 1.0 |
| A₀ = B₀ | 0.1 (field amplitude at r ≈ 0) |
| ω (electron) | 1.0 |
| R^phys | 191 fm |
| Expected H/Q | 1.6969 |
| Expected L/Q | 1.0  (= ω, exact identity) |
| Expected Q_CS | 1  (integer) |

### What's different from our sandbox v1

```text
Aspect              | Sandbox v1               | Benchmark v2
--------------------|--------------------------|-------------------------
Function count      | 2 (α, γ) winner          | 4 (V, A, Q, J)
                    | 4 ansatz also tried      |
Laplacian           | 3D: (2/r) prefactor      | 2D: (1/r) prefactor
Cubic term          | 4g·γ³ (single-field)     | λ·Q(Q²−J²), −λ·J(Q²−J²)
Cross-coupling      | α↔γ direct symmetric     | V↔Q (Coul), A↔J (Amp)
ω in radial eqs     | explicit ω² × field      | absent as written (see 4)
m_J term            | not present              | (1/2)m_J²·s in f(s)
Integrator          | RK45                     | RK4 / DOP853 / bvp_solve
Method              | IVP (shoot-from-0)       | shooting OR collocation
Domain r_max        | 10 (calib) / 15 (spec)   | 20–50 scaled units
Grid points         | 500–800 (linear)         | 2000–5000 (log spacing)
Tolerances          | rtol=1e-10, atol=1e-12   | convergence check at 0.1%
```

Implication: our sandbox-v1 winner A4+H1+Q2 (2-function, 3D Laplacian, 4g·γ³) **is structurally incorrect** relative to Werbos's actual radial ODE — even though it reproduces H/Q at electron calibration. The 0.30% calibration match was coincidental in form; the v2 spec is materially different.

## The Q=0 neutral chaoiton — Werbos's #1 ask

From his cover note: *"the mass of the lightest neutral chaoiton is the last missing piece for the dark-matter prediction. If it comes out MeV-scale or above, the standard freeze-out calculation immediately gives Ω_DM. If it's very small, we have a dark radiation candidate instead. Either way, it's a result worth publishing."*

```text
Setup               | Value / form
--------------------|--------------------------------------------------
Locked ansatz       | J(r) = −A(r),  Q(r) = −V(r)
Parameters          | g=1.0625, λ=1.0, ω=1.0  (same as electron)
Q_CS                | 0  (no topological charge)
Output              | m_χ = H_code  (since Q_CS = 0, no H/Q normalization)
Werbos prelim est.  | 0.01–0.05 scaled units  ≈ 0.003–0.015 MeV
                    |   (linearized locked ansatz, semi-analytic)
Quartic correction  | May push m_χ to MeV scale
Decision branch     | m_χ ≥ MeV → cold DM, Ω_χ h² freeze-out
                    | m_χ ≪ MeV → dark radiation candidate
Planck target       | Ω_c h² = 0.120 ± 0.001
```

The Brennan 2024 JHEP claim (axion = longitudinal mode of massive vector field) Werbos cites for axion-WIMP unification: the J-field's massive longitudinal mode is the axion in this picture, so axion-search null results constrain m_J directly. Worth following the cite but not load-bearing for the Q=0 computation itself.

## Open structural questions in the benchmark spec

```text
Question                                | Tension with benchmark
----------------------------------------|---------------------------------
Where is ω in the radial equations?     | Eqs as written contain no ω.
                                        | Spectrum paper says "only ω
                                        | changes" for lepton scan — so
                                        | ω must enter somewhere. Likely
                                        | through m_J(ω) or as an added
                                        | −ω² term we have to insert.
What formula gives m_J from g?          | Doc says "determined by g" but
                                        | no formula. Probably m_J = √g
                                        | or m_J² = g — to be verified
                                        | empirically against H/Q=1.6969.
Why 2D Laplacian (1/r)?                 | Justified by thin-torus reduction.
                                        | But the major radius R doesn't
                                        | appear in the ODEs — so where
                                        | does R enter the scaling?
V↔Q, A↔J source structure?              | Doesn't fall out of straight EL
                                        | on the L = −F·F − G·G + J·A − f
                                        | Lagrangian. Likely from the
                                        | sin/cos pairing in the time ansatz
                                        | acting through ∂²_t terms.
```

These four ambiguities are what need pinning down before sandbox v2 can faithfully reproduce Werbos's setup. Email Werbos with these specifics once v2 is built and the ω-insertion choice is forced.

## Planned sandbox v2 — three priorities

```text
Priority | Task                                     | Effort   | Gate
---------|------------------------------------------|----------|------
A        | Q=0 neutral chaoiton, locked ansatz      | 1–2 days | m_χ value
B        | 4-function ODE rebuild + lepton re-test  | 3–5 days | ω^2.22 ?
C        | 3-body proton bound state                | 1–2 wks  | binding/r
```

**Priority A — Q=0 neutral chaoiton** (deliverable to Werbos this week):

1. New script `m6_1_neutral_chaoiton.py` in `sandbox_v2/`
2. Implement the 4-function ODE with `J = -A`, `Q = -V` substitution → reduces to 2 free functions (A, V) with the constraint baked in
3. RK4 with shooting at r=0 (regularity BCs) outward to r_max=30 (start), confirm K_0 decay
4. Sweep initial amplitude A₀ to find the lowest-energy localized solution
5. Compute `m_χ = H` using the energy functional from §3 of benchmark
6. Report m_χ in scaled units AND in MeV (via R^phys = 191 fm calibration)
7. If MeV-scale: compute σv → Ω_χ h² freeze-out estimate
8. If sub-keV: flag as dark-radiation candidate

**Priority B — 4-function ODE rebuild + lepton spectrum re-test:**

1. New script `m6_1_full_ansatz.py` (NOT the old `m6_0_full_ansatz.py` which had Vector Laplacian and worsened scaling)
2. Implement the **benchmark's exact** 4-function ODE with 2D Laplacian
3. Two ω-handling branches to test:
   - 3a: **m_J²(ω)** — try m_J² = g·ω² (ω embedded in mass term)
   - 3b: **explicit −ω²** — add `−ω²·field` term on LHS of each radial eq
4. Run at ω=1.0; tune the ω-handling until H/Q reproduces 1.6969 to 0.5%
5. Run at ω={11.0, 13.0, 40.7}; check whether muon/pion/tau masses match Werbos's 4–6% gaps OR our previous 31–44% gaps
6. Outcome decides whether ω^2.22 scaling reproduces with the corrected ODE

**Priority C — 3-body proton bound state** (parked, post-A/B):

1. Compute chaoiton-chaoiton pair potential V(R) using the far-field formula `~exp(−2k·R)/R` (Yukawa-screened) from the bosonization-log far-field analysis
2. Set up classical 3-body bound-state problem with this V(R)
3. Find minimum-energy configuration, extract binding + r_rms
4. Compare to proton mass 938 MeV / charge radius 0.84 fm
5. This is harder and not Werbos's urgent ask — defer until A and B land

## Sawada / nuclear-detection adjacencies (notes, not action)

Werbos forwarded the bounced 2026-05-13 email to T. Sawada (Japan) claiming Ouroboros vindicates the 1989/2003 long-range nuclear-force anomalies. The far-field analysis in the bosonization log gives the J-field a Yukawa range ≈ 6.5 scaled units, with A-coherent (linear in mass number, **corrected from earlier A² claim**) enhancement matching Sawada's neutron-Pb observation. This is downstream of the Q≠0 chaoiton work and not part of sandbox v2's load-bearing path. Tracking in [[reference-sawada-anomaly]].

## Comparison summary — v1 results vs v2 plan

```text
Sandbox v1 result                         | v2 expectation
------------------------------------------|------------------------------
H/Q = 1.6918 (0.30% off)                  | should reproduce to <0.5%
                                          |   with correct 4-fn ODE
ω^2.04 scaling, 31–44% lepton gaps        | unknown — depends on how ω
                                          |   enters benchmark eqs
H/Q amplitude-invariant                   | still expected (scale-free)
Q=0 not yet attempted                     | PRIORITY A — Werbos's ask
3-body proton not attempted               | PRIORITY C — deferred
```

## Collaboration-norms note

Werbos's reply pattern across the multi-email burst:

- He's using DeepSeek + Claude in his loop; some of his "technical replies" are AI-drafted (he's transparent about this — *"I pass on comments from DeepSeek, not Claude, because MY Claude timed out again"*).
- He sees the OpenWave sandbox as the "tool that doesn't time out" — values the persistent state we get from Claude Code (Opus 4.7) over his chat-window-bounded sessions.
- He's prioritizing breadth (multiple parallel asks: Q=0, finer ω, proton, Sawada outreach) over depth on any one. We should keep depth-on-Q=0 as the first deliverable to anchor the collaboration.
- Asked us to send him "anything else" — Werbos values raw context. When sandbox v2 lands, the email back should include the script, the plot, the numbers, and the open structural questions from §4.

This collaboration norm matters for sequencing — deliver ONE clean result (the Q=0 mass) before opening the four-way correspondence on muon/tau/proton/Sawada. Otherwise the thread sprawls and load-bearing physics gets buried under email volume.

## SANDBOX v2 RESULTS, ROADBLOCKS, NEXT STEPS

### First-run setup (2026-05-18)

Built `sandbox_v2/m6_1_neutral_chaoiton.py` (main script) + `m6_1_investigate.py` (sweep harness). Implements the full 4-function radial ODE from the benchmark document §3, with the locked ansatz `J=−A, Q=−V` as initial conditions at `r=r_min`, integrated outward via `scipy.solve_ivp` (RK45 default; LSODA fallback for stiff cases).

Goal: extract `m_χ = H_code` at electron parameters `(g=1.0625, λ=1.0, ω=1.0, A₀=B₀=0.1)`, convert to MeV via `R^phys = 191 fm`, decide between cold-DM and dark-radiation branch.

### Investigation sweep — 4 batches, 22 parameter points

Each row in the sweep covers a hypothesis the structural questions in § "Open structural questions in the benchmark spec" surfaced:

```text
Batch | What it tests                              | Outcome
------|--------------------------------------------|---------------------------
S1×6  | Werbos's locked_sign=−1 vs flipped +1 at   | sign=−1: Bessel J_0
      | m_J²=0, across r_max ∈ {30, 60, 120}       |   oscillation; energy H
      |                                            |   grows LINEARLY with
      |                                            |   r_max (14→30→60). Sign
      |                                            |   changes 2→4→8 over the
      |                                            |   tail. NO localization.
      |                                            | sign=+1: I_0 explosion;
      |                                            |   fields reach 1e25 to
      |                                            |   1e103. Slope = +0.99.
S2×8  | Asymmetric seed satisfying V²−A²=−m_J²/λ   | All 8 FAIL stiffness
      | at r_min, both signs, m_J²∈{.01,.1,.5,1}   |   ("step size below
      |                                            |   spacing"). The cubic
      |                                            |   term blows up.
S3×4  | ω-insertion ∈ {explicit, in_m_J}           | sign=−1 explicit: still
      |                                            |   Bessel-like.
      |                                            | in_m_J: stiff fail.
      |                                            | sign=+1 explicit:
      |                                            |   trivial V = const.
S4×2  | 3D Laplacian (sandbox-v1 convention)       | sign=−1: Bessel-like
      |                                            |   (H=2.9, smaller due
      |                                            |   to extra 1/r damping).
      |                                            | sign=+1: I_0 explosion.
```

Output: `plots/m6_1_investigation_sweep.json` (22 rows, raw data) + `plots/m6_1_neutral_chaoiton_profile.png` (default-run field profiles).

### Diagnostic — why IVP can't find the K_0 decay mode

The benchmark equations + locked ansatz force the consistency constraint `V² − A² = −m_J²/λ` at every r where V ≠ 0 (algebraic consequence of equating the top eqs `Δ_r V = Q` with the bottom eqs `Δ_r Q = V + m_J²Q + λQ(Q²−J²)` under `Q = −V`). The seed `V(0) = A(0) = 0.1` violates this when `m_J² ≠ 0` → S2 stiffness.

When `m_J² = 0` the consistency is trivially satisfied. The equations reduce to:

```text
sign = −1:  V'' + V'/r + V = 0   →   Bessel J_0 (oscillatory)
sign = +1:  V'' + V'/r − V = 0   →   K_0 + I_0 mixture
```

The **regularity boundary condition at r=0** (`V'(0)=0`, V(0) finite) selects:

```text
For sign = −1:  J_0  (since Y_0 → −∞ at r=0)
For sign = +1:  I_0  (since K_0 → −∞ at r=0)
```

Neither is the K_0 modified-Bessel decay the benchmark §4 itself predicts at infinity. The K_0 mode requires V(r→0) → −∞, which the regular-r=0 BC rules out.

### Empirical conclusion

**IVP with the regularity-at-zero BC structurally cannot find the Q=0 chaoiton in this framework.** This is not a bug in our reproduction — it's a clean mismatch between the solver class and the problem class. The benchmark §10 itself recommends `bvp_solve` (shooting or collocation). We chose IVP first because it's simpler; sandbox v2's negative result is what forces the move to BVP for v3.

### Next steps for v3

```text
Step | Action                                      | Effort
-----|---------------------------------------------|--------
A ✅ | Capture v2 results into this doc            | done
B    | Build BVP solver: m6_1_bvp.py in sandbox_v2 | 0.5-1 day
     | using scipy.solve_bvp. Install scipy as    |
     | dependency (add to pyproject.toml).        |
C+F  | Email Werbos with: v2 sweep result +       | 1 hour
     | IMMEDIATE questions Q8-Q12 (incl. ask for  |   to draft
     | his Python source per Spectrum paper §Refs |
     | "Code available on request from QAGI LLC") |
```

The order matters: B first because the BVP result may CHANGE the questions for Werbos (e.g. if BVP succeeds and yields m_χ, the ask shifts from "what method?" to "we got X MeV, does that match yours?"). Doing B → revise email → send C+F.

### B done — BVP solver result (2026-05-18)

Added `scipy>=1.13,<2.0` to `pyproject.toml`. Built `sandbox_v2/m6_1_bvp.py` using `scipy.integrate.solve_bvp` with:

- 8-state vector `[V, V', A, A', Q, Q', J, J']`
- 4 BCs at r=r_min (regularity: V'=A'=Q'=J'=0)
- 4 BCs at r=r_max (decay: V=A=Q=J=0)
- Initial guess: locked-ansatz exponential envelope `V₀·exp(-r/L)` etc.

**Sweep over 180 parameter points** (2 signs × 5 m_J² × 3 decay-lengths × 6 (V₀,A₀) pairs):

```text
Outcome                           | Count | Interpretation
----------------------------------|-------|---------------------------
Trivial solution V=A=Q=J=0        |  148  | The only solution most of
                                  |       | the time. Linear system,
                                  |       | homogeneous BCs.
Numerical-noise non-trivial       |   30  | max|V| < 0.001, H ≈ 0;
  (max|V| < 0.001)                |       | not a real solution.
Significant non-trivial           |    2  | max|V| ≈ 3.7. Both at
  (max|V| > 0.001)                |       | sign=±1, m_J²=0.5,
                                  |       | V₀=1.0, A₀=0.1, L=10
Failed convergence                |    0  | All converged.
```

The 2 "significant" non-trivial cases have:

```text
Field profile                | Diagnostic
-----------------------------|------------------------------------
V oscillates around 0        | NOT K_0 decay
Q oscillates around 0        | Bessel-like, not decaying
A ≡ 0 to machine precision   | NOT the locked-ansatz solution
J ≡ 0 to machine precision   |
Locked-ansatz drift |Q-sV|=5 | Locked ansatz is BROKEN
Q_CS ≈ 0 (10⁻¹⁰⁰+)           | Q_CS=0 only because A=J=0
```

**r_max sensitivity test** (THE decisive check):

```text
r_max | m_χ found (MeV) | Meaning
------|-----------------|----------------------------------
  20  |     0           | trivial
  25  |  2788           | non-trivial
  30  | 10524           | non-trivial (5× the r_max=25 value)
  40  |     0           | trivial
  50  |     0           | trivial
  60  |  2053           | non-trivial (different from above)
```

A physical eigenstate would have an **r_max-invariant mass**. The fact that m_χ varies by 5× with r_max — and only exists at specific r_max values where Bessel zeros happen to coincide with our truncation — proves these are **r_max artifacts**, not real chaoitons.

### Definitive empirical conclusion

**The benchmark ODE + locked ansatz + Q_CS=0 + decay-at-infinity BC has NO physical chaoiton solution.**

Both IVP (sandbox v2 first attempt) and BVP (this attempt) yield only:

- Trivial V=A=Q=J=0  (the "no particle" state), or
- Bessel-like oscillations that don't decay (artifacts of finite r_max)

The "non-trivial" BVP solutions break the locked ansatz (|Q−sV| drift ≈ 5) AND have A=J=0, so they're not the locked-ansatz Q=0 chaoiton Werbos describes — they're a separate V-only sector.

### What this means for v3

Two interpretations of the negative result:

```text
Interpretation                          | Consequence
----------------------------------------|--------------------------
A. Our solver / setup is wrong          | Werbos's code will reveal
                                        |   the missing ingredient
B. Werbos's framework genuinely lacks   | Cold-DM prediction in
   a localized Q=0 chaoiton             |   DarkMatterv1 is overstated;
                                        |   collaboration value shifts
                                        |   to lepton spectrum work
```

Cannot distinguish A vs B without Werbos's actual code (request via Q12) OR a completely different ansatz (asymmetric J/A; ω-quantized; eigenvalue). **C+F email is the right next move.**

### IMMEDIATE-QUESTIONS (gates next steps)

```text
ID  | Question                                       | First surf
----|------------------------------------------------|-----------
Q8  | We tried IVP first (Bessel J_0 / I_0, no       | v2 (BVP)
    | decay), then BVP (only trivial + r_max-        |
    | artifact solutions). What numerical method     |
    | actually gave you the 0.003-0.015 MeV semi-    |
    | analytic estimate in DarkMatterv1 §5?          |
Q9  | Is the "semi-analytic" estimate based on the   | v2 (BVP)
    | LINEARIZED Bessel solution truncated at some   |
    | characteristic length, rather than a true      |
    | nonlinear localized chaoiton?                  |
Q10 | The locked ansatz J=±A, Q=±V seems             | v2 (BVP)
    | incompatible with Q²−J² ≠ 0 (needed for the    |
    | cubic to give localization). Is there an       |
    | asymmetric ansatz J(r)/A(r) ≠ const that you   |
    | actually use?                                  |
Q11 | Sign convention — DarkMatterv1 §5 says J=−A;   | v2 (BVP)
    | the K_0 decay structure requires J=+A. Which   |
    | is it?                                         |
Q12 | Could you share your Python source for the     | v2 (BVP)
    | Q=0 calculation? (Spectrum paper §Refs: "Code  |
    | available on request from QAGI LLC")           |
```

### STILL OPEN-QUESTIONS (load-bearing only)

```text
ID  | Question                                | First | Status at end of v2
    |                                         | surf  |
----|-----------------------------------------|-------|-------------------
Q2  | Discrete-spectrum mechanism — is there  | 0a    | ANSWERED EMPIRICALLY
    | a stability/quantization criterion      | §9.9  | (NEGATIVE) — v1 scan
    | beyond |A|+|J|<0.15 that selects        |       | sees all ω in [1,60]
    | ω={1,11,13,40.7} as DISCRETE            |       | localized equally;
    | eigenvalues? AND do the predicted       |       | reproduction shows
    | lepton masses match Werbos's claimed    |       | 31-44% lepton gaps
    | 4-6% gaps, or larger?  [merges old Q5]  |       | (not 4-6%); ω^2.22
    |                                         |       | = log(207)/log(11)
    |                                         |       | suggests fitted.
    |                                         |       | Sent to Werbos
    |                                         |       | 2026-05-17 — no
    |                                         |       | reply yet.
Q3  | Analytical derivation of ω = 2mc²/ℏ?    | 0a    | PARTIAL — 1.2% via
    |                                         | §9.9  | R^phys calibration,
    |                                         |       | not analytical
Q6  | QCD reconciliation (3-chaoiton proton)  | 0a    | OPEN — Priority C
    |                                         | §9.9  | in v2 plan,
    |                                         |       | uninvestigated
```

### RESOLVED-QUESTIONS

```text
ID  | Question                                | Resolution
----|-----------------------------------------|------------------------
Q1  | Physical distinction from Duda's LdGS   | BOTH ARE TOPOLOGICAL.
    | beyond topological-invariant choice?    | 0a §2.1, §2.4, §9.4 and
    |                                         | Werbos's v2 paper §3
    |                                         | confirm Ouroboros uses
    |                                         | Chern-Simons linking =
    |                                         | Hopf invariant. Duda
    |                                         | uses Brouwer degree. The
    |                                         | distinction is invariant-
    |                                         | type + field structure
    |                                         | (matrix M vs two-vector
    |                                         | A,J), NOT topology-vs-
    |                                         | no-topology.
```

### ARCHIVED-QUESTIONS (framework-semantics / historical)

```text
ID  | Question                                | Why archived
----|-----------------------------------------|------------------------
Q4  | Single-deeper-field vs two-field        | Werbos 2017 explains
    | ontology                                | two-field via toroidal-
    |                                         | poloidal mutual
    |                                         | confinement. Duda's
    |                                         | "single deeper field"
    |                                         | preference is aesthetic,
    |                                         | not falsifiable physics.
Q7  | Cold-fusion citation trail              | Historical/credibility
    |                                         | question; not framework
    |                                         | physics. Tracking
    |                                         | separately if needed.
```

Active question count at end of v2: 5 IMMEDIATE + 3 STILL OPEN = **8 load-bearing**.
Q5 merged into Q2 (same finding, theoretical + empirical angles).

### HARDEST-PIECES TRACKER

Stable across all checkpoints — nothing resolved v1 → v2:

```text
Hardest piece     | Source           | v1 | v2 | Status at end of v2
------------------|------------------|----|----|---------------------
V(M) potential    | Duda §III +      | ✓  | ✓  | UNRESOLVED — Duda
form              | Duda 2026-05-15  |    |    | admits "the most
                  | email            |    |    | difficult"; gates
                  |                  |    |    | both M5 and M6
f(J·J) form       | Werbos 2017      | ✓  | ✓  | UNRESOLVED — "any
                  | paper "any       |    |    | nonneg fn"; specific
                  | nonneg fn"       |    |    | form is calibration
                  |                  |    |    | choice
ω quantization    | Werbos Spectrum  | ✓  | ✓  | UNRESOLVED + v1
mechanism         | §6.1 "the key    |    |    | sharpened: scan
                  | open question"   |    |    | shows ω continuous,
                  |                  |    |    | not discrete
Q=0 chaoiton      | NEW — v2 BVP     | —  | ✓  | EMPIRICALLY NEGATIVE
existence under   | sweep            |    |    | under both IVP and
locked ansatz     |                  |    |    | BVP; pending Werbos
                  |                  |    |    | clarification (Q8-Q12).
                  |                  |    |    | If unresolved, the
                  |                  |    |    | DarkMatterv1 cold-DM
                  |                  |    |    | prediction is at risk.
```

Pattern: the framework's hardest pieces are **stable and unresolved**. Each new sandbox round surfaces a NEW hardest piece (Q=0 existence for v2) without resolving the prior ones. This is a load-bearing observation for the email to Werbos — we're not running out of questions; we're accumulating them.
