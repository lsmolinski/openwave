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
