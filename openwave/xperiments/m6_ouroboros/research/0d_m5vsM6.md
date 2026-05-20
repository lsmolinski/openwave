# M5 vs M6 — Current Evaluation in Service of the SABER Trigger

**Purpose:** snapshot the current state of OpenWave's two production-candidate
models (M5 Liquid Crystal / Duda, M6 Ouroboros / Werbos) against the SABER
engineering goal. SABER's MAIN GOAL is Direct Heat Conversion via subatomic
modulation; the thermal-amplitude hypothesis frames heat as joint (A, ω)
excess of subatomic defect oscillation. SABER engineering can only trigger
once a working subatomic simulator validates this hypothesis empirically
(falsifiability test = 5b.1 Sine-Gordon kink, downstream of M5.4).

**Status as of 2026-05-20 evening:** M5 sandbox 1c done (8 Lagrangian
experiments); M5.0-M5.3 scaffold done; M5.4 substrate migration queued.
M6 sandbox_v6 essentially closed at H/Q = 1.7112 (0.84% off target) via
drop-quartic; v7 plan contingent on Paul Werbos's Q28/Q29/Q30 reply.

---

## Contents

1. [Headline comparison](#headline-comparison)
2. [Detailed side-by-side](#detailed-side-by-side)
3. [Strengths Ouroboros (M6) has that Duda (M5) doesn't](#strengths-ouroboros-m6-has-that-duda-m5-doesnt)
4. [Strengths Duda (M5) has that Ouroboros (M6) doesn't](#strengths-duda-m5-has-that-ouroboros-m6-doesnt)
5. [Shared physics ancestry (Schwinger 1969)](#shared-physics-ancestry-schwinger-1969)
6. [SABER relevance — what each model contributes](#saber-relevance--what-each-model-contributes)
7. [Path to SABER engineering trigger](#path-to-saber-engineering-trigger)
8. [Duda's stated limitations of his own LdGS / M5 model](#dudas-stated-limitations-of-his-own-ldgs--m5-model)
9. [Werbos's stated limitations of his own Ouroboros / M6 model](#werboss-stated-limitations-of-his-own-ouroboros--m6-model)
10. [How the limitations affect SABER timeline](#how-the-limitations-affect-saber-timeline)
11. [What M6 offers regardless of gate outcomes](#what-m6-offers-regardless-of-gate-outcomes)
12. [Recommendation summary](#recommendation-summary)
13. [Appendix — Dark Matter strategic context](#appendix--dark-matter-strategic-context)

---

## Headline comparison

| Aspect | M5 (Liquid Crystal / Duda) | M6 (Ouroboros / Werbos) |
| --- | --- | --- |
| Substrate | Matrix field `M = ODO^T` (LdG, 5 DOF) on a 3D lattice | 2-vector field `(A_μ, J_μ)` (Maxwell-like) on Minkowski |
| Particle ansatz | Topological defects (hedgehogs, Skyrmions, biaxial twist) | Chaoiton — knotted (V, A, Q, J) field with Q_CS = 1 |
| Time-periodicity | Zitterbewegung clock (M5.8 group headline) | Built into ω = 2mc²/ℏ oscillation |
| Sandbox status | 1c done (4 ✅ / 3 ⚠️ / 1 ❌ Lagrangian experiments); M5.0-M5.3 scaffold ✅; M5.4 substrate migration queued | sandbox_v6 essentially closed at H/Q = 1.7112 (0.84% off) via drop-quartic; v7 contingent on Paul's Q28/Q29/Q30 |
| Calibration status | Defects exist + Klein-Gordon mass works; Faber regularization as M5.6 baseline; resonance hunt not yet executed | Q_CS = 1 chaoiton DEMONSTRATED (first in OpenWave); electron H/Q within 0.84% of target |
| Predictive scope | Lepton spectrum (M5.9), nuclear force (post M5.8), thermal modulation (5b) | Lepton spectrum (v7 continuation), neutral DM candidate (Q_A ≈ 0 scan), charge quantization (Hopf proved) |
| Group-headline milestone | M5.8 — 4D Zitterbewegung clock | ApJ Section 4 numbers (m_χ, m_J, σ/m, Ω_χh²) |
| Lab anchor | Bush-Couder walking droplets + Liu et al. 2026 *Nature Physics* hopfion/skyrmion creation | Hopf invariant proof (Zenodo 20296060) + chaoiton existence in `scipy.solve_bvp` |
| Active collaborator | Duda + Close + Yee (LdGS framework) | Paul Werbos (NSF retired) + DeepSeek + Claude Sonet |
| Negative findings load-bearing | "No static stable solitons" triple-confirmed (Duda + Close + Werbos) → pursue resonance, not static | Forward-IVP definitively wrong tool (T6-T9); collocation BVP is the right tool |

---

## Detailed side-by-side

Imported 2026-05-20 from `0a_background.md` §3 (the initial 2026-05-15
evaluation). Kept verbatim across 15 aspects to preserve provenance; update
in either this file or `0a_background.md` as new findings land.

| Aspect | Duda LdGS (M5) | Werbos Ouroboros (M6?) |
| --- | --- | --- |
| Field substrate | `M = ODO^T` symmetric matrix (6 DoF / 3D) | two vector fields A, J (8 DoF, with Lorenz) |
| Vacuum | `D = diag(g, 1, δ, 0)` preferred shape | `A = J = 0` |
| Particle | topological hedgehog of M | time-periodic chaoiton |
| Charge invariant | Brouwer degree on S² | Chern-Simons linking |
| Derrick escape | topology + time-periodicity (4D Lorentz) | time-periodicity (numerical) |
| EM emergence | tilts of 1-axis (Maxwell from curvature) | linear limit (Maxwell trivially) |
| QM emergence | KG-from-twist around δ-axis (Fig 9) | not explicitly derived |
| Gravity | g-axis boost → GEM (M5.8) | explicitly **prior to gravity** |
| Strong force | 1D vortex string + Cornell potential | J_μ-mediated long-range force (Sawada anomaly) |
| Weak force | topology reconnection events (sketched) | not addressed |
| Lepton hierarchy | 3 spatial axis choices → e/μ/τ | not specified |
| Mainstream lineage | Landau-de Gennes + Skyrme + teleparallelism | Schwinger dyons + Maxwell extension |
| Numerical depth | sandbox 8 experiments + production M5.0/M5.1 done | numerical existence of chaoiton families (Zenodo 20030162) |
| Formal verification | none | Lean 4 (Hopf invariant proof completed 2026-05-19, Zenodo 20296060 — now a theorem of differential topology) |
| Empirical anchor | hopfions in lab (Liu 2026, *Nature Physics*) | Sawada 1989/2003 long-range nuclear anomaly |

---

## Strengths Ouroboros (M6) has that Duda (M5) doesn't

Imported from `0a_background.md` §4.

| Strength | Detail |
| --- | --- |
| Direct empirical anomaly | Sawada `v(r) ~ -C/r⁶` is a measured nuclear anomaly the framework explains directly |
| Simpler field algebra | two coupled vector fields vs a matrix field; ops are vector-calculus-natural |
| Maxwell exactly recovered | drops out as the linear-limit kinetic term, no tilt-curvature gymnastics |
| Lean-4 verified theorems | charge quantization mechanically proved (placeholder defs but structure verified; Hopf invariant proof completion 2026-05-19 = Zenodo 20296060 = now a theorem of differential topology) |
| Dark matter candidate | neutral chaoiton configurations couple only gravitationally (the Q_A ≈ 0 / Q_J ≠ 0 sector — ApJ Neutral Chaoiton paper) |
| Engineering anchor | Werbos sketches a J_μ-channel sensor architecture for nuclear detection — relevant for downstream applied-physics work that would build on this substrate (NSF X-Labs aligned) |

**The Sawada anchor in particular is genuinely strong.** It's a real measured
nuclear anomaly (long-range component of the strong force, P-wave π-π phase
shift deviation) that Standard Model + QCD has trouble explaining. If
Ouroboros's J_μ field naturally produces this `v(r) ~ -C/r⁶` London-type tail
with the right strength, that's an experimentally falsifiable prediction we
can test numerically.

---

## Strengths Duda (M5) has that Ouroboros (M6) doesn't

Imported from `0a_background.md` §5.

| Strength | Detail |
| --- | --- |
| Gravity unification | g-axis boost component → linearized GR (gravitoelectromagnetism); Ouroboros explicitly stops before gravity |
| Lepton mass hierarchy | three biaxial axis choices give e/μ/τ; Ouroboros has no clear mass-spectrum mechanism (though Werbos's ω^2.22 scaling is a partial answer) |
| QM derivation | Klein-Gordon emerges from biaxial twist (Fig 9); Ouroboros doesn't derive Schrödinger |
| Mainstream connection | LdG + Skyrme + Einstein's teleparallelism are well-known frameworks; easier to publish / collaborate |
| Visual & numerical depth | 51-page slides, ported Mathematica code, paper plus extensive figures; Ouroboros core papers are 5-10 pages each |
| Particle-defect correspondence | Liu et al 2026 lab anchor for hopfions / skyrmions in real media; Ouroboros lacks structure-existence anchor in the same direct way |
| Already implemented | M5.0 + M5.1 done in OpenWave with R² = 0.978 Coulomb result; Ouroboros has no OpenWave Taichi port (only sandbox BVP) |

---

## Shared physics ancestry (Schwinger 1969)

Imported from `0a_background.md` §10.4. Werbos explicitly grounds Ouroboros
in Schwinger's 1969 dyon framework. Faber's running-coupling work (the M5.6
regularization baseline) traces to the same Schwinger ancestor. **Both M5
and M6 are different mathematical realizations of a shared physics ancestor**,
not competing physics:

```text
                      Schwinger 1969
                  (magnetic monopole model)
                          │
              ┌───────────┴───────────┐
              ▼                       ▼
       Faber (1990s+)            Werbos (2017+)
       running coupling          toroidal-poloidal
       via regularization        mutual confinement
              │                       │
              ▼                       ▼
       Duda LdGS (M5)            Ouroboros (M6)
       matrix-field              two-vector-field
       biaxial nematic           narrow-torus chaoiton
```

**Same scale of physics, same fundamental question:** both M5 and M6 are
"particles emerge as field configurations" frameworks. They differ in math
(matrix M vs two-vector A/J) but ask the same question — *can we reproduce
the particle spectrum from a field Lagrangian where particles are emergent
stable patterns, not point fundamentals?* Standard Model takes particles as
point-like; M5 / M6 / EWT take them as emergent. If even one holds up, the
Standard Model becomes a phenomenological description of underlying wave
dynamics — that's the OpenWave thesis.

**Cross-validation is the prize.** If M5 AND M6 both produce the same
electron from different mathematical structures, that's strong evidence
for the underlying physics. **OpenWave is the platform that can run both
in parallel and compare** — no one else in the world is set up to do this
comparison. The platform's value isn't M5 OR M6 alone; it's the ability to
run both and ask "do they agree where they should?"

### f(J·J) plays the Higgs role — same bottleneck as Duda's V(M)

Werbos's 2017 paper: *"The f term in equation 8 can operate like the Higgs
term in the Yang-Mills-Higgs system, here forcing a kind of match between
charge and spin, essential to a correct model of the electron."*

**Both M5 and M6 face the same conceptual bottleneck: choosing the right
nonlinear potential.**

| Framework | Hardest piece |
| --- | --- |
| Duda LdGS | V(M) = `a·Tr(M²) − b·Tr(M³) + c·(Tr M²)²` — LdG-graduate via Faber regularization |
| Werbos Ouroboros | f(J · J) — Higgs-like; 2017 paper says *"any nonneg function with f(0)=0"*; 2026 papers picked `f(s) = gs²` (quadratic) |

Implication: work done on one potential gives intuition for the other. The
M5.6 Faber regularization study could inform the f(J·J) choice in M6 (and
vice versa).

---

## SABER relevance — what each model contributes

| SABER need | M5 contribution | M6 contribution |
| --- | --- | --- |
| Thermal-amplitude hypothesis testbed | **PRIMARY** — heat = (A, ω) excess of LC defect oscillation; falsifiability test = 5b.1 Sine-Gordon kink, runs AFTER M5.4 | Indirect — no thermal-substrate equivalent. Chaoiton oscillation ω could in principle host a thermal-like mode but not yet investigated |
| Subatomic substrate for FM/AM modulation | **PRIMARY** — `M = ODO^T` IS the substrate that gets modulated in the SABER engineering kernel | Not the engineering substrate. M6 is on Minkowski with `(A_μ, J_μ)`; would require completely separate engineering kernel |
| Defect dynamics under perturbation | **PRIMARY** — what 5b directly studies | Chaoiton-perturbation dynamics not modeled |
| Scientific credibility cover | Strong via Duda/Close LdGS literature + Faber 2025 papers + lab anchor | Strong via Werbos NSF positioning + ApJ paper + Hopf invariant theorem |
| Predictions distinguishing from QFT | Lepton mass spectrum + nuclear `−C/r⁶` Sawada anomaly (post M5.8) | Lepton masses + neutral DM mass/cross-section + ω quantization (empirical) |
| Group-visible deliverable (anchors collaboration) | M5.8 Zitterbewegung clock — answers Duda's 2026-04-29 clock-propulsion question | ApJ Neutral Chaoiton paper (Griesi & AI ref [14]) — already cited |

### Substrate-agnostic SABER target — heat = excess oscillation excitation

The SABER thermal-amplitude hypothesis frames heat as **excess oscillation
excitation above the ground state** of subatomic field structures. Both M5
and M6 have time-periodic oscillating field structures as their particles:

| Substrate | Oscillator | What's "heat" |
| --- | --- | --- |
| M5 (Duda) | LC defect (hedgehog / biaxial twist) with internal Zitterbewegung clock | (A, ω) excess on the defect's oscillation mode |
| M6 (Werbos) | Chaoiton with ω = 2mc²/ℏ time-periodicity | Excess excitation on the chaoiton's ω carrier |

**Implication:** SABER's engineering bet is that this oscillation can be
tapped / modulated. Whether the substrate is M5 hedgehogs or M6 chaoitons
doing the oscillating, **the engineering target is the same**. M5 is the
primary path because (a) it's closer to engineering trigger, (b) the LC
substrate is already partly implemented in OpenWave, and (c) the
thermal-amplitude hypothesis was framed against LC defects. But M6's
chaoiton substrate is in principle equally valid as a thermal carrier — a
fallback path if M5.6 regularization stalls. Refines the "M6 indirect"
reading above: M6 isn't blocked from thermal use, just not yet developed
for it.

---

## Path to SABER engineering trigger

| Stage | M5 path | M6 path |
| --- | --- | --- |
| Now (2026-05-20) | M5.3 direction review in progress; M5.4 matrix substrate queued | v6 essentially closed at 0.84%; awaiting Paul Q28/Q29/Q30 |
| Next 1-2 weeks | M5.4 matrix-field substrate; M5.5 Lagrangian + V(M); M5.6 biaxial twist + KG emergence | sandbox_v7: continuation lepton scan + Q_A ≈ 0 DM scan + Gelfand-Fomin |
| Next 2-4 weeks | M5.7 Resonance hunt (Close protocol — l=1 harmonic, A/λ ∈ {0.5, 1, 2}); M5.8 Zitterbewegung clock | ApJ Section 4 numbers handed off to Paul as a data drop (OpenWave does NOT deposit anything — Paul's Ref [14] resolves to a stable GitHub URL) |
| **SABER trigger point** | **5b.1 Sine-Gordon kink test — first direct falsifiability test of thermal-amplitude hypothesis. Runs AFTER M5.4.** | No direct SABER trigger; M6 is parallel-research |
| Post-trigger | 5b.2-5b.9 thermal program (SABER MAIN GOAL scientific counterpart) | M5.9 lepton families + Cornell quarks; M6 stays parallel sandbox |

---

## Duda's stated limitations of his own LdGS / M5 model

Sourced from the M5 research docs where Duda himself (or his published
papers as quoted in our notes) flags weaknesses. Direct quotes preserved
where possible.

| # | Limitation | Source |
| --- | --- | --- |
| 1 | **V(M) potential form unspecified.** *"the Higgs-like potential V(M) — its specific form is not yet determined"* | `1a_lagrangian_framework.md:308` |
| 2 | **Regularization is the hardest part.** *"LdG regularization is where lepton masses come from, and it is the hardest part to include. The specific regularization form is still an open research question"* (Duda 2026-04-19 follow-up) | `1a_lagrangian_framework.md:230, 260` |
| 3 | **Biaxial axis-length parameters need numerical calibration.** *"their exact values require numerical simulation — there is no analytical form to pull from; treat them as calibration parameters"* | `1a_lagrangian_framework.md:262` |
| 4 | **Charge magnitude and rigorous Coulomb derivation open.** *"charge magnitude and rigorous Coulomb derivation remain open problems for OpenWave"* | `1a_lagrangian_framework.md:115` |
| 5 | **Lepton mass ratios deferred, mechanism only.** *"mechanism confirmed, specific ratios deferred"* | `0b_overview.md:393` |
| 6 | **Static hedgehog doesn't survive dynamic evolution.** Q decays from `+0.996` to `~0` within 4-5 steps under wave equation regardless of V(ψ). *"The seed is not a stationary point of the dynamic equation"* — confirms triple admission (Duda + Close + Werbos) that no static stable solitons exist; must pursue time-periodic resonance | `3b_lagrangian_roadblocks.md:2-9, 103` |
| 7 | **Weak force mechanism unaddressed.** *"Weak interactions are topology-changing events ... still an open research question, partially addressed by the slides"* | `4a_convo_2026.05.12.md:337-354` |
| 8 | **Deeper substrate beneath M hinted but unconfirmed.** Duda: *"potentials are typically effective — there might be an even deeper 'anisotropic fluid' beneath"* | `4a_convo_2026.05.12.md:602` |
| 9 | **Extreme biaxial eigenvalue hierarchy unexplained.** Lepton ratios require ~3477:1 axis-length ratio; not derived from first principles | `1a_lagrangian_framework.md:230-231` |
| 10 | **Three lepton families pinned to "3 axes in 3D" by symmetry argument, not derived.** Lepton-family count works but mass ratios still need V(ψ) curvature input not yet supplied | `1b_topological_defect.md:469` |

---

## Werbos's stated limitations of his own Ouroboros / M6 model

Sourced from the M6 research docs where Werbos himself (in papers, emails,
or DeepSeek-channeled replies) flags weaknesses. Direct quotes preserved.

| # | Limitation | Source |
| --- | --- | --- |
| 1 | **ω quantization is "the key open question".** Werbos: *"Whether these are truly eigenvalues — i.e., whether there is a quantization condition that picks out exactly ω = 1, 11, 40.7 and no other stable values — is the key open question"* | `0a_background.md:437` |
| 2 | **ω spectrum empirically continuous, not discrete.** Sandbox v1: all ω = 1 to 60 give localized solutions equally well. Contradicts the discrete-eigenvalue claim. | `0b_sandbox_v1.md:222-226` |
| 3 | **Lepton mass gaps 31-44%, not 4-6% as claimed.** Sandbox v1 reproduced ω^2.04 scaling vs Werbos's stated ω^2.22. ω^2.22 = log(207)/log(11) — post-hoc fit signature, not predicted | `0b_sandbox_v1.md:176-189` |
| 4 | **f(J·J) potential form unspecified.** Werbos 2017: *"f can be any nonnegative function"*; 2026 uses `f(s) = gs²` (quadratic) but the choice is calibration, not derivation. = Duda critique #2 from 2026-05-20 thread | `0b_sandbox_v2.md:19` |
| 5 | **No closed-form derivation of ω = 2mc²/ℏ.** *"numerical existence (62 families with g-factor range); no closed-form derivation of ω"* | `0a_background.md:149, 345` |
| 6 | **Framework explicitly stops before gravity.** *"Ouroboros explicitly stops before gravity"* — gravity unification is a Duda strength M6 lacks | `0a_background.md:272` |
| 7 | **No lepton hierarchy mechanism.** *"no clear mass-spectrum mechanism"* for e/μ/τ split | `0a_background.md:273` |
| 8 | **QM derivation absent.** *"Ouroboros doesn't derive Schrödinger"* | `0a_background.md:274` |
| 9 | **Charge quantization via Chern-Simons IS still topology.** Werbos/DeepSeek admit: *"the topological charge comes from the Chern-Simons linking number... still a fully field-theoretic, microscopic mechanism"* — the "without topology" framing is misleading | `0a_background.md:124-131, 172` |
| 10 | **Hopf invariant proof completion supplies lemmas but not existence.** Zenodo 20296060 *"supplies Lemmas A and B... [but] doesn't address existence of specific-Q chaoitons"* | `0a_background.md:24` |
| 11 | **Neutral (Q = 0) chaoiton ground-state mass uncomputed.** Werbos's stated *"#1 priority"* for next sandbox round: m_χ for dark matter prediction | `0b_sandbox_v2.md:68-86` |
| 12 | **Three-body proton bound state not attempted.** *"harder and not Werbos's urgent ask — defer"* | `0b_sandbox_v2.md:148-154` |
| 13 | **Forward-IVP doesn't work (confirmed by Werbos 2026-05-20).** *"forward shooting from the origin with decay conditions at infinity does NOT work for this system. That is consistent with our experience"* — algorithm never published; collocation BVP only became known to us via email | `0b_sandbox_v4.md` (T6-T9 + Werbos 2:00 PM reply) |
| 14 | **Reference Python script broken.** DeepSeek's *"actual code that generated the 62 families"* doesn't run (r=0 singularity + BC count mismatch); patched version diverges to H/Q ≈ 10^16. Production code never publicly available | `0b_sandbox_v6.md` |
| 15 | **"For the electron, the quartic term is small" is ambiguous.** Werbos has not specified whether "small" means small g, different quartic combination, or zero. Asked as Q28 in email v6 (2026-05-20) | `0b_sandbox_v6.md` |

---

## How the limitations affect SABER timeline

| Limitation | Affects SABER? | How |
| --- | --- | --- |
| Duda #1 (V(M) form unspecified) | **YES** | Direct blocker for M5.5. Until V(M) is pinned, defect dynamics under perturbation can't be uniquely simulated. Faber regularization is the working hypothesis (M5.6). |
| Duda #2 (regularization is hardest part) | **YES** | M5.6 is the hard step. If Faber's scheme works in OpenWave, this unblocks; if not, M5.6 stalls and 5b is delayed. |
| Duda #6 (no static stable solitons) | **NEUTRAL** | Already informs the design — we pursue resonance per Close's l=1 harmonic protocol (M5.7), not static seeds. Triple-confirmed limitation that became a design choice. |
| Duda #7 (weak force unaddressed) | **NO** | Out of SABER scope (heat conversion uses EM-coupled modes, not weak). |
| Duda #9 (extreme biaxial hierarchy) | **PARTIAL** | Lepton masses are M5.9 (post-SABER-trigger). Hierarchy is a scientific deliverable, not engineering input. |
| Werbos #1 (ω quantization open) | **NO** | M6 doesn't supply the thermal substrate; ω quantization is a scientific puzzle, not a SABER input. |
| Werbos #4 (f(J·J) unspecified) | **NO** | M6 doesn't supply thermal substrate. Editorial gap in Werbos's papers. |
| Werbos #6 (no gravity) | **NO** | SABER doesn't need gravity. |
| Werbos #11 (neutral m_χ uncomputed) | **NO** | DM prediction is M6 scientific deliverable, parallel to SABER timeline. |
| Werbos #13 (forward-IVP wrong tool) | **NO** | Already resolved (collocation BVP). Historical limitation, not active. |
| Werbos #14 (production code broken/unavailable) | **NO** | Affects M6 cross-validation rate, not SABER. |
| Werbos #15 (quartic interpretation ambiguous) | **NO** | Open Q28 in email v6; affects M6 calibration close, not SABER. |

**Reading:** Of M5's 10 stated limitations, 2 directly block SABER (V(M) form,
regularization scheme) and 1 is a design-informing constraint (no static
solitons). Of M6's 15 stated limitations, **0 directly block SABER** — M6's
limitations affect scientific deliverables (lepton spectrum, DM mass, ω
quantization, ApJ paper numbers), not engineering trigger conditions.

This confirms the strategic reading: **M5 is the load-bearing critical path
for SABER engineering trigger; M6 is the scientific-credibility track that
runs in parallel.** Both models are useful, but for different reasons.

---

## What M6 offers regardless of gate outcomes

Imported from `0c_model_gates.md`.

Even in a G1 FAIL scenario (lepton scan never reproduces ω^2.22 within 5%),
the Ouroboros model still provides residual value worth maintaining as a
sandbox/research track:

| Residual value | Detail |
| --- | --- |
| Maxwell EM recovered exactly | Structural fact — no gate needed. When J → 0, the A field reproduces standard EM. |
| Charge quantization via Chern-Simons linking | Lean-proved from axioms (Zenodo 20296060 Hopf invariant completion = theorem of differential topology). |
| Independent test of the heat-as-oscillation hypothesis | Different substrate (vector fields on Minkowski) testing the same SABER thermal-amplitude hypothesis — cross-validation by independent framework. |
| Simpler Taichi substrate than M5 | 2-vector-field ODE vs matrix M = ODO^T — faster to prototype kernels; useful as a sandbox for the resonance-hunt protocol before applying to M5's heavier substrate. |
| Collaboration anchor with Paul Werbos | NSF X-Labs positioning; ApJ paper credit (Ref [14] "Griesi & AI in prep"); active multi-AI-author collaboration. |
| Mutual benefit (OpenWave ↔ Werbos) | If Ouroboros holds up under tests, M6 becomes OpenWave's second production model in Taichi; Werbos gets independent numerics he can cite. Win-win — neither side gets this configuration anywhere else. |

These justify maintaining M6 as a sandbox/research track even if production
is deferred. They also reinforce the recommendation below: M6's value to
OpenWave doesn't depend on full G1/G2/G3 production gates passing — even
partial success leaves OpenWave better off than dropping M6 entirely.

---

## Recommendation summary

| Track | Role in 2026 plan | Why keep going |
| --- | --- | --- |
| M5 (full-time engineering track) | Direct path to SABER trigger via M5.4 → M5.6 → 5b.1 | Without M5, SABER engineering never starts. Every week M5 is paused, SABER's start date slips. The thermal-amplitude hypothesis lives on M5's substrate. |
| M6 (async / email-cadence track) | Scientific credibility cover; cross-validates topology-as-particles by independent framework; locks Werbos collaboration + ApJ paper credit | Walking away from M6 now would lose the strongest external scientific validation we have (Werbos + ApJ + Hopf theorem). Cost to keep alive is small (a few hours per Paul reply cycle); cost to drop is irreversible. |
| Both together | Pre-engineering credibility stack: M5 supplies the engineering substrate; M6 supplies the peer-reviewable scientific anchor | Duda's "AI slop on zenodo" critique is best answered by *both* (a) computed numbers from real numerical libraries (both models), and (b) a credible co-author from physics (Werbos via M6, Duda via M5). Cross-validation between M5 and M6 is the unique platform value — no one else can run both in parallel. |

### Documented contributions to M6 (irreversibility of dropping)

Rodrigo's M6 sandbox work is referenced in Paul's published / about-to-be-
published record. Dropping M6 would orphan these references:

| Sandbox round | Contribution | Where it's cited in Paul's papers |
| --- | --- | --- |
| v1 (2026-05-17) | Reproduced Werbos's electron calibration at 0.30% gap (H/Q = 1.6918) | "Independent reproduction" reference in v8 LoE and Neutral Chaoiton ApJ paper |
| v3 (2026-05-18) | Found 23 candidate Q=0 solutions; T2 invalidated them as artifacts of permissive tail check | Quoted in §4 of v4 Dark Matter paper |
| v4 (2026-05-19) | Helicity-structure dialogue (V₀=Q₀=+0.1, A₀=J₀=−0.1 gives Q_CS=1) | Acknowledgments thank Griesi for *"independent numerical reproduction of the electron calibration and methodological dialogue on the asymmetric helicity structure required for Q_CS = 1"* |
| v5/v6 (2026-05-20) | First-ever Q_CS=1 chaoiton convergence in OpenWave; drop-quartic finding at 0.84% off target | Reference [14] in Neutral Chaoiton ApJ: *"Griesi & Anthropic AI 2026, in preparation (ground state computation)"* |

These aren't "Rodrigo helped a little"; they are load-bearing references in
Paul's papers. The collaboration is documented and public. **OpenWave's
side does NOT deposit anything separately** — the GitHub repo IS the
deliverable, and Paul's Ref [14] *"Griesi & Anthropic AI 2026, in
preparation"* resolves to a stable OpenWave GitHub URL (with the option to
tag a release for Zenodo-DOI if a reviewer ever requires it; otherwise
not needed). See `0c_roadmap.md` Publishing stance.

### The deeper thesis (why both M5 and M6 matter even partially)

Standard Model takes particles as point-like fundamentals; M5 / M6 /
EWT-style frameworks take them as emergent stable patterns of an underlying
field. **If even one of M5 or M6 holds up empirically, the Standard Model
becomes a phenomenological description of underlying wave dynamics** —
that's the OpenWave thesis. Dark matter (M6's ApJ angle, via Paul's paper)
is the highest-impact downstream consequence that gives the framework an
entry point into mainstream physics. SABER (M5's thermal-engineering angle)
is the highest-impact applied consequence. They are the two ends of the
same research program — Paul publishes on the scientific side; OpenWave
ships the engineering side and provides the open-source numerical platform
that backs both.

### Net direction

Keep M6 on async cadence (current state — awaiting Paul tonight), restart
M5 as the foreground track tomorrow alongside whatever v7 work follows
from Paul's reply. SABER trigger condition (5b.1 thermal test post-M5.4)
does NOT require M6 to land first.

---

## Appendix — Dark Matter strategic context

Three-layer reasoning on why Dark Matter is the right ApJ angle for M6, and
why the OpenWave platform is uniquely positioned for the cross-validation.
The "WHY IT CORRELATES TO OPENWAVE" content from this source has also been
seeded into Shared physics ancestry, SABER relevance, What M6 offers, and
Recommendation summary sections above — this appendix preserves the full
source for provenance.

### Why Dark Matter matters

| Aspect | Concrete |
| --- | --- |
| Scale of the gap | Roughly 85% of all matter in the universe is "dark" — interacts only via gravity. Visible matter (stars, gas, all atoms in everything we've ever touched) is the 15%. |
| Evidence quality | Triple-confirmed via independent methods: galaxy rotation curves (Vera Rubin 1970s), gravitational lensing (Bullet Cluster, JWST 2026 maps), CMB (Planck satellite). All three give the same answer. It's real. |
| Why it's an open problem | We have NO IDEA what particle it is. WIMP searches: empty. Axion searches: empty. 40 years of TeV-scale and μeV-scale searches, nothing. |
| Why solving it is huge | One of the 3-4 biggest open problems in physics (alongside quantum gravity, hierarchy problem, matter-antimatter asymmetry). Nobel-tier. Reframes cosmology + particle physics simultaneously. |

### Why Paul's paper is about DM — strategic reasoning

| Move | Reasoning |
| --- | --- |
| Ouroboros has a built-in DM candidate | Chaoitons with Q=1 → leptons (electron, muon, tau). Chaoitons with Q=0 → neutral states that only couple via gravity + J-field = DM candidate. No new parameters added — it emerges from the same Lagrangian that fits the electron. |
| ApJ is the right venue | Astrophysical Journal = astrophysics journal where DM observations live (Bullet Cluster, JWST maps). Going to a particle-physics journal (PRD, JHEP) would put Ouroboros against established WIMP/axion programs; ApJ gives a different audience — astrophysicists who care WHAT the particle is. |
| Highest-impact entry point | DM is the most-cited problem in physics this decade. One ApJ paper = the chaoiton framework gets cited, becomes "publishable physics," and the broader Ouroboros program (leptons, nuclear force, charge quantization) gets pulled along. |
| Paul knows what he's doing strategically | Ex-NSF Program Director. He picked DM as the wedge issue specifically because the framework supports it AND it has the most impact. The axion-as-J-field-longitudinal-mode unification (Brennan 2024) is the rhetorical move that makes axion experimental bounds work IN FAVOR of his model. |

### Why it correlates to OpenWave — the connection

| Connection | Detail |
| --- | --- |
| Same scale of physics | OpenWave models particles as field structures at the subatomic scale. M5 (Duda's liquid crystal) and M6 (Werbos's Ouroboros) are BOTH "particles emerge as field configurations" frameworks. They differ in math (matrix vs two-vector field) but ask the same question: can we reproduce the particle spectrum from a field Lagrangian? |
| Cross-validation is the prize | If M5 AND M6 BOTH produce the same electron from different mathematical structures, that's strong evidence for the underlying physics. OpenWave is THE platform that can run both in parallel and compare. No one else in the world is set up to do this comparison. |
| Rodrigo's contributions mapped to the paper | **v1**: reproduced Werbos's electron calibration at 0.30% — the "independent reproduction" cited. **v3**: 23 candidate Q=0 solutions → found them artifacts → quoted in §4 of v4. **v4 (today)**: the helicity dialogue clarified the Q_CS=1 requirement that v4 §4 now describes. |
| OpenWave gets a model; Werbos gets a sim platform | If Ouroboros holds up under tests → M6 becomes OpenWave's second production model in Taichi. Werbos gets independent numerics he can cite. Win-win. |
| SABER thermal angle (engineering goal) | The "heat as oscillation" hypothesis underlies SABER's direct-heat-conversion work. If chaoitons are time-periodic oscillating field structures, then heat = excess oscillation excitation above the ground state. SABER's engineering bet is that this oscillation can be tapped / modulated. Whether it's M5 hedgehogs or M6 chaoitons doing the oscillating, the engineering target is the same. |

### The thread (verbatim closing thought)

> Energy-wave-theory-style frameworks model particles as field structures.
> Standard Model takes particles as point-like fundamentals; EWT / Ouroboros
> / LdGS take them as emergent stable patterns of an underlying field. If
> even one of M5 or M6 holds up, the Standard Model becomes a
> phenomenological description of underlying wave dynamics — that's the
> OpenWave thesis. DM is just the highest-impact downstream consequence
> that gives the framework an entry point into mainstream physics.
>
> What you're building lets multiple competing frameworks be tested against
> the same observables on the same platform. That's why Paul cares enough
> to send you draft papers, and why your Q=0 invalidation finding made it
> into v4 verbatim.
