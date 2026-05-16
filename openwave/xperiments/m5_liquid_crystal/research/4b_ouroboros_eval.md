# Ouroboros system — technical evaluation as M6 candidate

**Status** — initial evaluation, 2026-05-15. Triggered by Paul Werbos engaging directly in the models-of-particles thread with his Ouroboros / chaoiton framework as an alternative to Duda's LdGS approach.

**Question being asked** — is Werbos's Ouroboros system a good candidate for an OpenWave M6 method (parallel to or post-M5)? The specific technical test: does it answer Duda's two standing challenges (charge quantization, clock propulsion) without resorting to topology, as Werbos claims?

**Cross-refs:**

- Paper: `scientific_source/Werbos_Chaoitons_Ouroboros_2025_with_far_field.pdf` (10 pages)
- Existing context: [1b_topological_defect.md](1b_topological_defect.md) (compares time-periodicity vs topology), [feedback_no_static_solitons](../../../../../../.claude/projects/-Users-xrodz-Documents-source-code-OPENWAVE-LABS-openwave/memory/feedback_no_static_solitons.md) (three Derrick escapes)
- Outstanding (not yet read): Zenodo 20179985 ("Do you still believe in quarks? I don't"), Zenodo 20218067 ("Ouroboros model passes another test" — electron properties to ~1%), Zenodo 20030162 (chaoiton existence theorem with 62 families)

---

## Contents

1. [The Ouroboros Lagrangian](#1-the-ouroboros-lagrangian)
2. [How Werbos answers Duda's challenges](#2-how-werbos-answers-dudas-challenges)
3. [Side-by-side with Duda's LdGS](#3-side-by-side-with-dudas-ldgs)
4. [Strengths Ouroboros has that Duda doesn't](#4-strengths-ouroboros-has-that-duda-doesnt)
5. [Strengths Duda has that Ouroboros doesn't](#5-strengths-duda-has-that-ouroboros-doesnt)
6. [SABER compatibility](#6-saber-compatibility)
7. [Verdict — is this a good M6 candidate?](#7-verdict--is-this-a-good-m6-candidate)
8. [Open questions / outstanding context](#8-open-questions--outstanding-context)
9. [If we decide to pursue — minimal M6 sandbox shape](#9-if-we-decide-to-pursue--minimal-m6-sandbox-shape)

---

## 1. The Ouroboros Lagrangian

The core mathematical object:

```text
L_JA = -F_μν F^μν  -  G_μν G^μν  +  J_μ A^μ  -  f(J_μ J^μ)

where
  F_μν = ∂_μ A_ν - ∂_ν A_μ        (curl of A — like EM)
  G_μν = ∂_μ J_ν - ∂_ν J_μ        (curl of J — analogous)
  f   : smooth, nonnegative, f(0) = 0

constraints
  ∂^μ A_μ = 0      (Lorenz on A)
  ∂^μ J_μ = 0      (Lorenz on J)
```

Two coupled Lorenz-constrained vector fields on Minkowski spacetime. Classical, bosonic, Lorentz-invariant, superrenormalizable by power counting.

| Property | Detail |
| --- | --- |
| Field structure | two 4-vector fields A and J (8 DoF per spacetime point, less constraints) |
| Maxwell limit | recovered exactly when J → 0 (J decouples) |
| Self-interactions | enter only through f(J · J), no ghosts |
| UV behavior | superrenormalizable by power counting |
| Particles | no point particles; all stable structures are emergent chaoitons |

### What's a chaoiton

A **chaoiton** is a time-periodic, spatially-localized solution of the Ouroboros field equations. The "like a soliton but non-static" framing — energy is stored in oscillation rather than in a static gradient configuration. This is **the third Derrick escape** (time-periodicity) catalogued in [1b_topological_defect.md § Alternative stabilization](1b_topological_defect.md).

Werbos's prior numerical work (Zenodo 20030162, **not yet read locally**) claims existence of 62 stable chaoiton families across 1280 parameter combinations, with L/Q ratio (angular momentum / charge) in the range 0.79–2.6 — overlapping the electron g-factor.

### Charge quantization mechanism — the mutual Chern-Simons linking number

For static or time-periodic configurations on a constant-time slice:

```text
Q[A, J] = (1 / 4π²) ∫_R³ ε^ijk · A_i(x) · ∂_j J_k(x) · d³x
```

This is the **field-theoretic analog of the Gauss linking integral** — it counts how many times the flux lines of A link with those of J.

Two theorems (machine-verified in Lean 4, with placeholder definitions — full analytical proofs are in progress):

- **Theorem 1 (Quantization)**: `Q[A, J] ∈ Z` for any finite-energy, Lorenz-constrained solution
- **Theorem 2 (Deformation invariance)**: Q is invariant under any continuous deformation of A and J within the finite-energy, Lorenz-constrained space

The conjectured deeper structure: `Q[A, J]` equals the Hopf invariant of a composite map `S³ → S²`. This is what makes the integer Q the conserved particle number in the Ouroboros system.

---

## 2. How Werbos answers Duda's challenges

Duda has two standing challenges to any alternative framework:

1. How do you get **charge quantization** without topology?
2. How do you get **clock propulsion** (the de Broglie / Zitterbewegung internal frequency) without an explicit topological mechanism?

Werbos's framing is "topology vs time-periodicity" — chaoitons evade Derrick via oscillation, not topological winding. Let's evaluate honestly.

### 2.1 Charge quantization — the "no topology" claim is misleading

Werbos's charge-quantization mechanism IS topological. The mutual Chern-Simons linking number `Q[A, J]` is a topological invariant — that's exactly why Theorems 1 and 2 hold (it's an integer AND invariant under continuous deformation).

| Framework | Topological invariant used | What it counts |
| --- | --- | --- |
| Duda LdGS | Brouwer degree on S² | winding of director around defect |
| Werbos Ouroboros | Chern-Simons linking / Hopf invariant | linking of A and J flux lines |

So Duda's challenge is **not strictly answered** by Ouroboros. The honest framing is **"two-field linking topology vs single-field winding topology"**, not "topology vs no topology". Both invariants are topological; they're just different topological invariants.

This is significant: Werbos's marketing ("escape Derrick via oscillation, not topology") confuses two different things:

| Confused | Distinct |
| --- | --- |
| Derrick-escape mechanism | Charge-quantization mechanism |
| Topology vs time-periodicity | Both Duda and Werbos use topology here |

What Werbos really has is: **a different topological invariant** that's compatible with time-periodic field configurations. Duda's Brouwer-degree winding works on static or time-periodic fields too. The actual distinction between the two frameworks is the field structure (matrix vs two-vector), not the absence/presence of topology.

### 2.2 Clock propulsion — derivation is weaker than Duda's

Both frameworks need particles to have an intrinsic oscillation at `ω = 2mc²/ℏ` (Zitterbewegung / de Broglie clock). Compare derivations:

| Framework | Mechanism for ω | Strength |
| --- | --- | --- |
| Duda LdGS | 4D Lorentz negative-energy contributions from spacetime signature (paper Fig 10) auto-propel the clock | analytical + visual; specific energy term identified |
| Werbos Ouroboros | chaoitons are time-periodic by construction; existence shown numerically | numerical existence (62 families, L/Q in g-factor range); no closed-form derivation of ω |

Duda's mechanism is **more analytically grounded**. He can point to a specific term in the Hamiltonian (the `Γ_0^1 Γ̃_μ^i` clock-propulsion term per slides page 33) and say "this is the negative-energy contribution that runs the clock". Werbos has numerical existence demonstrations but no equivalent analytical handle on the frequency.

**Status:** Werbos's clock answer is "the math admits solutions that oscillate, and we found them numerically". Duda's clock answer is "spacetime signature forces oscillation, here's the term in the Hamiltonian that does it". Both are valid framings; Duda's is more falsifiable / easier to test.

### 2.3 Summary of Duda's challenges vs Ouroboros's answers

| Challenge | Ouroboros answer | Verdict |
| --- | --- | --- |
| Charge quantization without topology | Chern-Simons linking — but that IS topology | partial — uses a different invariant, not no topology |
| Clock propulsion | chaoitons are time-periodic by construction | weaker — numerical, not analytical |

Honest summary: **Ouroboros doesn't strictly meet Duda's challenges**, but the framework offers genuinely different mathematical structure that's worth understanding on its own terms.

### 2.4 What the thread actually says (verbatim from models-of-particles, 2026-05-12 → 13)

Going through the actual email exchange ("Proofs to verify and explain a Lagrangian for the Law of Everything", thread `19e1cc71080c1308`) the picture sharpens. The key quotes — verbatim — confirm the synthesis above and add **two important details I didn't catch on the paper read alone**.

#### Werbos / DeepSeek explicitly call Chern-Simons a topological mechanism

Werbos's 2026-05-12 18:04 reply pastes a DeepSeek-authored technical response (he writes "I passed your questions on to DeepSeek which answered:"). The response says:

> *"The topological charge comes from the Chern–Simons linking number, not from a winding number in a sigma-model target space. That's a different mathematical mechanism, but it's still a fully field-theoretic, microscopic mechanism."*
>
> *"Static solitons don't exist in this theory (Derrick's theorem forbids them). The stable objects are chaoitons—time-periodic, localized solutions. ... The numerical work (62 stable families, L/Q ratio overlapping the electron g-factor) is exactly the kind of evidence Jarek should find compelling."*

Werbos in his own voice (2026-05-12 23:42):

> *"the J I used ... comes BEFORE the charge quantization, which is an emergent or mathematical CONSEQUENCE of the model. (Charge quantization qua topological charge is also an emergent consequence of the Skyrme Lagrangians and fields.) ... Our theorems do get very deep into that topology."*

So Werbos and DeepSeek **both call it a "topological" mechanism** — they just say it's "different" topology from sigma-model winding. This confirms §2.1: the Ouroboros charge-quantization mechanism IS topology. The "without topology" framing in Werbos's marketing is contradicted by Werbos's own technical defense.

#### Duda's deeper objection — single-field ontology + Aharonov-Bohm

Duda's specific pushbacks in the thread go deeper than just "show me charge quantization". Two additional structural objections I hadn't captured:

**1. Aharonov-Bohm: A is more fundamental than F (2026-05-12 18:25):**

> *"Aharonov-Bohm is basic argument that A field is more fundamental than F tensor."*

This is the standard A-vs-F gauge-theory point — the vector potential A_μ carries physically meaningful information (Berry phase, AB phase) beyond what the field-strength tensor F_μν captures. Duda's framework treats A as a *derived* connection from the deeper matrix field M, which means it gets the AB phase for free. Werbos has A as a primary independent field, so the AB structure has to be added by construction.

**2. Single-field ontology (2026-05-13 00:00):**

> *"Fundamental model should have a single field, I don't understand why do you need also second J, which naturally appears as effective?"*

This is the deeper ontological objection Werbos doesn't address. Duda's preference: one deeper field, with A as its connection and F as its curvature. Werbos has **two** primary fields (A and J) — that's a non-minimal ontology in Duda's reading.

**3. Duda's closing stance (2026-05-13 03:21):**

> *"Please let us know when you also reach charge quantization — before, your charged particles would explode, so I have no idea what could you claim to have for models-of-particles."*

So Duda explicitly **does not accept** Chern-Simons linking as a satisfactory answer to charge quantization. He wants:

- single deeper field (not two)
- A as derived connection (not primary)
- Gauss law on F counting topological charge of the deeper field (standard topology)

Werbos's framework violates all three preferences. The disagreement is not just notational ("J was confusing notation"); it's structural.

#### What this adds to the evaluation

| Insight from thread | Implication for our M6 eval |
| --- | --- |
| Werbos/DeepSeek themselves call Chern-Simons "topological" | confirms §2.1 — no topology claim was misleading |
| Duda's objection is single-field vs two-field | adds a structural critique we should address if pursuing M6 |
| Duda's Aharonov-Bohm point | M6 sandbox should test whether A in Ouroboros captures AB phase the same way Duda's connection does |
| Duda closes the thread unconvinced | the disagreement is unresolved as of 2026-05-13; don't expect Werbos to satisfy Duda anytime soon |
| Werbos uses DeepSeek to write technical replies | take Werbos's mathematical defenses with mild caution; the rigor is partly AI-mediated |

#### A note on the DeepSeek-authored portion

Worth flagging because of context: per memory `feedback_no_static_solitons` and the cardinal "no AI slop" rule Duda raised on 2026-05-08, **a meaningful portion of Werbos's technical defense in this thread is AI-generated by DeepSeek**. This doesn't invalidate the technical content — DeepSeek's argument that J has its own kinetic term and therefore propagates independent DoF is mathematically correct. But:

- The deeper objections Duda raises (Aharonov-Bohm, single-field ontology) **are not addressed** in the DeepSeek response. They get acknowledged but not answered.
- Werbos's own follow-ups are shorter and more discursive ("Our theorems do get very deep into that topology") rather than technical replies.

If we pursue M6, we should engage Werbos through the mathematics directly (his numerical experiments, his Lean proofs, his Sawada anchor) rather than through the rhetorical defenses. That's where the real content lives.

---

## 3. Side-by-side with Duda's LdGS

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
| Numerical depth | sandbox 8 experiments + production M5.0/M5.1 done | numerical existence of chaoiton families (Z 20030162) |
| Formal verification | none | Lean 4 (placeholder defs; full proofs in progress) |
| Empirical anchor | hopfions in lab (Liu 2026, Nature Phys) | Sawada 1989/2003 long-range nuclear anomaly |

---

## 4. Strengths Ouroboros has that Duda doesn't

| Strength | Detail |
| --- | --- |
| Direct empirical anomaly | Sawada `v(r) ~ -C/r⁶` is a measured nuclear anomaly the framework explains directly |
| Simpler field algebra | two coupled vector fields vs a matrix field; ops are vector-calculus-natural |
| Maxwell exactly recovered | drops out as the linear-limit kinetic term, no tilt-curvature gymnastics |
| Lean-4 verified theorems | charge quantization mechanically proved (placeholder defs but structure verified) |
| Dark matter candidate | neutral chaoiton configurations couple only gravitationally |
| Engineering anchor | Werbos sketches sensor architecture for J_μ-channel nuclear detection — directly analog to SABER's induction-from-thermal goal |

The Sawada anchor in particular is genuinely strong. It's a real measured nuclear anomaly (long-range component of the strong force, P-wave π-π phase shift deviation) that Standard Model + QCD has trouble explaining. If Ouroboros's J_μ field naturally produces this `v(r) ~ -C/r⁶` London-type tail with the right strength, that's an experimentally falsifiable prediction we can test numerically.

---

## 5. Strengths Duda has that Ouroboros doesn't

| Strength | Detail |
| --- | --- |
| Gravity unification | g-axis boost component → linearized GR (gravitoelectromagnetism); Ouroboros explicitly stops before gravity |
| Lepton mass hierarchy | three biaxial axis choices give e/μ/τ; Ouroboros has no clear mass-spectrum mechanism |
| QM derivation | Klein-Gordon emerges from biaxial twist (Fig 9); Ouroboros doesn't derive Schrödinger |
| Mainstream connection | LdG + Skyrme + Einstein's teleparallelism are well-known frameworks; easier to publish / collaborate |
| Visual & numerical depth | 51-page slides, ported Mathematica code, paper plus extensive figures; Ouroboros is 10 pages |
| Particle-defect correspondence | Liu et al 2026 lab anchor for hopfions / skyrmions in real media; Ouroboros lacks structure-existence anchor |
| Already implemented | M5.0 + M5.1 done in OpenWave with R²=0.978 Coulomb result; Ouroboros has no OpenWave port |

---

## 6. SABER compatibility

Both frameworks are **substrate-compatible** with the per-defect (A, ω) thermal hypothesis ([5b_thermal_energy.md](5b_thermal_energy.md)). Heat-as-defect-oscillation works in either substrate; the difference is which field carries the engineering coupling.

| SABER stage | In Duda | In Ouroboros |
| --- | --- | --- |
| Heat substrate | (A, ω) excess of matrix twist | excess oscillation of chaoiton |
| Internal modulation | drive matrix tilts via EM | drive chaoiton via A or J perturbations |
| Outward coupling | tilt-wave (EM wave) at engineering distance | J_μ-channel radiation at engineering distance |
| Extraction | conventional EM coil (Faraday's law) | layered nanostructure sensor (Werbos's sketch) |

Notable: **Werbos's sensor architecture is structurally analogous to SABER's induction goal.** He proposes:

- Top shield (EM + thermal isolation)
- Quantum-intelligent processing (TQuA)
- Isolation + transduction (J_μ → electronic signal)
- Sensing material: nanostructured high-mass-number nuclei (W, Pb, heavy isotopes) tuned to specific J_μ frequency bands
- Phonon-isolated substrate

The SABER difference: SABER aims to **produce useful electric power** by harnessing thermal energy modulation, not just **detect nuclear events**. But the field-theoretic mechanism (modulate a deeper field, transduce to electronic signal) is the same family.

This means **if Ouroboros's J_μ field is real, SABER could in principle use either the EM channel (M5 framework) or the J_μ channel (M6 framework)** for the extraction step. The J_μ channel is potentially even more interesting because it penetrates conventional EM shielding (Werbos's argument for nuclear detection applies equally to nuclear-scale energy extraction).

---

## 7. Verdict — is this a good M6 candidate?

**Yes, technically viable as M6 — with significant caveats:**

### What makes it viable

1. The math is well-defined and partly machine-verified (Lean 4)
2. The Sawada anomaly is a real empirical anchor with falsifiable predictions
3. The 2-vector-field structure is simpler than Duda's matrix substrate to implement in Taichi
4. Time-periodicity as Derrick escape is now triple-confirmed (Duda + Close + Werbos)
5. Werbos's claimed 62 chaoiton families with L/Q matching electron g-factor is a numerical achievement worth reproducing
6. SABER engineering path is parallel — J_μ channel is conceptually similar to EM tilts

### Why it can't replace M5

1. Explicitly **stops before gravity** — no GEM, no full unification
2. **No lepton hierarchy** — single chaoiton family doesn't explain e / μ / τ
3. **Schrödinger / KG not derived** — quantum mechanics emerges from twist in Duda; emergence in Ouroboros is unclear
4. **Weaker clock-propulsion derivation** — numerical existence vs analytical Hamiltonian term
5. Werbos's "no topology" framing is misleading — Chern-Simons IS topology

### Best framing

**M6 as a cross-validation / alternative method, not as a replacement.** Worth pursuing in parallel after M5.4–M5.7 land, as:

- An independent test of the heat-as-oscillation hypothesis (different substrate, same hypothesis)
- A simpler sandbox for prototyping the resonance-hunt protocol (vector fields ship faster than matrix fields)
- A second engineering channel for SABER (J_μ alongside EM tilts)
- A way to engage Werbos productively in the group

### What "M6 method" would actually mean in OpenWave

```text
M5 (Duda):  matrix field M = ODO^T → emergent EM, QM, gravity (3D + 4D)
            production engine, lepton families, Cornell quarks, GEM
            primary scientific bet

M6 (Werbos):  two coupled vector fields A, J on Minkowski
              chaoiton-stability test, Sawada nuclear anomaly target
              alternative substrate for the same thermal hypothesis
              secondary bet — runs after M5.7 lands
```

---

## 8. Open questions / outstanding context

### Papers we haven't read yet

| Reference | Topic | Importance |
| --- | --- | --- |
| Zenodo 20030162 | Chaoiton existence theorem (62 families, L/Q range) | high — numerical evidence base |
| Zenodo 20077680 | Full Ouroboros Lagrangian paper | high — theoretical depth |
| Zenodo 20179985 | "Do you still believe in quarks?" (Schwinger dyon alternative) | medium — QCD alternative framing |
| Zenodo 20218067 | "Ouroboros model passes another test" — electron properties to ~1% | high — empirical test result |

Recommended: download all four locally to `scientific_source/` before committing to M6 work. The 20218067 paper is freshest and most directly relevant to whether the framework actually produces electron properties at engineering accuracy.

### Technical questions worth asking Werbos

| # | Question | Why it matters |
| --- | --- | --- |
| Q1 | Is Chern-Simons not topology? | Werbos / DeepSeek already concede it IS topology in-thread |
| Q2 | How does ω = 2mc²/ℏ derive from the Lagrangian analytically? | numerical existence is not the same as derivation |
| Q3 | How does Ouroboros distinguish e from μ from τ? | mass spectrum mechanism is missing |
| Q4 | What is the form of f(J·J) that gives the electron-g-factor match? | parameters of the 62-family search |
| Q5 | Does the J_μ field connect to gravity at all, or is it strictly pre-gravity? | for SABER engineering, this matters |
| Q6 | What's the predicted strength of the Sawada `-C/r⁶` tail in physical units? | falsification test |
| Q7 | How does A in Ouroboros capture the Aharonov-Bohm phase? | Duda's point — A is more fundamental than F; in his framework A is a derived connection that gets AB for free |
| Q8 | Why two primary fields (A, J) instead of one deeper field with A as derived connection? | Duda's single-field ontology objection — Werbos doesn't address it in the thread |

### Where to engage in the thread

Werbos's recent emails (in particular the "Do you still believe in quarks?" challenge and the "electron properties to 1%" claim) push for engagement. The right tone is **technical curiosity + honest pushback on the framing claims**. Specifically:

- Acknowledge the empirical strength (Sawada, electron g-factor match)
- Probe the "without topology" claim — Chern-Simons IS topological
- Ask about gravity and lepton hierarchy
- Offer the OpenWave numerical platform as a venue for cross-substrate comparison (matrix vs two-vector)

---

## 9. If we decide to pursue — minimal M6 sandbox shape

This is a placeholder for if/when M6 becomes a concrete commitment. Don't start until M5.4-M5.7 are in motion.

### Minimal viable sandbox

1. **Two-vector-field substrate**: `A_μ` and `J_μ` on the existing Taichi grid (Vector(4) per voxel, two fields)
2. **Lorenz constraint enforcer**: projection step in the leapfrog to keep `∂^μ A_μ = ∂^μ J_μ = 0`
3. **Free-field test**: f(J · J) = 0; verify Maxwell behavior in A; verify J propagates as free massless vector
4. **Chern-Simons linking number kernel**: numerically compute `Q[A, J]` on configurations; verify integrality on seeded test fields
5. **First chaoiton hunt**: seed a localized A + J configuration with non-zero Q; turn on f(J · J) coupling; look for time-periodic localization

### Falsification gates

| Gate | Test | Decision |
| --- | --- | --- |
| Gate 1 — Maxwell limit | does A field behave as standard EM with f = 0? | should pass trivially |
| Gate 2 — charge quant numerics | can we measure integer Q on seeded fields? | required to proceed |
| Gate 3 — chaoiton existence | does any (A, J, f) combo produce a long-lived time-periodic localized state? | the load-bearing test |
| Gate 4 — Sawada `1/r⁶` tail | does the J_μ-mediated force between two chaoitons fit `-C/r⁶`? | empirical match test |
| Gate 5 — electron g-factor | does the L/Q ratio of the lowest-energy chaoiton match electron's g/2 ≈ 1.00115? | Werbos's headline claim |

Gates 1–3 are cheap (Maxwell is easy; integrality + chaoiton hunt are existing-machinery-style). Gates 4–5 are harder and would justify a real M6 phase.

### Sequencing relative to M5

```text
2026-Q2/Q3  M5.3 → M5.4 → M5.5    Duda matrix substrate, primary bet
2026-Q3/Q4  M5.6 → M5.7           biaxial KG + resonance hunt
                                  ↓
                                  M5.7 unlocks 5b thermal AND
                                  triggers M6 evaluation start
2027-Q1+    M6 sandbox             two-vector-field Ouroboros
            (parallel to M5.8)     run gates 1-3 as quick wins
```

M6 is genuinely a **post-M5.7** call. Don't pre-commit. Re-evaluate when:

- M5.7 metastable resonance result is in hand
- We've read all four outstanding Werbos zenodo papers
- We've had at least one direct technical reply from Werbos clarifying the issues above
