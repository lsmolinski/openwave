# M5 / Liquid-Crystal — Question Tracker + Hardest Pieces

**Purpose:** single source of truth for open research questions on M5, the resolution chronology, and the long-running **hardest-pieces** status board. Mirrors the M6 tracker pattern. Update this doc when a question opens, gets sent to the group (Duda / Close / Yee), gets answered, or gets demoted. Source material: [`3b_lagrangian_roadblocks.md`](3b_lagrangian_roadblocks.md) (the group-questions Q1–Q4), [`4a_convo_2026.05.12.md`](4a_convo_2026.05.12.md) (Duda's 2026-05-14/15 substrate + eigenvalue replies + his own open items), and `NEPTUNYA/SABER/1_model_selection.md` (Duda's stated limitations of his own LdGS model — physics content reproduced here; the SABER engineering-timeline analysis stays private per the cardinal rule).

**Sister docs:**

| Doc | Purpose |
| --- | --- |
| [`0b_M5_roadmap.md`](0b_M5_roadmap.md) | Phase sequence (M5.0 → M5.9) + current state |
| [`3b_lagrangian_roadblocks.md`](3b_lagrangian_roadblocks.md) | M5.2 negative-result diagnosis + original group questions |
| [`4a_convo_2026.05.12.md`](4a_convo_2026.05.12.md) | Duda thread study — substrate decision, eigenvalue map, slides |
| [`1a_lagrangian_framework.md`](1a_lagrangian_framework.md) | Framework + the full email thread history |

**Last updated:** 2026-05-26 (**M5.5 CORE COMPLETE**). M5.4 matrix-field migration is done + merged (full 3×3 `M = ODO^T`; matrix triple buffer + `eigen_decompose` lynchpin + matrix seeders + `‖M−D‖_F`/`‖Ṁ‖_F` trackers; **M5.1 Coulomb reproduced on M**, R²=0.9704 relaxed + R²=0.9959 analytical page-18, both attractive). **M5.5 (paper Lagrangian + V(M)) core complete (2026-05-26)**: Duda's Eq.18 action ported to production (`compute_curvature_flux` + `evolve_M` + matrix Hamiltonian) — "Evolve PDE" is LIVE and **energy-conserving** (symplectic; GUI "bounded-not-bug" verified — a hedgehog *sloshes* under bounded, H-conserving nonlinear curvature dynamics, NOT a blow-up); `V(M)` Eq.13 LdG implemented (off by default); Faber `F²`=Skyrme-family-kinetic + Derrick-stability validated. **M5.5.5 (EM-from-tilts + Faber EM Lagrangian) folds into M5.6.** **Next phase = M5.6** (biaxial twist + KG emergence); the *questions* that gate it are Q7 (exact V(M) coeffs) + Q8 (Faber regularization exact form) — both Duda-open, neither blocking (Eq.13 baseline + Faber port proceed). **Added earlier this arc:** hardest-pieces board + Duda's stated limitations + M6 → M5 cross-pollination. **Q12** (Bell / Kochen-Specker vs M5's definite-orientation defect) — foundational note in `1b § Foundational stance`.

---

## Active count

```text
0 IMMEDIATE   M5.5 CORE COMPLETE (2026-05-26) — Eq.18 action live in
              production ("Evolve PDE"), energy-conserving (bounded-not-bug
              verified); V(M) Eq.13 implemented (off by default). M5.4
              substrate + Coulomb gate done + merged. Next phase = M5.6
              (biaxial twist + KG); Q7/Q8 gate it but neither blocks.

2 NEAR-TERM   Q7  exact V(M) form — now gates M5.6 (Duda's own open research
(gate M5.6)       question; Eq.13 LdG baseline implemented + off by default)
              Q8  Faber regularization exact form — gates M5.6
                  (mechanism validated: F²=Skyrme kinetic, Derrick-stable;
                  exact running-coupling scheme still open)

5 BACKGROUND  Q4  Liu et al. lab anchor — does it change sim priority?
(long-tail)   Q9  deeper substrate beneath M (anisotropic fluid?)
              Q10 weak-force clean SU(2) mechanism (gap; beta decay as
                  topology reconnection is a partial answer)
              Q11 Close's Eq.23 exact implementation form — gates M5.7
              Q12 Bell / Kochen-Specker vs M5's definite-orientation
                  defect (foundational; defense template + MERW bridge)

5 RESOLVED    Q1  initial-condition construction (no static soliton —
              triple-confirmed: Duda Fig.10 + Close + Werbos + M5.2)
              Q2  V(ψ) shape (wrong substrate — matrix M required)
              Q3  connection/curvature A-layer (load-bearing; paper §II)
              Q5  substrate: full M = ODO^T vs Q-tensor (full M)
              Q6  eigenvalue → physics mapping (1=EM, δ=QM, g=gravity)

Total: 12 questions (0 immediate, 7 open, 5 resolved). M5.5 core complete; M5.6 next.
```

---

## OPEN QUESTIONS

| ID | Question | Surfaced | Gates | Status |
| --- | --- | --- | --- | --- |
| Q7 | Exact V(M) form — Duda's Eq.13 Landau-de Gennes `V_LG = a Tr(M²) − b Tr(M³) + c (Tr(M²))²`, or "slightly different"? Duda calls this "the most difficult" part and an open research question. | Duda 2026-05-15 (`4a §12`) | M5.6 | OPEN. Eq.13 LdG + Eq.12 eigenvalue-preference **both implemented in production** (`V_M`/`dV_M`, off by default at `ldg_a/b/c=0`); exact `Λ=(1,δ,0)` coeffs still Duda-open. Not blocking — M5.6 turns it on with the biaxial vacuum. |
| Q8 | Faber regularization — exact form to "activate" V(M) and produce the running-coupling effect. Duda points to Faber's two papers as the starting scheme. | Duda 2026-04-19 / 05-15 | M5.6 | OPEN. **Mechanism validated (M5.5.3)**: Eq.18 `F²` is the Skyrme-family 4th-order kinetic → Derrick-stable (`E(L)=E_F/L+E_V·L³`, finite core+mass), no separate Skyrme term. Faber's exact running-coupling scheme (`reference_faber_regularization`) still open — port + adapt in M5.6. |
| Q11 | Close's Eq.23 exact implementation form — direct form vs vector-potential `s = ∇×A` vs divergence-cleaning projection; and should the resonance amplitude sweep span the full `m ∈ {−1, 0, +1}` dipole family? | Rodrigo → Close 2026-04-19 (`1a`) | M5.7 | OPEN, pending Close. Resonance-hunt protocol can proceed on best-current reading (`l=1`, `A/λ ∈ {0.5, 1, 2}`). |
| Q4 | Liu et al. *Nature Physics* 2026 (first direct laser creation of isolated hopfions + skyrmions) — does the lab confirmation change what the framework should simulate first (hedgehog before hopfion? simpler stabilizer before LdG?)? | `3b` (2026-05-12) | — | OPEN. Lab-existence anchor for topology-as-particles. Current plan keeps hedgehog-first; hopfion is post-M5.8 frontier (`project_particle_defect_correspondence`). |
| Q10 | Weak force — is there a clean matrix mechanism (SU(2) chiral) the way EM/gravity/strong have one? | `4a §7` | — | OPEN GAP. Partial answer in slides: beta decay appears as a topology-reconnection event (defect-class transition), not a force in the EM/gravity sense. No clean SU(2) mechanism yet. |
| Q9 | Is there a deeper substrate beneath the matrix field M (an "anisotropic fluid")? Duda hints the matrix may be effective, with something deeper below. | Duda 2026-05-15 (`4a §1, §12`) | — | OPEN / PARKED. Matches OpenWave's existing granule-level picture (matrix effective, granules deeper). Not actionable now. |
| Q12 | How does M5's definite-orientation defect (real spin axis + clock phase, deterministic local PDE) reconcile with Bell / Kochen-Specker contextuality? | Duda↔Hadley thread 2026-05-11; Close/Adenier thread 2026-05 | — | OPEN (foundational). Defense template: field evolution is local-realistic (cannot violate Bell); QM statistics enter at the measurement/collapse layer (Duda's Malus analogy + Feynman-ensemble), bridged by MERW `ρ=\|ψ\|²` (`4a §11b.4`). Close building the parallel Adenier / Eq.10-factoring defense. Full note in [`1b_topological_defect.md` § Foundational stance](1b_topological_defect.md#foundational-stance--m5-as-a-local-realistic-field-bell--kochen-specker). |

---

## RESOLVED QUESTIONS

| ID | Question | Resolution |
| --- | --- | --- |
| Q1 | Initial-condition construction for a stable single particle — relax against full Lagrangian (a), exact-soliton ODE BVP (b), or something else (c)? | RESOLVED — premise was wrong: **there is no static stable solution**. The particle is a **time-periodic resonance**. Triple-confirmed: Duda paper Fig.10 (4D Lorentz negative-energy auto-propels the clock) + Close email ("l=1 amplitudes, I doubt you'll find completely stable non-radiating solutions") + Werbos chaoiton claim + our own M5.2 empirics. See `feedback_no_static_solitons`. |
| Q2 | Is there a V(ψ) admitting a stable 3D hedgehog soliton (we tried φ⁴, biharmonic)? | RESOLVED — wrong substrate. The real potential is Eq.13 LdG `V_LG(M)` on the **matrix field M = ODO^T**, not any V(ψ) on Vector(3). None of our scalar V(ψ) probes were testing the actual proposed potential. Folds into Q5 + Q7. |
| Q3 | Is an explicit connection/curvature `A`-as-primary-field layer load-bearing for hedgehog stability, or the same physics in different clothing? | RESOLVED — load-bearing (paper §II). `A_μ = [M, ∂_μ M]` (Eq.19), `F_μν = ∂_μ A_ν − ∂_ν A_μ` (Eq.20). The Brouwer-degree winding integral we use is the right diagnostic, but the underlying field is the matrix M. |
| Q5 | Substrate: full real-symmetric matrix `M = ODO^T` (6 DoF) vs traceless Q-tensor (5 DoF)? | RESOLVED 2026-05-15 — **full M, no Q-tensor pivot**. Confirmed three ways: Duda's image 2 writes `M = ODO^T` with `D = diag(g,1,δ,0)`; the slides reaffirm it everywhere; Duda's follow-up used the same notation. Refactor green-lit (`4a §2`). |
| Q6 | Eigenvalue → physics mapping — what does each eigenvalue of D tag? | RESOLVED 2026-05-15 (Duda direct quotes, `4a §8`) — `1` = EM tilts (highest Lagrangian contribution); `δ ~ ℏ` = QM twists (quantum phase); `g ≫ 1` = gravity boosts (hedgehog = black hole); `0` (4D only) = null/time axis (clock propulsion). The 3 leptons come from 3D spatial axis-choice, NOT from the g eigenvalue. |

---

## HARDEST PIECES — status board

The long-running hardest-pieces tracker for M5 (mirrors M6's). These are the load-bearing unknowns the framework must resolve, distinct from the discrete Q-numbered questions above — a piece can persist across many phases.

| Hardest piece | Gates | Status (2026-05-26) | M6 lend? |
| --- | --- | --- | --- |
| Matrix substrate migration (Vector(3) ψ → `M = ODO^T`) | M5.4 | ✅ **COMPLETE (2026-05-26).** `ti.Matrix.field(3,3)` triple buffer + `eigen_decompose` lynchpin (director recovery 0.9995) + matrix seeders + `‖M−D‖_F`/`‖Ṁ‖_F` trackers. **M5.1 Coulomb reproduced on M** — R²=0.9704 relaxed + R²=0.9959 analytical page-18 cross-val, both attractive. Operators verified vs analytic. Rendering re-sourced from M (on-screen verified); live wave path retired (ψ engine dormant legacy). | — |
| **V(M) potential form** (Duda's #1 limitation, Q7) | M5.6 | 🔶 Eq.13 LdG + Eq.12 eigenvalue-preference **both implemented in production** (`V_M`/`dV_M`, off by default); exact `Λ=(1,δ,0)` coeffs = Q7 (Duda open). M5.6 turns it on. | ✅ strong (see below) |
| Regularization (Faber) to activate V(M) + running coupling (Duda #2, Q8) | M5.6 | 🔶 **Mechanism validated (M5.5.3)**: Eq.18 `F²` IS the Skyrme-family kinetic → Derrick-stable (finite core+mass), no separate term. Faber's exact running-coupling scheme still open — port + adapt in M5.6. | 🔶 partial |
| First metastable resonance (no static soliton) | M5.7 | Design locked (time-periodic per Close); awaiting M5.6. **M5.5.4 shipped the live energy-conserving Eq.18 evolution** the hunt runs on. | ✅ strong |
| KG-from-twist emergence | M5.6 | Closed-form available (slide p.32 + `5a §5`); **evolution machinery validated (M5.5.2 twist evolution, energy drift 0.21%)**. Biaxial-hedgehog port (incl. the z-axis disclination handling) not yet done. | — |
| Zitterbewegung clock / 4D negative-energy propulsion | M5.8 | Mechanism known (Fig.10); toy-model numerics available (slide p.33). **Externally validated 2026-05** as THE hard problem: group consensus that ALL macroscopic analogues fail at intrinsic non-damping oscillation (Couder droplets need a shaker; spinning needles damp like a pendulum) — only the elementary-particle clock is intrinsic + non-damping, so propulsion must be intrinsic. Validates the 4D negative-energy direction. See [9c § clock propulsion](9c_time_dynamics.md#the-de-broglie-clock-is-intrinsic--group-consensus-on-clock-propulsion-2026-05). | 🔶 partial |
| Biaxial eigenvalue hierarchy → lepton masses (Duda #3/#9/#10) | M5.9 | Calibration parameters; post-SABER-trigger. | — |
| Weak-force clean SU(2) mechanism (Duda #7, Q10) | — | GAP. Beta decay as topology reconnection is partial. Out of SABER scope. | — |
| Deeper substrate beneath M (Duda #8, Q9) | — | PARKED. Matches granule-level picture. | — |

### Duda's stated limitations of his own LdGS / M5 model

Physics limitations Duda (or his papers, as quoted in our notes) flags about his own framework. Sourced from the M5 research docs; reproduced from `NEPTUNYA/SABER/1_model_selection.md` (physics content — the SABER engineering-timeline analysis stays in that private doc).

| # | Limitation | Source |
| --- | --- | --- |
| 1 | **V(M) potential form unspecified.** *"the Higgs-like potential V(M) — its specific form is not yet determined"* | `1a:308` |
| 2 | **Regularization is the hardest part.** *"LdG regularization is where lepton masses come from, and it is the hardest part to include. The specific regularization form is still an open research question"* | `1a:230, 260` |
| 3 | **Biaxial axis-length parameters need numerical calibration** — *"their exact values require numerical simulation … treat them as calibration parameters"* | `1a:262` |
| 4 | **Charge magnitude + rigorous Coulomb derivation open** | `1a:115` |
| 5 | **Lepton mass ratios deferred, mechanism only** | `_overview:393` |
| 6 | **Static hedgehog doesn't survive dynamic evolution** — Q decays +0.996 → ~0 within 4–5 steps regardless of V(ψ); confirms the triple admission that no static stable solitons exist | `3b:2-9, 103` |
| 7 | **Weak force mechanism unaddressed** — *"topology-changing events … still an open research question"* | `4a:337-354` |
| 8 | **Deeper substrate beneath M hinted but unconfirmed** — *"there might be an even deeper 'anisotropic fluid' beneath"* | `4a:602` |
| 9 | **Extreme biaxial eigenvalue hierarchy unexplained** — ~3477:1 axis-length ratio not derived from first principles | `1a:230-231` |
| 10 | **Three lepton families pinned to "3 axes in 3D" by symmetry, not derived** — mass ratios still need V(ψ) curvature input | `1b:469` |

### Which M5 phase each limitation gates (physics roadmap)

SABER engineering-timeline impact is analyzed in the private `1_model_selection.md` (§ "How the limitations affect SABER timeline"); here we keep only the OpenWave-physics phase-gating. The chain `M5.5 → M5.6 → M5.7 → 9b` is the documented bridge to the SABER thermal counterpart (see [`9b_thermal_energy.md`](9b_thermal_energy.md)).

| Limitation | Gates | Reading |
| --- | --- | --- |
| #1 V(M) form | M5.5 | Direct blocker for the paper Lagrangian. Until V(M) is pinned, defect dynamics under perturbation aren't uniquely simulable. |
| #2 regularization | M5.6 | The hard step. Faber's scheme is the working hypothesis; if it ports cleanly, M5.6 unblocks. |
| #6 no static solitons | (design) | Already a design choice — pursue Close's `l=1` resonance protocol (M5.7), not static seeds. |
| #7 weak force | — | Out of all near-term scope. |
| #9 biaxial hierarchy | M5.9 | Lepton masses are post-SABER-trigger; a scientific deliverable, not an engineering input. |

---

## M6 → M5 cross-pollination — what the closed M6 sandbox can lend M5

M6 is closed (DM paper published, sandbox-only, no Taichi port), but its sandbox **methodology** transfers to M5. M6's scipy sandbox is far cheaper to prototype in than M5's Taichi matrix substrate, so several M5 unknowns can be de-risked in a scipy M5-sandbox *before* committing Taichi kernels — exactly the path M6 took. **This is sandbox methodology transfer, NOT an M6 production port.**

### The V(M) lead — M6 settled its potential the hard way, and the method transfers

M6's biggest relevance to M5's #1 hardest piece (V(M) form): M6 **settled its neutral-sector potential + geometry** as a concrete BVP. The neutral chaoiton ground state is the 3D spherical `l=1` p-wave ODE

```text
β'' + (2/r)β' − (2/r²)β − m_J²β + 4gβ³ = 0
```

found in **sandbox v9 phase 2** by treating `m_J` as a free eigenvalue with a proper `l=1` origin BC (`β ~ B0·r`) and `B0` fixed — *after* forward-IVP/shooting had failed. The transferable lesson is methodological, not the equation itself:

| M6 lesson (v9 phase 2 + v10) | M5 application |
| --- | --- |
| Forward-IVP/shooting does NOT find localized soliton ground states; **BVP with the mass parameter as a free eigenvalue does** | M5.6/M5.7 — solve the matrix-field defect radial profile as a BVP-eigenvalue problem, not by guessing V(M) and forward-stepping (note: M5.5.4's *forward* leapfrog confirms this — an unregularized hedgehog sloshes, energy-conserving, but does not settle into a localized profile without the M5.6 regularized BVP) |
| Correct `l=1` p-wave origin BC (`β ~ B0·r`) is what unlocked the ground state | M5's hedgehog is an `l=1`-class structure — the BC handling transfers directly |
| `η` geometry-conversion factor `(∫β²r dr)/(∫β²r² dr)` reconciles cylindrical vs spherical energy integrals | M5 will compare matrix-field energies across defect geometries — same conversion needed |
| Pohozaev-type virial identity as a self-consistency diagnostic (M6's `m_J/η = 1.21` family-invariant, Q47) | A scaling identity is a cheap correctness check on any M5 soliton profile |
| Sandbox-first: prototype in scipy, validate, *then* port to production | Prototype the M5 V(M) + matrix-defect-profile BVP in a scipy M5-sandbox before the Taichi M5.6+ kernels (the path M5.5 followed: `sandbox_v4` scipy/numpy → Taichi port) |

Canonical M6 numerical recipe (charged sector via 2-function IVP, neutral sector via the BVP above) is consolidated in `m6_ouroboros/research/0d_canonical.md` — the reference if we want to mirror the methodology.

### M6 strengths available to M5 regardless of M6's gate outcomes

Imported from `m6_ouroboros/research/0b_model_gates.md` § "What M6 offers regardless of gate outcomes" and `1_model_selection.md`:

| M6 strength | Relevance to M5 |
| --- | --- |
| Neutral chaoiton **ground state reachable** via BVP (v9 phase 2) | Proof-of-method that a topological soliton's ground state IS findable numerically — encouragement + recipe for M5.7's metastable resonance hunt |
| DM mass prediction delivered (`m_χ = 0.46 MeV`, `m_J = 0.62 MeV`, `C = 770 MeV·fm`) | Independent topology-as-particles result; cross-validates the shared "no static soliton / time-periodic" premise M5 also relies on |
| Lepton spectrum at PDG precision via 2-function reduction | Demonstrates the ansatz-reduction discipline (collapse DoF before solving) that M5.5 can borrow for the matrix field |
| Charge quantization via Chern-Simons / Hopf (Lean-proved) | Different topological flavor (Hopf-linking) from M5's Brouwer-winding — but confirms charge-from-topology is robust across frameworks |
| Methodology cross-pollination (η, `l=1` BC, virial diagnostic) | Listed above — the concrete transferable pieces |

**Caveat (load-bearing):** M6 is **not** the thermal substrate and does not supply M5's `(A, ω)` thermal degree of freedom — that lives on M5's matrix field. M6's contribution to M5 is **method + credibility cross-validation**, not physics substrate. A hypothetical M6-substrate thermal cross-check would require a full Taichi production port (ensemble + external-drive kernels), which is NO-GO / deferred indefinitely.

---

## Notes on scope

- This tracker covers **M5 physics/framework questions for the group** (Duda / Close / Yee) plus the hardest-pieces board. M5.4 implementation decisions (Taichi storage layout, eigen-kernel design, granule/glyph UX calls) are tracked as **roadmap tasks** in [`0b_M5_roadmap.md`](0b_M5_roadmap.md) and [`4b_rendering_features.md`](4b_rendering_features.md), not here — same split as M6 (questions here, sandbox tasks in the work log).
- M5 has **no active email round in flight** (unlike M6's intense Werbos/DeepSeek cadence). Duda's last substantive reply was 2026-05-15. **M5.5 (core complete 2026-05-26) now gives concrete results to seed a Q7/Q8 outreach**: the Eq.18 action runs in production + is energy-conserving, and the **V(M)-is-rotation-invariant finding** (V acts only on the eigenvalue/regularization sector, not the twist) sharpens the Q7 ask — Duda's exact `Λ=(1,δ,0)` Eq.13 coeffs + Faber's exact running-coupling scheme (Q8) become the natural questions once M5.6 needs them turned on. None blocks current work.
