# M5 / Liquid-Crystal — Question Tracker + Hardest Pieces

**Purpose:** single source of truth for **(a) what physical reality the M5 substrate reproduces and how (MODEL EMERGENCE)**, (b) open research questions on M5, (c) resolution chronology, and (d) the long-running **hardest-pieces** status board. Mirrors the M6 tracker pattern. Update this doc when a question opens, gets sent to the group (Duda / Close / Yee), gets answered, demoted, or when new emergent physics is verified. Source material consolidated here from [`1b_topological_defect.md`](1b_topological_defect.md) § STRATEGIC MAPPING + § OUTGOING-WAVE L+T + § Empirical validation, and [`4a_convo_2026.05.12.md`](4a_convo_2026.05.12.md) § 7 Force unification + § 8 Eigenvalue mapping — so this tracker reads stand-alone. The deeper *explanations* (what is a defect, what is a twist, derivation walk-throughs) stay in 1b / 4a / 5a; this doc carries only the emergent-reality catalog + numerical mechanism + status.

**Sister docs:**

| Doc | Purpose |
| --- | --- |
| [`0b_M5_roadmap.md`](0b_M5_roadmap.md) | Phase sequence (M5.0 → M5.9) + current state |
| [`3b_lagrangian_roadblocks.md`](3b_lagrangian_roadblocks.md) | M5.2 negative-result diagnosis + original group questions |
| [`4a_convo_2026.05.12.md`](4a_convo_2026.05.12.md) | Duda thread study — substrate decision, eigenvalue map, slides |
| [`1a_lagrangian_framework.md`](1a_lagrangian_framework.md) | Framework + the full email thread history |

**Last updated:** 2026-06-05 (**§ Empirical validation table added** as the document header — every headline M5 claim now has a runnable reproduction, re-verified on the 4×4 engine after the M5.8.1 promotion + ψ-retire). Previous: 2026-05-28 (**M5.7 COMPLETE** — free disperses→4D/M5.8, driven sustains→9b). M5.4 matrix-field migration done + merged; M5.5 (Eq.18 action + V(M)) core complete (2026-05-26). **M5.6 (biaxial twist + KG emergence + Maxwell) complete (2026-05-27)**: KG mass is **GEOMETRIC** (minimal coupling to the hedgehog connection `Â`, not an added `V_ψ` — Fig.9 reproduced); biaxial hedgehog dynamically sources its own twist + restoring mass (M5.8 clock seed); **Faber regularization ported** → mass pinned `E∝1/r₀` (the M5.9 calibration handle); **Maxwell recovered both routes** (hydro↔EM dictionary + Faber `R=Γ×Γ` closed field strength); production port done — biaxial seeder + analytic Cardano eigensolver (replaced the `ti.sym_eig`-on-biaxial bug) + V-on amplitude confinement + EM/thermal/clock viz (`∇·n̂`/`‖∇×n̂‖` meshes + E/B glyphs, glyph↔mesh alignment fixed). **M5.7 (resonance hunt) complete (2026-05-28):** the free vs driven split — a *free* defect disperses (M5.7.1 seeded + M5.7.2 intrinsic, both resolution-confirmed nulls ⇒ the metastable particle/clock is 4D, M5.8) but a *driven* defect sustains a bounded, frequency-selective `(A,ω)` excess (M5.7.3 9b.1 preview, EM-lever ~3× free at the resonant natural mode ⇒ the driven-thermal substrate for 9b). Carry-overs (5e two-defect, 5f dipole-viz, gauge-charge, ψ-cleanup, Eq.23 cross-check) migrated to M5.8. **Next phase = Rodrigo's call (M5.8 4D Zitterbewegung clock, or full 9b.1 thermal).** **Q7 (exact V(M) coeffs) + Q8 (Faber exact running-coupling) are no longer gates** — M5.6 shipped the Eq.13 `b=0` confinement interim + the validated Faber port; their exact Duda-open forms now feed M5.9 lepton calibration, non-blocking. **Q12** (Bell / Kochen-Specker vs M5's definite-orientation defect) — foundational note in `1b § Foundational stance`.

---

## Empirical validation — key findings → runnable reproductions

The model's verified emergent physics, each with a runnable reproduction script — on the production 4×4 engine where the finding is engine-level (charge, Coulomb, Eq.18, LdG), standalone numpy/sympy where the finding is analytic (Maxwell routes, KG-geometric, Faber mass), and Duda's own 1+1D toy for the clock. **All re-verified 2026-06-05** after the M5.8.1 4×4 promotion + ψ-retire; the engine-level checks reproduced their original 3×3 numbers, empirically confirming the promotion is physics-preserving:

| Key finding | Result | Runnable reproduction |
| --- | --- | --- |
| Topological charge is INTEGER + ADDITIVE (Q = director winding) | ✅ Q=±0.996 single hedgehogs; pair: ±1 around each core, 0.000 on the enclosing sphere | `m5_8_1_topo_charge.py` |
| Coulomb: attractive 1/d force between opposite defects | ✅ R²=0.9781, b<0 — matches the M5.1 reference (0.978) digit-for-digit | `m5_4_coulomb_matrix.py` |
| Coulomb — Duda's page-18 Mathematica cross-validation | ✅ R²=0.9959, attractive | `m5_4_coulomb_page18.py` |
| E = ∇·n̂ / B = ∇×n̂ → Maxwell (hydro↔EM dictionary) | ✅ PASS — ∇·B=0, Gauss charge identity, Faraday, Lorentz force | `m5_6_4a_hydro_em.py` |
| EM-from-tilts: Faber R=Γ×Γ → abelian Coulomb far-field | ✅ PASS — closed field strength, ‖R‖~1/r², running-coupling onset at r₀ | `m5_6_4b_faber_curvature_em.py` |
| KG mass is GEOMETRIC (biaxial twist, not an added V) | ✅ PASS — explicit mass cancels exactly; mass from the hedgehog connection Â | `m5_6_1_kg_operator_check.py` |
| The biaxial hedgehog SOURCES its own twist + restoring mass (`C_μν≠0` — the M5.8 clock seed) | ✅ ‖C‖~1/r² mass source confirmed; twist sourced by the background curvature + restoring force — stable, energy-conserving | `m5_6_2a_biaxial_hedgehog.py` + `m5_6_2b_biaxial_evolution.py` |
| Finite defect mass pinned by the core radius, E∝1/r₀ (Faber — the M5.9 lepton calibration handle) | ✅ E·r₀ constant (CV=0.0%); E₀=α_f·ℏc·π/(4r₀); r₀=2.2132 fm → 0.511 MeV e⁻; reproduced on Duda's matrix substrate | `m5_6_3a_faber_regularization.py` + `m5_6_3b_faber_on_M.py` |
| Eq.18 matrix evolution conserves energy | ✅ secular drift 2.15%→1.13%→0.03% as dt→0 (symplectic) | `m5_5_4_matrix_evolution_check.py` |
| LdG well confinement (V(M) Eq.13, b=0 amplitude well) | ✅ confines ~3.3× across k∈[0.5,25], no blow-up | `m5_6_5c_prod_scale.py` |
| Duda's 1+1D time-crystal kink (arXiv:2501.04036) — the clock | ✅ ω* = √(70/61) = 1.0712 reproduced EXACTLY; energy drift ~10⁻¹⁵; clock robust under +5% perturbation | `m5_8_0b_collective_clock.py` |
| Clock propulsion is SPONTANEOUS + bounded ("propelled by mass") | ✅ ω=0 is an unstable maximum → relaxation spins the clock up to ω*; every seed → ω* (no noise needed); E floor 2.1257 > 0 | `m5_8_0c0d_propulsion_robustness.py` |
| M5.7 resonance hunt: free defect disperses (expected — Derrick), DRIVEN defect sustains | ✅ free → dispersal (the empirical push to the 4D clock, M5.8); driven → bounded, frequency-selective (A,ω) excess, ~3× free at resonance | `m5_7_*.py` |

No particle stability yet — that is M5.8's job (the particle is a 4D time-crystal, not a static 3D lump; Derrick + Duda + the M5.7 null all agree). Everything above comes straight out of the medium topology.

---

## Active Count

```text
0 IMMEDIATE   M5.7 COMPLETE (2026-05-28). Free vs driven split: free
              defect DISPERSES (M5.7.1 seeded + M5.7.2 intrinsic, both NULL
              resolution-confirmed) ⇒ particle/clock is 4D (M5.8); DRIVEN
              defect SUSTAINS a bounded frequency-selective (A,ω) excess
              (M5.7.3 9b.1 preview, EM-lever ~3× free at resonance) ⇒ the
              driven-thermal substrate for 9b. Root cause of free dispersal:
              V confines amplitude not orientation (M5.6.5c). Carry-overs
              migrated to M5.8. NEXT PHASE = Rodrigo's call (M5.8 4D clock,
              or full 9b.1 thermal).

1 SELF-       Q11 Close's Eq.23 exact form — NO LONGER a hard gate (PLAN B
DETERMINE         = self-determination, likely Plan A): Eq.23 is in his
(was gate)        published paper (read + test all 3 forms ourselves, the
                  physical one preserves ∇·s=0); seeded-vs-intrinsic is
                  empirical (the sim answers it). Don't bottleneck on email.

6 BACKGROUND  Q4  Liu et al. lab anchor — does it change sim priority?
(long-tail)   Q7  exact V(M) coeffs — NO LONGER gates (M5.6 shipped Eq.13
                  b=0 confinement interim); exact Λ=(1,δ,0) form feeds M5.9
              Q8  Faber exact running-coupling — NO LONGER gates (M5.6
                  ported + validated the mechanism, mass pinned E∝1/r₀);
                  exact scheme refines M5.9 calibration
              Q9  deeper substrate beneath M (anisotropic fluid?)
              Q10 weak-force clean SU(2) mechanism (gap; beta decay as
                  topology reconnection is a partial answer)
              Q12 Bell / Kochen-Specker vs M5's definite-orientation
                  defect (foundational; defense template + MERW bridge)

5 RESOLVED    Q1  initial-condition construction (no static soliton —
              triple-confirmed: Duda Fig.10 + Close + Werbos + M5.2)
              Q2  V(ψ) shape (wrong substrate — matrix M required)
              Q3  connection/curvature A-layer (load-bearing; paper §II)
              Q5  substrate: full M = ODO^T vs Q-tensor (full M)
              Q6  eigenvalue → physics mapping (1=EM, δ=QM, g=gravity)

Total: 12 questions (0 immediate, 7 open, 5 resolved). M5.7 complete; M5.8 (4D clock) / 9b.1 (thermal) next — Rodrigo's call.
```

---

## OPEN QUESTIONS

| ID | Question | Surfaced | Gates | Status |
| --- | --- | --- | --- | --- |
| Q7 | Exact V(M) form — Duda's Eq.13 Landau-de Gennes `V_LG = a Tr(M²) − b Tr(M³) + c (Tr(M²))²`, or "slightly different"? Duda calls this "the most difficult" part and an open research question. | Duda 2026-05-15 (`4a §12`) | ~~M5.6~~ → M5.9 | OPEN — **no longer gating**. M5.6 (✅) shipped the working interim: the 3-term Eq.13 form has **no biaxial minimum** (`∂V/∂λ` collapses nonzero eigenvalues to one root), so a `b≠0` term erodes δ → uniaxial; M5.6.5c uses a **`b=0` amplitude well** `V=a·Tr(M²)+c·(Tr(M²))²` (min at `Tr(M²)=1+δ²`) that confines amplitude while leaving δ exactly flat (live in production, `LDG_STIFFNESS_K`). A fully biaxial-STABLE vacuum needs an extra invariant — that exact `Λ=(1,δ,0)` form is the Duda-open piece, feeds M5.9 lepton calibration. |
| Q8 | Faber regularization — exact form to "activate" V(M) and produce the running-coupling effect. Duda points to Faber's two papers as the starting scheme. | Duda 2026-04-19 / 05-15 | ~~M5.6~~ → M5.9 | OPEN — **no longer gating**. M5.6 (✅) **ported Faber's MTF** (`Universe` 11/2025/113): reproduced `E₀=α_f ℏc·π/(4r₀)`, mapped `Λ=q₀⁶/r₀⁴` onto `M=ODO^T` via spatially-melting eigenvalues, **mass pinned `E∝1/r₀`** (`r₀=2.2132 fm → 0.511 MeV`), and confirmed Faber `R=Γ×Γ` → homogeneous Maxwell + running-coupling onset at `r₀` (M5.6.4). The mass knob is now the physical `r₀` tied to `α_f`. The exact `arctan` profile re-derivation (vs imposed) is the remaining Duda-open piece → M5.9 BVP. |
| Q11 | Close's Eq.23 exact implementation form — direct form vs vector-potential `s = ∇×A` vs divergence-cleaning projection; and should the resonance amplitude sweep span the full `m ∈ {−1, 0, +1}` dipole family? | Rodrigo → Close 2026-04-19 (`1a`) | M5.7 | OPEN, pending Close. **PLAN B = self-determination (and likely Plan A — see below).** M5.7.1 already proceeded on best-current reading. We do NOT need to wait on Close: (1) Eq.23 is in his **published** paper (`research/theory/Equation-of-Everything.pdf`) — we can read it directly and **empirically test all three candidate forms ourselves** (the physical one preserves `∇·s = 0` + gives bounded, energy-conserving dynamics — a numerical discriminator, not an opinion); (2) the deeper "seeded resonance vs the defect's intrinsic oscillation" question is **empirical** — the simulation answers it (M5.7.1: seeded l=1 disperses; M5.7.2: does the intrinsic oscillation cohere?), Close's view is confirmation not gate. Close's framework is cited (published); the answers are ours to derive. M5 must not bottleneck on advisor email latency. |
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

## MODEL EMERGENCE — observed reality from the matrix field `M = ODO^T`

What known physical phenomena the M5 substrate reproduces from the Duda Eq.18 action, with the numerical mechanism + status + what's pending. **No model comparisons** (those live in the `1_model_selection.md`). **No conceptual explanations** of substrate / defect / twist / eigenvalues — those stay in `1b` / `4a` / `5a`. Just emergence + status.

**Substrate (recap, see Q6 RESOLVED).** Symmetric 3×3 `M(x) = O(x)·D·O^T(x)` on a 3D lattice; 5 DoF/voxel. Preferred shape `D = diag(g, 1, δ, 0)`. Eigenvalue → physics map (M5.6.4-verified):

| Eigenvalue | Physical mechanism | Phase |
| --- | --- | --- |
| `g ≫ 1` | gravity (boost; hedgehog ↔ black hole) | M5.8 (4D) |
| `1` | EM (tilts → Maxwell) | ✅ M5.6 |
| `δ ~ ℏ` | QM (twist → Klein-Gordon, mass is GEOMETRIC — Fig.9) | ✅ M5.6 |
| `0` (4D only) | time signature → Zitterbewegung clock propulsion (Fig.10) | M5.8 |

### PARTICLES (full spectrum)

| Phenomenon | M5 numerical mechanism | Status |
| --- | --- | --- |
| Quantized electric charge | Brouwer winding (S² → S² topology) on M's principal director; integer, conserved by topology | ✅ M5.1 (`compute_winding_number`) |
| Electron rest mass | Faber regularization `E₀ = α_f · ℏc · π / (4 r₀)`; physical knob `r₀ = 2.2132 fm → 0.511 MeV` | ✅ M5.6.3 (mass pinned `E ∝ 1/r₀`) |
| Lepton families `e / μ / τ` | 3 axis-choices of biaxial `Λ = (g, 1, δ)`; same Q, different masses | 🚧 M5.9 calibration |
| Neutrino | Closed vortex loop (distinct topology class from hedgehog) | 🚧 M5.9 / sub-eV |
| Quarks (u, d, s, …) | String endpoint — 1D topological vortex (not 0D hedgehog); Cornell `V(r) = −α/r + σ·r`, `σ ≈ 1 GeV/fm` | 🚧 M5.9 |
| Nucleons (p, n) | 3-string Y-configuration of quark endpoints | DEFERRED → 9d |
| Atomic orbital binding | EM Coulomb + de Broglie standing-wave quantization (NOT a separate force) | → 9d (cross-mass-class) |
| Composites (nuclei, atoms, hopfions) | Multi-defect / knot configurations | DEFERRED → 9d |
| Zitterbewegung clock `ω = 2mc²/ℏ` | Biaxial hedgehog dynamically sources its own twist (M5.6.2b: `ψ = 0` is NOT static = the clock seed); 4D Lorentz negative-energy from `0`-eigenvalue propels it | 🚧 M5.8 (GROUP HEADLINE; build spec `5a §10d` + reorganized roadmap M5.8 gates/breadth/carry-overs/unblocks; faithful kinetic per M5.6.5d) |
| Spin (real, transportable, conserved) | Director orientation = topological invariant; preserved by dynamics | ✅ inherent; empirical anchor: graphene spin-transport |
| First metastable resonance ("first long-lived M5 particle") | `l = 1` harmonic seed on biaxial substrate; extended-lifetime localization (NOT static stability — Q1 RESOLVED) | ✅ **M5.7 COMPLETE (2026-05-28) — the particle is 4D (M5.8), not 3D.** **M5.7.1 (seeded l=1) + M5.7.2 (intrinsic oscillation) = TWO confirmed NULLS:** both seeded and intrinsic 3D orientation energy disperse (apparent N=32 localization washed out at N=48 in both; intrinsic osc frequency 0.25→0.10/t = resolution-dependent, not a converged clock). **Root cause: V confines amplitude but NOT orientation** (M5.6.5c, rotation-invariant). Consistent with Derrick + Duda's framework: the stable particle is the **4D Zitterbewegung clock** (Fig.10), so 3D dispersing is *expected* — M5.7 empirically motivates M5.8 (4D) as necessary. Nuance: topological defect itself permanent (winding conserved); *driven* oscillation is a separate question — **M5.7.3 (9b.1 preview, 2026-05-28) confirmed the EM-like lever sustains a bounded, frequency-selective `(A,ω)` excess (~3× free, resonant at the natural mode), resolution-confirmed** → the driven-thermal substrate for 9b. Energy-conservation bankable (N=48 H-drift 0.01%). (`5a §5h/§5i`, `9b § M5.7.3`) |
| Lorentz contraction / relativistic kinematics | Wave-equation relativistic invariance of the matrix dynamics | ✅ structural; quantitative test = M5.8 |
| Annihilation (Q=+1 + Q=−1 → vacuum + photons) | Topology cancellation event releases stored field energy as outgoing EM waves | ✅ substrate; quantify in 9a |

**Mass → Zitterbewegung target table (M5.8 validation):**

| Particle | Target `ω = 2mc²/ℏ` | Notes |
| --- | --- | --- |
| Electron | 1.55 × 10²¹ rad/s | M5.8 exit gate (within 10%) |
| Muon | 3.21 × 10²³ rad/s | cross-particle test |
| Tau | 5.39 × 10²⁴ rad/s | cross-particle test |
| Neutrino | ~10¹⁵ rad/s | sub-eV (closed vortex loop) |
| Up quark | ~7 × 10²¹ rad/s | string endpoint |
| Down quark | ~1.4 × 10²² rad/s | string endpoint |
| Proton | 2.85 × 10²⁴ rad/s | 3-string Y-config (9d) |

### FORCES

| Phenomenon | M5 numerical mechanism | Status |
| --- | --- | --- |
| Coulomb `1/r` (electrostatic) | Director splay around hedgehog defect; `∇·n̂ ≈ 2/r` far field — **pure topology, no waves needed** | ✅ M5.1 R²=0.978 (1/d relaxation), M5.4 R²=0.9959 (analytical page-18 cross-val), both attractive |
| Maxwell (full EM) | 1-axis tilts of `M`; two independent routes verified — (a) hydro↔EM dictionary `∇×u ↔ B`, `ω×u ↔ E`; (b) Faber `R = Γ × Γ` closed field strength (Maurer–Cartan ⇒ `dR=0` ⇒ homogeneous Maxwell) | ✅ M5.6.4 both routes |
| Strong — confinement (linear) | 1D vortex string, Cornell `V = −α/r + σ·r` (Coulomb + linear); `σ ≈ 1 GeV/fm` per unit length | 🚧 M5.9 (string tension) |
| Strong — short-range (running coupling) | Faber `α_sol(d)`: `R` is non-abelian; Maxwell is long-range abelian limit; short-range = running coupling onset at `r₀` | ✅ M5.6.4 (mechanism); 🚧 M5.9 (calibration) |
| Gravity (linear, GEM) | g-axis boosts on 4D-extended `M`; gravito-EM has same form as Maxwell with `Γ⁰` as boost component | 🚧 M5.8 (4D + SO(1,3) Lorentz) |
| Weak force | PARTIAL — beta decay as topology-reconnection event (defect-class transition: neutron → split → proton + electron + neutrino). No clean SU(2) chiral mechanism yet. | gap → Q10 |
| Magnetism (intrinsic per-defect) | S¹-loop topology (vs S² for electric charge — different topology, same SO(3) parent); the T-component of the defect's outgoing wave (the L+T decomposition) | ✅ mechanism; ⌛ inertial-observability via frequency-downshift = Phase 4 falsifiable target |
| Magnetism (coherent macroscopic) | Coherent T-component summing across aligned defects (permanent magnets) or coherent motion (electromagnets) — reproduces both standard regimes by construction | ✅ structural |
| Atomic orbital "force" | NOT fundamental — EM Coulomb + de Broglie standing-wave quantization | → 9d |

### EM WAVES (9a — photon emergence)

| Phenomenon | M5 numerical mechanism | Status |
| --- | --- | --- |
| Maxwell wave-packet propagation | Tilt modes of `M` propagate at `c` (abelian linear limit of Faber `R`) | ✅ M5.6.4 |
| Outgoing-wave divergence + curl decomposition | Each defect emits a **divergence/radial** part (`∇·n̂` = electric splay — the longitudinal/Coulomb piece in EWT/M3 terms) + a **curl/circulation** part (`∇×n̂` = magnetic — the transverse/radiative piece); the Helmholtz split (M5-native div+curl = EWT's L+T). Shared single `ω = 2mc²/ℏ` Zitterbewegung clock at defect core | ✅ structural; per-defect quantification = 9b.7 |
| Running coupling vs distance | Faber `α_sol(d)`: `‖R‖ ~ 1/r²` far field; `‖R‖·r²` rolls off within `r₀` = onset of non-abelian short-range | ✅ M5.6.4 |
| Photon (quantized EM excitation) | Tilt-mode standing waves on the matrix substrate; quantization mechanism = the 9a research program | 🚧 9a |
| Annihilation → photon emission | Q=+1 + Q=−1 topological cancellation dumps stored field energy into outgoing EM waves | ✅ structural; quantify in 9a |
| Static field vs radiation distinction | Static EM = static topology (Frank elastic, no waves); radiation = dynamic tilt-mode propagation. Different mechanisms in the same field | ✅ structural |

### THERMAL (9b — heat / SABER counterpart)

> **Falsifiable hypothesis:** heat is the joint `(A, ω)` excess of subatomic defect oscillation above its ground state `(A₀, ω₀)`. The energy-dimensional quantity is `(A·ω)²` — `E_kin = ½m(Aω)²`; the bare `A·ω` has units of *velocity* (peak speed) not energy. Falsifiability test = 9b.1.

| Phenomenon | M5 numerical mechanism | Status |
| --- | --- | --- |
| Temperature observable | Joint `(A, ω)` excess; `A = ‖M − D‖_F` amplitude tracker (WAVE_MENU 2), `ω = ‖Ṁ‖_F` clock-rate tracker (WAVE_MENU 3) — already wired in production | 🚧 9b.0 (def + absolute-zero check) |
| Heat = per-defect kinetic energy | `E_kin = ½m(Aω)²` — single-particle quantum-mechanical degree of freedom, NOT ensemble-only statistical | hypothesis → 9b.1 falsifies/validates |
| Modulation lever (engineering trigger) | EM-wave tilt-wave drive on the `(A, ω)` state; AM/FM (5 sub-experiments 10.1a–e) | 🔶 **M5.7.3 PREVIEW (2026-05-28): the lever WORKS** — a continuous EM-like drive sustains the defect `(A,ω)` at ~3× free, bounded steady-state, frequency-selective (resonant at the natural mode); both FM (rings at f_d) + AM (3× amplitude); resolution-confirmed N=32→N=48. Field-theoretic basis for heat-as-driven-excess. Full 9b.1 (10.1a–e) is next. (`9b § M5.7.3`) |
| Soliton-breather correspondence | Compare 9b.1 amplitude oscillations to known Sine-Gordon / φ⁴ breather modes | 🚧 9b.2 |
| Multi-defect amplitude statistics | Seed N ∈ {10², 10³, 10⁴} defects, run to equilibrium → Boltzmann or Bose-Einstein distribution | 🚧 9b.3 |
| Specific heat `C_V(T)` | Derived from 9b.3 statistics; reproduce Dulong-Petit (high T) + Einstein-Debye (low T) + electronic scaling | 🚧 9b.4 (stiff prediction — wrong scaling = hypothesis falsified) |
| Blackbody spectrum | Thermalized-ensemble emission; Wien displacement law + Stefan-Boltzmann scaling | 🚧 9b.5 |
| Phase-coherence transition at low T | Defect-ensemble transition to coherent ground state (superconducting / BEC analog); reproduce critical temperatures for known materials | 🚧 9b.6 |
| Outgoing-wave thermal content | Thermal excess expressed in outgoing wave's amplitude/frequency/polarization at engineering distance from defect core — the addressable channel for indirect modulation | 🚧 9b.7-.8 |
| Heat ↔ magnetism wave-level coupling | Thermal excess scales the T-component of the outgoing wave — explains Curie temperature, magnetocaloric, pyromagnetic from one substrate | 🚧 9b.9 |

### What 9a / 9b / 9c share — same dynamics, three readings

| Reading | What it measures | Role of `ω` |
| --- | --- | --- |
| 9a (EM waves) | Outgoing wave at frequency ω = photon | ω = wave frequency |
| 9b (Thermal) | `(A, ω)` excess of defect oscillation | ω = thermal frequency |
| 9c (Time / clock) | Zitterbewegung `ω = 2mc²/ℏ` | ω = intrinsic local clock |

Shared infrastructure: `(A, ω)` trackers + outgoing-wave measurement + EM viz (`∇·n̂` / `‖∇×n̂‖` meshes + E/B glyphs). Build once, each stage reads a projection. See `0b_M5_roadmap.md § Parallel research stages` for the full link.

### CURRENT LIMITATIONS — what M5 does NOT yet do

| Limitation | Phase that addresses it |
| --- | --- |
| No FREE metastable resonance in 3D — M5.7.1 (seeded) + M5.7.2 (intrinsic) both disperse; the stable particle/clock is 4D (`5a §5h/§5i`) | M5.8 (4D Zitterbewegung clock — the actual stable particle). Driven oscillation (9b.1 lever) is the separate near-term test (M5.7.3 preview). |
| No Zitterbewegung frequency measurement | M5.8 (faithful kinetic on `O(x) ∈ SO(3)`, per M5.6.5d) |
| Lepton mass ratios not calibrated to PDG | M5.9 |
| Gravity emergence not run (still 3D) | M5.8 (4D + SO(1,3) Lorentz) |
| Quark string + Cornell `σ ≈ 1 GeV/fm` not validated | M5.9 |
| Photon quantization mechanism not formalized | 9a |
| Thermal hypothesis not falsified/validated | 9b.1 (SABER engineering trigger) |
| Weak force — no clean SU(2) chiral mechanism | Q10 (long-tail open) |
| Inertial observability of magnetism via frequency-downshift | Phase 4 falsifiable target |
| Heat ↔ magnetism wave-level cross-validation | 9b.9 |
| Foundational stance vs Bell / Kochen-Specker not formally settled | Q12 (foundational; defense template in `1b § Foundational stance`) |
| Composites (nuclei, atoms, hopfions) deferred | 9d (post-9b) |
| Exact biaxial-stable `V(M)` coeffs `Λ = (1, δ, 0)` | Q7 → M5.9 (M5.6 ships b=0 confinement interim) |
| Exact Faber running-coupling profile (vs imposed `arctan`) | Q8 → M5.9 (M5.6 ports the validated mechanism) |

---

## HARDEST PIECES — status board

The long-running hardest-pieces tracker for M5 (mirrors M6's). These are the load-bearing unknowns the framework must resolve, distinct from the discrete Q-numbered questions above — a piece can persist across many phases.

| Hardest piece | Gates | Status (2026-05-27) | M6 lend? |
| --- | --- | --- | --- |
| Matrix substrate migration (Vector(3) ψ → `M = ODO^T`) | M5.4 ✅ | ✅ **COMPLETE (2026-05-26).** `ti.Matrix.field(3,3)` triple buffer + `eigen_decompose` lynchpin (director recovery 0.9995) + matrix seeders + `‖M−D‖_F`/`‖Ṁ‖_F` trackers. **M5.1 Coulomb reproduced on M** — R²=0.9704 relaxed + R²=0.9959 analytical page-18 cross-val, both attractive. Operators verified vs analytic. Rendering re-sourced from M (on-screen verified); live wave path retired (ψ engine dormant legacy). | — |
| **V(M) potential form** (Duda's #1 limitation, Q7) | M5.6 ✅ | ✅ **M5.6 (2026-05-27).** Eq.13 LdG implemented + V-on confinement live (`LDG_STIFFNESS_K`). Key finding: 3-term form has **no biaxial minimum** → M5.6.5c uses a `b=0` amplitude well (confines `Tr(M²)`, leaves δ flat). Exact biaxial-stable `Λ=(1,δ,0)` form = Q7, Duda-open → M5.9. | ✅ strong (see below) |
| Regularization (Faber) to activate V(M) + running coupling (Duda #2, Q8) | M5.6 ✅ | ✅ **M5.6 (2026-05-27).** Faber MTF ported: `E₀=α_f ℏc·π/(4r₀)` reproduced, `Λ=q₀⁶/r₀⁴` mapped onto M, **mass pinned `E∝1/r₀`**, `R=Γ×Γ`→Maxwell + running-coupling onset at `r₀`. Mass knob is now physical `r₀`. Exact `arctan`-profile re-derivation → M5.9. | 🔶 partial |
| First metastable resonance (no static soliton) | M5.7 ✅ → M5.8 | ✅ **M5.7 COMPLETE (2026-05-28) — free vs driven split.** Free defect disperses (M5.7.1 seeded + M5.7.2 intrinsic, both NULL; root cause V confines amplitude not orientation, M5.6.5c) ⇒ the stable particle is the **4D Zitterbewegung clock (M5.8)** (3D dispersing is expected per Derrick — M5.7 empirically motivates the 4D extension). **Driven defect SUSTAINS** a bounded, frequency-selective `(A,ω)` excess (M5.7.3 9b.1 preview, EM-lever ~3× free at the resonant natural mode) ⇒ the driven-thermal substrate for 9b. Energy-conservation bankable (N=48 H-drift 0.01%). `5a §5h/§5i`, `9b § M5.7.3`. | ✅ strong (→ M5.8 / 9b) |
| KG-from-twist emergence | M5.6 ✅ | ✅ **DONE (2026-05-27).** KG mass is **GEOMETRIC** — minimal coupling to the hedgehog connection `Â`, the explicit mass term cancels, core regularization generates `mass²(r)=3r_c²/(2 reg²)` (Fig.9 reproduced, `5a §5a-c`). Biaxial hedgehog port + z-axis disclination handling done; defect dynamically sources its own twist (M5.8 clock seed). | — |
| Zitterbewegung clock / 4D negative-energy propulsion | M5.8 | Mechanism known (Fig.10); toy-model numerics available (slide p.33). **Externally validated 2026-05** as THE hard problem: group consensus that ALL macroscopic analogues fail at intrinsic non-damping oscillation (Couder droplets need a shaker; spinning needles damp like a pendulum) — only the elementary-particle clock is intrinsic + non-damping, so propulsion must be intrinsic. Validates the 4D negative-energy direction. See [9c § clock propulsion](9c_time_dynamics.md#the-de-broglie-clock-is-intrinsic--group-consensus-on-clock-propulsion-2026-05). | 🔶 partial |
| Biaxial eigenvalue hierarchy → lepton masses (Duda #3/#9/#10) | M5.9 | Calibration parameters; post-SABER-trigger. | — |
| Weak-force clean SU(2) mechanism (Duda #7, Q10) | — | GAP. Beta decay as topology reconnection is partial. Out of SABER scope. | — |
| Deeper substrate beneath M (Duda #8, Q9) | — | PARKED. Matches granule-level picture. | — |

### Duda's stated limitations of his own LdGS / M5 model

Physics limitations Duda (or his papers, as quoted in our notes) flags about his own framework. Sourced from the M5 research docs; reproduced from `NEPTUNYA/SABER/1_model_selection.md` (physics content — the SABER engineering-timeline analysis stays in that doc).

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

SABER engineering-timeline impact is analyzed in the `1_model_selection.md` (§ "How the limitations affect SABER timeline"); here we keep only the OpenWave-physics phase-gating. The chain `M5.5 → M5.6 → M5.7 → 9b` is the documented bridge to the SABER thermal counterpart (see [`9b_thermal_energy.md`](9b_thermal_energy.md)).

| Limitation | Gates | Reading |
| --- | --- | --- |
| #1 V(M) form | M5.5 | Direct blocker for the paper Lagrangian. Until V(M) is pinned, defect dynamics under perturbation aren't uniquely simulable. |
| #2 regularization | M5.6 ✅ | The hard step — **ported cleanly (2026-05-27)**. Faber's MTF mapped onto M, mass pinned `E∝1/r₀`, `R=Γ×Γ`→Maxwell. M5.6 unblocked; exact running-coupling form → M5.9. |
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
| Sandbox-first: prototype in scipy, validate, *then* port to production | Prototype the M5 V(M) + matrix-defect-profile BVP in a scipy M5-sandbox before the Taichi M5.6+ kernels (the path M5.5 followed: `sandbox_v5` scipy/numpy → Taichi port) |

Canonical M6 numerical recipe (charged sector via 2-function IVP, neutral sector via the BVP above) is consolidated in `m6_ouroboros/research/0d_canonical.md` — the reference if we want to mirror the methodology.

### M6 strengths available to M5 regardless of M6's gate outcomes

Imported from `m6_ouroboros/research/0b_model_gates.md` § "What M6 offers regardless of gate outcomes" and `1_model_selection.md`:

| M6 strength | Relevance to M5 |
| --- | --- |
| Neutral chaoiton **ground state reachable** via BVP (v9 phase 2) | Proof-of-method that a topological soliton's ground state IS findable numerically — encouragement + recipe for M5.7's metastable resonance hunt |
| DM mass prediction **published** (DM paper v2, Zenodo 20350105, 2026-05-22): `m_χ = 0.460 MeV`, `m_J = 0.618 MeV`, `C = 770 MeV·fm` (canonical g=1.0, B0=0.5; `0d_canonical.md`). Supersedes the earlier v8 scan numbers (0.998/1.033/6.7e-4) | Independent topology-as-particles result; cross-validates the shared "no static soliton / time-periodic" premise M5 also relies on |
| Lepton spectrum at PDG precision via 2-function reduction | Demonstrates the ansatz-reduction discipline (collapse DoF before solving) that M5.5 can borrow for the matrix field |
| Charge quantization via Chern-Simons / Hopf (Lean-proved) | Different topological flavor (Hopf-linking) from M5's Brouwer-winding — but confirms charge-from-topology is robust across frameworks |
| Methodology cross-pollination (η, `l=1` BC, virial diagnostic) | Listed above — the concrete transferable pieces |

**Caveat (load-bearing):** M6 is **not** the thermal substrate and does not supply M5's `(A, ω)` thermal degree of freedom — that lives on M5's matrix field. M6's contribution to M5 is **method + credibility cross-validation**, not physics substrate. A hypothetical M6-substrate thermal cross-check would require a full Taichi production port (ensemble + external-drive kernels), which is NO-GO / deferred indefinitely.

---

## Notes on scope

- This tracker covers **M5 physics/framework questions for the group** (Duda / Close / Yee) plus the hardest-pieces board. M5.4 implementation decisions (Taichi storage layout, eigen-kernel design, granule/glyph UX calls) are tracked as **roadmap tasks** in [`0b_M5_roadmap.md`](0b_M5_roadmap.md) and [`4b_rendering_features.md`](4b_rendering_features.md), not here — same split as M6 (questions here, sandbox tasks in the work log).
- M5 has **no active email round in flight** (unlike M6's intense Werbos/DeepSeek cadence). Duda's last substantive reply was 2026-05-15. **M5.6 (complete 2026-05-27) now gives a strong results bundle to seed a Q7/Q8 outreach**: the Eq.18 action runs in production + is energy-conserving, KG mass is geometric (Fig.9), Maxwell recovered both routes, Faber ported with mass pinned `E∝1/r₀`, and the **V(M)-is-rotation-invariant finding** (V acts only on the eigenvalue/regularization sector, not the twist) + the **3-term-LdG-has-no-biaxial-minimum** result sharpen the Q7 ask — Duda's exact `Λ=(1,δ,0)` Eq.13 coeffs + Faber's exact running-coupling scheme (Q8) are the natural questions, now feeding M5.9 calibration. **M5.7 (complete 2026-05-28) adds a further outreach hook** — the free-disperses/driven-sustains split (the particle is 4D; the EM lever sustains a bounded `(A,ω)` excess) is a clean result to share with Close (Q11) + Duda, no longer a question. None blocks current work (M5.7 complete; next = M5.8 or 9b.1, Rodrigo's call).
