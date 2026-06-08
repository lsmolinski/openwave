# 9b — Thermal Energy as Wave Dynamics (Hypothesis)

**Parallel research stage** — scientific counterpart to SABER's Direct Heat Conversion (SABER MAIN GOAL). This document holds the **physics hypothesis** and the **scientific test program** (9b.0 → 9b.9). Engineering implementations stay in the SABER repo-discipline rule (`feedback_repo_discipline`). 9b in this public repo validates the hypothesis numerically; SABER repo holds device-side engineering.

Unblocks once M5.7 (resonance hunt) lands; runs in parallel with M5.8 (Zitterbewegung) / M5.9 (lepton families). Empirical falsification (or validation) of the SABER thermal-amplitude hypothesis.

---

**Hypothesis**: Thermal energy is subatomic-based (thermal energy = particle excess angular momentum) rather than kinetic-based at the fundamental level, reframing heat entirely.

**Background**: Thermodynamics defines heat as molecular vibrations — mechanical movement, a classical phenomenon. OpenWave proposes that heat is a **quantum mechanics property**: a modulation of the fundamental energy wave's A/λ within a particle's surrounding standing wave structure. This is not bulk kinetic motion — it's a wave-state change at the subatomic scale.

**Research goal, Quantum Thermodynamics concept**: Investigate thermal energy's true nature — kinetic vs wave-based. The hypothesis that heat relates to standing wave dynamics within particle radius is testable in simulation by comparing predictions between the two models. If thermal energy is encoded in standing wave amplitude or frequency modulation rather than particle velocity, the wave-based model should produce different (and potentially more accurate) predictions for phenomena like specific heat, thermal conductivity, and blackbody radiation.

**Connection to wave resonance**: thermal energy may be **increased wave steepness** (A/λ ratio) beyond the particle's fundamental — storing more energy in standing waves. The particle naturally tends to return to its fundamental steepness, releasing excess as radiation (photons). Thermal equilibrium between two bodies is steepness exchange through resonance. Steepness conservation (A/λ = const) holds for isolated energy redistribution; external energy input (heating) or output (radiation) changes steepness. Photon emission is the resonance-mediated frequency readjustment described by Wolff's coupled-oscillator model. See [Wave Resonance](9c_time_dynamics.md#wave-resonance-as-the-fundamental-energy-exchange-mechanism) for the full analysis.

**Connection to topology**: in the M5 framework, the thermal-energy hypothesis evolves into a per-defect framing: each topological defect has an intrinsic ground-state oscillation amplitude (its Compton wavelength `A_0 ≈ ℏ/(mc)`) and any excess amplitude above `A_0` is the thermal degree of freedom. The wave-mechanical "increased steepness" picture and the "per-defect amplitude excess" picture are the same physics described at different abstraction levels. See [1b_topological_defect.md](1b_topological_defect.md) for the topological framing.

## Heat and Light: The Fundamental Conversion Cycle

Heat and electromagnetic radiation are in constant mutual conversion. This is observable:

- **Blackbody radiation**: any particle that isn't reflecting light is always emitting EM waves — from its stored thermal energy. This is universal — every object above 0 K radiates
- **Infrared emission**: heated objects emit IR radiation (captured by thermal cameras). The hotter the object, the shorter the peak wavelength (Wien's law)
- **Absorption**: EM waves hitting matter convert to heat (microwave → water heating, sunlight → surface warming)

Heat converts itself to light and light converts itself to heat. There must be a resonance / energy exchange relationship mediating this — Wolff's coupled-oscillator model describes exactly this: wave systems exchanging frequency (energy) through resonance until they reach equilibrium. The photon is not a separate entity "emitted" — it is the traveling wave disturbance produced when a particle's standing wave steepness relaxes toward fundamental.

## The Snap In / Snap Out Mechanism

The fundamental conversion in the universe is between two wave states:

- **Standing waves** = matter and heat (energy stored, localized, "snapped in")
- **Traveling waves** = movement and light (energy flowing, propagating, "snapped out")

The physics is the conversion **between** standing and traveling waves. Every energy process in the universe is one of these conversions:

- **Combustion**: chemical bonds (standing waves between atoms) snap out → traveling waves (heat + light + kinetic)
- **Nuclear**: nuclear bonds (standing waves between nucleons) snap out → massive traveling wave release
- **Photon emission**: thermal standing wave energy snaps out → traveling EM wave
- **Absorption**: traveling EM wave snaps in → standing wave energy increase (heating)
- **Matter creation**: traveling waves snap in → stable standing wave structure (particle)
- **Annihilation**: standing waves snap out → traveling waves (gamma photons)

Understanding this snap mechanism at the wave equation level is fundamental to understanding how energy moves between matter and radiation — and is the natural connection point between this Layer-6 research and the M5 topological-defect physics, where "snap" events become amplitude-transfer events at the per-defect level.

### Thermal-mechanics detail

A working hypothesis for the heat output domain: **the thermal degree of freedom is the joint amplitude+frequency `(A, ω)` excess of a topological defect's intrinsic Zitterbewegung above its ground-state values `(A₀, ω₀)`** — both components carry thermal content; they are coupled by **wave-steepness conservation** (`A/λ = const`), the same conservation principle articulated in [9c_time_dynamics.md](9c_time_dynamics.md). Each topological defect (= each particle in M5) has its own intrinsic ground-state oscillation `(A₀, ω₀)` with `A₀ ≈ ℏ/(mc)` (its Compton wavelength) and `ω₀ = 2mc²/ℏ` (the Zitterbewegung clock). Any excess on top of *either* component is what aggregates statistically into macroscopic thermodynamic temperature. This extends thermodynamics from ensemble-only statistical mechanics (where heat is a collective property requiring many particles) to a per-particle quantum-mechanical degree of freedom — the same way M5 makes other things "more fundamental, not less" by deriving collective phenomena from defect-level dynamics.

**Refinement — the ground-state oscillation is a constant-amplitude spin (2026-05-30).** The intrinsic Zitterbewegung is a *steady rotation* of the defect's internal frame (`ψ = ωt`, a constant-radius circular motion — like a point on a wheel), **not** a pendulum: `A₀ ≈ ℏ/(mc)` is the **fixed rotational radius** (a spin has amplitude, just constant), `ω₀ = 2mc²/ℏ` the rate. So the joint `(A, ω)` excess is **radius-growth (`A − A₀`, AM) + rate-shift (`ω − ω₀`, FM)** above that spinning baseline — observable as the launcher's WM2 "Thermal Amp" (`‖M−D‖`) and WM3 "Thermal Clock" (`‖Ṁ‖`). It spins rather than swings because `V(M)` is rotation-invariant — no restoring force on the orientation phase (M5.6.5c) — which is *also* why a free defect's excess disperses and the excess must be **driven** to be maintained (the M5.7.3 / 9b.1 framing below). Substrate-physics framing per OpenWave `0c §L7`.

Two consequences worth noting for the 9b observables: (1) **the radius `A` is the eigenvalue gap** — the twisting `b`–`c` block `R(ψ)·diag(δ,0)·R(ψ)ᵀ` makes the components oscillate with amplitude `(λ_b−λ_c)/2 = δ/2`, so thermal-`A` (WM2) is the *growth* of that `δ/2` excursion and the QM-`δ` axis sets its scale; (2) **the observable cycles at `2ω`** (apolar `n̂⊗n̂` doubling), the likely origin of `ω₀ = 2mc²/ℏ` vs the phase rate `mc²/ℏ` (confirm in M5.8). And **what oscillates is a collective, phase-locked mode of the orientation field localized on the defect** — not the defect-as-a-point, not independent voxels — so the 9b temperature observable (9b.0) and the `(A, ω)` modulation (9b.1) act on that whole defect-localized mode (and its outgoing wave, 9b.7), not on isolated voxels.

**Direct prediction at the boundary**: if all defects sit exactly at `(A₀, ω₀)` with zero thermal excess, the framework predicts the system temperature is exactly **0 Kelvin** (absolute zero). Conversely, the temperature observable is — by construction — the measurement of joint `(A, ω)` excess above ground state. This recovers absolute zero correctly and consistently with the third law of thermodynamics, BEC, and superconductivity (phase-coherent ground states are the zero-excess configuration). The framework being self-consistent at this boundary is a pre-numerical sanity check on the hypothesis itself — and one of the 9b sub-stages below (9b.0) verifies it numerically.

The closest precursor framing in OpenWave is the steepness-conservation / energy-starvation work in [9c_time_dynamics.md](9c_time_dynamics.md). The topological-defect deep dive in [1b_topological_defect.md](1b_topological_defect.md) covers the defect's ground-state oscillation; this Phase 7 pathway is what extends it to the excited (thermally-excited) regime.

**Connection to Time Dynamics**: the Zitterbewegung frequency `ω` IS the defect's intrinsic local clock — its rate of internal cycling. Modulating defect `ω` therefore corresponds to *locally engineering the rate of time at the subatomic scale*. The Phase 7 (A, ω) modulation experiments are simultaneously thermal-mechanics validation AND time-dynamics validation under different framings — same physical operation, two complementary scientific framings. See [9c_time_dynamics.md](9c_time_dynamics.md) for the time-dynamics-side treatment.

**Infrastructure foundation from M5.1 task 6**: the `lagrangian_engine.relax_director_step` kernel built for M5.1 task 6 (gradient-descent ground-state finder) is the **mathematical precursor** to 9b's thermal-modulation primitives, despite the very different physics intent. The relaxation uses tangent-projected diffusion `∂_τ n = ∇²n − (n·∇²n)n` to drag the field toward the minimum-energy configuration — i.e. the **γ → ∞ limit** of a damped wave equation `∂²_t ψ = c²·∇²ψ − γ·∂_t ψ` (pure-dissipation overdamped regime). 9b's thermal-modulation kernels will reuse the **same primitives** — the Laplacian stencil (`compute_laplacian`), the tangent projection (to preserve `|n|=1`), the soft-core pinning (to preserve topology) — but with a **physical γ** representing real damping (radiation, phonon coupling, EM-load impedance) instead of an artificial time-stepping parameter τ. So although task 6's relaxation has no physics interpretation in M5.1 (it's a ground-state finder for visualization + initial conditions), the kernel and its surrounding test infrastructure (`research/scripts/m5_1_relax.py`, `pin_centers/signs` plumbing in `SimulationState`, the `relax_field` helper that updates both `psi_am` AND `psi_prev_am` to preserve ψ̇=0) will all carry forward into 9b with minimal modification. The applied-tech implications of this math overlap — engineering γ as the DHC engineering target — are explored in SABER work outside this scope.

#### 9b sub-stages — Thermal mechanics from defect (A, ω) statistics

| Sub-phase | Layer / domain | Headline test |
| --- | --- | --- |
| **9b.0 / Temperature observable + absolute-zero sanity check** | Per-defect thermal (definition + boundary case) | Define the temperature observable as a measurement of joint `(A, ω)` excess above ground-state `(A₀, ω₀)`. Implement it as a numerical probe operating on the M5.7 metastable defect output. Verify the boundary case: a defect sitting at exactly `(A₀, ω₀)` with no perturbation reads `T = 0 K`. Verify monotonicity: small symmetric `(A − A₀, ω − ω₀)` excitation reads small positive `T`; larger excitation reads larger `T`. Cheapest pre-numerical sanity check on the framework — if 9b.0 fails, the temperature definition is wrong and all downstream 9b.x stages are meaningless. Establishes the temperature observable that 9b.1 → 9b.6 will all use |
| **9b.1 / Single-defect (A, ω) modulation** | Per-defect thermal | Take a single biaxial hedgehog from M5.7. **Bidirectional perturbation of both amplitude AND frequency, from the defect's *current* state** (not only from topological ground — a real defect's instantaneous state is `(A₀, ω₀)` *plus* current thermal excess; modulation perturbs from wherever it is). Five sub-experiments share the same kink + KG background: **(10.1a) AM-up** — bump A above current; measure relaxation, breather modes, emission spectrum. **(10.1b) AM-down** — pull A below current via destructive interference; measure whether the field refills the kink (a per-defect cooling event). **(10.1c) FM-up** — drive ω above current via tuned external wave; measure A response per wave-steepness conservation `A/λ = const`. **(10.1d) FM-down** — slow ω below current via tuned interference; measure A rise. **(10.1e) coupling cross-check** — validate `A/λ = const` under all four perturbation directions. Cheapest go/no-go on the joint (A, ω) thermal hypothesis |
| **9b.2 / Soliton-breather comparison** | Per-defect thermal | Compare 9b.1's amplitude oscillations to known soliton-breather modes (Sine-Gordon, φ⁴). Breather modes should match observed thermal excitation patterns. Cross-validation against established field-theory math |
| **9b.3 / Multi-defect amplitude statistics** | Defect ensemble thermal | Seed N defects (10², 10³, 10⁴) with varying initial amplitudes. Run to thermodynamic equilibrium. Extract amplitude distribution. Predict: should match Boltzmann (classical) or Bose-Einstein (quantized) statistics for the appropriate temperature definition |
| **9b.4 / Specific heat reproduction** | Defect ensemble thermal | From 9b.3 statistics, derive specific heat `C_V(T)`. Validate against experimental: Dulong-Petit at high T, Einstein-Debye at low T, electronic heat-capacity scaling for free electrons. Stiff prediction — wrong scaling = hypothesis falsified |
| **9b.5 / Blackbody spectrum** | Thermal-EM coupling | Emission spectrum from a thermalized defect ensemble. Validate Wien's displacement law + Stefan-Boltzmann scaling. The heat → light channel |
| **9b.6 / Phase-coherence transition** | Quantum thermal | At low T, do defect ensembles transition into phase-coherent ground states (analogous to superconductivity / BEC)? Reproduce critical temperatures for known materials |
| **9b.7 / Per-defect outgoing-wave thermal-character measurement** | Per-defect thermal-EM coupling | Take a single defect from 9b.1 at varied joint (A, ω) states. Measure the outgoing wave-field's amplitude / frequency / polarization at macroscopic distance from the defect core. Confirm thermal-excess content is also expressed in the outgoing wave (the dual-located thermal-degree-of-freedom prediction — heat lives both inside the defect AND in the spherical waves the defect propagates outward) |
| **9b.8 / Wave-modulation back-reaction** | Per-defect thermal | Drive resonance with the defect's outgoing wave-field at macroscopic distance (perturbation acts on the wave, NOT on the defect interior). Measure whether the defect's joint (A, ω) state responds via back-reaction. Cheapest go/no-go on whether outgoing-wave manipulation is physically viable as an alternative to direct defect modulation |
| **9b.9 / Heat-magnetism wave-level co-scaling** | Per-defect thermal-EM coupling | Cross-validation joining 9b.7 with the Phase 4 magnetic-emergence work. Confirm that thermal excess scales the latent magnetic-component (T-component) magnitude of the outgoing wave as predicted by the L+T-as-channels framing |

These 9b sub-stages are deliberately speculative — they assume the joint (A, ω) framework holds and cascade from there. If 9b.1 falsifies the hypothesis, the rest don't run as planned. But if 9b.1 confirms it, OpenWave has identified a new mechanism for thermal physics: heat as a single-particle quantum-mechanical phenomenon rather than an ensemble-only statistical one.

> **Why 9b.7 / 9b.8 / 9b.9 matter**: 9b.1's sub-experiments perturb the defect directly, but any real instrument operates at macroscopic distance from the defect, on the outgoing wave-field — never directly on the sub-fm-scale defect interior. 9b.7 verifies that the outgoing wave carries the thermal information; 9b.8 verifies that perturbing the outgoing wave back-reacts on the defect's state (closing the indirect-manipulation loop); 9b.9 cross-validates the per-defect heat-magnetism wave-level coupling that joins thermal (9b) to EM emergence (Phase 4) at the per-defect substrate. Together they upgrade the program from "we can perturb defects we can't reach" to "we can perturb something we *can* reach (the outgoing wave) and validly infer / predict the defect-level effect."

**9b.1 four-outcome decision matrix**:

| AM (10.1a + 10.1b) | FM (10.1c + 10.1d) + coupling (10.1e) | Interpretation |
| --- | --- | --- |
| ✅ | ✅ + `A/λ = const` validated | Joint (A, ω) hypothesis validated; proceed to 9b.2 → 9b.9 |
| ✅ | ❌ | Revise: heat is amplitude-only — update the physics framing accordingly |
| ❌ | ✅ | Revise: heat is frequency-only or coupling inverted; rethink the physics framing |
| ❌ | ❌ | Hypothesis falsified at the per-defect level; close out 9b as informative negative |

**Two-tier prerequisite distinction for 9b.1**:

| Test scope | Minimum prerequisites | Why |
| --- | --- | --- |
| **9b.1 internal validation** (tests the hypothesis at the field-theoretic level) | M5.7 (metastable defect) + M5.6 (KG wave dynamics from twist) | The "external wave" driving 10.1c/d can be a KG perturbation in the M5 medium. Sufficient to test the (A, ω) modulation hypothesis |
| **9b.1 full-lab-relevance** (tests with the same wave classes a physical drive would use — EM, magnetic, acoustic) | M5.7 + Phase 4 (EM emergence) + Phase 5 (gravity, optional for non-gravitational drives) | A physical drive uses real-world wave classes. To predict its effect on defect modulation, the simulation drive must be the same class as the lab drive |

The internal validation can run as soon as M5.7 lands; full lab-relevance results follow Phase 4. Heat-substrate validation (9b.x) and lever-class validation (Phase 4 / Phase 5) are independent prerequisites — both must reach completion before the combined physics can be predicted reliably.

#### ⚠️ M5.7 finding (2026-05-28) — what it means for the 9b `(A, ω)` source

M5.7.1 (seeded l=1 perturbation) and M5.7.2 (the defect's intrinsic oscillation) both showed that **a *freely-evolving* 3D defect DISPERSES its excess orientation energy** — it does not self-sustain a localized coherent oscillation. Root cause: the `V(M)` amplitude well confines `Tr(M²)` but is rotation-invariant, so it does NOT confine the director *orientation* (M5.6.5c); orientation/twist energy radiates away. The stable *intrinsic* clock requires the **4D** Lorentz mechanism (M5.8), not 3D dynamics. This sharpens — does **not** invalidate — the 9b thermal program, on two counts:

| Implication for 9b | Why |
| --- | --- |
| **The topological defect itself is permanent.** Winding is conserved — what disperses is the *excess* `(A, ω)` oscillation energy, not the particle. So there is always a defect to carry a thermal state. | Topology (Q1/Q3) — the defect cannot unwind. |
| **9b's `(A, ω)` excess must be DRIVEN/sustained, not free.** A free defect sheds its excess; a *continuously driven* one (the EM-modulation lever — exactly 9b.1's 10.1a–e) can hold a steady-state `(A, ω)`. **This is the right framing for heat anyway** — thermal energy is *maintained* excess (a bath/drive), not a one-shot free oscillation that radiates away. | A free-dispersal null does not answer the driven-steady-state question; 9b.1 *is* the driven test. |
| **The intrinsic-clock `ω` (for `T = 0 K` ground-state `ω₀`) comes from M5.8 (4D), not M5.7.** | M5.7 established 3D has no converged free clock frequency (it shifted 0.25→0.10/t with resolution). |

**Direct consequence — the M5.7.3 / "9b.1 preview" (driven-oscillation) experiment is the natural bridge:** apply a continuous EM-wave-like drive to the defect and test whether it sustains a steady-state `(A, ω)` (vs the free-evolution dispersal). If it does, the joint-`(A, ω)`-as-heat hypothesis has its substrate confirmed at the field-theoretic level and 9b.1's full sub-experiments (AM/FM 10.1a–e) follow.

#### ★ The MOLTEN-CLOCK bridge — the first in-data contact between the M5.8 clock and the heat hypothesis (N-6d, 2026-06-07)

The M5.8 N-6d chaos battery (`sandbox_vn/m5_8_2n_chaos_battery.py`) classified the saturated quartic breather as a **molten clock**: a persistent coherent fundamental (ω₁ — the N-1 attractor) riding on low-dimensional CHAOTIC dressing (λ_max ≈ +0.4–0.7 per τ, D₂ ≈ 2.7–3.0 — dynamical, not additive noise, by the control ladder). Put against this hypothesis:

| Heat hypothesis statement | N-6d measurement |
| --- | --- |
| Heat = the joint `(A, ω)` excess of defect oscillation ABOVE its ground state | The excited breather (H_dyn ≈ 2.8× H_rest) hosts its excess energy as the chaotic dressing AROUND the coherent clock |
| The excess is the per-defect kinetic content (temperature's substrate) | The dressing is quantifiable: λ_max / D₂ / broadband power fraction = a **temperature-like axis for a single defect** |
| Less excitation ⇒ less heat | The least-excited arm (R1-settled) is markedly less molten (D₂ 1.70 vs 2.7–3.0) — the first two points of a per-defect "melting curve" T(E) |

**Candidate dictionary (hypothesis to test, registered 2026-06-07): molten-ness IS the microscopic face of per-defect heat** — cold particle = regular clock, hot particle = molten clock; the M5.7.3 drive pumps exactly this axis. **N-6c CONFIRMED IT (2026-06-07, the deep-settle route): the clock REGULARIZES toward the ground state** — the coldest measured state (T = 2.36) reads near-regular (λ_max +0.110 ≈ the estimator bias floor, D₂ 1.68 vs 2.7–3.0 hot): the picture closes — a (near-)clean clock at rest + a measurable thermal (molten) sector when excited, the clock and its heat in one object. Bonus preparation-physics: a cold defect CANNOT be made by gentle kicking (the spontaneous release floors excitation at ~1.7× H_rest) — cooling requires DAMPING, and FULL coldness is forbidden by un-sittability itself (the self-start). The melting curve now has both anchors (cold λ+0.11/D₂1.68 ↔ hot λ+0.75/D₂3.0). 9b.1 inherits the battery (λ/D₂/broadband fraction) as its temperature observables alongside WM2/WM3. Honest caveats: stepper-stiffness numerical heating must be separated from physical molten-ness (dt-cross-check of λ/D₂ = the hardening item); single defect, natural units.

**The M5.7 arc, resolved (2026-05-28) — free disperses, driven sustains:**

| Sub-step | Question | Result |
| --- | --- | --- |
| M5.7.1 (seeded, free) | does an injected resonance persist? | **NULL** — disperses |
| M5.7.2 (intrinsic, free) | does the defect self-sustain a clock? | **NULL** — disperses → the clock is **4D (M5.8)** |
| M5.7.3 (driven) | can an EM lever hold an excited `(A, ω)`? | **YES** — 3.0×, bounded, frequency-selective |

So the **free** defect can't self-sustain (that's the 4D Zitterbewegung clock's job — M5.8), but a **driven** defect holds a frequency-tunable, bounded `(A, ω)` excess — which is **exactly what thermal energy is in the SABER hypothesis** (a *maintained* excess from a bath / drive) and **what the EM modulation lever does**. 9b.1 has its field-theoretic foundation, with the resonant modulation frequency identified.

#### ✅ M5.7.3 — 9b.1 driven-oscillation PREVIEW: the lever works (2026-05-28)

`sandbox_v7/m5_7_3_driven_oscillation.py`. A continuous EM-wave-like drive — a fixed localized
director-rotation forcing `F_drive(x,t) = A_drive·sin(2π f_d t)·D̂(x)`, `D̂ = w(x)·[G_y, M_bg]`
(the so(3) rotation tangent, symmetric+traceless = orientation-sector / V-flat, shell-localized to
the defect — an incident-EM-tilt-wave proxy, the same channel M5.6.4 maps to Maxwell) — added to the
Eq.18 acceleration. Metric: `A_core(t) = ⟨‖M − M_bg‖_F⟩` over the defect shell (the excitation
amplitude), sustained = mean over the last third; plus the core-director FFT response at `f_d`.

**Result — a driven defect sustains a bounded, frequency-selective `(A, ω)` excess (resolution-confirmed):**

| Drive `f_d` | gain vs free, N=32 | gain vs free, N=48 | response @ `f_d` | H-growth | plateau |
| --- | --- | --- | --- | --- | --- |
| free (0.0) | 1.0× (A=0.31) | 1.0× (A=0.275) | — | — | — |
| **0.10 (resonant)** | 2.7× | **3.0× (A=0.82)** | 59% | +1% (bounded) | ✅ rises then flattens 0.78→0.84→0.81 (even ticks down 0.837→0.813 at the end) |
| 0.20 | 1.6× | 1.5× | 94% | +0% | — |

- **Sustained + resolution-stable** — the gain *grew* 2.7→3.0× from N=32→N=48 (it did NOT wash out, unlike the free M5.7.1/.2 nulls). The lever genuinely holds the defect excited.
- **Genuine steady-state plateau** — `A_core` at `f_d=0.10` rises then flattens at ~0.82 (drive-in balancing dispersal-out), H-growth +1% ⇒ a *maintained excess*, not a runaway. This is the thermal picture: heat = maintained `(A, ω)` excess from a continuous drive.
- **Frequency-selective resonance** — 3.0× at `f_d=0.10` vs 1.5× at 0.20, peaking at the defect's **natural mode** (`f≈0.10/t`, exactly the M5.7.2 N=48 intrinsic frequency). The lever couples best on resonance — **this is the modulation frequency a real EM lever should use** (the 9b.1 FM handle).
- Both `(A, ω)` levers demonstrated: **ω** (the core director responds at the drive frequency — FM) and **A** (3× amplitude gain — AM).

**⇒ 9b's field-theoretic foundation is established:** the substrate supports a driven, frequency-tunable, bounded `(A, ω)` excess — the basis of heat-as-driven-excess and the entry point to the full 9b.1 program (the AM/FM 10.1a–e sub-experiments).

**Literature cross-ref (2026-06-05):** the IBM 2D discrete-time-crystal experiment (Switzer et al., arXiv:2501.18036, `theory/time_crystal_2d_dtc.pdf`) maps, in a driven quantum spin lattice, exactly the dial 9b cares about — **where the driven system ABSORBS drive energy (ergodic phase, heats toward thermalization) vs where absorption is BLOCKED (MBL/DTC coherence)**, controlled by kick angle + coupling + disorder. The quantum many-body analog of 9b.1's frequency-selective absorption: which knob gates a driven defect lattice between coherent (non-absorbing) and thermalizing (absorbing). Companion anchor: the **continuous** time crystal (Greilich et al., *Nat. Phys.* 20, 631 (2024), `theory/time_crystal_continuous.pdf`) **MELTS into chaotic auto-oscillation as temperature rises** (chaos parameter `K→1` at 17 K) — a defect-level picture of thermal decoherence: temperature converts a coherent intrinsic clock into aperiodic ticking, the coherent↔molten transition adjacent to 9b's `(A, ω)`-statistics view of heat. Physics cross-refs only; engineering implications → SABER DHC docs. **Honest caveats** (this is a *preview*, not 9b.1 proper): sandbox numpy (not production Taichi); the drive is an *idealized localized forcing pattern*, a proxy for a self-consistent incident EM wave (fuller treatment = 9b.1 / Phase 4); `A_core` is a shell-averaged excitation proxy, not yet the full 9b.0 thermodynamic temperature observable; frequency is in sandbox `/t`-units (physical Hz calibration = M5.8/M5.9); the closed Dirichlet box reflects energy, so the steady-state is a drive-vs-internal-dispersal balance (the right qualitative picture, but a true radiation sink is a refinement). Engineering-relevant framing (the lever + its resonant frequency) → SABER DHC docs.

#### Thermal-viz support — carried over from M5.6.5b (deferred 2026-05-27)

These rendering features were built-but-deferred during M5.6.5b (the EM/thermal/clock viz work, `4b Part 3`) because they only become *meaningful* once the 9b thermal program runs. The amplitude-`A` (`‖M−D‖_F`, WAVE_MENU 2) and clock-`ω` (`‖Ṁ‖_F`, WAVE_MENU 3) trackers are already wired; these three build on them. Pick them up here:

| Viz feature | What it shows | Supports |
| --- | --- | --- |
| **Joint-thermal `(A·ω)²`** color mode | the heat *energy* content in one view. Use `(A·ω)²`, **NOT** `A·ω` — for a defect oscillator `E_kin = ½m(Aω)²`, so `A·ω` alone has dimensions of *velocity* (peak speed = amplitude × frequency); the energy-dimensional quantity is the square. Product-of-squares of the two existing trackers. | 9b.0 (temperature observable made visible), 9b.1 |
| **Granule heat-map** | color each granule (point cloud at voxels) by local thermal `A` — a sparse per-defect "heat map" | 9b.3 (multi-defect amplitude statistics) |
| **Modulation-response** view | apply an EM-wave lever, render the `Δ(A·ω)²` shift before/after (or as a time-trace) — the SABER-method "modulation" picture made visible (physics framing only; device specifics stay in the SABER repo) | 9b.1 (AM/FM sub-experiments 10.1a–e), 9b.8 (wave back-reaction) |

Dynamic-charge note: the gauge-stable charge view (`|∇·n̂|` / topological winding, so the EM-charge view doesn't sign-flip under dynamics) is the *other* deferred M5.6.5b viz item — it's homed at **M5.7** (first sustained dynamic runs), not here.

Even if the framework needs refinement (first hypotheses rarely survive intact), the *framework for asking the question* — per-defect amplitude as a thermal degree of freedom alongside topology and mass — is what M5+ matter physics enables. **No other framework currently lets us pose the question this way.** That's the deeper value: even partial validation moves thermal physics from "ensemble-only statistical" to "single-defect quantum-mechanical".
