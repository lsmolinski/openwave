# OpenWave Models: the comparison table

## Why multiple models

OpenWave's mission is to build a platform where multiple candidate field-theoretic models are evaluated **numerically, side by side, in the same computational environment**. No single alternative framework can map the space of possibilities on its own: when several independent models are run against the same observables with the same pass/fail criteria, the comparison **triangulates** what is actually out there.

Features that survive across frameworks are likely load-bearing physics; features that only work in one framework, or only with hand-tuning, reveal themselves as such. A null result in one model is ambiguous (model wrong, or engine wrong?); a positive result in any model certifies the engine for all of them.

Anybody can contribute to building these numerical validations. Every cell in the table below is backed by a runnable script or a research document in this repository, and every claim is reproducible, refutable, and extendable under [Apache 2.0](LICENSE).

## Main use case: emergent science feeding new technology

OpenWave serves as a repository of **emergent science with concrete numerical validation status**, so that applied researchers can build on what holds and route around what does not.

> The goal of the OpenWave effort is to build working models on an open-source platform that anyone can use for applied research, with different models hosted in the same computational environment, supporting new tech development from emergent science.

## Validation status legend

| Icon | Meaning |
| --- | --- |
| ✅ | validated in-platform (runnable reproduction exists) |
| ⚠️ | partial, or validated with documented caveats |
| ❌ | tested and failed, or honest negative on record |
| 🔶 | in progress |
| 🚧 | planned, not yet tested in-platform |

A ❌ is a result, not an embarrassment: documented negatives (with the scripts that produced them) are part of the platform's value.

## COVERAGE MATRIX (Phenomenological Coverage)

Models:

- **M5: LC** (Liquid-Crystal topological defects, Jarek Duda, with Robert Close and Manfried Faber inputs);
- **M6: Ouroboros** (chaoiton framework, Paul Werbos);
- **M3: EWT** (Energy Wave Theory, Jeff Yee, built on Milo Wolff and Gabriel LaFreniere pioneer work).

**Column order:** models are sequenced by their validated + partial count (✅ + ⚠️ in the [Summary count](#summary-count) below): M5 (12), then M6 (9), then M3 (8). The order updates as validations land.

Script paths are relative to `openwave/xperiments/`. Rows are grouped by domain: particles, forces, waves + quantum emergence.

### PARTICLES (full spectrum)

| Criteria | Liquid Crystal (M5) | Ouroboros (M6) | EWT (M3) |
| --- | --- | --- | --- |
| Charge quantization | ✅<br>Topological winding number of the hedgehog defect (Gauss-Bonnet integer Q = ±1)<br>`m5_liquid_crystal/research/sandbox_v2/m5_1_winding.py` | ⚠️<br>Mutual Chern-Simons linking number of A and J flux lines; elementary charge within 0.6% (author's claim + Lean 4 artifacts, not yet re-derived in-platform)<br>`m6_ouroboros/research/0d_canonical.md` | ❌<br>Charge sign imposed via `cos(source_offset)`, not emergent from wave physics<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Electron rest energy (mass) | ✅<br>Hedgehog rest energy with Faber core regularization; mass pinned E ∝ 1/r₀, physical knob r₀ = 2.2132 fm → 0.511 MeV<br>`m5_liquid_crystal/research/sandbox_v6/m5_6_3b_faber_on_M.py` | ✅<br>Electron calibration H/Q = 1.6969 reproduced to 0.090% at g = 1.0 (canonical benchmark)<br>`m6_ouroboros/research/sandbox_v8/ouroboros_benchmark.py` | ⚠️<br>Wave-center standing-wave lock-in demonstrated; mass values come from EWT's analytic equations, not yet from in-sim dynamics<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Lepton mass spectrum (μ, τ) | 🚧<br>M5.9 target: three axis choices of the biaxial hierarchy 0 < δ ≪ 1 ≪ g as the three lepton families<br>`m5_liquid_crystal/research/0b_M5_roadmap.md` | ⚠️<br>Muon 0.80%, tau 6.47% at the chosen harmonic indices; caveat: our scan found no discrete-spectrum mechanism selecting those ω values (every ω in 1-60 equally localized)<br>`m6_ouroboros/research/sandbox_v8/ouroboros_benchmark.py` | ❌<br>K-selectivity not achieved: all K = 2..10 equally stable at perfect placement, K = 10 breaks worst under perturbation<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Neutrinos | 🚧<br>Closed vortex loop, a distinct topology class from the hedgehog; flavour oscillation as transitions between loop states; hopfions as excited-oscillation candidates (lab-existence anchor: Liu et al., Nature Physics 2026)<br>`m5_liquid_crystal/research/0b_question_tracker.md § PARTICLES` | 🚧<br>Not addressed; where active neutrinos fit is an open question in the framework<br>(none yet) | ⚠️<br>The neutrino is EWT's fundamental wave-center unit (postulated as the substrate; magic-number K = 2, 8, 10 stability is the open target)<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Quarks | 🚧<br>String endpoint of a 1D topological vortex (not a 0D hedgehog); Cornell confinement, M5.9 target<br>`m5_liquid_crystal/research/0b_M5_roadmap.md` | 🚧<br>Not directly addressed; a 3-chaoiton proton (Schwinger H-particle) is implicit in the dyon framing, not computed<br>(none yet) | 🚧<br>Not modeled in-sim<br>(none yet) |
| Antimatter + annihilation | ⚠️<br>Q → −Q under reflection; annihilation = topology-cancellation releasing stored field energy as outgoing waves (structural, not yet quantified; two-defect ± demo planned)<br>`m5_liquid_crystal/research/0b_question_tracker.md § PARTICLES` | 🚧<br>Q_CS = −1 positron analog identified, not yet computed numerically<br>(none yet) | ⚠️<br>Opposite-phase wave centers annihilate in-sim, with documented assists (0.5λ threshold, damping, velocity clamp)<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Magnetic moment μ + spin J | 🔶<br>ELECTRON-ID project in progress: 3×3 fixed-clock hedgehog, integrate μ and J, unit-free g-factor ≈ 2 target<br>`m5_liquid_crystal/research/0b_M5_roadmap.md § ELECTRON-ID` | 🚧<br>Spin from chaoiton field rotation is paper-level, not yet in-platform<br>(none yet) | 🚧<br>Bohr magneton listed as target, not yet attempted<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| de Broglie clock (Zitterbewegung) | ✅<br>Bounded, self-starting, frequency-rigid 3+1D time crystal: quadratic action fails, signed quartic saturates, ω₁ is a start-independent attractor; classified a molten clock that regularizes toward a cold ground state<br>`m5_liquid_crystal/research/sandbox_vn/m5_8_2h_omega_attractor.py`, `m5_8_2j_zbw_ratio.py` | ⚠️<br>Time-periodicity is built into the ansatz (e^{iωt}) rather than emergent; the L/Q = ω Noether identity reproduces exactly across the full scan<br>`m6_ouroboros/research/sandbox_v8/ouroboros_benchmark.py` | 🚧<br>EWT carries particle frequency but no in-sim clock-propulsion mechanism<br>(none yet) |
| Particle stability (Derrick escape) | ✅<br>Escape via time-periodic resonance: static solitons confirmed impossible (M5.2 negative), the saturated breather self-starts from exact rest and holds resolution-robustly (24³ → 48³)<br>`m5_liquid_crystal/research/sandbox_v8/m5_8_2g_spontaneity.py`, `sandbox_vn/m5_8_2i_dispersal_gate.py` | ✅<br>Escape via oscillation: true neutral ground state found by BVP (zero sign changes, K₁ exponential tail)<br>`m6_ouroboros/research/sandbox_v11/` | ⚠️<br>Standing-wave lock-in holds at perfect placement, fragile under perturbation; annihilation requires threshold + damping assists<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Dark matter candidate | 🚧<br>None claimed (hopfions are a forward-looking candidate, post-M5.8)<br>(none yet) | ✅<br>Neutral chaoiton: m_χ = 0.460 MeV with mediator m_J = 0.6184 MeV parameter-free via the exact scaling symmetry; canonical β(r) profile + dipole form factor independently computed in-platform<br>`m6_ouroboros/research/sandbox_v11/`, `sandbox_v11/dm_paper_supplement/` | 🚧<br>None claimed<br>(none yet) |
| Spin-½ statistics (720° belt trick) | ⚠️<br>The biaxial fundamental group is the quaternion group π₁ = Q₈, which supports spinor-like 720° behavior (a uniaxial director cannot); structural, the charge-class demo is pending<br>`m5_liquid_crystal/research/0b_M5_roadmap.md § NG-9` | 🚧<br>Not addressed<br>(none yet) | 🚧<br>Scalar field carries no spinor structure<br>(none yet) |
| Baryons (p, n) | 🚧<br>Nucleons as 3-string Y-configurations of quark endpoints (proton encloses charge into the hedgehog, the lighter baryon; the neutron must compensate the charge, the heavier); deferred to the 9d stage<br>`m5_liquid_crystal/research/9d_composite_particles.md` | 🚧<br>The 3-chaoiton proton (Schwinger H-particle) and the ≈ 0.84 fm proton radius remain author claims, not yet computed<br>(none yet) | ⚠️<br>K = 10 tetrahedron holds at perfect placement using the Combined Wolff-LaFreniere equation, breaks under perturbation<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Mesons (π, K) | 🚧<br>Pion as twist/reconnection of a vortex loop; kaon as a knotted vortex loop (strangeness = knotting); paper-level, lands with the 9d stage<br>`m5_liquid_crystal/research/9d_composite_particles.md` | ⚠️<br>Pion+ candidate at ω = 15.0 reproduces 139.57 MeV to 3.25% in the canonical benchmark<br>`m6_ouroboros/research/sandbox_v8/ouroboros_benchmark.py` | 🚧<br>Not modeled in-sim<br>(none yet) |

### FORCES

| Criteria | Liquid Crystal (M5) | Ouroboros (M6) | EWT (M3) |
| --- | --- | --- | --- |
| Electric force (Coulomb 1/r) | ✅<br>E(d) ~ 1/d between two hedgehogs from pure topology, R² = 0.978; reproduced on the production matrix field (R² = 0.97) + analytical page-18 cross-validation (R² = 0.996)<br>`m5_liquid_crystal/research/sandbox_v2/m5_1_coulomb.py`, `sandbox_v4/m5_4_coulomb_matrix.py` | 🚧<br>Static two-charge derivation exists at paper level; force-level Coulomb between chaoitons not yet tested in-platform<br>(none yet) | ❌<br>Sinc envelope barriers block far-field attraction/repulsion; signed envelope is a modeling choice<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Magnetic force | ⚠️<br>Per-defect magnetism = S¹-loop topology (vs S² for electric) carried by the curl/transverse component of the defect's outgoing wave; coherent macroscopic regimes (magnets, electromagnets) structural; quantitative observability is a pending falsifiable target<br>`m5_liquid_crystal/research/0b_question_tracker.md § FORCES` | ⚠️<br>Contained in the A_μ Maxwell sector by construction; no per-defect magnetic structure computed<br>`m6_ouroboros/research/0d_canonical.md` | 🚧<br>Scalar model carries no polarization structure to support magnetism<br>(none yet) |
| Strong force / confinement | ⚠️<br>Short-range mechanism verified: running-coupling onset at r₀ (non-abelian ‖R‖·r² roll-off, Maxwell as the abelian long-range limit); the linear Cornell confinement V(r) = −α/r + σr via 1D vortex string is the M5.9 target<br>`m5_liquid_crystal/research/sandbox_v6/m5_6_4b_faber_curvature_em.py` | 🚧<br>Sawada long-range nuclear anomaly v(r) ~ −C/r⁶ identified as falsifiability target, not yet tested<br>(none yet) | 🚧<br>Listed as end-game target<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Weak force | ⚠️<br>Sketched: beta decay as a topology-reconnection event (defect-class transition n → p + e + ν), with the W boson as a tilt-sector excitation of the connection; no clean SU(2) chiral mechanism yet (open question)<br>`m5_liquid_crystal/research/0b_question_tracker.md § FORCES` | 🚧<br>Not addressed<br>(none yet) | 🚧<br>Not addressed<br>(none yet) |
| Gravity | 🔶<br>The 4×4 time axis is the gravity slot: gravity enters only via the boost tilt (GEM ∝ (b·g)², exactly zero at zero boost, negative = the clock-fuel block); dynamical metric not implemented<br>`m5_liquid_crystal/research/sandbox_vn/m5_8_2q_delta_scaling.py` | 🚧<br>Not in the Lagrangian (the framework explicitly stops before gravity)<br>(none yet) | 🚧<br>Not modeled<br>(none yet) |

### WAVES + QUANTUM EMERGENCE

| Criteria | Liquid Crystal (M5) | Ouroboros (M6) | EWT (M3) |
| --- | --- | --- | --- |
| EM waves (Maxwell) | ✅<br>Maxwell recovered by two independent routes: the hydrodynamic dictionary (abelian) and Faber's curvature R = Γ×Γ; tilt modes propagate at c, with the divergence/curl (electric/magnetic) decomposition of each defect's outgoing wave<br>`m5_liquid_crystal/research/sandbox_v6/m5_6_4a_hydro_em.py`, `m5_6_4b_faber_curvature_em.py` | ✅<br>A_μ is the electromagnetic four-potential by construction; delocalized J-field wave modes coexist with solitons<br>`m6_ouroboros/research/0d_canonical.md` | ⚠️<br>Scalar wave propagation only (no polarization structure)<br>`m3_wolff_lafreniere/research/` |
| Quantum wave equation (Klein-Gordon) | ✅<br>Klein-Gordon emerges from the biaxial twist with GEOMETRIC mass (minimal coupling to the hedgehog connection; the explicit mass term cancels, core regularization generates it)<br>`m5_liquid_crystal/research/sandbox_v6/m5_6_1_kg_operator_check.py`, `m5_6_1b_twist_evolution.py` | 🚧<br>QM not derived; the classical field carries the e^{iωt} ansatz, quantum behavior is outside current scope<br>(none yet) | ⚠️<br>The scalar wave equation is the postulated substrate, not an emergent result<br>`m3_wolff_lafreniere/research/0_WAVE_EQUATION.md` |
| Orbital quantization (atomic structure) | 🚧<br>EM Coulomb + de Broglie standing-wave quantization (not a separate force); cross-mass-class binding deferred to 9d<br>`m5_liquid_crystal/research/9d_composite_particles.md` | 🚧<br>Not addressed<br>(none yet) | ⚠️<br>Standing-wave lock-in demonstrated: same-phase wave centers sit in energy wells at λ separation; selectivity fragile under perturbation<br>`m3_wolff_lafreniere/research/0_STATUS.md` |

### Summary count

| Status | Liquid Crystal (M5) | Ouroboros (M6) | EWT (M3) |
| --- | --- | --- | --- |
| ✅ validated in-platform | 7 | 4 | 0 |
| ⚠️ partial / with caveats | 5 | 5 | 8 |
| ❌ honest negative | 0 | 0 | 3 |
| 🔶 in progress | 2 | 0 | 0 |
| 🚧 planned / not yet | 7 | 12 | 10 |
| **Total criteria** | **21** | **21** | **21** |

The counts tell the platform story plainly: M5 has the widest validated coverage and the widest remaining program, M6 has fewer but sharper validated claims concentrated where its design aims, and M3 carries the honest negatives that motivated the later models (its documented failures are what pointed at topology and time-periodicity).

## Reading the table

The three frameworks escape Derrick's theorem three different ways (standing-wave interference, topology + time-periodic resonance, oscillation), and the table makes the triangulation visible: particle stability requires time-periodicity in every framework that achieves it, charge quantization only emerges where there is topology, and lepton mass spectra remain the open problem in all three. That convergence-and-divergence pattern is the platform's scientific product.

## MODEL PROFILE (decision-relevant attributes)

The coverage matrix scores phenomena; this companion table scores the model-level attributes a reader needs to weigh the columns: parameter economy, what formal artifacts back the claims, and what would falsify each model next.

| Attribute | Liquid Crystal (M5) | Ouroboros (M6) | EWT (M3) |
| --- | --- | --- | --- |
| Free parameters | δ (quantum phase), g (gravity/time axis), plus 1-2 potential (LdG) coefficients; the boost dressing b enters the clock sector. Calibration handles: Faber r₀ (mass), Coulomb units<br>`m5_liquid_crystal/research/4c_convo_2026.06.08.md` | 3 claimed (g, λ, ω); the neutral sector's exact scaling symmetry closes the (g, λ) plane, making m_J parameter-free (in-platform result)<br>`m6_ouroboros/research/0c_sandbox_v11.md` | EWT's analytic wave constants (amplitude, wavelength, density); in-sim runs add documented envelope/threshold choices per script<br>`m3_wolff_lafreniere/research/0a_equations.md` |
| Formal artifacts | Every claim backed by a runnable open script + research note; documented negatives (M5.2, the M5.7 nulls) preserved as results<br>`m5_liquid_crystal/research/10_summary_report.md § Reproduction` | Author's Lean 4 proof artifacts (linking number, mountain-pass existence, power counting) + our independent numerical reproduction of the canonical profile and benchmark<br>`m6_ouroboros/research/0d_canonical.md` | Runnable scripts + an explicit honest-blockers status doc<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Falsifiable near-term tests | Unit-free g-factor ≈ 2 from the fixed-clock electron (ELECTRON-ID, in progress); the absolute-ω calibration chain (Coulomb units + LdG-to-rest-energy)<br>`m5_liquid_crystal/research/0b_M5_roadmap.md § ELECTRON-ID` | Author's roadmap: NEGF vertex check, sub-MeV searches, six-peak Gaia-stream annual modulation of the J-field flux<br>`m6_ouroboros/research/0e_dm_paper_review.md` | K-selectivity under perturbation (currently failing, and documented as such)<br>`m3_wolff_lafreniere/research/0_STATUS.md` |
| Direct-detection compatibility (DM) | n/a (no DM candidate claimed) | Dipole-suppressed chaoiton-proton cross section claimed compatible with direct-detection bounds (the monopole coupling vanishes by angular-momentum orthogonality; the numerical chain is still being reconciled on the author's side)<br>`m6_ouroboros/research/0e_dm_paper_review.md`, `sandbox_v11/dm_paper_supplement/` | n/a |

## Per-model results of record

| Model | Deep dive |
| --- | --- |
| Liquid Crystal (M5) | [`m5_liquid_crystal/research/10_summary_report.md`](openwave/xperiments/m5_liquid_crystal/research/10_summary_report.md): the results-of-record; [`0b_M5_roadmap.md`](openwave/xperiments/m5_liquid_crystal/research/0b_M5_roadmap.md): full program; [`0b_question_tracker.md`](openwave/xperiments/m5_liquid_crystal/research/0b_question_tracker.md): emergence catalog + open questions |
| Ouroboros (M6) | [`m6_ouroboros/research/0d_canonical.md`](openwave/xperiments/m6_ouroboros/research/0d_canonical.md): canonical numerical specification |
| EWT (M3) | [`m3_wolff_lafreniere/research/0_STATUS.md`](openwave/xperiments/m3_wolff_lafreniere/research/0_STATUS.md): targets, achieved, honest blockers |

## Contributing a model or a validation

A new framework enters as a new `openwave/xperiments/<model>/` directory with its own `research/` folder; a new validation enters as a runnable script plus a research note documenting pass/fail against the shared criteria above. Negatives are as publishable here as positives. See [`CONTRIBUTING.md`](CONTRIBUTING.md) and [`SYS_ARCH.md`](SYS_ARCH.md) for structure, and open an issue to propose a new column.
