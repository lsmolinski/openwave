# PHASE 2: ENERGY LAYERS

The ENERGY LAYERS hierarchy, EWT particle terminology, dual/non-dual geometry, and concepts to investigate are documented in [0_OVERVIEW.md](0_OVERVIEW.md). Task checklists are tracked in [0_ROADMAP.md](0_ROADMAP.md).

This document contains Phase 2 strategy, context sources, and open questions.

> **🚧 TODO — Revision pending (from Phase 3 concept review)**: insert a **Layer 0: Vacuum at rest** before Layer 1, and refine Layer 1-3 semantics to align with the Lagrangian / topological framework. The current hierarchy assumes a pre-existing "base wave" (EWT view). The Lagrangian view treats the medium as a *static* ordered vacuum (ground state of V(ψ)) with waves arising as perturbations. See [3b_concept_review.md](../../m5_lagrangian_field/research/3b_concept_review.md) for the full rationale and proposed updated hierarchy table. Scope: update this doc + [0_OVERVIEW.md](0_OVERVIEW.md) + [README.md](../../../../README.md) once the M3 near-field work completes and the M5 / LAGRANGIAN-FIELD METHOD is under construction.

---

## WHERE WE STAND

Phase 1 exhaustively tested whether Coulomb force emerges from two fundamental particles (single WCs at K=1 or simplified K=10). Key results:

| Effect | Status | Mechanism |
| --- | --- | --- |
| Standing wave lock-in | ✅ Emerges | Sinc `cos(k·Δr+Δφ)` creates energy wells at λ/2 — same-phase WCs lock in |
| Annihilation | ✅ Emerges | Opposite-phase WCs: deepest well at r=0, barriers at λ/2 (positronium) |
| Charge discrimination | ✅ Exists | 2D flux: 100% same≠opposite at ALL separations (22/22) |
| Consistent Coulomb direction | ❌ Not from single WCs | Sinc flips force every λ/2 — intrinsic to monochromatic spherical interference |
| Correct Coulomb sign | ❌ Inverted | Isolated interaction: same→ATT, opposite→REP (180° off Coulomb) |

**Critical realization**: spherical 3D waves from single WCs (fundamental particles) produce lock-in and annihilation (near-field) but do NOT produce Coulomb alone. The sinc oscillation IS the correct near-field physics. Coulomb requires standalone particles (K≥10) with non-dual geometry + spin. A neutrino (K=1) is electrically neutral — we cannot expect Coulomb from a single WC.

Full Phase 1 results: [1z_results.md](1z_results.md)

---

## PHASE 2 STRATEGY

### 🔶 LAYERS 1-3: Particle Formation from Standing Waves (M3 Scalar, near-field)

**Goal**: Demonstrate stable standalone particle formation in M3 3D simulator. The near-field is purely longitudinal standing wave physics — M3 (scalar method) should be sufficient. We already simulate lock-in and annihilation; the missing piece is the K=10 electron tetrahedron stability.

**Why M3 is enough for near-field**: lock-in is a longitudinal phenomenon — sinc energy wells from spherical standing wave interference. No transverse component (spin) is needed for lock-in itself. M3 already demonstrates the right near-field physics. Far-field (Coulomb, magnetic) needs spin → needs vector field → that's Phase 2b on M4.

### 🚧 LAYER 4: Electro-Magnetism and Traveling Waves (M4 Vector, far-field)

**Goal**: Demonstrate how spin creates charge, Coulomb force, and magnetic force. This requires M4 (vector wave engine) because far-field forces emerge from the transverse wave component that spin creates — a scalar field cannot represent this.

**Why spin is the key**: Phase 1 showed single WCs (neutrinos) cannot produce Coulomb. The neutrino has no spin, no charge. The electron (K=10) has spin because its non-dual tetrahedral geometry forces continuous WC rotation. Spin creates L→T conversion, which:

- Generates a transverse out-wave (magnetic force)
- Promotes a phase shift from the spinning geometry (charge sign)
- Distorts the wavefront from pure spherical to toroidal/elliptical
- Creates the traveling wave pattern beyond particle radius that IS the electric field

**Context sources for spin development** (load these when starting Phase 2b):

EWT theory and equations:

- [EWT Forces: Electric](https://energywavetheory.com/forces/electric/) — charge as traveling wave amplitude
- [EWT Forces: Magnetic](https://energywavetheory.com/forces/magnetic/) — transverse wave from spin
- [EWT Forces: Unification](https://energywavetheory.com/forces/unification-of-forces/) — force hierarchy
- [EWT Electron](https://energywavetheory.com/subatomic-particles/electron/) — spin mechanism, L→T
- [EWT Constants: Fine Structure](https://energywavetheory.com/physics-constants/fine-structure-constant/) — α = L→T conversion ratio

Wave equation foundations:

- Wolff-original complex sinusoidal: `ψ = A · e^(iωt) · sin(kr)/r` — the imaginary component IS the transverse wave (see [0a_equations.md](0a_equations.md))
- LaFreniere radiation pressure model: `lafreniere/Gabriel_LaFreniere/sa_phaseshift.html` — core phase shift, diffractive lens, wave speed in compressed medium
- [LaFreniere Electron](http://www.interwaves.com/en/sa_electron.htm) — standing → traveling transition, near-field → far-field

Scientific papers (in `/scientific_source/`):

- `01. The Geometry of Spacetime and the Unification of Forces v2.3.pdf` — force unification framework
- `02. The Geometry of Particles and the Explanation of Their Creation and Decay v2.pdf` — dual/non-dual tetrahedra, particle geometry, spin
- `03. The Physics of SubAtomic Particles.pdf` — particle physics context
- `06. Constants and Equations - Waves.pdf` — primary constants reference
- Yee & Hauger wavelength shells: `references/Spin.pdf` — variable λ(r), standing → traveling transition

Additional references (in `/references/`):

- Smolinski: `references/MagnetismGravity_v4.pdf` — toroidal soliton geometry, r⁵ energy scaling, Energy Domain vs EMC Domain, non-linear wave equation `(∂²/∂t² - c²∇²)Ψ + k·Ψ³ = 0`
- Butto (2021): `references/` — vortex electron model, superfluid irrotational vortex, spin-½ as differential rotation (720° core / 360° boundary), Helmholtz theorems → toroidal geometry

Phase 1c vector wave research (validated infrastructure):

- [1c_vector_wave.md](1c_vector_wave.md) — Steps 1-2 validated: 3D base wave, L/T decomposition, WC out-wave with spin, sinc = correct K=1 physics
- Scripts: `sandbox_phase1_vector/step1_base_wave.py`, `step2_single_wc.py`, `step3_two_wc_force.py`

### 🚧 LAYER 5: Gravity Emergence

- Demonstrate gravity from spin energy deficit

### 🚧 LAYER 6: Emergent Waves

- **Goal**: Demonstrate Electro_magnetic Waves being generated and Thermal Energy.

### 🚧 LAYER 7: Composite Particles

- **Goal**: Demonstrate gravity from spin energy deficit and test composite particle (proton = 4e⁻ + 1e⁺, neutron = proton + e⁻ at center) formation.

---

## TOOLS AND TECHNIQUES TO EVALUATE

- Granule elliptical motion — vector amplitude decomposition (6 phasor numbers)
- Energy with λ(r) — variable wavelength energy equation `E = ρV(c·A/λ(r))²`
- Density variation ρ(r) — Smolinski density function, push-out operator
- Time dynamics — variable λ per voxel, local dt, time dilation
- Flux-based force computation — radiation pressure as primary Coulomb mechanism
- M4 vector wave engine — transverse displacement for magnetic force
- Taichi GPU acceleration — for 3D flux and large-grid simulations

---

## OPEN QUESTIONS

1. Does the electron tetrahedron stabilize on M3 (scalar) with variable λ(r), or does it require M4 vector forces?
2. Is the base wave required for lock-in, or do WCs alone create the standing wave structure?
3. Why is K=10 (1-3-6 tetrahedron) the most stable geometry? What makes K=2..9 unstable?
4. Can we simulate a dual tetrahedron (K=8 or K=20) and verify it produces zero net external amplitude (neutral)?
5. How exactly does non-dual spin promote the phase shift / wavefront distortion that creates charge? Is it the L→T conversion itself, the toroidal geometry, or the core compression?
6. Is flux or gradient the correct force computation for each regime (near-field gradient vs far-field flux)?
7. Does the Yee & Hauger λ profile (longest near core) or LaFreniere profile (shortest near core) produce the correct physics? Or do they describe different regions (core vs shells)?
8. Can the sinc lock-in + charge discrimination from 2D flux be unified into a single force law that gives strong force at near-field and Coulomb at far-field?
9. Does the standalone particle's far-field traveling wave naturally produce Coulomb, or is the composite particle level (proton) needed?
