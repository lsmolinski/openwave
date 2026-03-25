# PHASE 2: BUILDING BLOCKS

## EWT Particle Hierarchy

| Scale | Structure | Example | K | Geometry |
| --- | --- | --- | --- | --- |
| **Fundamental particle** | 1 wave center | Neutrino (electron-type) | 1 | Single WC, symmetric |
| **Standalone particle** | K wave centers locked at standing wave nodes | Electron e⁻ / Positron e⁺ | 10 | Tetrahedral 1-3-6 |
| **Composite particle** | Multiple standalone particles at nodes | Proton (4e⁻ + 1e⁺), Neutron | — | Tetrahedral |
| **Atom** | Composite particles + orbital electrons | Hydrogen, Helium | — | Nuclear + shells |
| **Molecule** | Bonded atoms | H₂O | — | Molecular geometry |

Not all standalone particles are stable. Stability depends on geometry: K=1 (neutrino) and K=10 (electron) are the most stable. K=2 through K=9 are temporary — they decay because the geometry doesn't support all WCs sitting at standing wave nodes simultaneously. The magic K numbers (1, 8, 10, 20, 28, 50) correspond to tetrahedral configurations with enhanced stability.

### Dual vs Non-Dual Geometry — The Origin of Charge

Stable standalone particles split into two families based on tetrahedral symmetry:

**Neutrino family (DUAL tetrahedra)** — two interlocked tetrahedra of opposite phase (matter + antimatter nodes). The opposite-phase halves destructively interfere externally → no net wave amplitude → **neutral charge**. The two halves spin in opposite directions → **spins cancel** → no net transverse wave → no magnetic force.

| K | Particle | Structure | Charge |
| --- | --- | --- | --- |
| 1 | Electron neutrino | Single WC (trivially symmetric) | Neutral |
| 8 | Muon neutrino | Dual 2-level tetrahedra (2×4 WCs) | Neutral |
| 20 | Tau neutrino | Dual 3-level tetrahedra (2×10 WCs) | Neutral |

**Electron family (NON-DUAL, single tetrahedron)** — all WCs on the same phase node. No internal cancellation → **net wave amplitude** → charge. The asymmetric geometry forces spin → L→T conversion → **transverse wave = magnetic force**.

| K | Particle | Structure | Charge |
| --- | --- | --- | --- |
| 10 | Electron / Positron | 3-level tetrahedron (1-3-6) | Charged (e⁻/e⁺) |
| 28 | Muon | Combined dual 2-level + 3-level | Charged |
| 50 | Tau | Higher dual structure | Charged |

The K values follow **tetrahedral numbers**: 1-level = 1, 2-level = 4, 3-level = 10, 4-level = 20. Duals are sums: 4+4=8, 10+10=20, 4+4+10+10=28. These match the **nuclear magic numbers** (2, 8, 20, 28, 50) — the same geometric stability principle operates at both subatomic and nuclear scales.

**Key insight for simulation**: it's not just about K being large enough for spin — it's about the **non-dual asymmetry**. K=8 (muon neutrino) has 8 WCs but is neutral because it's dual. K=10 (electron) has 10 WCs and is charged because it's non-dual. The electron is the **simplest non-dual tetrahedral number** — this is why it's the lightest charged particle.

**Electron vs Positron**: same K=10, same 1-3-6 geometry, but WCs occupy opposite standing wave nodes (λ/2 offset). Same structure, opposite phase = opposite charge. Both are non-dual single tetrahedra.

References: [Subatomic Particles](https://energywavetheory.com/subatomic-particles/), [Particle Creation & Decay](https://energywavetheory.com/subatomic-particles/particle-creation-and-decay/), [Neutrino](https://energywavetheory.com/subatomic-particles/neutrino/), [Electron](https://energywavetheory.com/subatomic-particles/electron/), [Proton](https://energywavetheory.com/subatomic-particles/proton/)

---

## WHERE WE STAND

Phase 1 exhaustively tested whether Coulomb force emerges from two fundamental particles (single WCs at K=1 or simplified K=10). Key results:

| Effect | Status | Mechanism |
| --- | --- | --- |
| Strong force (lock-in) | ✅ Emerges | Sinc `cos(k·Δr+Δφ)` creates energy wells at λ/2 — same-phase WCs lock in |
| Annihilation | ✅ Emerges | Opposite-phase WCs: deepest well at r=0, barriers at λ/2 (positronium) |
| Charge discrimination | ✅ Exists | 2D flux: 100% same≠opposite at ALL separations (22/22) |
| Consistent Coulomb direction | ❌ Not from single WCs | Sinc flips force every λ/2 — intrinsic to monochromatic spherical interference |
| Correct Coulomb sign | ❌ Inverted | Isolated interaction: same→ATT, opposite→REP (180° off Coulomb) |

**Critical realization**: spherical 3D waves from single WCs (fundamental particles) produce lock-in and annihilation (near-field) but do NOT produce Coulomb alone. The sinc oscillation IS the correct near-field physics. Coulomb requires something that only standalone particles (K≥10) with spin can produce. A neutrino (K=1) is electrically neutral — as its name literally implies. We cannot expect Coulomb from a single WC.

Full Phase 1 results: [01z_results.md](01z_results.md)

---

## KEY INSIGHT: THE FORCE HIERARCHY

**Invert the path.** Phase 1 tried to go directly from single WCs to Coulomb. That was wrong. Forces emerge in a hierarchy — each building on the previous:

```text
1. STRONG FORCE (fundamental, from wave interference)
   └── sinc lock-in: WCs lock at λ/2 nodes
   └── annihilation: opposite-phase WCs cancel waves out at r=0

2. PARTICLE FORMATION (from strong force lock-in)
   └── fundamental particles (K=1, neutrino) lock into standalone particles
   └── neutrino (K=1) → electron (K=10, tetrahedral 1-3-6)
   └── standing waves reorganize into concentric spherical nodes
   └── not all K values are stable: K=2..9 decay, K=10 is the most stable

3. SPIN EMERGENCE (from NON-DUAL standalone particle geometry)
   └── only non-dual tetrahedra (electron family: K=10, 28, 50) have net spin
   └── dual tetrahedra (neutrino family: K=1, 8, 20) have opposite spins that cancel
   └── in non-dual: WCs can't all sit at nodes simultaneously
   └── one WC goes off-node → displaces the next → continuous rotation = spin
   └── spin converts longitudinal → transverse wave energy (L→T at rate α)
   └── spin creates a NEW transverse out-wave (not present in the in-wave)

4. CHARGE + ELECTRIC FORCE (from non-dual spin, from standalone particle)
   └── charge ≠ phase. Phase → lock-in/annihilation (near-field only)
   └── charge = net wave amplitude at the first wavelength of the standalone particle
   └── dual particles: opposite-phase halves cancel → zero net amplitude → neutral
   └── non-dual particles: all same-phase → net amplitude → charged
   └── single WC (neutrino) = neutral, no spin, no charge — cannot produce Coulomb
   └── spin promotes the disturbance needed for Coulomb to emerge at far-field:
       phase shift, toroidal wavefront (not pure spherical), L→T asymmetry
   └── traveling waves beyond particle radius = electric force

5. MAGNETIC FORCE (from transverse wave created by spin)
   └── L→T conversion creates a new transverse wave that IS the magnetic force
   └── accurately modeled as the Bohr magneton
   └── requires spin alignment/coherence for net magnetic field

6. GRAVITATIONAL FORCE (from the spin energy deficit)
   └── L→T drainage creates a longitudinal amplitude deficit
   └── deficit = gravitational shading. Ratio ≈ α accumulated over K WCs
```

**Phase ≠ Charge:**

- **Phase** (source_offset 0 vs π): determines matter vs antimatter. Same phase → strong force lock-in. Opposite phase → annihilation. Antiparticles (e.g. positron) occupy opposite standing wave nodes (λ/2 offset from electron). This is NEAR-FIELD behavior, already validated in Phase 1
- **Charge**: emerges from the standalone particle's **geometry** (dual vs non-dual). Dual tetrahedra (neutrino family) have opposite-phase halves that cancel externally → neutral. Non-dual tetrahedra (electron family) have all WCs same-phase → net wave amplitude → charged. The non-dual geometry also forces spin (WCs can't all sit at nodes) → L→T conversion → transverse wave + phase shift + toroidal wavefront → Coulomb. This is why Phase 1 could not produce Coulomb from single WCs — a single WC is the trivially symmetric (dual) case, inherently neutral

---

## FUNDAMENTAL CONCEPTS TO TEST

### 1. Standing Wave Node Reorganization

A single WC doesn't create new energy — it reorganizes the base wave's standing wave nodes into concentric spherical rings, radially oriented to the WC. The amplitude concentrates toward center, but λ also increases (Yee & Hauger shells), so energy steepness A/λ is conserved. Energy stays constant; only the node geometry changes.

This reorganization IS the lock-in mechanism: spherical concentric nodes create the energy wells that trap other WCs at node positions.

### 2. Energy Steepness Conservation

```text
A(r) increases toward center   (amplitude concentrates)
λ(r) increases toward center   (Yee & Hauger: shell n=1 = 18λ₀ for K=10)
A/λ ≈ constant                 (steepness preserved)
E = ρV(c·A/λ)² ≈ constant     (energy conserved despite amplitude variation)
```

The only difference between disturbed and undisturbed regions is the standing wave structure (nodes/antinodes), not the total energy.

### 3. Scalar vs Vector Energy at Near-Field

3D spherical wave superposition always produces elliptical granule motion. Even at near-field, the energy field has vector character (longitudinal + transverse projections). Questions:

- Does lock-in work with scalar energy only, or does the vector structure matter?
- Does the elliptical motion affect lock-in stability?
- How do we extract scalar energy from the elliptical trajectories?

### 4. Charge from Dual vs Non-Dual Geometry

Charge is determined by the **symmetry** of the standalone particle, not by the phase of individual WCs:

**Non-dual (electron family, K=10, 28, 50)**: all WCs on the same phase node → out-waves constructively interfere externally → net wave amplitude = charge. The asymmetric geometry can't accommodate all WCs at nodes → forced continuous spin → L→T conversion. This spin:

- Converts longitudinal → transverse energy (L→T at rate α ≈ 1/137)
- Creates a transverse out-wave NOT present in the in-wave (the magnetic force)
- Promotes a phase shift (LaFreniere: core λ/2 shift from medium compression)
- Modifies the wavefront shape (no longer pure spherical — toroidal/elliptical)
- Produces the far-field traveling wave disturbance that IS the electric field

**Dual (neutrino family, K=1, 8, 20)**: two interlocked tetrahedra of opposite phase → out-waves destructively interfere externally → no net wave amplitude = neutral. Spins of opposite halves cancel → no net transverse wave → no magnetic force. Even K=20 (20 WCs, far more than the electron's 10) is neutral because of this cancellation.

The electron (K=10) is the **simplest non-dual tetrahedral number** — the lightest charged particle. This is why Phase 1 could not produce Coulomb from single WCs: K=1 is the trivially symmetric case, inherently neutral.

### 5. Flux vs Gradient Force

Phase 1 showed:

- `F = -∇E` (gradient): charge-blind when base wave present, wrong sign for isolated interaction
- `S = -c²·ψ·∇ψ` (flux/radiation pressure): 100% charge discrimination in 2D

LaFreniere explicitly uses radiation pressure, not energy gradient, for Coulomb. Evaluate which force mechanism applies at each regime.

### 6. Granule Elliptical Motion and Vector Force Decomposition

When 3D spherical waves superpose, medium granules trace elliptical trajectories — not linear oscillation. This ellipse is fundamental: it encodes both electric and magnetic field information at every point in space.

The ellipse at each grid point can be decomposed by projecting its vector amplitude onto a radial direction from the particle experiencing the force:

- **Longitudinal projection** (along radial): the electric field component. Gradient of longitudinal energy → electric force
- **Transverse projection** (perpendicular to radial): the magnetic field component. Gradient of transverse energy → magnetic force
- These two projections have a **90° relationship** — the right-hand rule and the E×B force relationship emerge from the orthogonal decomposition of the same elliptical motion

Force is a vector field: `F = -∇E_L - ∇E_T`. The energy must be decomposed into vector components (longitudinal and transverse) to compute the correct force direction. Scalar energy `E = ρV(fA)²` collapses the directional information — it gives correct 1/r² magnitude but loses the charge-dependent direction. Vector energy `E = E_L + E_T` preserves it.

Gravitational force is the residual: total energy deficit from L→T spin conversion, radially toward mass, always attractive, omnidirectional.

---

## PHASE 2 SUMMARY

Task checklists are tracked in [00_ROADMAP.md](00_ROADMAP.md). Below is the strategy and context for each sub-phase.

### 🔶 Phase 2a: Near-Field — Particle Formation from Standing Waves (M3 Scalar)

**Goal**: Demonstrate stable standalone particle formation in M3 3D simulator. The near-field is purely longitudinal standing wave physics — M3 (scalar method) should be sufficient. We already simulate lock-in and annihilation; the missing piece is the K=10 electron tetrahedron stability.

**Why M3 is enough for near-field**: strong force / lock-in is a longitudinal phenomenon — sinc energy wells from spherical standing wave interference. No transverse component (spin) is needed for lock-in itself. M3 already demonstrates the right near-field physics. Far-field (Coulomb, magnetic) needs spin → needs vector field → that's Phase 2b on M4.

### 🚧 Phase 2b: Far-Field — Electromagnetism and Traveling Waves (M4 Vector)

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

- Wolff-original complex sinusoidal: `ψ = A · e^(iωt) · sin(kr)/r` — the imaginary component IS the transverse wave (see [00a_equations.md](00a_equations.md))
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

- [01c_vector_wave.md](01c_vector_wave.md) — Steps 1-2 validated: 3D base wave, L/T decomposition, WC out-wave with spin, sinc = correct K=1 physics
- Scripts: `scripts_phase1_vector/step1_base_wave.py`, `step2_single_wc.py`, `step3_two_wc_force.py`

### 🚧 Phase 2c: Gravitational Force and Composite Particles

**Goal**: Demonstrate gravity from spin energy deficit and test composite particle (proton = 4e⁻ + 1e⁺, neutron = proton + e⁻ at center) formation.

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
5. Is flux or gradient the correct force computation for each regime (near-field gradient vs far-field flux)?
6. Does the Yee & Hauger λ profile (longest near core) or LaFreniere profile (shortest near core) produce the correct physics? Or do they describe different regions (core vs shells)?
7. Can the sinc lock-in + charge discrimination from 2D flux be unified into a single force law that gives strong force at near-field and Coulomb at far-field?
8. Does the standalone particle's far-field traveling wave naturally produce Coulomb, or is the composite particle level (proton) needed?
