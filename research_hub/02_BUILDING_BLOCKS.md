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

3. SPIN EMERGENCE (from standalone particle geometry)
   └── WCs can't all sit at nodes simultaneously in the tetrahedron
   └── one WC goes off-node → displaces the next → continuous rotation = spin
   └── spin converts longitudinal → transverse wave energy (L→T at rate α)
   └── spin creates a NEW transverse out-wave (not present in the in-wave)

4. CHARGE + ELECTRIC FORCE (from spin, from standalone particle)
   └── charge ≠ phase. Phase → lock-in/annihilation (near-field only)
   └── charge = wave amplitude at the first wavelength of the standalone particle
   └── spin promotes the disturbance needed for Coulomb to emerge at far-field:
       phase shift, toroidal wavefront (not pure spherical), L→T asymmetry
   └── single WC (neutrino) = neutral, no spin, no charge — cannot produce Coulomb
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

- **Phase** (source_offset 0 vs π): determines matter vs antimatter. Same phase → strong force lock-in. Opposite phase → annihilation. Antiparticles occupy opposite standing wave nodes (λ/2 offset). This is the NEAR-FIELD behavior, already validated in Phase 1
- **Charge**: emerges from the standalone particle's spin. The electron's 1-3-6 tetrahedral arrangement cannot have all WCs at nodes — the resulting spin creates L→T conversion, phase shift, and toroidal wavefront disturbance that produces far-field Coulomb. Charge = wave amplitude at the particle's first wavelength. Neutrino (K=1) has no spin, no charge — this is why Phase 1 could not produce Coulomb from single WCs

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

### 4. Charge from Standalone Particle Spin

The electron (K=10) has 10 same-phase WCs in a 1-3-6 tetrahedral arrangement. Not all WCs can sit at nodes simultaneously — the off-node WC displaces others in sequence, creating continuous spin. This spin:

- Converts longitudinal → transverse energy (L→T at rate α ≈ 1/137)
- Creates a transverse out-wave NOT present in the in-wave (the magnetic force)
- Promotes a phase shift (LaFreniere: core λ/2 shift from medium compression)
- Modifies the wavefront shape (no longer pure spherical — toroidal/elliptical)
- Produces the far-field traveling wave disturbance that IS the electric field

The fundamental particle (K=1, neutrino) = neutral, no spin, no charge. The standalone particle (K=10, electron) = charged, spinning, creates Coulomb. This is why we could not simulate Coulomb with single WCs in Phase 1 — a neutrino is literally neutral.

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

## ATTACK PLAN

### Phase 2a: Near-Field — Particle Formation from Standing Waves (M3 Scalar)

**Goal**: Demonstrate stable standalone particle formation in M3 3D simulator. The near-field is purely longitudinal standing wave physics — M3 (scalar method) should be sufficient. We already simulate lock-in and annihilation; the missing piece is the K=10 electron tetrahedron stability.

**Why M3 is enough for near-field**: strong force / lock-in is a longitudinal phenomenon — sinc energy wells from spherical standing wave interference. No transverse component (spin) is needed for lock-in itself. M3 already demonstrates the right near-field physics. Far-field (Coulomb, magnetic) needs spin → needs vector field → that's Phase 2b on M4.

- [ ] Stabilize K=10 electron tetrahedron on M3
  - Implement variable λ(r) from Yee & Hauger shells (non-uniform node spacing)
  - Test whether variable nodes accommodate the 15/45 non-node pairs (√3×λ/2, √2×λ/2)
  - Leapfrog integrator already in place, damping at 0.995
- [ ] Validate same-phase lock-in in 3D animation
  - Two fundamental particles (WCs), same phase, observe oscillation in energy wells
  - Measure well depth, oscillation period, escape energy
- [ ] Validate opposite-phase annihilation pathway
  - Two fundamental particles, opposite phase, observe attraction through barriers
  - Measure barrier heights at λ/2 intervals, compare with positronium lifetime data
- [ ] Test multi-WC lock-in: K=2, K=3, ..., K=10 progressive build-up
  - K=2 through K=9 should be unstable (temporary particles that decay)
  - K=10 should be the first fully stable standalone particle (electron)
  - Does the geometry self-organize or must it be prescribed?
  - Validates the EWT prediction: K=10 is special because the 1-3-6 tetrahedron is the simplest 3D geometry where all WCs can sit near nodes
- [ ] Characterize near-field → far-field transition boundary
  - Where do standing waves end and traveling waves begin?
  - How does this relate to standalone particle radius K²λ?
- [ ] Wrap M3 near-field validation
  - M3 has done its job once particle formation is demonstrated
  - No far-field analysis on M3 — Coulomb needs spin, which needs vector field (M4)
  - Rename xperiments

### Phase 2b: Far-Field — Electromagnetism and Traveling Waves (M4 Vector)

**Goal**: Demonstrate how spin creates charge, Coulomb force, and magnetic force. This requires M4 (vector wave engine) because far-field forces emerge from the transverse wave component that spin creates — a scalar field cannot represent this.

**Why spin is the key**: Phase 1 showed single WCs (neutrinos) cannot produce Coulomb. The neutrino has no spin, no charge. The electron (K=10) has spin because its tetrahedral geometry forces continuous WC rotation. Spin creates L→T conversion, which:

- Generates a transverse out-wave (magnetic force)
- Promotes a phase shift from the spinning geometry (charge sign)
- Distorts the wavefront from pure spherical to toroidal/elliptical
- Creates the traveling wave pattern beyond particle radius that IS the electric field

- [ ] Build / extend M4 vector wave engine with L→T spin conversion
  - Vector displacement: `ψ = (ψ_x, ψ_y, ψ_z)` per voxel
  - L→T conversion at standalone particle WCs (η = α ≈ 1/137)
  - Elliptical displacement trajectories (6 phasor numbers)
- [ ] Test standalone particle (K=10 electron) as Coulomb source on M4
  - Place one stabilized electron, measure far-field force on a test particle
  - Does the spinning tetrahedral geometry produce consistent Coulomb direction?
  - Compare with fundamental particle / neutrino (should be neutral — no Coulomb)
- [ ] Investigate 3D flux-based force
  - Extend Phase 1 2D flux test to full 3D spherical integration
  - Does solid-angle averaging smooth the sinc absolute direction?
  - Test at K=10 standalone particle scale (100λ radius)
- [ ] Investigate charge mechanism
  - How does spin + 1-3-6 geometry create the wavefront disturbance?
  - Role of LaFreniere core phase shift (λ/2 from compression)
  - Does variable λ(r) inside the standalone particle create the charge sign?
- [ ] Validate Coulomb
  - Direction: same charge repels, opposite attracts at ALL separations
  - Magnitude: 1/r² scaling
  - Symmetry: Newton's 3rd law (equal and opposite)
- [ ] Separate longitudinal (electric) and transverse (magnetic) force components
  - Electric: ∇E_L from traveling longitudinal waves beyond particle radius
  - Magnetic: ∇E_T from transverse wave generated by spin

### Phase 2c: Gravitational Force and Composite Particles

**Goal**: Demonstrate gravity from spin energy deficit and test composite particle (proton) formation.

- [ ] Gravitational force from spin deficit
  - L→T drainage creates longitudinal amplitude deficit
  - Test whether deficit accumulates over K WCs (factor of K×α?)
  - Compare with 10⁻⁴² EM-to-gravitational ratio
  - Validate against Smolinski's Scilab reference values
- [ ] Wave shading with standalone particle clusters
  - Multiple standalone particles → accumulated energy deficit → gravity
- [ ] Composite particle formation
  - Proton: 4 electrons + 1 positron at center (pentaquark, tetrahedral)
  - Neutron: proton + electron at center (neutralizes charge via destructive interference)
  - Test whether composite particles self-assemble from standalone particles

---

## TOOLS AND TECHNIQUES TO EVALUATE

- [ ] Granule elliptical motion — vector amplitude decomposition (6 phasor numbers)
- [ ] Energy with λ(r) — variable wavelength energy equation `E = ρV(c·A/λ(r))²`
- [ ] Density variation ρ(r) — Smolinski density function, push-out operator
- [ ] Time dynamics — variable λ per voxel, local dt, time dilation
- [ ] Flux-based force computation — radiation pressure as primary Coulomb mechanism
- [ ] M4 vector wave engine — transverse displacement for magnetic force
- [ ] Taichi GPU acceleration — for 3D flux and large-grid simulations

---

## OPEN QUESTIONS

1. Does the electron tetrahedron stabilize on M3 (scalar) with variable λ(r), or does it require M4 vector forces?
2. Is the base wave required for lock-in, or do WCs alone create the standing wave structure?
3. Why is K=10 (1-3-6 tetrahedron) the most stable geometry? What makes K=2..9 unstable?
4. How exactly does spin promote the phase shift / wavefront distortion that creates charge? Is it the L→T conversion itself, the toroidal geometry, or the core compression?
5. Is flux or gradient the correct force computation for each regime (near-field gradient vs far-field flux)?
6. Does the Yee & Hauger λ profile (longest near core) or LaFreniere profile (shortest near core) produce the correct physics? Or do they describe different regions (core vs shells)?
7. Can the sinc lock-in + charge discrimination from 2D flux be unified into a single force law that gives strong force at near-field and Coulomb at far-field?
8. Does the standalone particle's far-field traveling wave naturally produce Coulomb, or is the composite particle level (proton) needed?
