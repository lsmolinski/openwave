# PHASE 2: BUILDING BLOCKS

## WHERE WE STAND

Phase 1 exhaustively tested whether Coulomb force emerges from two single wave centers (K=1 or K=10). Key results:

| Effect | Status | Mechanism |
| --- | --- | --- |
| Strong force (lock-in) | ✅ Emerges | Sinc `cos(k·Δr+Δφ)` creates energy wells at λ/2 — same-phase WCs lock in |
| Annihilation | ✅ Emerges | Opposite-phase WCs: deepest well at r=0, barriers at λ/2 (positronium) |
| Charge discrimination | ✅ Exists | 2D flux: 100% same≠opposite at ALL separations (22/22) |
| Consistent Coulomb direction | ❌ Not from single WCs | Sinc flips force every λ/2 — intrinsic to monochromatic spherical interference |
| Correct Coulomb sign | ❌ Inverted | Isolated interaction: same→ATT, opposite→REP (180° off Coulomb) |

**Critical realization**: spherical 3D waves from single WCs produce lock-in and annihilation (near-field) but do NOT produce Coulomb alone. The sinc oscillation IS the correct near-field physics — it was never a bug. Something else is needed for far-field forces.

Full Phase 1 results: [01z_results.md](01z_results.md)

---

## KEY INSIGHT: THE FORCE HIERARCHY

**Invert the path.** Phase 1 tried to go directly from single WCs to Coulomb. That was wrong. Forces emerge in a hierarchy — each building on the previous:

```text
1. STRONG FORCE (fundamental, from wave interference)
   └── sinc lock-in: same-phase WCs lock at λ/2 nodes
   └── annihilation: opposite-phase WCs attract to r=0

2. PARTICLE FORMATION (from strong force lock-in)
   └── multiple WCs (K>1) locked in standing wave geometry
   └── neutrino (K=1) → electron (K=10, tetrahedral 1-3-6)
   └── standing waves reorganize into concentric spherical nodes

3. CHARGE EMERGENCE (from particle geometry, NOT from single WC phase)
   └── charge ≠ phase. Phase → lock-in/annihilation
   └── charge = property of the COMPOSITE particle (multi-WC arrangement)
   └── the specific geometry (tetrahedral, with spin) creates the wavefront
       disturbance that produces Coulomb

4. ELECTRIC FORCE (from charge, from particle)
   └── Coulomb emerges from the composite particle's wave disturbance
   └── not from individual WC phase — from the COMBINATION

5. MAGNETIC FORCE (from spin of the composite particle)
   └── L→T conversion drains longitudinal energy into transverse
   └── spin of multi-WC particle, not single WC

6. GRAVITATIONAL FORCE (from the spin energy deficit)
   └── L→T drainage creates a longitudinal amplitude deficit
   └── deficit = gravitational shading. Ratio ≈ α accumulated over K WCs
```

**Phase ≠ Charge:**

- **Phase** (source_offset 0 vs π): determines matter vs antimatter. Same phase → strong force lock-in. Opposite phase → annihilation. This is the NEAR-FIELD behavior, already validated in Phase 1
- **Charge**: emerges from the GEOMETRY of a multi-WC particle. The electron's 1-3-6 tetrahedral arrangement, combined with spin, creates a specific wave disturbance pattern that produces far-field Coulomb force. Charge is a compounding effect of same-phase WCs locked together

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

### 4. Charge from Multi-WC Geometry

The electron (K=10) has 10 same-phase WCs in a 1-3-6 tetrahedral arrangement. This specific geometry:

- Modifies the wavefront shape (no longer pure spherical — toroidal/elliptical from spin)
- Creates a phase shift (LaFreniere: core λ/2 shift from compression)
- Produces the far-field disturbance pattern that IS the electric field

Single WC (neutrino) = neutral, no spin, no charge. Multi-WC particle = charged, spinning, creates Coulomb.

### 5. Flux vs Gradient Force

Phase 1 showed:

- `F = -∇E` (gradient): charge-blind when base wave present, wrong sign for isolated interaction
- `S = -c²·ψ·∇ψ` (flux/radiation pressure): 100% charge discrimination in 2D

LaFreniere explicitly uses radiation pressure, not energy gradient, for Coulomb. Evaluate which force mechanism applies at each regime.

### 6. Granule Elliptical motion

We want to figure out the origin of fundamental forces.So first, force is a vector field and we know there is electric and magnetic forces and they have a 90 degree relationship.Gravitational force is radio towards large mass.We have the right hand rule for electromagnetic, electromagnetism and also the three finger pointing force relationship with electric field and magnetic field.We also have the fundamental elliptical track, a medium component, trajectory, traces.When three dimensional spherical waves superpose, this ellipse can be decomposed into vector amplitudes projected over a radial vector from a point that is experiencing the force.The ellipse projected creates the force fields and there is also the energy gradient needed to compute forces.So we need to decompose energy in a vector field as well, electrical energy, magnetic energy, vectors or scalar energy and the vector is the force.

---

## ATTACK PLAN

### Phase 2a: Near-Field — Particle Lock-In (M3/M4 3D)

**Goal**: Demonstrate stable particle formation in the 3D simulator. Start with the simplest case and build up.

- [ ] Stabilize K=10 electron tetrahedron on M3/M4
  - Implement variable λ(r) from Yee & Hauger shells (non-uniform node spacing)
  - Test whether variable nodes accommodate the 15/45 non-node pairs (√3×λ/2, √2×λ/2)
  - Leapfrog integrator already in place, damping at 0.995
- [ ] Validate same-phase lock-in in 3D animation
  - Two WCs, same phase, observe oscillation in energy wells
  - Measure well depth, oscillation period, escape energy
- [ ] Validate opposite-phase annihilation pathway
  - Two WCs, opposite phase, observe attraction through barriers
  - Measure barrier heights at λ/2 intervals, compare with positronium lifetime data
- [ ] Test multi-WC lock-in: K=2, K=3, ..., K=10 progressive build-up
  - At what K does the structure become stable?
  - Does the geometry self-organize or must it be prescribed?
- [ ] Characterize near-field → far-field transition boundary
  - Where do standing waves end and traveling waves begin?
  - How does this relate to particle radius K²λ?

### Phase 2b: Far-Field — Charge and Coulomb (Research Scripts + M3/M4)

**Goal**: Determine how charge and Coulomb force emerge from the composite particle.

- [ ] Investigate 3D flux-based force
  - Extend 2D flux test to full 3D spherical integration
  - Does solid-angle averaging smooth the sinc absolute direction?
  - Test at K=10 particle scale (100λ radius)
- [ ] Test composite particle (K=10) as Coulomb source
  - Place one stabilized K=10 electron, measure far-field force on a test WC
  - Does the tetrahedral geometry produce consistent Coulomb direction?
  - Compare with single WC (should show the difference)
- [ ] Investigate charge mechanism
  - How does the 1-3-6 geometry create the wavefront disturbance?
  - Role of LaFreniere core phase shift (λ/2 from compression)
  - Does variable λ(r) inside the particle create the charge sign?
- [ ] Validate Coulomb
  - Direction: same charge repels, opposite attracts at ALL separations
  - Magnitude: 1/r² scaling
  - Symmetry: Newton's 3rd law (equal and opposite)

### Phase 2c: Spin, Magnetic, and Gravitational Forces

**Goal**: Demonstrate the remaining forces from spin and energy deficit.

- [ ] Magnetic force from spin
  - L→T conversion at K=10 particle (η=α≈1/137)
  - Transverse energy gradient perpendicular to WC axis
  - Test spin alignment coherence → net magnetic field
- [ ] Gravitational force from spin deficit
  - L→T drainage creates longitudinal amplitude deficit
  - Test whether deficit accumulates over K WCs (factor of K×α?)
  - Compare with 10⁻⁴² EM-to-gravitational ratio
  - Validate against Smolinski's Scilab reference values
- [ ] Wave shading with particle clusters
  - Multiple particles → accumulated energy deficit → gravity

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

1. Does the electron tetrahedron stabilize with variable λ(r) alone, or does it need vector forces too?
2. Is the base wave required for lock-in, or do WCs alone create the standing wave structure?
3. How does the 1-3-6 geometry specifically produce charge? What is the minimal geometry that has charge?
4. Is flux or gradient the correct force computation for each regime (near-field vs far-field)?
5. Does the Yee & Hauger λ profile (longest near core) or LaFreniere profile (shortest near core) produce the correct physics? Or do they describe different regions (core vs shells)?
6. Can the sinc lock-in + charge discrimination from 2D flux be unified into a single force law that gives strong force at near-field and Coulomb at far-field?
