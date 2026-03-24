# CONCEPTUAL OVERVIEW

## GOALS: Matter Formation & Force Unification from Wave Interference

### PRIMARY (EMERGENT PARTICLES)

- Provide numerical evidence that all PARTICLES emerge from Energy Wave interference patterns in a spacetime medium — validating Energy Wave Theory (EWT) through simulation.

### SECONDARY (EMERGENT FORCES)

- Provide numerical evidence that all fundamental FORCES (electric, magnetic, gravitational) emerge from Energy Wave interference patterns in a spacetime medium — validating Energy Wave Theory (EWT) through simulation.

### TERTIARY (EMERGENT WAVES)

- Provide numerical evidence that ELECTROMAGNETIC WAVES and HEAT emerge from Energy Wave.

## CURRENT STATE

### NEAR-FIELD

| Effect | Status | Evidence |
| --- | --- | --- |
| Strong Force (lock-in) | ✅ Emerges | WC at Same phase: sinc nodes create energy wells at λ/2, WCs lock into these wells — this IS the strong force binding mechanism for quarks, nuclear binding, and particle formation |
| Anti-Particle Annihilation | ✅ Emerges | WC at Opposite phase: deepest well at r=0 (zero separation = complete wave cancellation), barriers at λ/2 explain positronium, high kinetic energy overcomes barriers → direct annihilation |
| Particle Formation | ⚠️ Concept works | Lock-in holds WCs, but M3/M4 K=10 tetrahedron unstable (needs variable λ for non-uniform nodes in equilibrium state) |

### Next to do

- solve K=10 tetrahedron stability on M3/M4
- rename M3 xperiments with near-field results.

---

### FAR-FIELD

| Effect | Status | Evidence |
| --- | --- | --- |
| Electric Force | | Charge emerges, Coulomb validation, ❌ spatial_term flips force direction |
| Magnetic Force | | Spin |
| Gravitational Force | | Shading |

### Candidates for testing

#### Energy Flux concept vs. F = -∇E

- F ∝ momentum flux: Poynting vector analog
- Coulomb force as radiation pressure from the traveling waves, not as energy gradient

## ⚠️ Key Hypothesis Under Investigation

- **BASE WAVE**:
  - The fundamental energy wave that fills all of space — the energy source for matter formation, force emergence, EM waves, and heat. Waves arriving from all directions create an isotropic field. In 1D, validated as a standing wave (`ψ = A₀·cos(kx)·cos(ωt)`) via Laplacian self-stabilization. In 3D, the isotropic superposition from all directions creates a uniform energy density field with vector displacement at every point. The base wave equation in 3D and its exact form are still to be determined — this is Step 1 of Phase 1c
  - TESTING: standing wave, quadrature, uniform, stochastic, laplacian, dual-channel

- **WAVE CENTERS (WCs)**:
  - Locations in 3D space where the base wave is elastically disturbed — the wave passes through and comes out changed (not just reflected). Phase 1b tested 10 disturbance models: passive (reflect, absorb, clamp, scatter) all failed; elastic L→T spin conversion is the only charge-sensitive mechanism found. The WC equation and exact disturbance mechanism are still being explored — this is Step 2 of Phase 1c
  - TESTING: (wave equation) in + out, partial standing, out wave only, A(r), λ(r)
  - sinc function: `ψ = A · cos (kr ± ωt ± φ) / kr`
  - TESTING: (interaction with base) additive reflection, absorber, elastic disturbance (amp, phase shift, spin)
  - effect: standing @near_field + traveling @far_field

- **3D SPACE & VECTOR WAVES**:
  - In 3D, spherical wave interference from the base wave + WC disturbance always promotes elliptical granule motion — even one WC + base wave is enough. Energy is not scalar in 3D — it requires a vector field with independent longitudinal (L) and transverse (T) components: `E = E_L + E_T`. Scalar models (M3/1D) collapse directional information into magnitude, producing correct 1/r² scaling but wrong force direction. Vector displacement is required to capture charge-dependent force — this is why Phase 1c must work in 3D

- **FORCE MECHANICS**:
  - One unified force `F = -∇E`: Force emerges from energy gradient
  - Two directions. `∇E_L` = electric force (radial/longitudinal gradient), `∇E_T` = magnetic force (perpendicular/transverse gradient). L/T defined relative to radial direction from the WC experiencing the force. Electric always present (one radial direction). Magnetic conditional — transverse has 360° freedom, cancels unless aligned/coherent (moving electrons, permanent magnets, spin coherence). Gravity = residual total energy deficit (omnidirectional, always attractive, weak)
  - `E = ρV·(c/λ·A)²`: Energy is a function of Α and λ (ρ might also be variable, c is absolute)
  - Wave character and superposition (from wave equations) determines A(r) and λ(r) for every point in space.
  - Wave oscillator as `f(r,t,φ)`: phase = spacial_term + temporal_phase + source_offset.
  - Major clarity gained on what produces what:
    - wave phase → annihilation vs. repulsion
    - standing wave / sinc → strong force / lock-in
    - UNKNOWN → electric / charges
    - spin → magnetic
    - drainage → gravitational

- **3D WAVE DISPLACEMENT: ELLIPTICAL MOTION, SPIN, HYPERBOLIC PATTERNS**:

  - granule elliptical motion, hyperbolic patterns
  - needs vector field with orthogonal components
  - 3D force orthogonality might connect to electromagnetism.
  - The ellipse at each point encodes L amplitude (semi-major), T amplitude (semi-minor), handedness (CW/CCW rotation direction = spin up/down = electron/positron), and ellipse plane orientation. Described by 6 phasor numbers (R_x, R_y, R_z, Φ_x, Φ_y, Φ_z). Complex sinusoids in QM encode this naturally: real = L, imaginary = T, |ψ|² = E_L + E_T. The imaginary unit i IS the 90° quadrature relationship between L and T components. Spin is the L→T wave transformation at the WC (Wolff: "the wave is spinning, not the particle"), with 720° spherical rotation required by 3D space geometry

- **NON-LINEAR WAVES**:
  - Variable λ(r) near WCs — wavelength changes with distance from the WC (Yee & Hauger shells `r_wavelength = 2Κλ - 2nλ`, WKB phase integral). The current energy formula `E = ρV(fA)²` uses constant f and can't see λ variation. Phase 1d implements `E = ρV(c·A/λ(r))²` where the `∇λ` term creates force from wavelength gradients. Converges with Phase 1c — variable λ breaks sinc periodicity while vector displacement provides charge sensitivity

- **ENERGY CONSERVATION & WAVE STEEPNESS**:
  - Wave steepness A/λ = constant for isolated energy redistribution. External energy input (heating) increases steepness. WC spin converts L→T while conserving `E_L + E_T = const`. The conversion ratio may be the fine-structure constant α

- **VISUALIZATION & RENDERING**:
  - Phase 1c focuses on raw math (numpy scripts) for an accelerated discovery phase — compute, sweep, analyze, document. No matplotlib animation, no 3D rendering. When the wave equations are validated and force emergence is confirmed, the equations get ported to OpenWave's M4 Taichi-based 3D rendering infrastructure (launcher, medium, wave_engine, particles, force & motion) for visualization and public demonstration

- **VALIDATION TARGETS** (the complete force unification validation set):
  - **Near-field (particle regime)**: same-phase lock-in — oscillatory force creates stable energy wells (quarks, orbital shells, bonding). Opposite-phase monotonic attraction → annihilation (wave cancellation at zero separation). Near-field to far-field transition boundary — clean crossover from lock-in to Coulomb. Note: the sinc force-flipping that's the PROBLEM in far-field is the DESIRED PHYSICS in near-field — at sub-wavelength distances the oscillating force direction creates the energy wells that lock particles together. The challenge is making it stop at the transition boundary
  - **Matter formation**: attraction + repulsion equilibrium holds particles together. Standing wave core forms at WC with correct radius (K²λ). Mass = energy in standing waves (E = mc²)
  - **Far-field (Coulomb regime)**: force direction emergent (not imposed ±1 — passes emergence test). 1/r² magnitude scaling vs Coulomb reference. Newton's 3rd law (equal and opposite). Opposite charge attracts, same charge repels — at ALL separations (no sinc flips). Energy conserved
  - **Magnetic**: transverse force component (∇E_T) reproduces magnetic field geometry. Spin alignment → coherent T field → detectable magnetism. Perpendicular to electric (∇E_L) at 90°
  - **Gravitational**: 10⁻⁴² EM-to-gravitational force ratio emerges from the model (shading / density deficit). Computed G matches Smoliński's Scilab reference values

## ⚠️ Emergence Criteria

Matter, forces need to emerge from wave physics (not be imposed):

1. The force must come from **wave interference** — the constructive/destructive pattern of actual oscillating waves
2. The charge sign must enter through the **phase relationship** between waves (how they interfere), not as a ±1 label on a smooth function
3. The force direction must depend on **both** charges' phases interacting, not on either one individually

## Force Regime Matrix

| Regime     | Same Phase                        | Opposite Phase                                |
| ---------- | --------------------------------- | --------------------------------------------- |
| Near-field | Lock-in (quarks, orbits, bonding) | Attraction → annihilation (wave cancellation) |
| Far-field  | Constructive → repulsion          | Destructive → attraction                      |

## MAIN CONCEPTS

Spacetime is an elastic fluid structure, vibrating in harmonic oscillations at extremely high frequencies (~10²⁵ Hz), that allows energy to transfer from one point to the next, forming an **energy field**. Through this fluid, waves travel in all directions and get reflected by wave centers that create disturbances. Those energy waves can become **standing** (forming matter) or remain **traveling** (transmitting all forces).

Everything emerges from disturbances in this wave field:

- **Particles (wave centers)**: localized disturbances that reflect and re-emit waves, changing the wave character and shifting energy spatially. The standing wave structure near a wave center IS the particle
- **Forces**: energy gradients created by wave interference between wave centers. These gradients generate all known forces — electric, magnetic, gravitational, strong, orbital — and hold matter together at every scale
- **Photons**: traveling wave packets — disturbances propagating through the medium. Unlike static force fields, photons are discrete packets of energy that travel through space and can apply force on distant particles upon absorption
- **Heat**: disturbances in standing wave form — thermal energy may be encoded in the amplitude or frequency modulation of standing waves within particle structure, rather than in bulk kinetic motion (see Thermal Energy Hypothesis)
- **Time**: not a fundamental dimension but an emergent property — the local rate of change determined by wave frequency (f = c/λ). Where λ is shorter, things change faster; where λ is longer, things change slower. Time dilation near massive bodies emerges from λ modulation by wave centers (see Time Dynamics)

**Simulation model**: The medium is modeled as a particle-based fluid simulation where:

- **Granules** = discretized fluid particles representing the medium
- **Granule motion** around equilibrium traces a circle (or ellipse under multi-source interference)
- **Vector displacement** captures both position and velocity of each granule
- **Cycle period** is related to wavelength — this is the source of time itself (time = rate of change in the medium)

All of these are proposed concepts within Energy Wave Theory. OpenWave's goal is to **numerically validate** them through computational analysis and rendered simulations — starting with force emergence from wave interference, then extending to matter formation, photon emission, and thermal dynamics.

## The Spacetime Medium & The Fundamental Energy Wave

An unknown but fluid-like substance permeates all of space and penetrates all matter. This medium vibrates in harmonic oscillations, forming waves that travel at the speed of light (c). The medium has measurable properties:

- Density: ρ = 3.86 × 10²² kg/m³ (38.6 qg/am³)
- Wave speed: c = 2.998 × 10⁸ m/s (0.3 am/rs)
- Fundamental amplitude: A₀ = 9.22 × 10⁻¹⁹ m (0.92 am)
- Fundamental frequency: f₀ = 1.05 × 10²⁵ Hz (0.0105 rHz)
- Fundamental wavelength: λ₀ = 2.85 × 10⁻¹⁷ m (28.5 am)

## Energy, Force, and Spacetime Curvature

The EWT energy equation: `E = ρ · V · (f · A)²`. In vector displacement (3D): `E = E_L + E_T = ρV(f·A_L)² + ρV(f·A_T)²`.

Force: `F = -∇E` — wherever wave interference creates spatial variation in energy, there is a force. Particles fall into low-energy valleys in the energy density landscape. Computing F from ∇E directly (not chain-rule expansion) means variable ρ(x), f(x), λ(x) are automatically captured.

This may be the mechanism behind "spacetime curvature" — not geometric bending of an abstract manifold, but a real energy density landscape sculpted by wave interference. The "curvature" is the shape of energy valleys and hills. Particles follow geodesics because they roll downhill in the energy field.

## Force Unification Hierarchy

One force `F = -∇E`, projected onto vector displacement components:

1. **Electric Force** (longitudinal, ∇E_L): constructive/destructive interference from WC phase offsets. Opposite charges (0 vs π) → destructive → energy well → attraction. Same charges → constructive → energy hill → repulsion. Always present (one radial direction)

2. **Magnetic Force** (transverse, ∇E_T): L→T spin conversion at WC creates transverse wave component (90° to longitudinal). Magnetic force from transverse energy gradient. Conditional — requires spin alignment/coherence (360° transverse freedom cancels unless aligned)

3. **Gravitational Force** (total deficit): WC clusters collectively drain energy from far field → amplitude "shadow" → net inward force. Smoliński: pressure deficit / buoyancy from medium density displacement. The 10⁻⁴² EM-to-gravitational ratio emerges from density hierarchy geometry. Two competing gravity models:

    - **Shading / directional attenuation**: WCs absorb part of in-wave, re-emit attenuated out-wave toward other particles → directional energy deficit → gravity. Likely requires M4 vector displacement (directional attenuation can't be represented in scalar M3)

    - **Buoyancy / medium density (Smoliński)**: Wave equation unchanged. WCs modify local ρ and f via λ modulation → granule velocity → pressure. Gravity emerges from E = ρV(fA)² where ρ and f become local variables. Works in M3 — only the energy calculation changes

For 1D sandbox: Phases 1-3 use constant ρ and f. Phase 5 adds option for position-dependent ρ(x) and f(x) computed from the wave field.

4. **Strong Force**: electric force at sub-wavelength distances (λ scale). Particles inside each other's standing wave radius experience concentrated energy from standing waves

5. **Orbital Force**: combination of electric + magnetic creating zero-amplitude nodes (orbitals) where electrons can be found — not classical planetary orbits

## Wave Centers and Matter

Wave centers elastically disturb the base wave — the wave passes through and comes out changed. Key properties:

- **Base wave interaction**: isotropic energy wave from all matter in the universe (Huygens wavelet principle) is disturbed at each WC
- **L→T spin conversion**: the disturbance transforms longitudinal → transverse (the physical meaning of spin). CW = electron, CCW = positron
- **Standing waves**: near the WC, in-waves and disturbed out-waves superpose → standing wave structure = the particle itself
- **Matter signature frequency**: every element has a unique frequency combination (like timbre). Discovering signatures enables wave-based matter manipulation
- **Matter formation**: force equilibrium (attraction + repulsion) assembles structures at every scale: subatomic → atoms → molecules → bulk matter

- Also refer to `../CLAUDE.md` file to search for any available context to the OpenWave project in a parent directory.

## SIMULATOR TOOLS

- 1D & 3D wave engines
- Wave Equations compute on GPU
- Multiple methods & options for wave equations, on scalar & vector fields
- Phasor superposition for exact analytical amplitude computation
- Variable λ logic (non-linear wave equations)
- 3D Rendering system
  - subatomic physics simulation with wave & particle 3D animation (leap-frog integrator)
  - flux mesh cross-section for variable display (color and magnitude via warp-mesh)
  - xperiment parametric design
  - data dashboard & instrumentation for data collection
