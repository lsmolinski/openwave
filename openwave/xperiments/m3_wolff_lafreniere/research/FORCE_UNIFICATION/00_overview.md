# CONCEPTUAL OVERVIEW

## GOALS: Force Unification from Wave Interference

### PRIMARY (EMERGENT FORCES)

- Provide numerical evidence that all fundamental FORCES (electric, magnetic, gravitational, strong) emerge from Energy Wave interference patterns in a spacetime medium — validating Energy Wave Theory (EWT) through simulation.

### SECONDARY (EMERGENT WAVES)

- Provide numerical evidence that ELECTROMAGNETIC WAVES and HEAT emerge from Energy Wave.

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

## The Spacetime Medium

An unknown but fluid-like substance permeates all of space and penetrates all matter. This medium vibrates in harmonic oscillations, forming longitudinal waves that travel through space at the speed of light (c). The medium has measurable properties:

- Density: ρ = 3.86 × 10²² kg/m³ (38.6 qg/am³ in simulation units)
- Wave speed: c = 2.998 × 10⁸ m/s (0.3 am/rs)
- Fundamental amplitude: A₀ = 9.22 × 10⁻¹⁹ m (0.92 am)
- Fundamental frequency: f₀ = 1.05 × 10²⁵ Hz (0.0105 rHz)
- Fundamental wavelength: λ₀ = 2.85 × 10⁻¹⁷ m (28.5 am)

## Wave Centers and Particle Formation

Wave centers are points where waves reflect and change character. Key properties:

- **Reflection**: An isotropic base energy wave is always present at any point in space — the result of reflections from all matter in the universe, reaching every point from all directions (Huygens wavelet principle)
- **Re-emission**: Every wave center reflects this base wave outward radially, changing its phase and character
- **Phase offset (source_offset)**: Each wave center has a phase offset that determines its charge signature. `cos(0) = +1` (positron), `cos(π) = -1` (electron). This is the origin of particle charge
- **Standing waves**: Near the wave center, reflected in-waves and emitted out-waves superpose to form standing waves — the structure of the particle itself

### Matter Formation from Standing Waves & Force Equilibrium

Once force unification is demonstrated, matter formation can be explained as the natural consequence of attraction/repulsion forces reaching equilibrium — keeping or separating particles to form the chain of material bonding. The same wave interference forces that create energy gradients would drive the assembly of increasingly complex structures: subatomic particles → electrons → protons/neutrons → atoms → molecules → bulk matter. Each level of organization emerges from force balance at the level below.

## Energy and Force

The EWT energy equation for a volume of medium:

```text
E = ρ · V · (f · A)²
```

Where ρ is medium density, V is volume, f is frequency, A is displacement amplitude.

Force is the negative gradient of energy density:

```text
F = -∇E
```

Force emerges from the **gradient of energy density** — wherever wave interference creates a spatial variation in energy, there is a force. Wave centers (which have mass) move toward lower energy configurations via F = ma. In essence, particles fall into low-energy valleys in the energy density landscape shaped by wave interference.

Computing F directly from ∇E (rather than expanding the chain rule) means that when ρ(x), f(x), or λ(x) become spatially variable in later phases, the force computation automatically captures all contributions without modification.

This may be the physical mechanism behind what Einstein described as "curvature in spacetime" — not a geometric bending of an abstract manifold, but a real energy density landscape sculpted by wave interference in the medium. The "curvature" is the shape of the energy valleys and hills created by constructive and destructive interference. Particles follow geodesics not because space is curved, but because they roll downhill in the energy density field.

## Force Unification Hierarchy

The theory proposes that all forces are manifestations of the same wave interference — ultimately, there may be **one single force** (F = -∇E, the gradient of energy). What we call electric, magnetic, and gravitational forces are the same energy gradient projected onto different directional components of the vector wave field. The force type depends on the angle/direction, not on separate mechanisms:

1. **Electric Force** (longitudinal): The fundamental force. Arises from constructive/destructive interference of energy waves between wave centers with different phase offsets. Opposite charges (0 vs π offset) create destructive interference between them → lower amplitude in the gap → particles fall into the energy well → attractive force. Same charges create constructive interference → higher amplitude between them → particles repelled from the high-energy zone → repulsive force

2. **Magnetic Force** (transverse): When longitudinal energy waves hit a spinning wave center, the energy required for spin is converted into a transverse wave component (90° to longitudinal). This reduces the longitudinal (electric) wave energy. Magnetic force emerges from the transverse wave interference

3. **Gravitational Force** (shading / push-out): When many wave centers cluster together (large bodies), they collectively absorb/scatter incoming energy waves. The region around the body has reduced longitudinal wave amplitude — a "shadow". Another body in this shadow region experiences a net force toward the first body because the energy density gradient points inward. Gravity is a residual effect of large-scale electric wave shading. An alternative but compatible formulation (Smoliński) treats gravity as a **pressure deficit (buoyancy)** in the medium: solitons (particles) displace medium density from their interior, creating a push-out effect. The 10⁻⁴² electromagnetic-to-gravitational force ratio emerges geometrically from a density hierarchy in the lattice packing (see Smoliński's Contributions).

4. **Strong Force**: Longitudinal electric force in short distances (λ scale). Particles inside each others standing wave radius will experience the strong force, a concentrated electrical force with stored energy from standing waves.

5. **Orbital Force**: a combination of electric and magnetic forces creating zero amplitude nodes (orbitals), where electrons can be found randomly, doesn't follow a classical planetary orbit.
