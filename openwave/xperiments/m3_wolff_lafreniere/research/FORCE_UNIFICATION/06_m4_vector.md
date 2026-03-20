# M4 VECTOR WAVE

## Elliptical Motion Geometry

![ellipse](images/granule_ellipse_small.gif)

### Displacement Always Traces an Ellipse

When each xyz component oscillates at the same ω:

```text
x(t) = R_x · cos(ωt + Φ_x)
y(t) = R_y · cos(ωt + Φ_y)
z(t) = R_z · cos(ωt + Φ_z)
```

The trajectory at any grid point is **always an ellipse** lying on some plane in 3D space. This is mathematically guaranteed — the same principle as polarization of light (linear, circular, elliptical are the only possible states for a monochromatic wave).

Degenerate cases:

- **All phases equal** (Φ_x = Φ_y = Φ_z): linear oscillation (ellipse collapsed to a line)
- **Phases differ by π/2, equal amplitudes**: circle
- **General case**: ellipse on a tilted plane

### 6 Numbers Fully Describe the Ellipse

The ellipse at each voxel is fully described by just 6 numbers: `R_x, R_y, R_z, Φ_x, Φ_y, Φ_z`.

These encode:

- **Shape:** the ratio of amplitudes (R_x, R_y, R_z) determines elongation
- **Orientation:** the phase differences (Φ_x - Φ_y, Φ_x - Φ_z) determine the tilt of the ellipse plane in 3D
- **Size:** the magnitudes of R_x, R_y, R_z

No other information is needed (plus the shared ω for temporal evolution).

### Drawing the Ellipse Without Time Sampling

Decompose the complex phasor vector into real and imaginary parts:

```text
a = (R_x·cos(Φ_x),  R_y·cos(Φ_y),  R_z·cos(Φ_z))
b = (R_x·sin(Φ_x),  R_y·sin(Φ_y),  R_z·sin(Φ_z))
```

These two vectors **are** the ellipse:

- **a** and **b** span the plane the ellipse lies on
- Normal to that plane: `a × b`
- Semi-major and semi-minor axes: diagonalize the 2×2 matrix formed by `a·a`, `a·b`, `b·b`

From 6 numbers → full geometric ellipse (center, plane, axes, orientation) — no time sampling needed.

## Extracting the 6 Phasor Values

### Method 1: Phasor Superposition (Analytical)

Sum complex phasors from all sources in one pass. Instant and exact.

### Method 2: Quadrature Sampling from Timestep Simulation

Sample displacement at two times separated by a quarter period (`t₀` and `t₀ + π/(2ω)`):

```text
s1 = displacement_x at t₀
s2 = displacement_x at t₀ + π/(2ω)

R_x = √(s1² + s2²)
Φ_x = atan2(-s2, s1) - ωt₀
```

Repeat for y and z components.

### Comparison

| Method | How | Cost |
| --- | --- | --- |
| Phasor superposition | Sum complex phasors from all sources, one pass | Instant, exact |
| Quadrature sampling | Run sim for at least a quarter period, sample twice | Requires simulation time, approximate if sources are moving |

The phasor route gives them analytically in one shot. The timestep route requires the sim to actually run and is only accurate if sources are stationary during the sampling window.

## M4/ POSSIBLE SOLUTION: Spin as the Key to Magnetic and Thermal Phenomena

**Spin** may be the unifying concept connecting magnetic force and thermal energy:

- **Magnetic moment from spin**: When longitudinal energy waves hit a wave center, part of the wave energy drives the wave center into rotation (spin). This rotational motion converts longitudinal wave energy into transverse wave emission — the magnetic component. The magnetic moment of a particle is directly tied to its spin rate and axis
- **Thermal energy from spin**: If thermal energy is wave-based (standing wave dynamics within the particle), then spin modulation — changes in spin rate or precession — could be the mechanism by which thermal energy is stored and transferred. Higher temperature = higher spin modulation = more energetic transverse wave emission
- **The connection**: Both magnetic fields and thermal energy may emerge from the same physical process (spin), manifesting differently depending on whether the spin is coherent (magnetic: aligned spins produce net transverse wave field) or incoherent (thermal: random spin axes produce isotropic energy distribution)

**Toroidal flow geometry (Smoliński)**: The spin/magnetic/thermal connection gains a concrete geometric picture from Smoliński's soliton model. The internal "engine" of a particle operates in the Energy Domain where wave energy circulates in **toroidal (doughnut-shaped) flows** governed by non-linear r⁵ scaling. This toroidal circulation is the physical mechanism behind spin:

- **Spin** = toroidal wave flow around the particle core, not an abstract quantum number
- **Magnetic moment** = natural dipole field from the circulating toroidal current (exactly as a current loop produces a magnetic field in classical EM)
- **Thermal energy** = modulation of the toroidal flow rate/precession — higher temperature means more energetic toroidal circulation
- **The m4 elliptical trajectories** observed in phasor analysis may be cross-sections of these toroidal flows — the ellipse semi-major axis captures the longitudinal (electric) component, the semi-minor axis captures the transverse circulation (magnetic/thermal)

**Vortex model of the electron (Butto, 2021)**: The toroidal flow picture gains further support from Butto's vortex theory (see [A New Theory for the Essence and Origin of Electron Spin](references/Spin.pdf), which models the electron as a **superfluid irrotational vortex** in the vacuum medium. Key insights for OpenWave:

- **Spin as differential rotation**: The vortex core rotates at twice the rate of the vortex boundary. The 1/2 spin quantum number is literally the ratio between core and boundary rotation rates — not an abstract quantum number but a measurable gap in rotational velocity. The core completes 720° while the boundary completes 360°
- **Planck constant = vortex angular momentum**: The circulation Γ = 2πrc is conserved in an irrotational vortex. Multiplied by the electron mass: Γ·m_e = 2πr_e·c·m_e = h (Planck's constant). This connects the fundamental quantum of action to the physical circulation of the vortex
- **Magnetic moment without g-factor**: The vortex model derives the electron magnetic moment directly from the circulating charge current μ = qcr, without needing the anomalous g-factor correction. The ratio μ/L = q/m_e emerges naturally from vortex hydrodynamics
- **Helmholtz vortex theorems apply**: (1) vortex strength is constant along its length, (2) a vortex filament cannot end in the fluid — it must form a closed loop or extend to boundaries, (3) an irrotational fluid stays irrotational. These constrain the topology of possible particle structures — closed vortex loops = toroidal geometry
- **Streamline curvature from differential rotation**: The curved streamlines of the vortex arise from the velocity differential between core and boundary — the faster-spinning core "drags" the surrounding medium, creating spiral arms. This is the physical origin of the wave pattern near the particle center

**Connection to OpenWave**: The vortex model suggests that our standing wave structure near the wave center may actually be a cross-section of a vortex flow pattern. The "standing wave" appearance could arise from the stable circulation pattern of the vortex, where the interference of in-waves and out-waves creates the stationary nodes we observe. In the m4 vector wave method, this vortex circulation should manifest as the toroidal component of the elliptical displacement trajectories.

This research direction requires the m4 vector wave method, as toroidal/vortex flows, spin, and transverse waves cannot be represented in the scalar m3 framework. The non-linear r⁵ scaling inside the Energy Domain also implies variable-wavelength wave equations (λ(r)) near the core, connecting to the [Yee & Hauger](references/Spin.pdf) and WKB research already planned for Phase 3.

## Why Scalar is Insufficient: Monopole = Longitudinal Only

**Key insight**: Wolff's scalar wave model may be scalar precisely because a **single WC (monopole) produces only longitudinal waves** — radial motion, linear oscillations, radially oriented. A monopole has no ellipses, no transverse component, no magnetic field. The scalar model correctly describes a single isolated particle.

**Magnetic fields emerge from multiple particles**: when two or more WCs combine their waves, the interference of radial waves from different directions creates **transverse displacement components** — the elliptical granule paths observed in M4. The transverse component IS the magnetic field. It doesn't exist for a monopole — it emerges from multi-particle wave superposition.

**Why Wolff needs complex numbers (i)**: in the scalar equation, the imaginary unit i encodes a 90° phase-shifted (quadrature) component — a mathematical proxy for the transverse wave that can't be geometrically represented in a scalar field. In a vector model, this transverse component is naturally encoded in the displacement direction itself. The complex number is an algebraic shortcut for what is fundamentally a geometric (vector) quantity.

**Implication for force unification**: there may be only **one force** — the gradient of energy (F = -∇E). What we call electric, magnetic, and gravitational forces are the same energy gradient projected onto different directions:

- **Longitudinal component** of the gradient → electric force
- **Transverse component** of the gradient → magnetic force
- **Residual/bulk component** (from density deficit) → gravitational force

The direction of the gradient depends on the vector structure of the wave field at each point. A scalar model collapses all directional information into magnitude — which is why it can produce force magnitude (1/r² scaling) but not force direction (charge sign).

## Per-Component Amplitude (Not Magnitude-Only)

Amplitude must be tracked **per spatial component** (A_x, A_y, A_z), not as a single magnitude |A|. The displacement vector already encodes this in M4 — but the current energy/force computation collapses it to |ψ|² = ψ_x² + ψ_y² + ψ_z² (magnitude squared), losing the directional information.

To recover force direction, energy and force may need to be computed **per component**:

```text
E_x = ρV(f·ψ_x_rms)²,  E_y = ρV(f·ψ_y_rms)²,  E_z = ρV(f·ψ_z_rms)²
F_x = -∂E_x/∂x,  F_y = -∂E_y/∂y,  F_z = -∂E_z/∂z
```

Or from vector operations on the full displacement field (divergence, curl, flux). The per-component approach preserves the directional structure that magnitude-only discards.

Soliton [Toroidal](03_additional.md#soliton-geometry-two-domains-and-toroidal-flows) relationship.

Lafreniere [magnetic field lines](03_additional.md#magnetic-fields-as-hyperboloid-wave-patterns) relationship.
