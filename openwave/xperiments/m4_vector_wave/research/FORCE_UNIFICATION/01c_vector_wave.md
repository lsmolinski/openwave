# PHASE 1c: Vector Wave Force

## ❌ The Problem

The main blocker is still the far-field oscillatory force: the sinc function sin(kr)/kr in the out-wave creates permanent nodes in the phasor RMS at all distances. Force direction flips every λ/2 of separation change, even where only smooth 1/r decay should exist. Confirmed in both 3D and 1D engines.

## The Base Wave (Fundamental Energy Wave)

The medium is not empty. A pre-existing **isotropic energy wave field** fills all of space — the fundamental longitudinal energy wave described by EWT. Its properties are known:

- Amplitude: A₀ = 9.22 × 10⁻¹⁹ m (0.92 am)
- Wavelength: λ₀ = 2.85 × 10⁻¹⁷ m (28.5 am)
- Frequency: f₀ = 1.05 × 10²⁵ Hz (0.0105 rHz)
- Density: ρ₀ = 3.86 × 10²² kg/m³ (38.6 qg/am³)
- Energy density: E₀ = ρV(fA₀)² — uniform everywhere (without WCs)

**What we don't know**: how granules oscillate (displace) in time. The base wave is NOT a uniform oscillation where everything goes up and down together (like water level). It must represent waves coming from **all directions** — an isotropic field with constant amplitude and energy density, but with granule displacements in **multiple phases** at each point.

**Possible nature**: a standing wave field everywhere — a fixed universal background. The medium stores energy in these standing waves as potential for everything: matter, forces, EM waves, heat, and even time itself (the displacement cycle at each point IS the local rate of change — the local clock).

**Base wave oscillation scale**: base wave λ = 28 am, while electron radius is 2800 am (100× λ). Particles don't "feel" base wave oscillations directly; they oscillate too fast at 10²⁵ Hz. Particles may only respond to the averaged-out RMS amplitudes.

**The base wave concept is not just context — it is the energy source.** The medium has stored energy (standing waves) as potential for matter, forces, EM waves, heat, and time. WCs redistribute this energy, and the redistribution pattern (affected by WC phase) creates the gradients that produce force.

The base wave (1b) provides the energy field that WCs disturb; vector displacement (1c) may carry the charge information that scalar magnitude discards; non-linear equations (1d) describe how the disturbance propagates. They may ultimately converge into a single solution.

**The open question**: can the reflected waves produce a far-field energy gradient that is smooth enough to avoid the sinc oscillation problem, while still being genuine wave interference?

### On the "Longitudinal" Assumption

The base wave is described as longitudinal in EWT literature, but **this is a hypothesis, not established fact**. In 1D simulation, we are forced to isolate a single wave mode (longitudinal displacement along x). This is a limitation of the dimensionality, not necessarily of the physics:

- **1D**: only longitudinal mode available — displacement along the propagation axis
- **2D**: two modes possible — longitudinal (along propagation) + transverse (perpendicular)
- **3D (reality)**: full elliptical motion — the granule displacement traces an ellipse in 3D space, characterized by: longitudinal amplitude, transverse amplitude, handedness (direction of granule motion around the elliptical track), and ellipse plane orientation

All of these properties can contribute to force direction. The unified force concept proposes that what we perceive as separate forces (electric, magnetic, gravitational) are actually **one 3D elliptical behavior** decomposed into components. At human scale, specific conditions make each component appear distinct — we named and described them as separate forces because that's how they manifest at our scale of mass, frequency, and inertial frame.

A 2D simulation could capture longitudinal + transverse, but **may not be sufficient** — the elliptical form can be oriented in multiple ways in 3D space, and this orientational freedom is likely essential for magnetic fields and spin. The 1D base wave work here is foundational (establishing the energy redistribution mechanism), but the full picture likely requires 3D (Phase 1c / M4).

#### Base Wave Contender Assessment (post-testing)

- **Quadrature** — confirmed as strongest base wave model. Flat energy, two orthogonal channels that map to L/T duality, traveling wave direction encodes handedness. In 1D the two-direction limitation (left/right) is a projection of 3D rotational freedom (CW/CCW elliptical). The quadrature captured something real: the two channels (cos/sin, real/imaginary, L/T) encode a **rotational degree of freedom** that 1D collapses into direction. This is exactly what complex sinusoids do in quantum mechanics — the imaginary unit i encodes a 90° phase relationship that represents physical rotation in 3D space. The quadrature base wave may be the correct 1D representation of the 3D wave field, with the two channels being L and T components. Carry forward into Phase 1c as the base wave model
- **Standing wave** — physically correct (Laplacian-validated), simplest model. Nodes exist at λ/2 but average out at particle scale (electron = 100λ). The node structure IS real at the base wave scale and may matter for sub-particle physics (neutrino K=1 has only 1λ core — nodes are significant at that scale). For force computation, the standing wave doesn't help with charge-dependent force — confirmed by all tests. The quadrature is a standing wave with a second orthogonal channel — it's a superset
- **Uniform + dual phase** — charge-blind in testing. Per-channel energy sum is symmetric regardless of which channel is boosted. Needs cross-channel coupling (like L→T) to break symmetry — which points back to the spin mechanism

## Wave Centers as Energy Redistributors

WCs do not emit waves. They do not inject energy. They create **disturbances** in the base wave field that **redistribute energy density**:

**Reflection / equilibrium (Jeff Yee on the mechanism)**:

> "Is the wave being reflected by a wave center, or is the wave center shifting to the point of equilibrium such that waves from the opposite side continue through at the same amplitude? The math for EWT would support #1 or #2, so it's hard for me to know if a wave center really does a reflection or it is just responding to the position of all waves (from all directions). In my writings, I favor #1, but it's still possible that those that believe a wave center is just the center point of where waves converge with equal amplitude."

Whether reflection (#1) or equilibrium positioning (#2), the result is the same: the base wave is disturbed, and the disturbance expands radially from the WC as a spherical wave of disturbed medium.

**Energy redistribution (globally conserved)**:

- **Near-field (inside particle radius, r < K²λ)**: energy is **concentrated** into the WC's own 3D spherical standing waves. These standing waves define the particle: its radius (K²λ), its mass (energy contained in standing waves, E = mc²), its identity. Core size depends on particle type via wave center count K: neutrino K=1 (1λ core), electron K=10 (100λ core). The particle's standing waves are radially oriented disturbances on the base wave — different from the base wave standing waves
- **Far-field (outside particle radius)**: energy is **drained** from the surrounding base wave field to supply the near-field concentration. This far-field energy deficit is the mechanism behind force and gravity

Both near-field and far-field are waves with oscillatory displacement (standing and traveling respectively). The disturbance decays with distance and restores to the undisturbed base wave far from the WC.

**How WC phase affects the far-field drainage**: the phase (source_offset: 0 = positron, π = electron) must affect the spatial pattern of the far-field energy drainage. But NOT via a simple ±1 sign multiplier — that was Phase 1a, and it's not emergent. The actual mechanism is **unknown and must be discovered**. This is the central open question of Phase 1b.

## Force Emergence from Energy Redistribution

1. WC1 disturbs the base wave → concentrates energy in its standing wave core → creates a far-field energy deficit (drainage) that radiates outward
2. This drainage reaches WC2's location and **disturbs WC2's standing waves** — warping the energy field around WC2
3. The warped energy field creates an **energy gradient at WC2's position**
4. F = -∇E → WC2 moves toward lower energy density → **force and motion**

The force direction depends on HOW WC1's drainage pattern interacts with WC2's standing waves — and this must depend on the phase relationship between the two WCs.

## Connections

**EWT / LaFreniere**: the base wave IS the "isotropic in-wave from all matter in the universe." LaFreniere's model describes WCs as reflecting incoming waves — the reflection creates local energy redistribution (standing waves near WC) and a far-field amplitude deficit.

**Gravitational shading**: the far-field amplitude deficit from energy redistribution IS the gravitational shading mechanism. WCs absorb/redirect base wave energy, creating a "shadow" in the far field. This connects to Smoliński's push-out / buoyancy model.

**WC disturbance scope**: the disturbance affects not just amplitude but also **wavelength λ** and potentially **density ρ** near the WC. This connects to the multi-variable energy gradient (∇A + ∇f + ∇ρ) and to non-linear wave equations (Phase 1d).

**Scalar base → vector emergence (hypothesis)**: the fundamental base wave might be scalar (longitudinal only) — but this is a hypothesis, not established. In 3D reality, granule displacement traces ellipses with longitudinal, transverse, handedness, and orientation components. Vector (transverse) waves may emerge from **spin** — the WC's toroidal wave rotation converts longitudinal to transverse. See [Base Wave Numerical Model](#on-the-longitudinal-assumption) for full dimensional analysis.

**Emergent wave hierarchy**:

- **Matter** = standing electromagnetic waves (concentrated base wave energy near WC)
- **Photons / EM waves** = traveling wave disturbances
- **Heat** = standing wave concentrated energy, related to spinning/magnetic momentum

**Dual-phase speculation**: the base wave may consist of two complementary phase modes. WCs lock onto one mode or the other (source_offset = 0 or π), creating the two charge states.

**Spin as longitudinal → transverse converter**: spin may convert longitudinal base wave energy into transverse wave components (720° spherical rotation = spin-1/2). The conversion ratio may be the fine-structure constant α. This connects to Smoliński's toroidal Energy Domain and Butto's vortex electron model.

## Open Questions

- What is the correct time-domain representation of the isotropic base wave in 1D?
- How does a WC's reflection/scattering produce the radial disturbance?
- What is the spatial structure of the far-field energy drainage?
- How does WC phase (0 vs π) modify the drainage pattern?
- Does the drainage itself oscillate (wave-like) or is it smooth?
- Can the drainage-drainage interaction between two WCs produce charge-dependent force direction?
- Does spin convert longitudinal to transverse at a fixed rate? Is the ratio α?

**Connection to spin**: "spin" is an observed quantum phenomenon never fully understood physically, and likely mislabeled as actual physical rotation of a particle. Wolff's insight: **"the wave is spinning, not the particle"** — spin is the transformation of the in-wave into the out-wave at the wave center, a property of 3D space (720° spherical rotation). Elastic disturbance may be the mechanism behind what we observe as spin — the WC changes the wave character as it passes through, and that change is what we measure as spin, magnetic moment, and charge.

**Status**: not ruled out — the physics is sound (λ variation creates energy gradients), but the current implementation framework can't capture it. Requires Phase 1d's variable-λ energy equation to test properly.

**Status**: most promising elastic model. Points strongly toward Phase 1c (full vector displacement with independent L and T components) as the path to resolve the remaining oscillation. The quadrature proxy demonstrates that L→T conversion produces charge sensitivity, but the full effect likely requires true two-component displacement, not phasor channel manipulation.

### Key Conclusions from Phase 1b

1. **Additive superposition cannot resolve the force direction problem** — regardless of base wave mode. The WC interference pattern dominates and produces the same sinc oscillation. The base wave only adds a constant phasor offset (Step 2b finding)

2. **Passive WC interactions fail in isotropic fields** — reflections cancel from all directions. Confirmed in M2 (3D, 12 experiments) and re-confirmed in 1D (weaker cancellation, same result). M2 prior art fully validated

3. **Elastic amplitude modulation is charge-blind** — symmetric scaling doesn't extract charge info from the wave. Elastic phase warping produces zero force — phasor rotation preserves RMS at every point (the energy formula uses constant λ, so it can't see phase changes)

4. **L→T spin conversion is the breakthrough mechanism** — it's the ONLY model that produces different behavior for opposite vs same charge. The charge-dependent direction of mode conversion (CW vs CCW) creates asymmetric energy landscapes. But the quadrature phasor proxy is limited — P² + Q² still produces sinc modulation. True independent L and T displacement components are needed

5. **The 1D scalar sandbox has a fundamental limitation** — it cannot represent true L/T conversion because it has only one displacement component. The phasor channels are mathematical, not physical. True L→T needs two independent displacement fields whose energies sum separately: `E = E_L + E_T`

6. **Phase warp (variable λ) needs variable λ in the energy equation** — the physics is real (λ varies near WCs), but `E = ρV(fA)²` with constant f can't see it. Phase 1d must implement `E = ρV(c·A/λ(r))²` where the `∇λ` term creates force from wavelength gradients

#### Path Forward: Phase 1c and 1d Converge

Both remaining paths are needed and will likely converge into one solution:

- **Phase 1c (vector displacement)**: independent L and T displacement fields → force from L/T balance via `E = E_L + E_T`. Makes elastic spin (Option G) fully distinguish charges without sinc

- **Phase 1d (variable λ)**: `E = ρV(c·A/λ(r))²` → force from `∇λ` (wavelength gradients). Yee & Hauger shells for λ(r), WKB phase integral. Makes elastic phase warp (Option F) produce force

**Recommended order**: Phase 1c first — the L→T spin is the only charge-sensitive mechanism discovered. It needs true vector displacement to work fully. Phase 1d (variable λ) can be layered on top of the vector engine, and the combined solution (variable λ + L→T spin in vector displacement) may be the unified answer.

**Phase 1b WRAP-UP tasks** (carry forward to Phase 1c/1d):

- [ ] Test energy redistribution: concentration near WC (r < K²λ), drainage in far field
- [ ] Determine how WC phase affects far-field drainage pattern (NOT via ±1 sign — mechanism must be discovered → L→T spin is the candidate)
- [ ] Test force emergence: drainage from WC1 disturbs WC2 → energy gradient → F = -∇E
- [ ] Validate against Coulomb reference (direction + 1/r² scaling)

### M2 Proposed Solution: Spin as L→T Conversion (14_spin_theory.md)

The M2 research proposed that the missing physics is **spin — the conversion of longitudinal waves to transverse waves** at the wave center. This is the key insight:

**Why L→T conversion breaks the symmetry:**

1. Incoming waves arrive from all directions — **pure longitudinal** (compression/rarefaction)
2. At the WC, some longitudinal amplitude converts to transverse amplitude
3. Outgoing waves are **mixed longitudinal + transverse** (reduced L, increased T)
4. The transverse component is a **NEW type of wave** that wasn't in the incoming field
5. Because in ≠ out (different wave character), the **symmetry that caused cancellation is broken**
6. Standing waves can now form from the L/T interference pattern

**Energy conservation**: `E ∝ ampL² + ampT² = constant`. The conversion redistributes energy between modes without creating or destroying it. The conversion ratio may be related to the **fine-structure constant α** (L/T coupling strength).

**Physical basis (Milo Wolff)**: spin is not the particle spinning like a ball — the **wave character is being transformed**. The in-wave must undergo a 720° phase shift (spherical rotation property of 3D space) to become the out-wave. This rotation is spin. Only two directions are possible (CW or CCW) → electron vs positron. Spin is a property of **3D space**, not of the particle itself — explaining why all charged particles have the same spin value.

**Connection to electromagnetism**: if L→T conversion is the mechanism, then:

- **Electric field** ← Longitudinal wave component
- **Magnetic field** ← Transverse wave component
- **90° geometric offset** because L and T are perpendicular by definition
- **Same phase** because they come from the same wave

This would explain why complex numbers are required in quantum mechanics: **real part = L, imaginary part = T**. The Schrödinger equation naturally encodes L/T duality. Probability `|ψ|² = ψ_real² + ψ_imaginary²` is the total energy from both components.

⚠️ **The M2 research points toward L→T conversion (spin) as the mechanism that breaks isotropic symmetry.** This connects directly to:

- The **quadrature model** — the two 90°-offset channels may represent L and T components. The base wave already encodes this duality
- **Complex sinusoids** — real = L, imaginary = T. This is why wave equations require complex numbers
- **Phase 1c (vector waves)** — full vector displacement naturally captures L and T as separate components
- The **"spinning in place" idea** — a standing wave that rotates phase instead of oscillating linearly requires two orthogonal components (L + T)

### ⚠️ Recommended Path Forward

Complete Step 2c testing (Options A/C/D) for validation — even if M2 suggests they'll fail, the 1D test has weaker cancellation and may reveal different behavior. Document results either way. Then evaluate whether to:

1. **Stay in 1D scalar** with the dual-channel (quadrature) model as a proxy for L/T — test if WCs can selectively disturb one channel
2. **Extend to 1D+1D (two-component displacement)** — add a transverse displacement track to the 1D engine, enabling L→T conversion at WCs without full 3D
3. **Move to Phase 1d (non-linear)** — variable λ(r) may break the sinc pattern independently of the L/T mechanism
4. **Move to Phase 1c (vector waves)** — full 3D vector displacement with spin, requires M4 engine

---

## Divergence & Curl

Divergence/curl/flux from M4 vector displacement; recovers charge sign from rotation direction

**Problem**: F = -∇(|ψ|²) uses scalar magnitude, which discards vector direction information. On-axis, vector reduces to scalar — no help for the standard test case.

**Opportunity**: vector displacement carries information beyond magnitude — ellipse rotation direction (handedness), divergence, curl, energy flux direction. These are **signed quantities** that could recover charge-phase information.

Force must be computed from a **different quantity** than |ψ|²:

- **Divergence** (∇·ψ): compression/rarefaction — scalar but signed
- **Curl** (∇×ψ): rotational displacement — vector, related to magnetic field
- **Energy flux** (ψ × ∂ψ/∂t or similar): directional energy flow
- **Per-component amplitude** (A_x, A_y, A_z separately): preserves directional structure

**One force, different directions**: F = -∇E is one force — electric (longitudinal), magnetic (transverse), gravitational (density deficit) are projections onto different components. Scalar collapses all directional information into magnitude — correct scaling (1/r²) but wrong direction (charge sign). See [04_magnetic_vector.md](04_magnetic_vector.md#why-scalar-is-insufficient-monopole--longitudinal-only) for full analysis.

**Connection to non-linear equations (Phase 1d)**: non-linear Ψ³ soliton, toroidal wave flows (r⁵), spin-as-vortex all require vector displacement. Phases 1d and 1c may converge.

Maybe this path needs 3D simulation, definitely 1c is not enough, but possibly 2D is not enough either. On the force unification concept, there is only ONE force, but at human scale (inertial frame / mass scale / frequencies that this scale of mass can experience), this single force appears at defined conditions that makes us perceive them as separate forces, so we named and describe them as so, but if they are a single 3D elliptical behavior that can be decomposed into 2 major amplitudes (90 degrees apart) and this elliptical form can be oriented in multiple orders in 3D space, this opens up the possibility. 2D simulation won't capture that.

Requires M4 vector wave method. Magnetic force arises from transverse wave interference. The elliptical displacement trajectories (6 phasor numbers per voxel) naturally encode both longitudinal (electric) and transverse (magnetic) components. Spin modeled as toroidal wave flow.

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

## On-Axis Vector Analysis

For two WCs on the x-axis, the displacement vectors at any point between them are anti-parallel (d̂₁ = +x̂, d̂₂ = -x̂). Only the x-component is nonzero, so |ψ_vec| reduces exactly to |ψ_scalar|. The vector energy |ψ_vec|² expands to:

```text
|ψ_vec|² = |ψ₁|² + |ψ₂|² + 2·|ψ₁|·|ψ₂|·cos(k·Δr + Δφ)·cos(θ_geo)
```

On-axis: cos(θ_geo) = -1 (constant), so F = -∇E is identical to scalar. No improvement for the standard on-axis test case. The vector advantage only manifests off-axis or with transverse components.

## Spin and Gravitational Shading

The energy for spinning comes FROM the longitudinal base wave. Spin reduces the longitudinal amplitude by converting it to transverse — this amplitude reduction may be the mechanism behind gravitational shading. The conversion ratio may be related to the **fine-structure constant α** (the coupling strength between electric/longitudinal and magnetic/transverse). If so, the gravitational force weakness (10⁻⁴² of electric) emerges from repeated application of this coupling.

## Permanent Magnets — Electron Alignment Persistence

Electrons are always spinning and generating magnetic fields, but in most materials these spin orientations are random and cancel out — no net magnetic field. When electrons move relative to protons (electric current), the electron-proton attraction creates resistance (reactive force) that aligns electron spins, producing the magnetic flux around a wire (electromagnet).

**Permanent magnets** (iron, neodymium, cobalt): these materials have atomic/molecular structures that allow electron spin alignment to persist. Once magnetized, their crystal structure delays the return to equilibrium — the electron alignment is "locked in" by the atomic lattice geometry. The stronger the locking mechanism, the more permanent the magnet. This is not a different force — it's the same spin-to-transverse-wave conversion, but with the alignment maintained by material structure rather than continuous electric current.

## Open Questions — Spin-Conversion Dynamics

- Does spin convert longitudinal to transverse at a fixed rate? Per spin cycle? Per unit time?
- Is the conversion ratio the fine-structure constant α (coupling strength)?
- Spin velocity → angular momentum → what determines it?
- Higher spin rate = stronger magnetic moment = more heat (thermal energy)?
- Magnetic moment = spatial distribution of transverse wave energy?
- Thermal energy = fluctuations in spin rate/axis (incoherent transverse emission)?
- What atomic lattice properties determine spin alignment persistence (permanent magnets)?

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

Soliton [Toroidal](00b_additional.md#soliton-geometry-two-domains-and-toroidal-flows) relationship.

Lafreniere [magnetic field lines](00b_additional.md#magnetic-fields-as-hyperboloid-wave-patterns) relationship.

## Other Possible Solutions

### Energy Flux (Radiation Pressure)

Force can also arise from **energy flux** — the directional flow of energy through the medium:

- **Energy density** (current: `F = -∇E`): energy *stored* per voxel. Force from the energy landscape shape
- **Energy flux** (`S = E · v_group`): energy *flowing* through a surface per unit time (W/m²). Has direction
- **Radiation pressure** (`P_rad = S/c`): force per unit area from wave momentum transfer (LaFreniere's mechanism)

For constant c: `∇P_rad = ∇E` — both approaches give the same force. They diverge when c varies spatially or waves are directional.

Energy flux could **naturally separate standing wave (near-field) from traveling wave (far-field) contributions** — standing waves have zero net flux, traveling waves have nonzero flux.

---

## Findings from Phase 1b Step 2c (hints for Phase 1c implementation)

Phase 1b tested L→T spin conversion (Option G) using the quadrature base wave's two phasor channels (P, Q) as a proxy for longitudinal and transverse components. The WC converts a fraction of L (P) into T (Q) with charge-dependent direction: `q = +1 → L→+T` (CW), `q = -1 → L→-T` (CCW).

**Result — the ONLY model that distinguishes charges**:

- Opposite charge: 12/24 attract, 12/24 repel (oscillates — sinc partially disrupted)
- Same charge: 24/24 unclear (both forces same direction — distinct from opposite)
- No other model (7 tested: 4 passive, 3 elastic) produced different behavior for opposite vs same charge

**Why the quadrature proxy is limited**: the two phasor channels P and Q are not truly independent — they combine into magnitude via `RMS = √(P² + Q²)/√2`. Converting energy between P and Q changes their individual values but the magnitude (and therefore energy) is still dominated by the sinc-modulated interference pattern. The L→T conversion shifts energy between channels, but `P² + Q²` doesn't fully decouple them.

**What Phase 1c needs**: true two-component displacement — independent L and T fields that contribute to energy **separately**:

```text
E_total = E_L + E_T = ρV(f·A_L)² + ρV(f·A_T)²
```

With independent components:

- L→T conversion at the WC reduces A_L and increases A_T locally
- The T component is NEW (not in the incoming field) → breaks isotropic cancellation (M2's core finding)
- The T component doesn't have the same spatial sinc structure as L → its contribution to the energy gradient is different
- Force direction depends on HOW the L/T balance changes with charge sign → charge-dependent force from wave physics, not imposed labels

**Implementation path**: extend the 1c engine with a second displacement field (ψ_T alongside ψ_L). Each has its own phasor (P_L, Q_L, P_T, Q_T). At WC positions, apply L→T conversion with charge-dependent direction. Energy computed from both: `E = E_L + E_T`. This is a minimal vector extension — 1c spatial domain but 2D displacement (longitudinal + transverse), sufficient to test the spin hypothesis without full 3D.

**Quadrature direction as charge/spin** — physical motivation in 3D? In 1c, the quadrature produces a traveling wave that flips direction with temporal offset sign (left ↔ right). In 3D, this maps to handedness — CW vs CCW rotation of the elliptical granule displacement. The two-direction limitation of 1c is a projection of the infinite orientational freedom in 3D. The quadrature model captured something real: the two channels (cos/sin, real/imaginary, L/T) encode a rotational degree of freedom that 1c collapses into direction. This is exactly what complex sinusoids do in quantum mechanics — the imaginary unit i encodes a 90° phase relationship that represents physical rotation in 3D space. The quadrature base wave may be the correct 1c representation of the 3D wave field, with the two channels being L and T components.

**Connection to M2 spin theory**: M2's `14_spin_theory.md` proposed exactly this mechanism and attempted implementation (`interact_wc_spinUP/DOWN`) but it "never worked correctly" in the 3D isotropic Laplacian field. The Phase 1b quadrature proxy test confirms the concept works (charge discrimination achieved) but the implementation needs true independent components. Phase 1c should revisit the M2 spin code with the improved understanding from Phase 1b testing.

**Convergence with Phase 1d**: Phases 1d and 1c may converge into one solution — variable λ(r) in the energy equation (1d) + vector displacement with L→T conversion (1c). The L→T spin creates charge-dependent asymmetry, while λ variation creates force from wavelength gradients. Both are needed: L→T alone still oscillates in scalar RMS, λ variation alone has no charge sensitivity. Together they could produce charge-dependent force from wavelength + mode gradients.
