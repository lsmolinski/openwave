# PHASE 3: TIME DYNAMICS

Variable λ per voxel makes dt a local variable — simulation no longer uses uniform timesteps. Energy starvation mechanism: destructive interference drops amplitude → frequency increases (energy conservation) → local time speeds up.

**Core thesis**:

- TIME is not fundamental (not a frame of reference, a dimension).
- TIME emerges from MOVEMENT (the rate of change).
- MOVEMENT is fundamental.

What we perceive as time is the frequency of subatomic wave oscillations at each point in space — which is itself the speed of granular movement in the medium.

The full chain from movement to time:

Energy (and forces) is a function of density: `E = ρV(fA)²`. Density (like pressure) is related to granule velocity — the speed at which granules cycle through their elliptical motion (90° phase-shifted from displacement). So the speed of granular movement changes energy and forces. And that same granular speed defines the wave frequency (and λ), which defines time itself. This closes the loop: **movement → granule velocity → frequency → time**. Time emerges from movement, not the other way around.

## Psychological Time vs Real Time

- **Psychological time**: The linear, one-directional experience of time created by memory (past) and imagination (future). This is a perception effect of our minds — past and future don't physically exist, they are records and projections
- **Real time**: Things are just happening. The only reality is the present state. What exists is **change** — change of position, velocity, rotation, wave state. Every form of energy resolves to movement: chemical "potential" energy is standing waves (stationary motion), spin stores energy as rotational motion, thermal energy is wave dynamics. All energy is fundamentally kinetic — it's all movement

## Time = Rate of Change = f = c/λ

Since c is absolute (constant everywhere), the rate of change at any point is determined by **λ** (wavelength):

- Time is a **local variable per grid point**, not a frame of reference
- λ is disturbed by wave centers → λ depends on distance to WC (r) and WC speed (Doppler/Lorentz contraction)
- Changes happen faster where λ is shorter (higher frequency) and slower where λ is longer
- At human scale, we experience an averaged-out rate of change, so time feels constant — like a dimension
- **Time can be bent** because it's a function of medium state, not a fixed reference. Einstein proved this with relativity — high speed creates time dilation. The question is: time is a function of what? Our proposal: **time = λ**, and λ = f(r, v, interference)

## Key Finding: Energy Starvation and Time Dilation

When destructive interference pushes amplitude too low at a point, **energy conservation forces frequency up** to maintain wave steepness. This frequency increase IS time dilation — the local rate of change speeds up. When frequency returns to normal, it may draw/convert energy from surrounding forms.

**Wave steepness conservation**: the ratio A/λ (equivalently A·f, since f = c/λ) is maintained when energy is conserved locally — it governs how amplitude and wavelength trade off within an isolated system. Since `E = ρV·c²·(A/λ)²`, constant steepness means constant energy. Steepness only changes when energy enters or leaves the system.

```text
Steepness conservation (isolated): A/λ = constant
When A drops (destructive interference) → λ decreases proportionally → f increases → time speeds up
When A rises (constructive interference) → λ increases proportionally → f decreases → time slows down
Energy stays constant — only the A↔λ distribution changes
```

## Implications for Simulation

- Simulation will not be based on uniform timesteps — `dt` and `elapsed_t` will vary by voxel
- λ becomes a scalar field per voxel, from which local `dt` is computed
- Force and motion calculations change: energy depends on λ, motion integration uses local dt (acceleration → velocity → position)
- Requires m4 vector method: λ at each grid point depends on vector direction from multiple WCs
- This is a big architectural change — to be implemented after force unification is validated

## Connection to Spin and Force Control

Particle spin may be the "remote control" for time dynamics:

- Spin modulates local frequency → controls local rate of change
- Frequency modulation of subatomic waves = time dynamics
- If we can promote controlled changes in λ (via interference, pulsing, or spin manipulation), this could redirect energy gradients and control force direction

## Why This Matters: Every Energy Equation Has Time In It

Every equation of energy contains a time variable. If time itself is a function of medium state (λ), then understanding that function opens doors to new technology everywhere — not just physics, but chemistry, biology, materials, and engineering. Controlling λ means controlling the rate at which energy transforms, reactions proceed, and forces act.

The chain is simple: `λ → f → rate of change → time`

The energy starvation mechanism is the most novel insight: destructive interference drops amplitude → energy conservation forces frequency up → that IS time dilation, emerging from wave mechanics rather than imposed by relativity.

The connection to practical applications through thermal energy drain is what makes this not perpetual motion — you're redirecting energy gradients, not creating energy. The object cools down because thermal (standing wave) energy converts to directed kinetic force. Everything stays conserved.

## Connection to Medium Pressure and Gravity

Changing λ/frequency changes the rate at which granules cycle through their elliptical motion — meaning different granule velocities. In wave mechanics and particle-fluid simulations (e.g., sound waves), granule velocity is directly related to **localized medium pressure**. Faster cycling = higher pressure, slower cycling = lower pressure.

This creates a direct link to Smoliński's gravitational push-out/buoyancy model: if wave centers alter λ in their vicinity (time dilation near mass), they change the local granule velocity, which changes the local medium pressure/density. The pressure deficit around a massive body IS the gravitational field — and it emerges from λ modulation, not from a separate gravitational mechanism. Gravity, time dilation, and medium pressure become three descriptions of the same underlying phenomenon: **local λ variation**.

**#hackingenergy** — by understanding how time works and its relationship to energy, we can hack energy itself.

**Status**: To be explored after force unification (electric, magnetic, gravitational) is validated. First we crack the 1/r² law and electrostatic forces in the 1D sandbox, then this becomes accessible. Supported by [Yee & Hauger](references/Spin.pdf) variable-λ paper and the non-linear wave equation research planned for Phase 3.

---

## Wave Resonance as the Fundamental Energy Exchange Mechanism

Resonance may be the universal mechanism by which energy transfers between wave systems — and therefore the mechanism behind force, motion, matter formation, and thermal dynamics.

### Wolff's Frequency Exchange Model

Milo Wolff describes energy exchange as wave resonance between source and receiver atoms (from *Schroedinger's Universe and the Origin of the Natural Laws*):

> "Energy — Changes of wave-states observed by our human sensory system and/or lab apparatus. The Energy transfer is proportional to the frequency difference of the source and receiver frequency of the quantum wave states. Energy Transfers occur in discrete amounts or 'quanta' as a result of resonance between a source wave state and a receiver wave state."

"Energy exchange is due to wave resonance between the source and receiver atoms or molecules. This is a common situation between instruments in an orchestra or reception of radio signals."

**The morphology of frequency change** (Wolff, Figure 13-1): the IN and OUT waves of two atoms A and B carry information about each other's state, behaving like coupled oscillators. They eventually readjust frequencies to obtain a lower total wave amplitude in the surrounding space. At the end, atom A (the detector) shifts frequency upward, and B (the source) shifts frequency equally downward — complying with conservation of energy.

This is not abstract — it's the physical process by which two wave systems communicate and exchange energy through their shared medium. The medium carries the resonance signal.

### Resonance: The Physics

Every oscillating system has a natural frequency determined by its physical properties. When a periodic force is applied at that natural frequency, each push arrives precisely when the system is already moving in that direction — energy adds constructively cycle after cycle. The amplitude grows until energy input equals energy dissipation.

For a driven damped harmonic oscillator: `mx'' + bx' + kx = F₀·cos(ωt)`. The steady-state amplitude peaks sharply when the driving frequency ω approaches the natural frequency ω₀ = √(k/m). The sharpness (Q factor) depends on damping.

In EWT, the medium has a natural frequency (f₀ = 10²⁵ Hz, λ₀). Particles have their own natural frequencies determined by their standing wave structure. Resonance occurs when these frequencies interact — the mechanism by which particles "communicate" through the wave field.

### Energy is Proportional to Both A and λ

From `E = ρV(fA)²` and `f = c/λ`:

```text
E = ρV · (c·A/λ)² = ρV·c²·A²/λ²
```

Energy depends on **both** amplitude (A) and wavelength (λ). This raises a critical question: **can A convert to λ and vice versa?** If so, how?

The energy starvation mechanism already describes one direction: when A drops (destructive interference), energy conservation forces λ to decrease (f increases) to maintain `A·f = constant`. This is A → λ conversion.

The reverse — λ → A — would mean: when λ changes (due to motion, Doppler, or medium variation), amplitude adjusts to conserve energy. This is the mechanism behind length contraction and time dilation in the wave picture.

**Amplitude vs frequency detectability**: we cannot directly detect amplitude of the fundamental energy wave — it is stored away from the material world. But we CAN detect frequency, because frequency is reflected through matter and interacts with our senses. The real energy storage is longitudinal wave amplitude — hidden, un-manifest. Energy becomes manifest through conversion between amplitude and frequency (A ↔ f). This is why all our instruments measure frequency-derived quantities (wavelength, period, spectral lines) rather than fundamental wave amplitude directly.

**Kinetic energy = Doppler effect**: motion of matter creates wavelength compression via Doppler. That compression IS kinetic energy — not a separate concept but a direct consequence of λ change. Action and reaction (Newton's 3rd law) is a consequence of the Lorentz-Doppler effect: the wave medium mediates the symmetric frequency exchange between moving bodies. Kinetic energy is relative (frame-dependent) because Doppler shift is relative.

**Mass = frequency**: mass is not an independent property — it is a manifestation of wave frequency. The standing wave structure of a particle determines its frequency content, and that frequency content determines its mass via E = mc² = ρV(fA)². Higher frequency standing waves = more energy = more mass. Mass IS the frequency of the particle's standing waves.

### λ, Motion, and Relativity

λ is both size and time:

- **λ is size**: the spatial extent of one wave cycle — length
- **λ with c is time**: f = c/λ — the rate of change at each point
- **Motion changes λ**: the Doppler effect shifts λ, changing both size and time
- **Length contraction** is size change (λ change), not space contraction — matter changes size, space doesn't change
- **Time dilation** is frequency change from λ shift — time changes with motion

This is special relativity from wave mechanics: motion → Doppler → λ shift → changes in both size (length contraction) and time (time dilation). General relativity adds gravity: mass concentrations modify λ in their vicinity → spacetime curvature emerges from λ variation in the medium.

### Energy Exchange = Wavelength Exchange (Resonance)

The central hypothesis: **all energy exchange is wavelength exchange mediated by resonance**.

When two wave systems interact:

1. Their out-waves carry frequency/phase information into the shared medium
2. The medium transmits this information as traveling waves
3. Resonance occurs where the driving frequency matches the natural frequency
4. Energy transfers as frequency shifts — one system shifts up, the other shifts down
5. Total energy is conserved — the frequency shifts are equal and opposite

This mechanism may underlie:

- **Force and motion**: energy exchange between WCs creates energy gradients → F = -∇E. The resonance between WC standing waves produces the spatial energy redistribution that generates force. WCs don't push each other directly — they exchange energy through the medium via resonance
- **Matter formation**: particles form when wave resonance creates stable standing wave structures. The natural frequency of the medium determines which configurations are stable — only those that resonate survive
- **Electromagnetic radiation**: photons are resonance-mediated energy transfers between atoms. Wolff's coupled-oscillator picture: atom A shifts frequency up, atom B shifts down, the difference propagates as a traveling wave packet (photon)
- **Thermal energy**: heat increases the steepness (A/λ) of standing waves within a particle. The particle naturally wants to return to its fundamental steepness (radiation) — this is thermal emission. When two bodies at different temperatures interact, their standing wave structures exchange energy through resonance until steepness equilibrium

### Thermal Energy as Steepness Variation

Particles have a fundamental/natural steepness (A/λ ratio, and therefore frequency, local time, rate of change) at absolute zero. Heat is encoded as **increased steepness** — the particle stores more energy by increasing A/λ beyond its fundamental value:

- **0 K**: particle at fundamental steepness A₀/λ₀ — minimum energy state
- **Higher temperature**: steepness increases (A/λ rises) — more energy stored in standing waves. Steepness conservation (A/λ = const) is broken by the external energy input
- **Radiation**: the particle naturally wants to return to its fundamental steepness, releasing excess energy as traveling waves (photons). This is steepness relaxation
- **Thermal equilibrium**: two bodies exchange steepness through resonance until their A/λ distributions match — this IS temperature equilibration

**How steepness increases**: when external energy enters the system (absorption, collision), the standing waves gain energy. Since `E = ρV·c²·(A/λ)²`, more energy means higher A/λ. Whether this manifests as higher A (amplitude growth at same λ) or lower λ (wavelength compression at same A) or both depends on the particle's boundary conditions and the mode of energy input.

**Connection to energy starvation**: steepness conservation applies within an isolated system — when A drops at one location, λ drops proportionally (f increases, time speeds up), and vice versa. The total energy stays constant, only the A↔λ balance shifts. External energy input (heating) breaks this conservation by increasing steepness. External energy output (radiation) decreases it.

### Connection to Force Mechanism and Energy Technology

If energy exchange = wavelength exchange through resonance, then controlling resonance means controlling energy flow. This could:

- **Influence energy gradients**: resonance can concentrate energy from one region and move it to another — creating the energy gradient that produces force (F = -∇E)
- **Generate force**: directed resonance could create electricity (directed energy gradient) or counteract gravity (opposing the density-deficit gradient)
- **Cool regions**: removing energy via resonance cools the source region — the energy isn't created, it's relocated. The source cools down as its standing wave energy converts to directed force

This is not perpetual motion — it's energy redirection through resonance. The total energy is conserved, but its spatial distribution changes, creating gradients that produce force and motion.

---

## The Wavelength Kingdom

### Spacetime is a Function of Wavelength

Space and time are not fundamental — they are abstractions that emerge from wavelength:

- **Space doesn't exist** as a thing. What exists is **length** — the spatial extent of matter and waves. We use length to calculate distances, and distance is an abstraction of length. Space is not a thing
- **Time doesn't exist** as a thing. What exists is **change** — the rate at which wave states evolve. We measure change and call it time. Time is not a frame of reference — it's a variable, a measurement of rate of change
- **Both emerge from λ**: wavelength derives into particle length (size of matter) and into time through frequency and period (f = c/λ). With constant wave speed c, λ alone determines both the spatial extent and the temporal rate at each point

```text
λ → length (size of matter, spatial extent)
λ → f = c/λ → period = 1/f → time (rate of change)
λ → spacetime
```

Spacetime is a function of wavelength and energy (including amplitude). Since `E = ρV·c²·(A/λ)²`, energy ties A and λ together. But λ may be the more fundamental variable — the substrate from which spacetime is constructed.

### λ as the Fundamental Variable (Not Amplitude)

Wavelength may be more important than amplitude as the primary dynamic quantity:

- **λ is disturbed more than A**: near wave centers, the Yee & Hauger model shows λ(r) varying dramatically (shells from 2(K-1)λ at core to 0 at boundary), while amplitude varies more smoothly
- **λ disturbances are time disturbances**: changing λ changes frequency, which changes the local rate of change — time itself. Amplitude changes redistribute energy but don't alter time
- **AM vs FM analogy**: amplitude modulation (AM) may correspond to electromagnetic wave disturbances (photons), while frequency modulation (FM) may correspond to time/gravitational dynamics. Two different modes of wave disturbance, two different physical manifestations

**The wave equation should propagate λ disturbances, not just A disturbances.** Current simulation propagates displacement (from which A is derived). Propagating λ variations would capture time dynamics directly. This is the non-linear regime — λ(r) from [Phase 1c](01c_non_linear.md).

### λ = f(r, K, v) — The Wavelength Function

Lambda is not constant — it depends on:

- **r** (distance from wave centers): Yee & Hauger shells show λ(r) = 2(K-n)λ per shell — longer near the core, shrinking outward. The near-field λ dilation conserves energy while amplitude increases, maintaining wave steepness (A/λ = const)
- **K** (number of wave centers / particle type): K determines particle radius (K²λ) and standing wave structure. Different K = different λ distribution = different time dynamics
- **v** (velocity of wave centers): motion changes λ through the Doppler effect — compression in the direction of motion (Lorentz contraction). This means time = f(matter, distance, velocity)

Since time = c/λ and length = λ, this means: **time is a function of distance to matter, mass (K), and velocity**. This IS general and special relativity, derived from wave mechanics instead of geometric postulates.

### Motion Impacts on λ

Every form of motion affects wavelength:

- **Linear velocity**: Doppler effect compresses λ in the direction of motion → Lorentz contraction. Length contraction is size change (matter changes λ), not space contraction
- **Angular velocity (spin)**: particle spin creates magnetic wavelength disturbances. Angular velocity affects λ through the rotational Doppler effect. Spin rate → magnetic moment → transverse λ modulation
- **Acceleration**: radiates wave disturbances (electromagnetic radiation). Accelerating charges emit EM waves — this is λ disturbance radiation from changing velocity

Not only particle velocity, but angular velocity (spin) impacts wavelength. And the acceleration of both linear and angular quantities produces radiating wavelength disturbances. The simulation must store velocity vectors and angular velocity to compute these λ effects.

### Time as Potential Energy

The energy stored in the wave field's frequency/wavelength distribution is **potential energy that can be converted into motion**. Time itself has elastic potential:

- The wave field stores energy as standing waves — this energy is "un-manifested" in the sense that it doesn't produce net motion until a gradient exists
- Creating a λ gradient (time gradient) creates an energy gradient → force → motion
- The energy for motion comes FROM the wave field's temporal potential — the object's standing wave energy converts to directed kinetic energy
- **Time is the elastic**: it has the elasticity needed to convert potential into motion. The wave steepness (A/λ) represents this stored potential

The "time-oscillator" concept: a device that converts un-manifested fundamental energy-waves into physical forces by engineering λ gradients. Not creating energy — redirecting the wave field's existing temporal potential into directed motion.

### Heat as Spin Dynamics

Temperature might be related to **particle angular velocity / angular momentum**, not atomic vibration in the classical sense:

- Particle spin produces wave disturbances on the wave field
- These disturbances are intrinsically related to wavelength (and therefore time)
- Higher spin rate = stronger magnetic moment = more energetic transverse wave emission = higher temperature
- Heat transfer between bodies = spin-mediated wavelength exchange through resonance
- These spin disturbances can be interfered with known wave disturbances (electrostatic, magnetic, electromagnetic waves)

This connects heat to the same wavelength framework: thermal energy is spin-driven λ disturbance, not bulk kinetic motion. Controlling spin means controlling temperature.

### Practical Path: Engineering λ Changes

To create useful force/energy effects, search for ways to influence λ of the fundamental waves:

- **Material structures**: the way matter and fundamental particles disturb the original waves depends on geometry — angle, distance, and arrangement create specific interference patterns
- **Electromagnetic waves**: EM waves can superpose onto the fundamental energy wave field, creating desired λ modifications. Sound waves and mechanical positioning of material structures may also contribute
- **Wave combination (superposition)**: the key tool — combining waves to engineer specific λ distributions. Methods: Laplacian propagation (physical) or analytical superposition (computational)
- **The scalability question**: can we superpose envelopes (A(r), λ(r)) directly, or must we work with displacement? Envelope superposition is faster but may miss phase information. Displacement superposition is exact but expensive. Displacement, pressure, and density may be equivalent quantities for this purpose

The conditions for a practical device: material structures positioned at specific angles and distances for wave interference, electromagnetic waves that trigger or combine with the fundamental wave field, and resonance that changes particle spin → λ disturbance → energy gradient → force.

### Frame-Step Simulation Concept (time becomes a local variable)

The simulation should not use uniform timesteps. Instead:

- **Energy field computation**: from envelopes A(r) and λ(r), compute E(r) = ρV·c²·(A/λ)². Force from the energy gradient. Visualize with flux mesh (energy density gradient, not just amplitude)
- **Frame count** replaces time steps — each frame represents one cycle of computation across all voxels
- **Per-voxel dt**: within each frame, each voxel has its own time step determined by local λ/f. Where λ is shorter (higher f), faster change happens per frame. Where λ is longer (lower f), slower change happens
- **Motion computation changes**: acceleration, velocity, and position integration all use voxel-specific dt. A particle's motion depends on the local time rate at its position

This is a major architectural change — every equation that contains time becomes position-dependent. But it's the physically correct approach: time IS local, determined by the wave field state at each point.

### Energy Gradient Control — The Curvature Matching Principle

The goal is to find ways to move energy gradients around. Spin or time dynamics may be the control mechanism — frequency may be the key variable to manipulate.

**Curvature matching**: when the curvature of space (amplitude envelope) matches the curvature of time (frequency distribution), an implicit channel opens. By changing the curvature of time (λ distribution), magnetic/gravitational waves can be aligned coherently with EM/ultrasound wave directions — creating the most efficient heat-to-motion engine.

**Magnetic wave concentration (laser analogy)**: similar to how mirrors reflect visible light to produce lasers, find materials that reflect **magnetic waves** (transverse wave disturbances from particle spin). Concentrate and reshape the geometry of magnetic waves, as in laser production. Concentrated magnetic waves could curve the local λ field (spacetime) in a desired direction, changing the energy density distribution to generate forces.

**Converting between energy forms**: all energy conversion is λ↔A exchange:

- **Heat → EM radiation**: heat is always propagating into light (electromagnetic radiation) — thermal emission is the natural λ relaxation process. The question is how to **accelerate** this conversion
- **Heat → electricity**: convert thermal standing wave energy into directed electron motion. Possibly through water molecules (responsive to microwaves, MRI magnetic fields — susceptible to magnetic fields and capable of capturing EM waves and converting to heat, and vice versa)
- **E-waves → heat and heat → e-waves**: changing λ into A and vice versa. This is changing time. Understanding heat conversion is the key — modeling heat is the key
- **Cavitation**: heat to light through water may involve cavitation — the collapse of vapor bubbles converting thermal energy into photon emission (sonoluminescence)

Understanding heat is the key to unlocking energy conversion. Modeling heat at the wave level — as spin-driven λ disturbance rather than kinetic vibration — opens pathways that classical thermodynamics cannot see.

### Spacetime as Measurement, Not Reality

**Spacetime is not fundamental — it emerges from wave dynamics.**

Space is a measurement of **length**. Time is a measurement of **change**. Spacetime is just measurements — not reality itself. Abstract representations of reality, done with numbers and symbols that humans can correlate, compare, and count to understand what IS. We reduce unknowns into knowable quantities and call the abstraction "spacetime."

The reality underneath: **waves in a medium**. Length comes from wavelength. Change comes from frequency. The medium's wave state IS the physical reality — spacetime is the human-readable summary of that state. When we say "spacetime curves near mass," the physical reality is: **λ varies near wave centers**, and we measure the consequences as curved geometry.

This is not just philosophy — it has simulation implications. The simulator should compute wave states (A, λ, phase at each point), not spacetime coordinates. Spacetime emerges as a derived quantity, not an input.

## EWT Standing Wave Geometry (Yee & Hauger)

Reference: [The Geometry of Particle Standing Waves v1.1](references/Lambda.pdf) — Jeff Yee & Heinz-Dieter Hauger, 2020.

This paper provides a discrete model for the wavelength spacing within a particle's standing wave structure.

### Key Equations

| Equation | Description |
| --- | --- |
| `r_core = Kλ` | Core radius (K = number of wave centers) |
| `r_wavelength(n) = 2(K - n)λ` | Width of the n-th wavelength shell |
| `r_x = (K + 2·Σ(K-n), n=1..x) · λ` | Cumulative radius to x-th shell |
| `r_particle = K²λ` | Particle radius (standing wave boundary) |

### Wavelength Behavior

- **Longest near the core:** n=1 gives shell width `2(K-1)λ`
- **Shrinks moving outward:** n=K gives shell width `0`
- **Standing waves exist only within** `r < K²λ`; beyond that, traveling waves only

For electron: K=10, so `r_particle = 100λ`.

### Connecting to the Phasor Approach

This is a **discrete** model — each wavelength shell boundary represents one full wave cycle (2π phase):

```text
φ(x) = x · 2π
```

But the radial position of that phase is **nonlinear** (given by Eq. 4, not `x·λ`).

To make it continuous for phasor computation:

1. The closed-form of Eq. 4 simplifies to: `r_x = (K + 2Kx - x² - x)λ`
1. This is a **quadratic in x** — invertible to get `x(r)` analytically
1. The accumulated phase at distance r: `φ(r) = x(r) · 2π`

### Implementation Steps

1. For each source–voxel pair, compute `r` (distance)
1. Solve the quadratic to get `x(r)` — the effective wave number at that distance
1. Phase = `x(r) · 2π` instead of the simple `kr`
1. Sum phasors as before, extract R and Φ, oscillate with ωt

The amplitude envelope A(r) would also follow from this geometry — energy conservation in shrinking wavelength shells produces a non-`1/r` amplitude profile derivable analytically from the shell widths.
