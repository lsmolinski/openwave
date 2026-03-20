# PHASE 6: Time Dynamics (M4 or new method)

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

When destructive interference pushes amplitude too low at a point, **energy conservation forces frequency up** to maintain wave steepness (amplitude × frequency must conserve energy). This frequency increase IS time dilation — the local rate of change speeds up. When frequency returns to normal, it may draw/convert energy from surrounding forms.

This connects amplitude, frequency, wavelength, and energy conservation into a single dynamic:

```text
Wave steepness conservation: A · f = constant (locally)
When A drops (destructive interference) → f increases → λ decreases → time speeds up
When A rises (constructive interference) → f decreases → λ increases → time slows down
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

## EWT Standing Wave Geometry (Yee & Hauger)

Reference: [The Geometry of Particle Standing Waves v1.1](references/Spin.pdf) — Jeff Yee & Heinz-Dieter Hauger, 2020.

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
