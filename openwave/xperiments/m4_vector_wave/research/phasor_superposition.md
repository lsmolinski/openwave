# Phasor Superposition for Analytical Wave Field Computation

Research notes on replacing per-timestep numerical wave superposition with analytical phasor precomputation.

## Problem

Given multiple wave sources emitting waves of the form:

```text
ψ = A(r) · cos(ωt ± kr + φ + π) · direction
```

Do we need to simulate the superposition of all waves from all sources at every timestep to get `displacement_am[i, j, k]`? Or can we derive a combined equation and just oscillate the temporal term?

## Phasor Superposition (Constant λ)

**Yes — analytically solvable**, as long as all sources share the same frequency ω.

### Key Insight

Each wave from source *n* at a grid point:

```text
ψₙ = Aₙ(rₙ) · cos(ωt - krₙ + φₙ)
```

Since they all share ω, the sum of N cosines at the same frequency is **always** a single cosine at that frequency:

```text
ψ_total = R · cos(ωt + Φ)
```

where **R** (resultant amplitude) and **Φ** (resultant phase) are computed once per grid point via phasor addition.

### Method

Represent each wave as a complex phasor (dropping the shared ωt):

```text
Pₙ = Aₙ(rₙ) · e^{i(-krₙ + φₙ)} · d̂ₙ
```

where `d̂ₙ` is the direction vector. Sum them:

```text
P_total = Σ Pₙ
```

Then at each grid point, extract the resultant vector phasor **per xyz component**:

- `R_x, R_y, R_z = |P_total_x|, |P_total_y|, |P_total_z|`
- `Φ_x, Φ_y, Φ_z = arg(P_total_x), arg(P_total_y), arg(P_total_z)`

### Two-Phase Computation

1. **Precomputation (once):** For each grid point, sum the phasors from all sources to get `R[i,j,k]` and `Φ[i,j,k]` (vector quantities, one per spatial component).
1. **Time evolution (cheap):** Evaluate per component — no per-source loop needed at runtime:

```text
displacement[i,j,k].x = R_x · cos(ωt + Φ_x)
displacement[i,j,k].y = R_y · cos(ωt + Φ_y)
displacement[i,j,k].z = R_z · cos(ωt + Φ_z)
```

The scalar amplitude at that point: `√(R_x² + R_y² + R_z²)`

### Caveats

- All sources must have the **same ω** (monochromatic). If frequencies differ, one phasor sum per distinct frequency is needed.
- Amplitudes A(r) must not change in time (no source motion or amplitude modulation).
- The wave equation must be **linear** (superposition principle holds).

## Analytical Amplitude (Replaces EMA-RMS Tracking)

The phasor approach gives amplitude **directly** — no RMS tracking needed.

The magnitude `|P_total|` *is* the steady-state amplitude at that point. It's exact, not an approximation converging over time like EMA-RMS.

Since the wave is `R · cos(ωt + Φ)`:

- Displacement oscillates between `-R` and `+R`
- RMS = `R / √2`

Both are known analytically from precomputation. **No temporal tracking required.**

This replaces both:

- The per-timestep multi-source superposition loop
- The EMA-RMS amplitude tracker (`amp_local_rms_am`)

## Extension: Variable Wavelength λ(r)

### Hypothesis

Wave centers may have longer wavelengths at the core (pushing amplitude up to conserve energy/steepness), with λ getting shorter as r increases until stabilizing at some distance.

### Phasor Approach with Variable λ

The phasor method **still works**. The only change is how phase is computed.

With constant λ, phase from source to point = `kr` (linear in distance).

With λ(r), since `k(r) = 2π/λ(r)`, the phase becomes an **integrated phase** (WKB / eikonal approach):

```text
φ(r) = ∫₀ʳ k(r') dr'
```

The phasor for source *n* becomes:

```text
Pₙ = Aₙ(rₙ) · e^{i(-φ(rₙ) + φₙ)} · d̂ₙ
```

Everything else is identical — sum phasors, extract R and Φ per component, oscillate with ωt.

### Considerations for Variable λ

1. **The phase integral:** For each source–voxel pair, compute `∫₀ʳ k(r') dr'`. If λ(r) has a closed-form expression, this integral may have an analytical solution. Otherwise, precompute it numerically once.
1. **Amplitude A(r):** If λ shortens with r to conserve energy/steepness, the amplitude envelope will differ from simple `1/r`. Derive from energy conservation — typically `A(r) ∝ 1/(r · √(k(r)))` (WKB amplitude correction).
1. **Still monochromatic:** As long as the temporal frequency ω remains the same everywhere (λ changing with r means the *spatial* structure changes, not the oscillation rate), phasor summation remains valid.

## EWT Standing Wave Geometry (Yee & Hauger)

Reference: *The Geometry of Particle Standing Waves v1.1* — Jeff Yee & Heinz-Dieter Hauger, 2020.

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

## Impact of Moving Wave Centers (Force & Motion)

If wave centers move due to force/motion computation, the phasor precomputation (R and Φ at each grid point) becomes invalid after each move and must be recomputed.

The recomputation is the same `O(grid_points × sources)` loop as the current per-timestep superposition. **If sources move every frame, there is no performance gain** over the current approach.

### When Phasors Win vs. Don't

| Scenario | Phasor benefit |
| --- | --- |
| Static sources, many render frames | Huge — precompute once, oscillate cheaply |
| Sources move every N frames | Moderate — amortize phasor cost over N frames |
| Sources move every frame | **None** — same cost as current approach |

### Multi-Rate Stepping (Potential Optimization)

Force/motion typically operates on a **slower timescale** than wave oscillation. If decoupled:

- **Outer loop (slow):** Compute forces, move sources, recompute phasor field
- **Inner loop (fast):** Oscillate `cos(ωt + Φ)` for visualization / sub-stepping

Phasors pay off during the inner loop. The benefit depends on whether the sim needs multiple wave oscillation steps per force update — if force and wave evaluation are 1:1, there's no gain.

### Benefits Regardless of Source Motion

Even with sources moving every frame, phasors still provide:

- **Exact amplitude** `R[i,j,k]` for free — no EMA-RMS convergence lag
- **Cleaner separation** of spatial structure from temporal oscillation
- **Direct inspectability** of amplitude/phase fields between force steps

## Elliptical Motion Geometry

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
