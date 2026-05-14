# ✅ PHASE 1: EXPLORING WAVE EQUATIONS

## Table of Contents

- [The Problem](#-the-problem)
- [Target Behavior (LaFreniere Reference)](#target-behavior-lafreniere-reference)
- [Critical Issues Summary](#critical-issues-summary)
- [Tested and Ruled Out (9/9)](#tested-and-ruled-out-99)
  - [Gradient Sampling / Gaussian Smoothing](#-gradient-sampling-radius-gaussian-smoothing)
  - [Smooth Envelope (Imposed Charge Sign)](#-smooth-envelope-interaction-imposed-charge-sign)
  - [Numerical Precision](#-numerical-precision)
  - [Inertia Filtering](#-inertia-filtering)
  - [Pressure/Velocity Gradient (90° Phase Shift)](#-pressurevelocity-gradient-90-phase-shift)
  - [Standing vs Traveling Wave Decomposition](#-standing-vs-traveling-wave-decomposition)
  - [Alternative Wave Equations (All 5 Forms)](#-alternative-wave-equations-all-5-forms)
- Major Paths
  - [PHASE 1a: Signed Disturbance (forced charge sign)](1a_signed.md#phase-1a-signed-disturbance-forced-charge-sign)
  - [PHASE 1b: Base Wave + WC Energy Redistribution](1b_base_disturbance.md#phase-1b-base-wave--wc-energy-redistribution)
  - [PHASE 1c: Non-Linear Wave Equations (1D)](01c_non_linear.md#phase-1c-non-linear-wave-equations)
  - [PHASE 1d: Vector Wave Force (M4 displacement direction)](01d_vector_wave.md#phase-1d-vector-wave-force)
- [Summary](#summary)
- [Other Possible Solutions](01d_vector_wave.md#other-possible-solutions)
- [1D Sandbox & Test Configuration](#1d-sandbox--test-configuration)

---

## ❌ The Problem

The core blocker in the force implementation is that the sinc function (`sin(kr)/kr`) in the out-wave introduces permanent oscillatory nodes in the phasor RMS at all distances, causing force direction to flip every λ/2 — even where only smooth 1/r decay is expected — confirmed in both 3D and 1D engines.

After exhausting all 10 linear wave approaches, it's been confirmed that the oscillatory force is a fundamental consequence of coherent wave interference — no linear operation on the superposed field can eliminate it while preserving charge-dependent direction — and non-wave signed approaches (Phase 1a) merely replicate Coulomb's law without emergent charge.

The next step is validating the energy field by plotting E(x) along the inter-particle axis at various separations to identify interference locations and verify that gradient direction and magnitude match the expected 1/r² Coulomb reference.

Force is computed as F = -∇E where E(x) = ρ·V·(f·A(x))², calculated directly from ∇E rather than chain-rule expansion so that future spatially-varying ρ(x), f(x), and λ(x) are automatically captured.

---

## Target Behavior (LaFreniere Reference)

**PRIMARY TARGET:** [Lafreniere Attraction](0a_equations.md#primary-reference-two-opposite-phase-wave-centers)

![alt text](images/wave_interference.gif)

This animation shows the complete wave interaction for two opposite-phase WCs:

- **Near-field** (~1λ around each WC): fixed standing wave rings — the particle structure itself. Standing wave forces dominate
- **Far-field** (beyond ~1λ): only traveling waves, 1/r amplitude falloff. This is the electrostatic (Coulomb) regime
- **Between WCs**: amplitude visibly reduced — destructive interference from opposing phases → lower energy zone → attraction. Counter-propagating traveling waves from each source meet at midpoint, forming a standing wave pattern there
- **Key for simulation**: standing wave region is sharply localized (~1λ, matches steep weight rolloff), far-field is cleanly traveling (no oscillatory artifacts), and the 1D envelope profile is what the phasor RMS should reproduce                   |

## Critical Issues Summary

1. ❌ **Far-field oscillatory force (MAIN BLOCKER)** — Force direction flips every λ/2 of separation change. The sinc function sin(kr)/kr creates permanent nodes in the phasor RMS. These sinc nodes dominate over the charge-phase signal, making force direction depend on separation modulo λ instead of charge. Confirmed in both 3D and 1D engines
2. ❌ **Near-field opposite-phase monotonic attraction** — Opposite-charge WCs should always attract (to annihilation), not oscillate like same-charge lock-in. Same root cause as #1
3. ✅ **1/r² force law scaling (RESOLVED)** — Force between two particles comes from the **interaction energy** E_int ∝ |Z₁|·|Z₂| ∝ 1/r, whose gradient is ∝ 1/r². Confirmed numerically. See [0c_challenges.md](0c_challenges.md#-5-the-1r-force-law-resolved-in-theory)
4. ❌ **Dual-treatment boundary** — Near-field needs raw oscillatory phasor (for lock-in physics), far-field needs smoothed envelope (for Coulomb). Unimplemented

---

## Tested and Ruled Out (9/9)

### ✅ Gradient Sampling Radius (Gaussian Smoothing)

Gaussian smoothing of the phasor RMS before computing ∇E was tested with σ = 0.25λ, 0.5λ, 1.0λ, 2.0λ. Results:

- **Same-charge direction**: improves with larger σ (8/17 → 13/17 at σ=2λ)
- **Opposite-charge direction**: degrades with larger σ (8/17 → 3/17 at σ=2λ)
- **Root cause**: smoothing extracts the 1/r envelope (always-negative gradient) but destroys the charge-phase signal encoded in the oscillation

**Conclusion**: Gradient smoothing / wider sampling cannot resolve the oscillatory force because it removes the charge signal along with the oscillation.

### ✅ Smooth Envelope Interaction (Imposed Charge Sign)

Compute force from the product of individual WC phasor magnitudes |Z₁|·|Z₂| (each smooth 1/kr) with charge sign imposed from source_offset difference:

- **17/17 correct direction** for both charge configurations (2λ to 10λ)
- **1/r² scaling confirmed** (constant ratio to Coulomb across all separations)
- **But**: charge sign is imposed by hand (`-1` for opposite, `+1` for same), not emergent from wave interference — Coulomb with extra steps, not force emergence

### ✅ Numerical Precision

The 1D sandbox uses numpy f64 (15 decimal digits). All physical constants are in simulation-friendly units (am, rs, qg) with magnitudes in O(0.01)–O(10). The sinc zeros are exact zeros of sin(kr), not near-cancellation artifacts. **Conclusion**: the oscillatory force is a real mathematical feature, not a floating-point artifact.

### ✅ Inertia Filtering

The phasor RMS already IS the time-averaged amplitude — it gives the exact analytical result without simulating oscillation cycles. The oscillatory problem is in the *spatial* sinc structure, not in temporal high-frequency artifacts.

### ✅ Pressure/Velocity Gradient (90° Phase Shift)

A 90° phase shift does not resolve the oscillatory force:

- **Velocity (∂ψ/∂t)**: velocity RMS = ω × displacement RMS. Since ω is constant, ∇(ω·RMS) = ω·∇(RMS) — identical gradient direction, zero benefit
- **Pressure (∝ -∂ψ/∂x)**: gives d/dr[sin(kr)/kr] instead of sin(kr)/kr — different zeros, but still oscillates with the same λ/2 period, just shifted by λ/4

No phase rotation or derivative of a periodic function removes its periodicity.

**Note**: granule velocity is still physically significant — it relates directly to medium density ρ (see [Time Dynamics](06_time_dynamics.md): faster cycling = higher pressure/density). This connection is explored in the Non-Linear Wave Equations section (Phase 1c).

### ✅ Standing vs Traveling Wave Decomposition

**Single WC analysis**: each component individually is smooth:

- **Standing** (in-wave): RMS = A·w(r)/kr — smooth, dies off with weight function
- **Traveling** (out-wave): RMS = A/kr — smooth 1/r everywhere

**Two WC analysis**: even with traveling-wave-only (smooth 1/kr per source), coherent superposition creates oscillatory interference:

```text
E_interaction ∝ cos(k(r₁ - r₂) - (φ₁ - φ₂)) / (kr₁ · kr₂)
```

**The oscillation is intrinsic to coherent wave interference**, not to standing vs traveling character.

### ✅ Alternative Wave Equations (All 5 Forms)

| Equation | Spatial zeros | Zero spacing | Force direction flips? |
| --- | --- | --- | --- |
| #1 Wolff | sin(kr)/r | λ/2 | Yes |
| #2 LaFreniere-Marcotte | sin(kr)/(kr) + (1-cos(kr))/(kr) | λ | Yes |
| #3 Phase-warped Marcotte | sin(x_c)/x_c | ~λ (warped near core) | Yes |
| #4 Combined Wolff-LF | sin(kr)/r + (1-cos(kr))/r | λ/2 | Yes |
| #5 Weighted (current) | (w+1)sin(kr)/(kr) + (w-1)cos(kr)/(kr) | depends on w(r) | Yes |

**Confirms that the oscillation is intrinsic to coherent wave interference** regardless of spatial function.

---

## Summary

Three remaining paths, all connected:

1. **Base wave + WC energy redistribution** ([PHASE 1b: Base Wave + WC Energy Redistribution](1b_base_disturbance.md#phase-1b-base-wave--wc-energy-redistribution)) — model the actual base wave and how WCs redistribute its energy. Central open question: how WC phase determines the drainage pattern
2. **Non-linear wave equations** ([PHASE 1c: Non-Linear Wave Equations (1D)](01c_non_linear.md#phase-1c-non-linear-wave-equations)) — variable λ(r), ρ(x), Ψ³; breaks sinc periodicity while keeping genuine wave interference
3. **Vector wave force** ([PHASE 1d: Vector Wave Force (M4 displacement direction)](01d_vector_wave.md#phase-1d-vector-wave-force)) — divergence/curl/flux from M4 vector displacement; recovers charge sign from rotation direction

---

## 1D Sandbox & Test Configuration

### Sandbox Status

✅ The 1D wave engine sandbox (`wave_engine_1D_v2.py`) is built and operational with:

- Weighted partial standing wave equation + phasor superposition
- Energy density (Joules) and force field (Newtons) panels
- Interactive controls: WC on/off toggles, separation slider, phase offset toggle
- Coulomb reference comparison with direction-match detection
- Force annotations at each WC position with attraction/repulsion labels

### Confirmed Issue

The same behavior seen in the 3D Taichi engine appears in the 1D sandbox:

- For both phase deltas (0° and 180°), force direction (attraction/repulsion) depends on WC separation distance
- Every λ/2 change in separation flips the force direction
- The wave interference between WCs creates a standing wave pattern in the phasor RMS that oscillates with separation
- At some separations the wave force matches Coulomb direction, at others it shows `WRONG DIRECTION`

This confirms the issue is in the wave equation / phasor physics, not in the 3D simulator implementation.

**Near-field opposite-phase issue**: Opposite-phase WCs at near-field separations should show **monotonic attraction** — always pulling together until annihilation. Instead, the simulator shows the same oscillatory lock-in behavior as same-phase WCs. The sinc node structure is overriding the charge-phase signal.

### Test Configuration

Standard test setup for reproducing LaFreniere reference animations:

- 2 wave centers separated by ~5-10λ
- Same charge test: both source_offset = π (or both = 0)
- Opposite charge test: source_offset = 0 and π
- Weighted partial standing wave with transition = 1.25λ
- Visualize: displacement + phasor RMS overlay, energy density, force field
- Compare 1D profiles against animation cross-sections
- Compute force from phasor RMS gradient at each WC position

**Two-regime tests:**

- Near-field (Test A): 2 particles within 1-2λ — observe oscillatory lock-in, measure stability, test Verlet integrator and f64 precision
- Far-field (Test B): 2 particles at 5λ, 10λ, 15λ, 20λ — measure force vs distance, compare against Coulomb 1/r²
