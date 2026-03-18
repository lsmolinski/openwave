# 1D WAVE SANDBOX RESEARCH

✅ The 1D wave engine sandbox (`wave_engine_1D_v2.py`) is built and operational with:

- Weighted partial standing wave equation + phasor superposition
- Energy density (Joules) and force field (Newtons) panels
- Interactive controls: WC on/off toggles, separation slider, phase offset toggle
- Coulomb reference comparison with direction-match detection
- Force annotations at each WC position with attraction/repulsion labels

## ❌ Confirmed: Far-field oscillatory force issue reproduced in 1D

The same behavior seen in the 3D Taichi engine appears in the 1D sandbox:

- For both phase deltas (0° and 180°), force direction (attraction/repulsion) depends on WC separation distance
- Every λ/2 change in separation flips the force direction
- The wave interference between WCs creates a standing wave pattern in the phasor RMS that oscillates with separation
- At some separations the wave force matches Coulomb direction, at others it shows `WRONG DIRECTION`

This confirms the issue is in the wave equation / phasor physics, not in the 3D simulator implementation. The sinc oscillation in the out-wave spatial function is the root cause.

**Near-field opposite-phase issue**: Opposite-phase (opposite charge) WCs at near-field separations should show **monotonic attraction** — always pulling together until annihilation (at zero separation, waves cancel completely and both WCs cease to exist). Instead, the simulator shows the same oscillatory lock-in behavior as same-phase WCs — alternating attraction/repulsion every λ/2. This is incorrect: the lock-in oscillation should only occur for same-phase pairs (where it creates stable bonding positions). For opposite-phase pairs, the charge-phase difference (π) should ensure destructive interference always dominates between them, producing consistent attraction at all separations. The sinc node structure is overriding the charge-phase signal — another symptom of the same root cause.

## Analysis: What Wave Equation Reproduces This?

**PRIMARY TARGET:** [Lafreniere Attraction](02_equations.md#primary-reference-two-opposite-phase-wave-centers)

![alt text](images/wave_interference.gif)

This animation shows the complete wave interaction for two opposite-phase WCs:

- **Near-field** (~1λ around each WC): fixed standing wave rings — the particle structure itself. Standing wave forces dominate
- **Far-field** (beyond ~1λ): only traveling waves, 1/r amplitude falloff. This is the electrostatic (Coulomb) regime
- **Between WCs**: amplitude visibly reduced — destructive interference from opposing phases → lower energy zone → attraction. Counter-propagating traveling waves from each source meet at midpoint, forming a standing wave pattern there
- **Key for simulation**: standing wave region is sharply localized (~1λ, matches steep weight rolloff), far-field is cleanly traveling (no oscillatory artifacts), and the 1D envelope profile is what the phasor RMS should reproduce

## Critical Issues Summary — Electrostatic Force Validation

1. ❌ **Far-field oscillatory force (MAIN BLOCKER)** — Force direction flips every λ/2 of separation change. The sinc function sin(kr)/kr in the out-wave creates permanent nodes in the phasor RMS at all distances. These sinc nodes dominate over the charge-phase signal: source_offset (0 or π) should determine force direction via phasor P addition (same charge: cos(φ₁)=cos(φ₂) → P adds → larger amplitude) or cancellation (opposite charge: cos(φ₁)=-cos(φ₂) → P cancels → smaller amplitude). The spatial structure of where constructive/destructive interference occurs relative to each particle creates the gradient → force. But the sinc zeros override this charge-phase structure, making force direction depend on separation modulo λ instead of charge. Confirmed in both 3D and 1D engines
2. ❌ **Near-field opposite-phase monotonic attraction** — Opposite-charge WCs should always attract (to annihilation), not oscillate like same-charge lock-in. Currently both phase configurations show the same oscillatory behavior — same root cause as #1
3. ❌ **Gradient sampling radius** — Current gradient uses 1 grid unit. A real WC has spatial extent (~1λ). Wider gradient window (Gaussian smoothing with σ ~ λ) may smooth out the sinc oscillation and reveal the underlying Coulomb trend. **Tested**: Gaussian smoothing (σ = 0.25λ to 2λ) helps same-charge direction but hurts opposite-charge — smoothing removes the oscillation but also destroys the charge-dependent sign information encoded in the oscillation phase
4. ✅ **1/r² force law scaling (RESOLVED)** — The 1/r³ concern was based on single-source self-energy (E ∝ A² ∝ 1/r², gradient ∝ 1/r³). But force between two particles comes from the **interaction energy** E_int ∝ |Z₁|·|Z₂| ∝ 1/r, whose gradient is ∝ 1/r². Confirmed numerically: constant ratio to Coulomb across 2λ-10λ. See [04_challenges.md](04_challenges.md#-5-the-1r-force-law-resolved-in-theory)
5. ❌ **Dual-treatment boundary** — Near-field needs raw oscillatory phasor (for lock-in physics), far-field needs smoothed envelope (for Coulomb). The weight function transition boundary should serve double duty, but this is unimplemented and untested

Issues 1, 2, and 3 are connected — the oscillation is fundamental to the sinc spatial structure and cannot be removed by numerical techniques (smoothing, wider gradient) without also destroying the charge-phase signal that determines force direction.

## Test Configuration Notes

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

All task checklists are tracked in [01_plan.md](01_plan.md) ROADMAP (Phases 1 and 2).

---

### CANDIDATE DRIVERS OF ENERGY

## ✅ Force Computation: F = -∇E (Implemented)

Force is computed as the negative gradient of energy density:

```text
E(x) = ρ · V · (f · A(x))²      energy per grid point
F(x) = -∇E(x)                    force = negative energy gradient
```

This approach (rather than the expanded chain-rule form) ensures that when ρ(x), f(x), or λ(x) become spatially variable, ∇E automatically captures all contributions without changing the force computation logic. The same pattern applies across 1D sandbox, M3, and M4 engines.

## ✅ TESTED: Gradient Sampling Radius (Gaussian Smoothing)

Gaussian smoothing of the phasor RMS before computing ∇E was tested with σ = 0.25λ, 0.5λ, 1.0λ, 2.0λ. Results:

- **Same-charge direction**: improves with larger σ (8/17 → 13/17 at σ=2λ)
- **Opposite-charge direction**: degrades with larger σ (8/17 → 3/17 at σ=2λ)
- **Root cause**: smoothing extracts the 1/r envelope (always-negative gradient) but destroys the charge-phase signal encoded in the oscillation. The charge information lives in the *phase* of the sinc oscillation, not in the envelope — smoothing removes both

**Conclusion**: Gradient smoothing / wider sampling cannot resolve the oscillatory force because it removes the charge signal along with the oscillation.

## ✅ TESTED: Smooth Envelope Interaction (Imposed Charge Sign)

Compute force from the product of individual WC phasor magnitudes |Z₁|·|Z₂| (each smooth 1/kr) with charge sign imposed from source_offset difference:

- **17/17 correct direction** for both charge configurations (2λ to 10λ)
- **1/r² scaling confirmed** (constant ratio to Coulomb across all separations)
- **But**: charge sign is imposed by hand (`-1` for opposite, `+1` for same), not emergent from wave interference. This is equivalent to the old "analytical signed envelope" approach — Coulomb with extra steps, not force emergence

## POSSIBLE SOLUTION: Pressure/Velocity Gradient (90° Phase Shift)

Force may not respond to displacement amplitude directly — it may respond to **pressure** or **energy density**, which are related but phase-shifted quantities:

- **Displacement** (what we compute as ψ): the wave displaces the medium from equilibrium
- **Velocity** (∂ψ/∂t): 90° phase-shifted from displacement — granule velocity
- **Pressure/density**: related to velocity (or compression), also 90° phase-shifted from displacement

In fluid mechanics, force derives from pressure gradients, not displacement gradients. If the wave equation is fundamentally propagating **pressure or density variations** (not just displacement), then force may be **90° phase-shifted** from where the displacement gradient predicts it. This would change which points in space experience force and could resolve the oscillatory artifact — the pressure nodes are at different spatial positions than the displacement nodes.

The chain: `displacement → velocity (90° phase) → pressure/compression → energy → force`. If force follows pressure gradient rather than displacement gradient, the sinc zeros that create our oscillatory force artifact would shift by λ/4 — potentially breaking the exact cancellation pattern that makes charge-phase information unrecoverable.

The phasor already contains the 90° information (the Q component is sin(ωt), shifted from the P component cos(ωt)).

## POSSIBLE SOLUTION: Multi-Variable Energy Gradient

**Deeper question: Is amplitude the only variable that creates force?** The current approach treats `F = -∇E` where `E = ρV(fA)²`, assuming ρ, V, and f are constants — so only the amplitude gradient drives force. But the energy equation has **three variables** that can change spatially and create energy gradients:

- **A** (amplitude): currently the only variable — computed via phasor RMS
- **f / λ** (frequency / wavelength): if λ varies with position (time dynamics, variable-λ research), then ∇f contributes to ∇E even where ∇A = 0
- **ρ** (density): if medium density varies with position (Smoliński's buoyancy model), then ∇ρ contributes to ∇E independently of amplitude

The full energy gradient expands to: `∇E = ρV · [2fA² · ∇f + f² · 2A · ∇A] + V(fA)² · ∇ρ`

This means force can arise from amplitude curvature (current model), frequency curvature (time dynamics), OR density curvature (buoyancy) — or any combination. The fact that we can't reproduce clean electrostatic force from amplitude gradient alone may be a signal that **the force equation is incomplete** — we may need to include ∇f and/or ∇ρ terms even for the electrostatic case. Wave interference between WCs doesn't just change amplitude — it also changes the effective local frequency and medium density through the mechanisms described in the [Time Dynamics](07_time_dynamics.md) and [Smoliński sections](03_additional.md#smolińskis-contributions-bcc-lattice-geometric-framework).

Computing F = -∇E means these additional variables are automatically included when implemented — no force logic changes needed.

## POSSIBLE SOLUTION: Numerical Precision and Inertia Filtering

Two additional factors that may be masking or distorting the force signal:

- **Floating-point precision**: f32 arithmetic on very small force values (10⁻³ N scale at subatomic distances) may introduce numerical artifacts. The sinc function has near-zero values between nodes where precision loss is worst. Testing with f64 in the 1D sandbox (numpy default) should rule this out
- **Inertia as low-pass filter**: Real particles don't respond to 10²⁵ Hz oscillations — their mass acts as a low-pass filter, averaging out rapid force fluctuations. The force we compute should be the **time-averaged** force over at least one wave period, not the instantaneous gradient. The phasor RMS already time-averages amplitude — but if force depends on pressure (which is phase-shifted), the time-averaging may need to account for the phase relationship between displacement and pressure
