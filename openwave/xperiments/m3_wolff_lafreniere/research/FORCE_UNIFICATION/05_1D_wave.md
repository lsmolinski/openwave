# 1D WAVE SANDBOX RESEARCH

The 1D wave engine sandbox (`wave_engine_1D_v2.py`) is built and operational with:

- Weighted partial standing wave equation + phasor superposition
- Energy density (Joules) and force field (Newtons) panels
- Interactive controls: WC on/off toggles, separation slider, phase offset toggle
- Coulomb reference comparison with direction-match detection
- Force annotations at each WC position with attraction/repulsion labels

## Confirmed

Far-field oscillatory force issue reproduced in 1D.** The same behavior seen in the 3D Taichi engine appears in the 1D sandbox:

- For both phase deltas (0° and 180°), force direction (attraction/repulsion) depends on WC separation distance
- Every λ/2 change in separation flips the force direction
- The wave interference between WCs creates a standing wave pattern in the phasor RMS that oscillates with separation
- At some separations the wave force matches Coulomb direction, at others it shows `WRONG DIRECTION`

This confirms the issue is in the wave equation / phasor physics, not in the 3D simulator implementation. The sinc oscillation in the out-wave spatial function is the root cause.

**Near-field opposite-phase issue**: Opposite-phase (opposite charge) WCs at near-field separations should show **monotonic attraction** — always pulling together until annihilation (at zero separation, waves cancel completely and both WCs cease to exist). Instead, the simulator shows the same oscillatory lock-in behavior as same-phase WCs — alternating attraction/repulsion every λ/2. This is incorrect: the lock-in oscillation should only occur for same-phase pairs (where it creates stable bonding positions). For opposite-phase pairs, the charge-phase difference (π) should ensure destructive interference always dominates between them, producing consistent attraction at all separations. The sinc node structure is overriding the charge-phase signal — another symptom of the same root cause.

## Analysis: What Wave Equation Reproduces This?

**PRIMARY TARGET:** [Lafreniere Attraction](02_equations.md#primary-reference-two-opposite-phase-wave-centers)

Both animations share key features that constrain which wave equation form can reproduce them:

1. **Clear standing wave rings near each particle center** — concentric circles with fixed nodes. This requires a standing wave component near the core (rules out pure traveling wave forms)

2. **Traveling waves far from center** — the outer rings move outward. The wave transitions from standing to traveling character with distance

3. **Interference happens in the traveling wave region** — the force-producing pattern occurs where the traveling waves from both sources overlap between the particles. The standing wave region near each core remains relatively undisturbed

4. **The 1D envelope modulation** — the amplitude envelope (what the phasor RMS captures) shows a smooth modulation that increases (constructive) or decreases (destructive) between particles. This is the energy density landscape whose gradient produces force

## Critical Issues Summary — Electrostatic Force Validation

1. **Far-field oscillatory force (MAIN BLOCKER)** — Force direction flips every λ/2 of separation change. The sinc function sin(kr)/kr in the out-wave creates permanent nodes in the phasor RMS at all distances. These sinc nodes dominate over the charge-phase signal: source_offset (0 or π) should determine force direction via phasor P addition (same charge: cos(φ₁)=cos(φ₂) → P adds → larger amplitude) or cancellation (opposite charge: cos(φ₁)=-cos(φ₂) → P cancels → smaller amplitude). The spatial structure of where constructive/destructive interference occurs relative to each particle creates the gradient → force. But the sinc zeros override this charge-phase structure, making force direction depend on separation modulo λ instead of charge. Confirmed in both 3D and 1D engines
2. **Near-field opposite-phase monotonic attraction** — Opposite-charge WCs should always attract (to annihilation), not oscillate like same-charge lock-in. Currently both phase configurations show the same oscillatory behavior — same root cause as #1
3. **Gradient sampling radius** — Current gradient uses 1 grid unit. A real WC has spatial extent (~1λ). Wider gradient window may smooth out the sinc oscillation and reveal the underlying Coulomb trend. Untested
4. **Energy gradient ∇E vs amplitude gradient A·∇A** — Computing F = -∇E directly from E(x) = ρVf²A² may behave better numerically: A² is always positive, oscillations are compressed, and combined with wider gradient window may converge toward 1/r² more naturally. Mathematically identical but numerically different
5. **1/r² force law scaling** — Even if direction is fixed, F ∝ A·∇A with A ∝ 1/r gives F ∝ 1/r³, too steep. Need either a different envelope, a correction from the interference pattern between two sources, or a different force equation. Smoliński's "Degraded EMC Wall" concept may explain how the discrete sinc oscillation averages out to smooth 1/r² at macroscopic scales (see [03_additional.md](03_additional.md#the-degraded-emc-wall))
6. **Dual-treatment boundary** — Near-field needs raw oscillatory phasor (for lock-in physics), far-field needs smoothed envelope (for Coulomb). The weight function transition boundary should serve double duty, but this is unimplemented and untested

Issues 1, 2, 3, and 4 are likely connected — solving the gradient sampling or switching to ∇E may resolve the oscillatory behavior. Issue 5 may also resolve if the interference pattern between two sources modifies the effective amplitude scaling.

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

## POSSIBLE SOLUTION: Energy Gradient vs Amplitude Gradient

Current force: `F = -2ρVf² · A · ∇A` (from expanding `∇(A²)`)

Alternative: compute `F = -∇E` directly by first computing `E(x) = ρVf²A²`, then taking `∇E` numerically. While mathematically identical, the numerical behavior may differ:

- `A²` is always positive — no sign oscillation from sinc zeros
- `A²` oscillates at double frequency with different curvature
- Numerical gradient of the smoother `E(x)` function may average out sinc artifacts better than computing `A·∇A` from two oscillating quantities
- Combined with a wider gradient window, `∇E` on the squared envelope may converge toward 1/r² more naturally

This is a one-line change in the sandbox — worth testing early.

## POSSIBLE SOLUTION: Gradient Sampling Radius

The force gradient ∇A is currently computed per single grid spacing: `np.gradient(rms, dx)` uses central differences over 1 grid unit. This raises the question:

**Is 1 grid unit the physically correct scale for force computation?**

A real particle (wave center) doesn't experience force from a single point — it responds to the energy/amplitude curvature over a spatial region comparable to its own size (~1λ). The current single-grid-unit gradient captures the local slope of the sinc oscillation, which produces the oscillatory force artifacts.

Possible approaches to test:

- [ ] Compute gradient over a larger window (e.g., 1λ or 2λ) using weighted averaging — this would smooth out the sinc oscillation and extract the underlying 1/r envelope trend
- [ ] Use a Gaussian-weighted gradient kernel with σ ~ λ — physically motivated as the WC's "sensitivity radius"
- [ ] Compare force vs separation curves for different gradient sampling radii — which produces the smoothest 1/r² decay?
- [ ] This connects to the dual-treatment idea: near-field uses raw 1-grid gradient (captures lock-in wells), far-field uses smoothed gradient (captures Coulomb trend)

## POSSIBLE SOLUTION: Numerical Precision and Inertia Filtering

Two additional factors that may be masking or distorting the force signal:

- **Floating-point precision**: f32 arithmetic on very small force values (10⁻³ N scale at subatomic distances) may introduce numerical artifacts. The sinc function has near-zero values between nodes where precision loss is worst. Testing with f64 in the 1D sandbox (numpy default) should rule this out
- **Inertia as low-pass filter**: Real particles don't respond to 10²⁵ Hz oscillations — their mass acts as a low-pass filter, averaging out rapid force fluctuations. The force we compute should be the **time-averaged** force over at least one wave period, not the instantaneous gradient. The phasor RMS already time-averages amplitude — but if force depends on pressure (which is phase-shifted), the time-averaging may need to account for the phase relationship between displacement and pressure

## POSSIBLE SOLUTION: What Physical Quantity Drives Force?

**Deeper question: Is amplitude the only variable that creates force?** The current approach treats `F = -∇E` where `E = ρV(fA)²`, assuming ρ, V, and f are constants — so only the amplitude gradient drives force. But the energy equation has **three variables** that can change spatially and create energy gradients:

- **A** (amplitude): what we're currently computing via phasor RMS
- **f / λ** (frequency / wavelength): if λ varies with position (time dynamics, variable-λ research), then ∇f contributes to ∇E even where ∇A = 0
- **ρ** (density): if medium density varies with position (Smoliński's buoyancy model), then ∇ρ contributes to ∇E independently of amplitude

The full energy gradient expands to: `∇E = ρV · [2fA² · ∇f + f² · 2A · ∇A] + V(fA)² · ∇ρ`

This means force can arise from amplitude curvature (current model), frequency curvature (time dynamics), OR density curvature (buoyancy) — or any combination. The fact that we can't reproduce clean electrostatic force from amplitude gradient alone may be a signal that **the force equation is incomplete** — we may need to include ∇f and/or ∇ρ terms even for the electrostatic case. Wave interference between WCs doesn't just change amplitude — it also changes the effective local frequency and medium density through the mechanisms described in the [Time Dynamics](07_time_dynamics.md) and [Smoliński sections](03_additional.md#smolińskis-contributions-bcc-lattice-geometric-framework).

We're currently computing force from the amplitude gradient. But force may not respond to amplitude directly — it may respond to **pressure** or **energy density**, which are related but phase-shifted quantities:

- **Displacement** (what we compute as ψ): the wave displaces the medium from equilibrium
- **Velocity** (∂ψ/∂t): 90° phase-shifted from displacement — granule velocity
- **Pressure/density**: related to velocity (or compression), also 90° phase-shifted from displacement

In fluid mechanics, force derives from pressure gradients, not displacement gradients. If the wave equation is fundamentally propagating **pressure or density variations** (not just displacement), then force may be **90° phase-shifted** from where the displacement gradient predicts it. This would change which points in space experience force and could resolve the oscillatory artifact — the pressure nodes are at different spatial positions than the displacement nodes.

The chain: `displacement → velocity (90° phase) → pressure/compression → energy → force`. If force follows pressure gradient rather than displacement gradient, the sinc zeros that create our oscillatory force artifact would shift by λ/4 — potentially smoothing the pattern.
