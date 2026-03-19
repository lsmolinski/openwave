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

## Force Regime Matrix

| Regime     | Same Phase                        | Opposite Phase                                |
| ---------- | --------------------------------- | --------------------------------------------- |
| Near-field | Lock-in (quarks, orbits, bonding) | Attraction → annihilation (wave cancellation) |
| Far-field  | Constructive → repulsion          | Destructive → attraction                      |

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

## ✅ RULED OUT: Numerical Precision

- **Floating-point precision**: The 1D sandbox uses numpy f64 (15 decimal digits) — more than sufficient. All physical constants are already in simulation-friendly units (am, rs, qg) with magnitudes in the O(0.01)–O(10) range, avoiding large/small number cancellation. The sinc zeros are exact zeros of sin(kr), not near-cancellation artifacts. M3/M4 use f32 but with the same unit scaling. **Conclusion**: the oscillatory force is a real mathematical feature of the sinc spatial function, not a floating-point artifact

## ✅ RULED OUT: Inertia Filtering

- **Inertia as low-pass filter**: Real particles don't respond to 10²⁵ Hz oscillations — their mass acts as a low-pass filter, averaging out rapid force fluctuations. **But the phasor RMS already IS this filter** — it gives the exact analytical time-averaged amplitude without simulating any oscillation cycles. The force we compute from ∇E(RMS) is already the time-averaged force. The oscillatory problem is in the *spatial* sinc structure of the RMS envelope, not in temporal high-frequency artifacts that further averaging could remove

## ✅ RULED OUT: Pressure/Velocity Gradient (90° Phase Shift)

The chain: `displacement → velocity (90° phase) → pressure/compression → energy → force`. If force follows pressure gradient rather than displacement gradient, the sinc zeros that create our oscillatory force artifact would shift by λ/4 — potentially breaking the exact cancellation pattern that makes charge-phase information unrecoverable.

- **Displacement** (what we compute as ψ): the wave displaces the medium from equilibrium
- **Velocity** (∂ψ/∂t): 90° phase-shifted from displacement — granule velocity
- **Pressure/density**: related to velocity (or compression), also 90° phase-shifted from displacement

A 90° phase shift does not resolve the oscillatory force:

- **Velocity (∂ψ/∂t)**: velocity RMS = ω × displacement RMS. Since ω is constant, ∇(ω·RMS) = ω·∇(RMS) — identical gradient direction, zero benefit
- **Pressure (∝ -∂ψ/∂x)**: gives d/dr[sin(kr)/kr] instead of sin(kr)/kr — different zeros, but still oscillates with the same λ/2 period, just shifted by λ/4

No phase rotation or derivative of a periodic function removes its periodicity. The sinc zeros are λ/2-spaced, and any linear operation (derivative, phase shift, scaling) preserves that spacing. The problem would persist, just at different separation values.

**Note**: granule velocity is still physically significant — it relates directly to medium density ρ (see [Time Dynamics](07_time_dynamics.md): faster cycling = higher pressure/density). This connection is explored below in the Non-Linear Wave Equations section, where ∇ρ could contribute force independently of ∇A.

## ✅ RULED OUT: Standing vs Traveling Wave Decomposition

Decompose the phasor into standing-wave-only and traveling-wave-only components, test force from each separately.

**Single WC analysis**: each component individually is smooth:

- **Standing** (in-wave): RMS = A·w(r)/kr — smooth, dies off with weight function
- **Traveling** (out-wave): RMS = A/kr — smooth 1/r everywhere
- The sinc oscillation `√((w+1)²sin²(kr) + (w-1)²cos²(kr))/kr` only appears when combining both components of the same WC

**Two WC analysis**: even with traveling-wave-only (smooth 1/kr per source), coherent superposition creates oscillatory interference. The interaction energy between two traveling waves contains:

```text
E_interaction ∝ cos(k(r₁ - r₂) - (φ₁ - φ₂)) / (kr₁ · kr₂)
```

The cos() term oscillates with separation at spatial frequency k — this is the same λ/2 force flip. **The oscillation is intrinsic to coherent wave interference**, not to standing vs traveling character. This explains why the smooth envelope test (|Z₁|·|Z₂|) works — it avoids coherent superposition entirely by computing each source's smooth amplitude separately.

## [ ] Test Alternative Wave Equations (Quick Sweep)

The 1D engine implements 5 wave equation forms ([02_equations.md](02_equations.md)), each with different spatial functions. We've only tested the force sweep with #5 (weighted partial standing wave). The other forms have **different zero structures** that may change the oscillation behavior:

| Equation | Spatial zeros | Zero spacing | Notes |
| --- | --- | --- | --- |
| #1 Wolff | sin(kr)/r | λ/2 | Pure standing wave everywhere |
| #2 LaFreniere-Marcotte | sin(kr)/(kr) + (1-cos(kr))/(kr) | λ | Wider spacing, partially traveling |
| #3 Phase-warped Marcotte | sin(x_c)/x_c | ~λ (warped near core) | Core correction shifts zeros |
| #4 Combined Wolff-LF | sin(kr)/r + (1-cos(kr))/r | λ/2 | Wolff normalization |
| #5 Weighted (current) | (w+1)sin(kr)/(kr) + (w-1)cos(kr)/(kr) | depends on w(r) | Far-field smooth, near-field λ/2 |

**Hypothesis**: equations #2 and #3 (LaFreniere-Marcotte) have zeros at **λ spacing** instead of λ/2. The force oscillation period may differ, the oscillation depth may be shallower, and the interaction with charge-phase (cos(Δφ)) may produce better direction accuracy at certain separations. The phasor decomposition for all 5 forms is already implemented in `wave_engine_1D_v2.py`.

**Test**: run `sweep_force_vs_separation.py` with each equation form, compare direction accuracy and oscillation pattern. Quick to implement — just switch the equation selector and run.

**Expected outcome**: the oscillation is likely still present (coherent interference is fundamental), but comparing across equations may reveal which spatial function best preserves charge-direction information, and whether LaFreniere's reference animations correspond to a specific equation form.

## POSSIBLE SOLUTION: Non-Linear Wave Equations ([Phase 1b fallback](01_plan.md#phase-1b-non-linear-wave-equations-1d--details) / [Phase 4](01_plan.md#phase-4-non-linear-wave-equations-m3-1d--3d--details))

**Why non-linear?** All linear operations on the sinc function (smoothing, phase shift, derivatives, scaling) preserve its λ/2-periodic zero structure. Only a non-linear wave equation — where the spatial function is no longer a pure sinc — can break the oscillatory pattern that causes force direction errors.

**Smoliński's Non-linear Soliton Wave Equation** (MagnetismGravity_v2, Sec 6.1, Eq. 18-19): The soliton is governed by a non-linear wave equation with a cubic stabilizing term:

```text
(∂²/∂t² - c²∇²) Ψ(r,t) + F(Ψ, ε_G, |ε_M|, N_ν) = 0

where F = k(|ε_M|) · Ψ³    (NLS cubic non-linearity)
```

The Ψ³ term prevents wave dispersion and creates soliton stability. The coefficient k(|ε_M|) describes the non-linear elasticity of the medium, depending explicitly on the magnetic deficit. **Key insight**: solutions to the non-linear wave equation with cubic term are NOT pure sin(kr)/kr — the soliton spatial structure differs from a sinc, potentially breaking the periodic zero pattern that causes our oscillatory force. This is a well-studied form (Non-linear Schrödinger solitons) with known analytical solutions.

**Is amplitude the only variable that creates force?** The current approach treats `F = -∇E` where `E = ρV(fA)²`, assuming ρ, V, and f are constants — so only the amplitude gradient drives force. But the energy equation has **three variables** that can change spatially and create energy gradients:

- **A** (amplitude): current phasor RMS — carries sinc oscillation from sin(kr)/kr spatial structure
- **f / λ** (frequency / wavelength): if λ varies with position ([Yee & Hauger](07_time_dynamics.md#ewt-standing-wave-geometry-yee--hauger) discrete wavelength shells), then ∇f contributes to ∇E even where ∇A = 0. Crucially, variable λ(r) replaces the uniform kr phase with a WKB integral ∫k(r')dr', making node spacing non-uniform — this **breaks the sinc periodicity** that causes the oscillatory force. **Smoliński r⁵ decomposition** (MagnetismGravity_v4 Sec 7.2): E ∝ r⁵ = r³ (volume) × r¹ (amplitude: A ∝ r) × r¹ (frequency: f ∝ 1/r). Inside the soliton, amplitude scales linearly with distance from center and frequency scales inversely — both are concrete non-linear relationships that modify the wave equation near WCs
- **ρ** (density): if medium density varies with position ([Smoliński's buoyancy model](03_additional.md#smolińskis-contributions-bcc-lattice-geometric-framework)), then ∇ρ contributes to ∇E independently of amplitude. **Smoliński's explicit density function** (MagnetismGravity_v4 Eq. 32): `ρ(r) = ρ₀(1 - (r/r_ν)^k)^P · Θ(r_ν - r)` — packing density decreases from soliton boundary toward core, only within soliton radius r_ν (Heaviside cutoff). This is a concrete implementable ansatz. **Connection to granule velocity**: ρ is directly related to granule velocity — the speed at which granules cycle through their elliptical motion (∂ψ/∂t). Faster cycling = higher local density/pressure, slower cycling = lower density. Wave interference between WCs doesn't just change displacement amplitude — it changes granule velocities, which changes local ρ. This means ∇ρ may carry force information that ∇A alone cannot provide, and with a different spatial structure than the sinc envelope

This means force can arise from amplitude curvature (current model), frequency curvature (time dynamics), OR density curvature (buoyancy) — or any combination. The fact that we can't reproduce clean electrostatic force from amplitude gradient alone may be a signal that **the force equation is incomplete** — we may need to include ∇f and/or ∇ρ terms even for the electrostatic case. Wave interference between WCs doesn't just change amplitude — it also changes the effective local frequency and medium density through the mechanisms described in the [Time Dynamics](07_time_dynamics.md) and [Smoliński sections](03_additional.md#smolińskis-contributions-bcc-lattice-geometric-framework).

Computing F = -∇E means these additional variables are automatically included when implemented — no force logic changes needed.

**Smoliński's Push-out Operator** (MagnetismGravity_v4 Eq. 90) formalizes exactly this multi-variable approach:

```text
P̂Φ = -∇ · (η_stat / η_soliton) ∇Φ
```

This is F = -∇E but with the impedance ratio (statutory density / soliton density) as a spatially varying coefficient. The force is the gradient of a potential weighted by the local density mismatch — precisely what our ∇E captures when ρ(x) becomes a spatial variable.

**Connection to dual-treatment boundary**: Smoliński's Isotropy Operator I (MagnetismGravity_v4 Sec 17.7.1) acts as a **geometric low-pass filter** at the Degraded EMC Wall boundary, averaging angular energy into isotropic gravity. This maps to our weight function transition: inside the boundary → non-linear toroidal dynamics (Energy Domain, r⁵), outside → isotropic spherical push-out (EMC Domain, r³). The weight function could serve as this geometric filter.

**Scheduling**: This is a fallback if all Phase 1 linear approaches fail (expensive implementation: variable λ and ρ fields in 1D sandbox). If Phase 1 succeeds without it, remains as Phase 4 for 3D porting. See [01_plan.md](01_plan.md#phase-1b-non-linear-wave-equations-1d--details) for tasks.

## POSSIBLE SOLUTION: Energy Flux (Radiation Pressure)

Force can also arise from **energy flux** — the directional flow of energy through the medium:

- **Energy density** (current: `F = -∇E`): energy *stored* per voxel. Static snapshot. Force from the energy landscape shape (conservative)
- **Energy flux** (`S = E · v_group`): energy *flowing* through a surface per unit time (W/m²). Has direction — tells which way energy is moving
- **Radiation pressure** (`P_rad = S/c`): force per unit area from wave momentum transfer. This is LaFreniere's force mechanism — waves exert pressure on wave centers

For the current monochromatic case with constant c: `P_rad = S/c = E·c/c = E`, so `∇P_rad = ∇E` — both approaches give the same force. They diverge when:

- **c varies spatially** (variable medium density → Smoliński's buoyancy model)
- **Waves are directional** (flux has direction, density doesn't)
- **Standing vs traveling waves**: standing waves store energy but don't flow (zero flux, nonzero density). Traveling waves carry energy directionally (nonzero flux)

This last point is significant: energy flux could **naturally separate standing wave (near-field) from traveling wave (far-field) contributions** — standing waves have zero net flux, traveling waves have nonzero flux. This connects to the dual-treatment boundary idea: near-field force from energy density gradient (standing wave lock-in), far-field force from energy flux gradient (traveling wave Coulomb regime).
