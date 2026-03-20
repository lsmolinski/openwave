# PHASE 1b: Base Wave + WC Energy Redistribution

Model the actual base wave and how WCs redistribute its energy. Central open question: how WC phase determines the drainage pattern

## ❌ The Problem

The main blocker is still the far-field oscillatory force: the sinc function sin(kr)/kr in the out-wave creates permanent nodes in the phasor RMS at all distances. Force direction flips every λ/2 of separation change, even where only smooth 1/r decay should exist. Confirmed in both 3D and 1D engines.

**Force computation**: F = -∇E where E(x) = ρ·V·(f·A(x))². Computing from ∇E directly (not chain-rule expansion) ensures future variable ρ(x), f(x), λ(x) are automatically captured. The same pattern applies across 1D sandbox, M3, and M4 engines.

**All linear candidates exhausted (10/10 ruled out).** The oscillatory force is a fundamental consequence of coherent wave interference — no linear wave equation or operation on the superposed field can eliminate it while preserving charge-dependent direction. Signed non-wave approaches (Phase 1a) bypass interference entirely but are "Coulomb with extra steps" — charge sign is imposed, not emergent.

**Validation targets**: plot E(x) along the axis connecting two particles at various separations, identify constructive/destructive interference locations, verify gradient direction and 1/r² magnitude scaling against Coulomb reference.

All task checklists are tracked in [00_roadmap.md](00_roadmap.md).

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

**With zero WCs**: the sandbox should display this uniform energy field — constant energy density, constant amplitude, but displacement oscillating with wave character (not flat).

**Why fundamentally different**: all 9 previous candidates modeled WCs as wave SOURCES emitting into empty space. The base wave model changes the paradigm — WCs are disturbers of an existing field, not sources.

**Prior art**: M1 (granule method) implemented base wave as background waves from 8 universe vertices (A·cos(kr - ωt)·direction / 8, no 1/r falloff, `BASE_WAVE_TOGGLE`). M2 (grid/Laplace) used boundary wall oscillators. Both used additive superposition (base_wave + source_waves), which still produces oscillatory interference. The true base wave + disturbance model must go beyond additive superposition.

**The base wave concept is not just context — it is the energy source.** The medium has stored energy (standing waves) as potential for matter, forces, EM waves, heat, and time. WCs redistribute this energy, and the redistribution pattern (affected by WC phase) creates the gradients that produce force.

The base wave (1b) provides the energy field that WCs disturb; non-linear equations (1c) describe how the disturbance propagates; vector displacement (1d) may carry the charge information that scalar magnitude discards. They may ultimately converge into a single solution.

**The open question**: can the reflected waves produce a far-field energy gradient that is smooth enough to avoid the sinc oscillation problem, while still being genuine wave interference?

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

**WC disturbance scope**: the disturbance affects not just amplitude but also **wavelength λ** and potentially **density ρ** near the WC. This connects to the multi-variable energy gradient (∇A + ∇f + ∇ρ) and to non-linear wave equations (Phase 1c).

**Scalar base → vector emergence**: the fundamental base wave might be scalar (longitudinal only). Vector (transverse) waves emerge from **spin** — the WC's toroidal wave rotation converts longitudinal to transverse.

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
