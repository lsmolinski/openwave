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

---

## Base Wave Numerical Model — Analysis and Implementation Plan

### Prior Art: M1 and M2 Base Wave Implementations

**M1 (granule method)** modeled the base wave as background waves from 8 universe-corner vertices:

```text
ψ_base = Σ_{v=0}^{7}  (A₀/8) · cos(k·r_v - ωt) · (r̂_v)
```

Key properties: no 1/r falloff (full amplitude everywhere), amplitude split equally among 8 vertices for energy conservation, additive superposition onto WC source waves. Result: just more coherent waves on top of waves — still produces oscillatory interference when combined with sources.

**M2 (Laplacian grid)** used a fundamentally different approach:

1. Single pulse injected at center
2. Wave equation `∂²ψ/∂t² = c²∇²ψ` propagated via 6-point Laplacian + leap-frog integrator
3. Dirichlet BC (ψ=0 at walls) reflects waves back — emulating reflections from all matter in the universe
4. After many reflections, energy self-distributes into a quasi-uniform field

Key finding: **"get out of the way and let the waves stabilize themselves"** — the base wave emerged naturally from reflections without dynamic intervention. Stable over 20,000+ steps, ~100% energy conservation.

**What both approaches teach us**:

- M1's additive superposition (base + source) is exactly what equations #1–#6 already do — this cannot resolve the oscillatory force blocker
- M2's Laplacian propagation with reflecting boundaries is more physically honest — the base wave is not prescribed analytically, it emerges from wave dynamics
- M2 proved that a self-consistent standing wave field forms naturally from reflections — no external driving needed after initial energy injection
- Neither M1 nor M2 attempted WC-as-disturbance (both added WC waves on top of the base wave)

### What We Know vs What We Don't Know

**What we know**:

- The medium exists, is not empty — it has ρ₀, A₀, λ₀, f₀
- Waves come from all matter in the universe (all directions)
- Energy density is uniform (without WCs)
- λ₀ = 28.5 am, f₀ = 10²⁵ Hz, A₀ = 0.92 am

**What we don't know**:

- How granules oscillate in time at each point — are they all in phase? Random phase? Standing wave nodes?
- Is it a single coherent wave or a statistical/stochastic field?
- Does it have spatial structure (nodes) or is it spatially uniform?

### Candidate Base Wave Models for 1D Testing

The central question: what does the base wave look like in 1D? We test multiple candidates to discover which produces the right physics. Each is a selectable mode in the v3 wave engine.

#### Model A: Uniform Oscillation

```text
ψ_base(x,t) = A₀ · cos(ωt)
```

Every point oscillates together in phase. No spatial structure. Energy density is perfectly flat: `E = ρV(fA₀)²`. The simplest model — a uniform "hum."

- **Pro**: flat energy → zero gradient → no force from base wave alone (correct baseline)
- **Con**: not physically realistic — waves arriving from all directions cannot all be in phase everywhere. No wave character (no propagation, no interference)
- **Use case**: null baseline for comparison

#### Model B: Standing Wave (two counter-propagating waves)

In 1D, "waves from all directions" reduces to waves from left and right:

```text
ψ_base(x,t) = A₀/2 · cos(kx - ωt) + A₀/2 · cos(kx + ωt)
             = A₀ · cos(kx) · cos(ωt)
```

A standing wave with fixed nodes at every λ/2. Energy density has spatial structure: zero at nodes (`kx = nπ + π/2`), maximum at antinodes (`kx = nπ`). RMS envelope = `|A₀ · cos(kx)|/√2` — periodic, not flat.

- **Pro**: genuine wave physics, standing waves from counter-propagation (matches M2 behavior)
- **Con**: energy is NOT uniform — has zeros at nodes. WC placement relative to nodes matters. Node spacing = λ/2 matches the oscillatory force period we're trying to resolve (coincidence or clue?)
- **Use case**: test whether WC disturbance of a structured base wave produces different force behavior than WC-as-source. The node structure could interact with WC phase in interesting ways
- **Open question**: where are the nodes? The node positions depend on boundary conditions and are arbitrary in infinite space — but in a finite simulation domain, reflections fix them

#### Model C: Multi-Phase Superposition (stochastic isotropy)

Superpose N plane waves with random phase offsets to simulate isotropic arrival from many directions:

```text
ψ_base(x,t) = (A₀/√N) · Σᵢ [cos(kx + φᵢ - ωt) + cos(kx + φᵢ + ωt)]
            = (2A₀/√N) · Σᵢ cos(kx + φᵢ) · cos(ωt)
```

As N → ∞, the RMS envelope converges to a **spatially uniform** value (central limit theorem). The instantaneous displacement varies point-to-point but the time-averaged energy density is flat.

- **Pro**: physically motivated (waves from many sources = many random phases), uniform energy density in the statistical limit, genuine wave character at each point
- **Con**: requires many sources (N ≥ 50–100) for convergence, stochastic — different random seeds give different instantaneous patterns (but same RMS). Computationally heavier
- **Use case**: the "honest" isotropic model — if the base wave is truly from all matter in the universe, this is the closest 1D representation. Test whether WC disturbance of a stochastic field behaves differently from disturbance of a coherent field

#### Model D: Dual-Phase Standing Wave (complementary quadrature)

Two standing waves offset by a quarter wavelength (90° spatial phase):

```text
Channel 1:  ψ₁(x,t) = (A₀/√2) · cos(kx) · cos(ωt)
Channel 2:  ψ₂(x,t) = (A₀/√2) · sin(kx) · sin(ωt)
```

Each channel individually has nodes. But the nodes never coincide — where cos(kx) = 0, sin(kx) = ±1 and vice versa. The **energy sum** is:

```text
E_base(x) = ρVf² [A₁·cos(kx)]² + ρVf² [A₂·sin(kx)]²
           = ρVf² (A₀²/2) · [cos²(kx) + sin²(kx)]
           = ρVf² A₀²/2   ← FLAT!
```

Uniform energy density from two structured standing wave channels.

- **Pro**: elegant — each channel is a real standing wave with nodes and antinodes, but their energies complement to produce a perfectly flat energy field. Wave character preserved while energy is uniform
- **Con**: requires two independent "channels" that don't interfere with each other — is this physically justified? Why would the medium support two orthogonal wave modes?

**Additional test: π-apart dual phase** — two standing waves offset by half a wavelength (180° spatial phase):

```text
Channel 1:  ψ₁(x,t) = (A₀/√2) · cos(kx) · cos(ωt)
Channel 2:  ψ₂(x,t) = (A₀/√2) · cos(kx + π) · cos(ωt)
           = -(A₀/√2) · cos(kx) · cos(ωt)
```

Direct superposition: `ψ₁ + ψ₂ = 0` (total cancellation). But energy sum: `E₁ + E₂ = ρVf²A₀²cos²(kx)` — still has nodes, same as Model B. The π-apart case only works as separate energy channels (not field superposition), and even then doesn't produce flat energy. This contrasts with the 90°-apart case where energy complementarity is exact.

**Both angular offsets should be tested** to understand how the phase relationship between base wave channels affects energy uniformity and WC interaction.

- **Use case**: test whether WCs preferentially lock onto one channel (creating charge states?). The dual-phase speculation in the main doc — "WCs lock onto one mode or the other (source_offset = 0 or π), creating the two charge states" — directly maps to this model
- **Connection to charge**: if the base wave has two complementary modes, and WCs disturb one mode more than the other based on their phase, this could be the mechanism for charge-dependent force direction

#### Model E: Laplacian Propagation (M2 port to 1D)

Instead of an analytical formula, use the actual wave equation with reflecting boundaries:

```text
∂²ψ/∂t² = c² · ∂²ψ/∂x²
```

Initialize with a pulse (or boundary oscillators), let it ring and stabilize. The base wave emerges from simulation dynamics — no assumed spatial form.

- **Pro**: most physically honest — no assumptions about the base wave form. The wave equation itself determines what the field looks like. Matches M2's approach (which successfully self-stabilized). Naturally handles WC interaction via the same wave equation (WC = boundary condition or source term within the Laplacian domain)
- **Con**: architecturally different from analytical phasor computation — requires time-stepping, warmup period for stabilization, and a different code path. Cannot use phasor superposition (which assumes known analytical form). More complex to implement and analyze
- **Use case**: the "ground truth" model. If any analytical candidate (A–D) matches the Laplacian result, we know that candidate is correct. If none match, the Laplacian reveals what we're missing

### On the "Longitudinal" Assumption

The base wave is described as longitudinal in EWT literature, but **this is a hypothesis, not established fact**. In 1D simulation, we are forced to isolate a single wave mode (longitudinal displacement along x). This is a limitation of the dimensionality, not necessarily of the physics:

- **1D**: only longitudinal mode available — displacement along the propagation axis
- **2D**: two modes possible — longitudinal (along propagation) + transverse (perpendicular)
- **3D (reality)**: full elliptical motion — the granule displacement traces an ellipse in 3D space, characterized by: longitudinal amplitude, transverse amplitude, handedness (direction of granule motion around the elliptical track), and ellipse plane orientation

All of these properties can contribute to force direction. The unified force concept proposes that what we perceive as separate forces (electric, magnetic, gravitational) are actually **one 3D elliptical behavior** decomposed into components. At human scale, specific conditions make each component appear distinct — we named and described them as separate forces because that's how they manifest at our scale of mass, frequency, and inertial frame.

A 2D simulation could capture longitudinal + transverse, but **may not be sufficient** — the elliptical form can be oriented in multiple ways in 3D space, and this orientational freedom is likely essential for magnetic fields and spin. The 1D base wave work here is foundational (establishing the energy redistribution mechanism), but the full picture likely requires 3D (Phase 1d / M4).

### Implementation Plan for v3 Wave Engine

The v3 engine (`wave_engine_1D_v3.py`) replaces the v2 equation selector (#1–#6) with:

```text
BASE_WAVE_MODE:
  A = "uniform"        — ψ = A₀·cos(ωt), flat energy
  B = "standing"       — ψ = A₀·cos(kx)·cos(ωt), nodes at λ/2
  C = "stochastic"     — N random-phase standing waves, ~flat energy
  D = "dual_phase"     — two 90°-offset standing waves, flat energy (test π-offset too)
  E = "laplacian"      — time-stepped wave equation, reflecting BC

WC_INTERACTION_MODE:   (Phase 2 — after base wave is validated)
  [to be determined]   — how WCs disturb the base wave
```

**Step 1 — Base wave only (no WCs)**:

- Implement all 5 modes
- Visualize: displacement ψ(x,t), RMS envelope, energy density E(x)
- Verify: uniform energy for modes A, C, D; structured energy for B; emergent structure for E
- Compare: do any analytical modes (A–D) match the Laplacian result (E)?

**Step 2 — Insert WC disturbance**:

- Design WC interaction mechanism (how does a WC modify the base wave field?)
- Test energy redistribution: concentration near WC, drainage in far field
- Test phase dependence: does WC phase (0 vs π) affect the drainage pattern?
- Measure force: F = -∇E at WC positions

Each base wave mode may interact with WCs differently — that's the point of testing all of them.

### Recommended Implementation Order

1. **Model D (dual-phase)** first — the cleanest analytical model that gives uniform energy density while maintaining real wave structure. Good baseline for verifying the v3 engine works correctly
2. **Model E (Laplacian)** second — the "ground truth" time-stepped model for comparison against analytical candidates. Architecturally different (time-stepping vs phasor), so validates from a different angle
3. **Model B (standing wave)** third — shows what happens with nodes. WC disturbance interaction with node structure could be revealing, and the λ/2 node spacing matching the oscillatory force period may be significant
4. **Models A and C** as needed — A is the trivial null baseline, C is the computationally heavier stochastic model (useful if coherent models fail)
