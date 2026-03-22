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

---

### Candidate Base Wave Models for 1D Testing

The central question: what does the base wave look like in 1D? We test multiple candidates to discover which produces the right physics. Each is a selectable mode in the v3 wave engine.

## ✅ Uniform (`BASE_WAVE_MODE = "uniform"`)

```text
ψ_base(x,t) = A₀ · cos(ωt)
```

Every point oscillates together in phase. No spatial structure. Energy density is perfectly flat: `E = ρV(fA₀)²`. The simplest model — a uniform "hum."

- **Pro**: flat energy → zero gradient → no force from base wave alone (correct baseline). Simple and beautiful mathematical manifold
- **Con**: not physically realistic — waves arriving from all directions cannot all be in phase everywhere. No wave character (no propagation, no interference)
- **Use case**: null baseline for comparison

**Test result**: RMS is stable and uniform, energy field perfectly flat. Works exactly as expected.

**Idea for WC interaction**: add a second phase wave (π apart) — the two waves sum to zero energy, but WCs disturb one or the other depending on charge sign. This creates a zero-energy vacuum that WCs can asymmetrically perturb.

## ✅ Standing Wave (`BASE_WAVE_MODE = "standing"`)

In 1D, "waves from all directions" reduces to waves from left and right:

```text
ψ_base(x,t) = A₀/2 · cos(kx - ωt) + A₀/2 · cos(kx + ωt)
             = A₀ · cos(kx) · cos(ωt)
```

A standing wave with fixed nodes at every λ/2. Energy density has spatial structure: zero at nodes (`kx = nπ + π/2`), maximum at antinodes (`kx = nπ`). RMS envelope = `|A₀ · cos(kx)|/√2` — periodic, not flat.

- **Pro**: genuine wave physics, standing waves from counter-propagation. **Physically validated**: the Laplacian propagation model self-stabilizes into exactly this pattern — confirming that counter-propagating waves in a reflecting domain naturally produce a standing wave field
- **Con**: energy is NOT uniform — has zeros at nodes. WC placement relative to nodes matters. Node spacing = λ/2 matches the oscillatory force period we're trying to resolve (coincidence or clue?)
- **Use case**: the physically motivated analytical model. Computationally optimized replacement for the Laplacian (same result, no warmup needed)
- **Open question**: where are the nodes? The node positions depend on boundary conditions and are arbitrary in infinite space — but in a finite simulation domain, reflections fix them

**Test result**: clean standing wave pattern. Nodes and antinodes clearly visible in energy panel.

**Scale argument for node suppression**: the energy field is discrete at base wave scale (λ = 28 am), but the smallest stable particle (electron) has radius K²λ = 100λ = 2800 am. At that scale, the standing wave nodes are 100× smaller than the particle — the electron "sees" an averaged-out, effectively continuous energy field. Like pixels: when resolution is much finer than the object, discreteness disappears. This is Planck-scale granule motion — at electron scale, the base wave's node structure averages out.

---

## ✅ Hypothesis: Node-Locking Charge (standing wave model) — FALSIFIED (Step 2a)

Hypothesis: charge is not intrinsic to the WC — it emerges from the WC's position in the standing wave node lattice. Even-node WCs have one "charge," odd-node WCs have the opposite. The sinc λ/2 force oscillation IS the charge mechanism itself.

**Prediction**: force alternates every node step (λ/2) — odd separations attract, even separations repel.

### Test Implementation (`wave_engine_1D_v3_vnode.py`, `sweep_node_locking.py`)

- Base wave: standing mode `A₀·cos(kx)·cos(ωt)`
- Two WCs with phase = 0 (no imposed charge), snapped to node positions
- Combined phasor: base standing wave + WC weighted partial standing waves (eq #5)
- Swept separation from 1 to 30 nodes (0.5λ to 15λ)

### Results — FALSIFIED

```text
Overall: 7/30 match (23%) — not consistent
ODD  separations: 7/15 attract (as predicted), 8/15 repel (opposite)
EVEN separations: 0/15 repel, 15/15 unclear (both forces same direction)
Force does NOT alternate every λ/2 — 0/29 consecutive transitions alternate
```

**Actual pattern observed**: the force has period **2λ (4 nodes)**, not λ/2 (1 node):

- Odd separations of form 4n+3 (sep 3, 7, 11, 15, 19, 23, 27) → **always attract** (7/7)
- Odd separations of form 4n+1 (sep 1, 5, 9, 13, 17, 21, 25, 29) → **always repel** (8/8)
- Even separations (all) → **unclear** — both WCs pushed in same direction (net translation, not attraction/repulsion)

### Why It Failed

1. **The even/odd node parity is too simple** — the actual force pattern involves the interaction between the WC's own sinc wave structure (which has its own λ/2 periodicity) and the base wave's node lattice. The WC wave dominates over the base wave node structure
2. **Even separations don't produce repulsion** — they produce asymmetric forces where both WCs are pushed the same direction, suggesting the base wave gradient (not WC interaction) dominates at same-type nodes
3. **The 2λ periodicity** comes from the WC wave's sinc function interfering with the base wave — a beat pattern between two periodic structures, not a simple node-type alternation

### What We Learned

- Charge cannot be purely spatial (position-dependent) in this model — the WC's own wave structure imposes its own periodicity that overrides the base wave's node pattern
- The standing wave base field is still a valid contender for WC disturbance testing, but the force mechanism must come from WC wave interference, not from node position assignment
- The "unclear" even-separation result (net translation instead of attraction/repulsion) is interesting — may indicate that same-type node positions create a field asymmetry rather than a symmetric force

## ✅ Stochastic (`BASE_WAVE_MODE = "stochastic"`)

**Original model (monochromatic)** — N plane waves with random phase offsets at the SAME wavenumber k:

```text
ψ_base(x,t) = (2A₀/√N) · Σᵢ cos(kx + φᵢ) · cos(ωt)
```

**Bug found**: `Σ cos(kx + φᵢ)` with same k and random φᵢ collapses mathematically to `|S|·cos(kx + arg(S))` where S = Σ e^{iφᵢ} — a single complex number. The result is always a **single phase-shifted standing wave**, identical to the standing wave model. Random phases at the same wavenumber cannot produce spatial uniformity in 1D. This is a fundamental mathematical limitation, not a code bug.

**Fixed model (broadband)** — each source has its own k_i spread around the base k:

```text
k_i = k · (1 + σ · δᵢ),   δᵢ ∈ [-1, 1] random,   σ = STOCHASTIC_K_SPREAD
ω_i = c · k_i

ψ_base(x,t) = amp · Σᵢ cos(k_i·x + φᵢ) · cos(ω_i·t)
```

With different k_i, spatial patterns don't combine into a single sinusoid. Different ω_i means cross terms average to zero over time → quasi-uniform RMS. Parameter `STOCHASTIC_K_SPREAD` controls bandwidth: 0 = monochromatic (collapses to standing wave), 1 = full octave spread (maximum spatial uniformity).

- **Pro**: demonstrates broadband vs monochromatic wave field behavior. Energy field averages out with bandwidth, and roughness further suppresses at larger scales (electron sees smoothed-out fluctuations, like smaller pixels)
- **Con**: hard to find physical motivation — why would the base wave be broadband? EWT describes a single fundamental wavelength λ₀. The stochastic model may represent statistical noise in the medium, but not the fundamental base wave itself
- **Use case**: useful as mathematical demonstration of how bandwidth creates spatial uniformity. May become relevant if the medium has λ variation from non-linear effects (Phase 1c)

**Test result** (σ=1): chaotic displacement, but energy field averages to near-uniform. Higher N and broader σ → flatter energy. The pixel analogy holds: at particle scale (100λ), small variations disappear.

## ✅ Quadrature (`BASE_WAVE_MODE = "quadrature"`)

Two standing waves with configurable spatial and temporal offsets:

```text
Channel 1:  ψ₁(x,t) = A₀ · cos(kx) · cos(ωt)
Channel 2:  ψ₂(x,t) = A₀ · cos(kx + δs) · cos(ωt + δt)

DUAL_SPATIAL_OFFSET  = δs (default π/2)
DUAL_TEMPORAL_OFFSET = δt (default π/2)
```

Per-channel energy sum (independent channels, not field superposition):

```text
E_base(x) = ρVf² A₀² [cos²(kx) + cos²(kx + δs)] / 2
```

For δs = π/2: cos²(kx) + sin²(kx) = 1 → **FLAT energy** ✓
For δs = π: 2cos²(kx) → **nodes** (same as standing wave)

The **spatial offset** controls energy flatness. The **temporal offset** controls the combined displacement waveform without affecting energy:

- δt = +π/2 → traveling wave moving right: `A₀·cos(kx - ωt)`
- δt = -π/2 → traveling wave moving left: `A₀·cos(kx + ωt)`
- δt = 0 → standing wave (shifted): `A₀√2·cos(kx - π/4)·cos(ωt)`

- **Pro**: elegant — each channel is a real standing wave with nodes and antinodes, but their energies complement to produce a perfectly flat energy field. Wave character preserved while energy is uniform. This may be the physical meaning behind complex sinusoids in wave equations (Schrödinger, Wolff) — the real and imaginary components represent two 90°-offset standing wave channels. The quadrature structure may encode the longitudinal + transverse wave components, and possibly spin, in a 1D representation of 3D elliptical granule motion
- **Con**: requires two independent "channels" that don't interfere with each other — is this physically justified? Why would the medium support two orthogonal wave modes? In 1D, only two travel directions exist (left/right), but in 3D space, orientational freedom is much richer — the two-direction model may be an oversimplification

**Test result**: perfectly flat uniform energy (at δs = π/2). Displacement shows traveling wave that **inverts direction** when temporal offset sign flips.

**Traveling wave direction and charge/spin**: the direction flip (left ↔ right) with temporal offset sign is physically interesting — it could encode charge sign or spin handedness. WCs that "tap" into one direction vs the other would see different energy landscapes. The challenge: in 3D, there are infinitely many directions, not just two — so mapping this to real charge physics requires understanding how the 3D elliptical motion reduces to the 1D two-channel picture.

**π-offset variant** (δs = π): direct superposition cancels to zero, per-channel energy sum has nodes (same as standing wave). Only the 90°-apart case gives flat energy — confirming that quadrature (not arbitrary offset) is the special configuration. Testable in v3 by setting `DUAL_SPATIAL_OFFSET = π`.

**Connection to charge**: WCs may preferentially lock onto one channel (creating charge states). If the base wave has two complementary modes and WCs disturb one mode more than the other based on their phase, this could be the mechanism for charge-dependent force direction.

## ✅ Laplacian (`BASE_WAVE_MODE = "laplacian"`) — RETIRED

Time-stepped wave equation with reflecting boundaries (Dirichlet BC):

```text
∂²ψ/∂t² = c² · ∂²ψ/∂x²
ψ(x_min) = ψ(x_max) = 0
```

Initialized with Gaussian pulse or analytical standing wave. After warmup (20+ periods), the field self-stabilizes. Live RMS tracking updates every period via peak |ψ| measurement.

- **Pro**: most physically honest — no assumptions about the base wave form. The wave equation itself determines what the field looks like
- **Con**: architecturally different from analytical computation — requires time-stepping, warmup period, and separate code path

**Test result**: the Laplacian **resolves itself into a standing wave** — identical to the standing wave model. This provides **physical motivation for the analytical standing wave model**: it represents real Laplacian wave propagation in a computationally optimized analytical form. The standing wave isn't an arbitrary choice — it's what the wave equation naturally produces in a reflecting domain.

**Status: RETIRED** — the Laplacian mode validated that the standing wave is the physically correct 1D base wave form. It can now be replaced by the analytical standing wave model for all further work. Kept in v3 code for reference but not needed for Phase 1b Step 2.

## On the "Longitudinal" Assumption

The base wave is described as longitudinal in EWT literature, but **this is a hypothesis, not established fact**. In 1D simulation, we are forced to isolate a single wave mode (longitudinal displacement along x). This is a limitation of the dimensionality, not necessarily of the physics:

- **1D**: only longitudinal mode available — displacement along the propagation axis
- **2D**: two modes possible — longitudinal (along propagation) + transverse (perpendicular)
- **3D (reality)**: full elliptical motion — the granule displacement traces an ellipse in 3D space, characterized by: longitudinal amplitude, transverse amplitude, handedness (direction of granule motion around the elliptical track), and ellipse plane orientation

All of these properties can contribute to force direction. The unified force concept proposes that what we perceive as separate forces (electric, magnetic, gravitational) are actually **one 3D elliptical behavior** decomposed into components. At human scale, specific conditions make each component appear distinct — we named and described them as separate forces because that's how they manifest at our scale of mass, frequency, and inertial frame.

A 2D simulation could capture longitudinal + transverse, but **may not be sufficient** — the elliptical form can be oriented in multiple ways in 3D space, and this orientational freedom is likely essential for magnetic fields and spin. The 1D base wave work here is foundational (establishing the energy redistribution mechanism), but the full picture likely requires 3D (Phase 1d / M4).

---

## Implementation Plan for v3 Wave Engine

The v3 engine (`wave_engine_1D_v3.py`) replaces the v2 equation selector (#1–#6) with:

```text
BASE_WAVE_MODE:
  "uniform"    — ψ = A₀·cos(ωt), flat energy
  "standing"   — ψ = A₀·cos(kx)·cos(ωt), nodes at λ/2
  "stochastic" — N broadband random-phase waves, ~flat energy
  "quadrature" — two 90°-offset standing waves, flat energy
  "laplacian"  — time-stepped wave equation, reflecting BC (RETIRED)

WC_INTERACTION_MODE:   (Step 2 — after base wave is validated)
  [to be determined]   — how WCs disturb the base wave
```

## ✅ Step 1 — Base wave only (no WCs) — COMPLETE

All 5 modes implemented and tested:

- ✅ **uniform**: flat energy, flat RMS — works as null baseline
- ✅ **standing**: nodes at λ/2, structured energy — physically validated by Laplacian
- ✅ **stochastic**: monochromatic bug found and fixed (broadband), quasi-uniform energy at σ=1
- ✅ **quadrature**: perfectly flat energy (at 90° offset), traveling wave direction flips with sign
- ✅ **laplacian**: self-stabilizes to standing wave — validates standing wave model, RETIRED

**Key finding**: Laplacian resolves to standing wave → standing wave is the physically correct 1D base wave form. Quadrature produces flat energy via two complementary channels — most mathematically elegant and may encode real physics (complex sinusoids, spin, charge).

## Step 2 — WC Disturbance and Contender Selection

### ✅ Step 2a: Node-Locking Charge Hypothesis — FALSIFIED

Tested whether charge emerges from WC position in standing wave node lattice. Prediction: force alternates every λ/2, odd separations attract, even repel. Result: 7/30 match (23%). Actual pattern has 2λ periodicity from WC sinc interacting with base wave nodes. Even separations produce net translation, not repulsion. See [Node-Locking Hypothesis](#-hypothesis-node-locking-charge-standing-wave-model--falsified-step-2a) for full results.

---

### ✅ Step 2b: Migrate WC Disturbance to v3 — COMPLETED, additive model ruled out

Ported WC wave computation from v2 (equation #5 weighted partial standing wave + phasor superposition) to v3. Combined phasor: base wave phasor + WC phasors → total RMS → energy → force. Includes separation slider, phase toggle, Coulomb comparison, force annotations.

- ✅ Port WC class, separation slider, phase toggle, Coulomb comparison from v2
- ✅ Combined phasor: base wave phasor + WC phasors → total RMS → energy → force
- ✅ Added total energy display (ΣE) for energy conservation inspection

❌ **Finding: additive superposition does NOT change force behavior.** Adding WC waves on top of the base wave produces the same sinc oscillation as v2 (WCs in empty space). The base wave contributes a constant offset to the phasor, but the WC interference pattern — and therefore the force direction flips — is identical.

**Energy conservation test (Option B — normalized additive)**: scaling the combined RMS to preserve total energy (ΣE_combined = ΣE_base) conserves energy perfectly, but the normalization is a uniform scale factor. The spatial energy pattern is unchanged → gradients unchanged → forces unchanged → sinc flip unchanged. This confirms: energy-normalizing additive superposition is not sufficient. WCs need to **warp** the energy landscape, not just add waves on top.

⚠️ **Root cause**: additive superposition treats WCs as wave sources emitting into the field. More waves = more energy. The base wave is just "more field" that gets superposed — the WC interference pattern dominates regardless of the base wave's contribution. This is exactly the old model (equations #1–#6 from v2) running on top of a base wave, and it cannot resolve the force problem.

❌ **Key conclusion: additive superposition cannot work for Phase 1b, regardless of which base wave mode is underneath.** The WC wave pattern dominates and produces the same sinc oscillation for all 3 contenders (standing, quadrature, uniform). Simply adding a base wave on top does not change force behavior compared to v2 (WCs in empty space) — the base wave contributes a constant phasor offset that shifts the overall amplitude level but does not change the spatial interference pattern or force directions.

**Open question (partially answered)**: does the base wave change force behavior compared to v2? Additive: NO. But non-additive WC interaction (warping, disturbing, redistributing the base wave energy) is still untested — this is the core of Step 2c. The base wave may matter when WCs actually disturb it rather than just superpose on top of it.

---

### 🔶 Step 2c: Non-Additive WC Disturbance Models

Additive superposition (base + WC waves) is ruled out — it produces the same sinc oscillation regardless of base wave mode. WCs must warp the energy field non-additively (warping, disturbing, redistributing the base wave energy spatially).

For the base wave to matter, WCs must interact with it NON-additively — through reflection, scattering, or multiplicative modulation that changes the spatial structure of the energy landscape itself. Phase 1b requires WCs as energy **redistributors**, not energy **sources**.

**Per-contender WC interaction ideas** (to explore after migration):

- **Quadrature**: WCs tap into one channel direction (left/right traveling wave → charge sign?)
- **Standing wave**: WCs as boundary conditions or scattering centers within the standing wave field
- **Uniform**: add π-apart dual-channel wave, WCs disturb one channel depending on charge (→ Step 2c)

**WC disturbance mechanisms to investigate**: reflection/scattering of base wave, channel-selective disturbance, boundary condition at WC position, radial disturbance expansion, standing wave formation (nλ core, radius = K²λ), energy concentration near WC, energy deficit in far field

Three candidate approaches:

#### 🚧 Option A: Multiplicative (Energy Redistribution)

WC modifies base wave amplitude via a position-dependent multiplier:

```text
RMS_total(x) = RMS_base(x) · R(x)
R(x) > 1 near WC (energy concentration)
R(x) < 1 far from WC (energy drainage)
∫ R²(x) · E_base(x) dx = ∫ E_base(x) dx  (energy conservation)
```

The spatial pattern is fundamentally different from additive — no sinc, just a smooth concentration/depletion shape determined by R(x). R(x) could be derived from the WC's standing wave structure (K²λ radius, 1/r far-field decay). Energy-conserving by construction if R is properly normalized.

**Physical interpretation**: WC absorbs base wave energy from surrounding field and concentrates it into its own standing wave core. The absorption creates a far-field energy deficit — the drainage pattern that generates force on other WCs.

#### ❌ Option B: Normalized Additive — RULED OUT

Scale combined RMS to preserve total energy: `RMS_conserved = RMS_combined · √(ΣE_base / ΣE_combined)`. Energy conserves perfectly, but normalization is a uniform scale factor → spatial pattern unchanged → gradients unchanged → forces unchanged → sinc flip unchanged. Additive superposition with uniform normalization cannot change force behavior.

#### 🚧 Option C: Scattering (Reflection/Re-emission)

WC reflects the incoming base wave. Reflected wave interferes with the original base wave — but the reflected energy comes FROM the incident wave, not added on top. This is the LaFreniere model: WC reflects in-waves → out-waves. The interference between incident and reflected creates the energy redistribution.

**Key difference from additive**: in scattering, the incident wave is partially absorbed at the WC position (amplitude reduced), and the scattered wave carries that energy outward. Total energy conserved because scattered = incident - transmitted.

**Connection to Phase 1c**: scattering naturally produces variable λ near the WC (wavelength shifts from reflection), connecting to the non-linear wave equations path.

#### 🚧 Option D: Local Absorber (Boundary Condition)

WC acts as a boundary condition in the wave field — a point where displacement is constrained (e.g., pinned to zero, or to a specific amplitude). The wave field adjusts around the constraint, creating concentration and depletion patterns. Similar to how a fixed point in a vibrating string creates standing wave patterns around it.

**Physical interpretation**: the WC IS a point of wave reflection (Dirichlet BC = ψ=0 at WC position). Incoming waves reflect off the WC, creating standing waves near it and an amplitude shadow in the far field.

---

### 🚧 Step 2d: Dual-Channel Base Wave (π-apart)

Implement the dual-phase uniform model: two π-apart base waves that sum to zero energy. WCs disturb one phase or the other depending on charge sign. Requires non-additive WC interaction from Step 2c.

Dual-Channel Base Wave feeding WCs: uniform, standing, quadrature
maybe dual-chanel waves π apart, also applies to standing waves not only uniform waves. maybe there are two fundamental waves always canceling each other out but wave centers or matter can disturb that equilibrium and make the dual waves out of anti-phase and this is what manifests as energy

---

### 🚧 Step 2e: Deeper Physics Discussion and Contender Selection

**3 remaining contenders** (node-locking eliminated):

- **Quadrature** — strongest candidate (flat energy + directional charge encoding + spin/complex sinusoid connection). Traveling wave direction flips with temporal offset sign → possible charge/spin mechanism
- **Standing wave** — physically validated by Laplacian, simplest model, node structure may interact with WC phase at particle scale (100λ average-out argument)
- **Uniform + dual phase** — zero-energy vacuum from π-apart waves, WCs asymmetrically perturb one phase depending on charge

**Discussion topics**:

- Evaluate quadrature direction (left/right traveling wave) as charge/spin mechanism — does this have physical motivation in 3D? How does 1D two-direction map to 3D elliptical orientational freedom?
- Evaluate standing wave at particle scale (100λ) — does the node structure matter or average out?
- Compare force behavior across the 3 base wave contenders with non-additive WC interaction — which produces the best charge-dependent force direction? (Additive was identical across all 3 — non-additive may differentiate them)
- Compare non-additive WC interaction models (A, C, D) — which produces energy-conserving force with correct charge-dependent direction?
- Which combination of base wave mode + WC interaction mechanism best resolves the sinc oscillation blocker?

## M2 Research Prior-Art Findings
