# PHASE 1d: Non-Linear Wave Equations

Variable λ(r), ρ(x), Ψ³; breaks sinc periodicity while keeping genuine wave interference

**Rationale**: All linear operations on the sinc function preserve its λ/2 periodicity. Only a non-linear wave equation — where the spatial structure itself is no longer a pure sinc — can break the oscillatory pattern.

## Findings from Phase 1b Step 2c (hints for Phase 1d implementation)

Phase 1b tested elastic phase warping (Option F) — a charge-dependent phase shift `Δφ(r) = q · strength · δ(r)` applied to the base wave phasor at each WC. The phase warp is equivalent to local λ variation: `Δφ = ∫Δk(r')dr'`.

**Result**: near-zero forces. Phasor rotation preserves `|P² + Q²|` at every point → RMS unchanged → no energy gradient → no force via `F = -∇E`.

**Root cause**: the current energy formula `E = ρV(fA)²` uses constant f (and therefore constant λ). Phase warping changes the wave's spatial phase structure but the energy formula doesn't "see" λ variation — it only sees amplitude A. The phase disturbance is physically real (λ varies near WC), but the energy equation is blind to it.

**What Phase 1d needs to fix**: make λ a spatial variable in the energy equation: `E = ρV(c·A/λ(r))²`. When λ(r) varies near the WC (shorter = denser = more energy), the energy gradient becomes:

```text
∇E = ρVc² · ∇(A²/λ²) = ρVc² · [2A·∇A/λ² - 2A²·∇λ/λ³]
```

The second term `∇λ` is NEW — it creates force from wavelength gradients, independent of amplitude gradients. This is the term that Phase 1b's elastic phase warp couldn't produce because λ was constant in the energy formula.

**Implementation path**: use Yee & Hauger shells for λ(r) profile, WKB phase integral for phasor computation with variable k(r), and the expanded energy formula with λ(r). The phase warp code from Phase 1b (`_compute_elastic_phase_rms`) already computes the rotated phasor — just needs λ(r) in the energy formula to produce force from it.

**Convergence with Phase 1c**: Phases 1c and 1d converge. Phase 1c proved that spin alone (L→T conversion) creates magnetic force but NOT electric force (on-axis limitation). The electric force mechanism is variable λ(r) — this is Phase 1d's primary target. The vector infrastructure from 1c (Steps 1-2) is ready for integration.

---

## Findings from Phase 1c (inputs for Phase 1d)

Phase 1c Step 3 tested two WCs with spin at variable separations. Result: spin alone does NOT create charge-dependent Coulomb force. The T component is perpendicular to the WC axis → creates magnetic force, not electric. The sinc oscillation persists in the L component.

**Conclusion**: spin → magnetic force (Phase 4). Variable λ(r) → electric force (THIS phase).

---

## LaFreniere Phase Shift — Core Geometry Creates Charge

LaFreniere's analysis (`sa_phaseshift.html`) reveals that the electron core creates a **λ/2 phase shift** in waves passing through it:

- **Core diameter**: exactly 1λ (full wavelength)
- **Medium compression**: the core volume is 7x smaller than the first onion layer → medium is compressed → λ is shorter inside the core → frequency is higher (f = c/λ, c is constant)
- **Phase advance**: the compressed core has shorter λ → more phase cycles per unit distance → wave accumulates extra phase → emerges with a λ/2 phase shift — peaks become valleys
- **Charge sign**: this phase shift IS the charge. Electron and positron differ by this λ/2 shift (matter vs antimatter nodes in the standing wave)

**Clarification on "wave speed"**: LaFreniere describes waves as "faster" inside the core, using a sound wave analogy (sound is faster at higher air pressure). In EWT, **c is constant** — the speed of light is absolute and invariant. What LaFreniere calls "faster wave speed" is more precisely **higher frequency** from shorter wavelength: compressed medium → shorter λ(r) → higher f = c/λ → more phase cycles per unit distance. The wave doesn't travel faster — it oscillates faster in the compressed region. c stays constant throughout; only λ and ρ vary.

This is directly relevant to Phase 1d because:

1. The phase shift emerges from **variable λ(r)** inside the core (NOT variable c) → variable k(r) = 2π/λ(r)
2. The compressed medium means **variable ρ(r)** — higher density inside core
3. **c is constant** — only λ(r) and ρ(r) are spatial variables in `E = ρ(r)·V·(c·A/λ(r))²`
4. The ∇λ force term produces force from the wavelength gradient independently of amplitude oscillation

**Implementation implication**: the WC out-wave should use WKB phase integration `φ(r) = ∫k(r')dr'` instead of simple `kr`. The variable k(r) from Yee & Hauger shells encodes the phase shift naturally. The energy equation with variable λ(r) then "sees" the force from wavelength gradients that the constant-λ formula misses.

---

## Smoliński Insights for Phase 1d (from MagnetismGravity v4)

### r⁵ Energy Scaling — Three Physical Origins

Smoliński's `E ∝ r⁵` decomposition (Sec 7.2) provides the concrete non-linear relationships for our energy equation:

```text
E ∝ r⁵ = r³ × r¹ × r¹
         ↑     ↑     ↑
    volume  A∝r   f∝1/r
```

1. **Volumetric Occupancy (r³)**: 3D spatial extent of the WC within the BCC lattice
2. **Amplitude Scaling (r¹)**: charge = wave amplitude A in meters. Stability requires A to be geometrically coupled with radius: A ∝ r
3. **Resonance Frequency (r¹)**: intrinsic standing wave frequency scales inversely: f ∝ 1/r. As wavelength decreases, energy concentrates

These are the three spatial variables in our energy equation: V (volume), A (amplitude), and f = c/λ (frequency). All three vary with r near the WC core. The r⁵ scaling tells us exactly HOW they vary.

### Density Hierarchy — Three Levels

The transition from vacuum to gravitational force involves three density levels:

- **N_ν,max ≈ 10⁵⁴**: absolute geometric capacity of the BCC lattice (max packing)
- **N_ν,stat ≈ 10⁵²**: statutory background density (vacuum equilibrium, Eulerian dilution ~99.37%)
- **N_ν,eff ≈ 10⁴⁸**: effective gravitational density (matter-occupied region after push-out)

The 10⁻⁴² EM-to-gravitational ratio emerges from this hierarchy. For Phase 1d, the relevant insight is that ρ(r) varies from N_ν,stat (far field) to a displaced value near the WC core. This density variation contributes to ∇E through the ρ(r) term in `E = ρ(r)·V·(c·A/λ(r))²`.

### Soliton Push-Out — Gravity as Density Deficit

Gravity emerges as "volumetric density of the displacement deficit" — the WC displaces medium from its interior, creating a pressure deficit. This is not a separate force but a geometric consequence of the soliton structure within the high-density medium.

The Push-out Operator `P̂Φ = -∇·(η_stat/η_soliton)∇Φ` formalizes this: force from gradient of potential weighted by local density mismatch. This IS our `F = -∇E` with variable ρ(x) — Smoliński provides the explicit density profile.

### Geometric Soliton Core — Fine Structure Decomposition

The inverse fine structure constant decomposes as (Eq. 5):

```text
1/α = (4π³ + π² + π) - ε_M + Σε_G
       ↑ Geometric Core     ↑ Magnetic  ↑ Gravitational
       (matter/charge)       (spin)       (deficit × K)
```

This confirms the force hierarchy: the geometric core creates the baseline (electric), ε_M creates the magnetic correction (spin energy cost = 1/(N·π³)), and Σε_G = K_WC · ε_G creates the gravitational correction (accumulated deficit from K=10 wave centers).

### ε_M = 1/(8π⁷) — The Lattice Response Parameter

ε_M is NOT a fitting parameter — it's the intrinsic "stiffness" of the BCC vacuum lattice. It represents the structural resistance to deviation from ideal spherical wave symmetry. This single parameter bridges G, α, and lepton anomalous magnetic moments. For Phase 1d: ε_M defines the non-linear elasticity coefficient in the soliton wave equation.

**Smoliński's Non-linear Soliton Wave Equation** (MagnetismGravity_v2, Sec 6.1, Eq. 18-19):

```text
(∂²/∂t² - c²∇²) Ψ(r,t) + F(Ψ, ε_G, |ε_M|, N_ν) = 0

where F = k(|ε_M|) · Ψ³    (NLS cubic non-linearity)
```

The Ψ³ term prevents wave dispersion and creates soliton stability. The coefficient k(|ε_M|) describes the non-linear elasticity of the medium, depending explicitly on the magnetic deficit. Solutions are NOT pure sin(kr)/kr — the soliton spatial structure differs from a sinc, potentially breaking the periodic zero pattern. This is a well-studied form (Non-linear Schrödinger solitons) with known analytical solutions.

**Three energy gradient variables** — currently only A varies spatially; making ρ and f spatial variables turns E = ρV(fA)² into a multi-variable field:

- **A** (amplitude): current phasor RMS — carries sinc oscillation
- **f / λ** (frequency): Yee & Hauger discrete wavelength shells give λ(r) = 2(K-n)λ per shell. WKB phase integral replaces kr with ∫k(r')dr', breaking sinc periodicity. **Smoliński r⁵ decomposition**: E ∝ r⁵ = r³ (volume) × r¹ (A ∝ r) × r¹ (f ∝ 1/r)
- **ρ** (density): granule velocity determines local medium density/pressure. **Smoliński density function**: `ρ(r) = ρ₀(1 - (r/r_ν)^k)^P · Θ(r_ν - r)`. ∇ρ may carry force information that ∇A alone cannot provide

Computing F = -∇E means these additional variables are automatically included when implemented — no force logic changes needed.

**Smoliński's Push-out Operator** (Eq. 90): `P̂Φ = -∇·(η_stat/η_soliton)∇Φ` — force from gradient of potential weighted by local density mismatch, formalizing F = -∇E with variable ρ(x).

**Connection to dual-treatment boundary**: Smoliński's Isotropy Operator acts as a **geometric low-pass filter** at the Degraded EMC Wall boundary. Inside → non-linear toroidal dynamics (r⁵), outside → isotropic spherical push-out (r³).

**Connection to wave resonance**: variable λ(r) means different regions oscillate at different frequencies. Energy exchange between regions occurs through **wave resonance** — the mechanism by which A and λ convert into each other while conserving energy. Wolff's coupled-oscillator model describes how two wave systems readjust frequencies through their shared medium. This may be the physical process behind the non-linear energy redistribution that creates force gradients. See [Wave Resonance](06_time_dynamics.md#wave-resonance-as-the-fundamental-energy-exchange-mechanism) for the full analysis.

---

## Yee & Hauger — Variable λ(r) Profile (from Spin.pdf / Lambda.pdf)

The key equations for implementing variable wavelength in Phase 1d:

```text
r_core = Kλ                           (core radius, K = wave center count)
r_wavelength(n) = 2(K - n)λ           (width of n-th wavelength shell)
r_x = (K + 2·Σ(K-n), n=1..x) · λ     (cumulative radius to x-th shell)
r_particle = K²λ                       (particle radius = standing wave boundary)
```

For electron (K=10): core = 10λ, particle radius = 100λ. The wavelength shells are NOT uniform — they're widest near the core (n=1: width = 18λ) and shrink toward the boundary (n=K: width = 0). This non-uniform spacing means:

1. **k(r) varies with r**: wavenumber is NOT constant inside the particle
2. **Phasor phase uses WKB integral**: φ(r) = ∫₀ʳ k(r')dr' instead of simple kr
3. **The sinc function is replaced** by a non-uniform spatial structure
4. **Node positions are non-uniformly spaced** — this is why the tetrahedral electron can't sit on a uniform λ/2 lattice

For a single WC (K=1, neutrino): core = 1λ, particle radius = 1λ. Only one shell. The wavelength variation is minimal. This is consistent with the neutrino being the simplest soliton.

### ⚠️ λ(r) Direction: Yee & Hauger vs LaFreniere — Open Question

There is a **contradiction** between the two models on how λ varies near the WC:

| Model | λ near core | λ toward boundary | Mechanism |
| --- | --- | --- | --- |
| **Yee & Hauger** | LONGEST (shell n=1 = 18λ₀) | SHORTEST (shell n=K → 0) | Steepness conservation: A/λ = const, high A near core → long λ |
| **LaFreniere** | SHORTEST (compressed core) | returns to λ₀ | Medium compression: 7x smaller volume → higher density → shorter λ → higher f |
| **Smoliński** | f ∝ 1/r → longer λ at large r | f ∝ 1/r → shorter λ at small r | r⁵ decomposition: resonance frequency inversely proportional to radius |

Smoliński's f ∝ 1/r agrees with LaFreniere (shorter λ / higher f closer to core) but contradicts Yee & Hauger (longer λ near core). Possible resolutions:

1. **Different regions**: LaFreniere and Smoliński may describe INSIDE the core (r < Kλ, where medium is physically compressed). Yee & Hauger describe the standing wave SHELLS outside the core (Kλ < r < K²λ). Both could be correct in their respective domains
2. **Shell width ≠ local wavelength**: Yee & Hauger's "wavelength shell width" may represent the spatial extent of one wave cycle, not the local wavelength of the medium. A wide shell with high amplitude could still have short local oscillation wavelength
3. **One model is wrong**: the phase shift and force behavior in our simulation will reveal which λ(r) profile produces correct physics

**Phase 1d must test both profiles** — implement the Yee & Hauger shell model AND a LaFreniere-style compression model, and compare the force behavior. The correct profile should produce the λ/2 phase shift that creates charge sign and the ∇λ force term that produces Coulomb force.

---

## Butto Vortex Electron (reference only — future Phase 4)

Butto's model (2021) describes the electron as a superfluid irrotational vortex:

- **Spin-½ = differential rotation**: core rotates at 2x boundary rate (720° core / 360° boundary)
- **Planck constant = vortex angular momentum**: Γ·m_e = h
- **Magnetic moment without g-factor**: μ = qcr from vortex hydrodynamics
- **Helmholtz theorems**: vortex filaments form closed loops → toroidal geometry

Relevance to Phase 1d: the vortex model provides a physical picture for HOW the medium moves inside the soliton core — toroidal circulation that naturally creates the density variation and wavelength compression. The differential rotation (core 2x boundary) may explain why the core has a different effective λ than the surrounding shells. This connects to Smoliński's Energy Domain (toroidal, r⁵) vs EMC Domain (spherical, r³).

Not directly implementable in Phase 1d math scripts, but provides physical intuition for the variable-λ profile.
