# PHASE 1d: Non-Linear Wave Equations

Variable λ(r), ρ(x), Ψ³; breaks sinc periodicity while keeping genuine wave interference

**Rationale**: All linear operations on the sinc function preserve its λ/2 periodicity. Only a non-linear wave equation — where the spatial structure itself is no longer a pure sinc — can break the oscillatory pattern.

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

**Convergence with Phase 1d**: Phases 1c and 1d may converge into one solution — vector displacement with L→T conversion (1c) + variable λ(r) in the energy equation (1d). The λ gradient creates force from wavelength variation, while the L→T spin creates charge-dependent asymmetry. Both are needed: λ variation alone has no charge sensitivity, L→T alone still oscillates in scalar RMS. Together they could produce charge-dependent force from wavelength + mode gradients.
