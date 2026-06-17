# ✅ Phase 1d: Non-Linear Wave Equations

Variable λ(r), ρ(x), Ψ³; breaks sinc periodicity while keeping genuine wave interference

**Rationale**: All linear operations on the sinc function preserve its λ/2 periodicity. Only a non-linear wave equation — where the spatial structure itself is no longer a pure sinc — can break the oscillatory pattern.

---

## 1D Variable-λ Test Results (`step1d_variable_lambda.py`)

Tested two WCs on 1D axis with Yee & Hauger λ(r) profile and variable-λ energy equation `E = ρV(c·A/λ(r))²`. Compared K=1 (neutrino) and K=10 (electron) across separations 1.5× to 4× particle radius, four phase configurations.

### K=1 (neutrino): No λ variation → identical to constant-λ

- λ(r) = λ₀ everywhere (K=1 has only 1 shell, no variation)
- WKB phase = constant-k phase (Δφ = 0)
- Force directions: ALL MIXED, var-λ = const-λ exactly
- **Correct**: neutrino is neutral, sinc persists at K=1

### K=10 (electron): ∇λ term IS active but charge-blind

- λ(r) varies dramatically: core λ = 10λ₀, first shell = 19.5λ₀, shrinking to λ₀ at boundary
- WKB phase deficit: Δφ = -182π at particle boundary (massive non-linear phase accumulation)
- **var-λ and const-λ DIVERGE** at far-field separations (250λ, 400λ) — the ∇λ term reverses force direction compared to constant-λ. This confirms the ∇λ mechanism works.
- **BUT: force is charge-blind** — ALL four phase configs (same/opposite) get the same direction at each separation. The λ(r) profile depends on K (structure), not on source_offset (charge). The ∇λ contribution creates the same force for same and opposite phase.

### Conclusion: ∇λ works but is charge-blind in 1D

| Finding | Status |
| --- | --- |
| ∇λ changes force vs constant-λ | ✅ confirmed at K=10 |
| K=1 vs K=10 behave differently | ✅ confirmed (K=1 no effect, K=10 significant) |
| Charge-dependent force from ∇λ | ❌ charge-blind (same direction for all phases) |
| 1D sinc still dominates charge physics | ❌ sinc oscillation persists along axis |

---

## ⚠️ Critical Insight: 3D Spherical Integration May Break Sinc

The 1D test only captures the ON-AXIS force. But 3D spherical wave interference creates energy gradients in ALL directions (visible in M1 simulation screenshots — elliptical/hyperbolic interference patterns extending off-axis between two WCs).

The force on a WC is `F = -∇E` integrated over the **full 3D volume**, not just the connecting axis. Off-axis contributions have different oscillation periods:

```text
On-axis (θ=0):   path difference = d          → oscillation period λ/2
At angle θ:      path difference = d·cos(θ)   → period λ/(2·cos(θ)) — LONGER
Perpendicular:   path difference ≈ 0          → no oscillation
```

When integrating -∇E over all solid angles, these different periods **average out the sinc flips**. The on-axis λ/2 oscillation is counteracted by non-oscillating off-axis contributions. The net force could be smooth even though the on-axis slice oscillates.

This connects to:

- **LaFreniere's diffractive lens effect**: standing wave zones between particles focus energy from a wide solid angle onto the axis, amplifying the net force through 3D constructive interference
- **Smoliński's Degraded EMC Wall**: the geometric low-pass filter that smooths oscillatory structure from the discrete lattice into isotropic gravity — the 3D spherical averaging IS this filter
- **Phase 1c on-axis limitation**: the L→T spin test (Step 3) also failed because it only measured on-axis force. The T component creates off-axis gradients that project onto the axis differently

**The Coulomb force may be an inherently 3D phenomenon** — it emerges from the spherical integration of energy gradients, not from any 1D slice. The sinc oscillation IS real on-axis, but the NET force after 3D integration may be smooth and charge-dependent.

### 2D Cross-Section Test (`step1d_variable_lambda_2d.py`)

Tested with 1024² grid — K=1 at 64 vox/λ (excellent), K=10 at 4 vox/λ (adequate).

**Result: ALL separations are CHARGE-BLIND.** Same-phase and opposite-phase always get the same force direction. 2D off-axis averaging did NOT produce charge-dependent force.

⚠️ Critical observation: the force on WC1 depends on **WC1's own phase relative to the base wave** — NOT on WC2's phase. For example at K=10, sep=105λ: `same(0,0)=REP, same(π,π)=ATT, opp(0,π)=REP, opp(π,0)=ATT` — the direction tracks WC1's phase (0→REP, π→ATT), ignoring WC2 entirely.

**Root cause**: the base wave × WC interaction terms in E_total overwhelm the WC × WC cross term. The charge-dependent force lives in the **interaction energy** `E_int ∝ Re(P_wc1* · P_wc2)`, but `F = -∇E_total` is dominated by the much larger `∇(P_base · P_wc1)` and `∇(P_base · P_wc2)` terms. The WC-WC signal is buried.

### What We've Systematically Eliminated

| Mechanism | Test | Result |
| --- | --- | --- |
| Spin alone (L→T) | Phase 1c Step 3 (3D) | Magnetic only, no Coulomb |
| Variable λ(r) alone | Phase 1d 1D | ∇λ active but charge-blind |
| 2D off-axis averaging | Phase 1d 2D | Charge-blind, base wave dominates |
| Variable λ + 2D | Phase 1d 2D | Still charge-blind |

### What Remains: The Interaction Energy

The smooth envelope model from Phase 1 already showed that the WC-WC interaction energy `E_int ∝ |Z₁|·|Z₂|·cos(k·Δr + Δφ)` produces correct force direction (17/17) when extracted directly — but with charge imposed via ±1 label. The challenge has always been making this EMERGENT.

The interaction cross term DOES depend on relative phase (Δφ = source_offset difference). The problem is:

1. The cos(k·Δr) factor creates sinc oscillation in the interaction energy itself
2. `F = -∇E_total` mixes the interaction signal with the much larger base wave × WC terms

### Isolated Interaction Force Test (`step1d_analytical_force.py`)

Isolated the WC×WC cross term: `E_int = E(wc1+wc2) - E(wc1_alone) - E(wc2_alone)`, removing base wave and self-energy contamination. 2D grid, no base wave, only two WC out-waves.

**K=1 results**: PERFECTLY charge-dependent, PERFECTLY consistent:

```text
same (0,0): 10/10 ATT — CONSISTENT
same (π,π): 10/10 ATT — CONSISTENT
opp  (0,π): 10/10 REP — CONSISTENT
opp  (π,0): 10/10 REP — CONSISTENT
```

**K=10 results**: charge-dependent at ALL separations. Mostly consistent (7/9 ATT for same, 7/9 REP for opposite), with flips only at 105λ and 115λ (inside particle radius, near-field).

### ⚠️ Critical Assessment: What We Actually Showed

**The good:**

- The WC×WC interaction cross term IS inherently charge-dependent — same vs opposite phase produce opposite force directions
- 2D integration (off-axis contributions) produces consistent direction — no sinc flips at most separations
- The interaction term contains the charge-dependent physics we're looking for

**The problems:**

1. **Same charge ATTRACTS (180° off Coulomb)**: same phase (0,0) → ATT, opposite phase (0,π) → REP. This is the OPPOSITE of Coulomb (same should repel, opposite should attract). This may be the strong force / lock-in capture mechanism (LaFreniere: "external radiation pressure produces an attraction effect" for same-phase lock-in), not the Coulomb force.

2. **We manually subtracted self-energies**: `E_int = E(1+2) - E(1) - E(2)` is a mathematical decomposition, analogous to renormalization. The real wave field includes ALL terms. A particle in the real universe doesn't get to subtract its own self-energy — it feels `F = -∇E_total`. The question is: what mechanism in nature performs this "subtraction"? Possible answer: the WC's own field IS symmetric around itself → ∇E_self = 0 at its center by symmetry → F_self = 0 naturally. The problem was the base wave × WC term overwhelming the signal, not the self-energy.

3. **Base wave × WC dominates in E_total**: when we include the base wave, the WC-base interaction overwhelms the WC-WC interaction. The force on WC1 depends on WC1's phase relative to the base wave, not on WC2's charge. This means the Coulomb force is REAL but much weaker than the WC-base interaction at our simulation scale.

### Honest Status of Force Emergence

| Effect | Status | Evidence |
| --- | --- | --- |
| Strong force lock-in | ✅ Emerges | Sinc nodes create energy wells at λ/2 for same-phase WCs |
| Annihilation | ✅ Emerges | Opposite phase: deepest well at r=0, barriers at λ/2 explain positronium |
| Particle formation | ✅ Concept works | Lock-in holds WCs, but M3 tetrahedron unstable (needs variable λ for non-uniform nodes) |
| Charge-dependent interaction | ✅ Exists in cross term | Isolated E_int shows 100% charge-dependent direction (K=1) |
| Coulomb direction (correct sign) | ❌ Wrong sign | Same→ATT, opposite→REP. 180° opposite to Coulomb. May be strong force capture, not Coulomb |
| Coulomb from E_total (no subtraction) | ❌ Not yet | Base wave × WC terms overwhelm the WC×WC cross term in F = -∇E_total |

### Next: Standing Wave vs Traveling Wave Decomposition

The M3 wave engine already decomposes each WC's out-wave into standing + traveling components (weighted partial standing wave equation):

```text
ψ = A · [w·sin(kr+ωt+φ) + sin(kr-ωt-φ)] / kr
         ↑ in-wave (w→1 near)     ↑ out-wave (always present)

Phasor decomposition:
  C_n = (w+1)·sin(kr)/kr    ← STANDING component (dominates near WC)
  S_n = (w-1)·cos(kr)/kr    ← TRAVELING component (dominates far from WC)
```

**LaFreniere's key distinction**: standing waves → strong force (lock-in, capture). Traveling waves → electrostatic force (Coulomb). These are different force mechanisms from different wave components.

The standing wave component has sinc nodes `sin(kr)/kr` → lock-in at λ/2. The traveling wave component `sin(kr-ωt)/kr` carries energy flux outward → radiation pressure → Coulomb.

**Hypothesis**: if we compute force from ONLY the traveling wave component of two WC interactions, we may get the correct Coulomb behavior (same repels, opposite attracts) without sinc flips.

### In-Wave / Out-Wave Decomposition Results (`step1d_standing_traveling.py`)

Properly decomposed the weighted partial standing wave into in-wave (`e^{-ikr}`, inward) and out-wave (`e^{+ikr}`, outward) phasors. Computed isolated interaction energy and force from each component separately, plus the full wave.

**Result: ALL three components (in-wave, out-wave, full) produce the SAME λ/2 oscillation pattern.** The direction flips every 0.5λ of separation — alternating between "COULOMB CORRECT" and "INVERTED" at each half-wavelength step. No component is free of the sinc.

| Component | Force magnitude (7λ sep) | Charge-dependent? | Consistent? |
| --- | --- | --- | --- |
| In-wave | 8.9×10⁻¹² (negligible at far-field) | Yes | MIXED (λ/2 flips) |
| Out-wave | 8.6×10⁻⁶ (dominant) | Yes | MIXED (λ/2 flips) |
| Full | 1.6×10⁻⁵ | Yes | MIXED (λ/2 flips) |

**Definitive finding**: the sinc oscillation is in the OUT-WAVE itself. The pure traveling wave `sin(kr-ωt)/kr` has `sin(kr)/kr` spatial structure — the `1/kr` sinc envelope IS the traveling wave's radial profile for a monochromatic spherical source. There is no wave decomposition that removes the sinc — it's intrinsic to spherical wave propagation at a single frequency.

### Where We Stand: The Sinc Is the Physics

The sinc oscillation cos(k·Δr + Δφ) in the interaction energy is NOT a numerical artifact, model limitation, or decomposition failure. It is the **fundamental mathematical consequence of two coherent monochromatic spherical waves interfering in space.** Every approach that uses:

1. Single frequency ω (monochromatic)
2. Spherical wave propagation (1/r decay)
3. Coherent superposition (amplitudes add)

...will ALWAYS produce cos(k·Δr) in the interaction term, regardless of:

- Variable λ(r) (changes node positions but not the far-field periodicity)
- L→T spin conversion (creates perpendicular force, not radial)
- 2D/3D integration (sinc persists off-axis too for isolated interaction)
- Standing/traveling decomposition (both have sinc structure)
- Self-energy subtraction (reveals the cross term, which oscillates)

### ⚠️ What Must Be Different for Coulomb

For the Coulomb force to emerge without sinc flips, the interaction must violate at least one of the three conditions above:

1. **Non-monochromatic**: if the out-wave has a wavelength SPREAD (not single λ), the cos(k·Δr) terms at different k values partially cancel, smoothing the oscillation. This connects to the Yee & Hauger shells — if different shells emit at different effective wavelengths, the far-field is broadband. This is also the stochastic base wave concept from Phase 1b.

2. **Non-coherent superposition**: if particles are not phase-locked to each other (random phase jitter, thermal fluctuations), the cos(Δφ) term averages to zero for same-phase but survives for opposite-phase through a different mechanism.

3. **Not pure spherical propagation**: if the WC modifies the wavefront shape (e.g., toroidal/elliptical from spin, diffractive lens focusing), the 1/kr sinc might not apply.

4. **Force from energy FLUX, not energy DENSITY**: LaFreniere's radiation pressure model computes force from wave momentum transfer (Poynting-like flux), not from -∇E. Flux-based force may behave differently from gradient-based force.

5. **The force is statistical/averaged**: maybe the Coulomb force IS the time-averaged or ensemble-averaged sinc, where the averaging window is the particle radius K²λ. At K=10 (100λ radius), averaging over 100 wavelengths would capture ~200 sinc oscillation cycles, potentially averaging them out.

These are the remaining avenues. Each represents a fundamentally different approach from what we've tried.

---

### Avenue #1: Broadband from Yee & Hauger Shells (`step1d_broadband.py`) — ❌ ELIMINATED

Hypothesis: each shell emits at different λ → far-field is broadband → cos(k_n·Δr) at different k_n partially cancel.

K=10 has 9 shells at λ/λ₀ = [18, 16, 14, 12, 10, 8, 6, 4, 2]. Each contributes an out-wave at its own wavenumber.

**Result**: ALL K values MIXED. Multiple cosines at different frequencies create complex beating, not convergence. K=10 broadband: 28/50 ATT, 22/50 REP — still MIXED.

### Avenue #4: Flux-Based Force — Radiation Pressure

#### 1D Flux Test (`step1d_flux_force.py`) — ⚠️ PARTIAL SUCCESS

Hypothesis: force from energy FLUX `S = -c²·ψ·∂ψ/∂x` (wave momentum transfer) instead of energy gradient `F = -∇E`.

Time-averaged over one full period (32 samples). Two flux definitions tested: spatial `S = -c²·ψ·∂ψ/∂x` and temporal `S = ψ·∂ψ/∂t`.

**K=1**: MIXED for all methods.

**K=10**: Flux(spatial) same charge = **50/50 REP → CONSISTENT → COULOMB CORRECT** — first ever consistent result for any configuration at K=10. But opposite charge = 49/50 REP (should be ATT) — the baseline outward push from radiation pressure dominates over the charge-dependent correction in 1D.

| Method | K=10 same (expect REP) | K=10 opp (expect ATT) |
| --- | --- | --- |
| Gradient (-∇E) | 25/50 MIXED | 25/50 MIXED |
| Flux (spatial) | **50/50 REP CONSISTENT** | 49/50 REP (wrong) |
| Flux (temporal) | 49/50 REP | 48/50 REP |

The flux mechanism extracts consistent direction where the gradient cannot — but only for same charge. Opposite charge needs the 2D off-axis inward pressure (LaFreniere's diffractive lens).

#### 2D Flux Test (`step1d_flux_force_2d.py`) — ⚠️ KEY FINDING

Full 2D time-domain computation: weighted partial standing wave on 512² grid (32 vox/λ), time-averaged flux through a circle surrounding WC1.

**Breakthrough: 100% charge discrimination at ALL separations.** Same and opposite phase ALWAYS produce OPPOSITE force directions at every separation tested (22/22):

```text
sep=2.00λ:  same=ATT  opp=REP  → OPPOSITE
sep=2.50λ:  same=REP  opp=ATT  → OPPOSITE  ← COULOMB CORRECT
sep=3.00λ:  same=ATT  opp=REP  → OPPOSITE
sep=3.50λ:  same=REP  opp=ATT  → OPPOSITE  ← COULOMB CORRECT
...every separation: same ≠ opp
```

**The sinc oscillation persists** — the ABSOLUTE direction flips every λ/2. But the RELATIVE direction (same vs opposite) is ALWAYS opposite. At half-integer λ separations (2.5, 3.5, 4.5...): Coulomb correct (same=REP, opp=ATT). At integer λ (2, 3, 4...): inverted.

**Key properties of the 2D flux result:**

- **Charge discrimination: 100%** — same and opposite ALWAYS get opposite directions
- **Sinc persists**: absolute direction flips every λ/2 for both
- **Equilibrium points**: at quarter-λ offsets, forces are ~10x smaller — lock-in equilibrium (LaFreniere "capture")
- **Force decays with distance**: ~1/r trend under the oscillation
- **Coulomb correct at half-integer λ separations**: same=REP, opp=ATT

This IS the wave interference physics: the sinc oscillation produces alternating lock-in and Coulomb zones. The full 2D radiation pressure correctly captures the charge discrimination that the energy gradient approach misses.

### Avenue #5: Statistical Averaging Over Particle Radius (`step1d_averaged_force.py`) — ❌ ELIMINATED

Hypothesis: average the sinc interaction force over K²λ (particle radius). K=1 → 1λ window → no smoothing → neutral. K=10 → 100λ window → 200 sinc cycles → smooth → Coulomb.

**Result**: ALL K values (1, 2, 3, 5, 10) give 25/50 ATT, 25/50 REP — perfect 50/50 split. The sinc is perfectly symmetric: each half-cycle of attraction is exactly balanced by repulsion. Averaging over many cycles gives ~zero, not a net direction.

### Phase 1 — Mechanism Status After All Tests

| # | Avenue | Script | Result |
| --- | --- | --- | --- |
| — | Spin L→T (1c) | step3_two_wc_force.py | Magnetic only, not electric |
| — | Variable λ (1d, 1D) | step1d_variable_lambda.py | ∇λ active but charge-blind |
| — | Variable λ (1d, 2D) | step1d_variable_lambda_2d.py | Charge-blind, base wave dominates |
| — | Variable λ (1d, 3D) | step1d_variable_lambda_3d.py | K=10 resolution too low |
| — | Isolated interaction | step1d_analytical_force.py | Charge-dependent but wrong sign (ATT for same) |
| — | In/out wave decomp | step1d_standing_traveling.py | Sinc in ALL components |
| 1 | Broadband shells | step1d_broadband.py | ❌ Multiple cosines still oscillate |
| 4 | 1D flux | step1d_flux_force.py | ⚠️ Same-charge consistent at K=10 |
| 4 | **2D flux** | **step1d_flux_force_2d.py** | **⚠️ 100% charge discrimination, sinc persists** |
| 5 | Statistical averaging | step1d_averaged_force.py | ❌ Sinc symmetric → zero |

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

### Connection to Medium Pressure and Gravity — λ Unifies All Three

Changing λ/frequency changes the rate at which granules cycle through their elliptical motion — meaning different granule velocities. In wave mechanics and particle-fluid simulations (e.g., sound waves), granule velocity is directly related to **localized medium pressure**. Faster cycling = higher pressure, slower cycling = lower pressure.

This creates a direct link to Smoliński's gravitational push-out/buoyancy model: if wave centers alter λ in their vicinity (time dilation near mass), they change the local granule velocity, which changes the local medium pressure/density. The pressure deficit around a massive body IS the gravitational field — and it emerges from λ modulation, not from a separate gravitational mechanism. Gravity, time dilation, and medium pressure become three descriptions of the same underlying phenomenon: **local λ variation**.

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

**Connection to wave resonance**: variable λ(r) means different regions oscillate at different frequencies. Energy exchange between regions occurs through **wave resonance** — the mechanism by which A and λ convert into each other while conserving energy. Wolff's coupled-oscillator model describes how two wave systems readjust frequencies through their shared medium. This may be the physical process behind the non-linear energy redistribution that creates force gradients.

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
