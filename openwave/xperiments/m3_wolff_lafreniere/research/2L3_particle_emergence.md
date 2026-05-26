# LAYERS 1-3: PARTICLE EMERGENCE

Near-field standing wave physics — demonstrating that particles emerge from wave interference.

Task checklists in [0_ROADMAP.md](0_ROADMAP.md). Energy layers hierarchy in [0_OVERVIEW.md](0_OVERVIEW.md).

---

## WAVE EQUATIONS TESTED

Five wave equations tested in M3 for lock-in, annihilation, and K=10 tetrahedron stability. All are exact solutions to the 3D wave equation in spherical coordinates.

### Combined Wolff-LaFreniere (current best)

Phase convention: `spatial ± (temporal + offset)` → `(kr ± (ωt + φ))`

```text
Standard form (product):
  ψ(r,t) = 2A · sin(kr/2) · cos(kr/2 - (ωt + φ)) / r

  sin(kr/2) / r       = spatial envelope (sinc zeros at λ, 2λ... = particle structure)
  cos(kr/2 - (ωt+φ))  = traveling wave at half-wavenumber (force field propagation)
  2A                   = standing wave doubling (in + out constructive)

Phasor (expand product back to Phase + Quadrature):
  sin(kr/2)·cos(kr/2-(ωt+φ)) = ½sin(kr)·cos(ωt+φ) + ½(1-cos(kr))·sin(ωt+φ)
  C_n = A·sin(kr)/r            (Phase:      zeros at λ/2, λ, 3λ/2...)
  S_n = A·(1-cos(kr))/r        (Quadrature: zeros at λ, 2λ, 3λ...)
  |phasor| = 2·|sin(kr/2)| / r (using 1-cos=2sin²(x/2))

Center limit (r→0):
  ψ(0,t) = A·k·cos(ωt+φ)
```

The product form is a **modulated traveling wave**: the particle (envelope) shapes the outgoing force field (traveling oscillation). One formula, one physical picture.

**Key property**: sinc zeros at **λ intervals** (not λ/2). The λ/2 lock-in distance is at the **peak** of the envelope (sin(π/2)=1), not at a node. All 45 K=10 tetrahedron WC pairs sit within [0, λ) where there are no zeros.

**Quadrature term**: C_n and S_n have DIFFERENT spatial functions (sin(kr)/r vs (1-cos(kr))/r). The ratio |S_n/C_n| = |tan(kr/2)| varies with distance — encoding position-dependent elliptical oscillation. At λ/4: circular. At λ/2: purely quadrature. The product form absorbs this into the `kr/2` spatial phase of the traveling wave — the ellipticity is encoded as position-dependent phase, not lost.

### Equation comparison

| Property | Wolff-original | Combined W-L | Weighted PSW |
| --- | --- | --- | --- |
| Center amplitude | k√2 ≈ 0.67 | k ≈ 0.48 | 2.0 |
| Envelope zeros at | λ/2, λ, 3λ/2... | **λ, 2λ, 3λ...** | λ/2, λ, 3λ/2... |
| Well spacing | λ/2 | **λ** | λ/2 |
| First barrier | λ/2 from WC | **λ from WC** | λ/2 from WC |
| Quadrature | same spatial fn | **different fns** | weight-dependent |
| K=10 at perfect placement | unstable (zero at λ/2) | **holds** (peak at λ/2) | unstable (zero at λ/2) |

### What IS grounded in wave physics

- **`sin(kr ± ωt) / r`** — exact solution to 3D wave equation in spherical coordinates. Textbook.
- **Superposition principle** — adding solutions gives another solution. Textbook.
- **In-wave + out-wave** — Wolff's Space Resonance model. Huygens-Fresnel principle.
- **Standing wave at core** — the particle IS the standing wave. Central claim of Wolff/EWT.
- **Sinc lock-in** — energy wells from wave nodes. Direct consequence of wave equation.
- **Annihilation** — opposite-phase cancellation. Direct consequence of superposition.

### Why analytical equations, not PDE solvers

Validated in 1D sandbox (`wave_engine_1D_v3.py`, Laplacian mode + absorber disturbance):

**PDE/Laplacian solvers cannot simulate wave-center disturbance in a wave field.** A wave center modeled as a Dirichlet boundary condition (ψ=0 reflective point) is invisible to the isotropic field — the wave passes through without creating standing waves around it. The PDE solver treats the boundary as a passive constraint, not an active disturbance that reorganizes the wave structure.

This approach also cannot simulate the in-wave + out-wave superposition that creates the standing wave particle. The in-wave (Huygens reconstruction from the universe) and out-wave (re-emitted from the WC) must be explicitly modeled as counter-propagating analytical solutions — the PDE solver only propagates what's already in the field.

The Combined Wolff-LaFreniere analytical solution solves this: Wolff's phase term `sin(kr)/r` (in-wave/out-wave standing wave) combined with LaFreniere's quadrature term `(1-cos(kr))/r` (traveling wave component) creates the standing wave field with energy wells that trap other standing waves at wave nodes — the mechanism for particle formation from wave interference.

### What's a modeling choice (not derived)

1. **The wave equation form** — Combined W-L selected empirically (best K=10 stability at perfect placement). Physical justification for this specific combination not yet established.
1. **Monochromatic** (single frequency) — root cause of sinc barrier issue at far-field.
1. **The 1.25λ transition distance** for signed envelope blend — modeling parameter.

---

## SIGNED ENVELOPE — Energy & Force

Force computed from `F = -∇E` where `E = ρ·V·(f·A)²` and A is the signed envelope.

### Combined W-L envelope (hybrid near/far)

```text
Near-field:  charge_sign × A_eff × 2|sin(kr/2)| / r   (sinc wells at λ)
Far-field:   charge_sign × A_eff × 2 / r               (smooth, no barriers)
Transition:  blended via weight = 1 / (1 + (r/1.25λ)^8)
Center:      charge_sign × A_eff × k                   (≈ 0.48, finite)
```

### Charge sign — imposed, not emergent

`charge_sign = cos(source_offset)` gives ±1 (imposed). Same sign → constructive → repulsion. Opposite → destructive → attraction. True Coulomb from spin deferred to Layer 4 / M4.

---

## AMPLITUDE PROFILE: Sinc Steepness & Energy Conservation

```text
E_shell ∝ r² × (sin(kr)/kr)² = sin²(kr)/k²
```

Shell surface area (r²) exactly cancels 1/r² from sinc → each shell has roughly equal energy. The sinc function `sin(kr)/kr` (spherical Bessel j₀) stays finite at origin (→ 1). Variable λ(r) from Yee & Hauger would make profile less steep. Deferred.

---

## WHAT WE DEMONSTRATED

### Lock-in (same-phase WCs)

- 2 WCs at λ separation oscillate around equilibrium — visually demonstrated
- K=10 (1-3-6 tetrahedron) holds with Combined W-L at perfect placement
- Validated: rotation-invariant, translation-invariant, phase-invariant (e⁻ and e⁺ both stable)

### Annihilation (opposite-phase WCs)

- Requires INIT_VELOCITY (momentum injection) — phasor cross-term barriers block approach, signed envelope has no far-field attraction
- Annihilation threshold at 1.2λ (standing wave core overlap)
- Both head-on (annihilation1) and diagonal (annihilation2) demonstrated
- Far-field Coulomb attraction deferred to Layer 4 / M4

### K=10 tetrahedron — perfect placement only

- Combined W-L: all 45 pair distances in [0, λ) → maximum envelope → holds
- Other equations: sinc zero at λ/2 → pairs at λ/2 have zero lock-in force → unstable

---

## WHAT WE DID NOT ACHIEVE — K-Selectivity

### The critical test

Perturb WC positions by ±10-30% of λ. If K=10 recovers and K=2..9 don't → genuine wave physics selectivity.

### Result: K-selectivity NOT achieved

K=10 breaks WORSE than lower K under perturbation, with both energy sources. Higher resolution (26 voxels/λ) didn't help.

| Energy source | Perfect placement | With perturbation |
| --- | --- | --- |
| Total phasor | all K stable | K=10 breaks worst |
| Signed envelope | all K stable | K=10 breaks worst |

### Root cause analysis

1. **Signed envelope: shallow equilibria** — wells at D=λ are where other WC's envelope is zero. No deep restoring force. Perturbation pushes WCs out irreversibly.

1. **Total phasor: oscillating cross-term** — genuine interference wells but also barriers. Under perturbation, 45 competing pair forces amplify displacement instead of restoring.

1. **More pairs = more instability** — K=10 has 45 pairs (competing forces in many directions) vs K=2 has 1 pair. Simpler structures are more robust.

1. **Point-sampling limitation** — force from gradient at a single grid point. Collective 10-WC standing wave may need volume integration, not point sampling.

1. **Wells too broad or too deep** — Combined W-L envelope has no zeros in [0, λ), making the entire range a single basin. No fine structure to discriminate K=10 geometry from others.

---

## OPEN QUESTIONS

1. **Variable λ(r)**: does non-uniform wavelength break the symmetry that makes all K stable? Different K geometries have different multi-body wavelength profiles
1. **Non-linear terms**: Ψ³ soliton stabilizer (Smoliński) — K-dependent well structure?
1. **Volume-integrated force**: integrate ∇E over sphere around WC instead of point-sampling
1. **M4 vector field**: L→T spin converts longitudinal to transverse. Spin geometry IS K-dependent — might provide missing selectivity
1. **Self-organization**: scatter K WCs randomly, let only stable K survive (long sim times)
1. **Broadband waves**: multiple frequencies create K-dependent resonances?
1. **Different wave equations**: the energy well structure depends critically on the spatial function. More equations need testing.

---

## M3 SIMULATION PIPELINE

Scripts at `openwave/xperiments/m3_wolff_lafreniere/`. Launched from `_launcher.py`.

### Per-frame computation

1. **`wave_engine.propagate_wave()`** — displacement (Combined W-L), phasor RMS, signed envelope, EMA-RMS, energy
1. **`force_motion.compute_force_vector()`** — `F = -∇E` weighted multi-shell gradient (R=3, 1/d²). Energy switchable: signed envelope or phasor
1. **`force_motion.integrate_motion_leapfrog()`** — symplectic, velocity clamp at c, configurable damping
1. **`force_motion.detect_annihilation()`** — opposite-phase pairs within 1.2λ

### Xperiments

| Xperiment | Setup | Status |
| --- | --- | --- |
| annihilation1 | 2 WC, opposite phase, head-on | Working (INIT_VELOCITY + threshold) |
| annihilation2 | 2 WC, opposite phase, diagonal | Working (INIT_VELOCITY + threshold) |
| tetra_electron | K=10 electron, 1-3-6 | Stable at perfect placement only |
| tetra_positron | K=10 positron, rotated | Stable at perfect placement only |
| formation | K=2..10 with perturbation | K-selectivity NOT achieved |
| tetras | Multi-particle interaction | Pending |

### Features added

- **INIT_VELOCITY**: per-xperiment initial velocity (annihilation momentum injection)
- **VELOCITY_DAMPING**: per-xperiment (1.0 = no damping for annihilation, 0.995 for lock-in)
- **Velocity clamp to c** (0.3 am/rs): both Euler and leapfrog integrators
- **Annihilation threshold**: 1.2λ (standing wave core overlap)
- **Tetrahedron rotation/translation**: grid-independence validation
- **K=2..10 geometry generator**: with configurable perturbation
- **Multiple signed envelopes**: per wave equation (Wolff, Combined W-L, LaFreniere, Weighted PSW)

### The sinc barrier problem (Phase 1 conclusion, confirmed)

Coherent monochromatic wave interference ALWAYS produces `cos(k·Δr)` oscillation in phasor cross-terms. No single-frequency wave equation avoids this. Smooth far-field force requires imposed charge sign (not emergent). Carry-over approaches for Layer 4: 3D flux, variable λ(r), non-linear Ψ³, K=10 scale averaging.
