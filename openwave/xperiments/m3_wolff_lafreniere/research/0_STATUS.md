# PROJECT STATUS

Quick reference for OpenWave research progress, blockers, and next steps.

---

## 🎯 WHAT PROOFS WE NEED

The end-game targets — phenomena that must emerge from wave physics:

- **PARTICLE STABILITY** (across rotation/translation/phase)
  - Magic numbers stability: K = 2, 8, 10
  - Reference: <https://energywavetheory.com/subatomic-particles/neutrino>
- **PARTICLE ANNIHILATION**
  - Barriers that can be transposed with enough momentum
  - Not super-luminal
- **ATTRACTION / REPULSION** (the fundamental forces)
  - Electric Coulomb (Charge)
  - Magnetic Bohr Magneton (Spin)
  - Gravitational Constant
  - Strong / Gluonic
  - Orbital Force

---

## ✅ WHERE WE ARE — Current State (M3 Scalar, Layers 1-3)

### Achieved

- **Lock-in**: same-phase WCs oscillate in energy wells at λ separation
- **Annihilation**: opposite-phase WCs annihilate (caveats below)
- **K=10 tetrahedron** holds at perfect placement using Combined Wolff-LaFreniere equation
- **PDE vs analytical validated**: Laplacian/Dirichlet BC can't simulate WC disturbance; analytical solutions required

**Annihilation caveats** (not purely force-driven):

- Requires 0.5λ threshold parameter
- Only works on the WOLFF-LAFRENIERE envelope form
- Needs: INIT_VELOCITY + damping + velocity_clamp
- Energy/force gradient sample radius matters
- Annihilation detection vs. test same-position wave cancellation

### NOT Achieved (the open blockers)

- **K-selectivity**: all K=2..10 are equally stable at perfect placement; K=10 actually breaks WORST under perturbation (more pairs = more competing forces). The energy landscape (both phasor and envelope) doesn't discriminate K=10 from simpler geometries
- **Emergent charge sign**: imposed via `cos(source_offset)`, not from wave physics
- **Far-field Coulomb**: sinc barriers block force-driven attraction/repulsion; signed envelope is a modeling choice
- **Perturbation-robust stability**: Combined W-L has shallow equilibria, no fine structure to discriminate K=10 from others

> **Honest conclusion**: M3 scalar with monochromatic waves can demonstrate the MECHANISM (lock-in from standing waves) but cannot yet demonstrate the SELECTIVITY (why K=10 and not K=3). That selectivity likely requires the missing physics — variable λ, spin, or non-linear terms.

### Wave Equations Tested (5 total, all produce sinc)

| # | Equation | Envelope zeros | K=10 at perfect | Under perturbation |
| --- | --- | --- | --- | --- |
| 1 | Wolff-Original | λ/2 | unstable (zero at λ/2) | - |
| 2 | LaFreniere-Marcotte | λ | holds | breaks |
| 3 | Phase-warped Marcotte | ~λ | holds | breaks |
| 4 | **Combined W-L** (current best) | λ | **holds** | **breaks worst** |
| 5 | Weighted PSW | λ/2 | unstable | - |

---

## 🎯 MAIN GOAL: Crack the Wave Equation (Base + Disturbance)

**Hypothesis**: Combined W-L equation creates energy wells too deep / too stable — explaining lock-in but not selectivity.

**Goal**: derive the perfect equation that fits all our requirements:

- Base + disturbance, WC reflection vs. pure longitudinal, phase shift
- Transition standing vs. traveling waves
- Timestep instability handling
- Insights from de Broglie, Schrödinger

### Focus: Force from Granule Motion, Handedness from Elliptical Displacement

Review of granule motion and amplitude:

- Close to center and crossing center
- Amplitude decline
- Lambda change with radius
- The inverse square law
- Sinc function standing wave lock-in

---

## 🔬 OPEN AVENUES FOR K-SELECTIVITY

The 7 main physics directions to explore (from `2L3_particle_emergence.md`):

1. **Variable λ(r)** — Yee & Hauger shells make node spacing non-uniform. Could break the symmetry that makes all K equally stable. Different K geometries have different multi-body wavelength profiles
   - L→T conversion causes the slight difference in symmetry
   - L→T spin only at K=10 layer, so no Coulomb at K=1
1. **Non-linear ψ³** — Smolinski (Lukasz) soliton stabilizer. Changes spatial structure from pure sinc → K-dependent well structure. **Bridge to Lagrangian framework** (comes from valid quartic potential)
1. **Volume-integrated force** — integrate ∇E over a sphere around each WC instead of point-sampling the gradient at a single grid point
1. **M4 vector field** — L→T spin conversion is K-dependent; spin geometry IS K=10 specific
   - Ellipse + spin, toroidal/vortex/spiral structures
   - Lukasz density/pressure variations
1. **Self-organization** — scatter K WCs randomly, let only stable K survive (long sim times)
1. **Broadband waves** — multiple frequencies create K-dependent resonances?
1. **Different wave equations** — the energy well structure depends critically on the spatial function. **Could be resolved by Lagrangian derivation** (first-principles instead of empirical)

---

## 🛠️ NUMERICAL & IMPLEMENTATION TWEAKS

Engineering improvements (not physics changes):

- **Energy computation**: from phasor vs. envelope (workaround), energy flux
- **Force & motion numerical errors**: timestep surfaces these
  - `sample_radius = 1`, force integration (sphere ∇E)
  - `velocity_damping = 1`
- **Perturbations**, `init_velocity` tuning
- **Higher resolutions**
- **Other sources of numerical instability**

---

## 📊 INSTRUMENTATION (Future Work)

At the end we'll need instrumentation to collect numerical evidence:

- Measure barrier heights at λ/2 → compare with positronium
- Measure well depth, oscillation period, escape energy
- Wrap-up: collect proof, rename xperiments, review EWT phases

**Open question**: maybe we don't have much to compare against — only electron radius, energy, and mass are well-measured at this scale.
