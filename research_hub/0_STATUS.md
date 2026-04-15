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

## 🧪 LAGRANGIAN EVALUATION (Sub-Project) — see [2a_lagrangian_eval.md](2a_lagrangian_eval.md)

A new research thread sparked by an email exchange in the "Models of Particles" group with two prominent researchers. The goal: evaluate whether a Lagrangian formulation could derive the correct wave equation from first principles instead of testing equations empirically.

### The Real Payoff

If Duda's LdG Lagrangian (or Close's elastic solid Lagrangian) is correct, it would tell us:

- **The exact wave equation** (no more testing 5 candidates empirically)
- **The exact energy functional** (possibly different from `ρV(fA)²`)
- **Why charge is quantized** (topological invariant of the Lagrangian's symmetry group)
- **What conservation laws hold** (and which don't, under which conditions)

In short: it would replace our empirical wave-equation search with a first-principles derivation, and answer Duda's "show me your Lagrangian" challenge in a way that connects OpenWave to established physics formalisms (QFT, GR, Standard Model).

### Key People

- **Jarek Duda** (Jagiellonian University) — inventor of ANS data compression (zstd, JPEG XL, LZFSE). Works on liquid crystal particle analogs and topological field theory. Argues OpenWave needs (1) a Lagrangian, (2) topological charge quantization to prevent the electron from "exploding"
- **Robert Close** (Clark College, retired) — author of "Plane Wave Solutions to a Proposed 'Equation of Everything'" (Foundations of Physics, 2025). Derives the Dirac equation from classical wave mechanics in an ideal elastic solid. Has a constructed Lagrangian with classical interpretation
- **Yves Couder** (deceased) and team — bouncing droplet experiments showing orbit quantization, diffraction, and tunneling in classical wave-particle systems. The closest physical analog to what OpenWave simulates
- **John Bush** (MIT) — leads the modern walking droplet research program. His "Pilot-wave hydrodynamics" (Annual Review of Fluid Mechanics, 2015) is the canonical review of the field

### Duda's Three Challenges

1. **"Show me your Lagrangian"** — we test wave equations empirically; need to derive them from a variational principle
1. **"Without charge quantization your electron explodes"** — our `cos(source_offset)` is imposed; topology forces integer charges naturally (Gauss-Bonnet)
1. **"Particles as topological defects, not standing waves?"** — defects (hedgehogs) protected by topology vs. our standing waves protected by interference

### Two Lagrangian Candidates

| Source | Lagrangian type | Field | Quantization mechanism |
| --- | --- | --- | --- |
| **Duda** | Landau-de Gennes + Skyrme kinetic | Director field M(x) | Topological (winding numbers) |
| **Close** | Elastic solid spin density | Vector spin density | Nonlinear (fixed amplitudes) |

Both may be needed: Duda's topology handles charge quantization and far-field Coulomb; Close's nonlinearity handles amplitude quantization and the Dirac formalism.

### Key Insights from Duda's Replies

1. **Both topological AND standing wave quantization needed** — topology for charge, standing waves for orbit quantization (Couder droplets)
1. **Time crystals** explain WHY particles oscillate (mass-driven, not assumed)
1. **Coulomb needs regularization** → running coupling effect
1. **Klein-Gordon is effective**, not fundamental — a deeper model should average to it
1. **LdG + Skyrme kinetic term** — open question on Higgs-like potential
1. **Three lepton families** from biaxial hedgehog confirmed
1. **Sine-Gordon equation** as the entry point — soliton kinks with pair creation, special relativity
1. **Zitterbewegung** — electron trembling motion at 2mc²/ℏ, must emerge naturally
1. **Faber's 4D approach** automatically produces Zitterbewegung from curvature coupling

### Quick Tests Planned (scripts_phase2_lagrangian/)

8 numpy research scripts, no M4 refactor needed. Decision strategy: run all tests first, then select the winning equation/approach to implement in M4.

| # | Test | What it validates |
| --- | --- | --- |
| 1 | **Sine-Gordon 1D solitons** | Kink stability, pair creation/annihilation, Lorentz contraction |
| 2 | **Hedgehog energy vs distance** | Clean 1/r Coulomb from topology (no sinc) |
| 3 | **Topological charge quantization** | Winding number returns integers under perturbation |
| 4 | **Klein-Gordon from twist dynamics** | Massive dispersion ω² = c²k² + m² |
| 5 | **Lagrangian derivation** | Can our Combined W-L come from a valid Lagrangian? |
| 6 | **Three lepton families** | Biaxial hedgehog → 3 mass scales (e/μ/τ) |
| 7 | **Close's nonlinear vector wave eq** | Particle emergence from spherical harmonic seed |
| 8 | **Smolinski's non-linear Ψ³** | Direct K-selectivity test — does `-k·Ψ³` make K=10 stable? |

**Recommended order**: Test 1 (Sine-Gordon, build intuition) → Tests 2+3 (Coulomb from topology + winding number) → Test 5 (math) → Test 8 (Smolinski Ψ³ K-selectivity on familiar scalar setup) → Test 4 (Klein-Gordon dynamics) → Test 6 (lepton families, most complex) → Test 7 (Close's invitation).

#### Practical Implication

If we implement Test 2 (hedgehog energy vs distance) and seed the field with WL's in+out structure (instead of just a static hedgehog), we might get the best of both worlds:

- Static hedgehog topology → integer winding, charge quantization, Coulomb 1/r
- WL standing wave dynamics → lock-in wells, near-field interactions, energy

This could be the unification: not "topology vs waves" but "topology IS the geometry, waves ARE the dynamics."

### Why This Matters for K-Selectivity

- Avenue **#2 (non-linear ψ³)** gets a major boost — bridges directly to Duda's framework via quartic potential Lagrangian
- Avenue **#7 (different wave equations)** could be resolved entirely — Lagrangian derives the correct equation from first principles instead of empirical testing
- Tests 1+2 could potentially address the **far-field Coulomb blocker** and **charge quantization** simultaneously, while avenues 1-6 remain open for the **near-field K-selectivity** problem

### Architectural Implications

- Lagrangian model needs a **vector field → M4** (not M3 scalar) — director field requires 3+ components per voxel
- Needs a **background vacuum field → M2-like philosophy** (not M3's "WCs emit into void")
- Future **M5 method** would combine: M2 (background field) + M3 (near-field results) + M4 (vector infrastructure) + Duda's topology + time crystal dynamics + Close's Lagrangian

---

## 📊 INSTRUMENTATION (Future Work)

At the end we'll need instrumentation to collect numerical evidence:

- Measure barrier heights at λ/2 → compare with positronium
- Measure well depth, oscillation period, escape energy
- Wrap-up: collect proof, rename xperiments, review EWT phases

**Open question**: maybe we don't have much to compare against — only electron radius, energy, and mass are well-measured at this scale.
