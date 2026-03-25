# PHASE 1: RESULTS, CONCLUSIONS & CARRY-OVER

## What We Set Out To Do

Provide numerical evidence that the electric (Coulomb) force emerges from wave interference — specifically, that two wave centers in a base wave field produce charge-dependent force direction (same repels, opposite attracts) through wave physics alone, without imposing charge as a ±1 label.

---

## What Emerged From Waves (confirmed)

### Strong Force / Lock-In

The sinc oscillation `cos(k·Δr + Δφ)` in the interaction energy between two monochromatic spherical WCs creates energy wells at half-wavelength intervals. Same-phase WCs lock into these wells — this IS the strong force binding mechanism for quarks, nuclear binding, and particle formation.

- **Evidence**: all 1D, 2D, and 3D tests consistently show oscillatory attraction/repulsion at λ/2 intervals
- **Lock-in wells**: same phase WCs have wells at Δr = λ/2, 3λ/2, 5λ/2 ... (odd multiples)
- **Equilibrium points**: quarter-λ offsets where force crosses zero (LaFreniere's "capture")
- **Jeff Yee confirmation**: "I think it's a good thing that we're getting the physics we want for near-field, because this is the key to explaining particle formation. Simulations based on real math/physics that show this is actually the real breakthrough."

#### ❌ Electron Tetrahedral Stability (WIP)

M3 repulsion4 test: K=10 tetrahedral arrangement (1-3-6 geometry) is unstable. 15/45 WC pairs sit at non-node distances (√3×λ/2, √2×λ/2). The uniform λ/2 sinc lattice is incompatible with tetrahedral geometry. Needs variable λ(r) (non-uniform node spacing) and/or M4 vector forces.

### Annihilation

Opposite-phase WCs have their deepest energy well at Δr = 0 (zero separation = complete wave cancellation). Barriers at λ/2 intervals can temporarily trap opposite-phase pairs (positronium-like bound states, lifetime ~125 ps to ~142 ns). High kinetic energy overcomes barriers → direct annihilation.

- **Evidence**: interaction energy `E_int ∝ -cos(k·Δr)` for opposite phase → minimum at Δr = 0
- **Positronium analog**: barriers at λ/2, 3λ/2 create metastable states before annihilation
- **Consistent with experiment**: electron-positron annihilation is well-observed

### Charge Discrimination (from 2D flux)

The 2D radiation pressure (flux-based force `S = -c²·ψ·∇ψ`) produces **100% charge discrimination**: same and opposite phase ALWAYS get opposite force directions at every separation tested (22/22). The sinc oscillation determines the absolute direction, but the relative direction (same vs opposite) is always opposite.

- **Evidence**: `step1d_flux_force_2d.py` — 22/22 separations show same ≠ opposite
- **Coulomb-correct zones**: at half-integer λ separations (2.5λ, 3.5λ, 4.5λ...) the directions match Coulomb (same=REP, opp=ATT)
- **Force decays**: ~1/r trend under the oscillation (consistent with Coulomb scaling)

### Neutrino Neutrality (K=1)

Per EWT (Jeff Yee): spin only occurs at K≥10 (electron scale). A single WC (K=1, neutrino) has no spin, no charge, no L→T conversion. The sinc oscillation at K=1 IS the correct physics — lock-in for particle formation, no far-field Coulomb. Neutrinos are neutral.

- **No observational data** exists for neutrino forces — we can't validate K=1 against experiment
- **The sinc "blocker" was never a bug** — it's the correct K=1 physics

### Spin → Magnetic Force

L→T spin conversion at a WC creates a transverse energy component perpendicular to the WC axis. This produces force in the transverse direction (magnetic), NOT in the radial direction (electric). CW and CCW spin produce identical energy for a single WC; spin sign matters only for WC-WC interactions.

- **Evidence**: Phase 1c Step 3 — all 4 spin configs show MIXED radial force. The T component is perpendicular to the connecting axis.
- **Physical η = α ≈ 1/137** (fine structure constant): the fraction of longitudinal energy converted to transverse by spin

---

## What Did NOT Emerge (the remaining challenge)

### Far-Field Coulomb (consistent direction at all separations)

The sinc oscillation `cos(k·Δr)` flips the absolute force direction every λ/2 of separation change. We could not produce consistent Coulomb direction (same always REP, opposite always ATT) at ALL separations.

**Systematically eliminated mechanisms:**

| Approach | Phase | Result | Why it failed |
| --- | --- | --- | --- |
| All 5 wave equations | 1 | ❌ All produce sinc | Oscillation intrinsic to coherent interference |
| Gaussian smoothing | 1 | ❌ Destroys charge signal | Removes oscillation AND charge info |
| Signed disturbance | 1a | ❌ Charge imposed, not emergent | ±1 label on smooth potential |
| 10 WC disturbance models | 1b | ❌ Only L→T distinguishes charges | Passive models fail in isotropic fields |
| Spin L→T conversion | 1c | ❌ Creates magnetic, not electric | T perpendicular to WC axis |
| Variable λ(r) energy | 1d | ❌ Charge-blind | λ depends on K (structure), not phase (charge) |
| 2D/3D integration (gradient) | 1d | ❌ Base wave dominates | WC-base interaction >> WC-WC interaction |
| Statistical averaging (K²λ) | 1d | ❌ Sinc symmetric → zero | Equal positive/negative areas per cycle |
| Broadband shells | 1d | ❌ Multiple cosines oscillate | Beating pattern, not convergence |
| In-wave / out-wave decomp | 1d | ❌ Sinc in ALL components | Out-wave has sinc spatial structure |

### Correct Coulomb Sign

The isolated interaction energy (E_int = cross term) gives same→ATT, opposite→REP. This is 180° opposite to Coulomb. This may be the strong force capture mechanism (same-phase attraction via radiation pressure), not Coulomb.

---

## Fundamental Insight: The Sinc Is the Physics

The sinc oscillation `cos(k·Δr + Δφ)` is NOT a numerical artifact or model limitation. It is the **fundamental mathematical consequence** of two coherent monochromatic spherical waves interfering in space. Any model that uses:

1. Single frequency ω (monochromatic)
2. Spherical wave propagation (1/r decay)
3. Coherent superposition (amplitudes add)

...will ALWAYS produce `cos(k·Δr)` in the interaction energy. This is proven math — no variable λ, no spin, no decomposition, no dimensional trick changes it.

**The sinc simultaneously creates:**

- Strong force lock-in (same phase wells at λ/2)
- Annihilation pathway (opposite phase well at r=0)
- Charge discrimination (same ≠ opposite in both gradient and flux)
- Near-field particle formation (equilibrium at capture points)

**The sinc does NOT create:**

- Consistent far-field Coulomb direction at all separations
- The correct Coulomb sign (same=REP, opp=ATT) in the energy gradient approach

---

## Key Theoretical Insights Discovered

### LaFreniere Phase Shift

The electron core (1λ diameter) creates a λ/2 phase shift from medium compression (7x smaller volume → shorter λ inside → more phase cycles → phase advance). This phase shift — not spin — is what creates charge sign. The mechanism is variable λ(r) inside the core, where c stays constant but λ compresses.

### Yee & Hauger vs LaFreniere λ(r) Contradiction

Yee & Hauger shells: λ is LONGEST near core (shell n=1 = 18λ₀ for K=10), shrinking outward. LaFreniere: λ is SHORTEST inside core (compressed). Smoliński (f∝1/r): agrees with LaFreniere. This contradiction is unresolved — Phase 2 should test both profiles.

### Force from Flux vs Gradient

`F = -∇E` (energy gradient) is charge-blind when base wave is present and produces wrong sign for isolated interaction. `S = -c²·ψ·∇ψ` (energy flux / radiation pressure) produces 100% charge discrimination in 2D. These are different physical quantities:

- Gradient: WHERE energy is stored (scalar, no direction preference)
- Flux: WHERE energy is FLOWING (vector, knows wave propagation direction)

LaFreniere explicitly uses radiation pressure (flux), not energy gradient, for the Coulomb mechanism.

### Wave Center Hierarchy (K)

- K=1 (neutrino): no spin, neutral, sinc = lock-in physics only
- K=10 (electron): spin (L→T at rate α), charged, tetrahedral 1-3-6 arrangement
- Spin only at K≥10: off-node WC repositioning creates continuous rotation
- All observable force validation must be at K≥10 (no neutrino force data exists)

### Connection: λ → Pressure → Gravity

Variable λ/frequency → variable granule velocity → variable medium pressure/density. The pressure deficit around a massive body IS the gravitational field — emerging from λ modulation, not a separate mechanism. Gravity, time dilation, and medium pressure are three descriptions of the same phenomenon: local λ variation.

---

## Spin, Magnetism, and Thermal Energy

**Spin** may be the unifying concept:

- **Magnetic moment from spin**: L→T conversion at WC creates transverse wave emission. Coherent spin alignment = net transverse field = magnetism
- **Thermal energy from spin**: incoherent spin fluctuations (random T directions) = heat. Higher temperature = more energetic L→T cycling
- **Gravitational shading from spin**: spin converts L→T, reducing longitudinal amplitude. This L reduction is the gravitational "shadow" — the mechanism behind shading/push-out gravity

The fine-structure constant α may be the L→T conversion ratio. If gravity emerges from repeated application of this coupling, the 10⁻⁴² EM-to-gravitational ratio emerges naturally.

### Toroidal Flow Geometry (Smoliński)

The internal "engine" of a particle operates in the Energy Domain with toroidal (doughnut-shaped) wave flows governed by r⁵ scaling:

- **Spin** = toroidal wave flow around the particle core
- **Magnetic moment** = natural dipole from circulating toroidal current
- **Thermal energy** = modulation of toroidal flow rate/precession
- The elliptical displacement trajectories may be cross-sections of toroidal flows

### Vortex Electron Model (Butto, 2021)

The electron as a superfluid irrotational vortex — supporting the toroidal flow picture:

- **Spin-½ = differential rotation**: core rotates at 2× boundary rate (720° core / 360° boundary)
- **Planck constant = vortex angular momentum**: Γ·m_e = h
- **Magnetic moment without g-factor**: μ = qcr from vortex hydrodynamics
- **Helmholtz theorems**: vortex filaments form closed loops → toroidal geometry

---

## Carry-Over: Knowledge & Ideas for Phase 2

### The Most Promising Lead: 2D Flux

The `step1d_flux_force_2d.py` result is the closest we've come to Coulomb: 100% charge discrimination, ~1/r scaling, Coulomb-correct at half-integer λ separations. The remaining issue is the λ/2 absolute direction oscillation. Ideas:

1. **3D spherical flux integration** — the 2D test used a circle. Full 3D sphere may smooth the remaining oscillation through solid-angle averaging (untested)
2. **Flux computed on the M3/M4 Taichi engine** — the time-domain computation with GPU acceleration can handle larger grids and longer time averages
3. **Flux at K=10 particle scale** — test whether the 100λ particle radius naturally averages the flux oscillation (different from averaging the force, which failed)

### Unsolved: The Coulomb Sign Problem

The isolated interaction cross term gives same→ATT (should be REP). The 2D flux gives correct sign at half-integer λ separations but inverted at integer λ. The absolute sign may depend on:

- Whether we're measuring the force at a standing wave node or antinode
- The convention for source_offset (0=positron, π=electron) and whether it maps correctly to "charge"
- Whether the LaFreniere λ/2 core phase shift provides the missing π offset that flips the sign

### Infrastructure Ready for Phase 2

- **3D vector base wave** (step1_base_wave.py): 64³ grid, 200 Fibonacci sources, validated energy/L-T/speckle
- **WC out-wave with L+T** (step2_single_wc.py): sinc envelope, WKB phase, spin direction, energy concentration
- **Variable λ(r) profile** (Yee & Hauger shells): continuous inversion from discrete shells, WKB phase integral
- **2D flux computation** (step1d_flux_force_2d.py): time-domain wave field, radiation pressure measurement
- **Leapfrog integrator** in M3 and M4 force_motion.py: symplectic, with damping
- **K parameter**: tested K=1 (neutrino) through K=10 (electron)

### Deferred Mechanisms (not yet tested)

- Variable ρ(r) from Smoliński density function
- Ψ³ cubic non-linearity (NLS soliton stabilizer)
- Non-coherent superposition (thermal phase jitter)
- Non-spherical propagation (toroidal from spin, diffractive lens)
- Combined vector displacement + variable λ + flux force
- Energy flux (Poynting-like) in full 3D with vector displacement

### External Input Expected

Jeff Yee is bringing Dieter Hauger (co-author of wavelength shells paper) into the conversation. Hauger may have insights on:

- The standing → traveling wave transition mechanism
- Whether the λ(r) profile inside the core compresses or stretches
- How the phase shift produces charge sign

## Other Force Computation Approaches

Beyond `F = -∇E_total`, vector displacement enables alternative force quantities:

- **Divergence** (∇·ψ): compression/rarefaction — scalar but signed
- **Curl** (∇×ψ): rotational displacement — vector, related to magnetic field
- **Energy flux** (ψ × ∂ψ/∂t): directional energy flow (radiation pressure, LaFreniere's mechanism)
- **Per-component force**: `F_x = -∂E_x/∂x`, `F_y = -∂E_y/∂y`, `F_z = -∂E_z/∂z`

Standing waves have zero net flux, traveling waves have nonzero flux — flux naturally separates near-field (standing, lock-in) from far-field (traveling, Coulomb).

All Maxwell's equations have divergence and curl components.
