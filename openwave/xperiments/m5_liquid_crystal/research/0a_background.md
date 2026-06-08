# The Quest for the Wave Equation vs. Lagrangian Mechanics

In M3 (Wolff-LaFreniere) we've been searching for the perfect wave equation **empirically** — testing 5 candidates, tweaking parameters, looking for the one that satisfies all our requirements (lock-in, annihilation, K-selectivity, far-field Coulomb). This is the **pedagogical / engineering view** of physics: start from a physical system, derive the wave equation, study its solutions.

The Lagrangian approach **inverts this entire process**.

## The Inverted Hierarchy

In modern theoretical physics, almost no one writes a wave equation directly. Particle physicists ALWAYS start from a Lagrangian — because it gives them Lorentz invariance, conservation laws, gauge symmetries, and quantization for free. The equation of motion is just a consequence — they derive it via Euler-Lagrange and then forget about it.

```text
Pedagogical / engineering view        Modern theoretical physics view
(what we've been doing in M3)         (what Duda is asking us to do)
─────────────────────────────         ─────────────────────────────────
                                      
  Physical system                       Lagrangian (chosen from symmetries)
       ↓                                          ↓
  Newton's laws                              Action principle
       ↓                                          ↓
  Wave equation  ← THE GOAL              Euler-Lagrange equation
       ↓                                          ↓
  Solutions                              Wave equation  ← just a step
       ↓                                          ↓
  Lagrangian (afterthought)                   Solutions
```

### Why This Matters for OpenWave

For our purposes — deriving the *correct* wave equation for OpenWave — the bottom-up view (Lagrangian first) is the right one. That's exactly what Duda is asking us to do.

The Wikipedia "Wave equation" page is fine reference material for the linear free-wave case, but it won't help us find the right potential V(ψ) or vector structure. For that, we need the field theory pages.

### What This Means in Practice

| Aspect | M3 approach (engineering) | Lagrangian approach (theoretical) |
| --- | --- | --- |
| Starting point | Pick a candidate wave equation | Pick a Lagrangian (from symmetries) |
| Method | Test, tweak, compare | Derive Euler-Lagrange, solve |
| Energy | Postulate `E = ρV(fA)²` | Falls out of Hamiltonian `H = T + V` |
| Conservation laws | Check manually | Guaranteed by Noether's theorem |
| Charge quantization | Imposed via `cos(source_offset)` | Topological invariant of symmetry group |
| Wave equation | Empirically chosen (5 candidates tested) | Derived from first principles |
| Solutions | Phasor superposition (analytical) | Field evolution (PDE, possibly nonlinear) |
| What we control | The equation directly | The Lagrangian (more abstract, more powerful) |

### The Real News: We've Been Doing It Backwards

This is the most important realization from the Lagrangian eval: **we've been searching for the wave equation as the goal, when it should be the consequence.**

The M3 quest — testing Wolff, LaFreniere, Phase-warped Marcotte, Combined W-L, Weighted PSW — is the pedagogical approach. Each candidate is a guess, and we evaluate it empirically. After 5 candidates and dozens of tweaks, the best we have (Combined W-L) explains lock-in but breaks under perturbation, and we have no principled way to know if a sixth candidate would be better.

The Lagrangian approach replaces this empirical search with a derivation:

1. **Start with symmetries** (Lorentz invariance, gauge invariance, topological structure)
2. **Write the most general Lagrangian** consistent with those symmetries
3. **Choose the potential V(ψ)** based on the physics we want (mass, solitons, charge quantization)
4. **Derive the equation of motion** via Euler-Lagrange — automatically
5. **Conservation laws** come for free via Noether's theorem
6. **Solutions** are then either analytical or numerical

Steps 1-3 are conceptual choices (what symmetries does spacetime have?). Steps 4-6 are mechanical math. The wave equation is no longer a thing we search for — it's a thing we calculate.

### Suggested Learning Path (Lagrangian-First View)

If we want to study the Lagrangian view of wave equations in order:

1. **[Action (physics)](https://en.wikipedia.org/wiki/Action_(physics))** — what an action functional is
1. **[Euler-Lagrange equation](https://en.wikipedia.org/wiki/Euler%E2%80%93Lagrange_equation)** — the math machinery
1. **[Lagrangian (field theory)](https://en.wikipedia.org/wiki/Lagrangian_(field_theory))** — apply it to fields ψ(x,t)
1. **[Klein-Gordon equation](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation)** — first nontrivial example (free-wave + mass)
1. **[Sine-Gordon equation](https://en.wikipedia.org/wiki/Sine-Gordon_equation)** — first nonlinear example (solitons!)
1. **[Landau-de Gennes theory](https://en.wikipedia.org/wiki/Landau%E2%80%93de_Gennes_theory)** — Duda's specific framework

This sequence takes you from the math foundation up to Duda's framework in 6 steps — and at the end, you have all the tools needed to write down the OpenWave Lagrangian and derive its wave equation.

### Note: Why Wikipedia Hides This

Wikipedia organizes physics by **equation/topic**, not by **formalism**. The "Wave equation" page focuses on the PDE itself and its mathematical properties (d'Alembert's formula, characteristic curves, solution methods, boundary conditions). The Lagrangian/Hamiltonian treatments live on the dedicated pages listed above.

There's also a historical reason: the wave equation was derived in **1747 by d'Alembert** for a vibrating string using Newton's laws directly — *before* Lagrangian field theory existed in any usable form. Lagrange's analytical mechanics (1788) and Hamilton's principle (1834) came later, and field-theoretic Lagrangians weren't fully developed until the late 19th and early 20th century (Maxwell, then formalized by Hilbert, Noether, and the QFT pioneers).

So the Wikipedia "Wave equation" page reflects the **historical/pedagogical** path. Most readers (engineers, undergraduates, mathematicians) want the mechanical derivation that connects to physical systems they care about. The Lagrangian view is more abstract and assumes you already understand the equation — it's a graduate-level perspective.

For OpenWave's research goals, we need the graduate-level perspective.

## M4's Elliptical Granule Motion and the Director Field

### The Physical Basis of OpenWave

Before comparing M4 to Duda's director field, it's worth stating the physical picture OpenWave is built on:

> OpenWave simulates **Planck-scale granules oscillating in elliptical trajectories around equilibrium positions**. Every voxel in the grid represents a granule whose displacement from rest, over time, traces an ellipse. Waves propagate by these elliptical motions coupling between neighbors.

That is the intuition underneath everything — M1 tracks the granules directly; M2 solves the coarse-grained PDE of granule displacement; M3 computes the analytical standing-wave superposition; M4 resolves the full vector trajectory of each granule's ellipse. The Lagrangian / topological framework rides on top of the same underlying picture — it just adds the observation that *orientation information is already encoded in those ellipses*.

### What M4 Tracks (the 6-phasor)

At each voxel, M4 carries six numbers that fully describe the local granule's elliptical oscillation:

- 3 amplitudes (A_x, A_y, A_z)
- 3 phases (φ_x, φ_y, φ_z)

Together, these six parameters specify the full elliptical trajectory — its shape (eccentricity, size), its orientation in 3D (tilt, rotation plane), and its timing (where the granule sits on the orbit at t = 0). This is a **dynamical** quantity — it describes what the granule *does over time*.

### What the Director Field Is (Duda's Framework)

In the Lagrangian / topological framework (Duda, liquid crystals), each point in space carries a unit vector **n(x)** — a local orientation, like a tiny compass needle. The director field is a **static** quantity: at each point, "which way is up?" The topology is then: how does this orientation wind around a defect?

### The Bridge: Ellipses Carry Orientation Information

An ellipse is not just a dynamic object — it has a natural **geometric orientation** built into its shape:

- The **major axis** defines one direction in space
- The **normal to the orbital plane** defines another direction (perpendicular to where the granule moves)
- The **handedness** (clockwise vs. counter-clockwise traversal) is a ± sign — chirality information

So the M4 ellipse **already contains** director-like structure. M4 just tracks the full dynamical ellipse; extracting the director is a matter of reading off one of its natural axes.

### Side-by-Side Comparison

| M4 ellipse (dynamical) | Director n(x) (static) |
| --- | --- |
| 6-phasor: 3 amplitudes + 3 phases per voxel | Unit vector per voxel (3 components, normalized) |
| Describes the full elliptical trajectory | Describes just the local orientation |
| Contains timing (phase on orbit) | No timing — purely geometric |
| Contains shape (eccentricity, axis lengths) | No shape information |
| Contains handedness (CW / CCW spin of the orbit) | No handedness — unit vectors have no chirality on their own |
| Major-axis direction ⊂ full state | Exactly the major-axis direction |

**If we take just the major-axis direction of the M4 ellipse at each voxel, that *is* a director field.** The winding of that director around a WC computes the topological charge. The chirality (CW/CCW) of the ellipse carries spin information on top.

### Close's "Spin Density Vector" Is Exactly This

In Robert Close's 2025 *Foundations of Physics* paper ("Plane Wave Solutions to a Proposed 'Equation of Everything'"), the primary field is a **spin density vector** — the local axis of rotational motion in the elastic medium. That is, the normal to the plane of elastic motion at each point.

**Close's field is essentially M4's ellipse orientation, extracted as a director.** This is exactly why Experiment 7 (Close's nonlinear vector wave equation) is plausible to implement in M5: the raw data is already there in M4's representation — we'd just be evolving it under a different (nonlinear) PDE.

### The Implication: M4 Is a Superset, Not a Peer, of Duda's Framework

This is important:

- **M4 tracks more than Duda's framework needs.** The 6-phasor contains orientation (director info), chirality (spin info), *and* the full dynamical phase/amplitude information
- **Duda's director field is a projection of M4's state.** Take the major-axis direction at each voxel → you have a director field. Discard the rest
- **Topology is already implicit in M4.** The topological charge around a WC is a function of the ellipse orientations on a surrounding sphere — we already have the data; we just don't compute the winding number yet

### Why This Matters for M5

The direct-line implication for M5 design ([2a_path_to_m5.md](2a_path_to_m5.md)):

1. **No new vector infrastructure needed.** M4's (ψ_x, ψ_y, ψ_z) displacement and 6-phasor trackers already supply everything the director needs
2. **The winding-number tracker is a new kernel, not a new data structure.** It reads from existing ellipse data
3. **Experiment 7 (Close) maps cleanly onto M4's data layout** — same vector field, different time-evolution law
4. **Chirality handling is free.** M4's CW/CCW distinction (from the phase structure of the ellipse) already encodes spin — which is what Duda calls "topological charge + orientation"

The take-away for intuition: **OpenWave has always been about granules moving in ellipses. The Lagrangian / topological framework just asks: "and what integer winds around each WC?" — which is a measurement on top of what we already compute, not a replacement for it.**

---

## Liquid Crystal Primer — the underlying medium of M5

This section can be read first if you're unfamiliar with what a liquid
crystal actually is. The M5 framework treats the vacuum as
liquid-crystal-like; understanding the physical analog makes Duda's
mathematical choices feel less arbitrary.

### The physical material

A **liquid crystal** is a state of matter that sits between a liquid (no
positional order, fully fluid) and a crystal (rigid 3D lattice positional
order). It has:

| Property | Liquid | Liquid crystal | Crystal |
| --- | --- | --- | --- |
| Molecular positions | random | random | lattice |
| Molecular orientations | random | **ordered** (aligned along an axis) | ordered |
| Flow | yes | yes | no |

The defining feature: **molecules align along a preferred direction** (the
"director") even though they're not fixed in space. Think of cigar-shaped
molecules in a soap bubble — they slosh around freely but their long axes
point mostly the same way. LCD displays exploit this — applying voltage
rotates the director, changing how light passes through.

### The mathematical object — the order tensor

Because molecules are aligned but free to slide, you describe the local
ordering with a tensor (not just a vector). For uniaxial molecules with
axis along unit vector `n̂`:

```text
Q_ij = S × (n̂_i n̂_j - δ_ij/3)
```

where `S` is the scalar order parameter (0 = isotropic liquid; 1 =
perfectly aligned). For biaxial molecules (lower symmetry — two distinct
transverse axes), you need a richer tensor:

```text
M = O · D · O^T
```

where:

| Symbol | Meaning |
| --- | --- |
| `O` | 3×3 rotation matrix specifying local orientation of the principal axes |
| `D = diag(λ₁, λ₂, λ₃)` | diagonal of "axis lengths" — the order parameter eigenvalues |
| `M` | symmetric traceless 3×3 matrix (5 degrees of freedom) |

This `M` is the **Landau–de Gennes (LdG) tensor**. The number of DOF (5
per spatial point) is why M5's matrix-field substrate is "richer" than
M6's two 4-vector fields (8 DOF with Lorenz constraints reducing further).

### Why this is useful for particles — topological defects

In an ordered medium with a director field, you can't always smoothly
deform the orientation from one place to another. Sometimes the alignment
breaks down at a point — that's a **topological defect**:

| Defect type | What it looks like | Topological invariant |
| --- | --- | --- |
| **Hedgehog** | director points radially outward (or inward) at a 3D singularity | Brouwer degree ±1 (integer winding around S²) |
| **Disclination** | director rotates by 180° or 360° around a 1D line | line winding number |
| **Skyrmion** | localized topological texture (no singularity, but stable winding) | integer Hopf invariant |
| **Biaxial twist** | the two transverse axes rotate around a structure | richer biaxial-symmetry invariants |

These defects are **topologically stable** — you can't smoothly remove
them without cutting the field. That's the analogy Duda exploits: **a
particle is a topological defect in an underlying liquid-crystal-like
field, and its quantum numbers (integer charge, integer spin) are the
topological invariants of that defect.**

### How M5 maps onto this

| Liquid crystal concept | M5 (Duda's framework) interpretation |
| --- | --- |
| The medium (LC) | The vacuum — a Landau–de Gennes tensor field `M(x, t)` on 3D space, time-evolved |
| Vacuum state | `D = diag(g, 1, δ, 0)` — preferred biaxial ground state shape |
| Director rotation | The matrix `O(x, t)` |
| Hedgehog defect | Charged lepton (electron) — Brouwer degree ±1 = electric charge ±1 |
| Skyrmion / Hopfion | Candidate for excited states / neutrinos (Liu et al. 2026 lab anchor) |
| Director oscillation | Klein-Gordon-like wave dynamics around vacuum → particle mass |
| Standing-wave interference | Orbit quantization (atomic shells) |
| Biaxial-axis hierarchy | Three lepton families (e/μ/τ map to choices of which axis dominates) |
| LC distortion energy (Frank elastic) | Mass + gradient-energy terms in the Lagrangian |
| Higgs-like potential V(M) | Forces M near the vacuum shape; sets symmetry breaking |

### Why it's experimentally interesting

| Real-world LC connection | Why it matters |
| --- | --- |
| LCDs use nematic liquid crystals; defects are visible under polarized light | Visualization analog — actual lab pictures of these defects exist |
| Bush & Couder walking droplets reproduce pilot-wave behavior | Hydrodynamic analog of quantum mechanics in classical fluids |
| Liu et al. *Nature Physics* 2026 — first direct laser creation of isolated hopfions and skyrmions in real LC | Lab anchor: the topological structures M5 hypothesizes for particles have been made in matter |
| Schwinger 1969 + Faber regularization | Mainstream physics ancestry of the same math (the Schwinger family tree) |

### One-line summary

A liquid crystal is matter that has orientation order but no position
order — molecules align without locking into a lattice. Duda's M5 takes
the math of biaxial liquid crystals (the `M = ODO^T` tensor field) and
applies it to the vacuum, treating particles as topological defects in
the local orientation. Lab analogs (LCD displays, hopfion experiments)
make the framework empirically grounded in a way M6 is not.
