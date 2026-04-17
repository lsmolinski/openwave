# THE WAVE EQUATION

The central technical question of OpenWave: what is the correct wave equation that produces particle emergence, force unification, and quantization from first principles?

---

## M2 (PDE Solver) and the Lagrangian Connection

### The Question

> M2 PDE wave equation solver is pretty close to the Lagrange and Euler-Lagrange equations, right? But as we documented, our M2 method / PDE was not enough to show WC disturbances — a stabilized wave field governed by a PDE solver creates an isotropic standing-wave-like background field that does not get disturbed by WCs acting as Dirichlet boundary conditions (ψ=0). How does this new approach change that?

This question gets at the heart of why M2 failed and how the Lagrangian approach fixes it.

---

## M2 IS Already Lagrangian-Based

M2's Laplacian PDE solver evolves the standard wave equation:

```text
∂²ψ/∂t² = c²∇²ψ
```

This IS the Euler-Lagrange equation of the simplest possible field Lagrangian:

```text
L = ½(∂ψ/∂t)² - ½c²(∇ψ)²
                ↓
        Euler-Lagrange
                ↓
       ∂²ψ/∂t² = c²∇²ψ  ✓
```

So M2 is already doing Lagrangian field theory — the **free wave Lagrangian** with no potential V(ψ). Every Taichi tick is solving an Euler-Lagrange equation.

---

## Why M2's Lagrangian Couldn't Show WC Disturbance

Three structural limitations prevent the WC from disturbing the field:

### 1. The Lagrangian Has No Potential (V = 0)

```text
M2's Lagrangian:    L = ½(∂ψ/∂t)² - ½c²(∇ψ)²              ← LINEAR
                                                              (free wave)

Nonlinear option:   L = ½(∂ψ/∂t)² - ½c²(∇ψ)² - V(ψ)       ← NONLINEAR
                                                              (Sine-Gordon, phi-4, LdG)
```

**Without V(ψ), there are no solitons.** Linear wave equations have only two kinds of solutions: traveling waves (which disperse) and standing waves (which fill the entire domain). There's no localized stable structure. Any disturbance you create just spreads out.

The nonlinear potential is what creates **topologically protected solitons** — localized field configurations that cannot dissipate because the math forbids it.

### 2. Scalar Field Has No Topology

A scalar ψ has only magnitude. There's no "direction" to wind around a center. With a scalar field, you cannot define a winding number — there's no hedgehog topology possible.

A vector field n(x) (or matrix field M(x)) DOES have direction at each point. You can ask "how does the direction wind as I go around this center?" → integer winding number → topological protection.

**M2's scalar choice eliminated topology from the start.**

### 3. Dirichlet BC Makes the WC Passive, Not Active

This is the most subtle point. When M2 imposes ψ=0 at a WC location:

```text
M2 approach:        WC = constraint (ψ=0 forced at one grid point)
                    The field adjusts AROUND the constraint passively
                    Wave passes through, doesn't reorganize
                    The WC is "outside" the field, imposed externally

New approach:       WC = part of the field configuration itself
                    The field n(x) winds around the WC by topology
                    The WC is built INTO the field, not imposed on it
                    The surrounding field MUST organize around it
```

**Analogy**: in M2, the WC is like putting a thumbtack on a stretched membrane — the membrane just sags around the tack but is otherwise unaffected. In Duda's approach, the WC is more like a knot in a rope — the knot IS part of the rope; you can't have the rope without the knot pattern.

---

## How the New Approach Changes Each Limitation

| Problem in M2 | Lagrangian fix |
| --- | --- |
| Linear Lagrangian, no potential | Add V(ψ) — Sine-Gordon `1-cos(ψ)`, phi-4 `(ψ²-1)²`, or full LdG → creates solitons |
| Scalar field, no topology | Use vector field n(x) (M4) or matrix field M(x) → enables winding numbers |
| WC as Dirichlet BC (passive) | WC IS a topological defect of the vector field (active, built into the field configuration) |

---

## The Concrete Difference for Test 1 (Hedgehog Energy vs Distance)

**In M2**: we'd initialize the field as zero, drop a WC as a Dirichlet point, evolve, and watch the wave do nothing interesting around the WC.

**In the Lagrangian/topological approach**: we initialize the entire field as

```text
n(x) = (x - x_center) / |x - x_center|
```

— a hedgehog. Every voxel in the field already has a defined direction pointing radially outward from the center. Then we relax the field to minimize the Frank elastic energy:

```text
H = K · Σ |∇n|²
```

The defect is unavoidable — you literally cannot smoothly define n at the center because the radial directions all converge there. The surrounding field MUST organize itself to be consistent with the central winding.

The energy required to maintain this topological structure spreads outward as the 1/r Coulomb tail — automatically, without any imposed boundary condition.

---

## M2's Insight Was Correct

M2's "let waves stabilize themselves" philosophy WAS the right idea — you can't impose disturbances externally; the field has to host them naturally. M2 just used the wrong field type (scalar) and the wrong Lagrangian (linear).

The Lagrangian/topological approach keeps M2's philosophy but upgrades:

- **Linear → nonlinear** (V(ψ) ≠ 0, creates solitons)
- **Scalar → vector** (enables topology)
- **BC-imposed WCs → topological-defect WCs** (built into the field)

This is why the new approach is the natural evolution of M2 — same philosophy (field self-organization), upgraded ingredients (nonlinear potential + vector field + topological defects).

---

## M3: The Quest for the Wave Equation vs. Lagrangian Mechanics

In M3 (Wolff-LaFreniere) we've been searching for the perfect wave equation **empirically** — testing 5 candidates, tweaking parameters, looking for the one that satisfies all our requirements (lock-in, annihilation, K-selectivity, far-field Coulomb). This is the **pedagogical / engineering view** of physics: start from a physical system, derive the wave equation, study its solutions.

The Lagrangian approach **inverts this entire process**.

### The Inverted Hierarchy

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

The Wikipedia "Wave equation" page is fine reference material for the linear free wave case, but it won't help us find the right potential V(ψ) or vector structure. For that, we need the field theory pages.

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
1. **[Klein-Gordon equation](https://en.wikipedia.org/wiki/Klein%E2%80%93Gordon_equation)** — first nontrivial example (free wave + mass)
1. **[Sine-Gordon equation](https://en.wikipedia.org/wiki/Sine-Gordon_equation)** — first nonlinear example (solitons!)
1. **[Landau-de Gennes theory](https://en.wikipedia.org/wiki/Landau%E2%80%93de_Gennes_theory)** — Duda's specific framework

This sequence takes you from the math foundation up to Duda's framework in 6 steps — and at the end, you have all the tools needed to write down the OpenWave Lagrangian and derive its wave equation.

### Note: Why Wikipedia Hides This

Wikipedia organizes physics by **equation/topic**, not by **formalism**. The "Wave equation" page focuses on the PDE itself and its mathematical properties (d'Alembert's formula, characteristic curves, solution methods, boundary conditions). The Lagrangian/Hamiltonian treatments live on the dedicated pages listed above.

There's also a historical reason: the wave equation was derived in **1747 by d'Alembert** for a vibrating string using Newton's laws directly — *before* Lagrangian field theory existed in any usable form. Lagrange's analytical mechanics (1788) and Hamilton's principle (1834) came later, and field-theoretic Lagrangians weren't fully developed until the late 19th and early 20th century (Maxwell, then formalized by Hilbert, Noether, and the QFT pioneers).

So the Wikipedia "Wave equation" page reflects the **historical/pedagogical** path. Most readers (engineers, undergraduates, mathematicians) want the mechanical derivation that connects to physical systems they care about. The Lagrangian view is more abstract and assumes you already understand the equation — it's a graduate-level perspective.

For OpenWave's research goals, we need the graduate-level perspective.

---

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

The direct-line implication for M5 design ([3c_path_to_m5.md](3c_path_to_m5.md)):

1. **No new vector infrastructure needed.** M4's (ψ_x, ψ_y, ψ_z) displacement and 6-phasor trackers already supply everything the director needs
2. **The winding-number tracker is a new kernel, not a new data structure.** It reads from existing ellipse data
3. **Experiment 7 (Close) maps cleanly onto M4's data layout** — same vector field, different time-evolution law
4. **Chirality handling is free.** M4's CW/CCW distinction (from the phase structure of the ellipse) already encodes spin — which is what Duda calls "topological charge + orientation"

The take-away for intuition: **OpenWave has always been about granules moving in ellipses. The Lagrangian / topological framework just asks: "and what integer winds around each WC?" — which is a measurement on top of what we already compute, not a replacement for it.**
