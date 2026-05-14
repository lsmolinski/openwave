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
