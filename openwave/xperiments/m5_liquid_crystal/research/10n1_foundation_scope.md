# N1 + N2 foundation build: scope (the precision-safe closed-loop neutrino sim)

## What this is

The detailed scope for the **foundation phase (N1 + N2)** of the neutrino N-program ([`10a_neutrino_oscillations.md`](10a_neutrino_oscillations.md) § "Sub-tasks: phase-wired"), to be EXECUTED on Rodrigo's "go N1+N2". N1 (the numerical method) and N2 (the closed-vortex-loop sim) are planned and run **together** as one build , the buildable foundation, no Duda dependency. They produce a precision-safe closed-loop simulation with the mixing observables, ready for the N3 parameter search.

> **WORKFLOW NOTE , this build HOLDS GitHub (Rodrigo's call).** On "go N1+N2" the full task workflow runs (PLAN -> EXECUTE -> FINISH -> REVIEW) with the standard LOCAL artifacts (scripts in `sandbox_v10/`, a findings doc, checkpoints, a terminal TASK REVIEW). BUT **all GitHub #236 posting is HELD until the whole N-program (N1-N5) is finished**: #236 stays `In progress` throughout, NO status change and NO review comment is posted at the end of N1+N2, everything stays local until the N-program completes and the distilled writeup is published once. The resume-ping + mobile completion ping still apply (they are not GitHub).

## Goal of the foundation (N1 + N2)

A **precision-safe closed-vortex-loop simulation** in `sandbox_v10` that: (a) carries the physical regime `g ~ 1e10`, `delta ~ 1e-10` without precision loss (N1); (b) seeds the three neutrino flavour states as closed-loop configurations and evolves the SO(3) rotation among them (N2); (c) exposes the observables the N3 search will read to compare against the [#199](https://github.com/openwave-labs/openwave/issues/199) PMNS scorecard. The foundation does NOT derive the angles , that is N3 (search) + N4 (`theta_13`). It builds the machine they need.

## N1 , the numerical method (precision-safe + Duda's convention)

### The problem

`M = O * diag(g, 1, delta, 0) * O^T` with `g ~ 1e10` and `delta ~ 1e-10` spans a ~`1e20` dynamic range. Energy expressions (`Tr(M^2)`, the commutators `[M_mu, M_nu]`, the LdG potential) mix `g^2 ~ 1e20` and `delta^2 ~ 1e-20` in the same sum -> catastrophic cancellation in f32/f64.

### The approach (two complementary techniques)

1. **Non-dimensionalization.** `g` is a FIXED background scale (the rest-mass / boost axis). Measure energies in units of `g` (or `g^2`) so the DYNAMICAL quantities (the spatial structure, the SO(3) rotation) are O(1) and `g` enters only as a fixed multiplicative constant , removing the `1e20` range from the evolved quantities.
2. **Perturbative-`delta` expansion.** `delta ~ 1e-10` is tiny, so expand the field + energy in powers of `delta` around `delta = 0`: `M(delta) = M_0 + delta * M_1 + O(delta^2)`, computed order by order. The `delta = 0` field `M_0` is the SO(3)-symmetric (TBM) structure; the O(`delta`) term `M_1` is the leading SO(3)-BREAKING. **This is not just a numerical trick , it is physically aligned with the `theta_13` hypothesis** (`theta_13`, the SO(3)-breaking, is hypothesised to be the leading O(`delta`) effect; [`10a`](10a_neutrino_oscillations.md) § "the connecting hypothesis"). So N1's method directly sets up N4's test.

### Duda's convention

Adopt `D = diag(g, 1, delta, 0)` with `g` at INDEX 0 and the Minkowski metric `eta = diag(-1, 1, 1, 1)` (g and the minus both at index 0), per Duda's 2026-06-21 correction. (The existing engine uses the equivalent index-3 ordering; `sandbox_v10` uses index-0 to match Duda and avoid the recurring confusion.)

### N1 deliverable

A documented numerical formulation (in code) that: (1) represents `M` at `g ~ 1e10`, `delta ~ 1e-10` WITHOUT precision loss , verified by a **cancellation test** (a quantity that should be ~0 by cancellation comes out ~0, not ~`1e10` noise); (2) computes the order-by-order (`delta^0`, `delta^1`) energies/forces; (3) uses Duda's index-0 convention; (4) is validated against a known limit (`delta = 0` reproduces the SO(3)-symmetric baseline; the O(`delta`) correction is controlled).

## N2 , the closed-vortex-loop sim (sandbox_v10)

### The object

A neutrino = a CLOSED LOOP of topological vortex (Duda 2026-06-21; could be physically large, ~6.2 pm, Nature s41586-024-08479-6). In the M5 LdG substrate a vortex is a disclination line (the director winds by `2*pi` around it); a closed loop = that line forms a closed curve in 3D.

### What N2 builds

1. **Closed-loop seeder.** Seed `M(x) = O(x) * D * O(x)^T` for a closed disclination loop: the director frame `O(x)` winds around the loop. New geometry , the existing seeders do point-like hedgehogs + straight disclinations; a closed loop is the new primitive. Reuse the LdG potential + the evolution kernels where possible; build the loop topology fresh.
2. **The three flavour states.** The three neutrino flavours as three loop/frame configurations related by the SO(3) rotation (per [#199](https://github.com/openwave-labs/openwave/issues/199): oscillation = a single SO(3) rotation in flavour space). N2 seeds the three configs.
3. **The SO(3) rotation dynamics.** Evolve / rotate the frame among the three flavour configurations (the oscillation), in the N1 precision-safe formulation.
4. **The mixing observables.** The readout N3/N4 will turn into the mixing angles: the overlap structure between flavour states (the projection of one flavour onto another under the SO(3) rotation + the O(`delta`) breaking). N2 EXPOSES these; N3/N4 MEASURE the angles.

### N2 deliverable

A documented `sandbox_v10` simulation that seeds the three closed-loop flavour states, evolves the SO(3) rotation among them in the N1 formulation, and exposes the mixing observables , a working, precision-safe machine the N3 search drives. It does NOT yet claim any angle value (that is N3/N4).

### Honest engineering unknowns (N2)

- Seeding a STABLE closed disclination loop in the 4x4 substrate is the main new engineering (loop tension, core regularization, stability under evolution). May need iteration.
- The exact mapping (three flavours <-> three loop configs; mixing angles <-> the SO(3) + O(`delta`) geometry) is partly what N3/N4 explore. N2 builds the most physically-motivated machinery and exposes the observables, without hard-coding the answer.

## Joint Definition of Done (N1 + N2)

The foundation is DONE when:

1. **N1**: `M` at `g ~ 1e10`, `delta ~ 1e-10` is precision-safe (the cancellation test passes), order-by-order in `delta`, in Duda's index-0 convention.
2. **N2**: the three closed-loop flavour states seed STABLY, the SO(3) rotation evolves among them, and the mixing observables are exposed and finite.
3. A **findings doc** records the formulation, the loop construction, the observables, every parameter, and the validation , reproducible by a third party.
4. The machine is READY to be DRIVEN by N3 (the search varies `g`, `delta`, the potential; reads the observables).

NOT in scope for N1+N2 (resist scope drift): finding the parameters that reproduce TBM = N3; reproducing the angles = the TBM gate; `theta_13` = N4.

## Artifacts

| Artifact | Where |
| --- | --- |
| N1 numerical method + N2 loop-sim scripts | `sandbox_v10/` (naming finalized at go-time, proposed `n1_*.py` / `n2_*.py`; adjustable) |
| Findings doc (the foundation record) | `sandbox_v10/` (e.g. `n_foundation_findings.md`) , LOCAL |
| Checkpoints (each sub-result on arrival) | `sandbox_v10/checkpoints/` |
| This scope | `research/10n1_foundation_scope.md` |

## Risks

- **N1**: a wrong factorization silently reintroduces the precision loss , the cancellation test is the guard.
- **N2**: a stable closed disclination loop in the 4x4 substrate is non-trivial new geometry (the main engineering risk); budget iteration.
- **Scope discipline**: N1+N2 build the MACHINE only; do not drift into the N3 search or the angle measurement (separate, gated phases).

## Cross-refs

[`10a_neutrino_oscillations.md`](10a_neutrino_oscillations.md) (master plan + the N-program wiring + the TBM gate) · issue [#236](https://github.com/openwave-labs/openwave/issues/236) (HELD until the N-program finishes) · [#199](https://github.com/openwave-labs/openwave/issues/199) (the SO(3)/TBM scorecard = the target) · the engine for reuse: [`../medium.py`](../medium.py), [`../engine2_pde.py`](../engine2_pde.py) (seeders [`../engine1_seeds.py`](../engine1_seeds.py), observables [`../engine3_observables.py`](../engine3_observables.py)).

---
_Local scope doc. N1+N2 are planned together, run together on Rodrigo's "go N1+N2", and held local (no GitHub #236 posting) until the whole N-program finishes. Then the distilled writeup is published once._
