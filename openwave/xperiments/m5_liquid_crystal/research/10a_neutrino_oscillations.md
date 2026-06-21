# Neutrino oscillation parameters from topological-vortex dynamics: build plan (#236)

## Status

Issue [#236](https://github.com/openwave-labs/openwave/issues/236) (OpenWave board: **In progress**). Spun off on Dr. Duda's 2026-06-20 recommendation to "do less but more rigorously": pivot the first serious, article-targeted deliverable from charged-lepton masses to **neutrino oscillations**, because the geometry is simpler (only topological vortices, no point-like defects) and the payoff is higher (the oscillation parameters are a guessing game today, so a derivation from the deeper theory would be genuinely novel).

> **Round 3 (2026-06-21) sharpened the scope , see § "Round 3" below.** The deliverable narrows to the **4 PMNS mixing parameters** (`theta_12`, `theta_23`, `theta_13`, `delta_CP`); absolute masses (`Delta m^2`) are DEFERRED. Our own [#199](https://github.com/openwave-labs/openwave/issues/199) SO(3)/TBM result already gives 3 of the 4; the live target is **`theta_13` = 8.5deg, the SO(3)-breaking**, plus rigorizing the other 3 from the closed-vortex-loop field theory. The neutrino object = a **closed topological-vortex loop**. The parameter+potential search is a genuine open problem WE drive (Duda confirms even he lacks the exact values).

This doc is the LOCAL plan. It captures Duda's round-2 + round-3 replies verbatim, evaluates them, scopes the neutrino task, and records the (now active-search) blocker. Nothing is published to the issue or committed until the local work is solid (the same discipline as [`9a_lepton_mass_planning.md`](9a_lepton_mass_planning.md)).

| Related work | Link | Relationship |
| --- | --- | --- |
| Charged-lepton masses (PARKED) | [`9a_lepton_mass_planning.md`](9a_lepton_mass_planning.md) · issue [#200](https://github.com/openwave-labs/openwave/issues/200) | sibling; PARKS behind the SAME blocker (g, delta + the potential). The full production-engine build + Duda's round-1 review live there |
| PMNS mixing structure (precursor, CLOSED) | issue [#199](https://github.com/openwave-labs/openwave/issues/199) | established the rotation-group structure (neutrino oscillation = SO(3); quarks = SU(3)/CKM) + the `delta_CP = 180deg` prediction. This task does the DYNAMICAL field-theoretic derivation #199 did not |
| Effective Dirac / Standard-Model Lagrangian | issue [#197](https://github.com/openwave-labs/openwave/issues/197) | the QED-Lagrangian connection Duda's parameter argument rests on (g <-> the `mc^2` rest-mass term, delta <-> the quantum-phase Dirac term) |
| Scorecard row | [`../../../../MODELS.md`](../../../../MODELS.md) "neutrino sector" | update when this task lands |

## Dr. Duda's round-2 reply (verbatim, 2026-06-20, models-of-particles list)

```text
Dear Rodrigo,

Thank you for this initiative, but to convince mainstream we need really serious
simulation, e.g. to raport in article.

E.g. just looked into .../sandbox_v9/m5_9_1_lepton_mass_law.py and see

DELTA = 0.3 , while we have discussed it should be ~10^-10 to recreate QED contribution
of quantum phase ... what might be reason of huge values you get:
   [image: L_QED = -hbar*c psi-bar gamma^mu D_mu psi - m*c^2 psi-bar psi - F_mu_nu F^mu_nu / 4]

Also we need g, which should be ~10^10, I completely don't see.

Or for D=diag(g,1,delta,0) I see eta = diag(1,1,1,-1) - minus needs to be where g is:
eta = diag(-1,1,1,1) instead.

Or I see Faber potential V = (1-s^2)^3/r0^4) , while here we need for tensor like in
https://en.wikipedia.org/wiki/Landau-de_Gennes_theory ... but details are still to be found.

So maybe it would be better to do less, but more rigorously - e.g. neutrino oscillations
from one side seem the simplest for complete simulations (only topological vortices: no
point-like), from the other oscillation parameters are a guessing game now - this could be
first serious derivation from a deeper theory.

However, first of all they need the parameters (especially g, delta), but also potential -
what seems crucial for regularization in the center of e.g. topological vortex here.

Best,
Jarek
```

The attached image is the QED Lagrangian `L_QED = -hbar*c psi-bar gamma^mu D_mu psi - m*c^2 psi-bar psi - F_mu_nu F^mu_nu / 4`: the dominant rest-mass term `m*c^2 psi-bar psi` and the tiny quantum-phase (Dirac kinetic) term `psi-bar gamma^mu D_mu psi`. The point: g encodes the rest-mass scale (huge), delta encodes the quantum-phase correction (tiny).

## Evaluation of Duda's round-2 points

| Point | Assessment | Verdict |
| --- | --- | --- |
| `delta` should be ~`1e-10`, not 0.3 ("to recreate QED contribution of quantum phase") | Decisive. delta encodes the quantum-phase (Dirac) term, a ~`1e-10` correction relative to the rest mass. At delta=0.3 we made a `1e-10` effect a 30% structure -> huge values, no hierarchy. Applies to the toy AND the production build (both used delta=0.30). | ✅ valid, the core error |
| `g` should be ~`1e10`, "I completely don't see it" | Same regime error. g is the rest-mass (`mc^2`) scale, so it must be huge; we used `LC_G = 8.0`. | ✅ valid |
| `eta` minus must be on the g axis | Convention. Duda writes `D = diag(g,1,delta,0)` (g first) -> `eta = diag(-1,1,1,1)`. Our engine puts g at index 3 with `eta = diag(1,1,1,-1)`, so the minus IS on g (self-consistent); only the axis ORDERING differs. He inferred this from the toy (which has no g axis at all). Worth confirming, not an error in the production build. | ✅ clarify (notation, not a bug) |
| Faber `V = (1-s^2)^3/r0^4` too simple; need an LdG TENSOR potential | He reviewed `sandbox_v9/m5_9_1_lepton_mass_law.py` (the original numpy TOY). The production build [`sandbox_v9/m5_9_4_engine_lepton.py`](sandbox_v9/m5_9_4_engine_lepton.py) already uses the LdG tensor potential `a*Tr(M^2) - b*Tr(M^3) + c*(Tr M^2)^2`. The exact form that regularizes the vortex CORE is still the open piece ("details to be found"). | ✅ partly addressed (point him to the production build), partly open |
| "Do less but more rigorously" -> neutrino oscillations | Strategic. Topological vortices only (no point-like) = simpler complete simulation; oscillation parameters are unmeasured-from-theory -> a derivation is novel + article-worthy. ADOPTED (this task). | ✅ adopt |
| "First of all they need the parameters (g, delta), but also potential" | The shared blocker. Neutrinos and the parked #200 both need g, delta, and the vortex-core regularizing potential, in numbers. | ✅ the precondition |

**Net.** The load-bearing correction is the PARAMETER REGIME: g ~ `1e10` (rest-mass scale), delta ~ `1e-10` (quantum-phase correction), a ~`1e20` separation we had compressed to a factor of ~27 (8 vs 0.3). That is almost certainly why we got huge values and no hierarchy. Duda answered exactly the question we asked (what g, delta to use), and redirected the first rigorous deliverable to the cleaner neutrino target.

## The numerical reality (a real obstacle, worth telling Duda)

`g ~ 1e10` and `delta ~ 1e-10` in the same tensor is a ~`1e20` dynamic range. That exceeds f32 (the engine's current precision) and stresses f64; energy products would lose all precision. A serious simulation at this regime needs one of: **non-dimensionalization** that factors the scales out of the dynamics, a **perturbative delta expansion** (around delta=0), or extended precision. This shapes the whole build and partly explains why this "serious simulation" has not been done. It is itself a finding to raise with Duda (does he expect a non-dimensionalized formulation, or perturbative delta?).

## Round 3 (2026-06-21): the refined scope , 4 PMNS parameters, closed vortex loops, masses deferred

Duda's round-3 reply refines #236 into a sharp, article-grade target and confirms the parameter search is a genuine open problem (not his to hand over).

### Dr. Duda's round-3 reply (verbatim, 2026-06-21, models-of-particles list)

```text
Dear Rodrigo,

Thank you, indeed, while there are suggestions for delta~10^-10, g~10^10, finding the exact
ones turned out surprisingly difficult, even worse for details of potential - you should
start here, finding these parameters/details for agreements.

In .../sandbox_v9/m5_9_4_engine_lepton.py I still see eta = diag(1,1,1,-1), while -1 should
be where g is, and I see you also use diag(g,1,delta,0) - so should be in first coordinate:
eta = diag(-1,1,1,1).

Regarding neutrino mass, here we assume neutrino as closed loops of topological vortex - could
be huge, like this measured 6.2 pm in https://www.nature.com/articles/s41586-024-08479-6
However, mass of such loop should depend approximately linearly on its length - should be
rather seen as density: mass/length.
This way we can understand why oscillation between neutrinos of different masses do not violate
energy conservation - by also changing length of the loop.
I suspect neutrino masses used in oscillation formulation might be these densities per length,
but it also might be more complicated.

So I would leave masses for later, and first focus on the basic 4 parameters - if writing
convincing article able to pass peer review, this already would be huge: [image]

Best regards,
Jarek
```

The attached image is OUR OWN [#199](https://github.com/openwave-labs/openwave/issues/199) result (the SO(3)/TBM PMNS scorecard vs NuFIT 6.0, reproduced below). Duda is pointing us back at our own validated deliverable and saying: make it rigorous and finish it.

### Evaluation of round-3 points

| Point | Assessment | Impact |
| --- | --- | --- |
| "finding the exact g, delta turned out surprisingly difficult, even worse for the potential , you should start here" | The blocker is NOT "Duda hands us numbers". Even he lacks the exact parameters/potential. It is a genuine open SEARCH we drive: find the `(g, delta, potential)` that AGREE with data. | the blocker reframes from "wait for Duda" to "an active parameter + potential SEARCH, the first work" |
| eta should be `diag(-1,1,1,1)` for `diag(g,1,delta,0)` | CONFIRMED a notation bug in `m5_9_4`'s docstring: it writes `diag(g,1,delta,0)` (g first) but the engine orders `diag(1,delta,0,g)` with g AND the Minkowski minus both at index 3 (verified in `medium.py` + `signed_dot4`). Physics is correct; the docstring misleads. | fix the `m5_9_4` docstring; adopt ONE convention , Duda's `diag(g,1,delta,0)` + `eta=diag(-1,1,1,1)` , for the new neutrino code |
| neutrino = **closed loop** of topological vortex (could be ~6.2 pm, Nature s41586-024-08479-6) | The object model for #236: a CLOSED vortex loop, not point-like, not open. | defines the geometry of the neutrino simulation |
| mass ~ loop **length** => mass is a **density** (mass/length); oscillation = the loop CHANGES length (so different-mass flavours conserve energy) | The mass / oscillation mechanism. The "neutrino masses" in oscillation formulas may be per-length densities. | the mass approach (deferred); explains energy conservation across oscillation |
| "leave masses for later, focus on the basic 4 parameters" + the image (our #199 scorecard) | NARROWS the deliverable to the 4 PMNS MIXING parameters; absolute masses (`Delta m^2`) deferred to a later phase. | the scope cut that makes #236 tractable and article-grade |

### The deliverable, narrowed: the 4 PMNS parameters (our #199 scorecard is the target)

The image is our OWN [#199](https://github.com/openwave-labs/openwave/issues/199) result , the SO(3)/TBM prediction vs NuFIT 6.0. It already nails 3 of the 4:

| Parameter | SO(3)/TBM (#199 predicted) | NuFIT 6.0 (measured, NO) | Status |
| --- | --- | --- | --- |
| `theta_12` solar | 35.26deg (sin^2 = 1/3) | 33.7deg | ✅ tri-maximal |
| `theta_23` atmospheric | 45deg (maximal) | straddles 45deg | ✅ maximal |
| `delta_CP` | 180deg | ~177deg (1sigma [148,215]) | ✅ consistent |
| `theta_13` reactor | 0deg | 8.5deg | ⚠️ **THE SO(3)-BREAKING , the one open item** |

So #236 is the **rigorous completion of #199**: (a) DERIVE the TBM structure from the closed-vortex-loop FIELD THEORY (vs #199's group-structure argument), and (b) crack the one open item, `theta_13 = 8.5deg`, the SO(3)-breaking. An article on these 4 (the 3 that work + the `theta_13` derivation) "would already be huge" (Duda). Absolute masses are deferred.

### The connecting hypothesis (the central scientific bet of #236)

The SO(3)-breaking that lifts `theta_13` from 0 to 8.5deg is a natural home for the SMALL parameter `delta` (~`1e-10`, the quantum-phase term). If `delta` breaks SO(3) by the right amount to give `theta_13`, it UNIFIES Duda's two threads (the hard parameter search + the `theta_13` target) into one question. This is the hypothesis #236 is built to test.

## The target: neutrino oscillation parameters from topological-vortex dynamics

### Goal (round-3 scope)

Derive the four PMNS neutrino MIXING parameters (`theta_12`, `theta_23`, `theta_13`, `delta_CP`) from a field-theoretic simulation of CLOSED topological-vortex loops, rigorously completing [#199](https://github.com/openwave-labs/openwave/issues/199): reproduce the 3 it already gives (`theta_12` tri-maximal, `theta_23` maximal, `delta_CP = 180deg`) from the loop dynamics, and DERIVE the one open parameter, `theta_13 = 8.5deg` (the SO(3)-breaking, the central target, hypothesised to come from the small `delta`). Absolute masses (`Delta m^2`) are DEFERRED (the mass-as-loop-length-density model, a later phase). Compare to NuFIT 6.0. The bar is a peer-review-grade article.

### Approach (the rigorous pipeline)

1. **Parameter + potential search FIRST** (Duda: "start here"): an active search for the `(g, delta, vortex-core LdG tensor potential)` that reproduce the #199 scorecard. The exact values are unknown even to Duda , this is our work, not a wait.
2. Closed topological-vortex **loops** as the neutrino object (not point-like, not open); the three flavours = three loop configurations.
3. Oscillation = the SO(3) rotation among the loop configurations (#199), realized dynamically in the field.
4. Reproduce the TBM 3 (`theta_12`, `theta_23`, `delta_CP`) from the loop dynamics (rigorize #199's group-structure argument into a field-theoretic derivation).
5. **Derive `theta_13`** (the SO(3)-breaking) , test whether the small `delta` sources it (the central hypothesis).
6. Numerics: a non-dimensionalized or perturbative-`delta` formulation for the ~`1e20` dynamic range; adopt Duda's convention `diag(g,1,delta,0)` + `eta=diag(-1,1,1,1)` consistently.
7. Document every parameter + configuration + comparison to article standard.

### Definition of done

A documented closed-vortex-loop simulation that (a) reproduces the TBM 3 (`theta_12`, `theta_23`, `delta_CP`) from the loop dynamics, (b) derives `theta_13` (or honestly reports how far off + why), (c) searches and documents the `(g, delta, potential)` that give agreement, (d) compares to NuFIT 6.0, (e) is reproducible by a third party. Absolute masses deferred. The bar is peer-review-grade.

### Sub-tasks: phase-wired (the TBM milestone is the explicit gate)

The work splits into a buildable foundation (no Duda dependency, start NOW), the parameter search that must pass the **TBM gate**, then the `theta_13` crux, then the article. `Depends on` = the hard `blocked by` wiring.

| ID | Sub-task | Phase | Depends on | Effort | Risk |
| --- | --- | --- | --- | --- | --- |
| **N1** | Numerical method: non-dimensionalized / perturbative-`delta` formulation for the ~`1e20` range; adopt `diag(g,1,delta,0)` + `eta=diag(-1,1,1,1)` | A foundation (buildable NOW) | - | ~1-2 wk | LOW-MED |
| **N2** | Closed-vortex-loop sim: the 3 flavour loop configs + the SO(3) rotation dynamics (new `sandbox_v10`) | A foundation | N1 | ~1-2 wk | MED |
| **N3** | Parameter + potential SEARCH: scan `(g, delta, vortex-core LdG tensor potential)` for agreement with the #199 scorecard, guided by the `delta`->`theta_13` hypothesis | B the search | N1, N2 | weeks-months | HIGH |
| ★ **GATE** | **TBM milestone** (criteria below) | gate | N3 | - | - |
| **N4** | Derive `theta_13` (the SO(3)-breaking): test whether the small `delta` (or another mechanism) lifts it 0 -> 8.5deg | C the crux | **TBM GATE** | unknown | HIGH |
| **N5** | Compare to NuFIT 6.0 + the peer-review article (params, potential, configs, the TBM derivation + `theta_13`) | D the deliverable | TBM GATE (TBM-only article possible) + N4 (full) | ~2-4 wk | MED |
| **N6** | Absolute masses (`Delta m^2`) via the mass-as-loop-length-density model | E DEFERRED | N4 | later phase | - |

**★ The TBM gate (the explicit milestone).** N3 has found a `(g, delta, potential)` that REPRODUCES the 3 TBM angles , `theta_12` tri-maximal (sin^2 = 1/3), `theta_23` maximal (45deg), `delta_CP = 180deg` , from the closed-loop FIELD DYNAMICS, not the group argument.

- **Why a gate, not just a step:** it is a real, publishable result on its own (the first field-theoretic derivation of TBM from vortex loops), AND it is the precondition for N4 , there is no point chasing `theta_13` until the loop dynamics demonstrably yield the TBM baseline.
- **PASS:** the 3 angles emerge from the simulation within tolerance -> proceed to N4; N5 may begin a TBM-result article.
- **FAIL (honest):** if no `(g, delta, potential)` reproduces TBM, that is itself a reportable finding (the closed-loop SO(3) dynamics do not yield TBM -> a different structure is needed), and it redirects the effort rather than wasting it.

**Dependency spine:** `N1 -> N2 -> N3 -> [TBM GATE] -> N4 -> N5` ; `N6` deferred behind N4. N1+N2 (the foundation) need no Duda input and are the first move; N3 (the search) is where Duda's intuition narrows the space.

**N1 + N2 are scoped in detail in [`10n1_foundation_scope.md`](10n1_foundation_scope.md)** (planned + run together on "go N1+N2"; all GitHub #236 posting HELD until the whole N-program finishes).

> **✅ N1 + N2 EXECUTED 2026-06-21 (LOCAL, #236 HELD).** Foundation record: [`sandbox_v10/n_foundation_findings.md`](sandbox_v10/n_foundation_findings.md). **N1** PASS , the perturbative-`delta` method recovers the `theta_13` channel to 9.4e-16 where naive f64 returns exactly 0 (the breaking sits 4+ orders below the f64 floor set by the g-sector). **N2** PASS , the closed disclination loop seeds (line tension `dE/dL=+6.74`, bare loop shrinks), the mixing-observable pipeline reproduces the #199 TBM angles exactly (35.264 / 45.000 / 0.000), and the `delta`->`theta_13` channel is exposed. **Central tension for N3:** at O(1) loop coupling `theta_13=8.5 deg` needs `delta~0.15` OR a ~1.5e9 resonant enhancement, NOT `delta~1e-10` , N3 must resolve whether the mixing sees the bare `delta` or `theta_13` rides a near-degenerate gap.

### Effort + timeline (three clocks, do not conflate)

| Clock | Reality |
| --- | --- |
| **Compute (runtime)** | NOT the gate. Per run ~minutes-to-tens-of-minutes on the M4 (the existing 4x4 runs were 1-5 min); a full N3 scan is hours-to-overnight. |
| **Engineering (N1, N2, N5)** | Bounded: ~`6-9` focused-work weeks total across the buildable parts. |
| **Science (N3 + N4)** | UNBOUNDED , the open research (find the parameters; derive `theta_13`). This sets the calendar, not compute or engineering. Could converge in weeks, take months, or hit a wall (the FAIL path is still a result). |

Realistic read: a multi-month effort gated by two open research questions (N3, N4), with the **TBM gate** as a real intermediate, publishable win partway through. For external comms: "a months-scale derivation effort with a defined intermediate milestone," never a dated deliverable.

### Risks / honest unknowns

- The parameter + potential search is the hard FIRST problem , Duda lacks the exact values too, so it is OUR open search, not a wait. Nothing rigorous runs until we find a `(g, delta, potential)` that reproduces the scorecard.
- The central hypothesis (the small `delta` sources `theta_13`, the SO(3)-breaking) may fail , `theta_13` could need a different mechanism.
- The `1e20` dynamic range is a genuine numerical obstacle (non-dimensionalization or perturbative delta required).
- Whether the closed-loop vortex dynamics actually yield the measured mixing angles is the open scientific question (the point of the task).
- Article standard means full reproducibility + documented parameters/potential/configurations.

### Stage / tracking

OpenWave issue [#236](https://github.com/openwave-labs/openwave/issues/236) (board: **In progress**; #200 moved to **Next**). Multi-session. Total invisibility upheld (public OpenWave physics only). Research body lands here (this plan + new scripts + a findings doc); the issue gets the distilled writeup updated as the work converges.

## Open questions for Dr. Duda (round 3, content bullets for the reply, not a draft)

The blocker is now an active search (Duda lacks the exact values), so these are leads to ask AND search ourselves:

- The parameter + potential search: any leads on the `(g, delta, vortex-core LdG tensor potential)` region that reproduces the #199 scorecard? (We will drive the search; his intuition narrows it.)
- The closed loop: the loop topology, and how the three flavours map to loop configurations; is oscillation literally the loop CHANGING LENGTH (the mass-density picture)?
- `theta_13` (the SO(3)-breaking): is it expected to come from the small `delta`, or another source? (the central hypothesis)
- Confirm the scope: the 4 PMNS mixing parameters first, absolute masses (as loop-length densities) later.
- The eta convention: we will ADOPT his `diag(g,1,delta,0)` + `eta=diag(-1,1,1,1)` for the new neutrino code, and fix the `m5_9_4` docstring to match the engine's actual index-3 ordering (the physics there is already correct; only the docstring's `diag(g,1,delta,0)` notation misleads).

## GitHub issue (#236, filed 2026-06-20)

Filed as [#236](https://github.com/openwave-labs/openwave/issues/236) on `openwave-labs/openwave` (labels `theory` + `help wanted`), board status **In progress**; #200 moved to **Next**. The issue carries the distilled TASK PLANNING (Goal + Approach + DoD + Sub-tasks); links to #199 (precursor), #200 (sibling), #197 (Dirac); the blocker (parameters + potential) stated in the body. No SABER references (total invisibility). The detailed plan stays here in `10a`; the issue's distilled writeup is updated as the work converges.

---
_Local plan (research/10a_neutrino_oscillations.md; research body in research/, scripts to a new sandbox when the parameters are locked). Iterate here; publish the distilled version to the new issue only after it is solid and committed. Cross-refs: [`9a_lepton_mass_planning.md`](9a_lepton_mass_planning.md), #199, #200, #197, [`../../../../MODELS.md`](../../../../MODELS.md)._
