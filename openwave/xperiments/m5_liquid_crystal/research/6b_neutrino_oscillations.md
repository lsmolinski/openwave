# Neutrino oscillation parameters from topological-vortex dynamics: build plan (#236)

## Status

Issue [#236](https://github.com/openwave-labs/openwave/issues/236) (OpenWave board: **In progress**). Spun off on Dr. Duda's 2026-06-20 recommendation to "do less but more rigorously": pivot the first serious, article-targeted deliverable from charged-lepton masses to **neutrino oscillations**, because the geometry is simpler (only topological vortices, no point-like defects) and the payoff is higher (the oscillation parameters are a guessing game today, so a derivation from the deeper theory would be genuinely novel).

This doc is the LOCAL plan. It captures Duda's round-2 reply verbatim, evaluates it, scopes the neutrino task, and records the shared blocker. Nothing is published to the issue or committed until the local work is solid (the same discipline as [`6a_lepton_mass_planning.md`](6a_lepton_mass_planning.md)).

| Related work | Link | Relationship |
| --- | --- | --- |
| Charged-lepton masses (PARKED) | [`6a_lepton_mass_planning.md`](6a_lepton_mass_planning.md) · issue [#200](https://github.com/openwave-labs/openwave/issues/200) | sibling; PARKS behind the SAME blocker (g, delta + the potential). The full production-engine build + Duda's round-1 review live there |
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

## The target: neutrino oscillation parameters from topological-vortex dynamics

### Goal

Derive the neutrino oscillation parameters, the mixing angles (`theta_12`, `theta_13`, `theta_23`), the mass-squared differences (`Delta m^2_21`, `Delta m^2_31`), and the CP phase `delta_CP`, from a complete field-theoretic simulation of topological vortices in the M5 substrate, at the CORRECT `(g, delta)` regime with a proper LdG tensor core-regularizing potential, and compare to NuFIT 6.0. Build on #199's SO(3) rotation structure (oscillation = a 3-axis flavour rotation) and its `delta_CP = 180deg` prediction.

### Approach (the rigorous pipeline)

1. Topological vortices only (no point-like defects), the simplest complete-simulation geometry (Duda's reason for the pivot).
2. Set `(g, delta)` to the regime Duda specifies (g ~ `1e10`, delta ~ `1e-10`), handled by a non-dimensionalized or perturbative-delta formulation (the `1e20` dynamic range cannot be carried naively).
3. Regularize the vortex core with the LdG tensor potential Duda specifies (the open "details to be found").
4. Realize the three flavour states as vortex configurations; the oscillation = the SO(3) rotation dynamics (#199).
5. Derive the mixing angles + mass-squared differences from the vortex energetics/dynamics; compare to NuFIT 6.0; the `delta_CP = 180deg` (#199) is the cross-check.
6. Document every parameter + the configurations + the comparison, to article standard (Duda's explicit bar).

### Definition of done

A documented topological-vortex simulation that (a) runs at the locked `(g, delta)` regime with the locked potential, (b) sets up the three flavour vortex states, (c) derives the oscillation parameters (or, honestly, how far off and why), (d) compares to NuFIT 6.0 with the `delta_CP = 180deg` cross-check, (e) is reproducible by a third party. Honest scope: matching every parameter in one pass is not guaranteed; the convincing deliverable is the rigorous, documented derivation framework.

### Sub-tasks (sequence)

1. **BLOCKER first**: lock `g`, `delta`, and the vortex-core LdG tensor potential with Duda (numbers + form). [precondition]
2. Non-dimensionalized / perturbative-delta formulation handling the `1e20` dynamic range. [the numerical method]
3. Topological-vortex setup for the three flavour states (single-SO(3) per #199). [geometry]
4. Oscillation dynamics -> the mixing angles + mass-squared differences. [the result]
5. Compare to NuFIT 6.0 + the `delta_CP = 180deg` cross-check (#199). [validation]
6. Documentation to article standard (parameters, potential, configurations, comparison). [Duda's bar]

### Risks / honest unknowns

- The blocker is real: `g`, `delta`, and the potential form are still unspecified in numbers; nothing rigorous runs until Duda pins them.
- The `1e20` dynamic range is a genuine numerical obstacle (non-dimensionalization or perturbative delta required).
- Whether topological-vortex dynamics actually yield the measured oscillation parameters is the open scientific question (the point of the task).
- Article standard means full reproducibility + documented parameters/potential/configurations.

### Stage / tracking

OpenWave issue [#236](https://github.com/openwave-labs/openwave/issues/236) (board: **In progress**; #200 moved to **Next**). Multi-session. Total invisibility upheld (public OpenWave physics only). Research body lands here (this plan + new scripts + a findings doc); the issue gets the distilled writeup updated as the work converges.

## Open questions for Dr. Duda (round 2, content bullets for Rodrigo's reply, not a draft)

- The precise `g`, `delta`: the values/scaling and the derivation, why g ~ `1e10`, delta ~ `1e-10` (is the constraint `g*delta ~ 1`, or `g/delta ~ 1e20`?). Confirm the mapping to the QED Lagrangian: g <-> the `mc^2` rest-mass term, delta <-> the quantum-phase Dirac term.
- The vortex-core potential: the specific LdG TENSOR form that regularizes the center (beyond the general Landau-de Gennes reference), the "details to be found".
- The numerics: given the `1e20` dynamic range, does he expect a non-dimensionalized formulation or a perturbative `delta` expansion?
- The neutrino setup: confirm the topological-vortex configuration for the three flavours, and that the oscillation is the single-SO(3) rotation dynamics (#199); do the mixing angles come from the vortex geometry or the potential?
- Validation: NuFIT 6.0 as the target (angles + `Delta m^2` + `delta_CP`); confirm `delta_CP = 180deg` (#199) is the cross-check.
- Point him to the PRODUCTION build ([`sandbox_v9/m5_9_4_engine_lepton.py`](sandbox_v9/m5_9_4_engine_lepton.py)), not the `m5_9_1` toy he opened: it already uses the LdG tensor potential and the 4x4 g-axis with the Minkowski signature on g.

## GitHub issue (#236, filed 2026-06-20)

Filed as [#236](https://github.com/openwave-labs/openwave/issues/236) on `openwave-labs/openwave` (labels `theory` + `help wanted`), board status **In progress**; #200 moved to **Next**. The issue carries the distilled TASK PLANNING (Goal + Approach + DoD + Sub-tasks); links to #199 (precursor), #200 (sibling), #197 (Dirac); the blocker (parameters + potential) stated in the body. No SABER references (total invisibility). The detailed plan stays here in `6b`; the issue's distilled writeup is updated as the work converges.

---
_Local plan (research/6b_neutrino_oscillations.md; research body in research/, scripts to a new sandbox when the parameters are locked). Iterate here; publish the distilled version to the new issue only after it is solid and committed. Cross-refs: [`6a_lepton_mass_planning.md`](6a_lepton_mass_planning.md), #199, #200, #197, [`../../../../MODELS.md`](../../../../MODELS.md)._
