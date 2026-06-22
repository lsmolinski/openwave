# N4c , bullet-proofing the neutrino findings against cold-read peer review (the plan)

> **LOCAL / HELD.** The triage of the two AI peer reviews ([`10c_AI_reviewers.md`](10c_AI_reviewers.md)) and the
> plan for the N4c response round, OpenWave issue [#236](https://github.com/openwave-labs/openwave/issues/236).
> Executed results + adjusted findings: [`10e_findings_N4c.md`](10e_findings_N4c.md). Reviewed doc (now renamed):
> [`10b_findings_N4b.md`](10b_findings_N4b.md). Master plan: [`10a_neutrino_oscillations.md`](10a_neutrino_oscillations.md).

## Why this round

Two independent AI reviewers, with **zero** context on M5/OpenWave/#199/#197, reviewed the consolidated findings.
Both, independently, converged on the same central critique and the same single decisive test. About **70% of
both reviews is valid and survives full context** (it is the same ground Duda will probe), and ~30% is a
framework/scope gap (they reviewed a complete, SM-embedded, quantum theory of neutrinos; the document is a
focused mixing-parameter derivation in an alternative substrate framework that defers masses by design). N4c
acts on the valid 70% before any exchange with Duda, so the doc is honest to exactly the critiques he will raise.

## The one decisive test (both reviews, independently)

> **Is the magic-crossing tilt `alpha*` energetically selected (`dE/dalpha = 0` at `alpha*`), or is it just the
> tilt where trimaximal mixing is read off?**

If the loop energy is stationary at `alpha* = 46.94 deg`, then `theta_12 = 35.26` is a genuine prediction and the
substrate does real work. If the energy is flat in `alpha` or minimized elsewhere, then `theta_12` is also a fit,
and the honest summary becomes "an SO(3)/mu-tau ansatz reproduces TBM as group theory requires, with M5 dynamics
playing no discernible role in the angles." This is the linchpin; it is run **first**, and the result determines
how `theta_12` is framed in 10e.

## Triage , real gaps to act on (survive full context)

| Critique | Raised by | Verdict | N4c action |
| --- | --- | --- | --- |
| `alpha*` not shown energy-selected -> `theta_12` may be a fit | both | valid, decisive | `n4c_alpha_energy.py`: `E(alpha)` + `dE/dalpha` at `alpha*` |
| Scorecard triple-counts ONE mu-tau assumption as 3 predictions | both | valid; our own Q3 contradicts the headline | reframe scorecard (below) |
| `delta_CP` SIGN not predicted: `E(+chi)=E(-chi)` (our N4b-3) | Claude | valid, consistent with our data | downgrade to \|`delta_CP`\|=90; sign open |
| Mass ratios 1:1.15:1.68 checkable NOW, look compressed | Claude | valid, near-term falsifier | `n4c_mass_ratio.py`: ratio vs observed |
| "Hessian" is a large-displacement Gram/overlap; `U=eigvecs` asserted | both | valid | fix terminology; label the bridge a postulate |
| Loop instability (`dE/dL>0`) -> mass matrix on non-stationary configs | both | valid, elevate | promote from Q7; tie to the `alpha`-energy work |
| `delta` + potential independence = "substrate dropped out" of angles | both | valid, uncomfortable | reframe robustness honestly (angles are symmetry-generic) |

## Triage , honesty / framing fixes (doc edits in 10e)

| Fix | Change |
| --- | --- |
| Headline scorecard | stop calling the mu-tau consequences "predictions" (reframe below) |
| NuFIT comparison has no error bars | add 1sigma/3sigma; state `theta_12` ~2sigma off, `delta_CP` consistent |
| "DUNE/HK converging toward 270" | soften: "being tested; some global fits disfavor maximal" |
| "`theta_13` ~ O(10 deg) natural" | keep the "fit" label clean; naturalness is plausibility, not prediction |

## Triage , context gaps (reviewer reviewed a bigger claim than the doc makes)

| Reviewer expectation | Why it is a framework/scope mismatch | Legit residual kept |
| --- | --- | --- |
| Connect to SM (SU(2)_L, Higgs, seesaw) | the framework derives particles as substrate defects, not SM excitations; no Yukawa/seesaw by construction | justify the overlap -> measured-mixing bridge |
| Quantum theory of oscillation `P(a->b)` | the framework is a classical nonlinear field; the field IS the physical object | same bridge question |
| "Why 3 flavors?" | the 3 = SO(3) orientations is the [#199](https://github.com/openwave-labs/openwave/issues/199) starting point; the 3-generation origin is a separate program piece | state plainly the doc does not derive 3 |
| `g~1e10`/`delta~1e-10` ad hoc | the scales are the [#197](https://github.com/openwave-labs/openwave/issues/197) rest-mass / quantum-phase axes; and N4b-4 shows the mixing is independent of them | cite #197; note not load-bearing for mixing |
| "TBM is outdated" | TBM + `theta_13` corrections (TM1/TM2) is a current framework; N4 IS the post-TBM step | minor; the real core (`theta_13` fitted) already owned |

## The honest scorecard 10e must land

| Parameter | Honest status |
| --- | --- |
| `theta_23` = 45, `theta_13` = 0 (baseline) | CONSEQUENCES of the mu-tau-mirror INPUT (one assumption, via Harrison-Scott) |
| \|`delta_CP`\| = 90 | CONSEQUENCE of the same mu-tau-reflection INPUT; SIGN not predicted (handedness-degenerate) |
| `theta_12` = 35.26 | the ONE candidate genuine prediction (magic crossing) , conditional on `alpha*` being energy-selected |
| `theta_13` = 8.56 | one-parameter FIT (`g_chiral`) |

So: one structural assumption + one conditional prediction + one fit. Materially weaker than "3 predicted", but
defensible and exactly what Duda would accept.

## N4c run plan

| # | Step | Output |
| --- | --- | --- |
| 1 | `n4c_alpha_energy.py` , the decisive test: `E(alpha)` (self-energy + tight-binding trace + ground state) over the tilt; locate stationary points; compare to `alpha*` | `n4c_alpha_energy_summary.json` + panel |
| 2 | `n4c_mass_ratio.py` , gate eigenvalues -> splitting ratio under both natural maps (eigval = m, eigval = m^2) vs observed `Dm31^2/Dm21^2 ~ 33.6` | `n4c_mass_ratio_summary.json` + panel |
| 3 | `n4c_scorecard.py` , the honest scorecard figure (predictions vs imposed-symmetry consequences vs fit; NuFIT with error bars) | `n4c_scorecard.png` |
| 4 | Write `10e_findings_N4c.md` , same structure as the reviewed doc (sections, images, links); honest reframe; the two new results; bridge + stability elevated; framework + scope note | 10e |
| 5 | Rename `10b_findings.md` -> `10b_findings_N4b.md`; forward-pointer to 10e; update the 10a cross-ref | rename + edits |

Self-imposed honesty rule for this round: report the `alpha`-energy result as it comes out, including a negative
(no dynamical selection). A negative is a real result and is recorded plainly, not spun.

## Cross-refs

[`10c_AI_reviewers.md`](10c_AI_reviewers.md) (the verbatim reviews) · [`10e_findings_N4c.md`](10e_findings_N4c.md)
(executed) · [`10b_findings_N4b.md`](10b_findings_N4b.md) (reviewed) ·
[`10a_neutrino_oscillations.md`](10a_neutrino_oscillations.md) (master plan) ·
[#236](https://github.com/openwave-labs/openwave/issues/236) (HELD).
