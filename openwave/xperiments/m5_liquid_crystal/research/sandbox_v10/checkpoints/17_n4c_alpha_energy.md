# N4c checkpoint 17 , the decisive alpha-energy test (theta_12 prediction-vs-fit)

The single most important review point (both cold readers, independently): is the magic-crossing tilt `alpha*`
energetically SELECTED (`dE/dalpha = 0` there), making `theta_12 = 35.26` a real prediction, or just the tilt
where trimaximal mixing is read off (a fit)?

## Result (`n4c_alpha_energy.py`)

| quantity | value |
| --- | --- |
| magic crossing `alpha*` | 0.81934 rad = **46.945 deg** (theta12=35.2644, theta23=45, theta13=0, TBM err 0) |
| `E_self` flatness (peak-to-peak / mean) | **0.09%** (= discretization noise floor) |
| `E_self` global argmin | **8.6 deg** (small-tilt edge -> degenerate `alpha->0`, three loops coincide) |
| `dE_self/dalpha` at `alpha*` | +3.26e2 (rel 2.76 to characteristic slope) -> NOT stationary at `alpha*` |
| `trace M_ab` feature | broad **MAX at 48.15 deg** (1.2 deg from `alpha*`), 2.9% variation |
| `lambda_min`, `lambda_max` features | 18-21 deg away from `alpha*` (no coincidence) |

## Verdict (honest)

**NEGATIVE on dynamical selection, with a weak structural coincidence.**

- The substrate self-energy `E_self(alpha)` is FLAT to 0.09% -> the substrate is essentially INDIFFERENT to the
  tilt; its true global minimum is at the degenerate small-tilt edge (`alpha->0`, three loops coincide, mixing
  undefined). The loops do NOT relax to `alpha*` by energy minimization. **This directly confirms the reviewers'
  "substrate dropped out of the angles" concern at the self-energy level.**
- The ONE energetic signal is a broad, shallow MAXIMUM of the tight-binding trace `tr(M_ab)` (2.9% variation)
  within ~1 deg of `alpha*`: the three loop displacements are maximally mutually distinct near the magic tilt.
  Suggestive, but a maximum (not a minimum) and weak.
- So `theta_12 = 35.26` is PINNED by ONE scalar magic condition on the energy-overlap matrix (a derived
  geometric locus, NOT a free continuous fit), coinciding with a shallow trace extremum, but it is NOT a
  variationally-selected ground state. It sits BETWEEN "pure fit" and "dynamical prediction": geometrically
  determined, conditional on the mu-tau arrangement.

## Consequence for the scorecard

The honest scorecard downgrades `theta_12` from "PREDICTION" to "geometrically PINNED (magic condition),
conditional on the mu-tau arrangement; not energetically selected." Combined with the mu-tau inputs
(`theta_23`, `theta_13=0`, \|`delta_CP`\|=90 are consequences), the real content is: 1 imposed symmetry + 1
geometrically-pinned angle (`theta_12`) + 1 free coupling (`theta_13`) + `delta_CP` magnitude (sign open).

Artifacts: `n4c_alpha_energy.py`, `n4c_alpha_energy_summary.json`, `n4c_alpha_energy.png`. Next: mass-ratio
falsifier (`n4c_mass_ratio.py`).
