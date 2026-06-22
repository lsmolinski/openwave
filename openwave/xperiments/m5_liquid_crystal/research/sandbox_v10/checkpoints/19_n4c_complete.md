# N4c checkpoint 19 , N4c COMPLETE (review response + honest reframe)

The cold-read peer review ([`../10c_AI_reviewers.md`](../10c_AI_reviewers.md)) is triaged
([`../10d_bullet_proof_run.md`](../10d_bullet_proof_run.md)) and answered ([`../10e_findings_N4c.md`](../10e_findings_N4c.md)).

## Two decisive runs

| run | result |
| --- | --- |
| N4c-1 `n4c_alpha_energy.py` | `theta_12` is NOT energy-selected: `E_self` flat to 0.09% (substrate indifferent to tilt), global min at degenerate edge; only a shallow `tr(M_ab)` max within ~1 deg of `alpha*`. `theta_12` is geometrically PINNED (magic), not a dynamical prediction. Confirms "substrate dropped out of the angles". |
| N4c-2 `n4c_mass_ratio.py` | gate spectrum 1:1.148:1.682 -> splitting ratio 5.76 (lam=m) / 4.62 (lam=m^2) vs observed 33.6 -> ~5-7x too compressed. Tension FLAG (the eigenvalue->mass map is the deferred N6 question). |
| N4c-3 `n4c_scorecard.py` | honest provenance + pull plot: `theta_12` +2.3 sigma, `theta_23` +1.7 sigma, `theta_13` set (0), `delta_CP` consistent (huge range). |

## Honest scorecard (now canonical, in 10e)

ONE imposed mu-tau symmetry (generating `theta_23`=45, `theta_13`=0, \|`delta_CP`\|=90) + ONE geometrically
pinned angle (`theta_12`, NOT energy-selected) + ONE free coupling (`theta_13`) + an UNDETERMINED CP sign. Mass
spectrum ~6x too compressed vs data. Materially weaker than the prior "3 predicted", and it survives a hostile read.

## Docs

- `10c_AI_reviewers.md` (verbatim reviews), `10d_bullet_proof_run.md` (triage + plan), `10e_findings_N4c.md`
  (executed response, canonical framing).
- `10b_findings.md` -> renamed `10b_findings_N4b.md` (forward-pointer to 10e added; N0-N4b detail retained).
- `10a` related-work table updated (10e canonical, 10b_N4b detail, 10c/10d added).

## Pending (unchanged, deferred by design)

- N5 (peer-review article), N6 (absolute masses , must resolve the N4c-2 mass-ratio tension).
- All #236 GitHub posting still HELD until the N-program finishes. Nothing posted/committed (git is the user's).

Artifacts all < 135 KB; nothing > 1 MB.
