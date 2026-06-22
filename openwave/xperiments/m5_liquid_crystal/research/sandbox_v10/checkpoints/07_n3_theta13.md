# N3 checkpoint 07 , S3 the theta13 crux + central-tension RESOLUTION

`n3_theta13.py` (summary `n3_theta13_summary.json`). The crux resolved, with a sharper (and more honest)
answer than the original connecting hypothesis.

## The three results

| Step | Result |
| --- | --- |
| **S3a (negative control)** | a mu-tau-SYMMETRIC biaxial delta-structure gives **theta13 = 0 EXACTLY** (max 8.4e-14 deg across delta in [1e-3, 1]). The bare SO(3)-breaking delta does NOT source theta13. |
| **S3b (the real source)** | theta13 turns on ONLY with an explicit mu-tau ASYMMETRY eps (secondary delta-twist on the mu loop but not tau). Bilinear **theta13 ~ G*delta*eps**, G ~ 126 deg (clean linear regime delta<=0.1). **8.5 deg needs delta*eps ~ 0.07** (e.g. delta~0.1, eps~0.6) = an O(0.1) breaking. |
| **S3c (resonance)** | at delta=1e-10, tuning a near-degenerate mass gap peaks theta13 at only **1.3e-8 deg**; reaching 8.5 deg needs gap/spectrum ~ 3.3e-10 (absurd fine-tuning). Resonance is ruled out. |

## Resolution of the central tension (the foundation's open question)

The foundation surfaced: theta13=8.5 deg needs delta~0.15 OR a ~1.5e9 resonance, NOT delta~1e-10. S3 sharpens
this into a STRUCTURAL statement:

> A mu-tau-SYMMETRIC delta (any size) gives theta13 = 0 exactly. theta13 = 8.5 deg is a SEPARATE O(0.1)
> mu-tau-ASYMMETRY, not the bare quantum-phase delta ~ 1e-10. The connecting hypothesis ("delta sources
> theta13") is too simple: theta13 needs an explicit mu-tau-breaking, distinct from the SO(3)-breaking delta.

Plausible deeper home: the chiral / CP sector. The symmetric delta-biaxiality is CP-even and gives theta13=0;
a CHIRAL ingredient (handed secondary screw, opposite on the loops) is what breaks mu-tau. This hints
theta13 and delta_CP share a chiral origin , a lead for N4.

## Honest input/output ledger (for the writeup)

| Observable | status | origin |
| --- | --- | --- |
| theta12 = 35.26 (trimaximal) | PREDICTION | magic crossing of the loop overlaps (S2) |
| theta23 = 45, theta13 = 0 (TBM baseline) | PREDICTION | mu-tau mirror symmetry (S2) |
| theta13 = 8.5 | FIT (one free param) | requires an O(0.1) mu-tau asymmetry eps; NOT predicted (free), exactly as in real TBM-based flavour models |
| delta_CP | open | real mass matrix -> delta_CP in {0, 180}; a genuine value needs a complex/chiral coupling (S4 + N4) |

## Numerical note
The theta13(delta,eps) map is clean+monotone for delta<=0.1; at large delta*eps it is non-perturbative
(the TBM baseline drifts, one outlier cell delta=0.2/eps=0.4=32 deg) , the bilinear gain is taken from the
linear regime. Documented; not load-bearing for the conclusion (order-of-magnitude robust).

## Next: S4 `n3_scorecard.py`
Assemble the full NuFIT 6.0 scorecard (theta12/theta23/theta13/delta_CP), state delta_CP=180 in the real
CP-conserving limit (consistent with #199), generate the 3 key plots (magic-crossing TBM gate, theta13(delta,
eps) map, mass spectrum), write the master findings doc `n3_findings.md`.
