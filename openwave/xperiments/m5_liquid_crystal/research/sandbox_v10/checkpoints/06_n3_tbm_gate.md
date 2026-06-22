# N3 checkpoint 06 , ★ TBM GATE PASS (TBM angles from closed-loop field dynamics)

`n3_search.py` (summary `n3_search_summary.json`). The ★ TBM milestone (10a § "Sub-tasks", the explicit
gate) is PASSED: the 3 tribimaximal angles emerge from the closed-loop LC field theory, not the group
argument.

## The gate

| Quantity | Value |
| --- | --- |
| magic crossing `alpha*` (reference geom, n=48) | 0.8205 rad = **47.01 deg** |
| theta12 at alpha* | **35.2644 deg** (trimaximal, sin^2=1/3) |
| theta23 at alpha* | **45.0000 deg** (maximal) |
| theta13 at alpha* | **0.00 deg** |
| max err vs exact TBM | **0.0000 deg** -> PASS (tol 0.05) |
| mass spectrum (eigenvalues) | [1.692e3, 1.942e3, 2.851e3] (arb units) |

## Robustness (81 geometries: n in {40,48,56} x R_loop {7,9,11} x core_vox {1.5,2,2.5} x kappa {0,0.25,1})

| Metric | Result |
| --- | --- |
| magic crossing found | **81/81** |
| gate PASS (<0.05 deg) | **81/81** (TBM err median 0.0000, max 0.0081 deg) |
| alpha* range | [31.97, 63.89] deg, median 50.28 deg |
| nearest candidate magic angle | arccos(1/sqrt3)=54.7 deg (tetrahedral), off by 0.078 rad (4.5 deg) |

So the TBM gate is a ROBUST locus (every geometry has a magic crossing where TBM emerges), not a
fine-tuned coincidence. alpha* drifts with the geometry (not a single universal angle), so it is a
geometry-dependent OUTPUT, near but not exactly the tetrahedral angle.

## Honest decomposition (what is input vs output)

| Angle | Origin | input or output |
| --- | --- | --- |
| theta23 = 45, theta13 = 0 | the mu-tau MIRROR arrangement (e on axis, mu/tau = Rx(+-alpha) mirror pair) | geometric INPUT (a symmetry choice) |
| theta12 = 35.264 (trimaximal) | the MAGIC crossing of the loop-overlap matrix at alpha* | genuine OUTPUT of the field theory |

This is a field-theoretic REALIZATION of TBM: the mu-tau + magic structure (the de-risk theorem) is
carried by the closed-loop LC overlaps, with the magic condition fixing a definite geometric alpha*. The
novelty vs #199 (pure SO(3) group argument): a concrete field computation with a definite magic angle and
a mass spectrum. Be honest in the writeup that mu-tau is an assumed arrangement, magic is the derived
condition.

## Next: S3 `n3_theta13.py` , the crux (resolve the central tension)
theta13 = 0 holds by the EXACT mu-tau mirror. To source theta13 = 8.5 deg, break mu-tau. Connecting
hypothesis (10a): the LC delta (index-2 eigenvalue) sources it. Physical mechanism to test: make the loops
BIAXIAL (eigenvalues 1, delta, 0 with a secondary delta-twist director, Duda's D=diag(g,1,delta,0)); a
mu-tau-odd secondary winding gives theta13 proportional to delta. D4 sets the gain (~46 deg per unit
mu-tau breaking -> need effective breaking ~0.18 for 8.5 deg). Central tension: does delta=1e-10 produce
that via a RESONANT near-degenerate mass gap, or does theta13 need an O(0.1) effective breaking? S3 maps it.
