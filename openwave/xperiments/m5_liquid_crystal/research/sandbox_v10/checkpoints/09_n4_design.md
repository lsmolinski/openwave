# N4 design , theta13 + delta_CP from a CHIRAL loop coupling (fit -> prediction)

Continues N3 (TBM gate passed; theta13 = a fit O(0.1) mu-tau asymmetry, chiral in origin). N4 bet: the
SAME chirality sources theta13 AND a real delta_CP != 180, both from a CHIRAL coupling , the CP-odd
(Lifshitz / cholesteric) term the real symmetric N3 mass matrix was missing.

## The mechanism

Real symmetric mass matrix -> U real -> delta_CP in {0,180}, theta13 needs a free mu-tau asymmetry (N3).
A chiral coupling makes the matrix COMPLEX HERMITIAN:

```text
  M_H = M_real_sym  +  i * g_chiral * C
  C_ab = INT [ <dx dM_a, dy dM_b> - <dy dM_a, dx dM_b> + (yz) + (zx) ]   antisymmetric, reflection-ODD
  C antisymmetric & real -> i*C Hermitian; C is the curl-like (Lifshitz) chiral overlap (handedness).
  U = complex eigenvectors of M_H -> theta13 = |U_e3|^2, delta_CP via Jarlskog J = Im(U_e1 U_mu2 U*_e2 U*_mu1)
```

The chirality is sourced by the secondary delta-twist SCREW (same handedness on the loops = a global
loop handedness). With g_chiral=0 (or no screw) -> recover N3 (real, delta_CP in {0,180}).

## What N4 tests

| Q | test |
| --- | --- |
| does the chiral term give a genuine delta_CP (!= 0,180)? | scan g_chiral / screw chi, compute Jarlskog -> delta_CP |
| does it ALSO source theta13 (without a separate mu-tau asymmetry)? | measure theta13 vs the chiral knob at the magic gate |
| can ONE handedness give theta13~8.5 AND delta_CP~212 together? | 2D scan; look for a simultaneous match |
| is the handedness TOPOLOGICAL (predicted, not fit)? | relate the screw to a writhe/linking integer (if a match exists) |

## Honest fallback
If the chiral term gives delta_CP but theta13 still needs the separate mu-tau asymmetry, that is the
reportable result (theta13 and delta_CP have related-but-distinct chiral origins). Either way N4 adds the
CP sector the scorecard was missing. NuFIT 6.0 NO target: theta13=8.56, delta_CP=212.

## File: n4_chiral.py (sandbox_v10). Convention index-0. Headless f64, 16-core. LOCAL (#236 HELD).
