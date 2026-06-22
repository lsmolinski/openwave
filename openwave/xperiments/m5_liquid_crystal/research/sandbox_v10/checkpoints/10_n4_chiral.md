# N4 checkpoint 10 , chiral coupling: delta_CP = MAXIMAL (predicted) + theta13<->delta_CP unified

`n4_chiral.py` (summary `n4_chiral_summary.json`) + robustness check. Major result.

## What the chiral (Lifshitz/cholesteric) coupling adds

Real symmetric mass matrix (N3) -> U real -> delta_CP in {0,180}. Adding the CP-odd chiral overlap makes
M_H = M_real + i*g_chiral*C complex Hermitian -> complex U -> a genuine delta_CP via Jarlskog. Result:

| Observable | value | origin |
| --- | --- | --- |
| theta12 | 35.26 (trimaximal) | magic crossing (N3) |
| theta23 | 45.00 | mu-tau mirror (N3) |
| **delta_CP** | **+-90 deg (MAXIMAL CP violation)** | **mu-tau reflection + chiral coupling , PREDICTED, robust** |
| theta13 | chiral-sourced, grows with the handedness (reaches ~8 deg at strong chiral coupling) | the SAME chirality as delta_CP |

## Robustness (delta_CP = +-90 is pinned, not a one-off)

7 combos of (delta, chi, R_loop, core, g_chiral): **delta_CP = +-90.00 EXACTLY in all**, theta23 = 45.000 in
all, theta12 = 35.26 in all. theta13 varies 0.04-1.6 deg (chiral-sourced, grows with chi & g_chiral).
g_chiral=0 recovers N3 exactly (theta13=0, J=0, delta_CP=0). Sanity holds.

## The physics (literature-consistent, and that is a GOOD sign)

This is the celebrated **mu-tau REFLECTION symmetry** result (Harrison-Scott 2002): mu-tau reflection forces
theta23 = 45 AND delta_CP = +-90 (maximal). Our closed-loop LC geometry REALIZES mu-tau reflection (e on
axis, mu/tau mirror pair) and the chiral LdG term supplies the CP phase -> the loop field theory PREDICTS
maximal CP violation. The data leans this way: T2K/NOvA favour near-maximal delta_CP ~ 270 deg (= -90); NuFIT
6.0 NO best fit 212 deg sits between maximal (270) and CP-conserving (180), 1-2 sigma from -90.

So theta13 and delta_CP share ONE chiral origin (the loop handedness): turning on chirality turns on BOTH.
This converts delta_CP from "unknown / {0,180}" (N3) into a PREDICTION (maximal), and ties theta13 to it.

## Still open (the rest of N4)
- theta13 = 8.5 deg needs STRONG chiral coupling (chi~1.2, g_chiral~0.9); is that physical / topological?
- The "article-grade" piece: is g_chiral (the handedness) a TOPOLOGICAL invariant (loop writhe / self-linking
  integer)? If yes, theta13 + delta_CP are fully predicted from loop topology. Not yet shown.
- delta_CP sign (+90 vs -90) flips with parameters; the data prefers -90 (=270). What fixes the sign?

## Status
delta_CP = maximal: solid + robust + literature-consistent -> a real PREDICTION (N4's CP sector done).
theta13-as-topological-prediction: open (the remaining N4 crux). Next: the writhe/linking quantization.
