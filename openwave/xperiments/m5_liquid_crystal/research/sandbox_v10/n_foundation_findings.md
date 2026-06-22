# N1 + N2 foundation , findings (the precision-safe closed-loop neutrino machine)

> **LOCAL / HELD.** This is the foundation phase (N1 + N2) of the neutrino N-program, OpenWave
> issue [#236](https://github.com/openwave-labs/openwave/issues/236). Per
> [`../10n1_foundation_scope.md`](../10n1_foundation_scope.md) § WORKFLOW NOTE, all GitHub #236
> posting is HELD until the WHOLE N-program (N1-N5) finishes; #236 stays `In progress`, nothing
> is posted at the end of N1+N2. This doc is the local foundation record.

Scope: [`../10n1_foundation_scope.md`](../10n1_foundation_scope.md) · master plan:
[`../10a_neutrino_oscillations.md`](../10a_neutrino_oscillations.md) · the SO(3)/TBM target:
[#199](https://github.com/openwave-labs/openwave/issues/199) · engine reused:
[`../../medium.py`](../../medium.py), [`../../engine2_pde.py`](../../engine2_pde.py),
[`../../engine1_seeds.py`](../../engine1_seeds.py), [`../../engine3_observables.py`](../../engine3_observables.py).

## Headline

The foundation machine is built and self-validated. Two results carry the program forward, and
one honest tension is surfaced for the next phase:

| # | Result | Status |
| --- | --- | --- |
| 1 | The SO(3)-breaking (theta_13) channel is computable to **machine precision** via the perturbative-delta method; the naive f64 path **annihilates it** (returns exactly 0). | ✅ N1 PASS |
| 2 | The closed disclination loop seeds + carries a measurable **line tension** (bare loop shrinks); the mixing-observable pipeline **reproduces the #199 TBM angles exactly** and **exposes the delta -> theta_13 channel**. | ✅ N2 PASS |
| 3 | At O(1) loop coupling, theta_13 ~ 1e-8 deg for delta~1e-10 vs the observed **8.5 deg** , a ~**1.5e9** gap. N3 must resolve whether delta is really 1e-10 in the mixing or theta_13 rides a resonant near-degeneracy. | ⚠️ open, for N3 |

## N0 , the M5 engine connection (VERIFIED, not transcribed)

The N1/N2 scripts are standalone numpy , they re-implement the M5 LdG energy functional in f64
(the engine is f32 and cannot carry the g~1e10 / delta~1e-10 range; it also uses the index-3
time axis vs Duda's index-0; and the closed loop is new geometry). `n0_engine_equivalence.py`
proves that re-implementation is faithful, at the engine's toy scales (g=8, delta=0.30):

| Test | Result | Verdict |
| --- | --- | --- |
| **A , fidelity** | engine `compute_energyH_density_M` = `5.41607e-6` vs numpy port = `5.41626e-6`, rel diff **3.4e-5** | port reproduces the engine to f32 rounding -> the LC connection is verified ✅ |
| **B , convention** | signed curvature + LdG potential: engine index-3 vs Duda index-0 agree **bit-identical (rel 0.0)** | index-3 and index-0 are the SAME physics (relabel `sigma=[3,0,1,2]`) ✅ |

Consequence (also the convention-housekeeping answer): the engine's index-3 ordering is not
"wrong" , it is physically identical to Duda's index-0. The fix for peer-review readability is a
convention note + a consistency audit (catch any docstring that contradicts the code, the
m5_9_4 bug class), NOT a risky rewrite of the validated engine.

## N1 , the precision-safe numerical method

### The problem (Duda 2026-06-21, made concrete)

`M = O . diag(g,1,delta,0) . O^T`, g ~ 1e10, delta ~ 1e-10. The energy functional (read from
`engine3_observables.compute_energyH_density_M` + `engine2_pde.V_M`) mixes g^2 ~ 1e20 with
delta^2 ~ 1e-20 in one float sum. f64 ULP at the g-sector swamps everything below it, and the
**entire delta-sector , the SO(3)-breaking, the theta_13 channel , lives there**.

### The method

| Technique | What it does |
| --- | --- |
| Non-dimensionalization | g is a fixed background scale; the energy grades cleanly in powers of g (E_0 ~ g^4, E_1 ~ g^3, E_2 ~ g^2 in bare units) |
| Perturbative-delta | `D(delta)=D_0+delta.D_1` linear => `M=M_0+delta.M_1` exact => the polynomial energy expands exactly via the delta-graded commutators `C=C_0+delta.C_1+delta^2.C_2`, each order summed separately, never against the 1e39-scale E_0 |
| Convention | Duda index-0, `D=diag(g,1,delta,0)`, `eta=diag(-1,1,1,1)` (g and the metric minus both at index 0) |

### The cancellation test (3-way, PASS)

E_1 = the O(delta) SO(3)-breaking coefficient, computed three ways on a tiny grid (n=6):

| Method | E_1 | rel err vs reference |
| --- | --- | --- |
| naive f64 finite difference | `0.000e+00` (annihilated) | `1.0` |
| graded perturbative (the method) | `7.4761884549e+28` | **`9.4e-16`** |
| mpmath 50-digit reference | `7.4761884549e+28` | , (gold truth) |

Why naive dies (full grid n=20): the delta-graded orders are `E_0:E_1:E_2 = 3.16e39 : 2.99e29 :
3.78e20` (each ~g apart). The breaking's energy signal `delta_phys*|E_1| = 2.99e19` sits below
the ULP of E_0 (`6.04e23`) , so any single-accumulator f64/f128 derivative underflows to 0.

**Payoff**: E_1 IS the theta_13 channel ([`../10a`](../10a_neutrino_oscillations.md) § "the
connecting hypothesis"). N1 proves it is numerically reachable; N4 (theta_13) is on solid ground.

## N2 , the closed-loop sim + mixing observables

### (A) closed-loop seeder + line tension

A circular disclination loop (radius R, z=0 plane): the director winds by q*psi (q=1/2, the
stable nematic disclination) around the loop core, blended to the background far out. M =
block-diag(g [time, index 0], M_spatial [indices 1,2,3]), M_spatial uniaxial with principal axis
= director. Signed energy E(R) over R = 5..13 vox gives a clean monotone **line tension
dE/dL = +6.74** -> the bare loop has a collapse force and shrinks. Stabilization (a twist, a
boost-dressing GEM dip cf the #200 boost-GEM result, or a balancing mechanism) is the honest
open engineering problem for N3+. (Absolute E carries the bulk LdG vacuum offset; the slope is
the offset-independent physical quantity.)

### (B) observable pipeline + the delta -> theta_13 channel

The PDG angle extraction (`sin^2 th13=|U_e3|^2`, `sin^2 th12=|U_e2|^2/(1-|U_e3|^2)`,
`sin^2 th23=|U_mu3|^2/(1-|U_e3|^2)`) feeding on the tribimaximal matrix reproduces the #199
symmetric limit exactly: **theta12=35.264, theta23=45.000, theta13=0.000 deg**. The reactor
channel is then activated by a two-level (degeneracy-lifting) model
`tan(2 theta13)=2 delta V_loop/gap` (theta13 -> 0 as delta -> 0). The channel is finite and
exposed; the loop couplings `V_loop, gap` are the inputs the N3 search computes.

### ⚠️ The central tension (for N3)

With O(1) loop coupling, theta_13 is linear in delta: ~1e-8 deg at delta~1e-10. The observed
8.5 deg needs delta ~ 0.15 OR a coupling enhancement ~1.5e9 (a near-degenerate gap that
resonantly amplifies a tiny delta). **This is the sharpest question handed to N3**: is the
effective delta in the mixing ~1e-10, or is theta_13 sourced by a resonance? It also connects
back to Duda's regime (delta~1e-10 was his quantum-phase scale) , the mixing may not see the
bare delta directly.

### How N1 and N2 connect

The bare loop is undressed (g in the decoupled time block, not varied) -> no 1e20 range -> N1
is not yet needed. N1 becomes ESSENTIAL once the loop is boost-dressed (the dynamical,
theta_13-relevant regime): the g-axis mixes into the varying spatial structure and the
delta-sector drops below the f64 floor. N1 recovers it to machine precision; N2 built the loop +
the observable that reads it. Together = the machine N3 drives.

## Parameters (reproducibility)

| Parameter | Value | Where |
| --- | --- | --- |
| g (index-0 eigenvalue) | 1.0e10 | `n1` G_SCALE / `n2` |
| delta_phys (index-2) | 1.0e-10 | both |
| eta | diag(-1,1,1,1) (index-0 minus) | `n1` ETA_DIAG / SIGN_MAT |
| LdG (a,b,c) | (-2.0, 0.0, 1.0) representative; exact = N3 (Duda Q7) | `n1` LDG_A/B/C |
| curvature prefactor c^2 | 1.0 (bare units) | `n1` CURV_C2 |
| N1 grid (full / tiny) | 20 / 6 ; mpmath dps=50 | `n1` |
| N2 grid / loop radii | 40 / [5,7,9,11,13] vox ; q=1/2 | `n2` |
| N2 loop coupling V_loop, gap | 1.0, 1.0 (placeholders -> N3) | `n2` |

## Artifacts (all LOCAL, all < 100 KB , no raw data > 1 MB produced, nothing to delete)

| Artifact | Regenerate |
| --- | --- |
| `n0_engine_equivalence.py` + `_summary.json` | `python3 n0_engine_equivalence.py` (imports the live engine, CPU) |
| `n1_precision_method.py` + `_summary.json` + `.png` | `python3 n1_precision_method.py` |
| `n2_closed_loop.py` + `_summary.json` + `.png` | `python3 n2_closed_loop.py` |
| `checkpoints/00_design.md`, `01_n1.md`, `02_n2.md` | progress log (this run) |

## Definition of Done (from [`../10n1`](../10n1_foundation_scope.md)) , check

| DoD item | Status |
| --- | --- |
| N1 precision-safe (cancellation test passes), order-by-order, index-0 | ✅ |
| N2 loop seeds, energy/line-tension measured, mixing observables exposed + finite | ✅ |
| Findings doc reproducible by a third party | ✅ (this doc) |
| Machine READY for N3 to drive (vary g, delta, potential; read observables) | ✅ |

NOT in scope (correctly deferred): the parameter search (N3), reproducing the angles = the TBM
gate, theta_13 = N4, masses = N6.

## Next (N3) , ✅ EXECUTED 2026-06-21 (LOCAL)

Drive this machine: search `g, delta`, the LdG potential, the loop geometry + dressing, to (a) resolve
the central tension and (b) reproduce the three TBM angles from loop dynamics = the ★ TBM gate
([`../10a`](../10a_neutrino_oscillations.md) § "Sub-tasks: phase-wired").

**N3 record: [`n3_findings.md`](n3_findings.md).** ★ TBM GATE PASS , the 3 TBM angles emerge from the
closed-loop field dynamics at the magic crossing alpha* (theta12=35.26 trimaximal, theta23=45, theta13=0),
robust across 81/81 geometries. The central tension is RESOLVED + sharpened: a mu-tau-symmetric delta gives
theta13 = 0 EXACTLY (any delta) , theta13 = 8.5 deg needs an explicit O(0.1) mu-tau ASYMMETRY, NOT the bare
delta ~ 1e-10 (resonance ruled out). The connecting hypothesis is too simple; theta13 is a separate mu-tau-
breaking, plausibly chiral (the theta13 <-> delta_CP link). GitHub #236 stays HELD.
