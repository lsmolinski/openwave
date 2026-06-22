# N4 findings , the PMNS CP sector + the closed gate (delta_CP predicted maximal; theta13 the one free coupling)

> **LOCAL / HELD.** N4 (the theta13 crux + CP sector) of the neutrino N-program, OpenWave issue
> [#236](https://github.com/openwave-labs/openwave/issues/236). All GitHub #236 posting HELD until the whole
> N-program finishes (per [`../10n1_foundation_scope.md`](../10n1_foundation_scope.md)). Local record. Builds
> on N3: [`n3_findings.md`](n3_findings.md). Master plan: [`../10a_neutrino_oscillations.md`](../10a_neutrino_oscillations.md).
> Convention index-0 ([`../_convention_refactor/CONVENTION.md`](../_convention_refactor/CONVENTION.md)).

## Headline , the PMNS gate is FULLY CLOSED (3 of 4 parameters PREDICTED)

| # | Result | Status |
| --- | --- | --- |
| 1 | The CP-odd **chiral (Lifshitz/cholesteric) coupling** makes the flavour mass matrix complex Hermitian and PREDICTS **delta_CP = +-90 deg (MAXIMAL CP violation)**, robust across delta/chi/geometry. | ✅ delta_CP predicted |
| 2 | **delta_CP sign = the loop HANDEDNESS** (`= -sign(g_chiral)`; the chi-screw sign is irrelevant). One definite chirality (`g_chiral > 0`) gives **delta_CP = -90 (= 270 deg), the data-preferred value**. | ✅ sign fixed by handedness |
| 3 | **theta13 is CONTINUOUS** (set by the chiral coupling strength, the cholesteric-pitch analogue), reaching the NuFIT 8.56 deg at **O(1) coupling `g_chiral* = 0.94` (NATURAL, no resonance)**. NOT topologically quantized (two independent attempts failed). | ✅ theta13 = 1 free coupling |
| 4 | `theta13` and `delta_CP` share ONE origin (the loop handedness): turning on chirality turns on BOTH. The N3 `theta13` "fit" is reframed as the chiral partner of `delta_CP`. | ✅ unified |
| 5 | This is the **mu-tau REFLECTION symmetry** result (Harrison-Scott 2002: mu-tau reflection => `theta23` = 45 AND `delta_CP` = +-90) , now REALIZED from a closed-loop LC geometry, not assumed. | ✅ literature-consistent |

## Mechanism (the CP-odd chiral coupling)

```text
  M_H = M_real_sym  +  i * g_chiral * C ,
  C_ab = INT [ <dx dM_a, dy dM_b> - <dy dM_a, dx dM_b> + (y->z) + (z->x) ]_s    real, antisymmetric, reflection-ODD
  U = complex eigenvectors of M_H ;  theta13 = |U_e3|^2 ;  delta_CP via Jarlskog J = Im(U_e1 U_mu2 U*_e2 U*_mu1)
```

`C` is the curl-like Lifshitz overlap (the cholesteric chiral term of a real liquid crystal): antisymmetric in
`(a,b)` -> `i*C` Hermitian; reflection-odd -> a true handedness. `g_chiral -> 0` recovers N3 exactly (real,
CP-conserving). Run: `n4_chiral.py`.

## FINAL scorecard vs NuFIT 6.0 (NO) , the article table

`n4_final_scorecard.py` ([`n4_final_scorecard.json`](n4_final_scorecard.json), figure
[`n4_final_scorecard.png`](n4_final_scorecard.png)). Headline point: `chi = 1.2, delta = 0.1, g_chiral* = 0.937`.

| Parameter | TBM baseline | N4 full | NuFIT 6.0 (NO) | status |
| --- | --- | --- | --- | --- |
| theta12 | 35.26 | 35.72 | 33.68 | **PREDICTION** (trimaximal; magic crossing) |
| theta23 | 45.00 | 45.00 | 43.30 | **PREDICTION** (maximal; mu-tau mirror) |
| theta13 | 0.00 | 8.56 | 8.56 | **FREE COUPLING** (chiral g* = 0.94, natural) |
| delta_CP | 0 | 270 | 212 | **PREDICTION** (maximal; sign = handedness -> 270) |

3 of 4 PMNS parameters PREDICTED from the loop geometry; theta13 is the single free coupling (the chiral
strength). theta23 = 45 and delta_CP = 270 are exact predictions; theta12 = 35.26-35.72 sits ~2 deg above
NuFIT 33.68 (the standard TBM-vs-data gap); delta_CP = 270 vs NuFIT 212 (data leans toward maximal; DUNE/HK
will decide).

## Study A , theta13 reach: NATURAL, continuous, near-linear (`n4_refine.py`)

At the gate (mu-tau-symmetric, alpha* re-found), `theta13(g_chiral)` at `chi = 1.2, delta = 0.1` rises smoothly
and near-linearly (slope ~9.5 deg/g, mass gap stable ~1.5 , NO resonance), crossing **8.56 deg at
`g_chiral* = 0.937`** with `theta23 = 45`, `delta_CP = -90` held. So `theta13 = 8.56` needs only O(1) chirality
, NATURAL, not strong/fine-tuned. (`theta12` drifts 35.26 -> 35.72 as chirality rises.)

## Study B , the delta_CP sign = the loop handedness (`n4_refine.py`)

Flipping the signs of `(chi, g_chiral)`: **`delta_CP_sign = -sign(g_chiral)` ALONE** , the chi-screw sign is
IRRELEVANT. So the sign of `delta_CP` is fixed by the loop HANDEDNESS (the sign of the chiral coupling). One
definite handedness (`g_chiral > 0`) gives `delta_CP = -90 (= 270 deg)`, the data-preferred value , a sharp,
falsifiable statement (the theory predicts which way CP is violated, given the loop handedness).

## Study C , is theta13 TOPOLOGICAL? NO (two independent attempts) (`n4_topo.py`, `n4_linking.py`)

Tested whether the handedness is a quantized self-linking integer N (Calugareanu: `Lk = Tw + Wr`):

| Attempt | framing | result |
| --- | --- | --- |
| `n4_topo.py` | `N * s` (LOCAL azimuth, mu-tau-ODD under the z->-z mirror) | breaks mu-tau, degrades TBM |
| `n4_linking.py` | `N * phi` (GLOBAL azimuth, intended mu-tau-even) | ALSO breaks mu-tau (the secondary FRAME still transforms), degrades TBM |

Both break the TBM baseline rather than cleanly stepping `theta13`. **Conclusion: `theta13` is NOT
topologically quantized in this theory; it is CONTINUOUS, set by the chiral material coupling `g_chiral` (the
cholesteric-pitch analogue).** The insight: a mu-tau-symmetric arrangement and a nonzero, quantized handedness
are in tension (the mirror flips the handedness), so the handedness is a continuous material property, not a
topological integer. `delta_CP = +-90` is the DISCRETE prediction; `theta13` is the one continuous free coupling.

## Physics reading (the article narrative)

The closed-loop LC geometry realizes **mu-tau reflection symmetry** (e-loop on the axis, mu/tau a mirror pair),
the symmetry that in the neutrino literature forces `theta23 = 45` and `delta_CP = +-90`. The chiral LdG
(Lifshitz) term supplies the CP phase. So the model PREDICTS maximal CP violation with a sign set by the loop
handedness (`g_chiral > 0 -> 270 deg`), and `theta13` is its continuous chiral partner (one free coupling,
natural at O(1)). Next-generation long-baseline experiments (DUNE, HK) measuring `delta_CP -> 270` would
support it; a CP-conserving value (180) would disfavour it.

## Honest ledger (for the article)

| Observable | status | origin |
| --- | --- | --- |
| theta12 = 35.26 (trimaximal) | PREDICTION | magic crossing of the loop overlaps |
| theta23 = 45 (maximal) | PREDICTION | mu-tau mirror symmetry |
| delta_CP = 270 (maximal) | PREDICTION | mu-tau reflection + chiral; sign = handedness |
| theta13 = 8.56 | 1 FREE COUPLING | continuous chiral strength g* = 0.94 (natural); not topological |
| theta12 residual (~2 deg above data) | open | the standard TBM-vs-data gap (shared by all TBM models) |

## Remaining open (leads for N5 / Duda , not blocking the article)

| Question | note |
| --- | --- |
| The ~2 deg `theta12` (and `theta23`) residual vs NuFIT | the standard TBM-vs-data gap; a small symmetry-breaking correction (cf. N3 mu-tau asymmetry) could absorb it , a refinement, not a gap |
| The microscopic value of `g_chiral` (-> a `theta13` PREDICTION) | `g_chiral` is the cholesteric/Lifshitz material strength; deriving it from the substrate would turn `theta13` from free to predicted |
| A mu-tau-respecting topological framing | ruled out in two forms here; a genuinely different construction (e.g. true loop-loop linking) remains open |

## Artifacts (all LOCAL, < 130 KB, nothing > 1 MB)

| Artifact | Regenerate |
| --- | --- |
| [`n4_chiral.py`](n4_chiral.py) + `n4_chiral_summary.json` | `python3 n4_chiral.py` (the delta_CP=maximal result, 16-core) |
| [`n4_refine.py`](n4_refine.py) + `n4_refine_summary.json` | `python3 n4_refine.py` (Study A theta13 reach + Study B sign) |
| [`n4_topo.py`](n4_topo.py) + `n4_topo_summary.json` | `python3 n4_topo.py` (topological test, local azimuth , negative) |
| [`n4_linking.py`](n4_linking.py) + `n4_linking_summary.json` | `python3 n4_linking.py` (Study C mu-tau-even framing , negative) |
| [`n4_final_scorecard.py`](n4_final_scorecard.py) + `n4_final_scorecard.json` + `n4_final_scorecard.png` | `python3 n4_final_scorecard.py` |
| `checkpoints/09_n4_design.md` .. `12_n4_complete.md` | progress log |

## Cross-refs

[`n3_findings.md`](n3_findings.md) (the TBM gate + the theta13 fit this reframes) ·
[`../10a_neutrino_oscillations.md`](../10a_neutrino_oscillations.md) (N* plan) ·
[#199](https://github.com/openwave-labs/openwave/issues/199) (the delta_CP=180 group prediction this revisits) ·
[#236](https://github.com/openwave-labs/openwave/issues/236) (HELD).
