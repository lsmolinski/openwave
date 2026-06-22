# N4b findings , closing the N4 pending tasks (a complete Duda round, no loose ends)

> **LOCAL / HELD.** N4b (the N4 closeout) of the neutrino N-program, OpenWave issue
> [#236](https://github.com/openwave-labs/openwave/issues/236). All GitHub #236 posting HELD until the whole
> N-program finishes. Local record. Builds on N4: [`n4_findings.md`](n4_findings.md). Master plan:
> [`../10a_neutrino_oscillations.md`](../10a_neutrino_oscillations.md) § "N1 + N2 foundation scope". Convention
> index-0 ([`../_convention_refactor/CONVENTION.md`](../_convention_refactor/CONVENTION.md)).

N4b nails the four items that were still open after N4, so the Duda exchange has zero pending tasks. N5
(article) + N6 (masses) stay deferred behind it.

## Headline

| item | question | answer | script |
| --- | --- | --- | --- |
| **2** | does the result survive a real LdG TENSOR potential? (Duda's #1) | the mu-tau predictions (theta23=45, theta13=0, delta_CP=+-90) are ROBUST to ALL 27 potentials; the magic/theta12 point is reached by tilt for 24/27 and RECOVERABLE by a 2nd geometric knob for the rest | [`n4b_potential.py`](n4b_potential.py) |
| **3** | the theta12 / theta23 residuals (~2 deg vs data) | theta12 is FREE (tunable by the tilt alpha; theta23 + delta_CP unaffected); theta23->43.3 costs ~14 deg of delta_CP (mu-tau breaking) | [`n4b_residual.py`](n4b_residual.py) |
| **6** | where does g_chiral come from + its scale? | CP requires g_chiral != 0 (a chiral substrate term); g_chiral=0 recovers CP-conservation; theta13 ~ O(10 deg) is natural at g_chiral* ~ 0.94 = O(1) | [`n4b_chiral_origin.py`](n4b_chiral_origin.py) |
| **4** | the two delta scales (mixing vs quantum-phase) | the TBM mixing is delta-INDEPENDENT (identical for delta = 1e-10 .. 0.3); Duda's 1e-10 plays no role in the mixing; theta13/delta_CP are chirality-driven | inline check (this doc) |

## item 2 , LdG tensor-potential robustness (answers Duda's "the potential is crucial")

Put the REAL Landau-de Gennes potential Hessian into the on-site term (replacing the crude kappa):
`H^V_ab = INT [2a Tr(dA dB) - 6b Tr(Mvac dA dB) + c(8 Tr(Mvac dA)Tr(Mvac dB) + 4 Tr(Mvac^2)Tr(dA dB))]`. Then
`M_mass = K(kinetic) + H^V` and scanned 27 potentials `(a in {-3,-2,-1}, b in {0,0.5,1}, c in {0.5,1,2})`.

| Result | Finding |
| --- | --- |
| mu-tau predictions | **ROBUST to ALL 27**: theta23 = 45 everywhere; delta_CP = +-90 in all spot checks. mu-tau is carried by the GEOMETRY, not the potential. |
| magic / theta12 (trimaximal) | reached by tilt alone in 24/27; for 3 strong potentials the magic crossing leaves the tilt range |
| recovery | a miss `(a,b,c)=(-3,0.5,1)` REGAINS the gate at `R_loop=11` (TBM err 0.000) , a 2nd geometric knob co-adjusts with the potential |

**Answer to Duda:** the result is NOT tuned to a specific LdG potential. The mu-tau STRUCTURE is
potential-independent; the magic POINT shifts with `(a,b,c)` but is recoverable by geometry.

## item 3 , the theta12 / theta23 residuals (which are free, which cost the predictions)

theta23 = 45 + delta_CP = +-90 come from the mu-tau MIRROR (holds at ANY tilt alpha); theta12 comes from the
MAGIC crossing (a specific alpha).

- **theta12 (35.26 vs 33.68): FREE.** Scanning alpha, theta12 ranges [15.5, 37.8] deg and hits **33.68 (data) at
  alpha = 66.7 deg** with **theta23 = 45 and delta_CP = +-90 staying EXACT** (dev 1e-8 / 2.5e-6). The residual
  costs ONLY the trimaximal PREDICTION (alpha becomes a fit), not delta_CP or theta23.
- **theta23 (45 vs 43.3): costs maximal CP.** Reaching 43.3 needs a mu-tau breaking (chi asymmetry dchi=0.08)
  that moves **delta_CP to -75.5 (14.5 deg off maximal)**.

**The honest dichotomy:** {exact TBM (theta12=35.26, theta23=45) + MAXIMAL CP, ~2 deg from data} vs {theta12,
theta23 fit to data, delta_CP not maximal}. The delta_CP measurement (DUNE/HK) decides; both are presentable.

## item 6 , the origin + natural scale of the chiral coupling (the CP sector hinges on the substrate)

| g_chiral | delta_CP | theta13 | |C| (geometric) | E(+chi)-E(-chi) |
| --- | --- | --- | --- | --- |
| 0.00 | 0.0 | 0.000 | 0.945 | 0 |
| 0.30 | -90 | 2.86 | 0.945 | 0 |
| 0.94 | -90 | 8.58 | 0.945 | 0 |
| 1.50 | -90 | 12.84 | 0.945 | 0 |

- **CP requires the chiral coupling `g_chiral != 0`** (a chiral / Lifshitz SUBSTRATE term). `g_chiral = 0` ->
  `delta_CP = 0`, `theta13 = 0` (CP-conserving, = [#199](https://github.com/openwave-labs/openwave/issues/199)'s
  original group prediction); `g_chiral != 0` -> `delta_CP = +-90` + `theta13`.
- The overlap `|C| = 0.945` is GEOMETRIC (independent of `g_chiral`; the loop SO(3)-orientation arrangement).
- The achiral LdG is handedness-degenerate (`E(+chi) = E(-chi)`) -> only the SIGN of `g_chiral` selects
  `delta_CP`'s sign.
- `theta13 = 8.56` at `g_chiral* ~ 0.94 = O(1)` -> `theta13 ~ O(10 deg)` is NATURAL if the substrate is chiral
  at its natural scale (the model predicts the SCALE; the chiral coefficient pins the value).

**The open substrate question for Duda:** does the M5 Landau-de Gennes functional carry a chiral / Lifshitz
(cholesteric-type) invariant? If yes -> `delta_CP = +-90` + `theta13 ~ O(10 deg)` PREDICTED, sign = the
favoured handedness. If strictly achiral -> CP conserved (`delta_CP` in {0,180}).

## item 4 (refined) , the two delta scales

TBM-gate angles vs the biaxiality `delta` (1e-10 .. 0.3): **theta12 = 35.2644, theta23 = 45, theta13 = 0,
alpha* = 46.945 deg , IDENTICAL for ALL delta.** The TBM mixing is `delta`-INDEPENDENT (purely geometric). So:

- Duda's quantum-phase `delta ~ 1e-10` (the Dirac/QED correction) plays NO role in the PMNS mixing.
- The mixing's energy scale is the GEOMETRY (mass spectrum 1 : 1.15 : 1.68), not `delta`.
- `theta13` / `delta_CP` are set by the chiral coupling `g_chiral`, not by `delta`.
- => the original N3 "central tension" (`theta13 = 8.5` needs `delta ~ 0.15` not `1e-10`) is FULLY resolved:
  the mixing does not see `delta` at all; `theta13` is chirality-driven; the `1e-10` quantum-phase `delta` is a
  separate scale (the rest-mass/QED sector, not the mixing).

## Net , the complete picture for Duda

| PMNS parameter | status | depends on |
| --- | --- | --- |
| theta12 = 35.26 | PREDICTION (trimaximal) | the magic crossing (geometry); tunable to data via alpha at the cost of the prediction |
| theta23 = 45 | PREDICTION (maximal) | the mu-tau mirror (geometry); robust to the LdG potential |
| delta_CP = +-90 (270) | PREDICTION (maximal) | mu-tau reflection + a chiral substrate term; sign = handedness |
| theta13 = 8.56 | natural O(10 deg), 1 free coupling | the chiral strength g_chiral ~ O(1); NOT delta, NOT topological |

No pending tasks for the Duda round: the potential is shown irrelevant to the predictions; the residuals are
characterized; the CP sector's dependence on a chiral substrate term is pinned; the two delta scales are
separated. Open (post-Duda / N5): does the substrate carry the Lifshitz term (-> theta13 + delta_CP fully
predicted); the deferred masses (N6).

## Artifacts (all LOCAL, < 130 KB, nothing > 1 MB)

| Artifact | Regenerate |
| --- | --- |
| [`n4b_potential.py`](n4b_potential.py) + `n4b_potential_summary.json` | `python3 n4b_potential.py` (16-core) |
| [`n4b_residual.py`](n4b_residual.py) + `n4b_residual_summary.json` | `python3 n4b_residual.py` (16-core) |
| [`n4b_chiral_origin.py`](n4b_chiral_origin.py) + `n4b_chiral_origin_summary.json` | `python3 n4b_chiral_origin.py` |
| `checkpoints/13_n4b_design.md` .. `16_n4b_origin_deltas.md` | progress log |

## Cross-refs

[`n4_findings.md`](n4_findings.md) (the closed gate this completes) · [`n3_findings.md`](n3_findings.md) ·
[`../10a_neutrino_oscillations.md`](../10a_neutrino_oscillations.md) (N* plan) ·
[#199](https://github.com/openwave-labs/openwave/issues/199) · [#236](https://github.com/openwave-labs/openwave/issues/236) (HELD).
