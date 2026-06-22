# N4b checkpoint 15 , item 3: the theta12/theta23 residuals (which are free, which cost predictions)

`n4b_residual.py` (`n4b_residual_summary.json`). Key structural fact: theta23=45 + delta_CP=+-90 come from
the mu-tau MIRROR (holds at ANY tilt alpha); theta12 comes from the MAGIC crossing (a specific alpha).

| Part | Result |
| --- | --- |
| **A , theta12 residual** | TUNABLE by the tilt alpha: scanning alpha, theta12 ranges [15.5, 37.8] deg and crosses **33.68 (data) at alpha=66.7 deg** (magic alpha*=63.5), with **theta23 = 45 and delta_CP = +-90 staying EXACT** (max dev 1e-8 / 2.5e-6). So the theta12 residual costs ONLY the trimaximal PREDICTION (alpha becomes a fit), NOT delta_CP or theta23. |
| **B , theta23 residual** | matching NuFIT 43.3 needs a mu-tau breaking (chi asymmetry dchi=0.08); it moves **delta_CP to -75.5 (14.5 deg off maximal)**. So theta23's fix COSTS the maximal-CP prediction. |

## The honest dichotomy (for the article)
- **{exact TBM (theta12=35.26, theta23=45) + MAXIMAL CP}** , the clean predictive option, ~2 deg from data.
- **{theta12, theta23 fit to data}** , absorbs the residual but loses maximal CP (delta_CP ~ -75 not -90).
- The **delta_CP measurement (DUNE/HK)** decides between them. A near-maximal delta_CP (~270) favours exact
  TBM; a value pulled toward 255 would favour the slightly-broken fit. Either is presentable + falsifiable.

The ~2 deg residual is the standard TBM-vs-data gap, now characterized precisely: theta12 is free, theta23
trades against maximal CP.
