# N4c checkpoint 18 , the mass-ratio near-term falsifier

Reviewer (Claude) flagged: the spectrum 1:1.15:1.68 is checkable against data NOW, and looks compressed vs the
observed splitting hierarchy. Tested in `n4c_mass_ratio.py`.

## Result

| quantity | value |
| --- | --- |
| gate eigenvalues (at `alpha*`=46.945 deg) | [1690.97, 1940.99, 2844.86] = **1 : 1.148 : 1.682** |
| observed `Dm31^2/Dm21^2` (NuFIT 6.0 NO) | 2.513e-3 / 7.49e-5 = **33.55** |
| Map A (`lambda = m`): `(l3^2-l1^2)/(l2^2-l1^2)` | **5.76** -> off by **x5.8** |
| Map B (`lambda = m^2`): `(l3-l1)/(l2-l1)` | **4.62** -> off by **x7.3** |

## Verdict

The gate spectrum is QUASI-DEGENERATE (1 : 1.15 : 1.68); the data needs a splitting ratio ~34, i.e. a much more
spread spectrum. Under both natural eigenvalue->mass maps the loop splittings are ~5-7x too compressed -> a real
near-term TENSION.

**FLAG, not a falsification:** the eigenvalue->mass map is undefined here (the deferred N6 question, Duda's
mass-as-loop-length-density picture). A successful mass model (N6) must either spread the spectrum or use a
nonlinear map. Worth stating plainly before building further (the reviewer's point stands).

Artifacts: `n4c_mass_ratio.py`, `n4c_mass_ratio_summary.json`, `n4c_mass_ratio.png`. Next: honest scorecard
figure (`n4c_scorecard.py`) + write `10e_findings_N4c.md`.
