# N4b checkpoint 16 , item 6 (chiral origin) + item 4 (two delta scales) , N4b COMPLETE

## item 6 , origin + natural scale of g_chiral (`n4b_chiral_origin.py`)

| g_chiral | delta_CP | theta13 | |C| (geometric) | E(+chi)-E(-chi) |
| --- | --- | --- | --- | --- |
| 0.00 | 0.0 | 0.000 | 0.945 | 0 |
| 0.30 | -90 | 2.86 | 0.945 | 0 |
| 0.94 | -90 | 8.58 | 0.945 | 0 |
| 1.50 | -90 | 12.84 | 0.945 | 0 |

- **CP requires the chiral coupling g_chiral != 0.** g_chiral=0 -> delta_CP=0, theta13=0 (CP-conserving, the
  real N3 limit, = #199's original prediction). g_chiral!=0 -> delta_CP=+-90 + theta13. CP hinges on a chiral
  (Lifshitz) SUBSTRATE term.
- The overlap |C| = 0.945 is GEOMETRIC (independent of g_chiral; set by the loop SO(3)-orientation arrangement).
- The +chi/-chi screw configs are achiral-LdG DEGENERATE (E(+chi)=E(-chi)) -> the achiral substrate does NOT
  prefer a handedness; only the SIGN of g_chiral selects delta_CP's sign.
- theta13 = 8.56 at g_chiral* ~ 0.94 = O(1) -> theta13 ~ O(10 deg) is NATURAL if the substrate is chiral at
  its natural scale. **Open question for Duda: does the M5 LdG carry a Lifshitz/cholesteric invariant?**

## item 4 (refined) , the two delta scales (inline check)

TBM-gate angles vs delta (the biaxiality eigenvalue), delta = 1e-10 .. 0.3: **theta12=35.2644, theta23=45,
theta13=0, alpha*=46.945 deg , IDENTICAL for ALL delta.** The TBM mixing is delta-INDEPENDENT (purely
geometric). So:

- Duda's quantum-phase delta ~ 1e-10 (the Dirac/QED correction) plays NO role in the PMNS mixing.
- The mixing's energy scale is the GEOMETRY (mass spectrum 1:1.15:1.68), not delta.
- theta13 / delta_CP are set by the chiral coupling g_chiral, not by delta.
- => the original N3 "central tension" (theta13=8.5 needs delta~0.15 not 1e-10) is FULLY resolved: the mixing
  doesn't see delta at all; theta13 is chirality-driven; the 1e-10 quantum-phase delta is a separate scale.

## N4b COMPLETE
items 2 (potential robustness), 3 (residuals), 6 (chiral origin), 4 (two delta scales) all closed. Duda round
has no pending tasks. Next: N5 (article) + N6 (masses), both deferred. Consolidate into n4b_findings.md.
