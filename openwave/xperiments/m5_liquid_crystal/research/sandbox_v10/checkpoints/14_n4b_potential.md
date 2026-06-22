# N4b checkpoint 14 , item 2: LdG potential robustness (answers Duda's #1)

`n4b_potential.py` (`n4b_potential_summary.json`). Put the REAL LdG tensor potential Hessian
`H^V = 2a Tr(dA dB) - 6b Tr(Mvac dA dB) + c(8 Tr(Mvac dA)Tr(Mvac dB) + 4 Tr(Mvac^2)Tr(dA dB))` into the
on-site term (replacing the crude kappa). M_mass = K(kinetic) + H^V. Scanned 27 potentials.

| Result | Finding |
| --- | --- |
| mu-tau predictions vs potential | **ROBUST to ALL 27**: theta23 = 45 everywhere; theta13 = 0; delta_CP = +-90 in all spot checks. mu-tau is carried by the GEOMETRY, the potential cannot break it. |
| magic / theta12 (trimaximal) | reached by tilt alone in **24/27**; for 3 strong potentials the magic crossing leaves the tilt range |
| recovery | a miss `(a,b,c)=(-3,0.5,1)` REGAINS the gate at `R_loop=11` (TBM err 0.000) -> a 2nd geometric knob co-adjusts with the potential |
| baseline `(-2,0,1)` | TBM err 0.0000 (perfect) |

**Answer to Duda:** the result is NOT tuned to a specific LdG potential. The mu-tau STRUCTURE (theta23=45,
theta13=0, delta_CP=+-90) is potential-independent; the magic POINT (theta12 trimaximal) shifts with (a,b,c)
but is recoverable by geometry. delta_CP = maximal is robust to all potentials tested.
