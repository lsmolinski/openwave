# Mixing-integral normalization resolution (2026-06-08)

Response to Werbos/DeepSeek's 2026-06-08 question: the raw mixing integral
M_JA(0) computed on the canonical profile is 271x the DM-paper draft value, which
plugged into `g_JA_eff = M_JA / V_Yukawa` gives a mixing fraction > 1 (unphysical)
and a cross section above LZ. This file records why that is a formula artifact, not
a physics exclusion. All values reproduce from `v11_canonical_beta_profile.csv`.

## Canonical anchors

| Quantity | Value |
| --- | --- |
| hbar*c | 197.3269804 MeV*fm |
| m_e | 0.5109989461 MeV |
| R_phys = hbar*c/m_e | 386.1593 fm |
| m_J | 0.618427 MeV |
| proton charge radius | 0.84 fm |

## Amplitude-INVARIANT spatial quantities (these carry the physics)

Both are ratios in which the overall beta scale cancels exactly, so the 8x
amplitude difference (peak 0.697 vs draft 0.085) and the 271x integral difference
do NOT propagate into sigma_p.

| Quantity | Value | Note |
| --- | --- | --- |
| beta_peak | 0.69707 @ r = 687.8 fm | physical amplitude, fixed by a = m_J/(2 sqrt(g)) |
| origin slope B0 | 1.2949e-03 /fm (0.50005 natural) | l=1 regularity beta ~ B0 r |
| beta(0.84 fm) | 1.0877e-03 | = B0 * 0.84 |
| ratio beta(0.84)/beta_peak | 1.5605e-03 | amplitude cancels |
| **S_vert = ratio^2** | **2.4350e-06** | draft used 1.54e-4 (wrong peak) -> 63x too large |
| F(q_gal = 3e-4 MeV) | 0.9999985 | normalized F(0)=1; ~1 confirms spatial-only suppression |
| F(q = 3 MeV) | -6.140e-05 | strongly suppressed at LZ scales |

## Raw (NON-invariant) g_JA_eff -- the bug, for the record

`g_JA_eff = M_JA / V_Yukawa` is linear in the beta amplitude (M_JA is proportional
to beta; V_Yukawa is a fixed geometric volume), so it is NOT a valid mixing
fraction. With the physical amplitude it explodes:

| Quantity | Value |
| --- | --- |
| int beta r^2 dr | 9.1343e+08 fm^3  (nat 15.8627) |
| M_JA(0) = 4*pi*int | 1.1479e+10 fm^3  (nat 199.336) |
| V_Yukawa = (4pi/3)(hbar c/m_J)^3 | 1.3608e+08 fm^3 |
| g_JA_eff RAW | 84.35  (>>1, unphysical) |
| ratio vs draft M_JA (4.23e7) | 271x  =  8.21x amplitude  *  33.0x extent |

## Fix (Werbos side): make g_JA_eff amplitude-invariant

Option A -- Lagrangian-level: the J-A mixing is set by the fixed J^mu A_mu coupling
and m_J^2, independent of the soliton amplitude; beta enters only via F(q) and S_vert.

Option B -- normalize M_JA by the soliton's own J-norm sqrt(N_J), N_J = int beta^2 r^2 dr
= 2.4522e+08 fm^3 (scales as amp^2), instead of the external V_Yukawa, so the scale cancels.

## Open structural point (l=1)

The chaoiton is l=1 (p-wave). M_JA(0) = 4*pi*int beta r^2 dr implicitly projects
onto l=0; for a pure l=1 field the angular integral of Y_1 over the sphere vanishes,
so the true zero-momentum monopole coupling is zero -- arguably the meaning of
"neutral." If so, M_JA(0) is the wrong object and the surviving coupling is
dipole/higher (automatically suppressed), with S_vert as the spatial stand-in.

Regenerate: `python regenerate_supplement.py`
