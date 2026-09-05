# M5.32 term catalog (seeded by R14, 2026-09-05)

The catalog owed since PAUSE RECORD 2 ([task doc](../tasks/m5_32_task_details.md)): every term of the frozen basis B+ of R14-A with its measured behavior on the witness families ([`m5_32_r14_a_rows.json`](../data/m5_32_r14_a_rows.json), the certified sym stencil, g = 8, delta = 0.3 unless stated). Energy orientation: a positive coefficient is the energy-positive sense (the certified action is 4 I1 + V4). Columns: the hedgehog TAIL plateau (energy per unit r on the far shells; nonzero = an L-linear divergence) on two relaxed fields; the static cost per area of the (1,2) TWIST sheet and of the (2,3) ZIGZAG sheet at width 1.5 (an entry marked artifact vanishes in the continuum, certified by slab refinement); the dressed-PAIR slope E(16) - E(8) on the R3 ansatz (positive = attraction) and E(24) - E(10) on the relaxed lambda = 0 pairs (g = 32 fields); the like-charge COULOMB read E(12) - E(24) on the certified 3x3 pairs (positive = repulsion; the certified 4 I1 gives +0.28); the omega^2 coefficient on the vacuum (2,3) clock and on the hedgehog's local clock (positive = stable inertia); the sign of the term's own UV quadratic form (PSD, NSD, indefinite over the seven wave directions of R14-A; two-derivative terms only). Classes: C0/C1 the F x F quadratics, C2 the field-dependent insertions, C4 two-derivative, C5/C6 quartic.

| Term | Class | Tail n48 | Tail n32 | Twist 1.5 | Zigzag 1.5 | Pair ansatz | Pair relaxed | Coulomb | q vac clock | q hedgehog clock | UV form |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `I1` | C0C1 | 0.0746 | 0.0769 | 0 | 0 | -18.9 | -271 | 0.0692 | 0 | 70.4 | indefinite |
| `I2` | C0C1 | 0.115 | 0.128 | 0 | 0 | -37.8 | -756 | 0.0431 | 0 | 0 | indefinite |
| `I3` | C0C1 | 0.066 | 0.0704 | 0 | 0 | -18.9 | -308 | 0.0453 | 0 | 32.5 | indefinite |
| `I4` | C0C1 | 0.0712 | 0.0919 | 0 | 0 | -28.3 | -383 | 0.0337 | 0 | 2.05 | indefinite |
| `I5` | C0C1 | 0.0625 | 0.0854 | 0 | 0 | -28.4 | -451 | 0.00988 | 0 | 0 | indefinite |
| `I6` | C0C1 | 0.136 | 0.214 | 0 | 0 | -75.7 | -1.05e+03 | -0.00357 | 0 | 0 | indefinite |
| `E1` | eps | 0 | 0 | 0 | 0 | -0.000368 | -9.03 | 0 | 0 | 0 | PSD |
| `E2` | eps | 0 | 0 | 0 | 0 | -4.44e-05 | -0.45 | 0 | 0 | 0 | PSD |
| `E3` | eps | 0 | 0 | 0 | 0 | 4.93e-05 | -2.39 | 0 | 0 | 0 | PSD |
| `I1_h` | C2 | 0.0746 | 0.0769 | 0 | 0 | -19.4 | -561 | 0.0692 | 0 | 70.4 | indefinite |
| `J1` | C2 | -0.00424 | -0.00379 | 0 | 0 | 3.49 | 1.5e+03 | -0.00348 | 0 | -13.6 | indefinite |
| `J2` | C2 | -0.0112 | -0.0113 | 0 | 0 | 0.564 | -7.41e+04 | -0.0222 | 0 | -34.4 | indefinite |
| `Pgrad` | C2 | 0 | 0 | 0 | 0 | 6.06 | 2.81 | 0 | 0 | 0 | NSD |
| `K_T` | C4 | 0 | 0 | 0 | 0 | -424 | -3.11e+03 | 0 | 0 | 0 | PSD |
| `T1` | C4 | 41.2 | 39 | 0.652 | 0.03 | 423 | 3.46e+03 | -160 | 38.9 | 1.92e+04 | indefinite |
| `T2` | C4 | 0.000119 | 0.00722 | 0 | 0 | 0 | -0.00116 | 3.21e-12 | 0 | 0 | PSD |
| `T3` | C4 | 39.1 | 39.5 | 0 | 0.015 | 291 | 1.93e+03 | -55.2 | 0 | 0 | indefinite |
| `T4` | C4 | 0.00162 | -0.158 | 0 | 0 | -2e-14 | 0.121 | 3.76e-07 | 0 | 0 | indefinite |
| `K_lambda` | C4 | 0.0805 | 0.239 | 0 | 0.015 | 0 | 0.689 | -12 | 0 | 0 | PSD |
| `R_etaMeta` | C4 | -17.9 | -19.2 | 0 | 0 | 649 | 5.5e+03 | -5.54 | 0 | 0 | indefinite |
| `R_hcov` | C4 | -19.3 | -21.9 | 0 | 0 | 81.3 | 173 | -8.39 | 0 | 0 | indefinite |
| `K_P_h` | C4 | 4.84e+03 | 3.63e+03 | 0.725 (artifact) | 36.1 | -242 | 3.53e+07 | -7.27e+04 | 4.2e+04 | 2.05e+07 | PSD (exact; the stored form carries a 1e-4 finite-amplitude artifact, R14-A audit) |
| `Q_I1sq` | C5 | 5.8e-06 | 3.75e-06 | 0 | 0 | -3.56 | -631 | 5.77e-05 | 0 | 0.0112 |  |
| `C6a` | C6 | 0.575 | 0.585 | 0.284 | 0.0006 | -147 | -9.09e+03 | -126 | 0 | 292 |  |
| `C6b` | C6 | 0.32 | 0.316 | 0.284 | 0.0006 | -71 | -7.79e+03 | -112 | 0 | 17.5 |  |

Definitions (hashes with their owners: `m5_32_lagrangian.py`, `m5_32_terms_ext.py`, `m5_32_r7_a_kt_form.py`, `m5_32_r8_a_quartics.py`, `m5_32_r14_terms.py`, `m5_32_r14_a_lp.py`):

| Term | Definition |
| --- | --- |
| `I1` | I1 = sum_{mu<nu} eta^mu eta^nu <F_mu nu, F_mu nu>_eta, <F,G>_eta = tr(eta F eta G^T), F_mu nu = A_mu eta A_nu - A_nu eta A_mu; = (1/2) F_{mu nu a b} F^{mu nu a b} |
| `I2` | I2 = F_{mu nu a b} F^{a b mu nu} (derivative pair of one F contracted with the internal pair of the other; mixed pairs use delta) |
| `I3` | I3 = F_{mu nu a b} F^{mu a nu b} (eta on the mu-mu and b-b pairs, delta on the mixed nu-a pairs) |
| `I4` | I4 = R_{nu a} R^{nu a} = sum eta_nu eta_a R[nu,a]^2, R[nu,a] = sum_mu F[mu nu a mu] |
| `I5` | I5 = R_{nu a} R^{a nu} = sum_{nu a} R[nu,a] R[a,nu], R[nu,a] = sum_mu F[mu nu a mu] |
| `I6` | I6 = R^2, R = sum_nu R[nu,nu] = sum_{mu nu} F[mu nu nu mu] |
| `E1` | E1 = F[pqrs] F[xyxy] eps[pqrs] (eps[F1.d0,F1.d1,F1.i2,F1.i3] (F2.d0-F2.i2) (F2.d1-F2.i3); slots 0,1 derivative, 2,3 internal; one Levi-Civita symbol; the R0 rule extended: eps-d eta, eps-i delta) |
| `E2` | E2 = F[pqrx] F[syxy] eps[pqrs] (eps[F1.d0,F1.d1,F1.i2,F2.d0] (F1.i3-F2.i2) (F2.d1-F2.i3); slots 0,1 derivative, 2,3 internal; one Levi-Civita symbol; the R0 rule extended: eps-d eta, eps-i delta) |
| `E3` | E3 = F[pqrx] F[xysy] eps[pqrs] (eps[F1.d0,F1.d1,F1.i2,F2.i2] (F1.i3-F2.d0) (F2.d1-F2.i3); slots 0,1 derivative, 2,3 internal; one Levi-Civita symbol; the R0 rule extended: eps-d eta, eps-i delta) |
| `I1_h` | I1_h = sum_{mu<nu} eta^mu eta^nu tr(h_cov F h_cov F^T), h_cov = eta + 2 (eta u)(eta u)^T, u = timelike unit eigenvector of M eta (contravariant, u^T eta u = -1); = I1_frob in the vacuum eigenframe |
| `J1` | J1 = sum_{mu<nu} eta^mu eta^nu tr(F eta M eta F eta M eta) (contravariant F, M joined by eta; tr((FM)^2)-type) |
| `J2` | J2 = sum_{mu<nu} eta^mu eta^nu tr(F eta F eta M eta M eta) (tr(F^2 M^2)-type) |
| `Pgrad` | Pgrad = sum_mu eta^{mu mu} q(d_mu P_t, d_mu P_t), P_t = u u^T eta, q(X,Y) = sum_ab eta_a eta_b X[a,b] Y[a,b], d_mu u by first-order eigenvector perturbation of M eta along A_mu (the search knob -kappa is the coefficient) |
| `K_T` | K_T (R7): (1/2) sum_mu eta^mumu [tr(h A_mu h A_mu) - tr(eta A_mu eta A_mu)], h = h_cov |
| `T1` | eta^{mu nu} tr(A_mu eta A_nu eta) |
| `T2` | eta^{mu nu} tr(A_mu eta) tr(A_nu eta) |
| `T3` | div_b eta^{bd} div_d, div^b = sum_mu (A_mu)^{mu b} |
| `T4` | sum_{mu nu} (A_mu)^{mu nu} tr(A_nu eta) |
| `K_lambda` | K_lambda: E = (1/2) sum_a [sum_i (d_i lambda_a)^2 + omega^2 (d_t lambda_a)^2], lambda_a eigenvalues of N = M eta, d_mu lambda_a = (v_a^T eta A_mu eta v_a)/(v_a^T eta v_a) |
| `R_etaMeta` | R_G with G = etaMeta: sum_mu nu G_cd [(A_mu)^(nu c)(A_nu)^(mu d) - (A_mu)^(mu c)(A_nu)^(nu d)], mixed pairs by delta, G covariant (0,2); E-density = c_R R_G |
| `R_hcov` | R_G with G = hcov: sum_mu nu G_cd [(A_mu)^(nu c)(A_nu)^(mu d) - (A_mu)^(mu c)(A_nu)^(nu d)], mixed pairs by delta, G covariant (0,2); E-density = c_R R_G |
| `K_P_h` | K_P^h (H-adjoint form, the author's 2026-09-05 correction): E = (1/2)[sum_i tr(Om_i H Om_i^T H^-1) + omega^2 tr(Om_0 H Om_0^T H^-1)], Om_mu = P A_mu eta P, H = eta + 2 (eta u)(eta u)^T; = the Frobenius norm of the projected jet in the eta-orthonormal eigenbasis (PSD everywhere) |
| `Q_I1sq` | (F.F)^2 = I1^2 (class C5) |
| `C6a` | [sum_mu eta^mumu tr(d_mu M eta d^mu M eta)]^2 (C6) |
| `C6b` | sum_{mu nu} eta^mumu eta^nunu [tr(d_mu M eta d_nu M eta)]^2 (C6) |

Readings that the catalog carries (from the rungs, all audited): every F-built quadratic (C0/C1/C2) has a finite tail and a zero sheet cost (planar flatness, R13-W S3); the two-derivative entrants K_P^h, T1, T3 charge the hedgehog tail linearly in L (the hedgehog's (2,3)-frame connection and its divergence are charged); K_lambda and the R_G's are core-local or cancel on the tail only field by field; only T1..T4 and the C6 quartics pay the twist sheet, only K_lambda, K_P^h, T's and quartics pay the zigzag sheet; the certified 4 I1's omega^2 coefficient is negative on the hedgehog boost tangents (R14-A audit), which K_T repairs; no combination of the two-derivative class satisfies every row at once (R14-A, exact certificate), and the full basis has no witness below coefficient norm 100.
