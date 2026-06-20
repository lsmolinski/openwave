# Checkpoint 01 - signed energy ledger works; clock run needs fixing

Run 1 of m5_9_4 (grid 47^3, 76s). Result:

## WIN: the boost-GEM negative gravity contribution is real and measured

| config | H_euclid | H_signed | boost_GEM |
| --- | --- | --- | --- |
| seed_undressed (b*=0) | 4.68e-05 | 4.68e-05 | 0.00e+00 |
| seed_dressed (b*=0.13) | 7.36e-05 | 4.33e-05 | 3.03e-05 |
| minimized (grad flow) | 7.31e-05 | 4.35e-05 | 2.96e-05 |

- Undressed: signed == euclid EXACTLY (boost_GEM=0). The Minkowski signature is a no-op
  when there are no `(alpha,3)` components. Sanity check passes.
- Dressed: signed curvature (2.26e-05) < euclid (5.29e-05). The boost puts weight on the
  time axis, and under eta=diag(1,1,1,-1) those components LOWER the energy by 3.03e-05.
  **This is exactly the negative gravity contribution Duda said the 3x3 could not see.** ✅
- Gradient-flow minimization moves H_signed <1% => the production dressed seed already sits
  near a local energy minimum (Euclidean). The signed-energy descent that would DEEPEN the
  boost dip is the regularization-limited direction (Duda's open Higgs problem, sub-task 4).

## BUG: clock-activation run blew up

`H_clock_avg = 1.48e+04` vs static 4.3e-05 => 9 orders of magnitude runaway. Causes:
1. timestep too large: I used DT_SCALE 0.05; the validated headless template uses 0.007 for
   the constrained integrator. Fix: dt_rs_b = 0.007 * ... + a blow-up guard that bails.
2. deeper physics: the KINETIC term must ALSO be Minkowski-signed. Clock oscillation along
   the time axis (the (alpha,3) velocity from the boost+kick) should contribute NEGATIVELY
   (Duda's "negative contributions from oscillations"). Currently kinetic = 1/2||Mdot||^2
   (Frobenius, positive). Add kin_signed = 1/2*signed_dot4(Mdot,Mdot); the (alpha,3)/(3,3)
   velocity comps then lower the energy = the oscillation term.

Honest note on "~21%": I pre-baked a "platform measured ~21%" comparison. Do NOT assert that
without checking what 21% actually refers to. Reframe: MEASURE the clock's energy budget
(signed kinetic share, sign) and report what it is.

## Next
- Fix energy_decompose: add signed kinetic. Fix clock dt -> 0.007 + guard. Re-run.
- Then sub-task 4 (Higgs/confiner search) + sub-task 5 (the three configs -> masses).
