# Checkpoint 02 - sub-tasks 1-3 DONE (foundation + signed ledger + minimization)

Run 2 of m5_9_4 (grid 47^3, 46s, clock now stable). m5_9_4_results.json (17K) + ledger.png (70K).

## What is now established (HONEST)

| config | H_euclid | H_signed | boost_GEM | clock_neg |
| --- | --- | --- | --- | --- |
| seed_undressed (b*=0) | 4.68e-05 | 4.68e-05 | 0.00 | 0.00 |
| seed_dressed (b*=0.13) | 7.36e-05 | 4.33e-05 | 3.03e-05 | 2.5e-15 |
| minimized (grad flow) | 7.31e-05 | 4.35e-05 | 2.96e-05 | 0.00 |
| clock_time_avg | - | 9.89e-05 | - | (kin_s 1.83e-05, clk_neg 3.56e-06) |

## Sub-task status

- [1] ✅ production 4x4 engine stood up headless; dressed hedgehog seeded; electron rest
      energy read by the signed ledger; e_scale calibration = 1.176e+04 MeV/unit (on the
      static minimum). All parameters dumped to JSON (Duda's ask).
- [2] ✅ negative-energy instrumentation: signed curvature (boost-GEM, negative GRAVITY) +
      signed kinetic (clock_negative, negative OSCILLATION). Both measured. Boost-GEM is the
      dominant negative term (3.03e-05); clock_negative is small (3.56e-06) at b*=0.13 because
      the production clock sits in the spatial biaxial (1,2) plane, not the time axis. To grow
      it, increase boost-clock coupling (a spectrum knob).
- [3] ✅ true minimization (overdamped gradient flow, M_prev<-M each step) - relaxes the Faber
      ansatz (Duda's point); confirms the dressed seed is near a local Euclidean min (<1%).
      The signed-energy descent that DEEPENS the boost dip is regularization-limited = Duda's
      open Higgs problem (sub-task 4).

## Honest caveats to carry into the writeup
- The clock-active H_signed (9.89e-05) is HIGHER than static (4.35e-05): the constrained clock
  run pumps curvature (the kick excites spatial modes), so it is not a pure ground-state clock.
  => calibrate the rest mass on the STATIC minimized signed energy (robust); treat the clock as
  the measured internal DOF (the #220 face), not the calibration anchor.
- "~21% clock reduction" was NOT reproduced and was never asserted; the measured clock_negative
  is 3.56e-06 (small). Report the measured number, not a target.

## Next (sub-tasks 4 + 5, script m5_9_5)
- [4] Higgs/confiner search on the REAL engine: scan r0, measure signed E(r0) with/without a
      fixed-density confiner; does an interior minimum (scale-selection) appear, now with the
      4x4 boost? (the toy m5_9_3 showed yes for 3x3).
- [5a] #220 BONUS: measure the dynamical clock rate omega(r0) across a resolvable r0 range,
       confirm omega ~ 1/r0 (mass-scale law). This closes #220's residual (the V-on dynamical
       readout now exists).
- [5b] #200 spectrum: at the selected scale + biaxial eigen-structure, read masses, compare to
       m_mu/m_e=206.8, m_tau/m_e=3477. Honest: the discrete hierarchy ORIGIN is expected to stay
       input (the toy showed E~Lambda^3 with hierarchy as Yukawa input); report how far off.
