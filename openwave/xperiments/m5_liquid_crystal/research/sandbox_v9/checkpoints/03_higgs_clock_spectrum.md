# Checkpoint 03 - sub-tasks 4 + 5 DONE (Higgs/confiner, clock scaling, spectrum)

m5_9_5 (grid 47^3). Two runs: 400-step clock (under-resolved), 1200-step clock (separates but
blows up). m5_9_5_results.json (2.1K) + png (59K). All HONEST results below.

## (A) Higgs / confiner scale-selection [sub-task 4]

| r0_vox | curv_signed | V (LdG) | conf_int (1-s^2)^3 |
| --- | --- | --- | --- |
| 1.88 | 2.49e-05 | 1.05e-05 | 1.64e+01 |
| 2.82 | 2.26e-05 | 2.07e-05 | 5.52e+01 |
| 6.11 | 2.03e-05 | 1.77e-04 | 5.52e+02 |

- conf_int ~ r0^2.99 (R2=1.000): the fixed-density confiner ~ r0^3 EXACTLY (matches the m5_9_3 toy). ✅
- curv_signed ~ r0^-0.17 (nearly FLAT), V grows ~r0^2.4 => total static energy is V-DOMINATED and
  INCREASES with r0 => the minimum sits at the SMALL-r0 grid floor (collapse), not an interior scale.
- A confiner ~r0^3 grows with r0 (penalizes large r0) so it does NOT arrest small-r0 collapse; only
  a marginal B<1 puts argmin interior. VERDICT: the production LdG potential does NOT robustly select
  a discrete finite core scale. This CONFIRMS Duda's "Higgs-like regularization details still to be
  found" - finding the V that selects a scale is the open problem, now pinned precisely.

## (B) #220 BONUS: dynamical clock omega(r0) - NOT cleanly closed (honest)

- 400 steps: omega=0.19638 IDENTICAL for all r0 = exactly FFT bin 1 (Dw=2pi/(400*0.08)). Under-resolved.
- 1200 steps: frequencies start to separate (0.079, 0.065...) but 3 of 4 runs BLEW UP, and survivors
  still pin near bin 1 (0.0654=2pi/(1200*0.08)). The clock period EXCEEDS even the longer stable window.
- DEEPER (from part A): H_signed INCREASES with r0 on the production engine => E*r0 != const => the
  production fixed-coefficient LdG V is NOT scale-covariant. So omega~1/r0 is a property of the
  scale-COVARIANT Faber family (#220's ANALYTIC result, CV 0.00%), which the production engine
  deliberately breaks to select a scale.
- VERDICT for #220: the dynamical readout is blocked by (clock period >> CFL-stable window) +
  (constrained integrator instability at long times). #220's residual stays OPEN dynamically; the
  analytic scale-covariance PRINCIPLE stands, and part A REFINES it (production V is non-covariant).
  This is the honest "bonus": a refinement, not a closure.

## (C) #200 spectrum [sub-task 5] - hierarchy origin stays input (honest)

| delta | H_signed | ratio (norm @0.30) |
| --- | --- | --- |
| 0.15 | 3.17e-04 | 7.32 |
| 0.30 | 4.33e-05 | 1.00 |
| 0.50 | 1.08e-03 | 24.99 |
| 0.70 | 6.87e-03 | 158.74 |

- H_signed is NON-monotonic in delta (min at the production delta=0.30, rises both ways). The delta
  knob spans 158x over delta in [0.15,0.70]; the tau needs 3477x. So delta is NOT the lepton-mass knob.
- VERDICT: varying the biaxial structure does not produce 206.8 / 3477. The discrete hierarchy ORIGIN
  stays Yukawa INPUT (consistent with the toy's E~Lambda^3 result). The real-engine upgrade adds the
  negative contributions but does NOT make the spectrum emerge. Honest, expected per the DoD.

## Net for the writeup (sub-task 6 next)
The serious build DELIVERS what Duda asked (4x4, negative gravity + oscillation terms, minimization,
documented params) and HONESTLY maps the frontier: (i) the boost-GEM negative gravity term is real and
dominant, (ii) the Higgs/confiner does not robustly select a scale on the production V (open), (iii) the
hierarchy origin stays input, (iv) the dynamical clock readout is stability-limited. This is a convincing
framework + an honest open-problem map, exactly the DoD deliverable.
