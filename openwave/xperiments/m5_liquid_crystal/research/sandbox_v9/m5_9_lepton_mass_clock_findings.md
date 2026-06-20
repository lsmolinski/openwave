# M5.9 lepton mass law + clock-scaling findings (#200, #220)

Two issues that the contributor reframes (onspotgithub, 2026-06-18) collapse into one biaxial-hedgehog study:

- **#200**: do the charged leptons follow a mass law from the biaxial hierarchy?
- **#220**: is the geometric-mean absolute-scale calibration a genuine principle, or fortunate at the electron only?

**Reading order (2026-06-20).** The study ran in two passes:

| Pass | Scripts | Field | Method | Status |
| --- | --- | --- | --- | --- |
| Serious 4x4 build (§ below, FIRST READ) | `m5_9_4_engine_lepton.py`, `m5_9_5_higgs_clock_spectrum.py` | production 4x4 `diag(g,1,delta,0)` (boost axis live) | signed (Minkowski) energy minimization | current state of the art |
| Toy 3x3 precursor (lower sections) | `m5_9_1`, `m5_9_2`, `m5_9_3` | 3x3 `diag(1,delta,0)`, fixed Faber ansatz | evaluate / scan fixed profile | superseded on field dimension + the negative terms; the confiner DIAGNOSIS it reached stays valid |

The toy correctly DIAGNOSED the missing ingredient (a core-volume confiner) but was too simplified to compute masses: a 3x3 field with no `g` (gravity/boost) axis and no clock cannot see the NEGATIVE energy contributions, and a fixed Faber profile is not a minimization. A second-pass review on the models-of-particles list (2026-06-20) flagged exactly this. The serious build below runs on the production 4x4 Taichi engine, adds the signed-energy ledger that makes the negative gravity + oscillation contributions explicit, minimizes (Faber as the starting point only), and documents every parameter.

## M5.9.4-5: the serious 4x4 build (FIRST READ)

### Why (the second-pass review, 2026-06-20)

| Point raised | What the toy did wrong | What the serious build does |
| --- | --- | --- |
| "require full 4x4 field allowing negative mass/energy contributions from gravity and oscillations" | 3x3, no `g` axis, no clock => only POSITIVE gradient + potential | production 4x4 `diag(g,1,delta,0)`, boost axis live + the de Broglie clock; a signed (Minkowski) ledger exposes the negative terms |
| "Faber's ansatz is only the starting point of energy minimization" | evaluated/scanned the FIXED profile `s(r)=r/sqrt(r^2+r0^2)` | gradient-flow minimization (overdamped, `M_prev<-M` each step) relaxes the seed |
| "the most difficult is the Higgs-like core regularization, details to be found" | confiner only as a 3x3 toy | confiner tested on the real 4x4 signed energy (result: open, see below) |
| "I have no idea what parameters are used; the files are much too simple" | numpy toys, `delta=0.3` hardcoded, no `g` | the production `TOPOLOGY_SEED` block, every parameter dumped to `m5_9_4_results.json` |

### The signed-energy ledger (the load-bearing new physics)

The engine's DISPLAY energy (`compute_energyH_density_M`) is Euclidean (Frobenius norm, always positive), so it cannot represent a negative gravity term. The engine's DYNAMICS are Minkowski-signed (`signed_dot4`, `eta = diag(1,1,1,-1)`): the `(alpha,3)`/`(3,alpha)` time-axis components enter with a MINUS sign. `m5_9_4` adds a signed-energy kernel mirroring the display one but using `signed_dot4` for BOTH the curvature and the kinetic term. Per configuration (grid 47^3, all parameters in the JSON):

| config | H_euclid | H_signed (physical) | boost_GEM (neg. gravity) | clock_neg (neg. oscill.) |
| --- | --- | --- | --- | --- |
| undressed seed (`b*=0`) | 4.68e-05 | 4.68e-05 | 0.00 | 0.00 |
| dressed seed (`b*=0.13`) | 7.36e-05 | 4.33e-05 | 3.03e-05 | ~0 (static) |
| minimized (gradient flow) | 7.31e-05 | 4.35e-05 | 2.96e-05 | 0.00 |
| clock-active (constrained) | - | 9.89e-05 | - | kin_signed 1.83e-05, clock_neg 3.56e-06 |

- `boost_GEM = curv_euclid - curv_signed` is the negative GRAVITY contribution: the boost puts weight on the time axis, and under `eta=diag(1,1,1,-1)` those components LOWER the energy. It is `3.03e-05` (~40% of the curvature) for the dressed soliton and EXACTLY 0 for the undressed control. This is precisely the term the 3x3 structurally could not see.
- `clock_negative = kin_euclid - kin_signed` is the negative OSCILLATION contribution: clock motion along the time axis lowers the signed kinetic energy. It is real but small (`3.56e-06`) at the production dressing `b*=0.13`, because the production clock sits mostly in the spatial biaxial `(1,2)` plane rather than the time axis. Growing it requires stronger boost-clock coupling (a spectrum knob).
- Gradient-flow minimization moves `H_signed` < 1% from the seed => the production dressed seed already sits near a local energy minimum; the signed-energy descent that would DEEPEN the boost dip is the regularization-limited direction (the open Higgs problem below).
- Honest calibration: the rest mass is anchored on the STATIC minimized signed energy (robust). The constrained clock run pumps curvature (the kick excites spatial modes), so its `H_signed` (9.89e-05) is not a pure ground-state clock and is reported as the measured internal DOF, not the calibration anchor.
- Resolution-robust: re-run at 63^3 (`LEPTON_TARGET_VOXELS=262144`, `m5_9_4_results_64.json`) gives `boost_GEM = 9.32e-05` dressed vs EXACTLY 0.00 undressed, ~half the Euclidean curvature at both grids. The negative gravity dip is a structural feature of the Minkowski signature acting on the seeded boost, not a discretization artifact. The script reads `LEPTON_TARGET_VOXELS` / `LEPTON_FLOW_STEPS` / `LEPTON_CLOCK_STEPS` / `LEPTON_OUT_SUFFIX` from the environment so any run is exactly reproducible.

### The Higgs / core-volume confiner on the real engine (sub-task 4)

`m5_9_5` part A scans `r0` and measures the signed energy + the confiner integral:

| Quantity | Result | Reading |
| --- | --- | --- |
| `conf_int = integral (1-s^2)^3` | `~ r0^2.99` (R2=1.000) | the fixed-density confiner `~ r0^3` EXACTLY (matches the toy) |
| `curv_signed(r0)` | `~ r0^-0.17` (nearly flat) | the signed curvature is nearly scale-independent on the production seed |
| `V(r0)` (production LdG) | grows `~ r0^2.4` | the fixed-coefficient LdG potential is V-dominated at large `r0` |
| total static `H_signed(r0)` | INCREASES with `r0` | minimum sits at the small-`r0` grid floor (collapse), not an interior scale |

**Verdict (sub-task 4).** The confiner `~r0^3` is confirmed, but it grows with `r0` (penalizes large cores) so it does not arrest the small-`r0` collapse the production V drives; only a marginal coefficient puts the minimum interior. The production LdG potential does NOT robustly select a discrete finite core scale. This pins down precisely the open "Higgs-like regularization, details to be found": the needed `V` is one that selects a scale, and the present one does not.

### #220 BONUS: the dynamical clock omega(r0) (sub-task 5)

`m5_9_5` part B runs the constrained clock and FFTs a field probe (inside the act mask) for the natural frequency:

- 400 steps: `omega` identical for all `r0`, pinned at FFT bin 1 (`Dw = 2pi/(n*dt)`) => under-resolved (clock period exceeds the window).
- 1200 steps: frequencies begin to separate, but 3 of 4 runs blow up and survivors still pin near bin 1; the clock period exceeds even the longer CFL-stable window.
- Deeper (from part A): `H_signed` increases with `r0` => `E*r0 != const` on the production engine => the fixed-coefficient LdG `V` is NOT scale-covariant.

**Verdict (#220).** The direct dynamical `omega(r0)` readout is blocked by (clock period >> CFL-stable window) + (constrained integrator instability at long times). #220's residual stays OPEN dynamically. The analytic scale-covariance PRINCIPLE (M5.9.2, `E*r0 = const`, CV 0.00%) stands, and this build REFINES it: that covariance is a property of the scale-covariant Faber family; the production engine deliberately breaks it (its fixed-coefficient `V` selects a scale), which is the desired behavior for a definite mass but means `omega ~ 1/r0` is a statement about the covariant family, not the production dynamics as-tuned.

### #200 spectrum: hierarchy origin stays input (sub-task 5)

`m5_9_5` part C scans the biaxial `delta`:

| delta | 0.15 | 0.30 | 0.50 | 0.70 |
| --- | --- | --- | --- | --- |
| `H_signed` (norm @0.30) | 7.32 | 1.00 | 24.99 | 158.74 |

`H_signed` is non-monotonic in `delta` (minimum at the production `delta=0.30`), spanning 158x over `delta in [0.15,0.70]`; the tau needs 3477x. So `delta` is not the lepton-mass knob, and the discrete hierarchy does not emerge from the biaxial structure on the real engine either. Consistent with the toy's `E~Lambda^3` verdict: the negative contributions are now correctly included, but the hierarchy ORIGIN (why three, why 206.8 / 3477) stays Yukawa input.

### Net for the serious build

The 4x4 build delivers what the second-pass review asked for: the full field with the negative gravity + oscillation contributions made explicit and measured, real energy minimization (Faber as the start only), and fully documented parameters. It HONESTLY maps the frontier: (i) the boost-GEM negative gravity term is real and dominant; (ii) the production Higgs/confiner does not robustly select a scale (open); (iii) the eigenvalue hierarchy origin stays input; (iv) the dynamical clock readout is stability-limited. A convincing documented framework plus a precise open-problem map, which is the honest deliverable (matching mu/tau in one pass was never guaranteed). Artifacts: `m5_9_4_results.json` (17K), `m5_9_4_ledger.png` (70K), `m5_9_5_results.json` (2.1K), `m5_9_5_results.png` (59K); all kept (< 1 MB).

---

## Toy 3x3 precursor (m5_9_1 / m5_9_2 / m5_9_3, superseded on field dimension + the negative terms)

Scripts (self-contained, numpy, headless): [`m5_9_1_lepton_mass_law.py`](m5_9_1_lepton_mass_law.py), [`m5_9_2_clock_scaling.py`](m5_9_2_clock_scaling.py). Both build on the M5.6.3b Faber-on-M machinery (`M = O D(s(r)) O^T`, curvature energy `integral ||[M_mu, M_nu]||^2`, Faber potential `V = (1-s^2)^3/r0^4`).

## #200 lepton mass law: NOT reproduced by the current functional

The contributor proposes the three leptons as ground states of three eigen-axes, with `E ~ Lambda^(3/2)` from a gradient-vs-regularization balance and the hierarchy `Lambda_tau:Lambda_mu:Lambda_e ~ 229:35:1` as input (`Lambda ~ m^(2/3)`). M5.9.1 introduces the eigenvalue amplitude `Lambda` on the order parameter (`D_full = Lambda * diag(1, delta, 0)`) and measures the V-on energy.

| Test | Result | Reading |
| --- | --- | --- |
| Baseline (`Lambda=1`) | `E*r0 = 18.74` const, CV 0.00% | Faber `E ~ 1/r0` mass-pinning reproduced |
| Amplitude sweep at fixed `r0` | `E_curv ~ Lambda^4.00` (R2=1.000), `E_V ~ Lambda^0` | the M5 curvature is **quartic** in the order-parameter amplitude, not the quadratic `~Lambda^2/r0` the dimensional argument assumes |
| Lepton hierarchy (`Lambda = m^(2/3)`) | `E_mu/E_e = 1.3e6` (measured 206.8), `E_tau/E_e = 2.4e9` (measured 3477) | amplitude-`Lambda` does not reproduce the mass ratios |
| `r0`-relaxation at fixed `Lambda` | `E(r0)` monotonic, **no interior minimum** | `r0` is a free scale-setting modulus, not a balance |

**Verdict.** The `E ~ Lambda^(3/2)` mass law does not emerge from the current M5 Faber functional. Two structural reasons:

1. The curvature energy `integral ||[M_mu, M_nu]||^2` is **quartic** in the order-parameter amplitude, so amplitude-scaling gives `E ~ Lambda^4`, not the quadratic `Lambda^2/r0` the contributor's dimensional balance assumed.
2. The Faber potential integrates to `E_V ~ 1/r0` (scale-balanced), so the total energy is `~1/r0` with `r0` a **free modulus**. There is no gradient-vs-regularization minimization that selects `r0(Lambda)` and yields the `Lambda^(3/2)` law.

This is not a falsification of the three-eigen-axes idea. It is a concrete finding on what the functional needs: the **mass-selection mechanism is genuinely missing**, exactly the crux the issue named ("the Higgs-like core regularization, details open"). To produce the mass law the regularization must be a **fixed-density confiner that grows with core volume** (a `V0 * r0^3` term, the Higgs-like vacuum), so the gradient-vs-confiner balance fixes `r0(Lambda)` and selects the discrete spectrum. The current Faber `V` pins the electron (given `r0`) but does not select the lepton scales.

### The confiner build (M5.9.3): scale-selection closed, mass law E ~ Lambda^3, hierarchy origin still open

[`m5_9_3_confiner.py`](m5_9_3_confiner.py) adds the identified core-volume confiner: `E_conf = B * integral (1-s^2)^3 d^3x ~ B * r0^3`, a fixed false-vacuum energy density `B` over the melted-core volume (the Faber potential is the same integrand times `1/r0^4`, hence scale-covariant; the confiner uses a FIXED `B`, breaking scale-covariance). It then minimizes `E(r0)` per eigenvalue amplitude `Lambda`.

| Result | Reading |
| --- | --- |
| Confiner OFF: `E(r0)` monotonic (no minimum) | the free-modulus gap (M5.9.1) |
| Confiner ON: interior minimum at `r0* = 1.0` | **the confiner SELECTS the core scale, the structural gap is CLOSED** |
| `E* ~ Lambda^3.00` (Skyrme curvature + confiner, R2=1.000), `r0* ~ Lambda^1.0` | the mass law is `E ~ Lambda^3`, NOT the contributor's assumed `Lambda^(3/2)` |
| Frank gradient + confiner | `E* ~ Lambda^2` but no interior minimum (both terms grow with `r0`, collapse): the Skyrme curvature is the stabilizing gradient |

**What the confiner achieves and what it does not.**

1. **Scale-selection: closed.** The confiner fixes the free modulus, an interior energy minimum (a selected core size) now exists. This is the structural piece the issue's crux named.
2. **The mass law is `E ~ Lambda^3`** (the natural M5 quartic Skyrme curvature balanced against the `r0^3` confiner), which **corrects** the contributor's assumed `Lambda^(3/2)`. To reproduce `E ~ m` the eigenvalue hierarchy must then be `Lambda ~ m^(1/3)` (`Lambda_mu/e = 5.9`, `Lambda_tau/e = 15.1`), not the assumed `35 : 229`.
3. **Reproducing the masses is near-tautological.** Once the law (`E ~ Lambda^p`) and the eigenvalue input (`Lambda ~ m^(1/p)`) are set, `E ~ m` holds by construction, so `m_tau/m_mu = 16.82` is recovered exactly under any consistent convention (the contributor's 99.5% and the `Lambda^3` route both give it). This is a consistency check, not a parameter-free prediction.
4. **The hierarchy origin stays open.** The three-axes picture gives the **three-ness** geometrically (three eigendirections in 3D), but the eigenvalue **values** (`1 : 5.9 : 15.1`) remain Yukawa-like input. The current biaxial `D = diag(1, 0.3, 0)` does not carry those ratios. Deriving the eigenvalue hierarchy from the biaxial geometry is the irreducible remaining piece.

**Net for #200:** the confiner closes the scale-selection gap and pins the mass law (`E ~ Lambda^3`); M5 can now **accommodate** the lepton masses (selection + a definite law), but does not yet **predict** them parameter-free, the eigenvalue hierarchy origin is the open frontier.

## #220 geo-mean calibration: a PRINCIPLE by scale-covariance (dynamical confirmation pending)

The contributor notes the geo-mean law `E*r0 = alpha*(pi/4)*hbar*c` is an analytic identity (mass cancels: e/mu/tau all give 1.1309 MeV*fm), so the law itself is not the test. What is tested is whether the **sim clock** scales with mass, `omega1(L) ~ m_L`, so the geo-mean recovers each lepton's Zitterbewegung rather than returning the electron's for all.

M5.9.2 establishes the verdict quantitatively, grounded in the M5.9.1 covariance result:

| Step | Result |
| --- | --- |
| The geo-mean law | analytic identity `E*r0 = 1.1309 MeV*fm` for e/mu/tau (mass-independent by construction) |
| V-on scale-covariance (M5.9.1) | `E*r0 = const` at **CV 0.00%** across `r0`: `r0` is the only length, so any characteristic frequency obeys `omega ~ 1/r0 ~ E ~ m` |
| Geo-mean recovery, COVARIANT clock | recovery ratio `omega_gm/omega_ZBW = 1.000` for every lepton (CV 0.00%) => **PRINCIPLE** |
| Geo-mean recovery, RIGID clock (the N-6a V=0 artifact) | recovers `e` exactly but misses `mu` by 207x and `tau` by 3477x => electron-only |

**Verdict.** The geo-mean calibration is a genuine **principle**, not an electron-only coincidence. With V on the Faber soliton is exactly scale-covariant (`E*r0 = const`, CV 0.00%), so the clock scale obeys `omega ~ 1/r0 ~ m` and the geo-mean recovers every lepton's ZBW. The N-6a frequency rigidity (`omega ~ H^0.033`) was a **V=0 scale-free artifact**: the dressing knob `R_W` moves the cloud, not the scale-free core clock, exactly as N-6a's own reframe states ("the true ZBW mass-family test must vary the CORE, unavailable in the V=0 stack").

**Honest residual.** This closes #220 at the level of the analytic identity plus exact static scale-covariance, not a direct dynamical clock measurement. A direct `omega1(r0)` readout with V on requires the **V-on dynamical stack** (the M5.9/NG track, flagged unbuilt by N-6a). The scale-covariance argument is strong (an exact symmetry of a relativistic field theory implies `omega ~ 1/r0`), but the dynamical confirmation is the natural next build, and it would also be the run that, with a confiner added (the #200 gap), tests mass selection and clock-scaling together.

## Net

- **#200**: the core-volume confiner was built (M5.9.3). It **closes the scale-selection gap** (an interior energy minimum now exists) and pins the mass law at **E ~ Lambda^3** (correcting the assumed `Lambda^(3/2)`). M5 can now accommodate the lepton masses (selection + a definite law), but reproducing them is near-tautological given the eigenvalue input; the **origin of the eigenvalue hierarchy** (the Yukawa values) is the irreducible open piece. Partially closed.
- **#220**: the geo-mean absolute-scale calibration is a **principle** (precise, self-consistent across the mass family by scale-covariance), not electron-only. The remaining refinement is the direct dynamical `omega1(r0)` confirmation on the V-on stack.

Artifacts: `m5_9_1_results.npz`, `m5_9_2_results.npz` (both < 4 KB, kept).
