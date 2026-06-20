# M5.9 lepton mass law + clock-scaling findings (#200, #220)

Two issues that the contributor reframes (onspotgithub, 2026-06-18) collapse into one V-on biaxial-hedgehog study:

- **#200**: do the charged leptons follow a mass law `E ~ Lambda^(3/2)` from the biaxial hierarchy?
- **#220**: is the geometric-mean absolute-scale calibration a genuine principle, or fortunate at the electron only?

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
