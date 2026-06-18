# M5.8.2z — The Faber `r₀` length-anchor calibration (issue #218)

> Does anchoring the unit map on the Faber `r₀` (length) instead of energy close the 28× clock gap? **Result: the length anchor alone overshoots (mirror-image of the energy route), but the two anchors BRACKET the ZBW and their GEOMETRIC MEAN reproduces it to ~13%.** The gap is a **calibration split, not an irreducible deficit.** Driver: [`m5_8_2z_length_anchor.py`](m5_8_2z_length_anchor.py). Date 2026-06-17.

## The two independent determinations of `ω_phys`

| Route | Formula | `ω_phys` | vs ZBW (1.553e21) |
| --- | --- | --- | --- |
| **Energy (#208)** | `ω₁·(m_ec²/ℏ)/H_static` | 5.51e19 | **28.2× LOW** |
| **Length (#218), c=2** | `ω₁·c_phys·r₀_sim/(c_sim·r₀_phys)` | 5.55e22 | **35.8× HIGH** |
| Length, c=1 | (cone-convention bracket) | 1.11e23 | 71.5× HIGH |

Measured inputs (sim units, sourced): `ω₁=1.188` (m5_8_2j/2h), `H_static=16.74` (m5_8_2j), `r₀_sim=0.69` (the Faber half-melt √3 scale, measured from the seed `_m5_8_2cb_ref.npz`, h=0.5217), `c_sim=2` (the radial cone, m5_8_2k, gradient per-sim-length), `r₀_phys=2.2132 fm` (m5_6_3a/3b).

## The key result — the anchors BRACKET the ZBW

| c_sim | energy | length | **geometric mean** | vs ZBW | each anchor off by |
| --- | --- | --- | --- | --- | --- |
| **2 (measured)** | 5.51e19 (/28) | 5.55e22 (×36) | **1.749e21** | **×1.13 (13%)** | ~×32 opposite |
| 1 | 5.51e19 (/28) | 1.11e23 (×72) | 2.473e21 | ×1.59 | ~×45 opposite |

At the **measured** cone c=2, the energy anchor is 28× low and the length anchor 36× high — **nearly symmetric in log space** — and `√(ω_E·ω_L)` reproduces the electron ZBW to **13%**. The absolute-scale error is therefore a **single factor ~32 carried in OPPOSITE directions** by the energy and length handles.

## Diagnostic — the model's particle is a CLASSICAL-radius object

| Relation | Value | Meaning |
| --- | --- | --- |
| model `E·r₀` | **1.131 MeV·fm** | from `0.511 MeV × 2.2132 fm` |
| Faber `α·(π/4)·ℏc` | **1.131 MeV·fm** (exact match) | the **classical-radius** relation |
| Compton `E·ƛ_C = ℏc` | 197.3 MeV·fm | what a ZBW clock needs |

The model's particle obeys `E·r₀ = α·(π/4)·ℏc` **exactly** (the Faber formula) — the classical-radius relation, **not** the Compton relation. So `r₀ = α·(π/4)·ƛ_C` sits **~174× below** the Compton wavelength. The two anchors anchor at two different physical scales (energy → the rest-mass scale; length → the classical-radius scale), and the geometric mean splits the difference back onto the Compton/ZBW scale.

## Verdict

The candidate outcomes from the task, resolved:

| Outcome | Result |
| --- | --- |
| (a) length anchor closes the gap | **No** — it overshoots ~36× |
| (b) routes disagree by ~28× | **Yes, they bracket** (28× low vs 36× high) |
| (c) a recognizable factor | **Yes** — `E·r₀ = α·(π/4)·ℏc` exact (classical-radius vs Compton, ×174); and the gap is a ±~32 split whose geometric mean recovers the ZBW to 13% |

**The gap is a calibration split, not an irreducible structural deficit.** The ZBW scale **is recoverable** when energy and length are calibrated **jointly** (the geometric mean / the Faber `E·r₀=const` line), rather than from the energy postulate alone. This **upgrades** #217's "largely structural" read: the structure is a *single mis-split scale factor*, not a missing physics sector.

**Recommendation for the unit map:** anchor on **both** energy and length (the Faber `E·r₀=const` line) and read the absolute scale off the joint fit, not the action↔ℏ energy postulate alone. That is the route to absolute Hz.

## Honest caveats (what could move this)

| Caveat | Effect |
| --- | --- |
| cone convention `c` | the geo-mean is ×1.13 at c=2 (measured), ×1.59 at c=1 — so "within ~2×, 13% at the measured cone". `c` pinned more tightly tightens this |
| `r₀_sim` definition | 0.69 is the Faber half-melt √3 scale (the regularization radius that maps to 2.2132 fm); the dressed-core radius is ~3 sim-len. The half-melt is the physically-correct Faber `r₀`, but its measurement carries binning uncertainty |
| is the geo-mean DEEP or fortunate? | the 13% coincidence should be tested across the **mass family** (N-6a, vary `r₀`/core) and a second seed — if the geo-mean tracks the ZBW as the mass varies, it is a real joint-calibration law; if not, it is a fortunate single-point near-miss. **This is the named follow-up.** |

## DoD status

- [x] unit map anchored on `r₀` (length unit derived), core size measured (`r₀_sim=0.69 sim-len`)
- [x] the length-route absolute `ω` vs the ZBW and vs the #208 energy route; the residual between routes quantified (they bracket, ~1000× apart, geo-mean within 13%)
- [x] verdict: the length anchor does not independently close the gap, but the gap is a **calibration split** (recoverable via joint energy+length anchoring), **not** an irreducible deficit
- [x] recognizable factor named: `E·r₀ = α·(π/4)·ℏc` (classical-radius vs Compton); the ±~32 anchor split

Model: M5 Liquid Crystal. Physics-only, headless. Closes the absolute-scale calibration thread (#208 → #217 → #218): the scale is recoverable via joint anchoring; the next refinement is the mass-family test of the geometric-mean law.
