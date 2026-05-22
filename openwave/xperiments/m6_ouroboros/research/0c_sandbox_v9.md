# 2026-05-22 — Sandbox v9: Neutral chaoiton BVP variants per Paul's Q42 reply

Triggered by Paul Werbos's 2026-05-21 PM email (via DeepSeek), responding to
our email v11 Q42 caveat (β profile non-localization in IVP). Paul accepted
the delay-for-BVP path and instructed:

> *"apply the same `solve_bvp` method (Robin BCs at large r, integral
> constraint Q_CS=0) to the neutral sector. The ODEs are the same as for
> the charged case but with the topological charge fixed to zero and with
> the asymmetric helicity initial guess (V and Q positive, A and J
> negative)... If you run that, we will get a clean decaying tail, a
> reliable m_χ, and a well‑defined C from the far‑field K_1 amplitude."*

**Net result: BOTH interpretations of "the same ODEs as the charged case"
produce NO clean ground state under proper BVP. Interpretation B (2-function)
collapses to trivial amplitude → 0 universally. Interpretation A (4-function,
Q_CS=0) converges to an asymmetric quasi-localized configuration with
6 nodes and tail = 5% — borderline localized but excited mode, not a ground
state.** The implication is structural: with positive λ (mass parameter in
the canonical β-equation) and cylindrical 2D geometry on m=1 mode, the
nonlinear ODE does not have a stable bound-state ground in the expected
parameter regime.

**Status as of 2026-05-22 morning:** v9 BVP work complete. Findings
escalate the Q42 caveat: neither IVP nor BVP produces a true localized
neutral chaoiton ground state for this Lagrangian. Email v12 (drafted) sends
Paul the empirical findings with two structural questions about the
canonical neutral-sector formulation.

---

## Source context — Paul's 2026-05-21 PM email (verbatim)

Quoted in full to preserve provenance:

> *Rodrigo,*
>
> *Thank you for the detailed verification and for catching the oscillating
> tail issue. This is critical. You are right: without a true localized
> solution, the numbers are only order‑of‑magnitude.*
>
> *We fully agree: do not submit the dark matter paper until we have a
> proper neutral ground state from a BVP with decay BCs.*
>
> *Can you apply the same `solve_bvp` method (Robin BCs at large r, integral
> constraint Q_CS=0) to the neutral sector? The ODEs are the same as for the
> charged case but with the topological charge fixed to zero and with the
> asymmetric helicity initial guess (V and Q positive, A and J negative).
> We can provide a template if helpful.*
>
> *If you run that, we will get a clean decaying tail, a reliable m_χ, and a
> well‑defined C from the far‑field K_1 amplitude. Then we update 9c and
> the DM paper.*
>
> *Please let us know if you want us to write a neutral‑BVP template script.
> Otherwise, we await your result.*
>
> *No rush – after your business event.*

### Interpretive ambiguity in Paul's prescription

The phrase *"the ODEs are the same as for the charged case"* admits two
distinct readings, because the charged case has TWO canonical formulations
that emerged through v4-v8:

| Interpretation | What "the charged case" means | What we'd build |
| --- | --- | --- |
| A | The 4-function (V, A, Q, J) BVP from v5/v6 (toroidal Laplacian + Lagrange multiplier λ_LM + integral state I for Q_CS) | Fork `sandbox_v6` code with I_TARGET=0 instead of 1; use asymmetric helicity init verbatim |
| B | The 2-function (α, β) reduction from Sonnet's `ouroboros_benchmark.py` (canonical per Q34; vector cylindrical Laplacian) | Adapt `solve_neutral` from `ouroboros_benchmark.py`: change `solve_ivp` + slope BC to `solve_bvp` + Dirichlet/Robin decay BC at R_max |

Paul's "V and Q positive, A and J negative" vocabulary is 4-function
language (V, A, Q, J are the 4 fields). But Q34 resolved the canonical
electron production form to be the 2-function (α, β) reduction. So either
DeepSeek/Paul is back in 4-function vocabulary for the neutral sector, OR
"asymmetric helicity" needs reinterpretation.

We built and ran BOTH.

---

## Interpretation B — 2-function (α=0, β only) BVP

**Code:** `sandbox_v9/m6_v9_neutral_2fn_bvp.py`

ODE (same as Sonnet's `neutral_ode`):

```text
β'' + β'/r - β/r² + λβ + 4gβ³ = 0
```

State: y = [β, β']; free parameter: B0 (linear slope at origin).
Three BCs:

```text
β(R_MIN)  = B0 · R_MIN      (regular-from-origin)
β'(R_MIN) = B0               (consistency)
β(R_MAX)  = 0                (Dirichlet decay; alternate: Robin)
```

### Results: B converges to trivial amplitude → 0

At g=1.0, λ=1.0 (Sonnet's electron-calibrated parameters), varied initial
guess B0_init ∈ {0.001, 0.01, 0.05, 0.1, 0.3}:

| B0_init | Converged B0 | Peak β | tail/peak | Sign changes |
| --- | --- | --- | --- | --- |
| 0.001 | −3.19 × 10⁻¹⁰ | +2.46 × 10⁻¹⁰ | 0 | 4 |
| 0.01 | +3.34 × 10⁻¹¹ | +2.65 × 10⁻¹¹ | 0 | 4 |
| 0.05 | +1.05 × 10⁻⁷ | +8.32 × 10⁻⁸ | 0 | 4 |
| 0.1 | −3.68 × 10⁻⁸ | +2.91 × 10⁻⁸ | 0 | 4 |
| 0.3 | +1.31 × 10⁻⁷ | +1.04 × 10⁻⁷ | 0 | 4 |

**All cases converge to amplitudes ~10⁻⁷–10⁻¹⁰ — effectively trivial.**
The BVP solver finds the linear Bessel J_1 modes (which satisfy
β(R_MAX)=0 trivially at small amplitude) and the nonlinear cubic correction
is negligible. Sign changes = 4 indicates we're landing on the 5th radial
mode of the linear Bessel equation, not a nonlinear soliton.

### Robin BC: same behavior

Tested Robin BC β'(R_MAX) + k·β(R_MAX) = 0 for k ∈ {0.5, 1.0, 2.0}:
all cases converge to trivial amplitudes; tail/peak ratio ≈ 0.17 at best.

### Negative-λ scan

Tested λ < 0 (where linear far-field gives K_1 exponential decay rather
than J_1 oscillation):

| λ | Converged B0 | Peak | tail/peak | m_χ (MeV) |
| --- | --- | --- | --- | --- |
| -0.1 | 7.3 × 10⁻⁸ | 7.5 × 10⁻⁸ | 1.9 × 10⁻¹⁷ | 0.194 |
| -0.5 | 3.7 × 10⁻⁸ | 3.4 × 10⁻⁸ | 3.2 × 10⁻¹⁹ | — |
| -1.0 | 2.2 × 10⁻⁸ | 2.0 × 10⁻⁸ | 1.5 × 10⁻¹⁹ | — |
| -2.0 | 1.1 × 10⁻⁸ | 1.0 × 10⁻⁸ | 1.2 × 10⁻¹⁹ | — |

**With λ < 0, tail/peak ratios are excellent (10⁻¹⁷ to 10⁻²¹) — clean
exponential decay** — but amplitudes are still trivial (10⁻⁸). So even
with the "correct" sign of effective mass, the cubic NLS in cylindrical
m=1 geometry doesn't have a stable soliton.

**Why:** this is the 2D Townes-soliton pathology. Cubic nonlinearity in
2D cylindrical geometry is critical — solitons exist mathematically but
are saddle points in energy (unstable). The BVP relaxation finds the
trivial (amplitude → 0) energy minimum, not the Townes-saddle solution.

### B verdict

The canonical 2-function neutral ODE does NOT support a stable localized
ground state. Sonnet's "lightest at λ=1, m_χ = 0.998 MeV" from IVP is
a windowed-integration value on an oscillating linear-Bessel-mode tail,
not a true bound state.

---

## Interpretation A — 4-function BVP with Q_CS=0

**Code:** `sandbox_v9/m6_v9_neutral_4fn_bvp.py` (forked from
`sandbox_v6/m6_v6_4fn_calibrated_bvp.py`)

Single change from v6: `I_TARGET = 0.0` instead of `1.0` (Q_CS=0 for
neutral instead of Q_CS=1 for charged). Same Lagrange-multiplier
formulation, same asymmetric helicity exp(-r) seed (V₀=+0.1, A₀=−0.1,
Q₀=+0.1, J₀=−0.1).

### Results: two qualitative regimes

| Init (ω, λ_LM) | Converged ω | Converged λ_LM | m_eff² | Peak amps | Nodes V/A/Q/J | tail @r≥8 | Status |
| --- | --- | --- | --- | --- | --- | --- | --- |
| (1.0, 0.1) | 2.506 | −0.241 | −5.78 | 0.10/0.06/0.56/0.01 | 15/12/15/22 | 0.127 | Highly excited |
| (0.1, 0.1) | 0.667 | +2.383 | +0.056 | 0.10/0.006/0.10/0.02 | 6/5/6/0 | 0.0514 | Borderline; excited mode |

**The 4-function BVP DOES find non-trivial solutions for Q_CS=0** — the
extra fields (V, A, Q, J) and the Lagrange multiplier coupling provide
enough freedom for the solver to converge to a non-zero configuration.
But:

- All converged solutions have many nodes (5-22) — none is a clean
  ground state with ≤4 nodes
- The "best" case (ω_init=0.1) has tail/peak = 5%, just over typical
  localization criteria
- Q_J normalization for neutral mass extraction is undefined under the
  Q_CS=0 framework (H/Q_CS → ∞ trivially); no clean m_χ comes out

### A verdict

The 4-function BVP produces SOMETHING for Q_CS=0 but it's not a clean
neutral chaoiton ground state. The asymmetric helicity init that gave
clean Q_CS=1 ground states in v5/v6 doesn't have an obvious Q_CS=0
counterpart — the natural Q_CS=0 configuration is A=J=0 (no field at
all) or trivially symmetric, neither of which is the "asymmetric
helicity" Paul prescribes.

---

## Combined finding (both interpretations agree on the structural issue)

| Approach | What we found | Implication |
| --- | --- | --- |
| 2-function BVP, Dirichlet, λ>0 | Trivial → 0 amplitude | Linear Bessel modes only, no nonlinear soliton |
| 2-function BVP, Dirichlet, λ<0 | Trivial → 0 amplitude with clean K_1 decay | Townes-soliton pathology: critical 2D cubic NLS has unstable saddle, BVP relaxation finds trivial |
| 4-function BVP, Q_CS=0, asymmetric helicity | Excited modes (6-22 nodes), tail 5-13% | No ≤4-node ground state in this ansatz space |
| Sonnet's IVP (reference) | β profile with oscillating tail; m_χ = 0.998 MeV at λ=1, B0=0.001 | Windowed-integration artifact, not a true ground state (confirmed Q42) |

**The neutral chaoiton "ground state" cited in 9c §8 + §9 criterion 7
(m_χ = 0.998 MeV) does not survive proper BVP scrutiny in either
interpretation.** The number is an artifact of integrating an oscillating
profile up to a finite cutoff R_max=12.

---

## Two structural questions for Paul (in email v12)

### Q43 — Sign convention or missing term?

In Sonnet's canonical `neutral_ode`:

```text
β'' + β'/r - β/r² + λβ + 4gβ³ = 0
```

The +λβ term acts as a positive mass-squared in standard Klein-Gordon
sign convention, which would normally allow K_1 exponential decay at
infinity. But the sign in this equation (when written as β'' + ... + λβ
= 0) means the linear far-field is **J_1 Bessel oscillatory**, not K_1
decaying. To get K_1 decay we'd need λβ on the OTHER side (or the equation
to be β'' + β'/r - β/r² - λβ + 4gβ³ = 0).

Either the canonical ODE has a sign error, OR our convention for the
"mass" is opposite from standard. Concrete question: should the equation
be:

```text
β'' + β'/r - β/r² - λβ + 4gβ³ = 0   (mass-like term as standard)
```

instead of

```text
β'' + β'/r - β/r² + λβ + 4gβ³ = 0   (Sonnet's current form)
```

### Q44 — Geometry: cylindrical 2D vs spherical 3D

The −1/r² term in Sonnet's β-equation is the 2D cylindrical m=1 azimuthal
component (∇²_⊥ for the φ-component of a vector field in 2D). In 2D, the
cubic NLS is critical and has only unstable (Townes) solitons.

For a STABLE neutral chaoiton in classical field theory, we typically
need 3D spherical geometry (where cubic NLS is subcritical and stable
solitons exist). Is the canonical neutral chaoiton:

| Option | Geometry | ODE |
| --- | --- | --- |
| (a) 2D cylindrical m=1 (Sonnet's current) | r·dr volume | β'' + β'/r - β/r² + ... |
| (b) 3D spherical s-wave | 4π r² dr volume | β'' + 2β'/r + ... |
| (c) 3D spherical p-wave (l=1) | r²·dr volume | β'' + 2β'/r - 2β/r² + ... |

(b) is the most natural for a stable 3D soliton. (c) is also viable.
(a) — Sonnet's current — is mathematically problematic for stable
solitons.

If the answer is (b) or (c), the canonical 2-function reduction needs
revision, and so do the charged-sector mass scaling laws (because the
electron H/Q at 0.090% gap from v8 step 2 was also on the 2D cylindrical
form).

---

## Email v12 plan

Three asks (drafted in next conversation segment):

| Ask | Content |
| --- | --- |
| 1 | Report v9 findings: both interpretations fail to produce a clean neutral ground state under proper BVP. m_χ in 9c §8 is IVP-windowed artifact. |
| 2 | Ask Q43 (sign convention) and Q44 (cylindrical-vs-spherical geometry) |
| 3 | Offer to implement either fix immediately if Paul/DeepSeek can confirm direction |

---

## DeepSeek's reply 2026-05-22 10:03 AM — both Q43 + Q44 confirmed

Paul forwarded DeepSeek's reply to email v12. Quoted verbatim:

> *"Rodrigo, Thank you for the detailed BVP attempts. Your Q43 and Q44
> are exactly right.*
>
> ***Q43 – Sign convention:** The linear term in the neutral ODE should
> have a MINUS sign for exponential decay. The correct form is:*
>     *β'' + (2/r)β' - (2/r²)β - m_J² β + 4g β³ = 0*
> *(Note the minus sign on m_J².)*
>
> ***Q44 – Geometry:** The correct geometry is 3D spherical, not 2D
> cylindrical. For l=1 (p-wave), the ODE is as above (coefficient 2/r
> for the first derivative, and centrifugal -2/r²). This is subcritical,
> so stable solitons should exist."*

DeepSeek attached `neutral_bvp_solver.py` (saved verbatim at
`sandbox_v9/neutral_bvp_solver.py`).

### DeepSeek's verbatim template — collapses to trivial

Template BCs: β'(R_MIN) = 0 + Robin at R_max. 2 BCs, 2 states, no free
parameters. With B0 scan as initial guess (B0 ∈ [0.001, 0.5]):

| B0_init | Converged β(0) | β'(0) | Note |
| --- | --- | --- | --- |
| All 9 values | 0.000000 | 0.000000 | All converge to trivial β ≡ 0 |

**Diagnosis:** β'(R_MIN) = 0 is the wrong BC for l=1 (p-wave). The
regular l=1 solution behaves as β(r) ~ B0·r near origin, so
β(0) = 0 and β'(0) = B0 ≠ 0. DeepSeek's BC forces β(r) ~ B0·r² instead,
which combined with trivial nonlinear effects collapses to β ≡ 0.

### Modification: l=1 BC + B0 fixed + m_J free → ground state found

`sandbox_v9/neutral_bvp_solver_mJ_free.py` fixes:

- Origin BC: β(R_MIN) = B0·R_MIN and β'(R_MIN) = B0 (proper l=1 regular)
- Far-field BC: Robin β'(R_max) + m_J·β(R_max) = 0 (unchanged)
- B0 specified, m_J FREE eigenvalue (per DeepSeek's offer)

**RESULT: ground-state family found** at g=1.0, R ∈ [0.02, 20]:

| B0 | Converged m_J | Peak β @ r | tail/peak | Sign changes | H/Q (3D) | Status |
| --- | --- | --- | --- | --- | --- | --- |
| 0.40 | +0.460 | 0.62 @ 1.99 | 3.7×10⁻⁵ | 0 | 1.694 | ✅ Ground state |
| 0.45 | +0.488 | 0.66 @ 1.88 | 2.0×10⁻⁵ | 0 | 1.906 | ✅ Ground state |
| 0.50 | +0.514 | 0.70 @ 1.78 | 1.1×10⁻⁵ | 0 | 2.117 | ✅ Ground state |
| 0.55 | +0.540 | 0.73 @ 1.70 | 6.3×10⁻⁶ | 0 | 2.329 | ✅ Ground state |
| 0.60 | +0.564 | 0.76 @ 1.62 | 3.7×10⁻⁶ | 0 | 2.541 | ✅ Ground state |
| 0.70 | −0.609 | 0.82 @ 1.50 | 3.6×10⁻⁵ | 1 | 2.964 | ⚠️ Excited |
| 0.80 | −0.651 | 0.88 @ 1.41 | 1.5×10⁻⁵ | 1 | 3.388 | ⚠️ Excited |

### Net: TRUE NONLINEAR LOCALIZED SOLITON FOUND

| Finding | Value |
| --- | --- |
| Sign-changes in ground state | **0** (true ground state, l=1 p-wave) |
| Far-field tail/peak ratio | **10⁻⁵ to 10⁻⁶** (clean K_1 decay) |
| m_J at B0=0.5 (representative) | +0.5145 natural units → **0.531 MeV** physical (× ℏc/R_phys) |
| Peak β location (B0=0.5) | r = 1.78 natural units → **340 fm** physical |
| H/Q (3D spherical r²·dr integration, B0=0.5) | 2.117 |
| Continuous family in B0 ∈ [0.40, 0.60] | All give clean ground states with m_J ∈ [0.46, 0.56] |
| Transition to excited branch | B0 between 0.60 and 0.70 (m_J flips sign) |

**The 3D spherical l=1 form is subcritical (as DeepSeek predicted), and
the nonlinear cubic NLS DOES support stable solitons.** Email v13 (next)
reports findings + asks how to pick the canonical (g, B0) point in the
family for the DM-paper m_χ extraction.

### Open question for email v13

The BVP has a 1-parameter family of ground states (parameterized by B0).
To pick THE canonical neutral chaoiton for the DM paper, we need an
additional constraint. Three candidates:

| Option | What to pin |
| --- | --- |
| (a) Fix m_J at electron-calibrated value (m_J = 1.0 in natural units, matching the charged sector's λ) | Search for the B0 in the family that gives m_J = 1.0 — but our scan shows B0=0.60 gives m_J=0.56, suggesting we may need a different g value, OR the electron's λ doesn't directly map to neutral's m_J |
| (b) Pick the lightest (smallest H/Q × m_e) | B0=0.40 in our scan: H/Q = 1.694, m_χ ≈ 1.694 × m_e = 0.866 MeV. But subject to 3D-spherical-vs-cylindrical normalization issue. |
| (c) Self-consistency: same (g, m_J) as the charged sector and find which B0 lands | Asks Paul/DeepSeek which is canonical |

---

## Cross-references

- `0b_M6_roadmap.md` — v9 sandbox added to sequence
- `0b_model_gates.md` — G2 status updated (now "BVP exhausted; no clean ground state")
- `0b_question_tracker.md` — Q43, Q44 added as IMMEDIATE
- `0c_sandbox_v8.md` — Q42 caveat that triggered v9 work
- `sandbox_v9/m6_v9_neutral_2fn_bvp.py` — Interpretation B code
- `sandbox_v9/m6_v9_neutral_4fn_bvp.py` — Interpretation A code (fork of v6)
