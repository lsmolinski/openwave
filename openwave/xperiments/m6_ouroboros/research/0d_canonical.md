# M6 / Ouroboros — Canonical Numerical Specification

**Purpose:** the authoritative reference doc consolidating the winning
numerical methodology across v1 → v10. If you want to reproduce M6's
lepton spectrum, neutral chaoiton ground state, and DM paper inputs from
scratch — this is the canonical specification. The per-iteration
`0c_sandbox_v*.md` files document the back-and-forth, what failed, and
why; this doc is just the canonical implementation.

Last updated: 2026-06-05 (post sandbox_v11 — exact scaling symmetry of the
neutral sector found; (g, λ) scan answered analytically + numerically; Q47
resolved as scaling identity; DM-paper "BVP inconsistency" closed).

---

## 1. The Lagrangian (Werbos Ouroboros)

Two coupled vector fields `A_μ, J_μ` on Minkowski spacetime with metric
`(−, +, +, +)`. The Lagrangian (v9 LoE paper §2):

```text
L = -¼ F_μν F^μν - ¼ G_μν G^μν + m_J² A_μ J^μ - f(J_μ J^μ)
```

where:

| Symbol | Meaning |
| --- | --- |
| `F_μν = ∂_μ A_ν − ∂_ν A_μ` | A-field curl tensor (Maxwell-like) |
| `G_μν = ∂_μ J_ν − ∂_ν J_μ` | J-field curl tensor |
| `m_J` | J-field mass parameter (mediator mass) |
| `f(J·J)` | Higgs-like potential. Canonical choice: `f(s) = (g/4) s²` for the electron sector; specific form per `0c_sandbox_v6.md` discussions |
| `g` | Coupling constant — canonical g = 1.0 for both charged + neutral after v8/v10 calibration |

In the linear limit (J → 0), the A-field reproduces Maxwell exactly
(structural fact). Nonlinear coupling produces stable time-periodic
solitons ("chaoitons") in two sectors:

| Sector | Topological charge | Examples |
| --- | --- | --- |
| Charged | Q_CS = 1 | Electron (ω=1), muon (ω=12.82), tau (ω=50.0), pion+ candidate (ω=15.0) |
| Neutral | Q_CS = 0 | Dark matter candidate (m_χ ≈ 0.46 MeV) |

The two sectors use **different geometric realizations** of the chaoiton
ansatz — this is critical and discovered empirically through v1 → v9
(see `0c_sandbox_v8.md` step 3 paper-math comparison).

---

## 2. Charged sector recipe (electrons, muons, taus)

### 2.1 Ansatz — 2-function (α, β) reduction in cylindrical geometry

The canonical charged-sector ansatz (per Q34 resolved 2026-05-21):

```text
A_0 = α(r) cos(ωt),     A = 0
J_0 = 0,                J = β(r) sin(ωt) φ̂
```

with **vector cylindrical Laplacian** (f'' + f'/r − f/r²) and **r dr volume
measure** (2D cylindrical projection of the 3D toroidal chaoiton). The
α-field is the toroidal scalar; β is the poloidal scalar.

### 2.2 ODE (charged sector)

Reduced equations after the ansatz substitution:

```text
α'' + α'/r - α/r²            + (ω² − m_J²) α   + 4g α (α² + β²) = 0
β'' + β'/r - β/r² + 2ωα     − m_J² β            + 4g β (α² + β²) = 0
```

The cross-coupling `2ωα` term in the β-equation is the chiral coupling
that produces Q_CS = 1. Sign convention: minus on the mass terms.

### 2.3 Boundary conditions

| Boundary | Condition | Why |
| --- | --- | --- |
| r → 0 | slope BCs: `α'(0) = 0`, `β'(0) = B0` for some B0 | Regular at origin for vector cylindrical Laplacian |
| r → ∞ | decay: α(r), β(r) → 0 | Localized soliton |
| r_inner offset | `r_inner = 0.02` in solve_ivp | Avoid 1/r singularity numerically |

### 2.4 Solver — forward IVP

Use `scipy.integrate.solve_ivp` with RK45 from r = 0.02 outward, starting
with α(r_inner), β(r_inner) seeded by exp-like profiles consistent with
the slope BCs.

**Why IVP works here (and DOESN'T for neutral):** the charged 2-function
system has well-conditioned forward integration from the origin slope.
Forward shooting + ω sweep finds the discrete lepton spectrum. This was
empirically validated in sandbox_v4–v9; sandbox_v8 confirmed the canonical
form via Paul's `ouroboros_benchmark.py`.

### 2.5 Calibration anchor (electron)

| Anchor | Value |
| --- | --- |
| ω (electron) | 1.0 (sets the unit) |
| g (electron) | 1.0000 (from sandbox_v8 step 2 scan: gives 0.090% gap to H/Q=1.6969 target) |
| H/Q (at g=1.0000, ω=1.0) | 1.6890 (matches target 1.6969 to 0.090%) |
| m_e | 0.5109989461 MeV (PDG) |
| R_phys | ℏc / m_e = 386.16 fm (electron Compton wavelength) |

The H/Q value calibrates the dimensionless ratio that maps to electron
mass. Other particles use the same anchor + frequency sweep.

### 2.6 Lepton spectrum (ω-sweep)

Same ODE + same g, just sweep ω. Three lowest stable Q_CS=1 modes match
the lepton families:

| Particle | ω | Predicted mass (MeV) | PDG mass (MeV) | Gap |
| --- | --- | --- | --- | --- |
| Electron | 1.00 | 0.511 (calibration anchor) | 0.5110 | 0.090% |
| Muon | 12.82 | 106.50 | 105.66 | 0.80% |
| Tau | 50.0 | 1661.88 | 1776.86 | 6.47% |
| Pion+ (candidate) | 15.0 | 144.11 | 139.57 | 3.25% |

Mass scaling is roughly H ~ ω^2.04 (empirical) — close to but not exactly
Werbos's heuristic ω^2.22 = log(207)/log(11). Discrete-spectrum mechanism
(why exactly these ω and no others) remains an open question (Q2).

### 2.7 Reference script

`sandbox_v8/ouroboros_benchmark.py` (Paul Werbos + Claude Sonnet 4.6).
Runnable first try. Implements both `electron_calibration_scan()` and
`lepton_spectrum_scan()` functions. Use this as the canonical reference
for the charged sector.

---

## 3. Neutral sector recipe (dark matter candidate)

### 3.1 Ansatz — 1-function (β only) in 3D spherical l=1 p-wave geometry

The charged sector's 2-function ansatz **does NOT** carry over to the
neutral sector. The canonical neutral ansatz (per Q43+Q44 resolved
2026-05-22):

```text
β(r) only — α = 0 (no charge to support);
3D spherical l=1 p-wave geometry, r² dr volume measure
```

The neutral chaoiton is a **structurally different** field theory from the
charged one — different Laplacian, different volume measure, different
dimensionality. v9 phase 1 confirmed empirically that "use the same ODE
as charged" doesn't produce a clean ground state (Townes-soliton pathology
in 2D cylindrical; excited-modes-only in 4-function).

### 3.2 ODE (neutral sector, l=1 p-wave)

DeepSeek's confirmed Q43+Q44 form:

```text
β'' + (2/r) β' - (2/r²) β - m_J² β + 4g β³ = 0
```

Note the **minus sign on m_J²** (Q43 — produces K_1 exponential decay,
not J_1 Bessel oscillation) and the **3D spherical Laplacian with l=1
centrifugal term** (Q44 — `(2/r)β'` + `−2β/r²` for p-wave; subcritical
NLS, supports stable solitons).

### 3.3 Boundary conditions

| Boundary | Condition | Why |
| --- | --- | --- |
| r → R_MIN (e.g., 0.02) | `β(R_MIN) = B0 · R_MIN` AND `β'(R_MIN) = B0` | l=1 regular: β ~ B0·r near origin |
| r → R_max (e.g., 20) | Robin: `β'(R_max) + m_J · β(R_max) = 0` | K_1 decay matching |

**Critical:** DeepSeek's verbatim template used `β'(R_MIN) = 0`, which
forces `β ~ B0·r²` and collapses to trivial β ≡ 0. The correct l=1
regularity class is `β ~ B0·r` (slope = B0, value = 0 at origin).

### 3.4 Solver — BVP with m_J as free eigenvalue

Use `scipy.integrate.solve_bvp` with:

- **B0 specified** (1-parameter family in B0)
- **m_J treated as free eigenvalue** (p[0]) — solver finds the value
  that supports the soliton
- 3 BCs: 2 at origin + 1 Robin at R_max
- Initial guess: β(r) ≈ B0 · r · exp(−m_J_init · r), with m_J_init = 1.0

### 3.5 Canonical point (per Paul/DeepSeek Q45+Q46 reply)

| Parameter | Value | Source |
| --- | --- | --- |
| g | 1.0 | Same as electron (Q45: "same Lagrangian parameters") |
| B0 | 0.5 | v9 phase 2 reference point; works across [0.40, 0.60] for ground state |

At this point: m_J_raw = +0.5145, peak β = 0.697 at r = 1.78, sign
changes = 0, tail/peak = 1.1×10⁻⁵ (clean K_1 decay), (H/Q)_spherical =
2.1173, Q_J_spherical = 4.2586, H_spherical = 9.0168.

### 3.6 Ground state validation

A clean ground state must satisfy:

| Diagnostic | Threshold | At canonical (g=1.0, B0=0.5) |
| --- | --- | --- |
| sign_changes | 0 (no nodes) | ✅ 0 |
| tail/peak ratio | < 10⁻⁴ (K_1 decay) | ✅ 1.1×10⁻⁵ |
| Peak β location | r ∈ [1, 3] (compact) | ✅ r = 1.78 |
| solve_bvp.status | 0 (converged) | ✅ 0 |

Family-invariant range: B0 ∈ [0.40, 0.60] all give clean ground states
at g = 1.0. Transition to excited branch (1 node, m_J < 0) at B0 between
0.6 and 0.7.

### 3.7 Exact scaling symmetry (v11 — closes the (g, λ) parameter space)

With λ ≡ −m_J², the substitution β(r) = a·u(x), x = m_J·r, a = m_J/(2√g)
maps the neutral ODE onto the **universal, parameter-free equation**

```text
u'' + (2/x) u' - (2/x²) u - u + u³ = 0
```

Every (g, λ < 0) ground state is an exact rescaling of one profile u(x).
Verified on an 18-point grid in `sandbox_v11/m6_v11_neutral_gl_scan.py`:

| Identity | Consequence |
| --- | --- |
| η = m_J·η_u, η_u = 0.826287 (universal) | m_J_phys = m_e/η_u = **0.6184 MeV for ALL (g, λ)** — parameter-free prediction; explains Q47's 1.21024 = 1/η_u invariance |
| H/Q = K_u·m_J², K_u = 7.99999 ≈ 8 | m_χ = (−λ)^{3/2}·η_u·K_u·m_e = 3.3778·(−λ)^{3/2} MeV — depends only on λ; candidate exact identity H = 8 m_J² Q (Q48) |
| amplitude a = m_J/(2√g) | g is a pure amplitude scale — moves NO mass observable, but does scale the J-A mixing integral (α_JN) since it is linear in β |
| B0 = u'(0)·m_J²/(2√g) | at fixed B0, m_J ∝ g^{1/4} — converts the (g, B0) chart of v9/v10 to the (g, λ) chart |

Canonical point in (g, λ) coordinates: **(g, λ) = (1.0, −0.264697)** —
identical to the v9/v10 point (g=1.0, B0=0.5). λ* is fixed by
m_χ = 0.460 MeV (any g); g = 1.0 by electron universality (Q45). The two
masses CANNOT independently fix (g, λ): m_J carries no parameter
information and m_χ only fixes λ.

---

## 4. DM paper inputs (v10 extraction recipe)

Per DeepSeek's 2026-05-22 7:05 PM reply on Q45+Q46:

### 4.1 Geometry conversion factor η

```text
η = (∫₀^R_max β² · r dr) / (∫₀^R_max β² · r² dr)
```

This is the ratio of cylindrical-measure norm to spherical-measure norm
of the same β profile. Converts H/Q across geometries.

At canonical (g=1.0, B0=0.5): **η = 0.4251**.

### 4.2 m_χ extraction

```text
m_χ = η × (H/Q)_spherical × m_e
    = 0.4251 × 2.1173 × 0.5110 MeV
    = 0.4599 MeV
```

### 4.3 m_J in physical units

```text
m_J_corrected (natural) = m_J_raw / η = 0.5145 / 0.4251 = 1.2102
m_J (physical) = m_J_corrected × m_e = 1.2102 × 0.5110 MeV = 0.6184 MeV
```

**Caveat:** m_J_corrected = m_J/η is empirically family-invariant at
1.21024 across both B0 ∈ [0.10, 0.60] and g ∈ [0.5, 1.6] at clean
ground-state convergence. DeepSeek's "1.0 target" doesn't land — looks
like a Pohozaev-type virial identity, not a free parameter. Q47 sent to
Paul/DeepSeek 2026-05-22 evening to confirm interpretation. The
deliverable numbers stand under either reading.

### 4.4 C extraction from far-field K_1 amplitude

For l=1 spherical mode, the far-field asymptotic is approximately:

```text
β(r) ~ C_natural · exp(-m_J · r) · (1 + 1/(m_J · r)) / r,  r → ∞
```

Fit the converged β(r) tail (r ≥ 8 in natural units) to this form via
nonlinear least squares. Extract C_natural, then convert:

```text
C (MeV·fm) = C_natural · R_phys · m_e
           = 3.903 · 386.16 · 0.5110 MeV·fm
           = 770 MeV·fm
```

Fit residual at canonical point: ~5% max relative error. Tightens with
larger `r_min_fit` if needed.

### 4.5 Final DM paper inputs

| Quantity | Value | How |
| --- | --- | --- |
| m_χ | **0.4599 MeV** | η × (H/Q) × m_e |
| m_J | **0.6184 MeV** | m_J_corrected × m_e (post-Q46 geometry) |
| C | **770 MeV·fm** | K_1 far-field fit + unit conversion |
| Sub-MeV DM | Yes | m_χ < m_e by ~10% (0.46 vs 0.51 MeV) |

---

## 5. Calibration anchors + unit conversions

Quick reference for converting between natural units and physical:

| Constant | Value | Use |
| --- | --- | --- |
| m_e | 0.5109989461 MeV | Electron mass (PDG); base anchor |
| ℏc | 197.3269804 MeV·fm | Quantum-relativistic unit conversion |
| R_phys | ℏc / m_e = 386.16 fm | Compton wavelength of electron; length unit |
| ω (electron) = 1 | Sets natural frequency unit | ω_phys = ω_natural × m_e c² / ℏ |

For mass extraction: any dimensionless H/Q (natural) × m_e gives MeV.
For length extraction: any dimensionless length × R_phys gives fm.

---

## 6. Things that DON'T work (briefly — for future-me / others)

| What | Doesn't work because | Reference |
| --- | --- | --- |
| Forward IVP for 4-function (V, A, Q, J) BVP | Generic forward shooting from value-BC origin doesn't reach Q_CS=1 chaoiton in any parameter region tested (6 attempts T6–T9) | `0c_sandbox_v4.md` |
| Forward IVP for neutral chaoiton | Produces oscillating-tail profile, not localized ground state (β/K_1 ratio diverges 10 orders in tail). 9c-cited m_χ = 0.998 MeV was windowed-integration artifact | `0c_sandbox_v8.md` Q42 |
| 2D cylindrical for neutral | Cubic NLS in cylindrical m=1 geometry is critical → only unstable Townes solitons; BVP relaxes to trivial β ≡ 0 | `0c_sandbox_v9.md` phase 1 |
| Sign `+m_J²β` in neutral ODE | Gives J_1 Bessel oscillatory far-field, not K_1 exponential decay | Q43, `0c_sandbox_v9.md` |
| DeepSeek's verbatim template `β'(R_MIN) = 0` | Forces β ~ B0·r² near origin (wrong l=1 regularity class); collapses to trivial | `0c_sandbox_v9.md` phase 2 |
| Linear BVP with B0 as free parameter | Amplitude collapses to ~0 (linear-eigenvalue degeneracy) | `0c_sandbox_v9.md` phase 2 |
| g-scan to force m_J_corrected = 1.0 | m_J/η is virial-invariant at 1.21 across B0 and g; brentq lands in excited modes | `0c_sandbox_v10.md` |
| v1's H/Q = 1.6918 as "independent reproduction" | Different geometric realization (3D spherical scalar Laplacian) than canonical (2D cylindrical vector Laplacian); numerical coincidence, not derivation | Q36, `0c_sandbox_v8.md` step 3 |

---

## 7. Reference scripts (canonical implementations)

| Script | Sector | Purpose |
| --- | --- | --- |
| `sandbox_v8/ouroboros_benchmark.py` | Charged | Lepton spectrum + electron calibration; runnable canonical impl |
| `sandbox_v9/neutral_bvp_solver_mJ_free.py` | Neutral | True ground-state BVP with proper l=1 BC + B0 fixed + m_J free |
| `sandbox_v10/m6_v10_canonical_neutral_chaoiton.py` | Neutral | DM paper input extraction (η + m_χ + m_J + C) per Paul/DeepSeek Q45+Q46 recipe |
| `sandbox_v11/m6_v11_neutral_gl_scan.py` | Neutral | (g, λ) scan + scaling symmetry + λ* locus + Technical-Note failure-mode reproduction + β(r) profile export (`v11_canonical_beta_profile.csv`) |
| `sandbox_v8/extract_mJ_C_mchi.py` | Neutral (legacy v8) | Earlier extraction using IVP; superseded by v10 (windowed-integration artifact per Q42) |

---

## 8. Open canonical-interpretation questions (live)

| Question | Status |
| --- | --- |
| Q47 — m_J_corrected = m_J/η family-invariance | ✅ ANALYTICALLY RESOLVED (v11): 1.21024 = 1/η_u, exact consequence of the scaling symmetry (§3.7); not tunable |
| Q48 (NEW) — H = 8 m_J² Q exact identity? | Holds to 5 digits across the (g, λ) grid; candidate Pohozaev/virial identity for the 3D l=1 cubic NLS; sent to Paul/DeepSeek for proof |
| Q49 (NEW) — α_JN mixing-integral definition | β(r) profile delivered (`v11_canonical_beta_profile.csv`); need exact M_JA(0) definition (integrand + measure + normalization) to compute |
| Q2 — discrete ω selection mechanism | Empirically validated (3 lowest modes = leptons within targets); analytic derivation deferred |
| Q3 — closed-form ω = 2mc²/ℏ | Not derived; calibration-only |
| Q6 — QCD reconciliation (3-chaoiton proton) | Future work; pion+ at 3.25% is suggestive |
| Q35 — active neutrinos in framework | Sterile-like? Open |

---

## 9. Cross-references

- `0a_background.md` — initial 2026-05-15 evaluation
- `0b_M6_roadmap.md` — sandbox sequence + current state
- `0b_model_gates.md` — G1/G2/G3 production criteria
- `0b_question_tracker.md` — live question + hardest-pieces tracker
- `0c_sandbox_v8.md` — Sonnet's canonical script + 5-step work program
- `0c_sandbox_v9.md` — neutral BVP variants A/B + phase 2 ground-state-found
- `0c_sandbox_v10.md` — Paul/DeepSeek Q45+Q46 recipe applied; DM paper inputs delivered
- `0c_sandbox_v11.md` — (g, λ) scan + scaling symmetry; DM-paper Technical-Note "BVP inconsistency" closed
- `theory/` — all Werbos corpus (LoE v1-v9c, Lean theorem, Numerical Benchmark, Hopf invariant proof)
- Werbos's 9c LoE paper (Zenodo 20330894) — cites Reference [17] = our GitHub repo
