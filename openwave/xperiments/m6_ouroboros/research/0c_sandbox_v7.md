# 2026-05-21 — Sandbox v7: Ground-state BVP per Paul's Q28/Q29/Q30 reply

Triggered by Paul's 2026-05-21 morning email answering email v6's three
questions. Paul's answers were clear and actionable; v7 implementation went
fast but hit a structural wall on mode-selector basins. **Net result: v7
explored 7 variants and could not reach the (V₀>0, Q₀>0, A₀<0, J₀<0, ≤4-node)
basin Paul described. v6.6 H_electron/Q = 1.7112 (0.84% off target) remains
the best result; new email v7 sent asking Paul for algorithm specifics.**

**Source:** Paul Werbos email 2026-05-21 morning (replying to email v6).

**Status as of 2026-05-21 morning:** v7 BLOCKED on mode-selector. Solver
landscape has multiple wrong-sign basins; warm-start + sign anchors all
fail to land Paul's prescribed ground state. Email v7 sent asking for
solver details. Calibration accuracy unchanged from v6 (0.84% off via
drop-quartic). M6 production gates G1/G2 remain "nearly unblocked"; v7
work parked until Paul replies.

---

## Source context — Paul's 2026-05-21 reply (verbatim Q28/Q29/Q30)

### Q28 — Which quartic IS canonical for the electron H?

> *"None. For the electron (the fundamental Q=1 chaoiton at ω≈1, g=1.0625,
> λ=1.0), the quartic term is indeed negligible. It contributes <0.1% to
> the energy. In our production code, we included the full quartic in the
> equations of motion (to get the correct field shapes), but when computing
> H/Q for the electron, the quartic contribution is so small that it can
> be dropped without affecting the calibration. Your 'no quartic' variant
> is therefore the correct effective Hamiltonian for matching H/Q = 1.6969."*
>
> *"For the muon and tau (higher ω), the quartic becomes significant. For
> the neutral chaoiton (Q=0), the quartic is essential. So the canonical
> form for H is the full expression (with quartic), but for the electron
> alone, the quartic term is a fine-tuning that can be ignored."*

| Particle | Quartic role in H |
| --- | --- |
| Electron (Q=1, ω≈1) | Negligible — drop for H/Q matching |
| Muon (Q=1, ω≈12.78) | Significant — keep full |
| Tau (Q=1, ω≈40.7) | Significant — keep full |
| Neutral chaoiton (Q=0) | Essential — keep full |

### Q29 — Zero crossings: ground state vs excited mode

> *"Our converged ground state for the electron has V and Q with exactly 4
> zero crossings (excluding r=0). The Lean stability specification (≤4 nodes)
> refers to the radial functions excluding the origin. If your v6.6 has 5
> crossings, you are on the first excited state. That explains why your H/Q
> is slightly higher (1.7112 vs 1.6969) — excited states have higher energy."*
>
> *"To get the true ground state, initialise with a single-node profile (e.g.,
> the reference profile I gave earlier, which decays monotonically after a
> single peak). Use that as the initial guess in solve_bvp with the integral
> constraint Q_CS=1. The solver should settle into the 4-node ground state,
> and H/Q will drop to 1.6969 (or even lower — your current 1.7112 is already
> very close)."*

### Q30 — Sign of Q at origin

> *"Our ground state has Q(0) = +0.1 (positive), same sign as V(0). The
> asymmetric helicity prescription is:*
>
> *V₀ > 0, Q₀ > 0, A₀ < 0, J₀ < 0*
>
> *Your v6.6 giving Q(0) negative suggests you are in a different basin —
> likely because your initial guess or the sign of the Lagrange multiplier
> term flipped the relative sign. To correct, set the initial profiles
> such that V(r) and Q(r) have the same sign (both positive) at r=0, and
> A(r) and J(r) have the same negative sign. The reference profile I gave
> earlier has exactly that."*

### Strategic summary

Paul's closing: *"Once you implement these, you should get H/Q = 1.6969
(±0.2%) and ω very close to 1.0. Then you can proceed to the lepton scan
(warm-start from the electron) and the neutral scan."* and *"We are very
close. No need to rush — take your time."*

This maps cleanly onto **Scenario B** from our v6 v7-plan: quartic
negligible for electron H, ≤4-node ground state, Q(0)>0 helicity. v7
implementation should be a mode-selector that lands the ground state.

---

## Today's tasks (planned 2026-05-21 morning)

Phased schedule drafted before starting work. Phase 1 was the critical-path
gate; Phases 2-4 contingent on Phase 1 success.

### Phase 1 — Electron ground state (~2-3 hrs)

| # | Task | Expected outcome |
| --- | --- | --- |
| 1 | Create `sandbox_v7/m6_v7_4fn_ground_state_bvp.py` — fork v6.6 | Same ODE structure (full quartic in EoM); dual H computation (H_full + H_electron) |
| 2 | Implement Paul's 8-point reference profile as warm-start | Single-peak monotonic-decay initial shape |
| 3 | Enforce sign convention V₀=+0.1, Q₀=+0.1, A₀=−0.1, J₀=−0.1 | Lands in correct basin (Q(0) > 0 at convergence) |
| 4 | Dual H computation: H_full (with quartic) + H_electron (no quartic) | Captures Paul's electron-specific simplification |
| 5 | Run convergence; verify nodes V/Q ≤ 4, Q(0) > 0, H_electron/Q within 0.2% of 1.6969 | Electron ground state locked |
| 6 | Fallback if doesn't land cleanly | Iterate on profile / add Q(0)>0 BC / node-penalty regularization |

### Phase 2 — Lepton continuation scan (~2-3 hrs)

| # | Task | Expected outcome |
| --- | --- | --- |
| 7 | Build continuation wrapper: warm-start from converged electron, nudge ω | Replaces cold-start (step 2 diagnostic failed) |
| 8 | Sweep ω ∈ [1, 50]; find stable modes at muon (≈12.78) and tau (≈40.7) | Lepton spectrum |
| 9 | Use FULL H (with quartic) per Paul Q28 for muon/tau | Different H rule per particle |
| 10 | Compute m_μ/m_e and m_τ/m_e; compare to 207, 3477 | **G1 closure** |

### Phase 3 — Q_A≈0 neutral chaoiton DM scan (~1-2 hrs)

| # | Task | Expected outcome |
| --- | --- | --- |
| 11 | Configure V₀, Q₀ → small ε; A₀, J₀ asymmetric | Q_A ≈ 0 sector |
| 12 | Sweep m_J² for stable bound state | Lands m_χ |
| 13 | Extract m_χ, m_J, σ/m | **G2 closure** — Section 4 numbers |
| 14 | Use FULL H (with quartic) per Paul Q28 for neutral | Different H rule |

### Phase 4 — Email Paul + doc sync (~1 hr)

| # | Task | Expected outcome |
| --- | --- | --- |
| 15 | Draft email v7 with ground state H/Q + lepton ratios + Section 4 numbers | Data drop to Paul |
| 16 | Create 0c_sandbox_v7.md | New sandbox doc |
| 17 | Update 0b_model_gates.md gate status | Gates update |
| 18 | Update 0b_M6_roadmap.md Sandbox v7: planned → DONE | Roadmap sync |

---

## v7 implementation recipe

Script: `sandbox_v7/m6_v7_4fn_ground_state_bvp.py`.

### State vector and ODE — unchanged from v6

Full v6 ODE with quartic in equations of motion (per Paul Q28 — quartic
needed in EoM for correct field shapes regardless of which H we use).

```text
dV/dr  = V'
dV'/dr = -V'/r + Q
dA/dr  = A'
dA'/dr = -A'/r + J + λ_LM·(J + 2·r·J')
dQ/dr  = Q'
dQ'/dr = -Q'/r + V + m_eff²·Q + λ_bench·Q·(Q²-J²)
dJ/dr  = J'
dJ'/dr = -J'/r + A − m_eff²·J − λ_bench·J·(Q²-J²) − λ_LM·(A + 2·r·A')
dI/dr  = r·(A·J' − J·A')        # Q_CS density
```

### Three H functionals computed simultaneously

Per Paul Q28: the canonical H is the full one (for muon/tau/neutral); for
the electron specifically, drop the quartic.

| Form | Use |
| --- | --- |
| `H_full` = (1/2)(V'²+A'²+Q'²+J'²) + (1/2)ω²(V²+A²+Q²+J²) − V·Q + A·J + (g/4)·((V²+Q²)² + (A²+J²)² + 2(V·A−Q·J)²) | Canonical (muon, tau, neutral) |
| `H_electron` = same as H_full minus the (g/4)·… quartic term | **Primary metric for electron** per Paul Q28 |
| `H_benchmark` = v5 Numerical-Benchmark form with (2π)² prefactor | Cross-check only |

### Mode-selector experiments (anchor options)

| Anchor | BCs | Free params | Purpose |
| --- | --- | --- | --- |
| `V` | V(R_MIN) = +V_norm (11 res) | ω, λ_LM (2) | v6 default anti-collapse |
| `Q` | Q(R_MIN) = +V_norm (11 res) | ω, λ_LM (2) | v7 single Q-pin |
| `VQ` + `m_J` free | V(R_MIN)=+V_norm AND Q(R_MIN)=+V_norm (12 res) | ω, λ_LM, m_J² (3) | v7 dual-anchor with m_J² as eigenvalue |
| `VQ` + `lam_bench` free | same | ω, λ_LM, λ_bench (3) | v7 dual-anchor with λ_bench as eigenvalue |

### Strict pass criterion (tightened from v6)

| Test | v6 | v7 |
| --- | --- | --- |
| H/Q precision | < 0.001 of 1.6969 | < 0.003 (0.2% per Paul) on **H_electron**/Q |
| Q_CS precision | < 0.01 of 1.0 | unchanged |
| Max nodes (excluding r=0) | ≤ 4 | strict ≤ 4 |
| Q(0) sign | not checked | **must be > 0** (Paul Q30) |
| A(0), J(0) signs | not checked | **must be < 0** (helicity) |

---

## Attempt log

7 variants tested in ~30 minutes. None landed Paul's prescribed (V₀>0, Q₀>0,
A₀<0, J₀<0, ≤4-node) ground state. Each anchor + warm-start combination
ended in a different wrong-sign basin.

| # | Config | ω | λ_LM | V₀ | A₀ | Q₀ | J₀ | nodes V/Q | H_e/Q | Verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 7.1 | V anchor + warm-start + λ_LM=12 (default) | −1.03 | 185 | +0.10 | −21.2 | +10.6 | +47.9 | 4 / 9 | 34033 | Catastrophic blowup |
| 7.2 | V anchor + no warm-start + λ_LM=0.1 | 1.016 | 12.2 | +0.10 | −0.82 | **−0.12** | −0.20 | 5 / 5 | **1.711** | v6.6 reproduction ✓ |
| 7.3 | V anchor + warm-start + λ_LM=0.1 | 1.01 | 4.2 | +0.10 | **+1.48** | −0.12 | **+0.05** | 5 / 5 | 3.60 | A, J flipped |
| 7.4 | Q anchor + warm-start | 1.59 | 4.3 | **−0.08** | −3.98 | +0.10 | −2.30 | 6 / 5 | 105 | V flipped |
| 7.5 | Q anchor + no warm-start | 1.02 | 6.0 | **−0.05** | −2.65 | +0.10 | **+2.37** | 5 / 5 | 11.6 | V, J flipped |
| 7.6 | VQ anchor + m_J² free | 10⁴⁸ | 10²⁴ | 10¹³ | 10⁵⁶ | 10²⁸ | 10²⁷ | many | 10²⁰⁰ | Singular Jacobian → blowup |
| 7.7 | VQ anchor + λ_bench free | 75 | 1.1 | +0.10 | **+16.3** | +0.10 | −0.79 | 0 / 1669 | 10¹⁰ | λ_bench ran to −360k |

**Reproduction check (7.2):** v6.6 reproduces exactly (ω=1.016, λ_LM=12.21,
Q(0)=−0.12, 5-node V/Q, H_electron/Q=1.7112). Confirms v7's v6.6-baseline
behavior is identical — basin landscape is the same.

**Catastrophic basin (7.1):** Warm-start + Paul's "biased" λ_LM init (12.0
from v6.6's converged value) immediately drove ω negative, λ_LM to 185,
fields to massive values. Suggests the initial gradient at (warm-start
fields, λ_LM=12) points away from any sensible solution.

**Sign-permutation patterns (7.3-7.5):** Each anchor choice fixes one
sign but flips another. The solver appears to have multiple local minima
where exactly one of (V, A, Q, J) has the "wrong" sign relative to Paul's
prescription. The Q_CS=1 integral BC is satisfied across all of these
(integrates whatever fields the solver produces, regardless of signs).

**3rd-eigenvalue failures (7.6, 7.7):** Releasing m_J² or λ_bench as a
third eigenvalue at the dual-pinned BC count caused the eigenvalue to run
unbounded — m_J² to 10⁶⁵, λ_bench to −360k. Both indicate the augmented
Jacobian has a singular direction at the initial guess.

---

## Phase 1 reassessment

### Diagnosis

The solver landscape has multiple wrong-sign basins. **None** of the 7
variants reach Paul's prescribed ground state. v6.6 (5-node "first excited
state" per Paul's Q29) at H_electron/Q = 1.7112 (0.84% off target) remains
the closest result we can produce with the current ODE + BVP formulation.

### What this likely means

| Possibility | Probability | Notes |
| --- | --- | --- |
| Our ODE has a subtle structural difference from Paul's production code | medium | Paul's reference script was broken; we can't diff against it. Reconstruction may be off by a sign or factor. |
| Lagrange-multiplier corrections need different coefficients | medium-high | We derived ΔA = J + λ·(J + 2r·J'), ΔJ = −λ·(A + 2r·A') from H' = H − λ·Q_CS via integration by parts. Alternative derivations exist; coefficient may be off. |
| Solver basin selection needs a continuation method | medium | Start at m_J²=0 (Klein-Gordon limit, well-understood), find solution, slowly homotope to m_J²=0.5. ~1-2 days of careful work. |
| Switch solver to `scipy.optimize.root` method='lm' | possible | Newton-Raphson with hand-coded Jacobian is the alternative Paul mentioned in his 2026-05-20 algorithm reply. ~3-4 hrs to implement; may be the right tool. |
| v6.6's 5-node mode IS the ground state and Paul miscounted | very low | He's confident: *"Our converged ground state for the electron has V and Q with exactly 4 zero crossings."* |

### Reading on today's plan

| Phase | Original estimate | Result |
| --- | --- | --- |
| Phase 1 (electron ground state) | 2-3 hrs | ❌ Hit a wall at ~30 min; 7 attempts, all wrong basin |
| Phase 2 (lepton continuation) | 2-3 hrs | ⛔ BLOCKED on Phase 1 |
| Phase 3 (DM Q_A≈0 scan) | 1-2 hrs | ⛔ BLOCKED on Phase 1 |
| Phase 4 (email + docs) | 1 hr | ✅ Pivots — different email content (asking Paul, not handing off) |

### Strategic options

| Option | Detail | Expected outcome |
| --- | --- | --- |
| **A — Email Paul with algorithm questions** (CHOSEN) | Send concise "we hit a basin problem" email with the 7 attempts summary + three specific algorithm questions. Paul explicitly said *"no need to rush, take your time"*. | He may know the trick immediately. Cheapest path to unblock. |
| B — Try continuation in m_J² | Start at m_J²=0 (Klein-Gordon limit), find solution, slowly homotope to m_J²=0.5. | ~1-2 days; uncertain landing |
| C — Switch to scipy.optimize.root + custom Jacobian | Build Newton-Raphson with finite-difference Jacobian. Paul mentioned as alternative method. | ~3-4 hrs to implement; may or may not reach ground state |
| D — Accept v6.6 H_electron/Q = 1.7112 (0.84% off) | "Good enough" for SABER's purposes (M5 is engineering path). Document v7 attempts, return to M5. | ApJ Section 4 numbers handed off with caveat that we're on excited mode; m_χ DM scan still feasible |

**Selected: Option A.** Paul is actively engaged, explicitly invited
follow-up, and one specific clarifying question saves multi-day exploration
that may not converge. Continue once his reply lands.

---

## Email v7 to Paul (sent 2026-05-21 — after Phase 1 wall)

Reply on the same chain after running 7 v7 variants and finding no path
to the prescribed ground state. Tone: technical, direct, no rushing.
Content:

| Bucket | Summary |
| --- | --- |
| Thanks | Paul's Q28/Q29/Q30 answers were clear and actionable — incorporated all three (drop quartic for electron H, ≤4 nodes ground state, V₀ Q₀ > 0 helicity). |
| Implementation | Built sandbox_v7 with the 8-point warm-start + multiple anchor schemes (V-pin, Q-pin, dual VQ-pin with m_J² or λ_bench as 3rd eigenvalue). |
| Basin findings | 7 variants tried; none reach (V₀>0, Q₀>0, A₀<0, J₀<0, ≤4 nodes). v6.6 (5-node first-excited, H_e/Q=1.7112) remains best result. Each anchor fixes one sign but flips another. Dual-anchor + 3rd eigenvalue → catastrophic basin runaway. |
| Three specific asks | (1) initial-guess derivative profile in production code (is the slope at r=R_MIN explicitly zero or some specific non-zero value?); (2) initial Lagrange multiplier λ value (near 0, near ω, or specific number?); (3) does the production solver use any explicit sign-pinning, OR a continuation method (e.g., homotopy from m_J²=0 → 0.5)? |
| Closing | We're not rushing per his guidance. One specific algorithmic clarification likely unblocks immediately. |

---

## Next steps

### Immediate — hard stop, awaiting Paul's algorithm clarification

Email v7 sent. Three specific asks about solver setup:

1. Initial-guess derivative profile shape
2. Initial Lagrange multiplier value
3. Sign-pinning constraint OR continuation method

### Tomorrow (or after Paul's reply)

Three branches depending on Paul's answer:

| If Paul says... | We do |
| --- | --- |
| "Use specific derivative profile X" | Update initial guess in v7; rerun |
| "λ starts at 0, not 12 or 0.1" | Try v7 with that init value |
| "Use continuation in m_J² from 0 to 0.5" | Build continuation wrapper; this becomes v7.x |
| Doesn't reply within 24-48 hrs | Default to Option C — switch to `scipy.optimize.root` with custom Jacobian |

### Parallel — SABER path stays unblocked

Per cardinal rule: M5 is SABER's primary engineering substrate. With v7
blocked on solver details, M5 work proceeds in foreground. M5.4 substrate
migration is queued and not affected by M6 v7's status. The v6.6 H/Q=1.7112
(0.84% off) result is already sufficient as the "M6 cross-validation" cover
for the OpenWave platform position — full ground state can land later.

---

## Files

| Path | Contents |
| --- | --- |
| `sandbox_v7/m6_v7_4fn_ground_state_bvp.py` | v7 implementation. Forks v6.6; dual H computation (H_full + H_electron); anchor options V/Q/VQ; optional 3rd free eigenvalue (m_J² or λ_bench). Status: 7 variants tested, all wrong basin. |
| `sandbox_v7/m6_v7_4fn_ground_state_bvp_results.json` | Last run's JSON output. |

---

## Question tracker (v7 net change — migrated to `0b_question_tracker.md`)

The full question + hardest-pieces tracker now lives in
**`0b_question_tracker.md`** as the single source of truth across all
sandbox iterations. The v7-era net change is preserved below for
archaeology; **update `0b_question_tracker.md` for any new changes**,
not this file.

| v7 status change | Detail |
| --- | --- |
| RESOLVED post-v7 | **Q28** (quartic interpretation), **Q29** (node count for ground state), **Q30** (Q(0) sign) — all closed by Paul's 2026-05-21 reply. |
| DEMOTED post-v7 | **Q26** (5% gap) → effectively resolved by drop-quartic + Q28; **Q27** (Q_CS_grid vs Q_CS_I disagreement) → excited-mode artifact, not a tunable. |
| NEW post-v7 | **Q31** (initial-guess derivative profile), **Q32** (initial Lagrange multiplier value), **Q33** (sign-pinning constraint vs continuation method) — all sent in email v7 evening of 2026-05-21. Block v7 ground-state convergence. |

Active count entering v8 / Paul's Q31/Q32/Q33 reply: **3 IMMEDIATE + 4 OPEN + 2 DEMOTED = 9 active questions.** See `0b_question_tracker.md` for the full list across all sandbox iterations.
