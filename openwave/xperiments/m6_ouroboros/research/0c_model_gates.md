# M6 / Ouroboros — Production Gates

Three gates must pass before committing M6 to full Taichi production. The
model is a credible scientific candidate; these gates determine whether it is
a credible *engineering* candidate alongside M5.

Last updated: 2026-05-19.

---

## Gate definitions

```text
Gate | Decision question                       | Unblocked by
-----|-----------------------------------------|-------------------------------
G1   | Does the lepton scan reproduce          | Run v4 ourselves (Priority A,
     | ω^2.22 scaling at <5% muon/tau          | ~1 day once Lean ODE confirmed).
     | gaps with the correct Lean ODE?         | Werbos's session shows 1.1%
     |                                         | muon gap at ω=12.78 — we need
     |                                         | to reproduce that ourselves.
G2   | What is the true neutral chaoiton       | Werbos's code (Q12 in email,
     | ground-state m_χ — our 0.508 MeV        | sent 2026-05-18). Our v3 gets
     | or Werbos's 0.003-0.015 MeV?            | 0.508 MeV; his prelim estimate
     |                                         | is 30× lighter. One of us has
     |                                         | a wrong energy functional.
G3   | Does a discrete ω selection             | Werbos reply on Q2 (email) OR
     | mechanism exist that picks              | failed: lepton masses are 3
     | ω = {1, 11, 40.7} from first            | fitted inputs, not predictions.
     | principles?                             | This is Werbos's own "key open
     |                                         | question" (Spectrum §6.1).
```

---

## Pass/fail criteria

```text
Gate | PASS criterion                          | FAIL consequence
-----|-----------------------------------------|-------------------------------
G1   | Muon gap <5%, tau gap <10%, reproduced  | M6 stays in sandbox phase
     | with our own code, not just Werbos's    | indefinitely. Cross-validation
     | bosonization session result             | value is not real if we can't
     |                                         | independently verify the claim.
G2   | m_χ confirmed to within 2× of           | DM prediction in DarkMatterv1
     | Werbos's prelim estimate with           | is at risk. Still worth building
     | same code base / energy functional      | M6 but without the DM headline.
G3   | A rigorous selection mechanism          | Lepton masses are fitted (like
     | (quantization condition) identified     | Standard Model, not a revolution).
     | OR its absence confirmed rigorously     | M6 still worth building as an
     |                                         | alternative substrate; weaker
     |                                         | scientific claim.
```

---

## Gate dependencies

```text
G2 → depends on Werbos's code (waiting reply as of 2026-05-18)
G1 → independent; can run sandbox v4 now with Lean ODE
G3 → parallel; may be resolved by Werbos reply or remain open

Recommended sequence:
  Run G1 first (days) → run G2 once code arrives (days) → G3 (weeks/open)
  Don't wait for G3 before production decision — it may never close analytically.
```

---

## GO / NO-GO decision

```text
Scenario                         | Decision
---------------------------------|------------------------------------------
G1 PASS + G2 PASS                | GO — scaffold M6 in Taichi, parallel to M5.
                                 | Start from Vector(4) × 2 substrate.
                                 | Timeline: post-M5.4 substrate migration.
G1 PASS + G2 FAIL                | CONDITIONAL GO — build M6 but without DM
                                 | headline. Lepton spectrum + nuclear force
                                 | are enough for scientific value.
G1 FAIL                          | HOLD — don't commit engineering time until
                                 | we can reproduce lepton gaps <5% ourselves.
G3 PASS (discrete mechanism)     | Upgrade to STRONG GO at any stage.
G3 FAIL (mechanism absent)       | Lepton masses are fitted inputs.
                                 | M6 is valid as alternative substrate;
                                 | weaker claim than Werbos markets.
```

---

## What M6 offers regardless of gate outcomes

Even in a G1 FAIL scenario, the Ouroboros model provides:

- Maxwell EM recovered exactly (structural fact — no gate needed)
- Charge quantization via Chern-Simons linking (Lean-proved from axioms)
- An independent test of the heat-as-oscillation hypothesis on a different substrate
- A simpler Taichi substrate than M5 (2-function ODE vs matrix M = ODO^T)
- A collaboration anchor with Paul Werbos (NSF X-Labs positioning)

These justify maintaining M6 as a sandbox/research track even if production is deferred.

---

## Gate status log

```text
Date       | Gate | Update
-----------|------|----------------------------------------------------------
2026-05-18 | G1   | PENDING — have correct Lean ODE; haven't run lepton scan
  morning  |      | with it yet. Werbos's bosonization session shows 1.1% muon
           |      | gap at ω=12.78 — promising. Running this is v4 Priority A.
2026-05-18 | G2   | PENDING — code request sent in email to Werbos. Our v3
  morning  |      | gets 0.508 MeV neutral m_χ; his estimate is 0.003-0.015 MeV.
           |      | 30× discrepancy. Likely energy functional difference.
2026-05-18 | G3   | OPEN — three sandbox rounds (v1, v2, v3) and 3 email exchanges
  morning  |      | have not produced a discrete ω selection mechanism. Werbos
           |      | himself flagged this as "the key open question."
-----------|------|----------------------------------------------------------
2026-05-18 | G1   | UPGRADED PRIORITY — Werbos's reply hypothesizes that
  19:39    |      | leptons are the 3 lowest stable modes in Q=1 sector.
           |      | Running v4 Track 1 lepton scan tests this DIRECTLY.
           |      | If our 3 lowest stable ω match leptons within 4-6%,
           |      | G3 is empirically resolved alongside G1.
2026-05-18 | G2   | CODE PROMISED — Werbos: "I'll share the Python source
  19:39    |      | for the Q=0 scan and the calibration run... Give me a
           |      | day or two to get the files to you." ETA: 2026-05-20.
           |      | Will resolve calibration localization + neutral ground
           |      | state + new Q14 (locked vs A=0 description conflict).
2026-05-18 | G3   | RECEIVING PATH FORWARD — Werbos acknowledges "empirical
  19:39    |      | for now"; analytic proof = "major mathematical
           |      | undertaking" deferred. But the v4 Track 1 lepton scan
           |      | provides an EMPIRICAL test: if 3 lowest stable modes =
           |      | leptons, the framework's discrete-spectrum claim has
           |      | numerical backing even without analytic proof. This
           |      | could effectively close G3 for production purposes.
-----------|------|----------------------------------------------------------
2026-05-19 | G1   | ATTEMPTED, BLOCKED. Built sandbox_v4/m6_v4_lepton_scan.py
  PM       |      | with Lean 2-fn ODE; quick ω-sweep found H decreases
           |      | with ω (opposite of Werbos's bosonization). ω=1.0
           |      | marginally fails localization. m_μ/m_e ratio target
           |      | 207, found 0.04. Diagnostic confirmed v3's
           |      | H_CODE_ELECTRON_CALIB=0.494 is wrong (real H ≈ 10.7
           |      | at calibration point). Root cause per Werbos's Q13:
           |      | 2-fn ODE is for NEUTRAL sector only; CHARGED scan
           |      | needs 4-fn (φ,α,ρ,β) toroidal ansatz. G1 cannot
           |      | progress until 4-fn ODE is implemented.
2026-05-19 | G2   | v3 RESULT INVALIDATED. Built sandbox_v4/m6_v4_t2_
  PM       |      | neutral_ground.py with wider B₀ ∈ [1e-5, 1.0] log
           |      | scan + golden-section refinement. NEGATIVE: H is
           |      | monotonic in B₀; no minimum exists. All "localized"
           |      | scan points violate Lean ≤4-node regularity (found
           |      | 13-29 nodes). Linear part is β''+β'/r-β/r²+λβ=0 =
           |      | J_1(√λr) Bessel — oscillatory, not bound. v3's 23
           |      | solutions + m_χ=0.508 MeV are artifacts of permissive
           |      | absolute-tail check. The cited DM result is not
           |      | scientifically defensible until canonical Q=0 ODE is
           |      | settled. Q14 now 3-way: locked (v2 failed) / A=0
           |      | (v3/T2 disproved) / Q_A≈0 (untested).
2026-05-19 | G3   | BLOCKED. Empirical-via-G1 path now blocked alongside
  PM       |      | G1 itself. Analytic-mechanism path still deferred
           |      | per Werbos's own admission. No movement until G1
           |      | unblocks (Werbos's code OR solo 4-fn).
-----------|------|----------------------------------------------------------
2026-05-19 | TRIG | Werbos sent DM paper draft (Dark_Matter_in_Universe
  10:51 AM |      | v0.txt) for ApJ review citing Griesi v3 m_χ=0.508 MeV
           |      | as numerical input to relic-abundance Ω_χh²≈0.1-0.2.
           |      | The cited result is the one T2 invalidated 30min
           |      | earlier. Rodrigo replied same morning with honest
           |      | feedback: artifacts finding, 3-way Q14 conflict,
           |      | 4-fn extraction work-in-progress, code request as
           |      | path to firm up m_χ before submission. Awaiting Paul's
           |      | response.
2026-05-19 | REPLY| Werbos replied (via DeepSeek) with 3 clarifications:
  1:49 PM  |      | (1) m_eff² = m_J² − ω² (ω absorbs into mass term in
           |      | the static benchmark ODE); (2) f(s) v5 ≡ benchmark
           |      | when m_J²=0, λ=4g; use benchmark as general form;
           |      | (3) Toroidal R cancels in H/Q ratio, set R=1.
           |      | Endorsed Q14 → canonical Q=0 is Q_A≈0, Q_J≠0
           |      | (both fields active, EM-neutral). Code arrives "as
           |      | soon as I can" — no committed date.
2026-05-19 | G2   | T6 first attempt at 4-fn benchmark ODE with
  PM       |      | m_eff² substitution. Built m6_v4_4fn.py with
           |      | value BCs V₀=A₀=Q₀=J₀=0.1, derivatives=0.
           |      | RK45 → stiff, LSODA → hangs, RK45+max_step+blowup
           |      | event → works numerically.
           |      | RESULT: ALL 40 (m_J², λ_bench) grid combos BLEW UP.
           |      | Reveals 4-fn ODE is eigenvalue/shooting problem;
           |      | bound state exists only at specific calibrated
           |      | (m_J², λ_bench) which we don't yet have. Q_CS
           |      | also negative in m_eff²>0 cases (winding sign).
           |      | Path forward: T6→A (BVP solver), then T7 (Newton
           |      | shoot), then ask Werbos for specific m_J² + λ_bench
           |      | values if neither converges.
2026-05-19 | G2   | T6→A BVP attempt. Built m6_v4_4fn_bvp.py with
  PM later |      | scipy.solve_bvp + eigenvalue parameter(s). RESULT:
           |      | PARTIAL. 1-eigenvalue mode (m_J² + V_norm)
           |      | converges to a non-trivial (V, Q) bound state, but
           |      | (A, J) collapses to zero spontaneously → Q_CS = 0,
           |      | not Q_CS = 1 target. The eigenvalue m_J² is also
           |      | r_max-dependent (0.09 @ r_max=12; 0.94 @ r_max=24)
           |      | suggesting Bessel-zero artifact, not true bound
           |      | state. 2-eigenvalue mode (m_J², λ_bench + V_norm,
           |      | A_norm) fails universally with singular Jacobian
           |      | across 36 scan combinations. Fundamental issue:
           |      | Q_CS=1 is a topological/integral constraint, not a
           |      | boundary condition — solve_bvp can't enforce it.
           |      | Moving to T7 (shooting with Q_CS=1 as residual).
2026-05-19 | G2   | T7 shooting attempt. Built m6_v4_4fn_shoot.py
  PM later |      | using r_blowup distance as continuous figure of
           |      | merit. Coarse 2D scan (m_J², λ) ∈ [0.5, 8] ×
           |      | [0.1, 10] at V₀=A₀=Q₀=J₀=0.1: max r_blowup = 6.44,
           |      | well below r_max=15 → no bound state in scan
           |      | window. m_J²=1.0 row flat (m_eff²=0 degenerate).
           |      | Amplitude shoot: r_blowup monotone decreasing in
           |      | A_0; no resonance/sweet spot.
           |      | DEFINITIVE: NO bound state in symmetric ansatz
           |      | V₀=A₀=Q₀=J₀ regime anywhere scanned. Werbos uses
           |      | asymmetric initial values OR different param
           |      | region we haven't covered.
           |      | Next: email Werbos with sharper specific ask —
           |      | what are the canonical (m_J², λ_bench, V₀, A₀,
           |      | Q₀, J₀) values at his electron calibration?
2026-05-19 | TRIG | Email sent ~5:00 PM with T6/T6→A/T7 findings +
  5:00 PM  |      | request for canonical (m_J², λ_bench, V₀-J₀)
           |      | values + shooting algorithm.
2026-05-19 | REPLY| Werbos reply (via DeepSeek): asymmetric helicity
  4:21 PM  |      | (NOTE: arrived BEFORE the 5 PM email — these
  (cross-  |      | crossed in transit). V₀=Q₀=+0.1, A₀=J₀=-0.1
   crossed)|      | gives Q_CS=1. m_J²≈0.5, λ_bench=1.0,
           |      | m_eff²=-0.5 (negative is correct, time-periodic).
           |      | "Shooting with decay rate at infinity as target"
           |      | — algorithm not detailed. Also v2 paper draft
           |      | sent for review.
2026-05-19 | G2   | T8: Re-ran 4fn IVP / shoot / BVP scripts with
  PM later |      | Werbos's asymmetric helicity prescription.
           |      | RESULT: helicity is NECESSARY but NOT SUFFICIENT.
           |      | IVP at (m_J²=0.5, λ=1): Q_CS goes from 0
           |      | (symmetric) → 4483 (asymmetric) — confirms
           |      | helicity prescription. But still not Q_CS=1.
           |      | Amplitude shoot: ~12% better than symmetric but
           |      | still monotone. BVP with 1-eigenvalue +
           |      | asymmetric init guess: collapses to (V,Q)
           |      | sub bound state at m_J²=5.51 regardless of init.
           |      | A,J drift to 0 — no BC holds them.
           |      | Q_CS=1 chaoiton requires: (a) topological
           |      | constraint Q_CS=1 as Newton residual, or
           |      | (b) tabulated init profile in right basin, or
           |      | (c) Werbos's actual code.
           |      | Review of v2 paper sent ~4:30 PM same day.
           |      | Awaiting Werbos's follow-up.
-----------|------|----------------------------------------------------------
2026-05-20 | G2   | T9: Newton-residual two-stage shooter.
  PM       |      | Built sandbox_v4/m6_v4_4fn_newton.py with
           |      | scipy.optimize.differential_evolution (Stage 1
           |      | global) + Nelder-Mead (Stage 2 local) over 6
           |      | variables (m_J², λ_bench, |V₀|, |A₀|, |Q₀|, |J₀|)
           |      | with helicity signs locked. Pre-scan: 2D grid
           |      | at Werbos's exact |amps|=0.1, 324 points
           |      | (m_J² ∈ [0.05, 10] log, λ ∈ [0.01, 10] log,
           |      | r_max=15). RESULT: 0/324 reach r_max with peak<5.
           |      | Max r_reached = 7.241 with peak=95 (catastrophic
           |      | blowup). Wider scan 576 points (m_J² ∈ [0.01, 50],
           |      | λ ∈ [0.001, 50], r_max=30): 0/576 reach r_max.
           |      | Max r_reached = 8.06 with peak=98. Stage 1 6-var
           |      | DE search (6528 evals): r_reached=8.4, Q_CS=58662
           |      | — delayed-blow-up regime found, NOT bound state.
           |      | DEFINITIVE NEGATIVE: forward-IVP from value-BC
           |      | origin does NOT reach the Q_CS=1 chaoiton in any
           |      | parameter region we have explored. Across 6
           |      | independent attempts (T6/T6→A 1-eig/T6→A 2-eig/
           |      | T7/T8/T9), no bound state found. Werbos's
           |      | shooting algorithm is not forward-IVP — needs
           |      | algorithmic clarification (most likely backward
           |      | integration from K_0 asymptotic, OR BVP with
           |      | Q_CS Lagrange multiplier, OR specific init
           |      | profile shape, NOT just BC values). G2 path
           |      | depends on Werbos's reply OR major BVP rework.
2026-05-20 | TRIG | New paper arrived (12:53 PM): "The Neutral
  PM       |      | Chaoiton: A Dark Matter Candidate from the
           |      | Ouroboros Lagrangian" — ApJ-targeted compact
           |      | dark matter paper, distinct from the broader
           |      | Dark_Matter_in_Universe_v4 (Zenodo 20298669).
           |      | Cites both prior depositions (v8 LoE @ 20313063
           |      | and v4 DM @ 20298669). Acknowledgments thank
           |      | Rodrigo Griesi for "independent numerical
           |      | reproduction of the electron calibration and
           |      | methodological dialogue on the asymmetric
           |      | helicity structure required for Q_CS = 1."
           |      | Reference [14]: "Griesi, R., & Anthropic AI
           |      | 2026, in preparation (ground state computation)."
           |      | Paul asked when to push a "TRULY new" Zenodo
           |      | version of either paper — holds off pending
           |      | our T9 / ground-state numbers. With T9 negative,
           |      | the m_χ / m_J / σ/m numbers depend on resolving
           |      | the shooting-algorithm question first.
-----------|------|----------------------------------------------------------
2026-05-20 | TRIG | Duda public reply on models-of-particles list
  ~1:15 PM |      | (responding to v8 LoE Zenodo 20313063 announce):
           |      | Four substantive technical critiques:
           |      | (1) G_μν has no microscopic definition — the G
           |      |     curl tensor in -G^μν G_μν is just symbolically
           |      |     given; J itself is a primary undefined field.
           |      |     Same "single-field ontology" objection Duda
           |      |     raised 2026-05-13 ("Fundamental model should
           |      |     have a single field").
           |      | (2) f(J_μ J^μ) is unspecified in the LoE paper.
           |      |     "is it square? Higgs?" The Numerical Benchmark
           |      |     sub-document does pin f(s) = ½ m_J² s + ¼ λ s²,
           |      |     but the standalone LoE paper leaves it open.
           |      | (3) The topological charge Q = (1/4π²)∫ε^μνρσ F G dx
           |      |     "clearly G is crucial here, you don't specify".
           |      |     For an electron, the paper asserts H/Q=1.6969 but
           |      |     doesn't show the field configuration (ansatz +
           |      |     energy minimization). T6-T9 confirm empirically
           |      |     this construction is not available via forward IVP.
           |      | (4) For Coulomb V~1/r between two charges, "you need
           |      |     to integrate Hamiltonian for two topological
           |      |     charges in distance" — Werbos asserts Yukawa
           |      |     far-field but doesn't derive two-chaoiton
           |      |     interaction. v5 §6 has the form V~-C·exp(-m_J r)/r
           |      |     but as ansatz, not derived from H[Φ_1, Φ_2].
           |      | Closing: "I still see LLM-generated word salad just
           |      | to look nice for the user" — tonally harsh but
           |      | technically the same construction critique. Consistent
           |      | with Duda's 2026-05-08 "AI slop on zenodo" warning.
           |      | Our internal read: (1) and (4) are structural; (2) is
           |      | editorial gap fixable in next LoE revision; (3) is
           |      | what T9 just demonstrated empirically — Werbos has
           |      | the algorithm offline but had not described it in
           |      | the canonical paper. Now-resolved-internally with
           |      | Paul's 2:00 PM algorithm reply (see below); still
           |      | unresolved in any public-facing Werbos document.
-----------|------|----------------------------------------------------------
2026-05-20 | REPLY| Werbos algorithm clarification (2:00 PM via
  ~2:00 PM |      | DeepSeek). Response to our T9 definitive-negative
           |      | + algorithmic ask. Key acknowledgments and content:
           |      | - "forward shooting from the origin with decay
           |      |   conditions at infinity does NOT work for this
           |      |   system. That is consistent with our experience."
           |      |   → confirms our T6-T9 negative result.
           |      | - Method: collocation / finite-difference BVP via
           |      |   `scipy.integrate.solve_bvp` or `scipy.optimize.root
           |      |   (method='lm')`. Nonlinear eigenvalue problem with
           |      |   Q_CS=1 as INTEGRAL CONSTRAINT (not BC).
           |      | - Two free eigenvalues: ω and Lagrange multiplier λ
           |      |   (from H' = H - λ·Q). Origin BCs are derivative=0
           |      |   only (V'(0)=A'(0)=Q'(0)=J'(0)=0); values V(0),
           |      |   A(0), Q(0), J(0) FREE.
           |      | - R_max BCs are ROBIN (V'+k·V=0 etc., not zero
           |      |   Dirichlet) matching K_0(κr) exponential decay.
           |      |   k = √(ω²-m_J²), initial approximation k = ω.
           |      | - Initial profile guess: V=+0.1·exp(-r), A=-0.1·exp(-r),
           |      |   Q=+0.1·exp(-r), J=-0.1·exp(-r). Werbos's "V(0)=0.1"
           |      |   was an init-profile amplitude, NOT a value BC.
           |      | - Post-convergence: Gelfand-Fomin conjugate-point
           |      |   test on second variation (per Lean theorem).
           |      | - DeepSeek offered to write Python; Paul deferred
           |      |   to us ("I suspect you want Claude involved").
           |      | NEXT: T10 implementation — solve_bvp with Paul's
           |      | recipe. Hold M5 return briefly, give M6 one more
           |      | focused attempt with the correct method.
-----------|------|----------------------------------------------------------
2026-05-20 | NEXT | T10 plan: m6_v4_4fn_lambda_bvp.py
  PM       |      | - State: (V, V', A, A', Q, Q', J, J', I) size 9
           |      |   (I = accumulated Q_CS integral for closure)
           |      | - Free params: ω, λ (Lagrange multiplier)
           |      | - ODEs Paul-as-written (no explicit λ-terms in A, J;
           |      |   if singular, retry with λ-corrections from
           |      |   δQ_CS/δA = ∂_r J formulation)
           |      | - BCs: V'=A'=Q'=J'=0 at r=0 (5 incl. I(0)=0);
           |      |   Robin V'+k·V=0 etc. at R_max (4); I(R_max)=1 (1).
           |      |   Add normalization (e.g. V(0)=0.1) if 1 BC short.
           |      | - Init: exponential decay profiles, ω=λ=1
           |      | - Grid: 200 points non-uniform on [R_MIN, 20-30]
           |      | - Acceptance: |H/Q-1.6969|<0.001, |Q_CS-1|<0.01,
           |      |   tail<0.05, nodes≤4, no blowup. Same T9 criterion.
2026-05-20 | TRIG | Question tracker reset for v5 (see 0b_sandbox_v5.md
  PM       |      | tail). Status counts entering v6:
           |      |   3 IMMEDIATE   Q22 (Q_CS normalization),
           |      |                 Q23 (H functional kinetic factor),
           |      |                 Q24 (sample converged profile)
           |      |                 → all 3 in Werbos email v4 (sent PM)
           |      |   1 ACTIVE     Q20 (Duda critique #3, half-addressed
           |      |                 by v5; v6 closes if Q22-24 land)
           |      |   5 OPEN       Q2, Q3, Q6, Q19 (Duda #2 editorial),
           |      |                 Q21 (Duda #4, future v7+)
           |      | Resolved by v5: Q1, Q9-Q15, Q16, Q17, Q25
           |      |   (Hopf invariant proof complete; charge quantization
           |      |    is now a theorem of differential topology).
           |      | Archived: Q4 (= Duda #1, single vs two-field
           |      |   ontology — unfalsifiable preference; if math
           |      |   matches observation, two fields is just a
           |      |   description, not a theory-killer).
-----------|------|----------------------------------------------------------
2026-05-20 | G2   | T10 PARTIAL SUCCESS. Built m6_v4_4fn_lambda_bvp.py.
  PM       |      | scipy.integrate.solve_bvp with:
           |      | - 9-state (8 fields + I = accumulated Q_CS integral)
           |      | - 2 free params (ω, λ_LM); Lagrange-multiplier
           |      |   corrections derived from H' = H - λ·Q_CS:
           |      |     ΔA = J + λ·(J + 2r·J')
           |      |     ΔJ = [unconstr RHS] - λ·(A + 2r·A')
           |      | - V(R_MIN) = 0.1 anti-collapse normalization
           |      | - Initial profile: exp(-r) shapes with A,J non-
           |      |   proportional (J on exp(-1.5r) to break Q_CS≡0
           |      |   degeneracy of symmetric init)
           |      | - 50000 max nodes, tol=5e-3, r_max=12.
           |      | RESULT: solve_bvp.status = 0 (CONVERGED), 28k final
           |      | nodes, max residual 1.3e-4. First-ever clean
           |      | convergence to a Q_CS=1 chaoiton via Werbos's actual
           |      | method.
           |      | Converged values:
           |      |   ω = 1.047    vs Werbos 1.0   (4.7% over)
           |      |   m_eff² = -0.596 vs Werbos -0.5 (19% over)
           |      |   λ_LM = -1.212  (first time pinned)
           |      |   Q_CS = 1.000   (exact via integral constraint)
           |      |   H/Q_from-I = 52.640
           |      |   H/Q_from-grid = 52.645 (0.01% agreement)
           |      | Did NOT match:
           |      |   H/Q_CS = 52.64 vs Werbos 1.6969 (31× off)
           |      |   A has 17 nodes (excited mode, not ground state)
           |      |   tail = 0.17 (slow decay)
           |      | V_norm scan {0.5, 0.3, 0.2, 0.15, 0.1, 0.05, 0.03,
           |      | 0.01}: H/Q_CS bottoms at 52.64 at V_norm=0.1; never
           |      | reaches 1.6969 anywhere. 31× ratio is STRUCTURAL,
           |      | most likely a Q_CS or H normalization convention
           |      | mismatch (Paul uses Chern-Simons (1/4π²)·∫ε A ∂J;
           |      | we use 2π·∫r·(A·J'-J·A')dr radial toroidal form).
           |      | NEXT: targeted Werbos email v4 — three specific
           |      | normalization questions + ask for sample converged
           |      | profile to anchor ground-state basin selection.
           |      | Major shift from T9 negative: forward-IVP-wrong-tool
           |      | confirmed solved; only normalization/profile match
           |      | remaining. M6 chaoiton existence empirically
           |      | demonstrated for the first time.
```

---

## Werbos reply update (2026-05-18 19:39)

The first round of gate-relevant clarification arrived. Key impacts:

- **G1**: Sharper test designed — v4 lepton scan can falsify or confirm
  Werbos's "lowest 3 stable modes = leptons" hypothesis directly.
- **G2**: Code coming in 1-2 days. Will close calibration + verify neutral
  ground state in one round.
- **G3**: Analytic proof deferred indefinitely (Werbos's own admission).
  Empirical G1 result becomes the effective gate.

**Refined GO/NO-GO**: with G1 + G2 both within reach this week, the
production decision can plausibly land by 2026-05-25.

---

## Status revision (2026-05-19, post-v4 first runs)

The 2026-05-18 "production decision by 2026-05-25" timeline is now too
optimistic. Both G1 and G2 hit unexpected blockers in v4 first runs:

```text
Gate | Was (2026-05-18)               | Is (2026-05-19)
-----|--------------------------------|--------------------------------
G1   | READY — run lepton scan w/     | BLOCKED — 2-fn ODE wrong tool
     | Lean 2-fn ODE; expect muon@    | for charged sector. Need 4-fn
     | 1.1% gap (Werbos bosonization) | benchmark ODE first.
G2   | CODE COMING — Werbos's source  | v3 RESULT INVALIDATED. T2 showed
     | will close calibration + m_χ   | 23 solutions are artifacts.
     | ground state in one round      | Existence in 2-fn form uncertain.
G3   | EMPIRICAL VIA G1               | BLOCKED with G1.
```

**Why the optimism was wrong**: we assumed the Lean 2-fn ODE was the
canonical form for both sectors. Per Werbos's Q13 reply, it is the
reduced form for the NEUTRAL sector only; the charged calibration +
lepton scan use the 4-function (φ, α, ρ, β) toroidal ansatz from his
pre-Lean code. The two ODEs are structurally different (toroidal l=0
vs vector l=1), not reductions of each other.

**Revised timeline**: production decision likely 2026-05-23 to
2026-05-27 if Werbos's code arrives this week; 2026-05-25 to
2026-05-30 if we have to implement the 4-fn ODE solo.

**Two unblock paths**:

```text
Path | Trigger                          | Resolves
-----|----------------------------------|---------------------------
A    | Werbos sends his Python source   | G1 + G2 + Q12 + Q14
     | (in response to our 2026-05-19   | in one shot
     | ask)                             |
B    | Solo: implement 4-fn benchmark   | G1 + G2 paths (with risk
     | ODE in sandbox_v4/, calibrate    | from 3 open gaps:
     | to H/Q=1.6969, then run scans    | ω-in-static, f(s) form, R)
```

The 2026-05-19 reply to Werbos was framed around Path A but acknowledges
Path B is the fallback if his code doesn't materialize.

---

## Status revision 2 (2026-05-19 PM, post-Werbos reply + T6)

Werbos resolved the 3 implementation gaps (m_eff² formula, f(s) mapping,
R cancels) but T6 first run revealed the deeper problem: the 4-fn ODE
is an **eigenvalue/shooting problem**, not solvable by plain IVP from
generic initial conditions. Scan of (m_J², λ_bench) at V₀=A₀=Q₀=J₀=0.1
across 40 grid points: ALL BLEW UP.

```text
Implementation reality       | Status
-----------------------------|---------------------------------------
Gaps from Werbos's reply     | RESOLVED in formula. Need ONE more
                             | piece: the specific m_J², λ_bench
                             | values at his calibrated point.
4-fn ODE is eigenvalue-like  | Generic IVP doesn't find bound state.
                             | Need BVP (solve_bvp with decay BC at
                             | r_max) or Newton-style shooting.
Path B has new sub-paths     | T6→A: BVP solver (2-4 hr, medium risk)
                             | T7: Newton shoot on (m_J², λ_bench)
                             | T8: ask Werbos for specific values
                             | (1-line email, lowest risk)
```

**Updated path sequence** (Rodrigo direction, 2026-05-19 PM):

```text
Step | Action                              | If succeeds        | If fails
-----|-------------------------------------|--------------------|-------------
1    | DONE — document T6 negative result  | proceed to (2)     | —
2    | T6→A: implement BVP solver          | use to calibrate   | → (3)
3    | T7: Newton shooting on              | use to calibrate   | → (4)
     | (m_J², λ_bench)                     |                    |
4    | Ask Werbos for specific calibration | resolves           | M6
     | values (m_J², λ_bench at canonical  | immediately        | parked
     | point)                              |                    |
```

The escalation principle: prove what we tried before going back to him.
