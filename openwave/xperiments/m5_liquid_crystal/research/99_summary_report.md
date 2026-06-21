# M5 Summary Report — Results & Canonical Implementation

**Purpose:** the results-of-record for the M5 Liquid-Crystal program (OpenWave). A narrative report of the 3+1D time-crystal field test, then every headline test as `topic | verdict + one-sentence finding + source script`, organized into the three validation bodies (**Liquid Crystal**, **Time Crystal**, **Zitterbewegung-clock existence**), then the **canonical implementation** (what works — the build that produced these results). The per-phase narrative lives in [`0b_M5_roadmap.md`](0b_M5_roadmap.md) and the hardest-pieces board in [`0b_question_tracker.md`](0b_question_tracker.md).

**Reading the verdict column:** the first line of each descriptive cell is the result — ✅ PASS / POSITIVE, ❌ NEGATIVE, ⚠️ NEGATIVE-INFORMATIVE or caveat. Scope throughout (unless noted): single defect, natural/lattice units, Duda's Eq.18-class matrix dynamics `M = ODOᵀ`, sandbox grids 24³/48³, f32 GPU with f64 numpy cross-checks.

**Orientation:** "M5" is OpenWave's implementation of Duda's liquid-crystal model; M5.7 / M5.8 are its phases (the resonance hunt; the Zitterbewegung clock). "N-1…N-6e" are the clock-measurement sub-experiments; the remaining labels (NG-, G-, Track C) index the project tracker and are not load-bearing for the results.

**Last updated:** 2026-06-08 (the full N-1…N-6e ZBW program, 2026-06-07).

---

## Report — the 3+1D time-crystal field test

Following exchange with Jarek Duda, this is the finished 3+1D field test of the time-crystal mechanism from arXiv:2501.04036. It is a results report — just what came out. Everything stated here is reproducible from the OpenWave repository (scripts plus per-experiment notes).

**Setup.** The field is the matrix `M = O D Oᵀ` with `D = diag(1, δ, 0, g)`, δ = 0.3, time as the fourth index. The dynamics is the Eq.18 signed Hamiltonian with the Minkowski (α,3) block. The only stable integrator found is a constrained spectral-projection scheme that keeps the positive-inertia directions per voxel and projects the momentum onto that subspace every step — every cheap positive-inertia kinetic tried has slow growing modes. All runs are headless on 24³ and 48³ grids (with a 63³ check on the runaway onset), f32 on GPU with f64 numpy cross-checks (they agree to about 1% everywhere checked), and a full-Euclidean twin ran as a kill-control at every stage.

**The quadratic action does not saturate at 3+1D.** The bare `−ΣF²` action runs away at a fixed phase, about two clock periods, and the onset is dt-invariant (it scales exactly with dt under halving), identical in f64, f32 and Metal, grid-independent across 24³ and 63³, and independent of the potential and the momentum clamp. That the runaway onsets at a fixed phase rather than a fixed time indicates a phase-locked (parametric) instability, not a fixed-clock numerical limit. This is the field-level version of the toy model's own requirement for a βR⁴ floor.

**A saturating quartic restores bounded dynamics.** Adding `V = u + β u²` on the signed spatial-curvature density `u = 2 Σ_(i<j) ⟨F_ij,F_ij⟩_s`, with β fixed by `2β|u_min| ≈ 3`, bounds the runaway: H breathes over the full horizon, the floor-bounce is real, and it is dt-converged and f64-confirmed.

**The signature is what drives it, not the nonlinearity.** With the same seed and kick, the Euclidean twin's clock decays while the Minkowski-plus-quartic clock grows by about 3×. The (α,3) block is doing the work.

**The saturated state self-starts.** A configuration damped to rest, restarted with momentum exactly zero, regrows kinetic energy from zero (to 5.76 by τ = 4000, H conserved to 0.5%), reproduced at half dt to four significant figures — a field starting at exactly zero momentum develops motion at conserved energy. The collective-coordinate analysis behind this finds the dressed minimum to be un-sittable — energetically minimal but carrying ghost breathing inertia (K_bb < 0) — which reads as the 3+1D face of the negative-energy auto-propulsion in Fig 10 of the paper.

**Its frequency is a property of the state.** The breathing fundamental is the same — within one FFT bin — whether the run is kicked, exactly unkicked, or seeded with small random jitter. One correction belongs here: an earlier reading claimed a strictly periodic ω₀ = 0.262 with an exact harmonic comb; a window-length sweep showed that was an FFT-window artifact — the peak moved with the window length and the comb was just bin arithmetic. The state carries a resolved fundamental plus a 2ω harmonic over a broadband, aperiodic background — a coherent comb, not a strictly periodic line (the background is classified below as low-dimensional chaos).

**The defect holds, and the holding is resolution-robust.** This is the M5.7 dispersal question asked again under the quartic with full backreaction. The structural alignment decays slower at 48³ than at 24³ at every matched time, settling toward a plateau well above the random baseline, and the result does not depend on the masking. This is the reverse of the M5.7 free-defect behaviour, where the wash-out strengthened with resolution.

**The frequency belongs to the core, not the dressing.** Across a family of states spanning a 2.6× range in rest energy (varying the dressing width), the fundamental is constant to within one bin — the exponent of ω against energy is about 0.03, not 1. So the naive ω ∝ mass law does not hold for the dressing energy; the frequency is set by the core, which the V = 0 frame cannot vary because it is scale-free.

**The absolute number, and where it falls short.** Under two explicit postulates — the rest energy mapped to 0.511 MeV and the lattice action unit mapped to ℏ — the clock runs at about 5.5×10¹⁹ rad/s, roughly 28× below the electron Zitterbewegung 2mc²/ℏ. Because the frequency is rigid against energy, that gap is structural rather than something closable by energy bookkeeping; it points at physics this V = 0 sandbox does not contain — the V-on/Faber-r₀ core, the faithful kinetic, or the action-to-ℏ postulate itself. This is reported as the first measurement, not as a pass or fail.

**Two more readings.** A chaos battery classifies the saturated state as a molten clock — a persistent coherent fundamental on low-dimensional chaotic dressing whose intensity grows with excitation — and the deep-settled cold state reads near-regular, so the clock appears to tick cleanly toward the ground state and melt as it heats. Full coldness is forbidden by the same un-sittability that makes it self-start. Separately, the rotation Noether charge of the clock kick is zero: the Θ-twist cancels over the hedgehog sphere, so the breathing channel carries no net frame angular momentum at this scale — though box torque on the 24³ grid bounds that measurement, so it is not read as a spin-½ statement either way.

**What is not claimed.** No particle-stability result; the state is a molten clock, not a strict single-line oscillator. The absolute frequency is 28× off. There is no intrinsic spin above the box-torque floor. It is a single defect in natural units, and the explicit stepper stiffens near the quartic floor so the very-late-time horizons are numerically limited.

**Standing.** The repository reproduces everything above. The 2+1D moving-defect (pilot-wave) test is the natural next step.

The detailed per-test evidence — every claim above plus the Liquid-Crystal substrate validations and the canonical build — follows.

---

## 1. Liquid-Crystal validations (the substrate is faithful)

These establish that the matrix substrate reproduces the established physics Duda's framework requires — charge, forces, EM, geometric mass. All re-verified on the production 4×4 engine after the M5.8.1 promotion (the reproduction of the original 3×3 numbers is itself proof the promotion is physics-preserving).

| Topic | Result · finding · source |
| --- | --- |
| Topological charge | ✅ POSITIVE<br>Charge is an integer, additive winding of the director: Q = ±0.996 for single hedgehogs, ±1 around each core of a pair and 0.000 on the enclosing sphere — quantization and additivity both emergent, not imposed.<br>`m5_8_1_topo_charge.py` |
| Coulomb force (1/d) | ✅ POSITIVE<br>Opposite-sign defects attract with a 1/d law fit at R² = 0.9781 (b < 0), matching the M5.1 reference (0.978) digit-for-digit on the matrix substrate.<br>`m5_4_coulomb_matrix.py` |
| Coulomb — Duda page-18 cross-check | ✅ POSITIVE<br>Duda's own Mathematica page-18 configuration reproduces an attractive 1/d fit at R² = 0.9959 — an independent geometry confirming the force law.<br>`m5_4_coulomb_page18.py` |
| Maxwell from tilts (hydro↔EM) | ✅ PASS<br>Defining E = ∇·n̂ and B = ∇×n̂ recovers the full Maxwell set — ∇·B = 0, the Gauss charge identity, Faraday, and the Lorentz force — to machine precision via the hydro↔EM dictionary.<br>`m5_6_4a_hydro_em.py` |
| Maxwell from tilts (Faber curvature) | ✅ PASS<br>The Faber field strength R = Γ×Γ is a closed 2-form (Maurer–Cartan ⇒ homogeneous Maxwell), giving abelian Coulomb ‖R‖ ~ 1/r² far-field with a running-coupling onset at the core radius r₀.<br>`m5_6_4b_faber_curvature_em.py` |
| KG mass is GEOMETRIC | ✅ POSITIVE<br>The Klein-Gordon mass emerges from minimal coupling to the hedgehog connection Â (Fig.9 reproduced) — the explicit mass term cancels exactly, so mass is geometry, not an added V_ψ.<br>`m5_6_1_kg_operator_check.py` |
| Biaxial hedgehog sources its own twist | ✅ POSITIVE<br>The biaxial defect dynamically generates C_μν ≠ 0 with a ‖C‖ ~ 1/r² restoring-mass profile — it cannot sit static, which is the seed of the M5.8 clock.<br>`m5_6_2a_biaxial_hedgehog.py` |
| Faber mass calibration (E ∝ 1/r₀) | ✅ POSITIVE<br>Porting Faber's MTF regularization pins the defect mass to the core radius, E·r₀ = const (CV = 0.0%), giving the physical knob that maps r₀ = 2.2132 fm → 0.511 MeV — the M5.9 calibration handle.<br>`m5_6_3a_faber_*` / `m5_6_3b_*` |
| Eq.18 matrix evolution conserves energy | ✅ PASS<br>The Duda Eq.18 leapfrog conserves energy: the drift falls 2.15% → 1.13% → 0.03% as dt → 0 (O(dt²)), confirming the production EOM conserves the Hamiltonian to integrator order.<br>`m5_5_4_matrix_evolution_check.py` |
| LdG potential confinement | ✅ PASS<br>The Landau–de Gennes V(M) confines the amplitude Tr(M²) ~3.3× across k ∈ [0.5, 25] with no blow-up — the well is stable across two decades of stiffness.<br>`m5_6_5c_prod_scale.py` |
| Free 3D defect disperses (M5.7) | ❌ NEGATIVE (expected — Derrick)<br>Both a seeded l=1 perturbation and the defect's own intrinsic oscillation disperse their orientation energy in pure 3D (localization washed out 32³→48³) — confirming a static/oscillatory 3D particle is forbidden and the stable object must be 4D.<br>`m5_7_1_l1_resonance_seed.py` · `m5_7_2_intrinsic_oscillation.py` |
| Driven defect sustains an (A,ω) excess (M5.7.3) | ✅ POSITIVE<br>A continuous EM-like drive holds a bounded, frequency-selective (A,ω) excess at ~3× the free baseline at the resonant mode (N=48-confirmed) — a resonant drive holds the defect's (A,ω) state in a maintained excess.<br>`m5_7_3_driven_oscillation.py` |

---

## 2. Time-Crystal validations (the mechanism works)

These are the M5.8 core: a time-periodic resonance — bounded, self-starting, signature-driven — exists in 3+1D. The arc is "quadratic fails → quartic saturates → the saturated state is a self-starting attractor."

| Topic | Result · finding · source |
| --- | --- |
| 1D toy-model clock (Duda 2501.04036) | ✅ POSITIVE<br>The 1+1D time crystal reproduces to 4–5 digits (ω\* = √(70/61) = 1.0712, E* < E(0)) with machine-precision energy conservation, and the GHOST finding — the ψ-sector is negative-kinetic ⇒ a constrained integrator is mandatory — set the rule for everything 4D.<br>`m5_8_0b_*` · `m5_8_0c0d_propulsion_robustness.py` |
| 4D Minkowski clock fuel is real | ✅ POSITIVE<br>The exact 4D Hamiltonian gives E(ω,b) = A(b) + ω²·C(b) with C(b) < 0 for every clock-plane × boost-axis once the time axis is dressed — and exactly zero undressed (the M5.7 null in functional form); the Euclidean flip kills it, so the signature IS the mechanism.<br>`m5_8_2a_4d_hamiltonian.py` |
| The quadratic action does NOT saturate | ❌ NEGATIVE (decisive)<br>The bare `−ΣF²` action runs away at a fixed phase τ ≈ 2 clock periods — dt-INVARIANT (ratio exactly 2.0 under dt-halving), f64 ≡ f32 ≡ Metal, 24³ ≡ 63³, V/clamp-independent — the field-level twin of the 1D toy's own requirement for a βR⁴ floor.<br>`m5_8_2cb_taichi_constrained.py` |
| The saturating quartic restores bounded dynamics | ✅ POSITIVE<br>`V_u = u + β·u²` on the signed spatial-curvature density (β = 1.558, data-anchored) bounds the runaway — H breathes over the full horizon (≥8 fundamental periods), floor-bounce real, dt-converged, f64-confirmed.<br>`m5_8_2d_quartic_saturation.py` |
| The signature IS the engine (kill-control) | ✅ POSITIVE<br>Same seed + kick: the Euclidean twin's clock DECAYS (1.4e-2 → 6.7e-3) while the Minkowski+quartic clock GROWS 3× (7.6e-3 → 2.4e-2) — the Minkowski signature, not generic nonlinearity, sustains the time crystal.<br>`m5_8_2e_invariant_matrix.py` (euclid mode) |
| Static reductions all closed (Track C) | ❌ NEGATIVE (decisive, ×3)<br>The global-clock BVP has no interior ghost minimum (C1), the static clock-phase twist is pure cost (C2-B), and the dressed minimum is UN-SITTABLE — K_bb(a*) = −67.6 < 0 with U″ > 0 (C3) — so the breathing state is irreducibly time-dependent: ḃ ≠ 0 IS the bounce.<br>`m5_8_2f_breathing_bvp.py` · `m5_8_2f2_localized_clock.py` · `m5_8_2f3_breather_orbit.py` |
| Spontaneity confirmed at field level | ✅ POSITIVE<br>A damped-settled configuration restarted with P = 0 EXACT regrows kinetic energy 0 → 5.76 by τ = 4000, H conserved to 0.5%, dt/2-converged to 4 significant digits — from exactly zero momentum the field develops motion at conserved H: the clock self-starts.<br>`m5_8_2g_spontaneity.py` |
| The invariant matrix is COMPLETE | ✅ PASS (all candidates classified)<br>Signed-u is the working invariant; Skyrme `β_E·u_E²` also saturates but damps the clock 10× (sign-blind); the M6-style amplitude quartic A4 is geometrically pinned (no lever arm — the fuel is curvature-sector); the covariant `𝒮+β𝒮²` is deferred-with-reason (cubic Legendre, static sector ≡ the validated quartic).<br>`m5_8_2e_invariant_matrix.py` · `m5_8_2l_invariant_completion.py` |
| The defect HOLDS, resolution-robust (N-2) | ✅ POSITIVE (G-2c-1)<br>Under the 4D quartic stack with full backreaction the free defect's structural alignment decays SLOWER at 48³ than 24³ at every matched t — the anti-M5.7 signature (its wash-out strengthened with resolution; this reverses) — settling to a plateau 0.8–0.88 ≫ the 0.5 random baseline, mask-insensitive.<br>`sandbox_vn/m5_8_2i_dispersal_gate.py` |

---

## 3. Zitterbewegung-clock existence + identification (N-1…N-6e)

The split that matters: **existence** (rows 1–6, all ✅ — a self-starting clock with a frequency that is a property of the state) vs **identification with the electron** (rows 7–10 — the calibration program, where the gaps are). "Does the ZBW clock exist?" — yes. "Is it the electron's?" — not yet, and the gap is localized.

| Topic | Result · finding · source |
| --- | --- |
| ω is an ATTRACTOR (N-1) | ✅ POSITIVE<br>Kicked, exactly-unkicked, and jittered starts all settle to the same breathing fundamental ω₁ ≈ 1.1 + a 2ω₁ harmonic (detrend-stable across filter windows) — the frequency is a property of the STATE, not the preparation.<br>`sandbox_vn/m5_8_2h_omega_attractor.py` |
| The "0.262 comb" was an FFT artifact (N-1) | ⚠️ NEGATIVE-INFORMATIVE (self-caught)<br>A window-length sweep showed the earlier "strictly periodic ω₀ = 0.262 + exact comb" was an FFT-window artifact (the peak moves with window length; the comb was bin arithmetic) — the real state is a resolved fundamental over a broadband, aperiodic background (classified below as low-dimensional chaos, N-6d).<br>`sandbox_vn/m5_8_2h_omega_attractor.py` (ref mode) |
| The breather is NOT a disclination-rod mode (N-1) | ✅ POSITIVE<br>The motion is near-isotropic (anisotropy A_z ≈ 1.0–1.3 vs the static rod baseline 5.58) and core-centered (3× core-shell concentration early), with no excess on the z-axis rod — the clock lives on the defect core, not the line.<br>`sandbox_vn/m5_8_2h_omega_attractor.py` (W4) |
| First ZBW ratio recorded (N-3) | ✅ MEASURED (no pass/fail)<br>ω₁/(2H_rest) = 0.0326–0.0343 across arms (5.4% spread, start-independent) in lattice-natural units — the first measurement of the ratio for the saturated breather, with the factor-2 apolar-doubling bookkeeping respected (no double-correction).<br>`sandbox_vn/m5_8_2j_zbw_ratio.py` |
| ω is RIGID across a mass family (N-6a) | ✅ POSITIVE (a property of the core)<br>Across a knob-gated dressing-width family spanning 2.6× rest / 3.9× state energy, ω₁ is constant within one FFT bin (exponent ω ∝ H^0.03 vs the naive law's 1) — the frequency belongs to the defect core, which the scale-free V=0 frame cannot vary.<br>`sandbox_vn/m5_8_2m_zbw_law.py` |
| Matter + light share a radial cone (N-4) | ✅ POSITIVE<br>The twist (KG) and both tilt (EM) channels hit the same maximal radial signal speed c(r̂)/2 = 1.0000 at N=48 and 64 — matter and light propagate on a common radial cone on the defect background (a necessary condition for emergent local Lorentz invariance, not a full demonstration); the tangential cone shapes differ per channel (one tilt polarization rank-1, non-propagating along one direction).<br>`sandbox_vn/m5_8_2k_tilt_cone.py` |
| Molten clock classification (N-6d) | ✅ POSITIVE (classified)<br>The saturated breather is a MOLTEN CLOCK — a persistent ω₁ comb riding low-dimensional chaotic dressing (λ_max +0.4–0.7/τ, D₂ 2.7–3.0, control-anchored) whose intensity grows with excitation — the DTC literature's middle "melting" phase: a coherent clock dressed by low-dimensional chaos, not destroyed by it.<br>`sandbox_vn/m5_8_2n_chaos_battery.py` |
| The clock regularizes toward the ground state (N-6c) | ✅ POSITIVE<br>The coldest measured state (deep-settle, T = 2.36) reads near-regular (λ_max +0.11 ≈ the estimator bias floor, D₂ 1.68 vs 2.7–3.0 hot): cold = ticking clock, hot = molten — and full coldness is FORBIDDEN by the un-sittability that drives the self-start (the clock and its excitation are one object).<br>`sandbox_vn/m5_8_2o_omega_of_E.py` |
| First absolute ω (N-6b) | ⚠️ NEGATIVE-INFORMATIVE (structural gap)<br>Under two explicit postulates (H_rest ↔ 0.511 MeV; lattice action ↔ ℏ) the clock runs at 5.5×10¹⁹ rad/s — ~28× below the electron ZBW 2m_ec²/ℏ = 1.55×10²¹ — and because ω is rigid the gap is STRUCTURAL (not closable by energy bookkeeping), pointing at the V-on/Faber-r₀ core, the faithful kinetic, or the action↔ℏ postulate.<br>`sandbox_vn/m5_8_2j_zbw_ratio.py` (N-6b addendum) |
| First spin readout (N-6e) | ⚠️ NEGATIVE-INFORMATIVE (clean null + bound)<br>The rotation Noether charge (gated ⟨P,Ṁ⟩/(2T) = 1.000000 exact) reads J = 0 on the clock kick — the Θ-twist cancels over the hedgehog sphere, so the breathing channel carries no net frame angular momentum; the measurement is bounded by box torque at 24³ (an ℏ/2-class result needs torque-free boundaries), and spin-½, if carried, lives in the polarized-seed / hopfion / far-field sectors.<br>`sandbox_vn/m5_8_2p_spin_readout.py` |

---

## 4. Canonical implementation — what works

The build that produced every result above. This is the load-bearing configuration; anything not listed here is either superseded (the ψ engine, retired) or backlog (the V-on core, the faithful kinetic).

### 4.1 The substrate + action

| Element | Canonical form |
| --- | --- |
| Field | Matrix `M = O D Oᵀ`, `O ∈ SO(3)` embedded 4×4 (time = index 3), `D = diag(1, δ, 0, g)`, δ = 0.3 — biaxial, NOT Vector(3) ψ (the M5.2 closed negative). |
| Action | Duda Eq.18 `ℒ = −Σ F_μν F^μν − V(M)`, `F = [∂_μM, ∂_νM]`; the 4D clock sector uses the signed (Minkowski `(α,3)`-block) flux `F → ηFη`. |
| Saturating invariant | `V_u = u + β·u²` on the SIGNED spatial-curvature density `u = 2Σ_(i<j)⟨F_ij,F_ij⟩_s`, β = 1.558 (anchored 2β·\|u_min(seed)\| ≈ 3) — the 1D βR⁴ floor, transplanted to the channel that runs. |
| Potential | LdG `V(M)` on the SPATIAL 3×3 block only (the g axis is decoupled — the M5.8.1 spatial-block restriction); Faber regularization for the r₀ mass scale. |
| Seed | Boost-dressed biaxial hedgehog `O = O_hh(x)·exp(b·w(r)·B)·R(ωt)`, w = exp(−(r/R_W)²), R_W = 3.5; the dressed minimum (NOT bare hedgehog + kick) is the constraint-manifold seed. |

### 4.2 The integrator (the only stable scheme found)

| Element | Canonical form |
| --- | --- |
| Kinetic | The faithful F-commutator metric at CC level; production ships the simple `½‖Ṁ‖²` (frequency-corrected only — the 5d diagnosis). |
| Stepper | The constrained spectral-projection integrator: per-voxel 10×10 inertia operator, keep positive-inertia directions (λ > 0.05·max\|λ\|), project P onto the kept subspace EVERY step, global-(α,3)-momentum clamp. Every cheap positive-inertia kinetic has slow growing modes — this is the one that holds. |
| Precision | f32 on GPU (Metal), f64 numpy probes as the independent cross-check — f32 ≡ f64 to ~1% everywhere checked; the deep-floor cascade (the late-time stiffness limit) is the known numerical bound, kept off the load-bearing windows by dt-convergence. |
| Guards | No free global dressing amplitude (the 2b ghost channel); bounded-energy monitor (mass can't go ≤ 0); momentum-aware det𝕂 zero-mode guard. |

### 4.3 The validated method discipline (earned, reusable)

| Rule | Why it exists |
| --- | --- |
| Surrogate GUIDES, direct quadrature DECIDES | Track C: spline/4-pt surrogates carry mesh artifacts; only machine-exact direct evaluation is load-bearing. |
| The dt-invariance discriminator | Real dynamics reproduces at matched τ under dt-halving; stepper-driven growth shifts with dt — the test that killed the quadratic-runaway-is-numerical hypothesis. |
| The exact 9-pt phase average | Densities have bounded harmonics (V,K ≤ e^{i8ψ}, Q ≤ e^{i16ψ}) ⇒ 9 uniform phase points average exactly; the 4-pt set is not shift-invariant under twist. |
| FFT-window discipline | No peak within ~2 bins of DC; combs at bin multiples are automatic; vary the detrend window as a dial; pair signed probes with rectifying envelopes (caught the 0.262 artifact). |
| The knob gate | Validate a parameter family actually changes the seed (fixed-domain invariant spread > 5%) BEFORE spending compute (caught the RC no-op). |
| Headless-first validation | Every M5.8 result was produced headless (gates + npz + plots); rendering gates nothing scientific (the 2026-06-07 decision). |

### 4.4 What is NOT in the canonical build (scoped out)

| Item | Status |
| --- | --- |
| The electron mass/ω calibration | NG-1 — needs the V-on/Faber-r₀ core (the 28× absolute-ω gap lives here); the ratio infrastructure exists, the units fix does not. |
| The faithful kinetic in production | NG-1 residual — the 5d diagnosis says it corrects the clock FREQUENCY only; re-measure ω₁ under it when production-scale clock runs need it. |
| Spin-½ as a measured charge | Open — the clock channel is J-neutral; spin-½ (if carried) needs the polarized-seed / hopfion / far-field sectors + torque-free boundaries. |
| Multi-defect / composites | 15a — everything above is single-defect; ensemble behavior is the many-body question. |
| Strict single-line periodicity | Not claimed — the state is a molten clock (a coherent comb on low-dimensional chaotic dressing), regularizing toward the cold ground state. |

---

## Reproduction

Every script above lives under `research/sandbox_v8/` (the M5.8.0–2g arc) or `research/sandbox_vn/` (the N-1…N-6e ladder). The full command table — prerequisites, run times, expected outputs — is in [`0b_M5_roadmap.md`](0b_M5_roadmap.md) § *Reproduction commands*. Each script carries its own PREREQUISITE + RESULTS docstring with the caveats inline.
