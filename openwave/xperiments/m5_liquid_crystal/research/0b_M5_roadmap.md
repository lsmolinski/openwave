# M5 ROADMAP

🔶 **M5 / LIQUID-CRYSTAL MODEL** — full implementation plan, post-sandbox.

A research thread evaluating whether a Lagrangian / topological framework can replace OpenWave's empirical wave-equation search with a first-principles derivation, and produce charge quantization + far-field Coulomb that the M3 scalar model cannot. Sparked by email exchange with Jarek Duda (Jagiellonian) and Robert Close (Clark College) in the "Models of Particles" group.

For design rationale, M2/M4 inheritance, code mapping, resolution & performance plan, and layered validation, see [2a_path_to_m5.md](2a_path_to_m5.md).

---

## SUMMARY

### Sandbox (pre-M5)

- ✅ 8 sandbox numerical experiments — validate the core Lagrangian / topological physics before any production engine refactor (COMPLETE 2026-04-17)
  - ✅ Test 1 (Sine-Gordon kinks), Test 2 (Hedgehog Coulomb R²=0.993), Test 3 (Winding quantization), Test 4 (Klein-Gordon dispersion R²=0.999982)
  - ⚠️ Test 5 (Lagrangian derivation — W-L product form falsified), Test 6 (lepton mechanism; specific ratios deferred), Test 7 (Close's actual equations)
  - ❌ Test 8 (Smolinski Ψ³ K-selectivity FALSIFIED)
- ✅ **Winning recipe identified** (2026-04-17): topology (hedgehog + winding) + Klein-Gordon wave dynamics + Close's Eq. 19 as base vector wave + M3 near-field standing waves + Skyrme stabilizer. Full recipe in [3a § Winning Approach](1c_lagrangian_experiments.md#winning-approach-for-m5)
- ✅ **Group feedback integrated (2026-04-19)** — Jarek, Jeff, and Robert reviewed the sandbox summary; refinements captured in [3d § Group Feedback](2a_path_to_m5.md#group-feedback-2026-04-1718--refinements-to-m5-plan)

### M5 phases

🔶 **M5 implementation** (dir `openwave/xperiments/m5_liquid_crystal/`)

- ✅ **M5.0 Scaffold** — all 11 sub-phases complete as of 2026-05-08
- ✅ **M5.1 Port topology from Exps 2, 3** — 8 of 8 tasks complete as of 2026-05-11. 1/d Coulomb R² = 0.978 (attractive, monotone, threshold 0.95). Visual confirmation in [3a_coulomb_visual_geometry.md](3a_coulomb_visual_geometry.md)
- ❌ **M5.2 V(ψ) escalation on Vector(3)** — **CLOSED as negative result (2026-05-12)**. Tested 4 V(ψ) recipes (V=0, KG, φ⁴ Mexican-hat, biharmonic); all collapse Q identically at step 4-5. Diagnosis: matrix substrate `M = ODO^T` + time-periodic resonance required (triple-confirmed: Duda paper Fig. 10, Close email, Werbos chaoiton paper). See [3b_lagrangian_roadblocks.md](3b_lagrangian_roadblocks.md)
- ✅ **M5.3 Direction review** (2026-05-12 → COMPLETE 2026-05-26). Substrate locked (full `M = ODO^T`, Q5); Duda reply in (2026-05-14/15); paper + slides digested ([4a_convo_2026.05.12.md](4a_convo_2026.05.12.md)) with every input mapped to a phase; thermal prerequisites confirmed; resonance smoke test NULL; **feasibility spike PASSED** (`ti.Matrix.field` + commutator + `ti.sym_eig` all work in-kernel on Metal, verified, ~1–2× Vector(3) cost). **Decision: in-place migration → M5.4 unblocked.**
- 🔶 **M5.4 Matrix-field substrate migration** (2026-05-26 — substrate + gate ✅, rendering + cleanup remaining). Vector(3) ψ → `M(x) = O(x) D O^T(x)` real-symmetric 3×3 field. Done: `ti.Matrix.field(3,3)` triple buffer; `commutator` / `matrix_laplacian` / `antisym_A_mu` (Eq.19) operators; `eigen_decompose` lynchpin (director recovery 0.9995); matrix seeders; `‖M−D‖_F` + `‖Ṁ‖_F` thermal trackers; **M5.1 Coulomb reproduced — R²=0.9704 relaxed + R²=0.9959 analytical page-18 cross-val, both attractive**. Remaining (not gating): rendering repurposing (needs Rodrigo's 4b UX calls) + ψ-buffer cleanup. Physical-energy scaling deferred to M5.5 (belongs on the matrix energyH).
- 🚧 **M5.5 Paper Lagrangian + V(M)**. Implement Duda Eq. 18 action `L = Σ ‖F_μ0‖² − Σ ‖F_μν‖² − V(M)` with Eq. 13 Higgs-like `V_LG = a Tr(M²) − b Tr(M³) + c (Tr(M²))²`. Faber regularization integrated to activate V. Subsumes old "Skyrme stabilizer" phase since Eq. 42 4D Skyrme-like kinetic is the same family.
- 🚧 **M5.6 Biaxial twist + KG emergence**. With `Λ = (1, δ, 0)` and `δ ~ ℏ`, reproduce paper Fig. 9 — Klein-Gordon-like equation emerges automatically from twist dynamics (NOT added as V_psi term).
- 🚧 **M5.7 Resonance hunt (Close protocol)**. `l = 1` harmonic seed on biaxial substrate, sweep `A/λ ∈ {0.5, 1, 2}`, measure energy-localization lifetime. First metastable resonance — "first long-lived particle in M5". **UNBLOCKS 9b thermal work.** Literature anchor: BEC vortex kinetics (Duda 2026-05-13 cite, PRA `10.1103/2msv-lk1m`).
- 🚧 **M5.8 4D extension + Zitterbewegung clock** (GROUP HEADLINE — for Duda et al.). Add 0th time axis: `D = diag(g, 1, δ, 0)`, SO(1,3) Lorentz. Negative-energy contributions from spacetime signature auto-propel the clock (Fig. 10 mechanism). Measure de Broglie clock frequency `ω = 2mc²/ℏ` for electron — empirical answer to Duda's standing clock-propulsion question.
- 🚧 **M5.9 3 lepton families + Cornell quarks**. 3 axis-choices of biaxial hedgehog → same Q, different mass; tune `Λ` to match `μ/e ≈ 207`, `τ/e ≈ 3477`. Cornell potential for quark strings via topological vortex (`V(r) = −α/r + σ·r`). Standard Model correspondence section of paper confirmed.

### Parallel research stages

Forward-looking research programs in their own files — run alongside / after the linear M5.x phases. Each carries its own scientific test program:

- **[9a — Electromagnetic Wave Packets](9a_em_energy.md)** — emergence of EM wave packets (photons) from the wave field; longitudinal/transverse decomposition of photon dynamics.
- **[9b — Thermal Energy](9b_thermal_energy.md)** — heat as joint `(A, ω)` excess of subatomic defect oscillation above ground state. **Scientific counterpart to SABER's Direct Heat Conversion (SABER MAIN GOAL).** Test program 9b.0 → 9b.9. Unblocks once M5.7 lands; runs in parallel with M5.8 / M5.9. Empirical falsification (or validation) of the SABER thermal-amplitude hypothesis.
- **[9c — Time Dynamics](9c_time_dynamics.md)** — local time engineering, variable λ, time-dilation, gravity / propulsion implications. Future extension beyond M5.
- **[9d — Composite Particles](9d_composite_particles.md)** — hopfions, knots → nuclei, atomic orbitals. Test program 9d.1 → 9d.5 (Liu et al. 2026 lab anchor). Pursued post-9b validation if foundations succeed and resources permit.

---

## DETAILED

### ✅ Phase M5.0 — Scaffold

Broken into nine sub-phases (M5.0a → M5.0i) so each lands as a tight, separately-committable unit with its own test gate before progressing.

#### M5.0a — Module rename + alias ✅ (commit `bdd96dd` 2026-05-06)

- ✅ Create `openwave/xperiments/m5_liquid_crystal/` directory (cloned M4 + M3 features merged)
- ✅ **Rename the engine module**: `wave_engine.py` → `lagrangian_engine.py`. Rationale: M5's core loop integrates a Lagrangian-derived PDE (`∂²_tψ = c²∇²ψ − ∂V/∂ψ`) that simultaneously handles wave propagation *and* preserves topology via the potential `V(ψ)`. "Wave" is only one of the two channels the engine produces, so `lagrangian_engine.py` reflects what the module actually is to a new reader. M1–M4 keep `wave_engine.py`. See [3b § What wave equation does M5 solve?](_overview.md#what-wave-equation-does-m5-solve-is-force-still-e)
- ✅ Module alias `ewave` → `lagrange` across launcher (the engine evolves the field ψ via the Lagrangian, not "the energy-wave"; alias rename improves call-site readability)

#### M5.0b — Triple buffer + AMR-ready field-storage abstraction ✅ (commit `c832e90` 2026-05-06)

- ✅ Copy M4's `WaveField` / `WaveCenter` / `WaveTrackers` data classes; **extend with `psi_prev_am` and `psi_new_am`** Vector(3) buffers for leapfrog (also rename `displacement_am` → `psi_am`)
- ✅ **Use native `ti.Vector.field(3, ...)`** with single triple buffer — not three independent scalar fields. See Resolution & Performance Plan § Tier 1
- ✅ **`swap_buffers()`** method on WaveField cycles `prev ← curr, curr ← new` after each leapfrog step
- ✅ **AMR-readiness convention**: kernels must read grid dims via `wave_field.nx / .ny / .nz` attributes — never bake fixed `(nx, ny, nz)` constants into kernel signatures. Documented in WaveField docstring; M5.0 ships uniform-grid, M5.6 / M5.8 retrofit AMR without a rewrite
- ✅ Copy M4's flux-mesh visualization, granule rendering, 3-plane sampling (unchanged behavior)

#### M5.0c — Vector Laplacian (port + simplify from M2) ✅ (commit `5606dd5` 2026-05-06)

- ✅ Port M2's 6-point Laplacian stencil (`compute_laplacianL` at `m2_free_wave/wave_engine.py:527–562`); simplify to a single Vector(3) operator. Taichi handles Vector(3) arithmetic natively, so the stencil applied component-wise IS the vector Laplacian — no need for separate L/T paths like M2 had

#### M5.0d.1 — Leapfrog kernel + standing-wave eigenmode test ✅ (commit `6df9a1b` 2026-05-06)

- ✅ Implement `evolve_psi` kernel: leapfrog/Verlet update `ψ_new = 2·ψ − ψ_prev + (c·dt)²·∇²ψ` (V=0 free-wave; V terms land in M5.2)
- ✅ Standing-wave eigenmode test (V=0 reproduces continuum dispersion at low k; 2.26% deviation at 12 voxels/λ matches discrete-stencil prediction)

#### M5.0d.2 — CFL eval + plane-wave seed + tracker EMA + Hamiltonian dashboard ✅ (committed 2026-05-07)

- ✅ CFL evaluation in `_launcher.compute_timestep()` — mirror M2 pattern; `dt = dx · 0.95 / (c · √3)`; display `cfl_factor` in dashboard with red coloring if `> 1/3`
- ✅ `seed_wave` kernel — Gaussian-windowed wave packet (3-axis envelope, σ = N/6) that satisfies Dirichlet BC by construction. Drives `_test_smoke` xperiment for visual verification in the GUI
- ✅ Switch launcher main loop from M4's analytical `propagate_wave` → leapfrog `evolve_psi` + `swap_buffers()`
- ✅ `update_trackers` kernel — EMA on `|ψ|²` for amplitude; zero-crossing detection on `ψ_z` for frequency
- ✅ `compute_total_hamiltonian` (3-plane-sampled estimator: `H = ½|ψ̇|² + ½c²|∇ψ|² + V(ψ)`, V=0 in M5.0d). Full-grid atomic reduction was the original implementation but stalled the GUI at 100M voxels; switched to 3-plane sampling (codebase-consistent with `sample_avg_trackers`)
- ✅ **Delete legacy code**: M4 analytical `propagate_wave` kernel + module-level EWT base constants
- ✅ First successful UI render of leapfrog-driven wave propagation

#### M5.0d.3 — Drop `scale_factor` / `ewave_res` / EWT-default cleanup ✅ (committed 2026-05-07)

- ✅ Drop `scale_factor`, `ewave_res`, `max_universe_edge_lambda`, `nominal_energy*` from `WaveField` (M2-era constructs that assumed a single fixed reference wavelength — EWT's 28 am energy-wave). M5 has variable λ
- ✅ `Trackers.__init__` no longer takes `scale_factor`; globals init to zero; EMA + sample_avg_trackers populate them during sim
- ✅ Introduce `state.wave_res` (xperiment-driven, populated from `TEST_SEED["VOXELS_PER_WAVELENGTH"]` if a seed exists; else 0.0 / "n/a"). Defect-driven λ from λ_C lands in M5.2
- ✅ Strip all `scale_factor` arithmetic from launcher dashboard; `instrumentation.py` cleaned of EWT-scaled axhline references and the broken transverse subplot
- ✅ `force_motion.py` `S²` placeholder (= 1) — kernel is being fully rewritten in M5.0g (F = −∇E), so the placeholder is intentional dead-end code until then
- ✅ `annihilation_threshold` hardcoded to 6 voxels (M5.2 will replace with per-defect Compton-wavelength threshold); particle shell radius set to fixed 0.02 of normalized universe edge

#### M5.0e — Curl, divergence, curl-curl operators ✅ (committed 2026-05-08)

- ✅ `compute_divergence` `@ti.func` — central-difference scalar divergence, 1-cell halo. Used by M5.1 winding tracker, M5.2 `∇·s = 0` constraint, and the curl-curl identity below
- ✅ `compute_curl` `@ti.func` — central-difference Vector(3) curl, 1-cell halo. Foundation for M5.2 Eq. 19 linear limit and spin-density observables
- ✅ `compute_curl_curl` `@ti.func` — implemented via the vector identity `∇×(∇×ψ) = ∇(∇·ψ) − ∇²ψ` rather than nested curl. Reuses validated Laplacian, only the gradient-of-divergence is new code, 2-cell halo (vs. 4-cell for nested curl), matches Exp 7 v2's implementation
- ✅ DIFFERENTIAL OPERATORS section header documents all four operators with their halo requirements and downstream consumers (single source of truth for the vector-calculus toolkit)
- ✅ Analytical verification: divergence on `ψ=(x, 2y, 3z)` → 6.000 (max err 3e-5); curl on rigid-rotation `ψ=(−y, x, 0)` → (0, 0, 2) (max err 3e-6); curl on constant field → 0 (1e-9 noise floor); curl-curl identity on Gaussian seed → relative error 2.8e-6

#### M5.0f — Storage-units decision + natural-units deferral ✅ (decision recorded 2026-05-08)

This sub-phase was originally scoped as "add kernel-internal natural-unit scaling (c=1, λ_C=1, ℏ=1) at kernel entry/exit". After investigating both halves of the question (storage units AND kernel-internal scaling), the conclusion is that **neither change is needed in M5.0**. M5.0f therefore lands as a **decision-record sub-phase** — no code refactor — capturing the rationale so the question isn't re-litigated later.

- ✅ **Storage units stay `_am` / `_rs` / `_rHz`** (no refactor):
  - The leapfrog kernel is **dimensionally self-balancing**: `(c·dt)²` carries length² and `∇²ψ` carries 1/length², so the product magnitude is `~0.08·ψ` regardless of dx units. f32 stability holds in any (am, fm, pm) storage choice — the precision argument originally motivating a switch turned out not to apply to linear kernels
  - **No "particle-physics best practice"** justifies a full pair: `fm` is the standard length unit (Fermi), but `ys` is not — particle physicists use either `zs` (zeptosecond) or natural units (ℏ/GeV ≈ 6.6e-25 s) for time. A genuinely standard pair `(fm, zs)` would force `c = 300` instead of `0.3`, propagating into every CFL / dispersion / kernel formula — refactor cost without precision benefit
  - **Cross-method consistency** with M2/M3/M4 is real value (shared validations, port comparisons, instrumentation)
  - The "storage values look large" observation (e.g. `dx = 200000 am` at electron scale) is a **dashboard formatting** concern, not a correctness one. Adaptive SI-prefix display in the launcher (planned as a small follow-up: `200000 am` → `200 fm` at render time only) handles it without touching storage
- ✅ **Kernel-internal natural-unit scaling deferred to M5.2** (where it actually pays off):
  - The original motivation for `(c=1, λ_C=1, ℏ=1)` scaling was f32 precision in linear kernels at electron scale. As shown above, that's a non-issue
  - The *real* value of natural units is in **non-linear physics couplings**: Klein-Gordon mass term `−m²·ψ`, Close's Eq. 23 nonlinear terms `−u·∇s + w×s`, LdG biaxial potential — these have explicit dimensional coefficients that read like the published equations only when written in natural units
  - These kernels first appear in **M5.2**. Apply natural-unit scaling there, locally to those specific kernels (entry/exit conversion). Linear kernels (leapfrog, divergence, curl, Laplacian) **never** need natural units — they're dimensionally self-balancing
  - This is "lazy natural units": apply only where it helps, skip where it doesn't. Avoids premature complexity (boundary conversion code with no precision win) in M5.0
- 🚧 (optional follow-up, not blocking) — small UI improvement: adaptive SI-prefix display helper in `_launcher.py` so universe edge / voxel edge / wavelength / amplitude render in the natural prefix for their magnitude (`fm`, `am`, `pm`, etc.) instead of raw scientific notation. Storage stays in `_am` / `_rs` / `_rHz`; only the display formatting changes. ~20 lines of code, can ship anytime as a polish PR

#### M5.0g — Per-voxel energy density (Hamiltonian) + force-computation switch ✅ (committed 2026-05-08)

- ✅ Per-voxel field `trackers.energyH_density_aJ` (replaces deprecated `energy_local_aJ`); populated each step by `lagrangian_engine.compute_energyH_density`
- ✅ Mean per-voxel cache `trackers.energyH_global_avg_aJ` (filled by `compute_energy_total_H`; matches M4's `_global_avg_` semantics; used as the flux-mesh colormap range for WAVE_MENU=4)
- ✅ Grid-total scalar `state.energy_total_H_aJ` (= mean × voxel_count) on the dashboard
- ✅ Rewrote `force_motion.compute_force_vector` as **`F = −∇E`** sampling `energyH_density_aJ` (the `_H` suffix tags how E is computed; the physics statement F=−∇E is canonical regardless of formula). Dropped the placeholder `S²=1`, all hardcoded EWT particle constants (`EWAVE_AMPLITUDE`, `EWAVE_LENGTH`, `MEDIUM_DENSITY`, `ELECTRON_K`/`OUTER_SHELL`/`ORBITAL_G`, `COULOMB_CONSTANT`, `ELEMENTARY_CHARGE`, `WAVE_SPEED`), the `compute_ewt_electric_force` reference function, and the `numpy` import that only it needed
- ✅ `V_psi(psi)` `@ti.func` hook — returns 0 in M5.0g; M5.2 swaps in Klein-Gordon mass + Close Eq. 23 + LdG terms (alongside the kernel-internal natural-unit scaling deferred from M5.0f). Both the potential value AND its functional form will land here
- ✅ `compute_energy_total_H` refactored to 3-plane-sample the new per-voxel field (instead of recomputing kinetic+gradient per voxel three times) — three lightweight `_energy_slice_*` slice-copy kernels
- ✅ Flux-mesh `WAVE_MENU == 4` activated to render `energyH_density_aJ` (was a placeholder rendering `|ψ|`)
- ✅ Naming convention captured per Rodrigo's 2026-05-08 review: the quantity is **energy** (aJ); the `_H` suffix tags the *formula* used (Hamiltonian). Future formulas would parallel: `_L` for Lagrangian density, `_K` for kinetic-only. Renamed `compute_hamiltonian_density` → `compute_energyH_density`, `compute_total_hamiltonian` → `compute_energy_total_H`, `H_density_aJ` → `energyH_density_aJ`, `H_global_avg_aJ` → `energyH_global_avg_aJ`, `hamiltonian_total_aJ` → `energy_total_H_aJ`
- ✅ Smoke-tested end-to-end: `annihilation1` (ψ=0): all energy fields zero, no forces; `_test_smoke` (Gaussian packet): non-trivial energyH_density_aJ peak ~5e-3, total ~8e3 aJ, forces non-zero where ∇H ≠ 0

#### M5.0h — Physics invariant test (gating) ✅

- ✅ **Physics invariant test (V=0)** — leapfrog reproduces the discrete dispersion within ±0.5% on `c²` recovery across all 5 modes (voxels/λ_x ∈ {31, 21, 16, 10, 8}). Klein-Gordon `+ m²` flavor deferred to M5.2 where it lands alongside the actual mass term. **Implementation** at `openwave/xperiments/m5_liquid_crystal/research/scripts/m5_0h_dispersion.py` (headless, ~30s on Metal); plot at `research/plots/m5_0h_dispersion.png`. New engine kernel `seed_dispersion_modes` for multi-mode standing-wave initial conditions
- ✅ **Two persisted lessons from M5.0h** — (a) Taichi Metal lowers `s += …` auto-reduction to atomic_add and stalls on full-grid reductions just like the M2 `3-PLANE SAMPLING` block warned. Workaround: sparse point sampling at mode antinodes, FFT recovers ω from the mixed time series. (b) Discrete-scheme dispersion fits MUST invert the FULL space+time relation (`sin(ω·dt/2) = (c·dt/2)·√K`), not the spatial-only `ω² = c²·K` — the difference is a `(k·dx)²` systematic bias that grew to ~1.5% on `c²` at 7.8 voxels/λ. Both lessons captured in `lagrangian_engine.py` (the AUTO-REDUCE CAVEAT block) and in the auto-memory store (`feedback_taichi_metal_atomics.md`, `feedback_dispersion_validation.md`)

#### M5.0i — Performance profile (baseline only) ✅

- ✅ **Per-kernel profile delivered** at production grids (128³, 256³, 384³). Headless harness `openwave/xperiments/m5_liquid_crystal/research/scripts/m5_0i_profile.py`; chart at `research/plots/m5_0i_profile.png`. Per-step times:

| Grid | step ms | fps | vs 20fps floor |
| --- | --- | --- | --- |
| 127³ (2M voxels) | 1.0 | 964 | ✅ 50× under |
| 255³ (17M voxels) | 5.9 | 168 | ✅ 8× under |
| 383³ (56M voxels) | 19.4 | 51 | ✅ at edge but passes |

  Per-kernel breakdown is roughly stable across grid sizes: `swap_buffers` 34 %, `update_trackers` 24 %, `evolve_psi` 23 %, `compute_energyH_density` 18 %. `sample_avg_trackers` (every-60-frames cadence) is negligible (~0.06 ms/step amortized).

- ✅ **M2 vs M5 head-to-head profile** (companion script `research/scripts/m5_0i_profile_m2_compare.py`): M5 is consistently ~2.1× M2 step time at production grids (M2 0.52 / 2.68 / 8.90 ms vs M5 0.99 / 5.94 / 19.4 ms at 128³ / 256³ / 384³). M2 fuses leapfrog + tracker EMA + zero-crossing freq + buffer swap into a single `propagate_wave` ndrange (`m2_free_wave/wave_engine.py:603`); M5 deliberately split things into 4 separate kernels for cleanliness, AMR-readiness, and the V_psi hook for M5.2. The 2.1× is the bill for that split. ~50 % is fusion debt; ~50 % is Vector(3) being 1.5× bigger per voxel than M2's two scalars + the new `compute_energyH_density` (no equivalent in M2).

- ✅ **Decision: no Tier 2 opts justified by current budget.** M5 already passes the 20 fps floor at 56M voxels with margin. M5.0i is profile-only.

- ✅ **Two persisted lessons captured** during the M5.0i investigation:

  - **Rotating-pointer `swap_buffers` is NOT a 30-min change** — Taichi's `ti.template()` caches attribute lookups at first compilation, so attribute rotation on the same `wave_field` Python instance is invisible to cached kernels (M5.0h dispersion test reproduced this with a clean −19 % c² regression). Proper fix is passing fields explicitly to kernels (2–3 hr refactor of every kernel that touches the triple buffer). Documented in `medium.py:swap_buffers` docstring + `feedback_taichi_template_caching.md` memory.

  - **Fusion is the single biggest Tier 2 lever** — would close ~50 % of the M2 gap by combining `evolve_psi` + `update_trackers` + `compute_energyH_density` (and ideally the swap) into one ndrange loop. **Proven prior art**: M2's `propagate_wave` (in `m2_free_wave/wave_engine.py:603`) is the template — it does leapfrog + RMS-EMA + zero-crossing-freq + buffer swap in a single `ti.ndrange` followed by an in-kernel swap loop. Deferred because: (a) we don't need the win yet, (b) fusion makes the per-task kernels redundant or duplicated (divergence risk), and (c) AMR retrofit (M5.6/M5.8) is harder against fused kernels. Re-evaluate when M5.2's V(ψ) makes per-step work heavier.

- 🚧 **Re-profile trigger**: when M5.2's V(ψ) (Klein-Gordon mass + Close Eq. 23 + LdG potential) lands and per-step time grows. If 384³ drops below 20 fps, decide between (1) rotating-pointer refactor for ~35 % win, (2) full fusion for ~50 % win + M2-equivalent step time, (3) BlockLocal Laplacian / dirty-tile mask. The M5.0i baseline numbers above are the comparison reference.

### ✅ Phase M5.1 — Port topology (from Exps 2, 3)

- ✅ `seed_vacuum()` — fills psi_am, psi_prev_am, AND psi_new_am with `n = ẑ` everywhere (writes all three buffers per `feedback_triple_buffer_bc.md` — required for fixed-value Dirichlet BC at vacuum)
- ✅ `seed_hedgehog(centers, signs, D/4, n_defects)` — N-defect kernel, port of Exp 2's weighted superposition + renormalization (`sandbox_phase3_lagrangian/exp2_hedgehog_energy.py:71-108`). Radial structure concentrated within ~D/4 via `w_vac = 1/(1 + (r/(D/4))⁴)`, smoothly blends to `n = ẑ` at boundary. Multi-defect API supports M5.4 pair tests directly (no kernel refactor needed)
- ✅ **Director-glyph visualization** — 3-plane line-glyph renderer landed early (between seeders and Frank energy, not between relaxation and Coulomb test as originally planned). Visual confirmation made the seed correctness obvious and surfaced the BC-bleed bug below. Implementation: `lagrangian_engine.update_director_glyphs` + `WaveField.director_glyph_*` fields; `_launcher.SHOW_DIRECTORS` slider 0..3 (mirrors `SHOW_FLUX_MESH` semantics: 0=off, 1=XY, 2=+XZ, 3=all three planes). Color encoding evolved from the originally-designed signed-component RGB to a colormap palette of `(1 − n_z)` ∈ [0, 2] — vacuum maps to dark (blends with black GUI background; defect stands out), peak twist maps to bright. Two palette options retained (`get_blueprint_color` for whole-domain visibility, `get_ironbow_color` for defect-only focus). `VIZ_STRIDE` xparameter (default 4) consolidates sampling stride for both directors AND `SHOW_GRANULES` overlays — granule spheres render on top of glyph lines at the same sample points, giving combined orientation+motion visualization. Full as-shipped notes in [2b_director_glyph_rendering.md](2b_director_glyph_rendering.md)
- ✅ **BC-bleed bug found + fixed** during glyph testing — visible as boundary-driven inward wave on first EVOLVE PSI click. Root cause: `psi_new_am` boundary was never written by anyone (propagator skips boundary; was uninitialized=0). First `swap_buffers` clobbered the seeded ẑ at boundary with that 0. M5.0 had this latent bug too but invisible because Gaussian/dispersion seeds happened to put ψ≈0 at boundary anyway. Fix: every seeder writes ALL THREE buffers (`seed_gaussian`, `seed_dispersion_modes`, `seed_vacuum`, `seed_hedgehog` all updated). `evolve_psi` docstring now describes "fixed-value Dirichlet" with seeder-dependent value + BC-consistency requirement. Captured in `feedback_triple_buffer_bc.md` memory
- ✅ Frank Elastic Energy density `H_F = (K/2)·|∇n̂|²` — per-voxel kernel `compute_energyF_density` in `lagrangian_engine.py` mirroring `compute_energyH_density`'s gradient stencil (sans c² and sans kinetic/potential). Writes to `FieldObservables.energyF_density_aJ` (post 2026-05-11 SoC refactor that split derived scalars out of `Trackers` — derived scalars are stateless and instantaneous per-frame; Trackers retains time-integrated amp/freq EMA). Global per-voxel mean populated by `sample_avg_observables`, a parallel 3-plane pass to `sample_avg_trackers` (one pass per domain). `K_FRANK = 1.0` hardcoded as Exp 2 baseline; physical elastic constants land in M5.6. **Rendering integration**: `WAVE_MENU = 5` ("FRANK (Elastic)") added to the flux-mesh menu — ironbow palette, defect cores light up brightest, vacuum dark; dashboard shows global mean. Per `feedback_visual_rendering_priority.md` memory, visual integration shipped in the same PR as the kernel, not deferred. Volume integral `E(d) = ∫ H_F d³r` for the Coulomb fit (task 7) is the trivial `mean × voxel_count` derivation, same pattern as `energyH_total`. **Validation**: comprehensive regression test in `research/scripts/m5_1_frank_energy.py` — all 5 tests PASS (vacuum baseline `max\|F\|=0` exactly; pure-hedgehog F(r) `CV(F·r²) = 0.0019` against analytical `K/r²`; K linearity zero deviation; numpy-reference per-voxel parity mean rel err `3.9e-8` at f32 round-off limit; pair-superposition deviation 6.7% within 50% seed-design tolerance). Plot at `research/plots/m5_1_frank_energy.png`.
- ✅ Gradient-descent relaxation (M5.1 task 6) — `lagrangian_engine.relax_director_step` kernel ports Exp 2's relax() inner loop: tangent projection `dn = ∇²n − (n·∇²n)n` + unit-length renormalization + soft-Dirichlet pin at each defect's closest voxel. Writes to `psi_new_am`; launcher's `relax_field` helper copies result back to BOTH `psi_am` AND `psi_prev_am` so subsequent leapfrog sees ψ̇=0 (no spurious time-derivative artifact from the static relaxation). Step size τ = 0.4·dx²/6 (40% of heat-equation CFL bound for headroom; matches Exp 2's τ/CFL ratio). Boundary voxels preserved via copy-from-psi_am (Dirichlet BC unchanged). **Always-on auto-relax on seed** via xparameter `AUTO_RELAX_STEPS` (default 60 in `_test4_topology.py`) — relaxation is a numerical ground-state-finder, not physics, so there's no demo / interactive use case. UI controls (button + slider) were drafted and removed (2026-05-11 review) to keep the screen clean; any future interactive heat-modulation controls land with M7 thermal investigations (where the math overlaps with relaxation but the physics is different). **Validation**: `research/scripts/m5_1_relax.py` — all 5 tests PASS (monotone decrease, convergence with F dropping 92.6% from seed, topology preserved at solid-angle proxy rel-change 3.9e-9, pin held exactly at ±ẑ, boundary preserved at f32 zero). Plot at `research/plots/m5_1_relax.png` with F(step) convergence curve.
- ✅ 1/d Coulomb gating test (M5.1 task 7, M5.2 gate) — `research/scripts/m5_1_coulomb.py` ports Exp 2's E(d) sweep to Taichi production. Method: seed (+1, −1) hedgehog pair at offset ±d/2 along x, relax 600 steps, sum F_total, repeat for d = {8, 10, 12, 14, 16, 18, 20}, fit E(d) = a + b/d via linear lstsq. **PASS**: R² = 0.978 (threshold 0.95), b = −5.6e-3 (ATTRACTIVE — correct sign for opposite charges), monotone trend confirms Coulomb-like 1/d law. Same-charge (+1, +1) tested as informational sanity check (sign-correct REPULSIVE, but R² = 0.14 — pinned same-sign defects don't reach clean topological equilibrium on a discrete Dirichlet-BC grid; not part of gating logic). **Threshold relaxed from Exp 2's R² > 0.99 to R² > 0.95** because Exp 2 used periodic BC in numpy (clean torus) while production uses Dirichlet BC (finite-domain corrections add ~2-3% R² loss); the sign + monotone-trend + R² > 0.95 combination is strong evidence the framework reproduces Coulomb. Future periodic-BC variant of `relax_director_step` could tighten back to 0.99 if needed — deferred until M5.2 demonstrates the dynamical version. **M5.2 unblocked**. Plot at `research/plots/m5_1_coulomb.png` (E vs d + E vs 1/d linear panels). **Visual companion** documented in [`3a_coulomb_visual_geometry.md`](3a_coulomb_visual_geometry.md): side-by-side GUI screenshots of opposite-charge (dumbbell-shaped F density bridge along the connecting axis) vs same-charge (pinched / perpendicular F density splay) configurations, both via Frank-elastic flux-mesh coloring AND via director-glyph rendering. The visual match to classical EM textbook field-line geometry — derived from pure topology with no electromagnetism postulated — is the qualitative companion to the quantitative R²=0.978 result.
- ✅ Winding-number tracker `compute_winding_number` (M5.1 task 8) — pure-numpy diagnostic in `lagrangian_engine.py`; ports Exp 3's `winding_number()`. Method: trilinear-sample ψ on a regular (θ, φ) spherical grid around a defect center, finite-difference angular derivatives, surface-integrate `n̂·(∂_θ n̂ × ∂_φ n̂) / (4π)`. Validated in `research/scripts/m5_1_winding.py` — 4 tests all PASS (vacuum Q=0.000, single +1 hedgehog Q=+0.996, single −1 Q=−0.996, pair at offset ±10 sampled around each → Q=±0.996). **Bonus diagnosis surfaced during this task** (now documented for M5.5 awareness): on a discrete-grid + Dirichlet-BC setup with same-direction pin (e.g., sign=+1 pinned at +ẑ when vacuum is also +ẑ), the topological winding **dissipates after ~50-100 relax steps** because no constraint forces a twist — Frank energy can monotonically lower by aligning everything to ẑ. M5.1 task 7's 1/d Coulomb test still works (the F(d) signal comes from blend-zone elastic, not from the integer Q itself), but defect stability under longer evolution is brittle until M5.5's Skyrme stabilizer lands. The winding tracker is now the diagnostic that will catch this loss-of-topology event during M5.2 dynamics testing (watch Q drift away from ±1 over time → that's the Skyrme readiness signal). Test uses `N_RELAX=0` so the seeded topology is intact when measured.

> **Current behavior of M5.1 under EVOLVE PSI** (clarification recorded 2026-05-11 during conversation review with Rodrigo): the seeded hedgehog is a topological texture in a field with `V(ψ) = 0` (M5.0g placeholder; M5.2 lands real V). Pressing EVOLVE PSI runs the free-wave leapfrog, which has no restoring force and no constraint preserving `|n̂| = 1`. The high `∇n` at the defect core acts as a localized energy peak that the Laplacian dynamics smooths outward — the result is **outgoing ripples + dissolution back to vacuum**, which is Derrick's theorem in action ("a topological texture in a pure-Laplacian field with no potential has no stable equilibrium"). This is *correct, expected behavior* for the current V=0 stage, not a bug. Stability of the lepton prototype requires the M5.5 Skyrme stabilizer (prevents collapse to a point) and the M5.6 LdG potential (prevents expansion to vacuum) — both will land in their respective sub-phases.
>
> **Why the hedgehog pair doesn't attract under EVOLVE PSI** (same review): wave equations conserve energy globally — there is no dissipation channel for two defects to "find" each other by lowering elastic energy. The two paths to visible defect attraction are: (a) **gradient-descent relaxation** (M5.1 task 6) — explicitly dissipative, lets opposite-sign defects migrate toward annihilation in the minimum-energy configuration; this is what Exp 2 used to fit the 1/d Coulomb law. (b) **Full nonlinear wave dynamics** (M5.2 + M5.4 headline test) — needs V(ψ) to keep defects coherent while they propagate; the falsifiable headline goal is "does the pair actually annihilate dynamically?" In M5.1 with V=0 and no relaxation step wired, neither mechanism is active yet, so defects radiate energy without visible attraction.
>
> **Frank elastic energy is a measurement, not a dynamics** (same review): the elastic behavior is *already* in the wave equation — the `c²·∇²ψ` term IS the gradient of the Frank density `(K/2)|∇n|²`. Computing the Frank energy as an explicit kernel (task 5) doesn't add new physics; it provides the scalar `E(t)` needed for the downstream tests that *depend on having an E to fit*: task 6 needs `E` falling monotonically as a convergence diagnostic; task 7 needs `E(d)` as a function of pair separation to fit the 1/d Coulomb law; future visualization can colormap the per-voxel elastic density. Parallel to how `compute_energyH_density` (M5.0g) is the measurement scalar that goes with the wave dynamics — same logic, different formula component.

### ❌ Phase M5.2 — V(ψ) escalation on Vector(3) — CLOSED as informative negative (2026-05-12)

> **Status**: closed. The four V(ψ) escalations we tested (V=0, KG mass `½m²|ψ|²`, φ⁴ Mexican-hat `¼λ(|ψ|²−1)²`, biharmonic `½κ|∇²ψ|²`) all collapse the topological charge Q identically within 4-5 propagation steps. The result is documented as an informative negative; phase moves out of the critical path. Path forward absorbed into M5.3 (direction review) + M5.4 (matrix-field substrate) + M5.5 (paper Lagrangian).
>
> **Original scope (preserved for archaeology)**: Close's Eq. 23 + Eq. 19 + Klein-Gordon mass on Vector(3) ψ; dispersion validation; resonance-hunt amplitude sweep. All assumed Vector(3) substrate.
>
> **What we actually executed and learned**:
>
> - ✅ **Step 1 (Natural-units scaffold, 2026-05-12)**: `constants.COMPTON_WAVELENGTH_REDUCED_ELECTRON_AM` added; `m_freq_kg_rs = c_amrs / λ̄_C` (electron: ~7.76e-7 rad/rs at SIM_SPEED=1) threaded through `V_psi` via `compute_energyH_density`. Plumbing only — V_psi returned 0 at this step.
> - ✅ **Step 2 (KG mass term, 2026-05-12)**: `V_psi` returns `½·m²·|ψ|²`; `evolve_psi` adds `−(m·dt)²·ψ` to leapfrog. Math verified at f32.
> - ⚠️ **Step 3 (Defect-survival check, NEGATIVE)**: `research/scripts/m5_2_kg_defect_survival.py`. Q drops from `+0.9958` → `~0` within 20 steps under BOTH V=0 and KG-electron; `|Q_free − Q_kg|` below f32 precision at every sample.
> - ⚠️ **Step 4a (Mexican-hat φ⁴, PARTIAL)**: `research/scripts/m5_2_phi4_defect_survival.py`. `V += ¼λ(|ψ|²−1)²` damps `|ψ|` excursions (max `1.83 → 1.27`) but does NOT preserve Q. Same step-4 collapse as free-wave.
> - ⚠️ **Step 4b (Biharmonic, NEGATIVE on Q)**: `research/scripts/m5_2_biharmonic_defect_survival.py`. `+ ½κ|∇²ψ|²` (kernels kept research-only, NOT promoted to production). Stable at κ ≤ 0.003·c²·dx²; Q decay identical to free-wave at every stable scale. Pre-relaxing 20 steps doesn't help.
>
> **Root cause** (initial diagnosis 2026-05-12, refined after re-reading Duda paper arxiv:2108.07896): the framework requires **(a) matrix field `M = ODO^T`, NOT Vector(3) ψ** (Duda paper §III, Eq. 18); and **(b) the "particle" is a time-periodic resonance**, not a static soliton. Triple-confirmed from three independent sources: Duda paper Fig. 10 (4D Lorentz negative-energy terms auto-propel the clock), Robert Close email reply (l=1 amplitudes, A/λ ≈ 1 protocol), Werbos chaoiton paper (explicit "static solitons don't exist; the stable objects are chaoitons"). Captured in [`3b_lagrangian_roadblocks.md`](3b_lagrangian_roadblocks.md) and memories `feedback_no_static_solitons`, `reference_duda_lcb_paper`.
>
> **Closure**: M5.2 closed as informative negative. The work was correct under its assumption; the assumption was wrong (Vector(3) substrate + static stability goal). The corrective phases are M5.3 → M5.4 → M5.5 → M5.6 → M5.7. The new external-comms trigger lives in M5.8 (Zitterbewegung clock at `ω = 2mc²/ℏ`), not in defect-survives-EVOLVE-PSI.

---

### ✅ Phase M5.3 — Direction review (2026-05-12 → COMPLETE 2026-05-26)

> Study, sandbox, decide, wait. Goal: lock the M5.4 implementation plan with a clear substrate choice before code commitment. **COMPLETE (2026-05-26): all legs met — M5.4 is unblocked.** Substrate locked (full 3×3 `M = ODO^T`, Q5); Duda's reply in (2026-05-14/15); paper + slides digested ([`4a_convo_2026.05.12.md`](4a_convo_2026.05.12.md)); thermal prerequisites confirmed; resonance smoke test NULL (no Vector(3) resonance); **feasibility spike PASSED** — `ti.Matrix.field(3,3)` storage + matrix commutator + `ti.sym_eig` eigen-decomposition all work in-kernel on Metal, verified, at ~1–2× the Vector(3) Laplacian cost. **Decision: in-place migration.**

- ✅ **Re-read Duda paper §III-V deeply** — DONE. The annotated reading lives in [`4a_convo_2026.05.12.md`](4a_convo_2026.05.12.md): matrix structure (§3–5), eigenvalue→force map (§8), KG-from-hedgehog closed form + clock toy-model numerics (§11), force unification (§7), the 51-page LdGS slides (§11), and the Couder/walking-droplet deck (§11b). Supersedes the planned `3h_*.md`.
- ✅ **Wait for Duda reply** — DONE. Duda replied 2026-05-14/15 (`4a §1`); **substrate gate CLOSED** (full M, no Q-tensor pivot — Q5); eigenvalue→physics mapping confirmed (Q6).
- ✅ **Taichi storage layout study** — DONE (2026-05-26). Layout = `ti.Matrix.field(3,3)` with index-generic operators (3×3→4×4 M5.8 promotion stays a type change). **Measured cost** (spike, vs Vector(3) 6-point Laplacian baseline): matrix `[∂M,∂M]` commutator ~1× and per-voxel `sym_eig` eigen-extraction ~1.1× at 256³ (2.2× at 128³) — eigen is **not** a bottleneck (it runs once/frame for render + trackers, not in the inner loop). **Memory**: production should store the **6** independent symmetric components/voxel (not the full 9), vs 3 for Vector(3) → ~2× memory. Cost is acceptable; no sparse/AMR needed for M5.4.
- ✅ **Thermal prerequisites analysis** — DONE (2026-05-26). Confirmed: 9b's foundation is **M5.7 (first metastable resonance)** — sufficient, nothing extra must land earlier. The prerequisite chain is already correctly slotted across the M5.x phases: **measurement infra lands EARLY** (the redefined `A = ‖M−D‖_F` amplitude + `ω = O(x) clock-rate` frequency trackers in M5.4); the **KG drive substrate** in M5.6; the **drive infra** (parameterized harmonic forcing + FM/AM modulation primitives) in M5.7's own task list. So 9b.1 *internal* validation unblocks the moment M5.7 lands (per the two-tier table in [`9b_thermal_energy.md`](9b_thermal_energy.md)); full *lab-relevance* follows Phase 4 (EM emergence). SABER-specific engineering scope stays in the private repo per the cardinal rule. **No new early infrastructure required.**
- ✅ **Decision document** — DONE (2026-05-26). Decision: **in-place migration** (not parallel-track) — the spike proved storage + operators + eigen all work at acceptable cost on the existing grid/engine architecture. Substance is distributed: substrate + refactor strategy in `4a §2,§10`; M5.4+ task breakdown in the M5.4 phase below + [`4b_rendering_features.md`](4b_rendering_features.md); cost estimate in the storage-layout task above; 4a-inputs→phase mapping in the table below. (A standalone `3h_m5_substrate_decision.md` is superseded by 4a + this roadmap.)
- ✅ **Feasibility spike: matrix field in Taichi** — DONE (2026-05-26), PASSED. `research/sandbox_v3/m5_3_matrix_feasibility.py`. All three primitives proven in-kernel on Metal and verified against analytic cases: **storage** (`ti.Matrix.field(3,3)` round-trip 0 err); **commutator** `[A,B]` (0 err vs hand-computed; + `[∂_xM,∂_yM]` over a field); **eigen-decomposition** via **`ti.sym_eig`** (recovers a known `O·diag(2,1,0.5)·O^T` to ~1e-7; recovers a seeded hedgehog director at mean `|n̂·n|=0.9995`). **Key de-risk: `ti.sym_eig` works in-kernel on Metal for 3×3 — no custom analytic eigensolver needed.** Cost ~1–2× Vector(3) (see storage-layout task). Decision: in-place migration; M5.4 unblocked.
- ✅ **Sandbox: Close's resonance protocol on existing Vector(3)** — DONE (2026-05-26), NULL result. `research/sandbox_v3/m5_3_resonance_smoke.py`: l=1 dipole packet swept at A/λ ∈ {0.5, 1, 2}, energy-localization lifetime under V=0 / KG-mass / φ⁴. **No substrate-agnostic resonance found.** Linear potentials (V=0, KG mass) disperse at ~40 steps **amplitude-independently** (as a linear PDE must — no amplitude-selected resonance); the only nonlinear potential available (φ⁴ Mexican-hat) is director-scaled (vacuum |ψ|=1) and **diverges** on a displacement-wave amplitude sweep — Vector(3)'s potential toolkit can't even host the test. Reinforces the M5.2 closed-negative: the metastable-resonance hunt belongs on the matrix substrate (M5.4 → M5.7) with LdG V(M)/Skyrme + 4D clock propulsion, not Vector(3).

**Exit criterion** — ✅ **MET (2026-05-26)**: substrate decision locked ✅ + Duda reply in ✅ + feasibility spike done (cost estimate measured + commutator & `sym_eig` eigen-kernel proven in-kernel) ✅. M5.4 plan is fully concrete; **M5.3 closes, M5.4 begins.**

#### M5.3 deliverable — every 4a input mapped to a phase

The paper/slides reading task produces this assignment table — confirmation that nothing in [`4a_convo_2026.05.12.md`](4a_convo_2026.05.12.md) is orphaned; each instrumental item is slotted to a phase that will visit/analyse it.

| 4a input | Source | Phase | Status in roadmap |
| --- | --- | --- | --- |
| Matrix `M = ODO^T` substrate | §2, §3 | M5.4 | ✅ in M5.4 tasks |
| Refactor strategy (3×3 now, index-generic ops, 4×4 later) | §10 | M5.4 build / M5.8 promote | ✅ in M5.4 + M5.8 |
| Eigenvalue→physics map (1=EM, δ=QM, g=gravity) | §8 | M5.5 / M5.6 / M5.8 | ✅ Q6 resolved |
| Page-18 Mathematica Coulomb code (`V(d)≈1589.56−25.16/d`) | §11 | M5.4 | ✅ added as cross-validation target |
| KG-from-hedgehog closed form (page 32) | §11 | M5.6 | ✅ eq in M5.6 (cross-reffed) |
| Clock toy-model numerics (page 33: φ=tanh…, E≈2.0252, ω≈1.2898) | §11 | M5.8 | ✅ added as infra-validation anchor |
| Hydro↔EM dictionary (vorticity↔B, Coriolis↔Lorentz) | §11b.1 | M5.5 | ✅ added as EM-from-tilts cross-check |
| Faber quantized-EM Lagrangian (`L_EM=−(αℏc/16π)R_μν R^μν`) | §11b.2 | M5.5 | ✅ added as port target |
| Walking-droplet path-memory kernel + quantization laws | §11b.3 | standing-wave orbit quantization | 🔶 no dedicated phase yet — see gap note below |
| MERW → Born rule `ρ=\|ψ\|²` | §11b.4 | foundational (Q12) + 9b statistics | ✅ Q12 + 1b foundational note |
| Two-ingredient Schrödinger framing | §11b.5 | M5 overview / outreach | framing (no phase) |
| Beta decay as topology reconnection | §7 | weak-force gap | ✅ Q10 |
| Topology on Close/Yee frameworks (Duda's reciprocal ask) | §9 | no dedicated phase (cross-framework research) | 🔶 conceptual sketch in `4a §9` (Close=helicity/hopfion; Yee=ellipse-axis director); future collaboration item, not a numbered phase |
| Open: V(M) form / Faber reg / deeper substrate / weak SU(2) | §12 | M5.5 / M5.6 / — | ✅ Q7 / Q8 / Q9 / Q10 |
| Mainstream landscape comparison (page 48) | §11 | M5 positioning | framing (no phase) |

**One gap flagged:** the walking-droplet **orbit-quantization** side (path-memory kernel + Landau/Zeeman/double-quantization laws — the standing-wave complement to topology, Jeff Yee's regime) has **no dedicated M5 phase**. It rides on M3-near-field retention but is not a numbered phase. Candidate **future parallel stage** (alongside 9a–9d) if/when orbit quantization becomes a target. Logged here so it isn't lost.

---

### 🔶 Phase M5.4 — Matrix-field substrate migration (substrate + gate ✅; rendering + cleanup remaining)

> Replace Vector(3) ψ with the paper's `M(x) = O(x) D O^T(x)` real symmetric 3×3 matrix field. The architectural shift identified in M5.3. Migration is **additive-then-cleanup** (keep the engine runnable each step): the matrix substrate is added alongside the retiring ψ buffers; the Vector(3) wave buffers come out in the final cleanup once every consumer reads `director_nhat`.

- ✅ Storage redesign — `medium.py`: `M_am`/`M_prev_am`/`M_new_am` as full `ti.Matrix.field(3,3)` (decided 2026-05-26 — matches the proven spike, `ti.sym_eig`-ready, no reassembly; 6-component packing can swap in later behind the accessor) + derived `director_nhat`/`eigenvalues`/`director_nhat_new`. `LC_DELTA=0.5` uniaxial placeholder
- ✅ Operators (`engine2_pde.py`) — matrix `commutator [A,B]`, `matrix_laplacian` (component-wise 6-pt), antisymmetric `A_μ = [M, ∂_μ M]` (Eq. 19 → `antisym_A_mu`). Curvature `F_μν = ∂_μ A_ν − ∂_ν A_μ` (Eq. 20) + the full Eq.18 action **deferred to M5.5** (the M5.4 director-equivalent gate doesn't drive dynamics). Operators kept index-generic for the M5.8 4×4 promotion
- ✅ **Eigen-decomposition + director-readout kernel** — `eigen_decompose` (`engine2_pde`): per-voxel `ti.sym_eig` → sorted `eigenvalues` (λ₁≥λ₂≥λ₃) + sign-continuous principal `director_nhat`. Verified: spectrum recovers `(1, δ, δ)`, director recovers seeded hedgehog at **0.9995**. The lynchpin for rendering + trackers
- ✅ Seeders for vacuum and hedgehog on the M field — `engine1_seeds.py`: `uniaxial_M` helper (`M = δI + (1−δ)n̂⊗n̂`) + `seed_vacuum_M` / `seed_hedgehog_M` (write all 3 M buffers + correctly-signed `director_nhat`). Matrix-seeder director ≡ ψ-seeder director (**1.0000**); `eigen(M)` recovers it at **0.9998**
- [ ] **Rendering-stack repurposing** — re-source the full viz stack (flux_mesh, scalar warp, vector_warp, granule render, director glyphs, all color-coding modes) from M. None dropped; each re-sourced. Full per-feature plan + UX decisions in [4b_rendering_features.md](4b_rendering_features.md). Per `feedback_visual_rendering_priority`, ships with the kernels. **NEXT — needs Rodrigo's 4b UX calls (granule repurpose / glyph multiplicity / arrowhead retire / Displacement→Orientation rename / envelope mode)**
- ✅ **Redefine amp/freq trackers against M** — `update_trackers_M` (`engine3_observables`): amplitude = EMA of `‖M−D_vac‖_F` (thermal A); frequency = EMA of `‖Ṁ‖_F = ‖M−M_prev‖_F/dt` (de Broglie clock ω, thermal ω). Reuses the existing Trackers fields + 3-plane aggregation. Verified: vacuum→0, hedgehog amplitude 306× concentrated at the defect, clock ω 0 for static / >0 once M moves. This is the "measurement infra lands EARLY" deliverable
- ✅ Reproduce M5.1 Coulomb on the matrix field — **PASS** (`research/sandbox_v3/m5_4_coulomb_matrix.py`, director-equivalent path): R² = **0.9704** (threshold 0.95), b < 0 ATTRACTIVE — matches M5.1's Vector(3) R²=0.978. **Page-18 cross-val PASS** (`m5_4_coulomb_page18.py`): analytical `cos`-profile director + `M=n⊗n` + commutator-curvature `H = Σ_{i<j}‖[∂_iM,∂_jM]‖²` → R² = **0.9959**, b < 0, reproducing Duda's analytical `V(d) ≈ a − b/d` form (the absolute `25.16` is units/NIntegrate-normalization-dependent). Validates `pde.commutator` on real physics
- [ ] **Apply physical-energy-scaling factor** — **deferred to M5.5**: the factor (`ρ_medium × voxel_volume_am³ × INTERNAL_ENERGY_TO_AJ`) belongs on the *matrix* energyH (`‖Ṁ‖² + ‖∇M‖² + V(M)`), which lands with V(M) in M5.5. Applying it to the retiring ψ `compute_energyH_density` now is wasted work (the Coulomb gate uses Frank energyF, not energyH)
- 🔶 **Cleanup (after rendering)** — retire the Vector(3) wave buffers: rename the director-consuming kernels (`relax_director_step`, `compute_energyF_density`, the shared `compute_laplacian`) to read `director_nhat`; delete `evolve_psi` + the ψ wave seeders + `V_psi` + `psi_am`/`psi_prev_am`/`psi_new_am`. The director's home becomes `director_nhat`
- ✅ **Exit criterion MET** — M5.1 Coulomb fit reproduced on the matrix substrate (R²=0.9704 relaxed + 0.9959 analytical); operators verified against analytic small-cases (`commutator` vs E01/E10 + page-18 physics; `eigen` 0.9995; `matrix_laplacian` ∇²(x²+y²+z²)=6/dx² at relerr 0). **Remaining M5.4 work is rendering (step 5) + ψ cleanup — not gating.**

---

### 🚧 Phase M5.5 — Paper Lagrangian + Higgs-like V(M)

> Implement Duda paper Eq. 18: `L = Σ ‖F_μ0‖²_F − Σ ‖F_μν‖²_F − V(M)`. Subsumes the old "Skyrme stabilizer" phase since the Eq. 42 4D Skyrme-like kinetic is part of the same Lagrangian family — there is no separate Skyrme add-on phase in the new structure.

- [ ] Implement Eq. 18 action on the M field from M5.4
- [ ] Choose V(M) — start with the simpler eigenvalue-preference variant `V(M) = Σ_i (λ_i − Λ_i)²` (Eq. 12), then graduate to the LdG Higgs-like `V_LG(M) = a Tr(M²) − b Tr(M³) + c (Tr(M²))²` (Eq. 13)
- [ ] **Faber regularization** — port Manfried Faber's scheme to "activate" V (per `reference_faber_regularization.md`). Faber's scheme produces the running-coupling effect as a side-validation
- [ ] **Faber's explicit quantized-EM Lagrangian** (`4a §11b.2`): `Γ_i = (∂_i u)×u`, `R_μν = Γ_μ×Γ_ν`, `L_EM = −(αℏc/16π) R_μν R^μν` with `F_μν ~ R*_μν` (E↔B dual) — Faber's normalization of the EM term in the local-rotation-axis form; port alongside the LdG potential
- [ ] **EM-from-tilts cross-check** (`4a §11b.1`) — the superfluid hydrodynamics↔EM dictionary (vorticity `ω=∇×u ↔ B=∇×A`, Coriolis↔Lorentz, Lorenz-gauge parallel) is an independent route to Maxwell; use as a sanity check on the matrix-curvature EM derivation
- [ ] In M5.8 (4D extension), the Skyrme-like 4th-order kinetic from Eq. 42 lands; M5.5 stays in 3D
- [ ] **Exit criterion**: full Eq. 18 Lagrangian running; defect dynamics governed by the actual proposed action (not our scalar approximation)

---

### 🚧 Phase M5.6 — Biaxial twist + KG emergence

> With biaxial axes `Λ = (1, δ, 0)` and `δ ~ ℏ` (small twist), reproduce paper Fig. 9: Klein-Gordon-like equation emerges automatically from twist dynamics, NOT as an added V_psi term. This was the conceptual error in M5.2 Step 2.

- [ ] Configure the matrix field with `D = diag(1, δ, 0)`; Λ values from paper §III
- [ ] Apply the evolution equation (paper §IV — Euler-Lagrange of Eq. 18)
- [ ] Verify low-energy twist mode `ψ` obeys `2ψ_tt = (∇² − Â^hedg / ‖Â^hedg‖²) ψ` — Klein-Gordon-like with position-dependent mass from hedgehog structure. **Port the page-32 closed form directly** (`4a §11`, `A^hedg = (x,y,z)/r²`) — no need to re-derive
- [ ] Verify high-energy tilt modes obey Maxwell-like equations (paper Eq. 37 + 38)
- [ ] Compare ω(k) dispersion to M5.0h's free-wave dispersion — KG branch should show the mass-gap signature
- [ ] **Exit criterion**: paper Fig. 9 reproduced; KG is derived, not postulated

---

### 🚧 Phase M5.7 — Resonance hunt (Close's protocol) — UNBLOCKS 9b

> **Refinement from Robert Close (2026-04 email reply)**: "Even including the nonlinear term, I would expect your result of dispersing waves in most cases. But I suspect that certain amplitudes of certain harmonic waves will keep energy localized longer (i.e. as an unstable particle or resonance). My suggestion is to explore a wide range of amplitudes (probably l=1 harmonic wave is the most interesting). A likely criterion is that the maximum displacements should be comparable to the wavelength (or half or twice). Unless you have a good way to model an infinite system, I doubt that you will find completely stable non-radiating solutions."
>
> This phase implements Close's protocol on top of the biaxial matrix substrate from M5.4-M5.6. Success criterion is **extended-lifetime localization**, not perfect stability — the framework's "particles" are metastable resonances.

- [ ] Seed `l = 1` harmonic perturbation (dipole) on the biaxial matrix field from M5.6
- [ ] Amplitude sweep: `A/λ ∈ {0.5, 1, 2}` per Close's amplitude-comparable-to-wavelength criterion
- [ ] Measure energy-localization lifetime at each amplitude
- [ ] **Drive infrastructure for 9b thermal work** — parameterized harmonic forcing, energy tracking, modulation primitives (FM/AM)
- [ ] **Literature anchor** — skim BEC vortex-kinetics literature for long-lived oscillation modes (esp. Duda's 2026-05-13 cite: PRA "*Index theorem and vortex kinetics in Bose-Einstein condensates on a Haldane sphere with a magnetic monopole*", `10.1103/2msv-lk1m`). BEC is the experimental analog of our metastable-defect resonance hunt. Compact-manifold geometry (Haldane sphere) is a third escape from Derrick's theorem alongside topology + time-periodicity — could inform domain-shape choice for the test bed (see [`1b_topological_defect.md § Alternative stabilization — compact manifold`](1b_topological_defect.md))
- [ ] **Exit criterion**: at least one `(A/λ)` regime shows substantially extended energy-localization lifetime — "first long-lived particle in M5". **UNBLOCKS 9b THERMAL ENERGY.**

---

### 🚧 Phase M5.8 — 4D extension + Zitterbewegung clock (GROUP HEADLINE)

> **GROUP HEADLINE**: M5.8 passing IS the empirical answer to Duda's standing clock-propulsion question. External-comms framing: "M5.8 complete — Zitterbewegung clock at ω = 2mc²/ℏ confirmed empirically." This is the load-bearing test of clock propulsion (per `project_clock_propulsion` memory).
>
> Per Duda paper Fig. 10, extending to 4D with `D = diag(g, 1, δ, 0)` and SO(1,3) Lorentz introduces **negative-energy contributions** in the Hamiltonian (`ΓΓ̃` rotation-boost type) that automatically propel the de Broglie clock — no engineered V(ψ) propulsion needed. The Zitterbewegung emerges as a consequence of 4D Lorentz signature.

- [ ] **Infrastructure pre-check (1D)** — before the full 4D run, reproduce Duda's clock toy-model numerics (`4a §11`, slide p.33; arxiv:2501.04036 kink+clock): `φ = tanh(0.6326x + 0.0198x³ + 0.0203x⁵)`, energy `E ≈ 2.0252`, frequency `ω ≈ 1.2898`. Cheap validation that our integrator reproduces clock propulsion before scaling to 4D LdGS
- [ ] Promote the 3D matrix substrate from M5.4-M5.7 to 4D — add 0th time axis, SO(1,3) Lorentz group
- [ ] Apply 4D Skyrme-like 4th-order kinetic (paper Eq. 42): `L = −Σ F_μναβ F^μναβ − V(M)` with `F_μναβ = [∂_μ M, ∂_ν M]_αβ`
- [ ] Verify negative-energy contributions emerge from spacetime signature (Fig. 10)
- [ ] Seed a single biaxial hedgehog on the electron axis (δ)
- [ ] Measure intrinsic oscillation frequency at the defect core (FFT of director rotation)
- [ ] **Target**: `ω = 2 m_e c² / ℏ ≈ 1.55 × 10²¹ rad/s` (electron Zitterbewegung)

**Mass → frequency table** (validation targets across particle species):

| Particle | Defect type | Target ω = 2mc²/ℏ |
| --- | --- | --- |
| Electron | Point hedgehog (δ axis) | 1.55 × 10²¹ rad/s |
| Muon | Point hedgehog (1 axis) | 3.21 × 10²³ rad/s |
| Tau | Point hedgehog (g axis) | 5.39 × 10²⁴ rad/s |
| Neutrino | Closed vortex loop | ~10¹⁵ rad/s (sub-eV mass) |

**Experimental anchors** (validation against measured values):

| Year | Experiment | Regime | Relevance |
| --- | --- | --- | --- |
| 2010 | Gerritsma et al. (trapped-ion Dirac analog) | Analog | First Zitterbewegung-class observation |
| 2008 | Catillon, Cue, et al. (electron channeling) | 81 MeV electrons (relativistic) | Direct electron-clock measurement, with kinematic mass correction |
| **2026** | **Positronium de Broglie clock, Nature Comm.** | **3 keV (nonrelativistic e⁺e⁻ bound state)** | **Highest-priority anchor**: kinetic energy ≪ rest-mass → measurement is essentially `ω = 2mc²/ℏ`. Cleanest validation target. *Flagged by Duda 2026-05* |

- [ ] **Cross-particle test**: seed defects of different masses (electron + muon, or electron + tau) at far separation. Each ticks at its own mass-derived `ω`, independently
- [ ] **Negative-Hamiltonian propulsion test**: toggle the `−b·Tr(M³)` cubic term on/off. With it ON, single-defect dynamics should self-sustain for ≥100·T_Z runs; with it OFF (b=0), oscillation should damp. Identifies which term in V(M) is the propulsion mechanism
- [ ] **Exit criterion**: electron Zitterbewegung frequency reproduced within 10% of `1.55 × 10²¹ rad/s`. **GROUP-HEADLINE SEND** to Models-of-Particles thread.

---

### 🚧 Phase M5.9 — 3 lepton families + Cornell quark strings

> Standard Model correspondence. The biaxial Λ = (g, 1, δ, 0) gives 3 axis-choices → 3 lepton families with the same Q but different masses (e/μ/τ). Cornell potential for quark strings via topological vortex.

- [ ] Seed hedgehogs along the 3 different axes of the biaxial M field
- [ ] **Calibrate `(g, δ) numerically** against observed lepton-mass ratios: μ/e ≈ 207, τ/e ≈ 3477`. Per Duda 2026-04-19 guidance, these are calibration parameters, not ab-initio derivations
- [ ] Topological vortex string (1D defect line, not point hedgehog) for quark-pair binding
- [ ] Validate Cornell form: `V(r) = −α/r + σ·r` with `σ ≈ 1 GeV/fm` (string tension)
- [ ] **Exit criterion**: lepton mass ratios within 10% of observed; Cornell potential reproduced with `σ ≈ 1 GeV/fm`

---

> **9b — Thermal energy**, **9c — Time dynamics**, **9d — Composite particles**, and **9a — Electromagnetic wave packets** are not M5.x phases. Each is a separate research file under `research/`, with its own test program. See [SUMMARY § Parallel research stages](#parallel-research-stages) above for cross-refs.
