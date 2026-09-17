# M4 Pipeline Engine — Implementation Roadmap

> **Status:** DRAFT (working document)
> **Scope:** The pipeline_engine as the realisation substrate for Enhanced EWT.
> **Audience:** Contributors implementing, testing, or extending M4.
> **Related:** `M4_engine_upgrade.md`, `M4_k_selectivity_Formalization.md`,
> `__M4_model_briefing.md`, manuscript v5.0.x (Zenodo), Yee's EWT corpus.

---

## 0. Executive Summary

The pipeline_engine is a **composition root** for M4 simulations: a list of
stateless processors executed against a shared `Context`, with explicit
`requires`/`provides` contracts, swappable sinks, and a typed `FeatureBag`
for runtime state. It exists because the previous monolithic `wave_engine.py`
made it impossible to express the actual physics of Enhanced EWT — every
mechanism was a hard-coded branch, every variant a new file, every comparison
a code duplication.

This document captures:

1. **Why** the previous implementation only *simulated* EWT and did not
   *realise* it (Section 1).
2. **What** the core conceptual commitments are that the new engine must
   honour (Sections 2–5).
3. **How** to build the engine (Block 1) and the physics (Block 2), with
   concrete, checkable work items (Sections 6–7).
4. **What** from the old implementation is obsolete and should not be ported
   (Section 8).
5. **Which** questions remain genuinely open and require author input or
   further research (Section 9).

Nothing in this document is a claim about nature. It is a **plan for a tool**,
written so that the tool can express the hypotheses we intend to test.

---

## 1. Why the previous implementation faked EWT

This section is diagnostic, not accusatory. The old `wave_engine.py` was an
honest attempt to translate EWT into code. It failed for structural reasons,
not for lack of effort.

### 1.1. Wave centres as pinned sources, not reflectors

Yee's EWT defines a wave centre as a **point at which incoming waves are
reflected to become outgoing waves**. The soliton is the standing-wave
interference of `Ψ_in` and `Ψ_out`, bounded by the particle radius.

The old implementation defined wave centres as **hard-pinned regions**:

```text
ψ = A · sin(ω t + offset) · r̂    inside a ball of radius R around each WC
```

This is a driven antenna, not a reflector. The field does not participate in
the WC's existence; the WC simply overwrites voxels. There is no reflection,
no `Ψ_in`/`Ψ_out` distinction, no conservation argument.

### 1.2. One field, not two

Yee's EWT distinguishes **longitudinal** (mass, charge) from **transverse**
(spin, magnetism), coupled through the fine-structure constant at the WC.
The old implementation had a single vector field `ψ` with no mode split.

Consequence: the fine-structure constant had to be **imposed** through
`cos(offset)` and the mass of the electron through an analytic formula.
Neither emerged from dynamics.

### 1.3. Wavelength structure ignored

Yee's standing-wave geometry prescribes *decreasing* wavelengths from the core:

```text
r_wavelength(n) = 2Kλ − 2nλ
r_x = (K + 2·Σ_{n=1..x}(K−n)) · λ
```

with a maximum radius `r_particle = K²λ`. The old implementation used a
single `base_wavelength` throughout. The electron's characteristic ring
structure — visible in the Lund stroboscope image — was absent.

### 1.4. Density profile fixed, not self-consistent

Enhanced EWT's gravitational sector is built on the **push-out** mechanism:
the soliton's energy density `ρ_E` displaces Elastic Medium Constituents
(EMC), creating a local deficit `ρ(r) < N_ν,stat`. The deficit is what the
outside vacuum presses against; the pressure gradient is what we call gravity.

The old implementation had `ρ(r)` as a **fixed profile relative to the domain
centre**. When a WC drifted, the well did not follow. There was no
self-consistency between `|Ψ|²` and `ρ(r)`.

### 1.5. Nonlinearity as an external potential, not a consequence of `c(ρ)`

M4.9 established the microscopic relation:

```text
v_phys(η) = a(η) · sqrt(k(η)/m₀)  ∝  η^{+1/2}
```

where `η = ρ/ρ₀`. The local wave speed depends on the local EMC density.
Because the soliton itself depletes EMC, `c²(r) = c₀²·ρ(r)/ρ₀`, and the
wave equation becomes *automatically nonlinear*:

```text
∂²Ψ/∂t² = c²(ρ) ∇²Ψ = c₀² (1 − β|Ψ|²) ∇²Ψ
```

The old implementation instead injected a Klein–Gordon-style potential
`V(ψ) = (c₁/4)u² − (c₂/6)u³` with `u = |ψ|²`. Mathematically this gives an
NLS soliton. Physically it is a prosthesis: the coupling constants `c₁`, `c₂`
were fitted, not derived.

### 1.6. Gravity as `−∇E`, not as `−∇ρ`

The old `compute_force_vector` used `F = −∇E_local`, where `E_local` was a
heuristic `ρ·V·(f·A)²` with `f` hard-coded to the base frequency. This
measured the *internal* energy gradient, not the *external* EMC density
gradient. The push-out picture requires `F_pressure = −α∇ρ(r)`, where `ρ(r)`
is the EMC packing density, not the soliton's energy density.

### 1.7. Summary table

| EWT primitive | Old implementation | Why it is a fake |
|---|---|---|
| WC as reflector | WC as hard pin | No reflection, no in/out |
| `Ψ_in` + `Ψ_out` | Single `ψ` | No distinction, no closure |
| Decreasing `λ(n)` | One `base_wavelength` | Ring structure absent |
| `Ψ_long` + `Ψ_trans` | One vector field | Spin not emergent |
| `α` from `\|Ψ_out\|/\|Ψ_in\|` | `α` as input | Not derived |
| `ρ(r)` self-consistent | Fixed profile | No push-out feedback |
| `c²(ρ)` nonlinearity | External `V(ψ)` | Coupling fitted |
| `F = −∇ρ` gravity | `F = −∇E` | Measures wrong gradient |
| EMC Wall as boundary | Dirichlet `ψ = 0` | No profile, no sign change |
| BCC lattice | Simple cubic stencil | Medium ≠ medium of theory |

---

## 2. Conceptual commitments

These are the principles the new engine and physics must respect. They are
*not* claims about nature; they are constraints on the tool.

### 2.1. Conservation is structural, not imposed

Energy conservation in EWT is a consequence of the **unitarity of reflection
at the WC** and the **Hermiticity of the wave equation**, not an external
constraint. If the WC reflection satisfies

```text
|Ψ_out|² + |Ψ_spin|² = |Ψ_in|²
```

and the wave equation derives from a real Hamiltonian density, then energy is
conserved by construction. Any deviation is a bug, not a feature.

The engine must therefore expose enough structure to *measure* the energy
budget and *verify* that `dE/dt + flux ≈ 0`.

### 2.2. The soliton is an open system in steady state

The medium carries an always-on base wave (Yee: "waves flow through all of
matter"). This wave supplies `Ψ_in` to the WC. The WC reflects part of it as
`Ψ_out` and converts part of it into transverse spin. In steady state:

```text
flux_in = flux_out + flux_spin
```

The soliton is a **dissipative structure** in the Prigogine sense — a local
concentration of energy maintained by continuous exchange with the medium,
not a closed system isolated from it.

This is a deliberate commitment to **variant B** (coupled dynamics), not
variant A (static background). It is harder to implement, but it is what the
manuscript describes.

### 2.3. The EMC density field is dynamic and self-consistent

`ρ(r)` is not an input. It is a solution of

```text
∂ρ/∂t = D ∇²ρ − γ(ρ − ρ₀) − β |Ψ|²
```

in the simplest model, or a richer evolution with inertia. In steady state,
`ρ(r)` and `|Ψ(r)|²` are mutually consistent. The EMC Wall — the local peak
`ρ > ρ₀` that reverses the sign of the nonlinearity — is a **consequence**
of the dynamics, not a parameter.

### 2.4. Scales are separable

Three distinct length scales appear in the theory:

| Scale | Value (electron) | Role |
|---|---|---|
| Soliton core | `K²λ = 100 λ_ν ≈ r_e` | Standing-wave boundary |
| EMC deficit range | `r_e/α ≈ 137 r_e` | Where `ρ(r) → N_stat` |
| Compton wavelength | `λ_C ≈ 2π · 137 r_e` | Quantum scattering scale |

**The core is simulated. The far-field deficit is treated analytically.**
No single simulation resolves both scales; the connection is made at the
boundary of the simulated domain.

### 2.5. Units are pluggable

No processor contains a dimensional constant. Every constant comes from a
single `UnitSystem` feature. Three implementations are anticipated:

- `NaturalUnitSystem` — `λ = 1`, `c = 1`. For research and algorithm testing.
- `OpenWaveUnitSystem` — `am`, `rs`, as in the legacy xparameters.
- `SIUnitSystem` — for output conversion and cross-checking against CODATA.

The same pipeline runs under any of them. Results are comparable after
conversion.

### 2.6. Topology and geometry are hypotheses, not assumptions

The 1-3-6 arrangement and the `n·λ` spacing of wave centres are **candidate**
configurations, not given truths. The engine must allow systematic variation
of:

- Topology (1-3-6, golden-angle, BCC, line, random).
- Spacing (exact `n·λ`, exact `(n+½)λ`, perturbed, swept).
- Coupling (none, push-out without wall, push-out with wall).

Only then can we say which of these factors is necessary, sufficient, or
irrelevant for stability.

---

## 3. The unit system

### 3.1. Contract

A `UnitSystem` is a feature providing at least:

```text
c            : wave speed in these units
lambda       : fundamental wavelength (λ_ν)
dx           : grid step in these units
dt           : time step in these units
rho_0        : statutory background density (N_ν,stat)
A_pi         : 4π³ + π² + π
eps_M        : 1 / (N_geom · π³)
N_geom       : 8π⁴ · (1 − ζ)
gamma        : 1 / eps_M
X_eff        : geometric dilution factor
N_nu_eff     : effective volume deficit
```

Plus conversion methods:

```text
to_physical_length(x)   -> metres
to_physical_time(t)     -> seconds
to_physical_energy(E)   -> joules
to_physical_density(r)  -> 1/m³
```

### 3.2. Implementations

**NaturalUnitSystem** (default for research):
```text
c      = 1
lambda = 1
dx     = 1 / K_grid     (e.g. 0.05 for 20 voxels/λ)
dt     = 0.9 · dx / c   (CFL-safe)
```
Advantages: f32-safe (all values near 1), `gamma` dimensionless, fast.

**OpenWaveUnitSystem**:
```text
c      = 0.3 am/rs
lambda = EWAVE_LENGTH / ATTOMETER   (≈ 28.5 am)
dx     = lambda / 12
dt     = 0.9 · dx / c
```
Advantages: comparable to legacy xparameters.

**SIUnitSystem**:
```text
c      = 299792458 m/s
lambda = 2.8540965e-17 m
dx     = r_e / 200
dt     = ...
```
Used primarily for output conversion.

### 3.3. Consequences for processors

A processor must never write `c = 1` or `gamma = 2.4e4` or `lambda = 28.5`.
It must always write:

```python
units = ctx.data.require(UnitSystem)
c = units.c
gamma = units.gamma
```

This rule is checkable by inspection: grep for numeric literals in
`physics/*.py`. If a processor contains a dimensional constant, it is a bug.

---

## 4. Universe scaling strategy

### 4.1. The scale problem

For K = 10 (electron), the core radius is `r_core = K²λ = 100 λ_ν`. With
`dx = 0.05` (20 voxels per λ_ν), the core occupies 2000 voxels in radius.
A full 3D simulation of the core is `(4000)³ ≈ 6.4 × 10¹⁰` voxels. This
is not tractable.

### 4.2. Three strategies

**(a) Full 3D simulation of the core** — ideal but infeasible at K = 10.
May become feasible for K = 1 (neutrino) where `r_core = 1λ_ν`.

**(b) 1D radial simulation** — exploits spherical symmetry. `r_core = 100λ_ν`
becomes 2000 voxels in 1D. Loses all angular structure (spin, 1-3-6 topology),
but is an excellent validation tool for radial dynamics.

**(c) 3D simulation of a truncated domain** — simulate `r ≤ r_domain` with
`r_domain ~ 50 λ_ν`, and impose an **analytic outer boundary condition** that
represents the far-field `1/r` tail.

**Recommendation:** Start with (c). Use (b) for validation. Attempt (a) only
for K = 1 and maybe K = 2.

### 4.3. The outer boundary

The analytic far-field is:

```text
ρ(r) → N_stat · (1 − r_s/r)   for r ≫ r_core
```

where `r_s = 2 G M / c²` is the gravitational radius. The boundary condition
at `r = r_domain` imposes this tail; it is not a Dirichlet `ψ = 0` and it is
not a reflecting wall.

Two additional boundary types are available for testing:

- **Absorbing (Sommerfeld)** — for isolating solitons from their environment.
- **Periodic** — for homogeneous-medium tests without boundaries.

### 4.4. The EMC wall

The EMC wall at `r_wall ~ r_e/α ≈ 137 r_e` is **outside the simulated
domain** in strategy (c). It enters only as a **modification of the analytic
boundary condition**, not as a simulated field. Its role — reversing the sign
of the nonlinearity to prevent outward leakage — is realised by the boundary,
not by a local peak in `ρ(r)`.

If later simulations extend the domain, the wall can be promoted to a
simulated feature without changing the contract.

---

## 5. Vacuum energy strategy

### 5.1. What the base wave is

The base wave is the always-on ground-state oscillation of the medium. In
natural units:

```text
Ψ_base(r, t) = A_0 · cos(k·r − ω·t)
```

with `k = 2π/λ = 2π` and `ω = 2π·c/λ = 2π`.

### 5.2. Is it a source or a carrier?

The base wave **carries** energy (it is a non-zero field), but it does not
**supply** energy in the sense of a reservoir that can be drawn down. In
steady state, the flux it carries in equals the flux it carries out of any
closed surface, so its net contribution to the energy budget of a soliton is
zero.

Implementation consequence: `E_base` is **excluded** from the conservation
check. The check is

```text
d/dt (E_soliton + E_deformation) + flux_through_boundary = 0
```

where `E_soliton = ∫ |∇Ψ_soliton|² dV` and `E_deformation = ∫ (ρ − ρ₀)² dV`.

### 5.3. Coupling to the soliton

The soliton interacts with the base wave only **at the WC**, through
reflection. There is no volume coupling. This preserves the locality of
the interaction and avoids the "pumping" problem: a volume-coupled base
wave would continuously inject energy everywhere, and the soliton would
either grow without bound or require an explicit damping term.

### 5.4. Open question

Whether the base wave is truly non-dissipative, or whether it slowly relaxes
toward an equilibrium (which would break conservation at cosmological
timescales), is an **author-gated question** and is not resolved by this
document. For in-simulation purposes, the base wave is treated as steady.

---

## 6. Block 1 — Engine implementation

Infrastructure only. No specific physics. Each item is a work unit.

### 1.0 — `UnitSystem` feature

- [ ] Define the `UnitSystem` protocol (fields + conversion methods).
- [ ] Implement `NaturalUnitSystem`, `OpenWaveUnitSystem`, `SIUnitSystem`.
- [ ] Add `UnitSystem` to `Pipeline.external_provides`.
- [ ] Wire `Runner.run(..., initial_features=[units])`.
- [ ] Add a linter rule: no dimensional literal in `physics/*.py`.
- [ ] Document in `README.md`.

### 1.1 — Multi-field FeatureBag

- [ ] Define `PsiBaseField`, `PsiLongField`, `PsiTransField`,
      `EMCDensityField`, `EMCFluxField` (optional) as separate features.
- [ ] Each field is a `ti.Vector.field(3, ti.f32, shape=grid)` with a
      triple-buffer variant for time integration.
- [ ] Ensure `FeatureBag` handles these as distinct types without aliasing.
- [ ] Add a test: allocate all five, write distinct values, read back.

### 1.2 — Trackers per voxel

- [ ] Define `TrackerFields` feature.
- [ ] Fields: `energy_long_local`, `energy_trans_local`, `amp_local`,
      `freq_local`, `rho_local`, `wall_proximity`.
- [ ] Implement `TrackersUpdate` processor in `Stage.MEASURE`.
- [ ] Implement 3-plane sampling (as in the legacy `sample_avg_trackers`)
      to compute global averages without full reductions.
- [ ] Add a test: uniform field → uniform tracker values.

### 1.3 — Multi-field evolution

- [ ] Generalise `LaplacianProcessor` to accept a field name in `__init__`.
- [ ] Generalise `LeapfrogProcessor` similarly.
- [ ] Define a coupling contract: `Ψ_long → Ψ_trans` conversion at WCs,
      with a tunable coefficient.
- [ ] Add a test: two coupled fields, no coupling coefficient → independent
      evolution; with coefficient → energy transfers.

### 1.4 — Source terms (additive, not overwrite)

- [ ] Define the contract: a source term **adds** to `psi_new`, never
      overwrites `psi_am`.
- [ ] Implement `SeedBaseWave` — seeds `Ψ_base` once.
- [ ] Implement `SourceTermInterface` — base class for additive sources.
- [ ] Add a test: two sources, superposition holds.

### 1.5 — Reflector interface

- [ ] Extend `WCState` with `reflect_coeff_long`, `reflect_coeff_trans`,
      `phase_shift`.
- [ ] These are *attributes*, not yet *behaviours*. The values are set by
      the experiment, not computed.
- [ ] Document the contract clearly: a reflector is a WC that satisfies
      unitarity on `Ψ_in`, `Ψ_out`, `Ψ_spin`.
- [ ] Add a test: reflection of a plane wave from a single reflector
      preserves energy.

### 1.6 — WC motion (drift)

- [ ] Extend `WCState` with `velocity`, `force`.
- [ ] Implement `WCMotionProcessor` in `Stage.POST_UPDATE`.
- [ ] Implement `WCDriftRule` contract: a callable that returns a force
      vector given local field values.
- [ ] Provide a default rule (`F = −∇ρ`) and a no-op (`F = 0`).
- [ ] Add a test: no-op drift → positions unchanged; default drift on a
      static `ρ` → WCs move to minimum.

### 1.7 — EMC Wall (analytic boundary)

- [ ] Define `EMCBoundary` feature: `r_wall`, `wall_profile`,
      `boundary_condition`.
- [ ] Implement `EMCBoundaryProcessor` in `Stage.POST_UPDATE`.
- [ ] Boundary condition modifies the far-field `1/r` tail. Not Dirichlet.
- [ ] Add a test: wave hitting the boundary is partially reflected,
      partially transmitted, energy conserved.

### 1.8 — Energy budget tracker

- [ ] Define `EnergyBudget` feature.
- [ ] Fields: `E_long`, `E_trans`, `E_deform`, `flux_boundary`, `dE_dt`.
- [ ] Implement `EnergyBudgetUpdate` in `Stage.MEASURE`.
- [ ] Add a test: harmonic oscillator in 1D → `dE/dt ≈ 0` to machine
      precision.

### 1.9 — Stability metrics

- [ ] Define `StabilityMetrics` feature.
- [ ] Fields: `sol_lifetime`, `localization`, `sphericity`, `freq_drift`,
      `wc_drift`.
- [ ] Implement `StabilityMetricsUpdate` in `Stage.MEASURE`.
- [ ] Add a test: static Gaussian → `localization` stays constant;
      spreading Gaussian → `localization` decreases.

### 1.10 — Experiment runner

- [ ] Implement `ExperimentRunner`:
      takes a list of configuration dicts, runs each, records summary.
- [ ] Output: CSV with one row per run, columns = config + final metrics.
- [ ] Support deterministic seeds per run.
- [ ] Add a test: 3-run sweep produces 3-row CSV.

### 1.11 — Geometric constants provider

- [ ] Define `GeometricConstants` feature.
- [ ] Implement `ComputeGeometry` lifecycle processor that imports from
      `m4_7_ewt_emergence_engine.py` and populates the feature.
- [ ] Wire as `external_provides` for pipelines that need it.
- [ ] Add a test: computed `A_pi`, `eps_M`, `N_geom` match the engine.

### 1.12 — Units and conversions

This is folded into 1.0. Listed separately only for traceability.

### 1.13 — Checkpoint / restart

- [ ] Define `CheckpointProcessor` (lifecycle + periodic).
- [ ] Serialise `FeatureBag` to disk (Taichi fields → numpy → npz).
- [ ] Deserialise on startup.
- [ ] Add a test: run 100 steps, checkpoint, run 100 more;
      run 200 from scratch; compare.

### 1.14 — Live monitor

- [ ] Adapt `live_monitor_viewer.py` for pipeline_engine.
- [ ] Panels: energy (long/trans/emc), boundary flux, WC drift, freq.
- [ ] Reads `live.json` written by `LiveJsonSink`.
- [ ] Add a test: launch monitor, run 100 steps, monitor updates.

### 1.15 — Research logging schema

- [ ] Define the output layout:
      `run_meta.json`, `timeseries.parquet` (or CSV), `summary.json`,
      `events.json`.
- [ ] `run_meta.json`: config, seed, code hash, unit system.
- [ ] `timeseries`: full time series per metric, with configurable cadence.
- [ ] `summary`: final metrics for sweep aggregation.
- [ ] `events`: annihilation, boundary hits, instability detected.
- [ ] Add a test: schema validates against a JSON schema.

### 1.16 — Universe scaling strategy

This is the decision documented in Section 4. Work items:

- [ ] Add `r_domain` and `r_core` to `UnitSystem`.
- [ ] Implement the analytic boundary condition (Section 4.3).
- [ ] Document the choice: strategy (c) as default.
- [ ] Add a test: for K = 1, simulation of the full core at 20 voxels/λ.

### 1.17 — Vacuum energy strategy

This is the decision documented in Section 5. Work items:

- [ ] Add `A_base` (base wave amplitude) to `UnitSystem`.
- [ ] Implement `SeedBaseWave` (in 1.4).
- [ ] Implement the conservation check **excluding** `E_base`.
- [ ] Add a test: base wave alone → `dE_soliton/dt = 0` trivially.

### 1.18 — Diagnostic hooks

- [ ] Define a stop-condition contract: `StopCondition` callable.
- [ ] Implement `DiagnosticProcessor` in `Stage.MEASURE`.
- [ ] Built-in conditions: `dE/dt > threshold`, `localization < threshold`,
      `sphericity < threshold`.
- [ ] Add a test: run with a forced violation → simulation stops early.

### 1.19 — Deterministic seeds

- [ ] Add `seed` to `RunContext` (already present).
- [ ] All random initialisation reads from `ctx.run.seed`.
- [ ] Add a test: two runs with same seed → bit-identical output.

### 1.20 — Parameter sweep DSL

- [ ] Define YAML/JSON schema for sweeps:
      `topologies`, `spacings`, `couplings`, `K`, `seeds`.
- [ ] Implement `SweepRunner` that consumes the schema and drives
      `ExperimentRunner`.
- [ ] Add a test: 2×2 sweep produces 4 runs.

### 1.21 — Artifact versioning

- [ ] Hash the configuration + code state.
- [ ] Store results in `output_dir / <hash> /`.
- [ ] Provide `list_runs()` and `load_run(hash)` utilities.
- [ ] Add a test: same config → same hash; different config → different.

---

## 7. Block 2 — Physics implementation

Each item is a *variant*. The engine (Block 1) is variant-agnostic; the
physics layer supplies the specific mechanisms.

### 2.0 — Soliton assembly

Before any specific mechanism, define how the pieces compose:

```text
Ψ_total = Ψ_base + Ψ_soliton
ρ(r)    = ρ₀ − β |Ψ_soliton|²        (initial guess)
c²(r)   = c₀² · ρ(r) / ρ₀
```

The soliton exists as a *fixed point* of the coupled system:
`Ψ_soliton` and `ρ` are mutually consistent. This composition is the
foundation; all variants below are refinements.

- [ ] Implement `SolitonAssembly` as a documented contract.
- [ ] Add a test: static `Ψ_soliton`, no coupling → no soliton, only
      dispersion.

### 2.1 — Base wave

- [ ] **B1a**: static base wave (`A_base` constant in time).
- [ ] **B1b**: oscillating base wave (time-dependent phase).
- [ ] **B1c**: base wave coupled to soliton (feedback).
- [ ] Recommended start: **B1b** — closest to EWT and analytically tractable.

### 2.2 — WC as reflector

- [ ] **B2a**: perfect reflection, no spin conversion.
- [ ] **B2b**: reflection with `α` conversion (`|Ψ_spin|² = α|Ψ_in|²`).
- [ ] **B2c**: geometry-dependent reflection (local `α`).
- [ ] Recommended start: **B2b** — this is what makes `α` emergent.

### 2.3 — Longitudinal ↔ transverse coupling

- [ ] **B3a**: conversion at WCs only.
- [ ] **B3b**: volume conversion proportional to `|Ψ_long|²`.
- [ ] **B3c**: with relaxation (`Ψ_trans → Ψ_long` possible).
- [ ] Recommended start: **B3a** — local, clean, matches Yee.

### 2.4 — EMC density dynamics

- [ ] **B4a**: instantaneous (`ρ = ρ₀ − β|Ψ|²`).
- [ ] **B4b**: relaxation dynamics (`∂ρ/∂t = D∇²ρ − γ(ρ−ρ₀) − β|Ψ|²`).
- [ ] **B4c**: inertial dynamics (full wave equation for `ρ`).
- [ ] Recommended start: **B4a** for tests; promote to **B4b** for
      self-consistent solitons.

### 2.5 — Wave speed modulation

- [ ] **B5a**: `c²(ρ) = c₀² · ρ/ρ₀`.
- [ ] **B5b**: power law `c² = c₀² · (ρ/ρ₀)^n`.
- [ ] **B5c**: anisotropic `c(ρ, ∇ρ)`.
- [ ] Recommended start: **B5a** — matches M4.9.

### 2.6 — Density-modulated nonlinearity

- [ ] **B6a**: `F = γ_nl · (1 − ρ/ρ₀) · |Ψ|² · Ψ`.
- [ ] **B6b**: explicit profile `mod(r)` instead of local `ρ`.
- [ ] **B6c**: nonlinearity in `c²(ρ)` instead of in `F`.
- [ ] Recommended start: **B6a** — matches manuscript Variant B.

### 2.7 — EMC Wall

- [ ] **B7a**: wall without peak (`ρ_wall < ρ₀`).
- [ ] **B7b**: wall with peak (`ρ_wall > ρ₀`).
- [ ] **B7c**: wall coupled to `|Ψ|²`.
- [ ] Recommended start: **B7b**, but note that at strategy (c) the wall
      is analytic (Section 4.4), so this is a *boundary condition* variant.

### 2.8 — WC motion rule

- [ ] **B8a**: `F = −∇ρ` (EMC density gradient).
- [ ] **B8b**: `F = −∇|Ψ|²` (energy gradient).
- [ ] **B8c**: `F = −∇(ρ + |Ψ|²)`.
- [ ] **B8d**: `F = 0` (control).
- [ ] Recommended start: **B8a** — matches push-out.

### 2.9 — WC topology

- [ ] **B9a**: `tetrahedron_10_locked` (r1 = 1λ, r2 = 2λ).
- [ ] **B9b**: `tetrahedron_10_unlocked` (legacy r1, r2).
- [ ] **B9c**: `golden_angle`.
- [ ] **B9d**: `bcc_lattice`.
- [ ] **B9e**: `line` (negative control).
- [ ] **B9f**: `random` (negative control).

### 2.10 — WC spacing

- [ ] **B10a**: sweep spacing at fixed topology.
- [ ] **B10b**: `n·λ` vs `(n+½)·λ`.
- [ ] **B10c**: perturbation ±10%, ±20%.

### 2.11 — K-selectivity

- [ ] **B11a**: sweep K = 2..12 at fixed topology, spacing, coupling.
- [ ] **B11b**: K × topology sweep.
- [ ] **B11c**: K × spacing sweep.
- [ ] **B11d**: full sweep.

### 2.12 — Energy conservation verification

- [ ] **B12a**: measure `dE/dt` for isolated soliton.
- [ ] **B12b**: measure boundary flux.
- [ ] **B12c**: compare stable vs unstable K.

### 2.13 — Emergent `α`

- [ ] **B13a**: measure `|Ψ_trans|² / |Ψ_long|²` near WCs.
- [ ] **B13b**: compare with `1/(8π⁷(1−ζ))`.
- [ ] **B13c**: verify `α` does not drift in time.

---

## 8. What is dropped from the old implementation

The following should **not** be ported. Each is listed with the reason.

| Old element | Reason for dropping |
|---|---|
| `interact_wc_dirichlet` | Hard pin, superseded by reflector (1.5, 2.2) |
| `interact_wc_neumann` | Hard pin variant, superseded |
| `interact_wc_soft` | Additive pin, at most one source variant, not the mechanism |
| `V_MODE = 1` (pure cubic) | Prosthesis, superseded by 2.6 |
| `V_MODE = 2` (quintic) | Prosthesis, superseded by 2.6 |
| `V_MODE = 3` (double-well) | Not relevant to EWT |
| `V_MODE = 4,5,6,7,9,10` | Simplified profiles, superseded by 2.4 |
| `energy_local_aJ` with hardcoded `base_frequency` | Superseded by `EnergyBudget` |
| `compute_force_vector` (`F = −∇E`) | Superseded by `F = −∇ρ` (2.8) |
| `DirichletBoundaryProcessor` (`ψ = 0`) | Superseded by analytic EMC boundary (1.7) |
| Simple cubic Laplacian | To be replaced by BCC stencil if needed |
| `seed_wave` modes 0, 1 | Kept as utilities, not central |
| `detect_annihilation` | Deferred until reflectors work |
| `select_voxels` | Already removed |

**Kept from the old implementation** (as utilities, not as core):

- Idea of flux mesh (rendering, Block 3).
- Idea of granule motion (rendering, Block 3).
- `constants.EWAVE_*` (absorbed into `OpenWaveUnitSystem`).
- `m4_7_ewt_emergence_engine.py` formulas (absorbed into `GeometricConstants`).

---

## 9. Open questions

These require author input or further research. They are **not** to be
resolved by inference.

### Q1. Is `r_wall = 137 r_e` a size or a range?

The manuscript says the deficit extends to `r_e/α ≈ 137 r_e`. This is
interpreted as a *range* — the distance over which `ρ(r) → N_stat` — not as
the size of the electron. Confirmation needed.

### Q2. Is the base wave truly steady?

If the base wave slowly relaxes (cosmological timescales), conservation is
only approximate in-simulation. If it is genuinely steady, conservation is
exact. The manuscript does not resolve this.

### Q3. Is the EMC wall a physical peak or an effective boundary?

At strategy (c), the wall is *outside* the simulated domain and enters as a
boundary condition. If the domain is later extended, the wall becomes a
simulated field. Which interpretation is preferred?

### Q4. What is the correct unit system for research?

Natural units are recommended for tractability. But the manuscript's
predictions are in SI. The choice affects how `γ` enters the dynamics.

### Q5. What is the correct definition of "stability"?

Several candidates: lifetime, localization, sphericity, frequency stability.
Which is primary? Or are all co-primary?

### Q6. Should the WC motion be continuous or discrete?

Yee's picture suggests continuous drift toward amplitude minima. But the
"lock-in" language suggests discrete jumps between nodes. Which is correct?

### Q7. Is `K = 10` a topological necessity or an energetic optimum?

The manuscript argues topology (winding number `Q = 10`). But if `Q` is not
conserved in the simulation (e.g., under continuous deformations that
change the count of wave centres), then `K = 10` must be selected
energetically. Which mechanism is operative?

### Q8. Does the soliton require spin to be stable?

Yee's picture includes spin from the start. If spin is required, then
`Ψ_trans` must be seeded along with `Ψ_long`. If not, the soliton can be
purely longitudinal. Which is correct?

---

## 10. Recommended execution order

**Phase A — Engine foundation (Block 1.0–1.3)**

1.0 (UnitSystem) → 1.1 (Multi-field) → 1.2 (Trackers) → 1.3 (Multi-field
evolution)

**Phase B — Physics interfaces (Block 1.4–1.7)**

1.4 (Source terms) → 1.5 (Reflector interface) → 1.6 (WC motion) → 1.7 (EMC
boundary)

**Phase C — Measurement (Block 1.8–1.9)**

1.8 (Energy budget) → 1.9 (Stability metrics)

**Phase D — Research infrastructure (Block 1.10–1.21)**

1.10 (Experiment runner) → 1.11 (Geometry provider) → 1.13 (Checkpoint) →
1.14 (Live monitor) → 1.15 (Logging schema) → 1.16 (Scaling) → 1.17 (Vacuum
energy) → 1.18 (Diagnostics) → 1.19 (Seeds) → 1.20 (Sweep DSL) → 1.21
(Artifacts)

**Phase E — Physics variants (Block 2)**

In the order 2.0 → 2.1 → 2.2 → 2.3 → 2.4 → 2.5 → 2.6 → 2.7 → 2.8 → 2.9 →
2.10 → 2.11 → 2.12 → 2.13.

**Phase F — Rendering (Block 3, out of scope here)**

Port the flux mesh, granule motion, and interactive controls from the old
launcher. Only after the physics is validated in headless mode.

---

## 11. Glossary

- **EMC** — Elastic Medium Constituent. The spherical unit of the vacuum lattice.
- **BCC** — Body-Centred Cubic. The lattice geometry of the EMC arrangement.
- **WC** — Wave Centre. A point that reflects incoming waves into outgoing waves.
- **`K`** — Number of wave centres in a soliton. `K = 1` neutrino, `K = 10` electron.
- **`K²λ`** — Maximum standing-wave radius of a soliton with `K` centres.
- **`λ_ν`** — Neutrino wavelength. The fundamental length scale in natural units.
- **`ρ_E`** — Energy density (high inside a soliton).
- **`ρ`** — EMC packing density (low inside a soliton).
- **`N_ν,stat`** — Statutory background EMC density (undisturbed vacuum).
- **`N_ν,eff`** — Effective EMC density inside the soliton.
- **Push-out** — The mechanism by which `ρ_E` displaces EMC, creating `ρ < N_stat`.
- **EMC Wall** — The local peak `ρ > N_stat` at the outer boundary of the soliton's deficit range.
- **`α`** — Fine-structure constant. In EWT, `|Ψ_out|/|Ψ_in|` at a WC.
- **`ε_M`** — Magnetic deficit. `1/(N_geom π³) ≈ 1/(8π⁷(1−ζ))`.
- **`A_π`** — Geometric core of the soliton. `4π³ + π² + π`.
- **`N_geom`** — Effective BCC stiffness. `8π⁴(1−ζ)`.
- **`γ`** — Nonlinear coupling. `1/ε_M`.
- **NESS** — Non-Equilibrium Steady State. The soliton's dynamical regime.
- **Reflector** — A WC that satisfies `|Ψ_out|² + |Ψ_spin|² = |Ψ_in|²`.
- **Feature** — A typed object stored in `FeatureBag`, keyed by its class.
- **Processor** — A stateless pipeline stage that reads and writes features.

---

## 12. How to use this document

This document is a **working roadmap**, not a specification. It records:

- **Why** the tool exists (Sections 1–2).
- **What** the tool must express (Sections 3–5).
- **How** to build it (Sections 6–7).
- **What** to skip (Section 8).
- **What** remains unresolved (Section 9).

Update it as decisions are made. Each work item in Blocks 1 and 2 should be
promoted to a `tasks/m4_<n>_task_details.md` when it is picked up, with
pre-registered pass/fail criteria. The roadmap row in `m4_roadmap.md`
references that task document.

When a work item is complete, mark the checkbox and add a one-line note in
Section 13 (Changelog).

---

## 13. Changelog

| Date | Change | Author |
|---|---|---|
| 2026-09-17 | Initial draft. | Lukasz Smolinski  |

---

## 14. OpenWave Compliance

### 14.1. What this document is, and what it is not

**Is** — a design rationale and a work plan for the pipeline_engine. It
explains *why* the tool exists, *what* it must express, and *how* the work is
organised.

**Is not** — a roadmap, a task document, or a findings note. It does not
replace `m4_roadmap.md`, `tasks/m4_<n>_task_details.md`, or `findings/`.
Those documents carry the criteria, the numbers, and the verdicts.

The correct flow for a work item is:

```text
this document  →  m4_roadmap.md row  →  tasks/m4_<n>_task_details.md
                (design intent)        (preview)               (the record)
                                       ↓
                              scripts/m4_<n>_*.py
                              data/m4_<n>_*.csv
                              plots/m4_<n>_*.png
                              findings/m4_<n>_*.md
```

When this document and a task document disagree, the task document wins. When
this document and the manuscript disagree, the manuscript wins. When the
manuscript and the model author disagree, the author wins.

### 14.2. TaskID mapping

The following TaskIDs are proposed for the roadmap. They are assigned in
creation order and are never reused. The list is a proposal; the M4 maintainer
assigns the final IDs when the rows are added to `m4_roadmap.md`.

**Block 1 — Engine**

| Proposed ID | Item | Depends on |
|---|---|---|
| M4.20 | UnitSystem feature | — |
| M4.21 | Multi-field FeatureBag | M4.20 |
| M4.22 | Trackers per voxel | M4.21 |
| M4.23 | Multi-field evolution | M4.21 |
| M4.24 | Source terms (additive) | M4.23 |
| M4.25 | Reflector interface | M4.23 |
| M4.26 | WC motion (drift) | M4.25 |
| M4.27 | EMC boundary (analytic) | M4.20 |
| M4.28 | Energy budget tracker | M4.22, M4.23 |
| M4.29 | Stability metrics | M4.22 |
| M4.30 | Experiment runner | M4.20 |
| M4.31 | Geometric constants provider | M4.20 |
| M4.32 | Checkpoint / restart | M4.21 |
| M4.33 | Live monitor | M4.22 |
| M4.34 | Research logging schema | M4.30 |
| M4.35 | Universe scaling strategy | M4.20 |
| M4.36 | Vacuum energy strategy | M4.20 |
| M4.37 | Diagnostic hooks | M4.29 |
| M4.38 | Deterministic seeds | M4.30 |
| M4.39 | Parameter sweep DSL | M4.30 |
| M4.40 | Artifact versioning | M4.30 |

**Block 2 — Physics**

| Proposed ID | Item | Depends on |
|---|---|---|
| M4.41 | Soliton assembly contract | M4.23, M4.24 |
| M4.42 | Base wave variants | M4.41 |
| M4.43 | WC reflector variants | M4.25, M4.42 |
| M4.44 | Longitudinal↔transverse coupling | M4.25 |
| M4.45 | EMC density dynamics | M4.21, M4.27 |
| M4.46 | Wave speed modulation | M4.45 |
| M4.47 | Density-modulated nonlinearity | M4.45, M4.46 |
| M4.48 | EMC Wall variants | M4.27, M4.45 |
| M4.49 | WC motion rule variants | M4.26, M4.45 |
| M4.50 | WC topology variants | M4.23 |
| M4.51 | WC spacing variants | M4.50 |
| M4.52 | K-selectivity sweep | M4.47–M4.51 |
| M4.53 | Energy conservation verification | M4.28, M4.52 |
| M4.54 | Emergent α | M4.44, M4.52 |

IDs `M4.1`–`M4.19` are already used or reserved by the existing roadmap
(M4.1 K-selectivity, M4.2 Coulomb, M4.3–M4.12 gravity and emergence work).
The proposed assignment continues the sequence without collision.

*End of document.*
