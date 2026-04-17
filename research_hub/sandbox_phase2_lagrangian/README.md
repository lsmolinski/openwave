# Phase 3 — Lagrangian Framework Numerical Experiments

Standalone numpy research scripts for evaluating Lagrangian-based wave equations.

> **📍 Status** (2026-04-17): all 8 sandbox experiments are complete. See [../3b_lagrangian_experiments.md](../3b_lagrangian_experiments.md) for results. The recommended-order workflow below is kept as historical record; no further sandbox work is planned (optional deferred: Exp 6.1 full biaxial Q-tensor). Next phase: M5 production implementation per [../3c_path_to_m5.md](../3c_path_to_m5.md).
>
> Note: directory name still `sandbox_phase2_lagrangian/` for backward compatibility (the content originated as a Phase-2 sub-project before being promoted to Phase 3 on 2026-04-17).

## Why "Sandbox"

OpenWave uses two layers of experimentation:

- **Sandbox** (`research_hub/sandbox_*/`) — quick numpy scripts to validate math, logic, and concepts. Pure exploration. No GPU, no Taichi, no production dependencies. This is where ideas are tested cheaply before committing engineering effort. If an experiment works in the sandbox, it graduates to the production engine.
- **Production** (`openwave/xperiments/m*/`) — Taichi-based 3D rendering and simulation on the official OpenWave platform. GPU-accelerated, full grid infrastructure, visualization, force & motion. This is the final product where validated equations run at scale.

**Sandbox = explore fast, fail cheap. Production = implement validated winners.**

## Purpose

Validate (or rule out) candidate Lagrangians before committing to any architecture change in M3/M4. These scripts test the core ideas quickly without requiring GPU acceleration or production engine refactors.

## Spec & Results

- **Spec**: [../3_LAGRANGIAN_FRAMEWORK.md](../3_LAGRANGIAN_FRAMEWORK.md) — full experiment specifications, hypotheses, success criteria
- **Results**: [../3b_lagrangian_experiments.md](../3b_lagrangian_experiments.md) — running log of results, comparisons to expected, conclusions

## File Naming Convention

```text
exp{N}_{short_name}.py
```

Examples:

- `exp1_sine_gordon_1d.py`
- `exp2_hedgehog_energy.py`
- `exp3_topological_charge.py`
- `exp4_klein_gordon.py`
- `exp5_lagrangian_derivation.py`
- `exp6_lepton_families.py`
- `exp7_close_vector_wave.py`
- `exp8_smolinski_psi3.py`

## Pattern (same as Phase 1 scripts)

Each script is self-contained:

- Pure numpy (no Taichi, no GPU)
- Hardcoded grid sizes (32³ to 128³ typically)
- Analytical or time-domain PDE evolution
- Output: numbers, plots (matplotlib), or both
- No production dependencies

Reference: `../sandbox_phase1_vector/` and `../sandbox_phase1_scalar/` for the established pattern.

## Recommended Order

1. **Experiment 1** — Sine-Gordon 1D (build intuition, simplest soliton)
2. **Experiments 2 + 3** — Hedgehog energy + topological charge (highest-value: Coulomb from topology)
3. **Experiment 5** — Lagrangian derivation (math, informs everything)
4. **Experiment 8** — Smolinski Ψ³ (K-selectivity on familiar scalar setup)
5. **Experiment 4** — Klein-Gordon dynamics
6. **Experiment 6** — Three lepton families (most complex)
7. **Experiment 7** — Close's nonlinear vector wave equation

## Workflow

For each experiment:

1. Fill in `Setup` section in `3b_lagrangian_experiments.md` with parameters/grid before running
2. Run script, capture output
3. Fill in `Results` and `Numerical Evidence` sections
4. Update `Comparison to Expected` table
5. Write `Conclusion` and `Next Steps`
6. Update `Status` in the Summary Dashboard at the top of `2b`
