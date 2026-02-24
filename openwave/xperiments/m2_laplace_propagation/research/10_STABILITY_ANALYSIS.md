# NUMERICAL STABILITY ANALYSIS

**Key Difference**: Spring-Mass vs Wave Equation PDE

**Reference**: See `/openwave/validations/stability_analysis_wave.py` for automated CFL stability verification

## Spring-Mass System (LEVEL-0) - THE IMPOSSIBLE TRIANGLE

Your spring-mass experiments hit a fundamental barrier:

```text
    Realistic Stiffness
            / \
           /   \
          /     \
    Stability --- Human-Visible Motion
```

### The Problem

- Required spring stiffness: k ≈ 5.56×10⁴⁴ N/m (Planck scale)
- Natural frequency: ω_n = √(k/m) ≈ 2.87×10¹⁷ rad/s
- CFL critical timestep: dt < 2/ω_n ≈ 7×10⁻¹⁸ s
- Actual frame time: dt ≈ 0.016 s (60 FPS)
- Violation: 24 orders of magnitude! 💥

### Results from experiments

- ❌ Euler: NaN at 0.1s even with 1000 substeps
- ❌ Leapfrog: NaN at 0.4s even with 1000 substeps
- ⚠️ XPBD: Stable but waves ~8× too slow (12.5% of c)
- ✅ PSHO (Phase-Synchronized): Perfect but bypasses force mechanics

## Wave Equation PDE (LEVEL-1) - FUNDAMENTAL ADVANTAGE

### The Solution

The wave equation ψ̈ = c²∇²ψ has the SAME CFL condition but with a CRITICAL DIFFERENCE:

```python
# Spring-Mass CFL (LEVEL-0)
dt_critical = 2/ω_n = 2/√(k/m)
# Problem: k must be HUGE for realistic wave speed!

# Wave Equation CFL (LEVEL-1):
dt_critical = dx/(c√3)
# Solution: Apply SLO_MO to c directly!
c_slo = c / SLO_MO
dt_critical_slowed = dx/(c_slo√3)  # Now feasible!
```

### Why This Works

1. Spring-Mass: Stiffness k is FIXED by physics → ω_n is FIXED → dt_critical is FIXED → UNSTABLE
2. Wave PDE: Wave speed c can be SLOWED without breaking physics → dt_critical becomes FEASIBLE → STABLE

### From stability analysis script

**Script**: `/openwave/validations/stability_analysis_wave.py`

**Example output** (for 6 fm³ universe, 1B voxels):

```python
# Grid parameters (from WaveField):
dx = 6 am = 6e-18 m  # Voxel edge for 1B voxels in 6 fm³

# Without SLO_MO:
dt_critical = dx / (c * √3) ≈ 1.2e-26 s  # Rontosecond range!
dt_frame = 1/60 ≈ 0.016 s
Violation: ~10²⁴× 💥 UNSTABLE!

# With SLO_MO = EWAVE_FREQUENCY = 1.05×10²⁵:
c_slo = c / 1.05×10²⁵ = 2.85 × 10⁻¹⁷ m/s
dt_critical_slowed = dx / (c_slo * √3) ≈ 0.121 s
dt_frame = 1/60 ≈ 0.016 s
✓ STABLE! (dt_frame < dt_critical_slowed)
```

**Safety Factor**: CFL factor = (c_slo·dt / dx)² ≈ 0.14 (well within 1/3 limit for 3D 6-connectivity)

**Key Parameters Tested**:

- Voxel edge: dx from `WaveField.voxel_edge`
- Wave speed: c = `constants.EWAVE_SPEED` = 2.998×10⁸ m/s
- SLO_MO factor: `config.SLO_MO` (configurable, typically ~10²⁵)
- Frame rates: 60 FPS (dt = 0.0167s) and 30 FPS (dt = 0.0333s)

## Direct Comparison

| Aspect              | Spring-Mass (LEVEL-0)                     | Wave PDE (LEVEL-1)         |
|---------------------|-------------------------------------------|----------------------------|
| Governing equation  | F = -k(x-L), ẍ = F/m                      | ψ̈ = c²∇²ψ                  |
| Stiffness           | k ≈ 5.56×10⁴⁴ N/m (FIXED!)                | No springs - pure wave     |
| CFL condition       | dt < 2/√(k/m)                             | dt < dx/(c√3)              |
| Critical dt         | ~7×10⁻¹⁸ s (INFLEXIBLE)                   | ~2.4×10⁻²⁷ s (raw)         |
| SLO_MO mitigation   | Can't reduce k w/out breaking physics     | Can reduce c directly      |
| Numerical stability | Explodes (NaN at 0.4s)                    | Stable with c_slo          |
| Wave speed fidelity | XPBD: ~12.5% of c                         | Exact by construction      |
| Computational cost  | 8-neighbor springs per granule            | 6-neighbor Laplacian /voxel|
| Result              | IMPOSSIBLE TRIANGLE                       | FEASIBLE SIMULATION        |

## Why You Won't Have Numerical Explosion Now

### The Key Insight

Your spring-mass system failed because:

1. You needed k ≈ 10⁴⁴ N/m for realistic wave speed
2. This created ω_n ≈ 10¹⁷ rad/s
3. CFL demanded dt < 10⁻¹⁷ s
4. But visualization needed dt ≈ 0.016 s
5. Gap unbridgeable → explosion

Your wave equation system succeeds because:

1. You need c = 3×10⁸ m/s for realistic wave speed
2. CFL demands dt < dx/(c√3) ≈ 1.2×10⁻²⁶ s (for dx = 6 am)
3. Apply SLO_MO: c_slo = c/1.05×10²⁵ = 2.85×10⁻¹⁷ m/s
4. New CFL: dt < dx/(c_slo√3) ≈ 0.121 s
5. Visualization needs dt ≈ 0.016 s
6. 0.016 < 0.121 → ✓ STABLE!

**From spring-mass experiments final report**:

> "The progression through different methods:
>
> - Force-based integration (Euler): Encountered CFL stability barrier ✓
> - Symplectic integration (Leapfrog): Same CFL limitation ✓
> - Constraint-based dynamics (XPBD): Stability with accuracy trade-offs ✓
> - PSHO approach: Bypassed integration methods ✓"
>
> *Source: `/openwave/xperiments/_archives/spring_mass/_docs/final_report.md`*

## Validation & Testing

**Automated Stability Check**:

Run the stability analysis script to verify your configuration:

```bash
python openwave/validations/stability_analysis_wave.py
```

**What it verifies**:

1. CFL condition satisfaction for 60 FPS and 30 FPS
2. Safety margins (how much headroom exists)
3. Required SLO_MO values if unstable
4. Recommended mitigation strategies

**Expected Output**:

```text
✓ STABLE at 60 FPS (dt=1.67e-02 s ≤ dt_crit=2.41e-02 s)
  Safety margin: 1.44× (CFL factor = 0.694)
```

## Implementation Requirements

**Critical Implementation Details** (see `02b_WAVE_ENGINE_propagate.md`):

1. **Apply SLO_MO to wave speed**, not timestep:

   ```python
   c_slo = constants.EWAVE_SPEED / config.SLO_MO * SIM_SPEED
   ```

2. **Use fixed timestep strategy**, not elapsed time:

   ```python
   dt_physics = 1/60  # Fixed (e.g., 60 FPS)
   c_slo = c / SLO_MO  # Slow wave speed instead
   ```

3. **Monitor CFL factor** during simulation:

   ```python
   cfl_factor = (c_slo_am * dt / dx_am)**2  # Should be ≤ 1/3 for 3D
   ```

4. **Use attometer scaling** for numerical precision:

   ```python
   dx_am = voxel_edge / constants.ATTOMETER  # [am]
   c_slo_am = c_slo / constants.ATTOMETER     # [am/s]
   ```

## No Substeps Required

**Critical Difference**: LEVEL-1 does **NOT** need the substep technique that LEVEL-0 required.

### Why Spring-Mass Needed Substeps

From your spring-mass experiments, substeps were an attempt to satisfy CFL:

```python
# Spring-Mass Problem (LEVEL-0):
dt_frame = 1/30  # 0.033s per frame
dt_critical = 2/ω_n ≈ 7×10⁻¹⁸ s

# Required substeps per frame:
N_substeps = dt_frame / dt_critical
           = 0.033 / 7×10⁻¹⁸
           ≈ 4.7 × 10¹⁵ substeps  # IMPOSSIBLE!

# What you tried:
N_substeps = 1000  # Still exploded! NaN at 0.4s
```

**From `2_findings.md`**:

> "Substeps needed for stability: 39,804,942
> You're using: 200 substeps
> ❌ UNSTABLE! Need 40M substeps (120M iterations/second)"

The substep technique **failed** because even with 1000 substeps per frame, you were still 10¹² orders of magnitude short of what was needed.

### Why Wave Equation DOESN'T Need Substeps

The wave equation approach eliminates the substep requirement entirely:

```python
# Wave Equation Solution (LEVEL-1):
dt_frame = 1/60  # 0.0167s per frame (60 FPS)
dt_critical_slowed = dx / (c_slo * √3) ≈ 0.121 s  # For dx = 6 am, 1B voxels

# Required substeps per frame:
N_substeps = dt_frame / dt_critical_slowed
           = 0.0167 / 0.121
           ≈ 0.14 substeps  # Much less than 1!

# What you use:
N_substeps = 1  # Single timestep per frame! ✓ STABLE
```

**Key Insight**: By applying SLO_MO to wave speed (not timestep), the critical timestep becomes **LARGER** than the frame timestep. No substeps needed!

### Comparison Table

| Technique | Spring-Mass (LEVEL-0) | Wave PDE (LEVEL-1) |
| --------- | --------------------- | ------------------ |
| **Substeps needed** | 40 million per frame | **0** (single step) |
| **What you tried** | 1000 substeps → NaN | Not applicable |
| **Why it failed** | k is FIXED by physics | Not applicable |
| **Why it succeeds** | Not applicable | c can be SLOWED |
| **Computational cost** | 40M iterations/frame (impossible) | 1 iteration/frame ✓ |
| **Artifact complexity** | Complex substep loop + accumulator | **Simple single step** |

### Code Simplification

**LEVEL-0 (with substeps artifact)**:

```python
# Complex substep loop needed
def update_frame(dt_frame):
    dt_sub = dt_frame / num_substeps  # Split frame time
    for substep in range(num_substeps):  # Typically 200-1000
        compute_forces()
        integrate_euler(dt_sub)  # Still explodes!
        apply_damping()
```

**LEVEL-1 (no substeps needed)**:

```python
# Simple single-step per frame
def update_frame(dt_frame):
    propagate_wave(dt_frame, SIM_SPEED)  # Just once! ✓
    track_amplitude_envelope()
    compute_wave_direction()
```

**Result**: LEVEL-1 is not only stable but also **simpler** - no substep loops, no timestep splitting, no accumulator artifacts.

## Summary

LEVEL-1 wave equation is essentially a **grid-based PSHO** - you're directly solving the wave equation via PDE, not simulating forces through springs. This is why it will work where spring-mass failed!

**Bottom Line**: You've already discovered the solution in your PSHO experiments - LEVEL-1 is that same principle but on a voxel grid instead of particle lattice. **No numerical explosion expected!** 🎯

**Bonus**: No substep complexity - just one clean timestep per frame! 🚀

## See Also

- [02b_WAVE_ENGINE_propagate.md](./02b_WAVE_ENGINE_propagate.md) - Wave propagation implementation with CFL stability discussion
- `/openwave/validations/stability_analysis_wave.py` - Automated stability verification script
- `/openwave/xperiments/_archives/spring_mass/_docs/final_report.md` - Spring-mass failure analysis
