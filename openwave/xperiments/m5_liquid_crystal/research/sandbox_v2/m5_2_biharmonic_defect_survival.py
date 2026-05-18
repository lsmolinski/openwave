"""
M5.2 Step 4b — Biharmonic Defect-Survival under EVOLVE PSI (headless)

Tests whether the 4th-derivative biharmonic term `+ ½κ|∇²ψ|²` preserves the
+1 hedgehog winding number under wave propagation, where φ⁴ Mexican-hat alone
(Step 4a) failed to. Three-way comparison:

  [A] V = 0                                      — free-wave (reference Derrick)
  [B] V = ½m²|ψ|² + ¼λ(|ψ|²−1)²                  — Step 4a Mexican-hat
  [C] V = ½m²|ψ|² + ¼λ(|ψ|²−1)² + ½κ|∇²ψ|²       — Step 4b adds biharmonic

The biharmonic term scales as E ~ 1/R under spatial dilation while gradient
energy scales as E ~ R; together they give a Derrick-balanced minimum at
finite R. Combined with the φ⁴ Mexican-hat (which holds |ψ| near 1), this is
the minimum-viable analog of the classical Skyrme model for 3D hedgehogs.

EXPECTED OUTCOME — IF SUCCESSFUL:
The biharmonic run [C] preserves Q over more propagation steps than the
free-wave [A] or Mexican-hat [B] runs. If Q decays at a similar rate to
[A]/[B], we'll know the specific cross-product structure of Faddeev-Skyrme
is required (Step 4b.2).

Pass criterion (observation):
  - Step at which |Q| drops below 0.5 — biharmonic should EXTEND this step
    significantly relative to free-wave and φ⁴ alone.

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.m5_2_biharmonic_defect_survival
"""

import sys
from pathlib import Path

import numpy as np
import taichi as ti

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from openwave.common import constants  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import medium  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import lagrangian_engine as lagrange  # noqa: E402

# ================================================================
# BIHARMONIC KERNELS — sandbox-only, not in production lagrangian_engine
# ================================================================
# The biharmonic stabilizer (M5.2 Step 4b experiment) produced a negative
# result, so the kernels live HERE rather than in production. If a future
# experiment validates a 4th-derivative term as part of the working V(ψ)
# recipe, promote these into lagrangian_engine.py at that point.
#
# Two-pass design: compute Laplacian into a scratch field, then compute the
# Laplacian-of-Laplacian inline. Reuses `lagrange.compute_laplacian` (the
# @ti.func used by evolve_psi) for the first pass.


@ti.kernel
def compute_psi_laplacian(
    wave_field: ti.template(),  # type: ignore
    lap_field: ti.template(),  # type: ignore  -- Vector(3) field with same shape as ψ
):
    """Fill lap_field[i,j,k] with ∇²ψ from wave_field.psi_am."""
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    for i, j, k in ti.ndrange((1, nx - 1), (1, ny - 1), (1, nz - 1)):
        lap_field[i, j, k] = lagrange.compute_laplacian(wave_field, i, j, k)


@ti.kernel
def apply_biharmonic_correction(
    wave_field: ti.template(),  # type: ignore
    lap_field: ti.template(),  # type: ignore  -- already-filled ∇²ψ
    dt_rs: ti.f32,  # type: ignore
    kappa_bi: ti.f32,  # type: ignore
):
    """Subtract (dt²·κ·∇⁴ψ) from psi_new_am. Requires lap_field already filled."""
    nx, ny, nz = wave_field.nx, wave_field.ny, wave_field.nz
    inv_dx2 = 1.0 / (wave_field.dx_am * wave_field.dx_am)
    kappa_dt2 = kappa_bi * dt_rs * dt_rs
    for i, j, k in ti.ndrange((2, nx - 2), (2, ny - 2), (2, nz - 2)):
        center = lap_field[i, j, k]
        face_sum = (
            lap_field[i + 1, j, k]
            + lap_field[i - 1, j, k]
            + lap_field[i, j + 1, k]
            + lap_field[i, j - 1, k]
            + lap_field[i, j, k + 1]
            + lap_field[i, j, k - 1]
        )
        bilap = (face_sum - 6.0 * center) * inv_dx2
        wave_field.psi_new_am[i, j, k] -= kappa_dt2 * bilap


# ================================================================
# CONFIG
# ================================================================
N = 65
UNIVERSE_EDGE_M = 1e-15
TARGET_VOXELS = N**3
DOMAIN_QUARTER_FRACTION = 0.20
WINDING_RADIUS = 5
N_RELAX = 20  # brief relax to smooth seed before propagation
N_PROPAGATE_STEPS = 400
SAMPLE_EVERY = 20


def seed_and_relax(wf, n_relax):
    centers = ti.field(dtype=ti.i32, shape=(1, 3))
    signs = ti.field(dtype=ti.i32, shape=(1,))
    centers[0, 0] = wf.nx // 2
    centers[0, 1] = wf.ny // 2
    centers[0, 2] = wf.nz // 2
    signs[0] = +1
    D_quarter = float(DOMAIN_QUARTER_FRACTION * max(wf.nx, wf.ny, wf.nz))
    lagrange.seed_hedgehog(wf, centers, signs, D_quarter, 1)
    cfl_bound = (wf.dx_am**2) / 6.0
    tau = 0.4 * cfl_bound
    for _ in range(n_relax):
        lagrange.relax_director_step(wf, tau, centers, signs, 1)
        wf.psi_am.copy_from(wf.psi_new_am)
        wf.psi_prev_am.copy_from(wf.psi_new_am)
    return wf.nx // 2, wf.ny // 2, wf.nz // 2


def propagate(wf, lap_field, c_amrs, dt_rs, m_freq_rs, lambda_phi4, kappa_bi, center):
    samples = []
    psi_np = wf.psi_am.to_numpy()
    Q0 = lagrange.compute_winding_number(psi_np, center, WINDING_RADIUS)
    n0 = np.linalg.norm(psi_np, axis=-1)
    samples.append((0, Q0, n0.mean(), n0.max()))

    for step in range(1, N_PROPAGATE_STEPS + 1):
        lagrange.evolve_psi(wf, c_amrs, dt_rs, m_freq_rs, lambda_phi4)
        if kappa_bi > 0:
            compute_psi_laplacian(wf, lap_field)
            apply_biharmonic_correction(wf, lap_field, dt_rs, kappa_bi)
        wf.swap_buffers()
        if step % SAMPLE_EVERY == 0 or step == N_PROPAGATE_STEPS or step <= 5:
            psi_np = wf.psi_am.to_numpy()
            if np.isnan(psi_np).any():
                samples.append((step, float("nan"), float("nan"), float("nan")))
                break
            Q = lagrange.compute_winding_number(psi_np, center, WINDING_RADIUS)
            n = np.linalg.norm(psi_np, axis=-1)
            samples.append((step, Q, n.mean(), n.max()))

    return samples


def print_table(label, samples):
    print(f"  {label}")
    print(f"    {'step':>6s}  {'Q':>10s}  {'<|ψ|>':>10s}  {'max|ψ|':>10s}")
    for step, Q, nmean, nmax in samples:
        Q_str = f"{Q:+10.4f}" if not np.isnan(Q) else "       NaN"
        nm_str = f"{nmean:10.5f}" if not np.isnan(nmean) else "       NaN"
        nx_str = f"{nmax:10.5f}" if not np.isnan(nmax) else "       NaN"
        print(f"    {step:>6d}  {Q_str}  {nm_str}  {nx_str}")


def first_decay_step(samples, threshold=0.5):
    for step, Q, _, _ in samples:
        if np.isnan(Q) or abs(Q) < threshold:
            return step
    return None


def main():
    ti.init(arch=ti.gpu)
    print("=" * 78)
    print("M5.2 Step 4b — Biharmonic Defect-Survival under EVOLVE PSI")
    print("=" * 78)

    wf = medium.WaveField([UNIVERSE_EDGE_M] * 3, TARGET_VOXELS)
    # Local scratch field for biharmonic (∇²ψ intermediate). Defined here in
    # the research script — not on FieldObservables — so production stays clean.
    lap_field = ti.Vector.field(3, dtype=ti.f32, shape=wf.grid_size)
    c_amrs = constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND
    cfl_safety = 0.95
    dt_rs = wf.dx_am * cfl_safety / (c_amrs * (3**0.5))
    m_freq_kg = c_amrs / constants.COMPTON_WAVELENGTH_REDUCED_ELECTRON_AM
    lambda_phi4 = 0.1 * (c_amrs / wf.dx_am) ** 2
    # Empirically: 0.01 blows up at step 100 due to nonlinear φ⁴ feedback;
    # 0.003 stable 200+ steps; 0.001 has 10× safety margin.
    kappa_bi = 0.001 * (c_amrs**2) * (wf.dx_am**2)

    print(f"Grid: {wf.nx}³  dx={wf.dx_am:.3f} am  dt={dt_rs:.4f} rs")
    print(f"c={c_amrs:.4f} am/rs   m_freq_kg(electron)={m_freq_kg:.4e} rad/rs")
    print(f"λ_φ⁴ = {lambda_phi4:.4e}   κ_bi = {kappa_bi:.4e}")
    print(f"λ·dt²  = {lambda_phi4 * dt_rs**2:.4f}")
    print(f"(c·dt/dx)² = {(c_amrs * dt_rs / wf.dx_am)**2:.4f}")
    print(f"κ·(π/dx)⁴·dt² = {kappa_bi * (np.pi/wf.dx_am)**4 * dt_rs**2:.4f}  (stability check)")
    print(f"propagate {N_PROPAGATE_STEPS} steps, relax {N_RELAX} steps first")
    print()

    # [A] V = 0
    print("[A] BASELINE — V = 0 (free-wave)")
    center = seed_and_relax(wf, N_RELAX)
    samples_A = propagate(wf, lap_field, c_amrs, dt_rs, 0.0, 0.0, 0.0, center)
    print_table("free Q(t):", samples_A)
    print()

    # [B] φ⁴ only
    print("[B] φ⁴ MEXICAN-HAT — V = ½m²|ψ|² + ¼λ(|ψ|² − 1)²")
    center = seed_and_relax(wf, N_RELAX)
    samples_B = propagate(wf, lap_field, c_amrs, dt_rs, m_freq_kg, lambda_phi4, 0.0, center)
    print_table("φ⁴ Q(t):", samples_B)
    print()

    # [C] φ⁴ + biharmonic
    print("[C] φ⁴ + BIHARMONIC — V = ½m²|ψ|² + ¼λ(|ψ|² − 1)² + ½κ|∇²ψ|²")
    center = seed_and_relax(wf, N_RELAX)
    samples_C = propagate(wf, lap_field, c_amrs, dt_rs, m_freq_kg, lambda_phi4, kappa_bi, center)
    print_table("φ⁴ + biharmonic Q(t):", samples_C)
    print()

    # Summary
    decay_A = first_decay_step(samples_A)
    decay_B = first_decay_step(samples_B)
    decay_C = first_decay_step(samples_C)
    print("=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print(f"  First step where |Q| < 0.5:")
    print(f"    [A] free-wave:        {decay_A}")
    print(f"    [B] φ⁴:               {decay_B}")
    print(f"    [C] φ⁴ + biharmonic:  {decay_C}")
    print()
    if decay_C is None:
        print("  ✅ BIHARMONIC PRESERVED Q for the entire run — defect-survival WIN")
    elif decay_A is not None and decay_C is not None and decay_C > 2 * decay_A:
        print(
            f"  ✅ BIHARMONIC EXTENDS lifetime by {decay_C / max(decay_A,1):.1f}× — significant improvement"
        )
    elif decay_A is not None and decay_C is not None and decay_C > decay_A:
        print(
            f"  ⚠️ Modest extension: {decay_C / max(decay_A,1):.2f}× — biharmonic helps but may need more"
        )
    else:
        print(
            "  ⚠️ Biharmonic did NOT meaningfully preserve Q — Step 4b.2 needed (full Faddeev-Skyrme)"
        )
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
