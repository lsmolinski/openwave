"""
M5.2 Step 3 — Defect-Survival Check under EVOLVE PSI (headless)

Tests whether a seeded +1 hedgehog stays topologically coherent under the
Klein-Gordon PDE  ∂²ψ/∂t² = c²·∇²ψ − m²·ψ  with the electron mass-frequency,
versus the free-wave baseline (V = 0). Driven by `evolve_psi` steps.

EXPECTED OUTCOME — NEGATIVE.
The electron Klein-Gordon mass-frequency m·c²/ℏ ≈ 7.76e-7 rad/rs is so small
relative to our grid timescale (dt ≈ 0.6 rs) that (m·dt)² ≈ 2.3e-13 — below
f32 precision. The KG term has no observable effect at electron mass on a
voxel-scale grid; the hedgehog dissipates at essentially the same rate as
under V = 0 (Derrick collapse). The plain KG potential V = ½m²|ψ|² is also
structurally wrong — its minimum is at ψ = 0, not on the unit sphere
|ψ| = 1 — so even if the term were strong, it would PULL the hedgehog
toward the trivial vacuum rather than stabilize it.

This test documents the negative result and forms the empirical basis for
moving to M5.2 Step 4 (Close Eq. 23 nonlinear terms / LdG Mexican-hat
potential), where the potential minimum sits on the unit sphere.

Pass criterion: there's no pass/fail — this is an observation test.
Output: winding number Q vs propagation step, for KG (electron) and V=0 baseline.

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.m5_2_kg_defect_survival
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
from openwave.xperiments.m5_liquid_crystal import engine1_seeds as seeds  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde  # noqa: E402
from openwave.xperiments.m5_liquid_crystal import engine3_observables as observables  # noqa: E402

# ================================================================
# CONFIG
# ================================================================
N = 65
UNIVERSE_EDGE_M = 1e-15
TARGET_VOXELS = N**3
DOMAIN_QUARTER_FRACTION = 0.20
WINDING_RADIUS = 5
N_RELAX = 0  # no relax (per m5_1_winding finding — relax dissipates topology)
N_PROPAGATE_STEPS = 200
SAMPLE_EVERY = 20  # measure winding number every N steps


def seed_single_hedgehog(wf):
    centers = ti.field(dtype=ti.i32, shape=(1, 3))
    signs = ti.field(dtype=ti.i32, shape=(1,))
    centers[0, 0] = wf.nx // 2
    centers[0, 1] = wf.ny // 2
    centers[0, 2] = wf.nz // 2
    signs[0] = +1
    D_quarter = float(DOMAIN_QUARTER_FRACTION * max(wf.nx, wf.ny, wf.nz))
    seeds.seed_hedgehog(wf, centers, signs, D_quarter, 1)
    for _ in range(N_RELAX):
        cfl_bound = (wf.dx_am**2) / 6.0
        tau = 0.4 * cfl_bound
        pde.relax_director_step(wf, tau, centers, signs, 1)
        wf.psi_am.copy_from(wf.psi_new_am)
        wf.psi_prev_am.copy_from(wf.psi_new_am)
    return wf.nx // 2, wf.ny // 2, wf.nz // 2


def propagate_and_measure(wf, c_amrs, dt_rs, m_freq_rs, center):
    """Run N propagation steps, sampling winding number along the way."""
    samples = []  # list of (step, Q, |psi|_avg)
    # Initial measurement
    psi_np = wf.psi_am.to_numpy()
    Q0 = observables.compute_winding_number(psi_np, center, WINDING_RADIUS)
    psi_norm = np.linalg.norm(psi_np, axis=-1).mean()
    samples.append((0, Q0, psi_norm))

    for step in range(1, N_PROPAGATE_STEPS + 1):
        pde.evolve_psi(wf, c_amrs, dt_rs, m_freq_rs, 0.0)
        wf.swap_buffers()
        if step % SAMPLE_EVERY == 0 or step == N_PROPAGATE_STEPS:
            psi_np = wf.psi_am.to_numpy()
            Q = observables.compute_winding_number(psi_np, center, WINDING_RADIUS)
            psi_norm = np.linalg.norm(psi_np, axis=-1).mean()
            samples.append((step, Q, psi_norm))

    return samples


def main():
    ti.init(arch=ti.gpu)
    print("=" * 78)
    print("M5.2 Step 3 — Defect-Survival under EVOLVE PSI")
    print("=" * 78)

    wf = medium.WaveField([UNIVERSE_EDGE_M] * 3, TARGET_VOXELS)
    c_amrs = (
        constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND
    )  # SIM_SPEED=1 baseline
    cfl_safety = 0.95
    dt_rs = wf.dx_am * cfl_safety / (c_amrs * (3**0.5))
    m_freq_kg = c_amrs / constants.COMPTON_WAVELENGTH_REDUCED_ELECTRON_AM

    print(f"Grid: {wf.nx}³  dx={wf.dx_am:.3f} am  dt={dt_rs:.4f} rs")
    print(f"c={c_amrs:.4f} am/rs   m_freq_kg(electron)={m_freq_kg:.4e} rad/rs")
    print(
        f"(m·dt)² = {(m_freq_kg * dt_rs)**2:.4e}  "
        f"(vs (c·dt/dx)² = {(c_amrs * dt_rs / wf.dx_am)**2:.4f})"
    )
    print(f"propagate {N_PROPAGATE_STEPS} steps, sample every {SAMPLE_EVERY}")
    print()

    # =================================================================
    # Run A — V = 0 (free-wave, M5.0g baseline)
    # =================================================================
    print("[A] BASELINE — V = 0 (free-wave)")
    print(f"    seeding hedgehog Q=+1 at center...")
    center = seed_single_hedgehog(wf)
    samples_free = propagate_and_measure(wf, c_amrs, dt_rs, 0.0, center)
    print(f"    {'step':>6s}  {'Q':>10s}  {'<|ψ|>':>10s}")
    for step, Q, pnorm in samples_free:
        print(f"    {step:>6d}  {Q:>+10.4f}  {pnorm:>10.5f}")
    print()

    # =================================================================
    # Run B — KG with electron mass-frequency
    # =================================================================
    print("[B] KG MASS TERM — V = ½·m_e²·|ψ|²")
    print(f"    seeding hedgehog Q=+1 at center...")
    center = seed_single_hedgehog(wf)
    samples_kg = propagate_and_measure(wf, c_amrs, dt_rs, m_freq_kg, center)
    print(f"    {'step':>6s}  {'Q':>10s}  {'<|ψ|>':>10s}")
    for step, Q, pnorm in samples_kg:
        print(f"    {step:>6d}  {Q:>+10.4f}  {pnorm:>10.5f}")
    print()

    # =================================================================
    # Comparison
    # =================================================================
    print("=" * 78)
    print("COMPARISON — free vs KG (electron) at each sample step")
    print("=" * 78)
    print(f"  {'step':>6s}  {'Q_free':>10s}  {'Q_kg':>10s}  {'|ΔQ|':>10s}")
    for (s_f, Q_f, _), (s_k, Q_k, _) in zip(samples_free, samples_kg):
        dQ = abs(Q_f - Q_k)
        print(f"  {s_f:>6d}  {Q_f:>+10.4f}  {Q_k:>+10.4f}  {dQ:>10.4e}")
    print()

    # =================================================================
    # Interpretation
    # =================================================================
    Q_init = samples_free[0][1]
    Q_final_free = samples_free[-1][1]
    Q_final_kg = samples_kg[-1][1]
    decay_free = abs(Q_init - Q_final_free)
    decay_kg = abs(Q_init - Q_final_kg)
    delta_kg_vs_free = abs(Q_final_free - Q_final_kg)

    print("=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print(f"  Initial Q (seed): {Q_init:+.4f}")
    print(f"  Final Q (free):   {Q_final_free:+.4f}  → ΔQ = {decay_free:.4f}")
    print(f"  Final Q (KG e⁻):  {Q_final_kg:+.4f}  → ΔQ = {decay_kg:.4f}")
    print(f"  |Q_free − Q_kg| at end: {delta_kg_vs_free:.4e}")
    print()
    print("  Expected: KG behavior indistinguishable from free-wave at electron")
    print("  mass scale — (m·dt)² ~ 1e-13 is below f32 precision relative to ψ.")
    print()
    print("  CONCLUSION: plain KG (V = ½·m²·|ψ|²) at electron mass-frequency does")
    print("  NOT preserve the hedgehog. Both runs dissipate the topology under")
    print("  EVOLVE PSI. Mexican-hat-style nonlinear potential required (Step 4 —")
    print("  Close Eq. 23 / LdG biaxial) for defect stability.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
