"""
M5.2 Step 4a — Mexican-hat φ⁴ Defect-Survival under EVOLVE PSI (headless)

Tests whether a seeded +1 hedgehog survives propagation when V(ψ) carries the
Mexican-hat term:

    V(ψ) = ½·m²·|ψ|² + ¼·λ·(|ψ|² − 1)²

The φ⁴ term moves the potential minimum from ψ = 0 (plain KG) to the unit
sphere |ψ| = 1 — the right shape for stabilizing a director field. Compared
side-by-side with the V = 0 baseline.

EXPECTED OUTCOME — PARTIAL SUCCESS.
The φ⁴ term should prevent |ψ| from "melting" away from 1 in vacuum (it locks
the field magnitude on the unit sphere). But it does NOT prevent the hedgehog
core from collapsing to a point: in 3D, kinetic energy ~ R and potential
energy ~ R³ both decrease monotonically as R → 0 (Derrick's theorem). To fully
stabilize against collapse, M5.2 Step 4b will add a 4th-derivative Skyrme
term — the classical Skyrme model is historically known to support stable 3D
hedgehog solitons.

Outcomes recorded per run:
  - Q(t) — winding number, samples every SAMPLE_EVERY steps
  - <|ψ|>(t) — mean field magnitude (φ⁴ should hold this near 1.0)
  - V(t) — mean potential energy density (gauges how much "stress" φ⁴ carries)

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.m5_2_phi4_defect_survival
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
# CONFIG
# ================================================================
N = 65
UNIVERSE_EDGE_M = 1e-15
TARGET_VOXELS = N**3
DOMAIN_QUARTER_FRACTION = 0.20
WINDING_RADIUS = 5
N_RELAX = 0
N_PROPAGATE_STEPS = 400  # longer than Step 3 — φ⁴ may slow the collapse
SAMPLE_EVERY = 40


def seed_single_hedgehog(wf):
    centers = ti.field(dtype=ti.i32, shape=(1, 3))
    signs = ti.field(dtype=ti.i32, shape=(1,))
    centers[0, 0] = wf.nx // 2
    centers[0, 1] = wf.ny // 2
    centers[0, 2] = wf.nz // 2
    signs[0] = +1
    D_quarter = float(DOMAIN_QUARTER_FRACTION * max(wf.nx, wf.ny, wf.nz))
    lagrange.seed_hedgehog(wf, centers, signs, D_quarter, 1)
    for _ in range(N_RELAX):
        cfl_bound = (wf.dx_am**2) / 6.0
        tau = 0.4 * cfl_bound
        lagrange.relax_director_step(wf, tau, centers, signs, 1)
        wf.psi_am.copy_from(wf.psi_new_am)
        wf.psi_prev_am.copy_from(wf.psi_new_am)
    return wf.nx // 2, wf.ny // 2, wf.nz // 2


def propagate_and_measure(wf, obs, c_amrs, dt_rs, m_freq_rs, lambda_phi4, center):
    samples = []  # list of (step, Q, |psi|_avg, V_avg)
    psi_np = wf.psi_am.to_numpy()
    Q0 = lagrange.compute_winding_number(psi_np, center, WINDING_RADIUS)
    psi_norm = np.linalg.norm(psi_np, axis=-1).mean()
    # Initial V contribution: just compute energyH and subtract kinetic+gradient
    # Simpler: compute (|ψ|²−1)² mean as a proxy for φ⁴ stress
    psi_sq = (psi_np ** 2).sum(axis=-1)
    v_phi4_proxy = (0.25 * lambda_phi4 * (psi_sq - 1.0) ** 2).mean()
    samples.append((0, Q0, psi_norm, v_phi4_proxy))

    for step in range(1, N_PROPAGATE_STEPS + 1):
        lagrange.evolve_psi(wf, c_amrs, dt_rs, m_freq_rs, lambda_phi4)
        wf.swap_buffers()
        if step % SAMPLE_EVERY == 0 or step == N_PROPAGATE_STEPS:
            psi_np = wf.psi_am.to_numpy()
            Q = lagrange.compute_winding_number(psi_np, center, WINDING_RADIUS)
            psi_norm = np.linalg.norm(psi_np, axis=-1).mean()
            psi_sq = (psi_np ** 2).sum(axis=-1)
            v_phi4 = (0.25 * lambda_phi4 * (psi_sq - 1.0) ** 2).mean()
            samples.append((step, Q, psi_norm, v_phi4))

    return samples


def print_table(label, samples):
    print(f"  {label}")
    print(f"    {'step':>6s}  {'Q':>10s}  {'<|ψ|>':>10s}  {'<V_φ⁴>':>12s}")
    for step, Q, pnorm, v in samples:
        print(f"    {step:>6d}  {Q:>+10.4f}  {pnorm:>10.5f}  {v:>12.5e}")


def main():
    ti.init(arch=ti.gpu)
    print("=" * 78)
    print("M5.2 Step 4a — Mexican-hat φ⁴ Defect-Survival under EVOLVE PSI")
    print("=" * 78)

    wf = medium.WaveField([UNIVERSE_EDGE_M] * 3, TARGET_VOXELS)
    obs = medium.FieldObservables(wf.grid_size)
    c_amrs = constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND
    cfl_safety = 0.95
    dt_rs = wf.dx_am * cfl_safety / (c_amrs * (3**0.5))
    m_freq_kg = c_amrs / constants.COMPTON_WAVELENGTH_REDUCED_ELECTRON_AM
    # Empirical: 1.0·(c/dx)² blows up around step 40 due to nonlinear |ψ|³
    # amplification; 0.3 is marginal; 0.1 is comfortable. Use 0.1 to match the
    # launcher default.
    lambda_phi4 = 0.1 * (c_amrs / wf.dx_am) ** 2

    print(f"Grid: {wf.nx}³  dx={wf.dx_am:.3f} am  dt={dt_rs:.4f} rs")
    print(f"c={c_amrs:.4f} am/rs   m_freq_kg(electron)={m_freq_kg:.4e} rad/rs")
    print(f"λ_φ⁴ = (c/dx)² = {lambda_phi4:.4e}  1/(am²·rs²)")
    print(f"λ·dt² = {lambda_phi4 * dt_rs**2:.4f}  (vs (c·dt/dx)² = "
          f"{(c_amrs * dt_rs / wf.dx_am)**2:.4f}; both < 4 for stability)")
    print(f"propagate {N_PROPAGATE_STEPS} steps, sample every {SAMPLE_EVERY}")
    print()

    # =================================================================
    # Run A — V = 0 (baseline; same as Step 3 [A])
    # =================================================================
    print("[A] BASELINE — V = 0 (free wave)")
    print(f"    seeding hedgehog Q=+1 at center...")
    center = seed_single_hedgehog(wf)
    samples_free = propagate_and_measure(wf, obs, c_amrs, dt_rs, 0.0, 0.0, center)
    print_table("free-wave Q(t):", samples_free)
    print()

    # =================================================================
    # Run B — Mexican-hat φ⁴ (with electron KG mass; KG invisible here)
    # =================================================================
    print("[B] MEXICAN-HAT φ⁴ — V = ½·m_e²·|ψ|² + ¼·λ·(|ψ|² − 1)²")
    print(f"    seeding hedgehog Q=+1 at center...")
    center = seed_single_hedgehog(wf)
    samples_phi4 = propagate_and_measure(wf, obs, c_amrs, dt_rs, m_freq_kg, lambda_phi4, center)
    print_table("φ⁴ Q(t):", samples_phi4)
    print()

    # =================================================================
    # Comparison
    # =================================================================
    print("=" * 78)
    print("COMPARISON — free vs Mexican-hat φ⁴")
    print("=" * 78)
    print(f"  {'step':>6s}  {'Q_free':>10s}  {'Q_φ⁴':>10s}  "
          f"{'<|ψ|>_free':>12s}  {'<|ψ|>_φ⁴':>12s}")
    for (s_f, Q_f, pn_f, _), (s_p, Q_p, pn_p, _) in zip(samples_free, samples_phi4):
        print(f"  {s_f:>6d}  {Q_f:>+10.4f}  {Q_p:>+10.4f}  "
              f"{pn_f:>12.5f}  {pn_p:>12.5f}")
    print()

    # =================================================================
    # Interpretation — find the first step where Q drops below 0.5
    # =================================================================
    def first_decay_step(samples, threshold=0.5):
        for step, Q, _, _ in samples:
            if abs(Q) < threshold:
                return step
        return None

    decay_free = first_decay_step(samples_free)
    decay_phi4 = first_decay_step(samples_phi4)
    psi_final_phi4 = samples_phi4[-1][2]
    psi_final_free = samples_free[-1][2]

    print("=" * 78)
    print("INTERPRETATION")
    print("=" * 78)
    print(f"  Q decay step (|Q| < 0.5):  free = {decay_free}  φ⁴ = {decay_phi4}")
    print(f"  Final <|ψ|>:               free = {psi_final_free:.4f}   "
          f"φ⁴ = {psi_final_phi4:.4f}")
    print()
    if decay_phi4 is not None and decay_free is not None:
        if decay_phi4 > decay_free:
            print("  ✅ φ⁴ EXTENDS defect lifetime — Mexican-hat shape helps slow Derrick collapse.")
        elif decay_phi4 < decay_free:
            print("  ⚠️ φ⁴ ACCELERATES decay — likely numerical (λ too large?).")
        else:
            print("  ⚠️ Identical decay — φ⁴ not measurably extending lifetime.")
    print()
    print(f"  <|ψ|> proximity to 1: φ⁴ should hold |ψ| near 1 (Mexican-hat minimum).")
    print(f"  Free baseline shows |ψ| drift; φ⁴ should hold tighter.")
    print()
    print("  NEXT (if φ⁴ alone is insufficient): Step 4b adds the 4th-derivative")
    print("  Skyrme term `+ ½κ |∇ψ × ∇ψ|²` to prevent core collapse — classical")
    print("  Skyrme model historically known to support 3D hedgehog solitons.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
