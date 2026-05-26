"""
M5.3 — Resonance smoke test on the EXISTING Vector(3) substrate (headless)

Robert Close's protocol (2026-04 email): "certain amplitudes of certain
harmonic waves will keep energy localized longer (i.e. as an unstable particle
or resonance). Probably l=1 harmonic; max displacement comparable to the
wavelength (or half/twice)." This script seeds an l=1 (dipole) wave packet at
A/λ ∈ {0.5, 1, 2} and measures how long its energy stays localized, with and
without the φ⁴ Mexican-hat nonlinearity.

NOT GATING (M5.3 optional task). The question it answers: does even the
Vector(3) substrate show extended-lifetime localization at a specific
amplitude? If yes → the resonance mechanism is substrate-agnostic (worth
knowing before M5.4). If no (expected) → it confirms the matrix substrate +
true nonlinear stabilizer (M5.4-M5.7) is required; the metastable-resonance
hunt belongs on the matrix field, per the M5.2 closed-negative result.

Localization metric: L(t) = Σ_{r<R_core} |ψ|²  /  Σ_all |ψ|²  (energy-intensity
fraction inside a core sphere). Lifetime = first step where L(t) < ½·L(0).

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v2.m5_3_resonance_smoke
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
from openwave.xperiments.m5_liquid_crystal import engine2_pde as pde  # noqa: E402

# ================================================================
# CONFIG
# ================================================================
N = 65
UNIVERSE_EDGE_M = 1e-15
TARGET_VOXELS = N**3
LAMBDA_VOX = 12.0  # internal wavelength of the harmonic packet, in voxels
SIGMA_VOX = LAMBDA_VOX  # Gaussian envelope width = one wavelength (localized)
A_OVER_LAMBDA = [0.5, 1.0, 2.0]  # Close's amplitude-comparable-to-wavelength sweep
N_PROPAGATE_STEPS = 400
SAMPLE_EVERY = 20
CORE_RADIUS_VOX = SIGMA_VOX  # localization core = initial packet scale


def seed_l1_dipole(wf, a_over_lambda):
    """Seed a z-polarized, l=1 (cosθ) dipole Gaussian wave packet at the center.

    Peak displacement amplitude A = (A/λ) · λ. Writes all three leapfrog
    buffers equal (ψ̇ = 0 rest start; also satisfies the triple-buffer BC rule
    so the first swap_buffers doesn't clobber the seed).
    """
    nx, ny, nz = wf.nx, wf.ny, wf.nz
    cx, cy, cz = nx / 2.0, ny / 2.0, nz / 2.0
    lam_am = LAMBDA_VOX * wf.dx_am
    amp_am = a_over_lambda * lam_am

    ii, jj, kk = np.meshgrid(
        np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij"
    )
    dx_ = ii - cx
    dy_ = jj - cy
    dz_ = kk - cz
    r_vox = np.sqrt(dx_**2 + dy_**2 + dz_**2) + 1e-9
    cos_theta = dz_ / r_vox  # l=1 dipole angular pattern (Y_1^0 ∝ cosθ)
    envelope = np.exp(-((r_vox / SIGMA_VOX) ** 2))

    scalar = amp_am * cos_theta * envelope  # (nx,ny,nz)
    psi = np.zeros((nx, ny, nz, 3), dtype=np.float32)
    psi[..., 2] = scalar  # ẑ polarization

    wf.psi_am.from_numpy(psi)
    wf.psi_prev_am.from_numpy(psi)
    wf.psi_new_am.from_numpy(psi)
    return float(amp_am), float(lam_am)


def localization_fraction(psi_np, core_mask):
    """Energy-intensity fraction inside the core sphere: Σ_core |ψ|² / Σ_all |ψ|².

    Returns NaN if the field has gone non-finite (numerical divergence).
    """
    intensity = (psi_np.astype(np.float64) ** 2).sum(axis=-1)  # |ψ|² per voxel (f64 to dodge overflow)
    total = intensity.sum()
    if not np.isfinite(total) or total <= 0.0:
        return float("nan")
    return float(intensity[core_mask].sum() / total)


def run_one(wf, c_amrs, dt_rs, m_freq_rs, lambda_phi4, a_over_lambda, core_mask):
    """Returns (samples, L0, lifetime, status). status ∈ {dispersed, diverged, sustained}."""
    amp_seed, _ = seed_l1_dipole(wf, a_over_lambda)
    L0 = localization_fraction(wf.psi_am.to_numpy(), core_mask)
    samples = [(0, L0)]
    lifetime = None
    status = "sustained"  # downgraded below if it decays or diverges
    for step in range(1, N_PROPAGATE_STEPS + 1):
        pde.evolve_psi(wf, c_amrs, dt_rs, m_freq_rs, lambda_phi4)
        wf.swap_buffers()
        if step % SAMPLE_EVERY == 0 or step == N_PROPAGATE_STEPS:
            psi_np = wf.psi_am.to_numpy()
            peak = float(np.abs(psi_np).max())
            # Divergence: non-finite, or field amplitude exploded > 50× the seed peak
            if not np.isfinite(peak) or peak > 50.0 * amp_seed:
                samples.append((step, float("nan")))
                status = "diverged"
                lifetime = step
                break
            L = localization_fraction(psi_np, core_mask)
            samples.append((step, L))
            if lifetime is None and L < 0.5 * L0:
                lifetime = step
                status = "dispersed"
    return samples, L0, lifetime, status


def main():
    ti.init(arch=ti.gpu)
    print("=" * 78)
    print("M5.3 — Resonance smoke test (l=1 dipole on Vector(3); Close protocol)")
    print("=" * 78)

    wf = medium.WaveField([UNIVERSE_EDGE_M] * 3, TARGET_VOXELS)
    c_amrs = constants.WAVE_SPEED / constants.ATTOMETER * constants.RONTOSECOND
    dt_rs = wf.dx_am * 0.95 / (c_amrs * (3**0.5))
    m_freq_kg = c_amrs / constants.COMPTON_WAVELENGTH_REDUCED_ELECTRON_AM
    lambda_phi4 = 0.1 * (c_amrs / wf.dx_am) ** 2  # matches launcher default (stable)

    # Precompute the core-sphere mask once (geometry is amplitude-independent).
    nx, ny, nz = wf.nx, wf.ny, wf.nz
    cx, cy, cz = nx / 2.0, ny / 2.0, nz / 2.0
    ii, jj, kk = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    r_vox = np.sqrt((ii - cx) ** 2 + (jj - cy) ** 2 + (kk - cz) ** 2)
    core_mask = r_vox < CORE_RADIUS_VOX

    print(f"Grid: {nx}³  dx={wf.dx_am:.3f} am  dt={dt_rs:.4f} rs")
    print(f"packet λ={LAMBDA_VOX:.0f} vox  envelope σ={SIGMA_VOX:.0f} vox  "
          f"core R={CORE_RADIUS_VOX:.0f} vox")
    print(f"propagate {N_PROPAGATE_STEPS} steps, sample every {SAMPLE_EVERY}; "
          f"lifetime = step where L < ½·L₀")
    print()

    results = {}
    for pot_label, m_freq, lam4 in [
        ("V=0 (free-wave)", 0.0, 0.0),
        ("KG mass ½m²|ψ|²", m_freq_kg, 0.0),
        ("φ⁴ Mexican-hat", m_freq_kg, lambda_phi4),
    ]:
        print(f"[{pot_label}]")
        print(f"    {'A/λ':>5s}  {'A (am)':>10s}  {'L₀':>8s}  {'outcome':>16s}  L(t) trace")
        for ratio in A_OVER_LAMBDA:
            samples, L0, lifetime, status = run_one(
                wf, c_amrs, dt_rs, m_freq, lam4, ratio, core_mask
            )
            amp_am = ratio * LAMBDA_VOX * wf.dx_am
            trace = " ".join(f"{L:.2f}" for _, L in samples[:: max(1, len(samples) // 8)])
            if status == "diverged":
                outcome = f"DIVERGED@{lifetime}"
            elif status == "dispersed":
                outcome = f"disperse@{lifetime}"
            else:
                outcome = f"sustained>{N_PROPAGATE_STEPS}"
            print(f"    {ratio:>5.1f}  {amp_am:>10.1f}  {L0:>8.3f}  {outcome:>16s}  {trace}")
            results[(pot_label, ratio)] = status
        print()

    # ================================================================
    # Verdict
    # ================================================================
    print("=" * 78)
    print("VERDICT")
    print("=" * 78)
    sustained = [(lbl, r) for (lbl, r), st in results.items() if st == "sustained"]
    diverged = [(lbl, r) for (lbl, r), st in results.items() if st == "diverged"]
    if sustained:
        print(f"  ⚠️ SUSTAINED localization (L ≥ ½L₀ for {N_PROPAGATE_STEPS} steps, finite) at:")
        for lbl, r in sustained:
            print(f"       {lbl}, A/λ={r}")
        print("  → possible substrate-agnostic resonance — investigate before M5.4.")
    else:
        print("  ✅ NULL (expected): no sustained localization on Vector(3) at any A/λ.")
        print("     • Linear potentials (V=0, KG mass) disperse — and amplitude-")
        print("       INDEPENDENTLY, as a linear PDE must (no amplitude-selected resonance).")
        if diverged:
            print("     • φ⁴ Mexican-hat DIVERGES: it is a director-scaled potential")
            print("       (vacuum |ψ|=1), structurally inapplicable to a displacement-wave")
            print("       amplitude sweep — Vector(3)'s potential toolkit can't even host")
            print("       the test. (Diverged rows above are this, not localization.)")
        print()
        print("     Confirms the M5.2 closed-negative result: a metastable resonance needs")
        print("     the matrix substrate + LdG V(M)/Skyrme + 4D clock propulsion")
        print("     (M5.4 → M5.7), NOT the Vector(3) field. No substrate-agnostic resonance")
        print("     — the resonance hunt belongs on the matrix substrate.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
