"""
Experiment 7: Close's Nonlinear Vector Wave Equation
(v2 — using Close's actual published equations)

Numerical Experiment 7 in the Lagrangian Framework Sub-Project.

Reference: Robert A. Close, "Plane Wave Solutions to a Proposed 'Equation
of Everything'", Foundations of Physics (2025) 55:27.

Hypothesis:
Close's nonlinear vector wave equation for spin density, seeded with a
spherical harmonic, should exhibit the particle-like dynamics Close
derives in his paper — plane-wave propagation in the linear limit, and
self-coupling (rotation/advection) in the nonlinear case.

Close's equations implemented here:
  Linear (Eq. 19):
      ∂²_t Q + c²·∇×∇×Q = 0
  Nonlinear "Equation of Everything" (Eq. 21), rewritten for Q:
      ∂²_t Q = −c²·∇×∇×Q − u·∇s + w×s
  where:
      s = ∂_t Q          (spin density; Eq. 17: s is the axial vector field)
      u = (1/2ρ)·∇×s      (medium velocity from s)
      w = (1/2)·∇×u       (angular velocity)

For plane waves, u·∇s and w×s vanish; for localized structures they
couple s to its own curl, which is Close's self-interaction.

Implementation notes (this is a faithful numerical port, NOT a proxy):
- Vector field Q(x, t) on 3D grid with 3 components per voxel
- Leapfrog for ∂²_t Q with nonlinear RHS
- s approximated at midpoint as (Q − Q_old)/dt
- curl, divergence, Laplacian all via central finite differences (np.gradient /
  np.roll)
- curl(curl(F)) computed via identity ∇(∇·F) − ∇²F for accuracy

Test plan:
1. Linear baseline (Eq. 19) with Y_l^m seeds → expect dispersion (no soliton
   without nonlinearity; should still behave as a transverse vector wave)
2. Full nonlinear (Eq. 21) with same seeds → looks for self-localization or
   rotational dynamics from the u·∇s − w×s coupling

Spec:    ../../1aa_lagrangian_framework.md  (Experiment 7)
Results: ../3b_lagrangian_experiments.md  (Experiment 7)
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# PARAMETERS
# ---------------------------------------------------------------------------

N = 48
DOMAIN = 16.0
DX = DOMAIN / (N - 1)
C = 1.0
RHO = 1.0                    # inertial density (Close's ρ)
DT = 0.02                    # CFL c·dt/dx = 0.06 — very stable
N_STEPS = 300                # t_final = 6
SEED_SIGMA = 2.0             # Gaussian envelope width
SEED_AMPLITUDE = 0.5

LM_MODES = [(1, 0), (2, 0), (2, 1)]
MODES = ["linear", "nonlinear"]

RESULTS_DIR = Path(__file__).parent / "exp7_results"


# ---------------------------------------------------------------------------
# GRID / DIFFERENTIAL OPERATORS
# ---------------------------------------------------------------------------

_X = _Y = _Z = _R = _THETA = _PHI = None


def make_grid():
    global _X, _Y, _Z, _R, _THETA, _PHI
    if _X is None:
        x = np.linspace(-DOMAIN / 2, DOMAIN / 2, N)
        _X, _Y, _Z = np.meshgrid(x, x, x, indexing="ij")
        _R = np.sqrt(_X**2 + _Y**2 + _Z**2)
        _THETA = np.arccos(np.clip(_Z / np.maximum(_R, 1e-12), -1.0, 1.0))
        _PHI = np.arctan2(_Y, _X)
    return _X, _Y, _Z, _R, _THETA, _PHI


def laplacian_vec(F, dx):
    """Component-wise 6-point Laplacian, periodic via np.roll."""
    lap = np.zeros_like(F)
    for i in range(3):
        f = F[..., i]
        lap[..., i] = (
            np.roll(f, 1, 0) + np.roll(f, -1, 0)
            + np.roll(f, 1, 1) + np.roll(f, -1, 1)
            + np.roll(f, 1, 2) + np.roll(f, -1, 2)
            - 6.0 * f
        ) / dx**2
    return lap


def divergence(F, dx):
    return (np.gradient(F[..., 0], dx, axis=0)
            + np.gradient(F[..., 1], dx, axis=1)
            + np.gradient(F[..., 2], dx, axis=2))


def grad_scalar(f, dx):
    return np.stack([
        np.gradient(f, dx, axis=0),
        np.gradient(f, dx, axis=1),
        np.gradient(f, dx, axis=2),
    ], axis=-1)


def curl(F, dx):
    """Vector curl via central differences."""
    out = np.zeros_like(F)
    out[..., 0] = (np.gradient(F[..., 2], dx, axis=1)
                   - np.gradient(F[..., 1], dx, axis=2))
    out[..., 1] = (np.gradient(F[..., 0], dx, axis=2)
                   - np.gradient(F[..., 2], dx, axis=0))
    out[..., 2] = (np.gradient(F[..., 1], dx, axis=0)
                   - np.gradient(F[..., 0], dx, axis=1))
    return out


def curl_curl(F, dx):
    """∇×∇×F  =  ∇(∇·F)  −  ∇²F   (vector calculus identity)."""
    return grad_scalar(divergence(F, dx), dx) - laplacian_vec(F, dx)


# ---------------------------------------------------------------------------
# INITIAL CONDITIONS
# ---------------------------------------------------------------------------

def spherical_harmonic(l, m, theta, phi):
    """Real-valued shape functions, amplitude ≤ 1."""
    if (l, m) == (1, 0):
        return np.cos(theta)
    if (l, m) == (2, 0):
        return 0.5 * (3.0 * np.cos(theta)**2 - 1.0)
    if (l, m) == (2, 1):
        return np.sin(theta) * np.cos(theta) * np.cos(phi)
    raise ValueError(f"(l, m) = ({l}, {m}) not implemented")


def seed_vector_sh(l, m):
    """Q_z(x, 0) = A · exp(−r²/(2σ²)) · Y_l^m(θ, φ); Q_x = Q_y = 0.
    Initial velocity s = ∂_t Q = 0 (stationary seed).
    """
    _, _, _, _, theta, phi = make_grid()
    r = np.sqrt(_X**2 + _Y**2 + _Z**2) if _X is not None else make_grid()[3]
    # Fix: always re-read r from cached grid
    _, _, _, r, _, _ = make_grid()
    Q = np.zeros((N, N, N, 3), dtype=np.float64)
    envelope = np.exp(-0.5 * (r / SEED_SIGMA)**2)
    Q[..., 2] = SEED_AMPLITUDE * envelope * spherical_harmonic(l, m, theta, phi)
    return Q


# ---------------------------------------------------------------------------
# EVOLUTION
# ---------------------------------------------------------------------------

def step_linear(Q, Q_old, dt, dx, c):
    """Eq. 19:  ∂²_t Q = −c²·∇×∇×Q.
    Leapfrog:  Q_new = 2Q − Q_old + dt²·accel.
    """
    cc = curl_curl(Q, dx)
    accel = -c**2 * cc
    return 2.0 * Q - Q_old + dt**2 * accel


def step_nonlinear(Q, Q_old, dt, dx, c, rho):
    """Eq. 21:  ∂²_t Q = −c²·∇×∇×Q  −  (u·∇)s  +  w×s
    with s = ∂_t Q,  u = (1/2ρ)·∇×s,  w = (1/2)·∇×u.
    """
    s = (Q - Q_old) / dt                 # backward time derivative
    cs = curl(s, dx)                     # ∇×s
    u = cs / (2.0 * rho)                 # velocity
    cu = curl(u, dx)                     # ∇×u
    w = 0.5 * cu                          # angular velocity

    # Advection term (u·∇) s  —  component-wise:  a_i = u_j · ∂_j s_i
    advection = np.zeros_like(s)
    for i in range(3):
        for j in range(3):
            advection[..., i] += u[..., j] * np.gradient(s[..., i], dx, axis=j)

    # Rotation w × s
    rotation = np.cross(w, s)

    cc = curl_curl(Q, dx)
    accel = -c**2 * cc - advection + rotation
    return 2.0 * Q - Q_old + dt**2 * accel


# ---------------------------------------------------------------------------
# DIAGNOSTICS
# ---------------------------------------------------------------------------

def energy_density(Q, s, dx, c):
    """Elastic-solid energy: ½ρ|s|²  +  ½μ|∇×Q|²  +  ½(λ+2μ)|∇·Q|².
    With μ = ρc² and taking incompressible limit (λ → ∞), the ∇·Q part
    usually vanishes; we keep it for completeness."""
    # kinetic
    kinetic = 0.5 * RHO * (s**2).sum(axis=-1)
    # shear / transverse gradient energy  (related to curl)
    curlQ = curl(Q, dx)
    shear = 0.5 * RHO * c**2 * (curlQ**2).sum(axis=-1)
    # compressional (proxy, if not incompressible)
    divQ = divergence(Q, dx)
    compression = 0.5 * RHO * c**2 * divQ**2   # assuming λ+2μ ≈ 2μ for simplicity
    return kinetic + shear + compression


def energy_in_ball(Q, s, radius, dx, c):
    _, _, _, r, _, _ = make_grid()
    H = energy_density(Q, s, dx, c)
    total = float(H.sum() * dx**3)
    inside = float(H[r < radius].sum() * dx**3)
    return inside / max(abs(total), 1e-30), total


def peak_Q(Q):
    return float(np.linalg.norm(Q, axis=-1).max())


# ---------------------------------------------------------------------------
# RUN ONE (mode, l, m)
# ---------------------------------------------------------------------------

def run_one(mode, l, m, n_steps=N_STEPS):
    Q = seed_vector_sh(l, m)
    Q_old = Q.copy()

    ball_radius = 2.0 * SEED_SIGMA
    step_fn = step_linear if mode == "linear" else step_nonlinear

    sample_stride = max(1, n_steps // 10)
    ts, energies, peaks, concs = [], [], [], []

    for step in range(n_steps + 1):
        s = (Q - Q_old) / DT
        if step % sample_stride == 0 or step == n_steps:
            Etot = float(energy_density(Q, s, DX, C).sum() * DX**3)
            conc, _ = energy_in_ball(Q, s, ball_radius, DX, C)
            ts.append(step * DT)
            energies.append(Etot)
            peaks.append(peak_Q(Q))
            concs.append(conc)
        if step < n_steps:
            try:
                if mode == "linear":
                    Q_new = step_linear(Q, Q_old, DT, DX, C)
                else:
                    Q_new = step_nonlinear(Q, Q_old, DT, DX, C, RHO)
            except FloatingPointError:
                print(f"  ** Instability at step {step}, aborting run **")
                break
            if np.any(~np.isfinite(Q_new)) or np.abs(Q_new).max() > 1e6:
                print(f"  ** Blow-up at step {step} (max|Q|={np.abs(Q_new).max():.2e}), aborting **")
                break
            Q_old = Q
            Q = Q_new

    return {
        "mode": mode, "l": l, "m": m,
        "t": np.array(ts),
        "energy": np.array(energies),
        "peak": np.array(peaks),
        "concentration": np.array(concs),
        "final_Q": Q,
    }


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Grid:  {N}³,  dx={DX:.4f},  domain=[{-DOMAIN/2}, {DOMAIN/2}]")
    print(f"Time:  dt={DT}, N_STEPS={N_STEPS}, t_final={N_STEPS*DT}")
    print(f"CFL:   c·dt/dx = {C*DT/DX:.3f}  (need < 1)")
    print(f"Physics: c={C}, ρ={RHO}")
    print(f"Seed: σ={SEED_SIGMA}, A={SEED_AMPLITUDE}")
    print(f"Equations implemented: Eq. 19 (linear) and Eq. 21 (nonlinear)")
    print()

    results = {}
    for mode in MODES:
        for (l, m) in LM_MODES:
            print(f"Running mode={mode:>10s}, Y_{l}^{m}...")
            r = run_one(mode, l, m)
            results[(mode, l, m)] = r
            print(f"  concentration {r['concentration'][0]:.4f} → "
                  f"{r['concentration'][-1]:.4f}")
            print(f"  peak |Q|      {r['peak'][0]:.4f} → {r['peak'][-1]:.4f}")
            drift = ((r["energy"].max() - r["energy"].min())
                     / max(abs(r["energy"][0]), 1e-30))
            print(f"  E drift       {drift:.2e}")
            print()

    # ----- PLOTS -----
    # 2 panels × 3 columns: (linear vs nonlinear) × (conc, peak, energy)
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    for row, mode in enumerate(MODES):
        for (l, m) in LM_MODES:
            r = results[(mode, l, m)]
            axs[row, 0].plot(r["t"], r["concentration"], "o-", ms=4,
                             label=f"Y_{l}^{m}")
            axs[row, 1].plot(r["t"], r["peak"], "o-", ms=4,
                             label=f"Y_{l}^{m}")
            axs[row, 2].plot(r["t"], r["energy"], "o-", ms=4,
                             label=f"Y_{l}^{m}")
        axs[row, 0].set_ylabel("energy-in-ball")
        axs[row, 1].set_ylabel("peak |Q|")
        axs[row, 2].set_ylabel("total H")
        for c in range(3):
            axs[row, c].set_title(f"{mode.capitalize()}")
            axs[row, c].grid(alpha=0.3)
        axs[row, 0].legend(fontsize=8)
    axs[1, 0].set_xlabel("t"); axs[1, 1].set_xlabel("t"); axs[1, 2].set_xlabel("t")
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "close_dynamics.png", dpi=110)
    plt.close()

    # Cross-section of |Q| for (nonlinear, Y_1^0) at final time
    r = results[("nonlinear", 1, 0)]
    Q_final = r["final_Q"]
    Q_mag = np.linalg.norm(Q_final, axis=-1)
    mid = N // 2
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))
    vmax = max(Q_mag.max(), 1e-10)
    for ax, (slc, name) in zip(axs, [
        (Q_mag[:, :, mid], "xy"),
        (Q_mag[:, mid, :], "xz"),
        (Q_mag[mid, :, :], "yz"),
    ]):
        im = ax.imshow(slc, origin="lower",
                       extent=[-DOMAIN/2, DOMAIN/2, -DOMAIN/2, DOMAIN/2],
                       vmin=0.0, vmax=vmax)
        ax.set_title(f"|Q| nonlinear Y_1^0 ({name} slice, t={r['t'][-1]})")
        plt.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "close_final_field.png", dpi=110)
    plt.close()

    print(f"Plots saved to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
