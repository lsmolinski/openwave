"""
M5.5.2 — first numerical evolution of the Eq.18 action (twist sector, 3D, headless)

Validates the action-derived dynamical evolution MACHINERY + energy conservation —
the first time Duda's Eq.18 Lagrangian is time-stepped numerically (M5.4 only did
the static energy). Derived from ℒ = T − U with the confirmed kinetic term
(5a §9):

    F_μ0 = 2[M_μ, Ṁ]           ⇒   T = 4 Σ_μ ‖[M_μ, Ṁ]‖²_F          (kinetic)
    F_μν = 2[M_μ, M_ν]         ⇒   U = 4 Σ_{μ<ν}‖[M_μ,M_ν]‖²_F + V   (potential)

DYNAMICAL DoF = the twist phase ψ(x,t): O(x,t) = O_bg(x)·Rx(ψ),  M = O D O^T.
With M_μ = M_μ^bg + ψ_μ M_ψ  (M_ψ = O_bg[Gx,D]O_bg^T), the twist is dynamical only
where the background varies (5a §9 degeneracy finding):

    P_μ  = [M_μ^bg, M_ψ]
    K(x) = 4 Σ_μ ‖P_μ‖²_F                                 (field-dependent kinetic metric)
    F̃_μν = C_μν − ψ_μ P_ν + ψ_ν P_μ ,   C_μν = [M_μ^bg, M_ν^bg]
    U    = 16 Σ_{μ<ν} ‖F̃_μν‖²_F
    EOM (V=0):  2K ψ_tt = Σ_μ ∂_μ J_μ ,  J_μ = −32 Σ_ν F̃_μν • P_ν

BACKGROUND: a SMOOTH, NON-SINGULAR director texture O_bg = Ry(γ(x)), γ a localized
3D bump → K ∝ ‖∇γ‖² (twist dynamical in the bump shell, vacuum elsewhere). This
validates the machinery cleanly. The KG MASS GAP needs the biaxial HEDGEHOG (C_μν≠0
from multiple generators) — which carries a z-axis DISCLINATION singularity → M5.6.

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v5.m5_5_2_twist_evolution
"""

import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

N = 32                      # cubic grid
L = 6.0                     # domain half-extent (so coords in [-L, L])
DELTA = 0.3                 # D = diag(1, δ, 0); twist lives on the δ axis
GAMMA0 = 0.8                # background tilt amplitude
R0, WID = 2.2, 1.2          # background bump: shell radius + width
DT = 0.02
N_STEPS = 1500

Gx = np.array([[0, 0, 0], [0, 0, -1], [0, 1, 0]], float)   # twist generator
Gy = np.array([[0, 0, 1], [0, 0, 0], [-1, 0, 0]], float)   # background tilt generator
D = np.diag([1.0, DELTA, 0.0])


def comm(A, B):
    return A @ B - B @ A


def frob2(M):
    """Frobenius² over the last two axes → scalar field."""
    return np.einsum("...ab,...ab->...", M, M)


def frob_inner(A, B):
    return np.einsum("...ab,...ab->...", A, B)


def matmul(A, B):
    return np.einsum("...ac,...cb->...ab", A, B)


def commf(A, B):
    return matmul(A, B) - matmul(B, A)


def central(f, axis, h):
    """Central difference of a scalar/tensor field along a spatial axis (interior; 0 at faces)."""
    out = np.zeros_like(f)
    sl_p = [slice(None)] * f.ndim
    sl_m = [slice(None)] * f.ndim
    sl_o = [slice(None)] * f.ndim
    sl_p[axis] = slice(2, None)
    sl_m[axis] = slice(0, -2)
    sl_o[axis] = slice(1, -1)
    out[tuple(sl_o)] = (f[tuple(sl_p)] - f[tuple(sl_m)]) / (2 * h)
    return out


def build_background():
    xs = np.linspace(-L, L, N)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    gamma = GAMMA0 * np.exp(-((R - R0) / WID) ** 2)          # smooth localized shell tilt
    c, s = np.cos(gamma), np.sin(gamma)
    O = np.zeros((N, N, N, 3, 3))
    O[..., 0, 0], O[..., 0, 2] = c, s                         # Ry(γ)
    O[..., 1, 1] = 1.0
    O[..., 2, 0], O[..., 2, 2] = -s, c
    OT = np.swapaxes(O, -1, -2)
    Mbg = matmul(matmul(O, np.broadcast_to(D, O.shape)), OT)
    Mpsi = matmul(matmul(O, np.broadcast_to(comm(Gx, D), O.shape)), OT)
    Mmu = [central(Mbg, ax, h) for ax in range(3)]            # M_μ^bg
    P = [commf(Mmu[ax], Mpsi) for ax in range(3)]             # P_μ = [M_μ^bg, M_ψ]
    K = 4.0 * sum(frob2(P[ax]) for ax in range(3))
    pairs = [(0, 1), (0, 2), (1, 2)]
    C = {p: commf(Mmu[p[0]], Mmu[p[1]]) for p in pairs}       # C_μν = [M_μ^bg, M_ν^bg]
    return dict(h=h, K=K, P=P, C=C, pairs=pairs, gamma=gamma, R=R)


def force_and_U(psi, bg):
    """Action-derived force Σ_μ ∂_μ J_μ and potential density U = 16 Σ‖F̃_μν‖²."""
    h, P, C, pairs = bg["h"], bg["P"], bg["C"], bg["pairs"]
    dpsi = [central(psi, ax, h) for ax in range(3)]
    Ftil = {}
    for (mu, nu) in pairs:
        Ftil[(mu, nu)] = (C[(mu, nu)]
                          - dpsi[mu][..., None, None] * P[nu]
                          + dpsi[nu][..., None, None] * P[mu])
    U = 16.0 * sum(frob2(Ftil[p]) for p in pairs)
    # J_μ = −32 Σ_ν F̃_μν • P_ν  (antisymmetric F̃: F̃_νμ = −F̃_μν)
    J = [np.zeros_like(psi) for _ in range(3)]
    J[0] = -32.0 * (frob_inner(Ftil[(0, 1)], P[1]) + frob_inner(Ftil[(0, 2)], P[2]))
    J[1] = -32.0 * (-frob_inner(Ftil[(0, 1)], P[0]) + frob_inner(Ftil[(1, 2)], P[2]))
    J[2] = -32.0 * (-frob_inner(Ftil[(0, 2)], P[0]) - frob_inner(Ftil[(1, 2)], P[1]))
    force = sum(central(J[ax], ax, h) for ax in range(3))
    return force, U


def main():
    print("=" * 70)
    print("M5.5.2 — numerical evolution of the Eq.18 action (twist sector, 3D)")
    print("=" * 70)
    bg = build_background()
    K = bg["K"]
    Kfloor = 1e-6 * K.max()
    active = K > Kfloor

    # ---- Stage 1: kinetic metric K(x) — dynamical twist confirmed ----
    print("\n[1] kinetic metric  K(x) = 4Σ‖[M_μ^bg, M_ψ]‖²_F")
    print(f"    K: max={K.max():.3e}  mean(active)={K[active].mean():.3e}  "
          f"active voxels={100*active.mean():.1f}% (∝ ‖∇γ‖², the bump shell)")
    print(f"    → {'PASS' if K.max() > 0 else 'FAIL'}  (twist is dynamical where the background varies)")

    # ---- Stage 2: leapfrog evolution — energy conservation + stability ----
    print("\n[2] leapfrog  2K ψ_tt = Σ_μ ∂_μ J_μ   (seed a localized twist packet)")
    xs = np.linspace(-L, L, N)
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    psi = 0.15 * np.exp(-((X - R0) ** 2 + Y**2 + Z**2) / 0.8)   # seed in the active shell
    psi *= active
    psi_prev = psi.copy()                                       # ψ̇=0 start
    inv2K = np.where(active, 1.0 / (2.0 * K + 1e-30), 0.0)

    def energy(p, p_prev):
        _, U = force_and_U(p, bg)
        pt = (p - p_prev) / DT
        return float(np.sum(K * pt**2 + U))

    H0 = energy(psi, psi_prev)
    Hs, amps = [], []
    for step in range(N_STEPS):
        force, _ = force_and_U(psi, bg)
        psi_tt = force * inv2K
        psi_new = (2 * psi - psi_prev + DT**2 * psi_tt) * active   # Dirichlet outside active
        psi_prev, psi = psi, psi_new
        if step % 50 == 0:
            Hs.append(energy(psi, psi_prev))
            amps.append(float(np.abs(psi).max()))

    Hs = np.array(Hs)
    drift = (Hs.max() - Hs.min()) / abs(H0) if H0 != 0 else 0.0
    blew_up = not np.isfinite(psi).all() or np.abs(psi).max() > 50 * 0.15
    print(f"    H0={H0:.4e}   energy drift over {N_STEPS} steps = {100*drift:.2f}%")
    print(f"    max|ψ|: seed=0.150 → final={np.abs(psi).max():.3f}  (stable, no blow-up: {not blew_up})")
    osc = np.std(amps) > 1e-4 * np.mean(amps) and not blew_up   # the packet evolves (not frozen)
    print(f"    twist packet evolves (not frozen): {osc}")
    ok = (K.max() > 0) and (drift < 0.05) and (not blew_up) and osc
    print(f"    → {'PASS' if ok else 'PARTIAL/FAIL'}")

    # ---- plot ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    PLOT = HERE / "plots" / "m5_5_2_twist_evolution.png"
    PLOT.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    mid = N // 2
    im0 = ax[0].imshow(K[:, :, mid].T, origin="lower"); ax[0].set_title("K(x) ∝ ‖∇γ‖²  (z=0 slice)")
    plt.colorbar(im0, ax=ax[0], fraction=0.046)
    ax[1].plot(np.arange(len(Hs)) * 50, Hs, "o-"); ax[1].axhline(H0, ls="--", c="k", alpha=0.5)
    ax[1].set_title(f"energy H (drift {100*drift:.2f}%)"); ax[1].set_xlabel("step"); ax[1].grid(alpha=0.3)
    ax[2].plot(np.arange(len(amps)) * 50, amps, "s-"); ax[2].set_title("max|ψ| (stability)")
    ax[2].set_xlabel("step"); ax[2].grid(alpha=0.3)
    plt.tight_layout(); plt.savefig(PLOT, dpi=110); plt.close()
    print(f"    plot → {PLOT}")

    print("\n" + "=" * 70)
    print(f"M5.5.2 MACHINERY: K(x)>0 ✓ | energy drift {100*drift:.2f}% | stable {not blew_up}")
    print("PASS — Eq.18 action evolves numerically + conserves energy" if ok else "PARTIAL/FAIL")
    print("(KG mass gap needs the biaxial hedgehog + disclination handling → M5.6)")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
