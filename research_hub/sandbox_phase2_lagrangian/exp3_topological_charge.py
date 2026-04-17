"""
Experiment 3: Topological Charge Quantization

Numerical Experiment 3 in the Lagrangian Framework Sub-Project.

Hypothesis:
The winding-number integral
    Q = (1/4π) ∮_S (∂_θ n × ∂_φ n) · n  dθ dφ
returns integers (±1 for hedgehog / anti-hedgehog, 0 for vacuum)
regardless of surface shape, surface radius, or smooth perturbation of n.

Tests:
1. Clean hedgehog        → Q = +1   (surface radius: 2, 3, 5, 7 grid units)
2. Clean anti-hedgehog   → Q = −1   (same radii)
3. Vacuum (n = ẑ)        → Q =  0
4. Perturbed hedgehog    → Q = +1   (noise levels: 5%, 10%, 20%, 50%)
5. Two-defect config     → Q_sphere_around_each = ±1, Q_sphere_around_both = 0

Spec:    ../3_LAGRANGIAN_FRAMEWORK.md  (Experiment 3)
Results: ../3b_lagrangian_experiments.md  (Experiment 3)
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# GRID / FIELD
# ---------------------------------------------------------------------------

N = 48
DOMAIN = 16.0
DX = DOMAIN / (N - 1)

RESULTS_DIR = Path(__file__).parent / "exp3_results"

_X = _Y = _Z = None


def make_grid():
    global _X, _Y, _Z
    if _X is None:
        x = np.linspace(-DOMAIN / 2, DOMAIN / 2, N)
        _X, _Y, _Z = np.meshgrid(x, x, x, indexing="ij")
    return _X, _Y, _Z


def build_hedgehog_field(center=(0.0, 0.0, 0.0), sign=+1):
    """3D director field of a clean hedgehog at `center`."""
    X, Y, Z = make_grid()
    dx, dy, dz = X - center[0], Y - center[1], Z - center[2]
    r = np.sqrt(dx**2 + dy**2 + dz**2)
    r_safe = np.maximum(r, 0.2 * DX)
    return np.stack([sign * dx / r_safe, sign * dy / r_safe, sign * dz / r_safe], axis=-1)


def build_vacuum_field():
    X, Y, _ = make_grid()
    n = np.zeros((N, N, N, 3), dtype=np.float64)
    n[..., 2] = 1.0
    return n


def build_two_hedgehogs(d=4.0):
    """Hedgehog at -d/2 (Q=+1) + anti-hedgehog at +d/2 (Q=-1), via weighted blend."""
    X, Y, Z = make_grid()
    c1 = np.array([-d / 2.0, 0.0, 0.0])
    c2 = np.array([+d / 2.0, 0.0, 0.0])
    dx1, dy1, dz1 = X - c1[0], Y - c1[1], Z - c1[2]
    dx2, dy2, dz2 = X - c2[0], Y - c2[1], Z - c2[2]
    r1 = np.sqrt(dx1**2 + dy1**2 + dz1**2)
    r2 = np.sqrt(dx2**2 + dy2**2 + dz2**2)
    r1s = np.maximum(r1, 0.2 * DX)
    r2s = np.maximum(r2, 0.2 * DX)
    n1 = np.stack([dx1 / r1s, dy1 / r1s, dz1 / r1s], axis=-1)
    n2 = np.stack([-dx2 / r2s, -dy2 / r2s, -dz2 / r2s], axis=-1)
    w1 = 1.0 / (r1 + 0.5)
    w2 = 1.0 / (r2 + 0.5)
    n = w1[..., None] * n1 + w2[..., None] * n2
    norm = np.linalg.norm(n, axis=-1, keepdims=True)
    return n / np.maximum(norm, 1e-12), c1, c2


def perturb_field(n, noise_amplitude, seed=0):
    """Add uniform random noise to each component, renormalize."""
    rng = np.random.default_rng(seed)
    noise = rng.uniform(-noise_amplitude, noise_amplitude, size=n.shape)
    n_pert = n + noise
    norm = np.linalg.norm(n_pert, axis=-1, keepdims=True)
    return n_pert / np.maximum(norm, 1e-12)


# ---------------------------------------------------------------------------
# WINDING NUMBER INTEGRAL ON A SPHERE
# ---------------------------------------------------------------------------

def trilinear_sample(n_field, x, y, z):
    """Trilinear interpolation of a vector field at (x,y,z) physical coordinates.

    Returns an array of shape (..., 3).
    """
    # Physical → grid index
    ix = (x + DOMAIN / 2) / DX
    iy = (y + DOMAIN / 2) / DX
    iz = (z + DOMAIN / 2) / DX

    i0 = np.floor(ix).astype(int)
    j0 = np.floor(iy).astype(int)
    k0 = np.floor(iz).astype(int)
    i0 = np.clip(i0, 0, N - 2)
    j0 = np.clip(j0, 0, N - 2)
    k0 = np.clip(k0, 0, N - 2)

    fx = ix - i0
    fy = iy - j0
    fz = iz - k0
    fx = np.clip(fx, 0.0, 1.0)
    fy = np.clip(fy, 0.0, 1.0)
    fz = np.clip(fz, 0.0, 1.0)

    out = np.zeros(x.shape + (3,), dtype=np.float64)
    for dxi in (0, 1):
        for dyi in (0, 1):
            for dzi in (0, 1):
                w = (fx if dxi else 1 - fx) * (fy if dyi else 1 - fy) * (fz if dzi else 1 - fz)
                out += w[..., None] * n_field[i0 + dxi, j0 + dyi, k0 + dzi, :]
    return out


def winding_number(n_field, center=(0.0, 0.0, 0.0), radius=2.0,
                   n_theta=64, n_phi=128):
    """Compute Q = (1/4π) ∫∫_S n · (∂_θ n × ∂_φ n) sinθ dθ dφ on sphere of radius R.

    Uses trilinear interpolation of the 3D field onto a spherical mesh,
    then finite differences on the (θ, φ) grid to form the surface integrand.
    """
    theta = np.linspace(1e-3, np.pi - 1e-3, n_theta)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    dth = theta[1] - theta[0]
    dph = phi[1] - phi[0]
    TH, PH = np.meshgrid(theta, phi, indexing="ij")
    sx = center[0] + radius * np.sin(TH) * np.cos(PH)
    sy = center[1] + radius * np.sin(TH) * np.sin(PH)
    sz = center[2] + radius * np.cos(TH)

    n_s = trilinear_sample(n_field, sx, sy, sz)  # shape (n_theta, n_phi, 3)

    # Finite differences on the (θ, φ) grid. Use central differences, wrap φ.
    dn_dtheta = np.zeros_like(n_s)
    dn_dphi = np.zeros_like(n_s)
    dn_dtheta[1:-1] = (n_s[2:] - n_s[:-2]) / (2 * dth)
    dn_dtheta[0] = (n_s[1] - n_s[0]) / dth
    dn_dtheta[-1] = (n_s[-1] - n_s[-2]) / dth
    dn_dphi[:, :] = (np.roll(n_s, -1, axis=1) - np.roll(n_s, 1, axis=1)) / (2 * dph)

    # Integrand: n · (∂_θ n × ∂_φ n)
    cross = np.cross(dn_dtheta, dn_dphi, axis=-1)
    integrand = (n_s * cross).sum(axis=-1)

    # Sphere measure: for unit-radius surface parameterization ∫∫ F(θ,φ) dθ dφ
    # (the sinθ factor is implicit in the cross product for this specific
    # form — n_s is unit, (∂_θ n × ∂_φ n) already contains the Jacobian.)
    return float(integrand.sum() * dth * dph / (4.0 * np.pi))


# ---------------------------------------------------------------------------
# TESTS
# ---------------------------------------------------------------------------

def test_clean_hedgehog():
    print("\n" + "=" * 72)
    print("TEST 1 — Clean hedgehog / anti-hedgehog / vacuum")
    print("=" * 72)
    n_h = build_hedgehog_field(sign=+1)
    n_a = build_hedgehog_field(sign=-1)
    n_v = build_vacuum_field()

    results = []
    for R in [2.0, 3.0, 5.0, 7.0]:
        qh = winding_number(n_h, radius=R)
        qa = winding_number(n_a, radius=R)
        qv = winding_number(n_v, radius=R)
        results.append((R, qh, qa, qv))
        print(f"  R={R:.1f}:  Q(hedgehog) = {qh:+.4f}   "
              f"Q(anti-hedgehog) = {qa:+.4f}   Q(vacuum) = {qv:+.4f}")
    return results


def test_perturbed_hedgehog():
    print("\n" + "=" * 72)
    print("TEST 2 — Perturbed hedgehog (noise on n, renormalized)")
    print("=" * 72)
    n_h = build_hedgehog_field(sign=+1)
    R = 3.0

    results = []
    for noise in [0.0, 0.05, 0.10, 0.20, 0.50, 1.0]:
        n_pert = perturb_field(n_h, noise, seed=42)
        q = winding_number(n_pert, radius=R)
        results.append((noise, q))
        print(f"  noise={noise:.2f}  Q = {q:+.4f}   "
              f"rounded → {int(round(q)):+d}   match? "
              f"{'✓' if abs(round(q) - 1) < 0.5 else '✗'}")
    return results


def test_two_hedgehogs():
    print("\n" + "=" * 72)
    print("TEST 3 — Two defects (Q=+1 and Q=−1) at separation d=4")
    print("=" * 72)
    n, c1, c2 = build_two_hedgehogs(d=4.0)

    # Sphere around defect 1 (should enclose Q=+1 only)
    q1 = winding_number(n, center=tuple(c1), radius=1.3)
    # Sphere around defect 2 (should enclose Q=-1 only)
    q2 = winding_number(n, center=tuple(c2), radius=1.3)
    # Large sphere enclosing both (should give total Q=0)
    q_both = winding_number(n, center=(0.0, 0.0, 0.0), radius=5.0)
    # Sphere far from both (should give Q=0)
    q_far = winding_number(n, center=(0.0, 0.0, 6.0), radius=1.0)

    print(f"  sphere around c1 (R=1.3):  Q = {q1:+.4f}   expected +1")
    print(f"  sphere around c2 (R=1.3):  Q = {q2:+.4f}   expected −1")
    print(f"  large sphere enclosing both (R=5):  Q = {q_both:+.4f}   expected  0")
    print(f"  sphere far from both (R=1):  Q = {q_far:+.4f}   expected  0")
    return q1, q2, q_both, q_far


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Grid {N}³,  domain=[{-DOMAIN/2}, {DOMAIN/2}],  dx={DX:.4f}")

    r1 = test_clean_hedgehog()
    r2 = test_perturbed_hedgehog()
    r3 = test_two_hedgehogs()

    # PLOTS
    # Plot 1: Q vs sphere radius for clean configs
    fig, ax = plt.subplots(figsize=(8, 4.5))
    Rs = [row[0] for row in r1]
    qh = [row[1] for row in r1]
    qa = [row[2] for row in r1]
    qv = [row[3] for row in r1]
    ax.plot(Rs, qh, "o-", label="hedgehog (expect +1)")
    ax.plot(Rs, qa, "s-", label="anti-hedgehog (expect −1)")
    ax.plot(Rs, qv, "^-", label="vacuum (expect 0)")
    ax.axhline(+1, ls=":", c="gray")
    ax.axhline(-1, ls=":", c="gray")
    ax.axhline(0, ls=":", c="gray")
    ax.set_xlabel("sphere radius R")
    ax.set_ylabel("Q")
    ax.set_title("Exp 3 — Winding number vs surface radius (clean configs)")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "winding_vs_radius.png", dpi=110)
    plt.close()

    # Plot 2: Q vs noise amplitude
    fig, ax = plt.subplots(figsize=(8, 4.5))
    noise_levels = [row[0] for row in r2]
    q_values = [row[1] for row in r2]
    ax.plot(noise_levels, q_values, "o-", ms=8, label="measured Q")
    ax.axhline(+1, ls="--", c="red", label="predicted Q=+1")
    ax.set_xlabel("noise amplitude")
    ax.set_ylabel("Q")
    ax.set_title("Exp 3 — Topological stability: Q vs noise")
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "winding_vs_noise.png", dpi=110)
    plt.close()

    print(f"\nPlots saved to: {RESULTS_DIR}")


if __name__ == "__main__":
    main()
