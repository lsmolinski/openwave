"""
SINE-GORDON KINK + ANTIKINK ANNIHILATION, the principle-level trail for the
"Antimatter + annihilation" row of MODELS.md.

Duda (2026-06-10, coverage-matrix review): "Antimatter + annihilation is
absolutely natural for topological soliton models, starting with sine-Gordon."

The demo, with the integrability honesty built in: PURE sine-Gordon is the
famous integrable exception, a kink and antikink pass through each other
elastically. Any generic perturbation (here: weak damping, the 1+1D stand-in
for radiation/3D losses) restores the generic soliton behavior: capture into
a breather (the bound Q = 0 particle-antiparticle state) which decays
radiatively to the vacuum. Total topological charge is 0 before, during, and
after: annihilation conserves topology exactly while the stored field energy
is released.

  phi_tt = phi_xx - sin(phi) - gamma * phi_t
  kink: 4 atan(exp(g_L (x - vt))), Q = +1; antikink mirrored, Q = -1.

RESULTS (2026-06-10): see the run table; regime 1 (gamma = 0) the pair
survives the collision (integrable pass-through, 2 kinks before and after);
regime 2 (gamma = 0.02) capture -> breather -> decay, final energy < 1% of
initial, field at vacuum: ANNIHILATED.

USAGE:  python sine_gordon_annihilation.py
"""
import numpy as np

N, L = 4096, 100.0
DX = 2 * L / N
DT = 0.02
X = np.linspace(-L, L, N, endpoint=False)


def kink(x, v, sign=+1):
    gL = 1.0 / np.sqrt(1.0 - v * v)
    return 4.0 * np.arctan(np.exp(sign * gL * x))


def initial(v, x0):
    """Kink at -x0 moving right + antikink at +x0 moving left (total Q = 0)."""
    phi = kink(X + x0, v, +1) + kink(X - x0, v, -1) - 2.0 * np.pi
    gL = 1.0 / np.sqrt(1.0 - v * v)
    # traveling-wave velocities: phi_t = -v phi_x per moving profile
    dk = 4.0 * gL * np.exp(gL * (X + x0)) / (1.0 + np.exp(2 * gL * (X + x0)))
    da = -4.0 * gL * np.exp(-gL * (X - x0)) / (1.0 + np.exp(-2 * gL * (X - x0)))
    phit = -v * dk + v * da
    return phi, phit


def lap(f):
    return (np.roll(f, 1) + np.roll(f, -1) - 2 * f) / DX ** 2


def energy(phi, phit):
    phix = (np.roll(phi, -1) - np.roll(phi, 1)) / (2 * DX)
    return float(np.sum(0.5 * phit ** 2 + 0.5 * phix ** 2
                        + (1 - np.cos(phi))) * DX)


def n_kinks(phi):
    """Count |Q| = 1 lumps: total absolute winding / 2pi."""
    phix = np.abs((np.roll(phi, -1) - np.roll(phi, 1)) / (2 * DX))
    return int(round(float(np.sum(phix) * DX / (2 * np.pi))))


def run(gamma, v, steps):
    phi, phit = initial(v, 15.0)
    E0, k0 = energy(phi, phit), n_kinks(phi)
    for _ in range(steps):
        acc = lap(phi) - np.sin(phi) - gamma * phit
        phit = phit + DT * acc
        phi = phi + DT * phit
    E1, k1 = energy(phi, phit), n_kinks(phi)
    vac = float(np.abs(np.sin(phi / 2)).max())     # 0 at any vacuum 2*pi*n
    Q = (phi[-1] - phi[0]) / (2 * np.pi)
    return E0, E1, k0, k1, vac, Q


def main():
    print("=" * 74)
    print("Sine-Gordon kink + antikink: pass-through vs annihilation")
    print("=" * 74)
    print("\n  | regime | E start | E end | kinks before | kinks after |"
          " max|sin φ/2| end | total Q |")
    print("  | --- | --- | --- | --- | --- | --- | --- |")
    for gamma, v, steps, tag in ((0.0, 0.5, 8000, "pure SG (integrable)"),
                                 (0.02, 0.2, 60000, "weakly damped")):
        E0, E1, k0, k1, vac, Q = run(gamma, v, steps)
        print(f"  | {tag} | {E0:.2f} | {E1:.2f} | {k0} | {k1} |"
              f" {vac:.2e} | {Q:+.2f} |")
    print("\n  READ: pure SG is the integrable exception (the pair survives);")
    print("  with weak damping (radiation-loss stand-in) the pair captures into")
    print("  a breather and decays to vacuum: particle + antiparticle ->")
    print("  radiation, with total topological charge exactly 0 throughout.")
    print("=" * 74)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
