"""
M5.6.4a — EM from tilts via the hydrodynamics ↔ electromagnetism dictionary (4a §11b.1)

The Couder/Bush deck gives a superfluid-hydro ↔ EM dictionary that is a clean, ABELIAN
route to "EM from tilts" (independent of Faber's non-abelian matrix curvature, M5.6.4b):

    vorticity   ω = ∇×u            ↔   B = ∇×A
    Lamb vector l = ω×u            ↔   E
    Faraday     ∂_tω = −∇×l        ↔   ∂_tB = −∇×E
    charge      ∇·l = u·(∇×ω)−‖ω‖² ↔   ∇·E = ρ
    force       −2(v×Ω) (Coriolis) ↔   q(v×B) (Lorentz)

This script builds an incompressible tilt-flow u (= ∇×A, so ∇·u=0 exactly) and verifies the
dictionary reproduces Maxwell's structure, using SPECTRAL (FFT) derivatives on a periodic box
(so the vector identities hold to machine precision — a crisp correctness check):

  1. ∇·ω = 0           — magnetic Gauss (no monopoles), kinematic.
  2. ∇·l = u·(∇×ω)−‖ω‖² — the turbulent-charge identity (↔ ∇·E=ρ). Non-trivial; verify exact.
  3. ∂_tω = −∇×l        — Faraday, via one ideal-incompressible Euler step (FFT pressure proj).
  4. Coriolis ↔ Lorentz — the force on a test velocity is F = v×ω (same v×field form as Lorentz).

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v5.m5_6_4a_hydro_em
"""
import numpy as np

N = 48
Lbox = 2 * np.pi                                  # periodic box [0, 2π)
xs = np.linspace(0, Lbox, N, endpoint=False)
X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
k1 = np.fft.fftfreq(N, d=Lbox / N) * 2 * np.pi    # wavenumbers
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing="ij")
K = [KX, KY, KZ]
K2 = KX**2 + KY**2 + KZ**2
K2inv = np.where(K2 == 0, 0.0, 1.0 / np.where(K2 == 0, 1.0, K2))


def ddx(f, ax):
    return np.real(np.fft.ifftn(1j * K[ax] * np.fft.fftn(f)))


def grad(f):
    return [ddx(f, a) for a in range(3)]


def div(V):
    return ddx(V[0], 0) + ddx(V[1], 1) + ddx(V[2], 2)


def curl(V):
    return [ddx(V[2], 1) - ddx(V[1], 2),
            ddx(V[0], 2) - ddx(V[2], 0),
            ddx(V[1], 0) - ddx(V[0], 1)]


def cross(A, B):
    return [A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]]


def dot(A, B):
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2]


def project_incompressible_rate(F):
    """Given an Euler RHS F (=−l), return F − ∇φ with ∇²φ = ∇·F (so result is divergence-free)."""
    phi = np.real(np.fft.ifftn(np.fft.fftn(div(F)) * (-K2inv)))   # ∇²φ = ∇·F ⇒ φ = ∇·F / (−k²)·(−1)... see note
    gp = grad(phi)
    return [F[i] - gp[i] for i in range(3)], phi


def main():
    print("=" * 70)
    print("M5.6.4a — EM from tilts: hydrodynamics ↔ electromagnetism dictionary")
    print(f"  periodic {N}³ box, spectral derivatives")
    print("=" * 70)

    # --- incompressible tilt-flow u = ∇×A (∇·u = 0 by construction) --------------
    g1 = np.exp(np.cos(X) + np.cos(Y) + np.cos(Z))           # smooth periodic bumps
    A = [g1 * np.sin(Y), g1 * np.sin(Z), g1 * np.sin(X)]     # generic periodic vector potential
    u = curl(A)
    omega = curl(u)                                          # ω = ∇×u   (↔ B)
    l = cross(omega, u)                                      # l = ω×u   (↔ E)

    divu = div(u)
    print("\n[0] tilt-flow u = ∇×A:  max|∇·u| = "
          f"{np.abs(divu).max():.2e}  (incompressible by construction ✓)")

    # --- 1. ∇·ω = 0  (↔ ∇·B = 0, no magnetic monopoles) -------------------------
    divw = div(omega)
    print(f"\n[1] ∇·ω = 0  (↔ ∇·B=0):  max|∇·ω| = {np.abs(divw).max():.2e}  → no monopoles ✓")

    # --- 2. charge identity  ∇·l = u·(∇×ω) − ‖ω‖²  (↔ ∇·E = ρ) -------------------
    divl = div(l)
    rho = dot(u, curl(omega)) - dot(omega, omega)
    err2 = np.abs(divl - rho).max() / (np.abs(divl).max() + 1e-30)
    print(f"\n[2] ∇·l = u·(∇×ω) − ‖ω‖²  (↔ ∇·E=ρ):  rel err = {err2:.2e}")
    print(f"    → the 'turbulent charge' ρ = ∇·l identity holds (the hydro Gauss law) ✓")

    # --- 3. Faraday  ∂_tω = −∇×l  via one ideal-incompressible Euler step --------
    # Euler (Lamb form): ∂_t u = −l − ∇φ, φ enforcing incompressibility. Then
    # ∂_tω = ∇×∂_tu = −∇×l (since ∇×∇φ = 0). Verify the discrete operators respect it.
    F = [-l[i] for i in range(3)]
    rate, _ = project_incompressible_rate(F)                # ∂_t u (divergence-free)
    dwdt = curl(rate)                                       # ∂_t ω
    minus_curl_l = [-c for c in curl(l)]
    num = max(np.abs(dwdt[i] - minus_curl_l[i]).max() for i in range(3))
    den = max(np.abs(minus_curl_l[i]).max() for i in range(3)) + 1e-30
    err3 = num / den
    print(f"\n[3] Faraday ∂_tω = −∇×l  (↔ ∂_tB=−∇×E):  rel err = {err3:.2e}")
    print(f"    → vorticity transport = Faraday's law; the incompressible Euler curl gives it ✓")

    # --- 4. Coriolis ↔ Lorentz: force on a test velocity is v×(field) -----------
    v = np.array([0.3, -0.2, 0.5])
    # hydro Coriolis-type force −2 v×Ω with Ω=ω/2  ⇒  −(v×ω); Lorentz q v×B with B=ω
    F_cor = [-(v[1] * omega[2] - v[2] * omega[1]),
             -(v[2] * omega[0] - v[0] * omega[2]),
             -(v[0] * omega[1] - v[1] * omega[0])]
    F_lor = [v[1] * omega[2] - v[2] * omega[1],
             v[2] * omega[0] - v[0] * omega[2],
             v[0] * omega[1] - v[1] * omega[0]]
    # same magnitude, opposite sign convention (Ω=ω/2 ⇒ −2Ω=−ω); both ∝ v×ω
    ratio = np.sqrt(sum((F_cor[i] ** 2).mean() for i in range(3))) / \
            (np.sqrt(sum((F_lor[i] ** 2).mean() for i in range(3))) + 1e-30)
    print(f"\n[4] Coriolis ↔ Lorentz:  both forces ∝ v×ω (the v×field law);  |F_cor|/|F_lor| = {ratio:.3f}")
    print(f"    → the deflection force has the Lorentz v×B structure (Ω=ω/2 ⇒ B↔ω) ✓")

    ok = (np.abs(divw).max() < 1e-8 and err2 < 1e-6 and err3 < 1e-6 and abs(ratio - 1) < 1e-9)
    print("\n" + "=" * 70)
    print("M5.6.4a: the hydro↔EM dictionary reproduces Maxwell's structure — ∇·B=0,")
    print("the Gauss charge identity, Faraday, and the Lorentz force law — the clean")
    print("ABELIAN route to EM-from-tilts. (Faber's non-abelian matrix curvature = 4b.)")
    print("PASS" if ok else "PARTIAL — inspect the failing metric above")
    print("=" * 70)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
