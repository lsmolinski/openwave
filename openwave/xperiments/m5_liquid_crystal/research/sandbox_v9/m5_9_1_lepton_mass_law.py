"""
M5.9.1 - Lepton mass law: does the V-on biaxial-hedgehog energy give E proportional to Lambda^(3/2)?

Tests the contributor (onspotgithub, 2026-06-18) proposal on the lepton-mass-spectrum
sector: the three charged leptons as ground states of three eigen-axes of the biaxial
Q-tensor, with a mass law E ~ Lambda^(3/2) from a gradient-vs-regularization balance and
the eigenvalue hierarchy Lambda_tau:Lambda_mu:Lambda_e ~ 229:35:1 taken as input
(Lambda ~ m^(2/3), so Lambda_mu/Lambda_e = (m_mu/m_e)^(2/3) = 34.97).

Machinery (M5.6.3b Faber-on-M): M = O D(s(r)) O^T with the Faber amplitude melt
s(r) = r/sqrt(r0^2+r^2), curvature energy E_curv = integral ||[M_mu, M_nu]||^2, and the
Faber potential V = (1-s^2)^3 / r0^4. Here we add an order-parameter eigenvalue amplitude
Lambda on the director eigenvalues: D_full = Lambda * diag(1, delta, 0).

We measure, per (Lambda, r0): E_curv, E_V, E_tot; and ask:
  A. baseline (Lambda=1): E*r0 ~ const  (reproduce the Faber 1/r0 mass-pinning)
  B. fixed r0: how does E scale with Lambda? (fit the exponent; the contributor needs 3/2)
  C. the lepton hierarchy: with Lambda_e:mu:tau = 1:34.97:229.5, do the E-ratios match the
     measured mass ratios m_mu/m_e = 206.77, m_tau/m_e = 3477.2?
  D. r0-relaxation: at fixed Lambda, is there an interior r0 that MINIMIZES E (the balance
     the contributor's dimensional argument assumes)?

Self-contained (numpy only). Headless.
    python m5_9_1_lepton_mass_law.py
"""
import numpy as np

DELTA = 0.3  # biaxiality eigenvalue spread (M5.6.3b default; the roadmap shows H is delta-flat)

# measured charged-lepton masses (MeV) and the contributor's Lambda = m^(2/3) hierarchy
M_E, M_MU, M_TAU = 0.510999, 105.658, 1776.86
LAM_E = 1.0
LAM_MU = (M_MU / M_E) ** (2.0 / 3.0)
LAM_TAU = (M_TAU / M_E) ** (2.0 / 3.0)


def matmul(A, B):
    return np.einsum("...ac,...cb->...ab", A, B)


def commf(A, B):
    return matmul(A, B) - matmul(B, A)


def frob2(M):
    return np.einsum("...ab,...ab->...", M, M)


def normalize(v, eps=1e-30):
    return v / (np.sqrt(np.einsum("...i,...i->...", v, v))[..., None] + eps)


def central(f, ax, h):
    out = np.zeros_like(f)
    sp_, sm, so = [slice(None)] * f.ndim, [slice(None)] * f.ndim, [slice(None)] * f.ndim
    sp_[ax], sm[ax], so[ax] = slice(2, None), slice(0, -2), slice(1, -1)
    out[tuple(so)] = (f[tuple(sp_)] - f[tuple(sm)]) / (2 * h)
    return out


def build_M(r0, lam=1.0, delta=DELTA, N=64, Lfac=9.0, melt=True):
    """M5.6.3b biaxial hedgehog with an order-parameter amplitude lam on the eigenvalues.
    D_full = lam*diag(1,delta,0); D_iso keeps the same trace; Faber melt s(r)."""
    D_full = lam * np.diag([1.0, delta, 0.0])
    D_iso = np.eye(3) * lam * (1.0 + delta) / 3.0
    L = Lfac * r0
    rhoc = 0.8 * r0
    xs = np.linspace(-L, L, N)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2)
    rho = np.sqrt(X**2 + Y**2)
    rhat = normalize(np.stack([X, Y, Z], -1) / np.sqrt(r**2 + r0**2)[..., None])
    azim = np.stack([-Y, X, np.zeros_like(X)], -1) / np.sqrt(rho**2 + 1e-12)[..., None]
    sm = np.clip(rho / rhoc, 0.0, 1.0)
    ePhi = azim * (sm * sm * (3.0 - 2.0 * sm))[..., None]
    ePhi = ePhi - np.einsum("...i,...i->...", ePhi, rhat)[..., None] * rhat
    eTheta = np.cross(ePhi, rhat)
    O = np.stack([rhat, eTheta, ePhi], axis=-1)
    s = (r / np.sqrt(r0**2 + r**2)) if melt else np.ones_like(r)
    Dfield = D_iso + s[..., None, None] * (D_full - D_iso)
    M = matmul(matmul(O, Dfield), np.swapaxes(O, -1, -2))
    return dict(M=M, s=s, h=h, r=r, rho=rho, L=L, N=N)


def energies(bg, r0):
    M, h, s = bg["M"], bg["h"], bg["s"]
    Mmu = [central(M, ax, h) for ax in range(3)]
    interior = tuple(slice(2, -2) for _ in range(3))
    curv = np.zeros(M.shape[:3])
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        curv += frob2(commf(Mmu[i], Mmu[j]))
    Vpot = (1.0 - s**2) ** 3 / r0**4
    E_curv = float(np.sum(curv[interior]) * h**3)
    E_V = float(np.sum(Vpot[interior]) * h**3)
    return E_curv, E_V


def fit_pow(x, y):
    """log-log slope p in y ~ x^p, plus R^2."""
    lx, ly = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    p, b = np.polyfit(lx, ly, 1)
    yhat = p * lx + b
    ss_res = np.sum((ly - yhat) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(p), float(r2)


def main():
    print("=" * 78)
    print("M5.9.1 - lepton mass law: E ~ Lambda^(3/2)? (V-on biaxial hedgehog)")
    print(f"  delta={DELTA}  hierarchy (Lambda=m^(2/3)): Lambda_mu/e={LAM_MU:.2f}  Lambda_tau/e={LAM_TAU:.2f}")
    print("=" * 78)

    # --- A. baseline: E*r0 ~ const at Lambda=1 (Faber 1/r0 mass-pinning) -------------
    print("\n[A] baseline (Lambda=1): E*r0 should be ~const (E ~ 1/r0)")
    print(f"    {'r0':>5} {'E_curv':>11} {'E_V':>11} {'E_tot':>11} {'E_tot*r0':>11}")
    prods = []
    for r0 in (0.8, 1.2, 1.6):
        ec, ev = energies(build_M(r0, lam=1.0), r0)
        et = ec + ev
        prods.append(et * r0)
        print(f"    {r0:>5.1f} {ec:>11.4f} {ev:>11.4f} {et:>11.4f} {et*r0:>11.4f}")
    cvA = np.std(prods) / np.mean(prods)
    print(f"    -> E*r0 CV = {100*cvA:.1f}%  (const if <15%): {cvA < 0.15}")

    # --- B. fixed r0: E vs Lambda exponent -------------------------------------------
    print("\n[B] amplitude sweep at fixed r0=1.0: E_curv, E_V, E_tot vs Lambda")
    lams = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0])
    ecs, evs, ets = [], [], []
    print(f"    {'Lambda':>7} {'E_curv':>12} {'E_V':>12} {'E_tot':>12}")
    for lam in lams:
        ec, ev = energies(build_M(1.0, lam=lam), 1.0)
        ecs.append(ec); evs.append(ev); ets.append(ec + ev)
        print(f"    {lam:>7.2f} {ec:>12.5f} {ev:>12.5f} {ec+ev:>12.5f}")
    p_c, r2c = fit_pow(lams, ecs)
    p_v, r2v = fit_pow(lams, np.maximum(evs, 1e-30))
    # E_tot dominated by whichever term; fit it directly
    p_t, r2t = fit_pow(lams, ets)
    print(f"    -> E_curv ~ Lambda^{p_c:.2f} (R2={r2c:.3f}) ; E_V ~ Lambda^{p_v:.2f} ; E_tot ~ Lambda^{p_t:.2f}")
    print(f"       contributor's mass law needs E ~ Lambda^1.5 ; pure amplitude scaling predicts curv ~ Lambda^4")

    # --- C. lepton hierarchy: do E-ratios match the mass ratios? ----------------------
    print("\n[C] lepton hierarchy (Lambda_e:mu:tau = 1:%.1f:%.1f), fixed r0=1.0" % (LAM_MU, LAM_TAU))
    Ee = sum(energies(build_M(1.0, lam=LAM_E), 1.0))
    Emu = sum(energies(build_M(1.0, lam=LAM_MU), 1.0))
    Etau = sum(energies(build_M(1.0, lam=LAM_TAU), 1.0))
    print(f"    E_e={Ee:.4f}  E_mu={Emu:.4f}  E_tau={Etau:.4f}")
    print(f"    E_mu/E_e = {Emu/Ee:.1f}   (measured m_mu/m_e = {M_MU/M_E:.1f})")
    print(f"    E_tau/E_e = {Etau/Ee:.1f}   (measured m_tau/m_e = {M_TAU/M_E:.1f})")
    print(f"    If E ~ Lambda^1.5 the ratios would be {LAM_MU**1.5:.1f} and {LAM_TAU**1.5:.1f} (= the masses, by construction)")

    # --- D. r0-relaxation: interior minimum at fixed Lambda? --------------------------
    print("\n[D] r0-relaxation at fixed Lambda=1: is there an INTERIOR E-minimum in r0?")
    r0s = np.array([0.5, 0.7, 1.0, 1.4, 2.0, 2.8])
    ets_r = [sum(energies(build_M(r0, lam=1.0), r0)) for r0 in r0s]
    for r0, e in zip(r0s, ets_r):
        print(f"    r0={r0:>4.1f}  E_tot={e:>10.4f}")
    imin = int(np.argmin(ets_r))
    interior_min = 0 < imin < len(r0s) - 1
    print(f"    -> argmin at r0={r0s[imin]:.1f} ; interior minimum: {interior_min}")
    print("       (monotonic E(r0) => r0 is a free scale-setting modulus, not a balance;")
    print("        the contributor's E ~ Lambda^(3/2) assumes a regularization that grows with r0)")

    print("\n" + "=" * 78)
    print("VERDICT (static structural):")
    mass_law = abs(p_c - 1.5) < 0.3
    print(f"  - E_curv exponent in Lambda: {p_c:.2f} (mass law wants 1.5) -> {'MATCH' if mass_law else 'NO MATCH'}")
    print(f"  - r0 interior minimum: {interior_min} (mass law's balance needs True)")
    print(f"  - baseline E~1/r0 reproduced: {cvA < 0.15}")
    print("=" * 78)
    np.savez("m5_9_1_results.npz", lams=lams, ecs=ecs, evs=evs, ets=ets,
             p_curv=p_c, p_V=p_v, p_tot=p_t, cvA=cvA, r0s=r0s, ets_r=ets_r,
             Ee=Ee, Emu=Emu, Etau=Etau, LAM_MU=LAM_MU, LAM_TAU=LAM_TAU)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
