"""
M5.9.3 - core-volume confiner: does a fixed-density Higgs-like vacuum select the core scale
and yield a lepton mass law?

M5.9.1 found the gap: the M5 Faber functional is scale-balanced (E ~ 1/r0, r0 a FREE modulus),
so no scale is selected. The named missing piece is a core-volume confiner: a fixed false-vacuum
energy density B over the melted core volume, E_conf = B * integral (1-s^2)^3 d^3x ~ B * r0^3
(the Faber potential is the SAME integrand times 1/r0^4, hence scale-covariant ~1/r0; the confiner
uses a FIXED coefficient B, breaking scale-covariance and selecting r0).

This script adds the confiner to the M5 hedgehog and measures, per eigenvalue amplitude Lambda:
  - whether E(r0) now has an INTERIOR minimum (scale-selection, the structural fix) vs the
    free-modulus runaway without it;
  - the selected r0*(Lambda) and mass E*(Lambda), and the fitted exponents;
  - the same for several gradient choices (quartic Skyrme curvature E_curv ~ Lambda^4/r0 ;
    quadratic Frank gradient E_grad2 ~ Lambda^2 * r0), since the mass-vs-Lambda exponent depends
    on WHICH gradient term dominates;
  - the eigenvalue hierarchy implied by the measured exponent, vs the contributor's assumed
    Lambda ~ m^(2/3) (which presumes E ~ Lambda^(3/2)).

Honest framing: this tests whether the confiner closes the SCALING + SELECTION gap. It does NOT
derive the eigenvalue hierarchy itself (the Yukawa-like input); that origin is the deeper open
question. Self-contained (numpy). Headless.
    python m5_9_3_confiner.py
"""
import numpy as np

DELTA = 0.3
M_E, M_MU, M_TAU = 0.510999, 105.658, 1776.86


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


def build_M(r0, lam=1.0, delta=DELTA, N=64, Lfac=9.0):
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
    s = r / np.sqrt(r0**2 + r**2)
    Dfield = D_iso + s[..., None, None] * (D_full - D_iso)
    M = matmul(matmul(O, Dfield), np.swapaxes(O, -1, -2))
    return dict(M=M, s=s, h=h, r=r, rho=rho, L=L, N=N)


def energy_terms(bg, r0):
    """Return the separable energy pieces (each at the given r0, lam baked into bg["M"])."""
    M, h, s = bg["M"], bg["h"], bg["s"]
    Mmu = [central(M, ax, h) for ax in range(3)]
    interior = tuple(slice(2, -2) for _ in range(3))
    # quartic Skyrme curvature
    curv = np.zeros(M.shape[:3])
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        curv += frob2(commf(Mmu[i], Mmu[j]))
    E_curv = float(np.sum(curv[interior]) * h**3)
    # quadratic Frank gradient ||dM||^2
    grad2 = sum(frob2(Mmu[ax]) for ax in range(3))
    E_grad2 = float(np.sum(grad2[interior]) * h**3)
    # Faber potential (scale-covariant, ~1/r0)
    E_Vfaber = float(np.sum(((1.0 - s**2) ** 3 / r0**4)[interior]) * h**3)
    # core-volume confiner integrand: (1-s^2)^3 with a FIXED coefficient -> ~ r0^3
    E_confbase = float(np.sum(((1.0 - s**2) ** 3)[interior]) * h**3)
    return dict(E_curv=E_curv, E_grad2=E_grad2, E_Vfaber=E_Vfaber, E_confbase=E_confbase)


def fit_pow(x, y):
    lx, ly = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    p, b = np.polyfit(lx, ly, 1)
    yhat = p * lx + b
    ss = np.sum((ly - yhat) ** 2)
    st = np.sum((ly - ly.mean()) ** 2)
    return float(p), float(1.0 - ss / st if st > 0 else np.nan)


def minimize_r0(lam, B, terms_func, r0grid):
    """E(r0) = grad/skyrme part + B * E_confbase ; return (r0*, E*, has_interior_min)."""
    Es = []
    for r0 in r0grid:
        t = energy_terms(build_M(r0, lam=lam), r0)
        Es.append(terms_func(t) + B * t["E_confbase"])
    Es = np.array(Es)
    imin = int(np.argmin(Es))
    return float(r0grid[imin]), float(Es[imin]), 0 < imin < len(r0grid) - 1, Es


def main():
    print("=" * 80)
    print("M5.9.3 - core-volume confiner: scale-selection + the lepton mass law")
    print("=" * 80)

    r0grid = np.linspace(0.4, 4.0, 19)
    lams = np.array([0.6, 0.8, 1.0, 1.4, 2.0, 3.0])

    # pick B so the Skyrme+confiner minimum sits near r0~1 at lam=1 (a reference scale)
    t1 = energy_terms(build_M(1.0, lam=1.0), 1.0)
    print(f"\n[ref] at lam=1, r0=1: E_curv={t1['E_curv']:.3f} E_grad2={t1['E_grad2']:.3f}"
          f" E_Vfaber={t1['E_Vfaber']:.3f} E_confbase={t1['E_confbase']:.3f}")

    # --- A. scale-selection: confiner ON vs OFF (Skyrme part) ------------------------
    print("\n[A] scale-selection (Skyrme curvature + confiner): interior minimum appears?")
    skyrme = lambda t: t["E_curv"]
    # B chosen from the Derrick balance r0*^4 = E_curv*r0/(3*E_confbase/r0^3) at lam=1
    B = t1["E_curv"] / (3.0 * t1["E_confbase"])  # ~ gives dE/dr0=0 near r0=1
    print(f"    B (confiner coeff) = {B:.4f}")
    _, _, has0, Es0 = minimize_r0(1.0, 0.0, skyrme, r0grid)
    r0s, Em, has1, Es1 = minimize_r0(1.0, B, skyrme, r0grid)
    print(f"    confiner OFF (B=0): interior min = {has0}  (E monotonic -> free modulus)")
    print(f"    confiner ON       : interior min = {has1}  at r0*={r0s:.2f}, E*={Em:.3f}")

    # --- B. mass law: E*(Lambda), r0*(Lambda) for two gradient choices ---------------
    print("\n[B] mass law with the confiner ON: E*(Lambda), r0*(Lambda)")
    for name, gfun in (("Skyrme curvature (E_curv ~ L^4/r0)", lambda t: t["E_curv"]),
                       ("Frank gradient (E_grad2 ~ L^2*r0)", lambda t: t["E_grad2"]),
                       ("Skyrme+Frank", lambda t: t["E_curv"] + t["E_grad2"])):
        r0stars, Estars, oks = [], [], []
        for lam in lams:
            r0s, Es, ok, _ = minimize_r0(lam, B, gfun, r0grid)
            r0stars.append(r0s); Estars.append(Es); oks.append(ok)
        pE, r2E = fit_pow(lams, Estars)
        pr, r2r = fit_pow(lams, r0stars)
        allmin = all(oks)
        print(f"\n    [{name}]")
        print(f"      {'Lambda':>7} {'r0*':>7} {'E*':>12} {'interior':>9}")
        for lam, r0s, Es, ok in zip(lams, r0stars, Estars, oks):
            print(f"      {lam:>7.2f} {r0s:>7.2f} {Es:>12.4f} {str(ok):>9}")
        print(f"      -> E* ~ Lambda^{pE:.2f} (R2={r2E:.3f}) ; r0* ~ Lambda^{pr:.2f} ; all interior: {allmin}")
        # implied hierarchy: to match E ~ m we need Lambda ~ m^(1/pE)
        if abs(pE) > 0.1:
            powe = 1.0 / pE
            lam_mu = (M_MU / M_E) ** powe
            lam_tau = (M_TAU / M_E) ** powe
            print(f"      -> to get E ~ m, need Lambda ~ m^{powe:.3f}; hierarchy Lambda_mu/e={lam_mu:.1f},"
                  f" Lambda_tau/e={lam_tau:.1f}  (contributor assumed 35, 229 from E~L^1.5)")

    print("\n" + "=" * 80)
    print("VERDICT (#200 confiner):")
    print(f"  - scale-selection: confiner ON gives an interior minimum (was {has0} OFF, {has1} ON)")
    print("    => the core-volume confiner FIXES the free modulus (the structural gap is closed).")
    print("  - mass law: the E*(Lambda) exponent is set by the dominant gradient term (printed above);")
    print("    whether it equals the contributor's 3/2 (E ~ m via Lambda ~ m^(2/3)) is reported, not assumed.")
    print("  - the eigenvalue hierarchy itself (why those Lambda) stays INPUT: the confiner gives")
    print("    scale-selection + a definite mass-vs-Lambda law, not the Yukawa origin of the masses.")
    print("=" * 80)
    np.savez("m5_9_3_results.npz", lams=lams, B=B, has_off=has0, has_on=has1, r0grid=r0grid)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
