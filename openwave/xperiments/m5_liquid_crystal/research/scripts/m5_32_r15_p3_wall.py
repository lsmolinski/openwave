"""M5.32 R15-P-iii: the reduced planar functional of L_P on the true degenerate vacuum (prediction (iii)).

Equations first.  Diagonal planar profiles M(z) = diag(g, 1, m2(z), m3(z)), g = 8, delta = 0.3,
N = M eta with spectrum (-g, 1, m2, m3); the exterior at rest is the degenerate vacuum (delta, delta).
    V4^dd(m2, m3) = W1 sum_{p=1..4} (tr N^p - C_p)^2,   C_p = (-g)^p + 1 + 2 delta^p,   W1 = 0.000724023879
    SPLIT          = (m2 - m3)^2 =: s^2
    K_P^23 static  = (c / 2)(m2'^2 + m3'^2)      [P23 = diag(0,0,1,1) on a diagonal profile, Om_z = diag(0,0,m2',m3') eta]
    K_P^23 inertia = c kappa_P s^2, kappa_P = 1     [a0 = G1 M - M G1 = s (e23 + e32) up to sign, Om_0 = a0 eta]
    I1 inertia     = 8 s^2 s'^2                     [the R14-D2 flank inertia; every certified curvature term is
                                                     static-zero on a planar profile]
Rotating-frame functional at fixed omega (E-orientation, the R14-D2 form):
    F[m2, m3] = int dz { (c/2)(m2'^2 + m3'^2) + V4^dd + mu s^2 - omega^2 [c s^2 + 8 s^2 s'^2] }
    V_eff(m2, m3) = V4^dd + (mu - omega^2 c) s^2      (the uniform part)
THEOREM (checked numerically here): V4^dd >= 0 with equality only at the vacuum spectrum, so for
omega^2 < mu / c every uniform state has V_eff > 0 = V_eff(exterior): no Maxwell crossing exists
below omega_c^2 = mu / c, and at omega_c the exterior loses stability along the split (its split
stiffness is exactly 4 mu, R15-M-a).  Pre-registered expectation: CONTINUOUS_ONSET at omega_c^2 = mu / (c kappa_P);
FIRST_ORDER_CROSSING would need a second minimum of V_eff below the exterior for omega < omega_c.
What DOES exist above onset: the Ising wall between the two split orientations +-s* (m2 <-> m3, a (2,3)
rotation by pi/2): its tension and 10-90 width are measured at omega^2 = 2 mu / c on the reduced line, and
the width of a static split perturbation decaying into the exterior is sqrt(c / (2 mu)) per component
(linearized along the split at fixed trace: m2' = s'/2, m3' = -s'/2, so (c/2)(m2'^2 + m3'^2) = (c/4) s'^2
against mu s^2: decay length (1/2) sqrt(c/mu), the author's (c_P / mu)^(1/2) scaling), read on the R15-M
hedgehog split profiles.  The same c/4 makes the fixed-omega s'^2 coefficient (c/4 - 8 omega^2 s^2): the
reduced functional is unbounded below wherever |s| > sqrt(c / (32 omega^2)) (the audit's correction of the
32-vs-16 factor in the first run's log line).  AUDIT NOTE (2026-09-05): with the off-diagonal (2,3) entry
free, the "Ising wall" between (a, b) and (b, a) is a Goldstone twist by the (2,3) rotation (a symmetry of
L_P), tension -> 0 on the infinite line: it is a saddle of the diagonal sector, not an object of the model.
Reduced-density certification: the profile placed on a 4 x 4 x n_z lattice slab and recomputed by the
registry (C15.lp_parts with a0 = G1 M - M G1) must agree with the reduced sums.

usage: python3 m5_32_r15_p3_wall.py
"""
import sys
ARGS = list(sys.argv[1:])
import os, json, glob
import numpy as np
from scipy.optimize import minimize
import m5_32_r15_common as C15

INS4, C13 = C15.INS4, C15.C13
log = C15.log
G, DELTA, W1 = C15.G, C15.DELTA, C15.W1
CP = [(-G) ** p + 1.0 + 2.0 * DELTA ** p for p in range(1, 5)]


def v4(m2, m3):
    lam = np.stack([np.full_like(m2, -G), np.full_like(m2, 1.0), m2, m3], axis=-1)
    return W1 * sum((np.sum(lam ** p, axis=-1) - CP[p - 1]) ** 2 for p in range(1, 5))


def v4_grad(m2, m3):
    lam = np.stack([np.full_like(m2, -G), np.full_like(m2, 1.0), m2, m3], axis=-1)
    g2 = np.zeros_like(m2); g3 = np.zeros_like(m3)
    for p in range(1, 5):
        tp = np.sum(lam ** p, axis=-1) - CP[p - 1]
        g2 += 2.0 * tp * p * m2 ** (p - 1)
        g3 += 2.0 * tp * p * m3 ** (p - 1)
    return W1 * g2, W1 * g3


def veff(m2, m3, mu, c, om):
    return v4(m2, m3) + (mu - om * om * c) * (m2 - m3) ** 2


def F_reduced(m2, m3, h, mu, c, om, grad=False):
    """finite-volume form: site terms V_eff h; link terms h [(c/2)(dm2^2 + dm3^2)/h^2 - 8 om^2 sbar^2 (ds/h)^2]."""
    s = m2 - m3
    site = veff(m2, m3, mu, c, om)
    dm2, dm3, ds = np.diff(m2), np.diff(m3), np.diff(s)
    sbar = 0.5 * (s[:-1] + s[1:])
    link = (c / 2.0) * (dm2 ** 2 + dm3 ** 2) / h - 8.0 * om * om * sbar ** 2 * ds ** 2 / h
    F = float(np.sum(site) * h + np.sum(link))
    if not grad:
        dens = np.zeros_like(m2); dens += site * h
        dens[:-1] += 0.5 * link; dens[1:] += 0.5 * link
        return F, dens
    g2v, g3v = v4_grad(m2, m3)
    gs_site = 2.0 * (mu - om * om * c) * s
    g2 = h * (g2v + gs_site); g3 = h * (g3v - gs_site)
    # link derivatives
    dl_ddm = c * dm2 / h; dl_ddm3 = c * dm3 / h
    g2[:-1] -= dl_ddm; g2[1:] += dl_ddm
    g3[:-1] -= dl_ddm3; g3[1:] += dl_ddm3
    dl_dds = -16.0 * om * om * sbar ** 2 * ds / h
    dl_dsbar = -16.0 * om * om * sbar * ds ** 2 / h
    gs = np.zeros_like(s)
    gs[:-1] += -dl_dds + 0.5 * dl_dsbar
    gs[1:] += dl_dds + 0.5 * dl_dsbar
    return F, g2 + gs, g3 - gs


def gradient_gate(mu=1e-2, c=1.0, om=0.05):
    rng = np.random.default_rng(3)
    n, h = 60, 0.5
    m2 = DELTA + 0.1 * rng.normal(size=n); m3 = DELTA + 0.1 * rng.normal(size=n)
    F, g2, g3 = F_reduced(m2, m3, h, mu, c, om, grad=True)
    d2, d3 = rng.normal(size=n), rng.normal(size=n)
    an = float(np.sum(g2 * d2 + g3 * d3))
    e = 1e-6
    fd = (F_reduced(m2 + e * d2, m3 + e * d3, h, mu, c, om)[0] - F_reduced(m2 - e * d2, m3 - e * d3, h, mu, c, om)[0]) / (2 * e)
    return {"analytic": an, "fd": fd, "rel": abs(an - fd) / abs(fd)}


def theorem_check(mu, c):
    """min over a wide plane grid of V_eff at omega^2 = 0.999 mu/c: must be >= 0 with the only zero at (delta, delta)."""
    grid = np.linspace(-2.0, 2.0, 801)
    M2, M3 = np.meshgrid(grid, grid, indexing="ij")
    om = np.sqrt(0.999 * mu / c) if c > 0 else 0.0
    V = veff(M2, M3, mu, c, om)
    i, j = np.unravel_index(np.argmin(V), V.shape)
    Vs = veff(M2, M3, mu, c, 0.0)
    zero_set = np.sum(Vs < 1e-9)
    return {"omega2_over_omega_c2": 0.999, "min_Veff": float(V[i, j]), "argmin": [float(grid[i]), float(grid[j])],
            "n_grid_points_with_V4dd_plus_split_below_1e-9": int(zero_set), "second_lowest_local_min_Veff": None}


def onset_scan(mu, c):
    """sigma*(omega): the split amplitude of the lowest V_eff minimum vs omega across omega_c."""
    if c == 0:
        return {"omega_c2": None, "note": "c = 0: no inertia, no onset"}
    omc2 = mu / c
    rows = []
    grid = np.linspace(-1.5, 1.5, 301)
    M2, M3 = np.meshgrid(grid, grid, indexing="ij")
    for f in (0.5, 0.9, 0.99, 1.0, 1.01, 1.1, 1.5, 2.0, 3.0):
        om = np.sqrt(f * omc2)
        V = veff(M2, M3, mu, c, om)
        i, j = np.unravel_index(np.argmin(V), V.shape)
        r = minimize(lambda x: float(veff(np.array(x[0]), np.array(x[1]), mu, c, om)), [grid[i], grid[j]], method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-14})
        rows.append({"omega2_over_omega_c2": f, "omega": float(om), "min": r.x.tolist(), "V_eff_min": float(r.fun), "split": float(abs(r.x[0] - r.x[1])), "sum": float(r.x[0] + r.x[1])})
    return {"omega_c2": omc2, "rows": rows}


def ising_wall(mu, c, f=2.0, L=None, n=None):
    """the wall between +-s* at omega^2 = f mu/c on the reduced line (both ends pinned)."""
    omc2 = mu / c
    om = np.sqrt(f * omc2)
    sc = onset_scan.__wrapped__ if hasattr(onset_scan, "__wrapped__") else None
    grid = np.linspace(-1.5, 1.5, 301)
    M2, M3 = np.meshgrid(grid, grid, indexing="ij")
    V = veff(M2, M3, mu, c, om)
    i, j = np.unravel_index(np.argmin(V), V.shape)
    r = minimize(lambda x: float(veff(np.array(x[0]), np.array(x[1]), mu, c, om)), [grid[i], grid[j]], method="Nelder-Mead", options={"xatol": 1e-12, "fatol": 1e-16})
    a, b = sorted(r.x, reverse=True)
    plus, minus = np.array([a, b]), np.array([b, a])
    wguess = np.sqrt(c / mu)
    L = L or 40.0 * wguess
    n = n or max(400, int(L / (wguess / 20)))
    h = L / n
    z = (np.arange(n) + 0.5) * h - L / 2
    prof = 0.5 * (1 + np.tanh(z / wguess))
    m2 = minus[0] + (plus[0] - minus[0]) * prof
    m3 = minus[1] + (plus[1] - minus[1]) * prof
    npin = 4
    x0 = np.concatenate([m2[npin:-npin], m3[npin:-npin]])

    def unpack(x):
        A = m2.copy(); B = m3.copy()
        A[npin:-npin] = x[:n - 2 * npin]; B[npin:-npin] = x[n - 2 * npin:]
        return A, B

    def fun(x):
        A, B = unpack(x)
        F, g2, g3 = F_reduced(A, B, h, mu, c, om, grad=True)
        return F, np.concatenate([g2[npin:-npin], g3[npin:-npin]])
    res = minimize(fun, x0, jac=True, method="L-BFGS-B", options={"maxiter": 50000, "maxfun": 200000, "ftol": 1e-16, "gtol": 1e-11})
    A, B = unpack(res.x)
    Fw, dens = F_reduced(A, B, h, mu, c, om)
    Funi = float(veff(np.array(a), np.array(b), mu, c, om)) * L
    s = A - B
    lo, hi = -0.8 * (a - b), 0.8 * (a - b)
    idx = np.where((s - lo) * (s - hi) < 0)[0]
    width = float(h * (idx.max() - idx.min())) if len(idx) > 1 else None
    return {"omega2_over_omega_c2": f, "omega": float(om), "plus": plus.tolist(), "minus": minus.tolist(), "s_star": float(a - b),
            "sigma": Fw - Funi, "width_10_90": width, "sqrt_c_over_mu": float(wguess), "width_over_sqrt_c_over_mu": (width / wguess) if width else None,
            "L": L, "n": n, "h": h, "success": bool(res.success), "nit": int(res.nit), "z": z.tolist(), "profile_m2": A.tolist(), "profile_m3": B.tolist()}


def slab_check(mu=1e-2, c=1.0, om=0.07):
    """the reduced densities against the registry on a 4 x 4 x nz slab."""
    nz, h = 64, 0.75
    z = (np.arange(nz) + 0.5) * h - nz * h / 2
    m2 = DELTA + 0.25 * np.exp(-z * z / 30.0) + 0.05 * np.exp(-(z - 6.0) ** 2 / 12.0)   # compact bumps: flat at both ends
    m3 = DELTA - 0.20 * np.exp(-z * z / 40.0)
    cfg = C15.cfg_dd(4, 4 * h, mu=mu, cP=c)
    cfg["n"] = 4; cfg["h"] = h
    M = np.zeros((4, 4, nz, 4, 4))
    M[..., 0, 0] = G; M[..., 1, 1] = 1.0
    M[..., 2, 2] = m2[None, None, :]; M[..., 3, 3] = m3[None, None, :]
    # interior cells only (the one-sided edges along z are excluded by summing cells 1..nz-2 in the reduced form)
    a0 = C13.a0_G1(M)
    pp = C15.lp_parts(M, cfg, a0)
    area = 16 * h * h
    s = m2 - m3
    # the sym stencil = the average of the forward and backward one-sided differences per cell
    def fb(f):
        fwd = np.zeros_like(f); bwd = np.zeros_like(f)
        fwd[:-1] = (f[1:] - f[:-1]) / h; bwd[1:] = (f[1:] - f[:-1]) / h
        return fwd, bwd
    f2, b2 = fb(m2); f3, b3 = fb(m3); fs, bs = fb(s)
    red = {"V4": float(np.sum(v4(m2, m3)) * h * area), "split": float(mu * np.sum(s * s) * h * area),
           "KP_static": float(c * 0.5 * 0.5 * np.sum(f2 ** 2 + b2 ** 2 + f3 ** 2 + b3 ** 2) * h * area),
           "KP_inertia_kappaP1": float(np.sum(s * s) * h * area), "I1_inertia": float(8.0 * 0.5 * np.sum(s * s * (fs ** 2 + bs ** 2)) * h * area)}
    lat = {"V4": pp["V4"], "split": pp["split"], "KP_static": pp["KP"], "KP_inertia_kappaP1": pp["kin_KP_raw"], "I1_inertia": pp["kin_I1"], "E_u": pp["E_u"]}
    rel = {k: abs(red[k] - lat[k]) / max(abs(lat[k]), 1e-300) for k in red}
    return {"reduced": red, "lattice": lat, "rel": rel, "note": "the sym stencil = the fwd/bwd average per cell; the reduced sums use the same form"}


def decay_length_from_m_fields():
    """the split profile along +x of the R15-M relaxed hedgehogs: exponential tail length vs (1/2) sqrt(c/mu)."""
    out = []
    for f in sorted(glob.glob(os.path.join(C15.CK, "m_hedgehog", "relax_*.json"))):
        rec = json.load(open(f))
        if rec["mu"] <= 0:
            continue
        e = rec["end"]
        r = np.array(e["axis_r"]); sp_ = np.array(e["axis_split"])
        L = rec["L"]
        sel = (r > 4.0) & (r < 0.5 * L - 2.5) & (sp_ > 1e-6)
        if np.sum(sel) < 4:
            out.append({"run": os.path.basename(f), "fit": None})
            continue
        slope, icpt = np.polyfit(r[sel], np.log(sp_[sel]), 1)
        out.append({"run": os.path.basename(f), "mu": rec["mu"], "cP": rec["cP"], "n": rec["n"], "L": L,
                    "tail_length": float(-1.0 / slope) if slope < 0 else None, "half_sqrt_c_over_mu": float(0.5 * np.sqrt(rec["cP"] / rec["mu"])) if rec["cP"] > 0 else None,
                    "fit_points": int(np.sum(sel)), "max_split": e["max_split"]})
    return out


def main():
    out = {"rung": "R15-P-iii", "CP": CP, "gradient_gate": gradient_gate()}
    log(f"gradient gate {out['gradient_gate']}")
    out["slab_check"] = slab_check()
    log(f"slab check rel {out['slab_check']['rel']}")
    out["grid"] = {}
    for mu in (1e-3, 1e-2, 1e-1):
        for c in (0.1, 1.0, 10.0):
            key = f"mu{mu:g}_c{c:g}"
            th = theorem_check(mu, c)
            sc = onset_scan(mu, c)
            wall = ising_wall(mu, c)
            verdict = "CONTINUOUS_ONSET" if th["min_Veff"] >= -1e-12 else "FIRST_ORDER_CROSSING"
            out["grid"][key] = {"mu": mu, "c": c, "theorem": th, "onset": sc, "ising_wall": {k: v for k, v in wall.items() if k not in ("z", "profile_m2", "profile_m3")}, "verdict": verdict}
            log(f"{key}: min V_eff below onset {th['min_Veff']:.3e} at {th['argmin']} -> {verdict}; omega_c^2 {sc['omega_c2']:.3e}; "
                f"split at 0.99/1.01/2 omega_c^2: {sc['rows'][2]['split']:.4f}/{sc['rows'][4]['split']:.4f}/{sc['rows'][7]['split']:.4f}; "
                f"Ising wall at 2 omega_c^2: sigma {wall['sigma']:.4e} width {wall['width_10_90']} = {wall['width_over_sqrt_c_over_mu']} sqrt(c/mu)")
            if key == "mu0.001_c1":
                out["wall_profile_mu0.001_c1"] = {k: wall[k] for k in ("z", "profile_m2", "profile_m3")}
    out["decay_length_M_fields"] = decay_length_from_m_fields()
    for d in out["decay_length_M_fields"]:
        log(f"decay length: {d}")
    json.dump(out, open(os.path.join(C15.DATA, "m5_32_r15_p3_wall.json"), "w"), indent=1)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    for key, gval in out["grid"].items():
        rows = gval["onset"]["rows"]
        ax[0].plot([r["omega2_over_omega_c2"] for r in rows], [r["split"] for r in rows], marker="o", label=key)
    ax[0].set_xlabel("omega^2 / omega_c^2"); ax[0].set_ylabel("split of the lowest V_eff minimum"); ax[0].set_title("the onset (continuous at omega_c^2 = mu/c)")
    ax[0].legend(fontsize=6); ax[0].grid(alpha=0.3)
    for gval in out["grid"].values():
        iw = gval["ising_wall"]
        if iw["width_10_90"] is None:
            continue
        ok = iw["sigma"] > 0
        ax[1].loglog([iw["sqrt_c_over_mu"]], [iw["width_10_90"]], "o" if ok else "x", color="C0" if ok else "C3")
    ws = [(gval["ising_wall"]["sqrt_c_over_mu"], gval["ising_wall"]["width_10_90"]) for gval in out["grid"].values() if gval["ising_wall"]["width_10_90"] and gval["ising_wall"]["sigma"] > 0]
    xx = np.array([1.0, 100.0])
    ax[1].loglog(xx, xx * np.median([w[1] / w[0] for w in ws]), "k--", label="width = 1.15 sqrt(c/mu) (the 4 bounded points, o)")
    ax[1].loglog([], [], "x", color="C3", label="runaway points (x): width = the line length")
    ax[1].set_xlabel("sqrt(c/mu)"); ax[1].set_ylabel("diagonal-sector wall 10-90 width at 2 omega_c^2"); ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)
    wp = out["wall_profile_mu0.001_c1"]
    ax[2].plot(wp["z"], np.array(wp["profile_m2"]) - np.array(wp["profile_m3"]), label="split m2 - m3")
    ax[2].plot(wp["z"], np.array(wp["profile_m2"]) + np.array(wp["profile_m3"]) - 2 * DELTA, label="trace m2 + m3 - 2 delta")
    ax[2].set_xlabel("z"); ax[2].set_title("diagonal-sector wall, mu 1e-3, c 1, 2 omega_c^2 (a Goldstone saddle: audit)"); ax[2].legend(); ax[2].grid(alpha=0.3)
    fig.suptitle("M5.32 R15-P-iii: the reduced planar functional of L_P on the degenerate vacuum")
    fig.tight_layout()
    fig.savefig(os.path.join(C15.PLOTS, "m5_32_r15_p3_wall.png"), dpi=110)
    log("wrote data/m5_32_r15_p3_wall.json + plots/m5_32_r15_p3_wall.png")


if __name__ == "__main__":
    main()
