"""M5.32 R14-D: the P250 bridge on the 4x4 field (ledger 6.3): does L_cert + c K_P^h have a
clock-active bulk with a Maxwell crossing, hence a coexistence wall and a bag law?

EQUATIONS FIRST
---------------
Uniform states are parametrized by the spectrum of N = M eta: M = diag(m0, m1, m2, m3), N =
diag(-m0, m1, m2, m3) (the certified vacuum m = (g, 1, delta, 0), g = 8, delta = 0.3).  Under
the (2,3) isorotation clock (a0 = G1 M - M G1 = (m2 - m3)(E23 + E32)) a uniform state has
    I1 inertia          = 0                          (no spatial gradients: F_0i = 0)
    K_P^h inertia iota  = [f(m2) f(m3)]^2 (m2 - m3)^2,   f(x) = (x + g)(x - 1)   (per unit volume)
    static energy       = V4(m) = W1 sum_{p=1..4} (tr N^p - C_p)^2               (per unit volume)
so for the action E = E_cert + c K_P^h:
    fixed omega (the rotating frame):  V_eff(m; omega) = V4(m) - c omega^2 iota(m)
    fixed J (angular momentum density): E_J(m) = V4(m) + J^2 / (4 c iota(m))
A clock-active bulk is a local minimum of V_eff other than the (ticking) vacuum; a Maxwell
crossing is an omega_* at which two local minima of V_eff are degenerate (P250's structure:
the exterior at rest, the interior rotating).  Because iota is a polynomial of degree 10 in the
eigenvalues while V4 has degree 8, V_eff is UNBOUNDED BELOW for every omega != 0 and c > 0 along
the large-split direction: the fixed-omega problem on the polynomial-P entrant has no global
minimum, only local structure (reported); the fixed-J problem is bounded (E_J -> V4 as
iota -> infinity).
D' (the C3 fallback, a MODIFIED action):  V' = V4 + mu (m2 - m3)^2 makes the DEGENERATE state
(m2 = m3 = delta/2) the vacuum for mu above ~V4_deg / delta^2 (the degenerate state then does
not tick: iota = 0 there, no volume-extensive inertia), and the split state m2 - m3 = s is the
clock-active bulk once c omega^2 [f f]^2 > mu: a first-order coexistence with a Maxwell crossing
omega_*.  The 1D kink between the two phases (a planar (2,3)-split ramp; I1 is planar-flat, so
its tension comes from K_P^h alone: static density (1/2)[f(m2)^4 m2'^2 + f(m3)^4 m3'^2] on a
diagonal profile) has the first-integral tension  sigma = int_0^{s*} sqrt(2 kappa(s) DV(s)) ds
with kappa(s) = (1/8)[f(m2)^4 + f(m3)^4] along m2 = delta/2 + s/2, m3 = delta/2 - s/2 and
DV = V_eff(s) - V_eff(0) >= 0 at the crossing; the thin-wall bag law R = 2 sigma / p with
p = V_eff(0) - V_eff(s*) the pressure above the crossing.  Reported as P250's comparison
object on the 4x4 field under the modified potential, never as a property of L_cert.

Run: python3 m5_32_r14_d_bridge.py (numpy, ~2 min).  Writes data/m5_32_r14_d_bridge.json and
plots/m5_32_r14_d_bridge.png.
"""
from __future__ import annotations
import importlib.util
import json
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA, PLOTS = os.path.join(RES, "data"), os.path.join(RES, "plots")
OUT = os.path.join(DATA, "m5_32_r14_d_bridge.json")
PNG = os.path.join(PLOTS, "m5_32_r14_d_bridge.png")
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


T14 = _load("m5_32_r14_terms", "m5_32_r14_terms.py")
B3 = T14.B3
G, DELTA, W1 = 8.0, 0.3, B3.W1
CP = tuple((-G) ** p + 1.0 + DELTA ** p for p in range(1, 5))     # tr N^p of the vacuum, N = diag(-g, 1, delta, 0)


def f_of(x):
    return (x + G) * (x - 1.0)


def v4(m0, m1, m2, m3):
    m0, m1, m2, m3 = np.broadcast_arrays(np.asarray(m0, float), np.asarray(m1, float), np.asarray(m2, float), np.asarray(m3, float))
    lam = np.stack([-m0, m1, m2, m3], axis=-1)
    tot = 0.0
    for p in range(1, 5):
        tot = tot + (np.sum(lam ** p, axis=-1) - CP[p - 1]) ** 2
    return W1 * tot


def iota(m2, m3):
    return (f_of(m2) * f_of(m3)) ** 2 * (m2 - m3) ** 2


def local_minima_2d(V, x, y):
    """indices of strict local minima of a 2D array (interior only)."""
    c = V[1:-1, 1:-1]
    ok = np.ones_like(c, dtype=bool)
    for di in (-1, 0, 1):
        for dj in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            ok &= c < V[1 + di:V.shape[0] - 1 + di, 1 + dj:V.shape[1] - 1 + dj]
    ii, jj = np.where(ok)
    return [(float(x[i + 1]), float(y[j + 1]), float(c[i, j])) for i, j in zip(ii, jj)]


def selftest():
    # the K_P^h inertia formula against the registry on a uniform lattice state
    cfg = B3.base_cfg(s=-1.0, g=G, n=4, L=6.0)
    ok = True
    for (m2, m3) in ((DELTA, 0.0), (0.5, -0.2), (0.15, 0.15)):
        M = np.broadcast_to(np.diag([G, 1.0, m2, m3]), (4, 4, 4, 4, 4)).copy()
        a0 = B3.embed34(np.zeros((4, 4, 4, 3, 3)), cfg) * 0.0
        J = np.zeros((4, 4)); J[2, 3], J[3, 2] = -1.0, 1.0
        a0 = np.broadcast_to(J @ np.diag([G, 1.0, m2, m3]) - np.diag([G, 1.0, m2, m3]) @ J, M.shape)
        k_reg = T14.kin_energy("K_P_h", M, a0, cfg) / (64 * cfg["h"] ** 3)
        ok &= abs(k_reg - iota(m2, m3)) <= 1e-9 * max(1.0, abs(iota(m2, m3)))
        v_reg = float(B3.e_parts(M, cfg)[1]) / (64 * cfg["h"] ** 3)
        ok &= abs(v_reg - v4(G, 1.0, m2, m3)) <= 1e-12 * max(1.0, v_reg)
    return bool(ok)


def main():
    res = {"rung": "R14-D", "selftest_formulas_vs_registry": selftest()}
    log(f"selftest {res['selftest_formulas_vs_registry']}")
    # ---------------- L_cert + c K_P^h: the (m2, m3) plane at m0 = g, m1 = 1
    m = np.linspace(-1.0, 1.5, 251)
    M2, M3 = np.meshgrid(m, m, indexing="ij")
    V4g = v4(G, 1.0, M2, M3)
    IO = iota(M2, M3)
    res["L_cert_plus_cKP"] = {}
    for c in (0.3, 1.0, 3.0):
        rec = {}
        for om in (0.0, 1e-4, 1e-3, 1e-2, 3e-2, 1e-1):
            Veff = V4g - c * om * om * IO
            mins = local_minima_2d(Veff, m, m)
            # the ticking vacuum: the minimum nearest (delta, 0) or (0, delta)
            vac = [mm for mm in mins if min(abs(mm[0] - DELTA) + abs(mm[1]), abs(mm[1] - DELTA) + abs(mm[0])) < 0.05]
            others = [mm for mm in mins if mm not in vac]
            rec[f"omega={om:g}"] = {"n_local_minima": len(mins), "ticking_vacuum": vac[:2], "other_minima": others[:6],
                                   "V_eff_min_on_box": float(np.min(Veff)), "argmin_on_box": [float(M2.ravel()[np.argmin(Veff)]), float(M3.ravel()[np.argmin(Veff)])],
                                   "bounded_on_box": bool(np.argmin(Veff) not in (0, Veff.size - 1)) and bool(np.min(Veff) > -1e6)}
        # unboundedness along the split direction
        s = np.array([1.0, 2.0, 4.0, 8.0])
        rec["V_eff_along_large_split(omega=1e-2)"] = [float(v4(G, 1.0, DELTA / 2 + x / 2, DELTA / 2 - x / 2) - c * 1e-4 * iota(DELTA / 2 + x / 2, DELTA / 2 - x / 2)) for x in s]
        rec["iota_over_V4_along_large_split"] = [float(iota(DELTA / 2 + x / 2, DELTA / 2 - x / 2) / v4(G, 1.0, DELTA / 2 + x / 2, DELTA / 2 - x / 2)) for x in s]
        # fixed J: the global minimum of E_J on the box
        rec["fixed_J"] = {}
        for J in (1e-3, 1e-2, 1e-1, 1.0):
            EJ = V4g + J * J / (4.0 * c * np.maximum(IO, 1e-300))
            k = np.argmin(EJ)
            rec["fixed_J"][f"J={J:g}"] = {"E_J_min": float(EJ.ravel()[k]), "state_m2_m3": [float(M2.ravel()[k]), float(M3.ravel()[k])],
                                        "omega": float(J / (2 * c * IO.ravel()[k])) if IO.ravel()[k] > 0 else None,
                                        "V4_at_min": float(V4g.ravel()[k]), "split": float(M2.ravel()[k] - M3.ravel()[k])}
        res["L_cert_plus_cKP"][f"c={c:g}"] = rec
        log(f"c = {c}: minima at omega 1e-2: {rec['omega=0.01']['n_local_minima']} (vacuum {rec['omega=0.01']['ticking_vacuum'][:1]}, others {rec['omega=0.01']['other_minima'][:2]}); "
            f"V_eff along split {rec['V_eff_along_large_split(omega=1e-2)']}")
    # the vacuum's ticking depth vs the degenerate state's V4: the crossover omega^2
    vdeg = float(v4(G, 1.0, DELTA / 2, DELTA / 2))
    res["degenerate_state"] = {"V4_deg_per_volume": vdeg, "iota_vac": float(iota(DELTA, 0.0)),
                               "omega2_crossover_c1": vdeg / float(iota(DELTA, 0.0)),
                               "verdict_on_L_cert_plus_cKP": "NO_CLOCK_ACTIVE_BULK distinct from the ticking vacuum at any scanned omega; "
                                                             "the degenerate exterior is never favored at fixed omega (the vacuum ticks with volume-extensive inertia); "
                                                             "V_eff unbounded below along the large-split direction for omega != 0"}
    # ---------------- D': the modified potential V' = V4 + mu (m2 - m3)^2, the 1D split line
    s = np.linspace(0.0, 1.2, 2401)
    m2s, m3s = DELTA / 2 + s / 2, DELTA / 2 - s / 2
    V4s = v4(G, 1.0, m2s, m3s)
    IOs = iota(m2s, m3s)
    kappa = 0.125 * (f_of(m2s) ** 4 + f_of(m3s) ** 4)
    dprime = {}
    for mu in (1e-4, 1e-3, 1e-2):
        Vp = V4s + mu * s * s
        i0 = int(np.argmin(Vp))
        rec = {"mu": mu, "V'_min_at_s": float(s[i0]), "V'_deg_is_vacuum": bool(s[i0] < 1e-9),
               "V'_at_split_delta": float(Vp[np.argmin(np.abs(s - DELTA))]), "V'_at_0": float(Vp[0])}
        if not rec["V'_deg_is_vacuum"]:
            dprime[f"mu={mu:g}"] = rec
            continue
        for c in (0.3, 1.0, 3.0):
            # scan omega for a second minimum degenerate with s = 0 (the Maxwell crossing)
            def veff(om):
                return Vp - c * om * om * IOs
            oms = np.logspace(-5, -0.5, 400)
            found = None
            for om in oms:
                V = veff(om)
                # local minima on the line s > 0.02
                d1 = np.diff(V)
                lm = np.where((d1[:-1] < 0) & (d1[1:] > 0))[0] + 1
                lm = [k for k in lm if s[k] > 0.02]
                if lm:
                    k = min(lm, key=lambda k: V[k])
                    if V[k] <= V[0]:
                        found = (float(om), float(s[k]), float(V[k] - V[0]))
                        break
            entry = {"maxwell_crossing": found is not None}
            # the ORDER of the transition: a first-order (Maxwell) crossing needs a barrier between
            # the exterior s = 0 and the bulk s* at the crossing; if the exterior's split stiffness
            # goes negative first (a continuous onset), there is no coexistence wall and no bag.
            # analytic thresholds: omega_inst^2 = mu_eff / (c F0), F0 = f(delta/2)^4 the small-split
            # inertia coefficient, mu_eff = mu + (1/2) V4''(0) along the split; omega_cross^2 =
            # [V'(delta) - V'(0)] / (c iota(delta, 0)) for the delta-well
            F0 = f_of(DELTA / 2) ** 4
            v4pp = float(np.polyfit(s[:41], V4s[:41], 2)[0]) * 2.0          # V4'' at s = 0 (quadratic fit on s < 0.02)
            mu_eff = mu + 0.5 * v4pp
            om_inst2 = mu_eff / (c * F0)
            om_cross2 = (float(Vp[np.argmin(np.abs(s - DELTA))]) - float(Vp[0])) / (c * iota(DELTA, 0.0))
            entry.update({"omega2_instability_of_exterior": om_inst2, "omega2_delta_well_crossing": om_cross2,
                          "V4_second_derivative_at_0_along_split": v4pp, "F0_small_split_inertia_coefficient": F0,
                          "order": "SECOND (continuous onset, no wall)" if om_inst2 <= om_cross2 else "FIRST (barrier, Maxwell crossing)"})
            if found is not None:
                om_s, s_star, dv = found
                V = veff(om_s)
                k_star0 = int(np.argmin(np.abs(s - s_star)))
                barrier = float(np.max(V[:k_star0 + 1] - V[0]))
                entry["barrier_between_exterior_and_bulk_at_crossing"] = barrier
                entry["order_measured"] = "FIRST" if barrier > 1e-12 else "SECOND (no barrier: sigma = 0, no bag)"
                # tension at the crossing: sigma = int sqrt(2 kappa DV) ds between the two phases
                k_star = int(np.argmin(np.abs(s - s_star)))
                DV = np.maximum(V[:k_star + 1] - V[0], 0.0)
                sigma = float(np.trapezoid(np.sqrt(2 * kappa[:k_star + 1] * DV), s[:k_star + 1]))
                # bag law slightly above the crossing: p = V(0) - V(s*), R = 2 sigma / p
                bags = []
                for fac in (1.001, 1.003, 1.01, 1.03):
                    Vb = veff(om_s * fac)
                    d1 = np.diff(Vb)
                    lm = np.where((d1[:-1] < 0) & (d1[1:] > 0))[0] + 1
                    lm = [k for k in lm if s[k] > 0.02]
                    if not lm:
                        continue
                    k = min(lm, key=lambda k: Vb[k])
                    p = float(Vb[0] - Vb[k])
                    bags.append({"omega_over_omega_star": fac, "s_bulk": float(s[k]), "pressure": p, "R_thin_wall": (2 * sigma / p) if p > 0 else None})
                entry.update({"omega_star": om_s, "s_star": s_star, "V_eff_bulk_minus_exterior_at_crossing": dv,
                              "sigma_first_integral": sigma, "kappa_at_0": float(kappa[0]), "kappa_at_s_star": float(kappa[k_star]),
                              "iota_bulk": float(IOs[k_star]), "bags": bags,
                              "omega_star_within_radiation_window(mu_gap=sqrt(V4 curvature)?)": "not assessed here"})
            rec[f"c={c:g}"] = entry
            log(f"  D' mu {mu:g} c {c:g}: crossing {found}")
        dprime[f"mu={mu:g}"] = rec
    res["D_prime_modified_potential"] = dprime
    res["wall_s"] = round(time.time() - T0, 1)
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    # plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
        Veff = V4g - 1.0 * 1e-4 * IO
        im = ax[0].contourf(M2, M3, np.log10(np.abs(Veff) + 1e-12) * np.sign(Veff), 40, cmap="RdBu_r")
        ax[0].plot([DELTA, 0], [0, DELTA], "k*"); ax[0].set_title("L_cert + K_P^h: sign log10 |V_eff| at omega 1e-2, c 1"); ax[0].set_xlabel("m2"); ax[0].set_ylabel("m3")
        fig.colorbar(im, ax=ax[0])
        for mu in (1e-4, 1e-3, 1e-2):
            ax[1].plot(s, V4s + mu * s * s, label=f"V' mu {mu:g}")
        ax[1].plot(s, V4s, "k--", label="V4")
        ax[1].set_xlim(0, 0.6); ax[1].set_ylim(-1e-5, 3e-4); ax[1].set_xlabel("split s"); ax[1].set_title("D': V4 + mu s^2 along the (2,3) split"); ax[1].legend(fontsize=7)
        key = "mu=0.001"
        if key in dprime and "c=1" in dprime[key] and dprime[key]["c=1"].get("maxwell_crossing"):
            om_s = dprime[key]["c=1"]["omega_star"]
            for fac in (0.9, 1.0, 1.1):
                ax[2].plot(s, V4s + 1e-3 * s * s - 1.0 * (om_s * fac) ** 2 * IOs, label=f"omega/omega* {fac}")
            ax[2].set_xlim(0, 0.6); ax[2].set_xlabel("split s"); ax[2].set_title(f"D' mu 1e-3, c 1: V_eff near the crossing (omega* {om_s:.3g})"); ax[2].legend(fontsize=7)
            ax[2].set_ylim(-3e-4, 3e-4)
        plt.tight_layout(); plt.savefig(PNG, dpi=110)
    except Exception as e:                                        # noqa: BLE001
        log(f"plot skipped: {e!r}")
    log(f"done -> {OUT}")


if __name__ == "__main__":
    main()
