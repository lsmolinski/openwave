"""M5.32 R15-P-ii: the K_P^23 tail on the relaxed L_P hedgehogs (the author's prediction (ii): finite tail,
L-exponent 0, against the certified K_P's 1.34 measured in R14-0).

Equations first.  On the R15-M end fields (n32 L48 and n48 L72 at the same spacing 1.5, mu in {0, 1e-2},
c_P in {0, 1}) the h^3-weighted static densities per cell
    e_KP(x) = (1/2) sum_i tr(Om_i^T eta Om_i eta),  Om_i = P23 A_i eta P23     (raw, unweighted by c_P)
    e_u(x)  = 4 sum_{i<j} <F_ij, F_ij>_eta
are shell-averaged in r and fitted as  e(r) ~ r^-q  over  6 <= r <= L/2 - 3  (log-log least squares); the
box-integrated energies E(L) give the L-exponent  p = ln[E(72) / E(48)] / ln(72 / 48)  for K_P^23 static, E_u,
V4^dd and the split term.  A Coulomb hedgehog has e_u ~ r^-4 (q = 4, p -> 0); an L-divergent tail has q <= 3
(p = 3 - q > 0).  Verdict per (mu, c_P): TAIL_FINITE if q_KP > 3 and p_KP < 0.3; TAIL_DIVERGENT otherwise.

usage: python3 m5_32_r15_p2_tail.py
"""
import sys
ARGS = list(sys.argv[1:])
import os, json, glob
import numpy as np
import m5_32_r15_common as C15
import m5_32_r15_m_hedgehog as R15M

INS4 = C15.INS4
log = C15.log
MH = os.path.join(C15.CK, "m_hedgehog")


def shell_profile(dens, r, L, nb=16):
    edges = np.linspace(0.0, L / 2, nb + 1)
    rc, ev = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (r >= a) & (r < b)
        if np.any(m):
            rc.append(0.5 * (a + b)); ev.append(float(np.mean(dens[m])))
    return np.array(rc), np.array(ev)


def fit_q(rc, ev, L):
    sel = (rc >= 6.0) & (rc <= L / 2 - 3.0) & (ev > 0)
    if np.sum(sel) < 3:
        return None
    slope = np.polyfit(np.log(rc[sel]), np.log(ev[sel]), 1)[0]
    return float(-slope)


def main():
    recs = {}
    for f in sorted(glob.glob(os.path.join(MH, "relax_*.json"))):
        rec = json.load(open(f))
        recs[(rec["n"], rec["L"], rec["mu"], rec["cP"])] = rec
    out = {"rung": "R15-P-ii", "R14_0_reference": "the certified K_P (K_P^h) static tail on the L_cert hedgehog: L-exponent 1.34 (R14-0)", "profiles": {}, "exponents": {}}
    for (n, L, mu, cP), rec in recs.items():
        cfg = C15.cfg_dd(n, L, mu=mu, cP=cP)
        M = np.load(rec["field"])
        X, Y, Z = INS4.coords(n, cfg["h"])
        r = np.sqrt(X * X + Y * Y + Z * Z)
        d = R15M.static_density(M, cfg)
        dkp = C15.kp23_static_density(M, cfg)
        key = f"n{n}_L{L:g}_mu{mu:g}_cP{cP:g}"
        prof = {}
        for name, dd in (("KP23_raw", dkp), ("E_u", d["E_u"]), ("V4", d["V4"])):
            rc, ev = shell_profile(dd / cfg["h"] ** 3, r, L)
            prof[name] = {"r": rc.tolist(), "density": ev.tolist(), "q": fit_q(rc, ev, L)}
        out["profiles"][key] = prof
        log(f"{key}: q_KP {prof['KP23_raw']['q']}, q_Eu {prof['E_u']['q']}, q_V4 {prof['V4']['q']}")
    for mu in (0.0, 1e-2):
        for cP in (0.0, 1.0):
            a = recs.get((32, 48.0, mu, cP)); b = recs.get((48, 72.0, mu, cP))
            if a is None or b is None:
                continue
            pa, pb = a["end"]["parts"], b["end"]["parts"]
            ex = {}
            for name in ("KP_static_raw", "E_u", "V4", "split"):
                va, vb = pa[name], pb[name]
                ex[name] = float(np.log(vb / va) / np.log(1.5)) if va > 0 and vb > 0 else None
                ex[name + "_values"] = [va, vb]
            qk = out["profiles"][f"n48_L72_mu{mu:g}_cP{cP:g}"]["KP23_raw"]["q"]
            pk = ex["KP_static_raw"]
            ver = "TAIL_FINITE" if (qk is not None and qk > 3.0 and pk is not None and pk < 0.3) else "TAIL_DIVERGENT"
            ex["verdict"] = ver
            out["exponents"][f"mu{mu:g}_cP{cP:g}"] = ex
            log(f"mu {mu:g} cP {cP:g}: L-exponents KP {pk}, E_u {ex['E_u']}, V4 {ex['V4']}, split {ex['split']}; q_KP(n48) {qk} -> {ver}")
    json.dump(out, open(os.path.join(C15.DATA, "m5_32_r15_p2_tail.json"), "w"), indent=1)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for key, prof in out["profiles"].items():
        ax[0].loglog(prof["KP23_raw"]["r"], np.maximum(prof["KP23_raw"]["density"], 1e-14), marker="o", ms=3, label=f"{key} q {prof['KP23_raw']['q'] if prof['KP23_raw']['q'] is None else round(prof['KP23_raw']['q'], 2)}")
        ax[1].loglog(prof["E_u"]["r"], np.maximum(prof["E_u"]["density"], 1e-14), marker="o", ms=3, label=f"{key} q {prof['E_u']['q'] if prof['E_u']['q'] is None else round(prof['E_u']['q'], 2)}")
    rr = np.array([6.0, 30.0])
    for a_ in ax:
        a_.loglog(rr, 1e-2 * (rr / 6.0) ** -4, "k:", label="r^-4"); a_.loglog(rr, 1e-2 * (rr / 6.0) ** -3, "k--", label="r^-3")
        a_.set_xlabel("r"); a_.grid(alpha=0.3); a_.legend(fontsize=6)
    ax[0].set_title("K_P^23 static density (raw), shell mean"); ax[1].set_title("E_u density, shell mean")
    fig.suptitle("M5.32 R15-P-ii: tails on the relaxed L_P hedgehogs")
    fig.tight_layout(); fig.savefig(os.path.join(C15.PLOTS, "m5_32_r15_p2_tail.png"), dpi=110)
    log("wrote data/m5_32_r15_p2_tail.json + plots/m5_32_r15_p2_tail.png")


if __name__ == "__main__":
    main()
