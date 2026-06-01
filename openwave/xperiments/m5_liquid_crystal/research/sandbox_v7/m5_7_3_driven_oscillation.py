"""
M5.7.3 — 9b.1 driven-oscillation preview: can an EM-like lever sustain the defect's (A, ω)?

M5.7.1 (seeded) + M5.7.2 (intrinsic) established that a *freely-evolving* 3D defect DISPERSES
its excess orientation energy — no free metastable clock (that's 4D / M5.8). But heat in the
SABER hypothesis is a *maintained* (A, ω) excess — sustained by a continuous drive (a bath / an
EM lever), not a one-shot free oscillation. So the question this previews (the heart of 9b.1):

    Does a continuous EM-wave-like drive hold the defect in an excited (A, ω) state
    (vs the free dispersal), and at what frequency does the lever couple best?

THE DRIVE (physically grounded).  From M5.6.4 EM = director tilts; from M5.7.1 "rotating the
director IS what an EM-wave lever does." So the drive is a continuous time-periodic director-
rotation forcing with a FIXED spatial pattern (an incident-wave's effect, localized to the
defect's textured shell):

    D̂(x)   = w(x) · [G_y, M_bg(x)]          (so(3) rotation tangent about ŷ, shell-localized)
    F_drive(x,t) = A_drive · sin(2π f_d t) · D̂(x)            added to the Eq.18 acceleration
    EOM:   M̈ = c²·div(G) − dV_M(M) + F_drive

`[G_y, M]` (antisymmetric·symmetric commutator) is **symmetric + traceless** → keeps M symmetric,
drives the ORIENTATION sector only (V is amplitude-only, so it stays V-flat). `w(x)` is a Gaussian
shell at r≈2.5r₀ (the active textured region), so the lever acts ON the defect, not the uniform
far field.

METRICS.
  • A_core(t) = ⟨‖M − M_bg‖_F⟩ over the core shell — the EXCITATION amplitude (how far the drive
    pushes the defect from rest). Free (A_drive=0): the M5.7.2 self-sourced-then-disperse curve.
    Driven: held elevated if the lever couples. **Sustained A_core = mean over the last third.**
  • Response: FFT the core director — does it ring at the drive f_d (forced response)?
  • Resonance: sweep f_d → the f_d that maximizes sustained A_core = the defect's natural mode
    (the frequency a real EM lever should use — the 9b.1 "FM" handle).
  • H(t) grows under drive (closed Dirichlet box has no energy sink — the box fills up; that is
    NOT a failure, just not the metric). We track it only to flag numerical blow-up.

BRIDGE TO SABER.  A sustained, frequency-selective driven (A, ω) is the field-theoretic basis of
"heat = maintained (A, ω) excess" + identifies the modulation frequency the EM lever should use.
(Physics framing only here; device/engineering specifics → SABER docs.)

USAGE:
    python -m openwave.xperiments.m5_liquid_crystal.research.sandbox_v7.m5_7_3_driven_oscillation
    M57_N=32 M57_STEPS=4000 M57_ADRIVE=0.02 python -m ...m5_7_3_driven_oscillation   # fast smoke
"""

import os

import numpy as np

# ----- substrate (matches _topo_biaxial1_von.py + M5.7.1/.2) ------------------------
DTYPE = np.float32
DELTA = 0.30
S2_STAR = 1.0 + DELTA**2
N = int(os.environ.get("M57_N", 44))
L = 5.0
RC = 0.9
RHOC = 0.9
C_WAVE = 0.3
C_COEF = 2.0
A_COEF = -2.0 * C_COEF * S2_STAR
DT = float(os.environ.get("M57_DT", 0.004))
N_STEPS = int(os.environ.get("M57_STEPS", 5000))
SAMPLE_EVERY = 5
R_CORE = 3.0 * RC
SHELL_PEAK = 2.5 * RC  # drive + measurement shell (active textured region)
SHELL_W = 0.7
A_DRIVE = float(os.environ.get("M57_ADRIVE", 0.02))  # peak drive acceleration (env-tunable)

G_Y = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], DTYPE)  # so(3) gen about ŷ


# ----- matrix helpers (numpy mirror of the @ti.func ops; identical to M5.7.1/.2) ----
def matmul(A, B):
    return np.matmul(A, B)


def commf(A, B):
    return np.matmul(A, B) - np.matmul(B, A)


def frob2(A):
    return np.einsum("...ab,...ab->...", A, A)


def tr(A):
    return np.einsum("...aa->...", A)


def central(f, axis, h):
    out = np.zeros_like(f)
    sp, sm, so = [slice(None)] * f.ndim, [slice(None)] * f.ndim, [slice(None)] * f.ndim
    sp[axis], sm[axis], so[axis] = slice(2, None), slice(0, -2), slice(1, -1)
    out[tuple(so)] = (f[tuple(sp)] - f[tuple(sm)]) / (2 * h)
    return out


def V_M(M, a, b, c):
    M2 = matmul(M, M)
    tr2 = tr(M2)
    return a * tr2 - b * tr(matmul(M2, M)) + c * tr2 * tr2


def dV_M(M, a, b, c):
    M2 = matmul(M, M)
    tr2 = tr(M2)[..., None, None]
    return 2.0 * a * M - 3.0 * b * M2 + 4.0 * c * tr2 * M


def curvature_force(M, h):
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    cc = lambda A, B: commf(commf(A, B), B)
    Gx = 8.0 * (cc(Mx, My) + cc(Mx, Mz))
    Gy = 8.0 * (cc(My, Mx) + cc(My, Mz))
    Gz = 8.0 * (cc(Mz, Mx) + cc(Mz, My))
    return central(Gx, 0, h) + central(Gy, 1, h) + central(Gz, 2, h)


def energy_total(M, M_prev, h):
    Mdot = (M - M_prev) / DT
    Mx, My, Mz = central(M, 0, h), central(M, 1, h), central(M, 2, h)
    curv = 4.0 * (frob2(commf(Mx, My)) + frob2(commf(Mx, Mz)) + frob2(commf(My, Mz)))
    return 0.5 * frob2(Mdot) + C_WAVE**2 * curv + V_M(M, A_COEF, 0.0, C_COEF)


# ----- biaxial hedgehog background (mirror of seed_biaxial_hedgehog_M) ---------------
def build_biaxial_M():
    xs = np.linspace(-L, L, N)
    h = xs[1] - xs[0]
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2) + 1e-9
    rho = np.sqrt(X**2 + Y**2) + 1e-9
    rhat = np.stack([X / r, Y / r, Z / r], -1)
    azim = np.stack([-Y / rho, X / rho, np.zeros_like(r)], -1)
    sdisc = np.clip(rho / RHOC, 0, 1)
    shrink = (sdisc**2 * (3 - 2 * sdisc))[..., None]
    ephi = azim * shrink
    ephi = ephi - np.einsum("...a,...a->...", ephi, rhat)[..., None] * rhat
    etheta = np.cross(ephi, rhat)
    srad = (r / np.sqrt(r**2 + RC**2))[..., None, None]
    diso = (1.0 + DELTA) / 3.0
    d0 = diso + srad[..., 0, 0] * (1.0 - diso)
    d1 = diso + srad[..., 0, 0] * (DELTA - diso)
    d2 = diso + srad[..., 0, 0] * (0.0 - diso)
    oprod = lambda v: v[..., :, None] * v[..., None, :]
    M = (
        d0[..., None, None] * oprod(rhat)
        + d1[..., None, None] * oprod(etheta)
        + d2[..., None, None] * oprod(ephi)
    )
    interior = np.zeros(r.shape, bool)
    interior[2:-2, 2:-2, 2:-2] = True
    active = (r > 2 * RC) & (rho > RHOC) & interior
    shell = active & (np.abs(r - SHELL_PEAK) < 2 * SHELL_W)  # core-shell for A_core + drive
    ix = int(round((SHELL_PEAK + L) / (2 * L) * (N - 1)))
    ic = (N - 1) // 2
    return dict(
        M=M.astype(DTYPE),
        h=DTYPE(h),
        r=r.astype(DTYPE),
        active=active,
        shell=shell,
        probe=(ix, ic, ic),
    )


def build_drive_pattern(bg):
    """D̂(x) = w(x)·[G_y, M_bg], symmetric+traceless, Gaussian-shell localized, peak-normalized."""
    M, r, active = bg["M"], bg["r"], bg["active"]
    Dh = commf(np.broadcast_to(G_Y, M.shape), M)  # [G_y, M_bg] (symmetric, traceless)
    w = np.exp(-((r - SHELL_PEAK) ** 2) / (2 * SHELL_W**2))[..., None, None]
    Dh = (Dh * w) * active[..., None, None]
    peak = np.sqrt(frob2(Dh)).max() + 1e-30
    return (Dh / peak).astype(DTYPE)  # peak ‖D̂‖_F = 1 over the shell


def principal_director(m3):
    w, v = np.linalg.eigh(m3.astype(np.float64))
    return v[:, np.argmax(w)]


# ----- lockstep evolve: free + driven fields, F_drive added per field ---------------
def run_lockstep(bg, freqs, n_steps):
    """freqs: list of drive frequencies f_d (0.0 = free baseline). Records per field:
    A_core(t) = ⟨‖M−M_bg‖_F⟩ over the shell, H_total(t), core director probe n̂(t)."""
    h, active, shell, probe = bg["h"], bg["active"], bg["shell"], bg["probe"]
    mask = active[..., None, None]
    Mbg = bg["M"]
    Dh = build_drive_pattern(bg)
    series = [dict(A=[], H=[], ndir=[]) for _ in freqs]
    cur = [Mbg.copy() for _ in freqs]
    prev = [Mbg.copy() for _ in freqs]
    tsteps = []

    def sample():
        for fi in range(len(freqs)):
            dM = cur[fi] - Mbg
            series[fi]["A"].append(float(np.sqrt(frob2(dM))[shell].mean()))
            series[fi]["H"].append(float((energy_total(cur[fi], prev[fi], h) * active).sum()))
            series[fi]["ndir"].append(principal_director(cur[fi][probe]))
        tsteps.append(len(tsteps) * SAMPLE_EVERY)

    sample()
    for step in range(1, n_steps + 1):
        t = step * DT
        for fi, fd in enumerate(freqs):
            f = C_WAVE**2 * curvature_force(cur[fi], h) - dV_M(cur[fi], A_COEF, 0.0, C_COEF)
            if fd > 0.0:
                f = f + (A_DRIVE * np.sin(2 * np.pi * fd * t)) * Dh  # the EM-like lever
            new = np.where(mask, 2 * cur[fi] - prev[fi] + DT**2 * f, Mbg)
            prev[fi], cur[fi] = cur[fi], new
        if step % SAMPLE_EVERY == 0:
            sample()
            if step % 1000 == 0:
                As = " ".join(f"f{freqs[i]}:A={series[i]['A'][-1]:.3f}" for i in range(len(freqs)))
                print(f"      step {step:>5}/{n_steps}  {As}", flush=True)
    finite = all(np.isfinite(c).all() for c in cur)
    out = []
    for fi, fd in enumerate(freqs):
        A = np.array(series[fi]["A"])
        H = np.array(series[fi]["H"])
        ndir = np.array(series[fi]["ndir"])
        tail = slice(2 * len(A) // 3, None)  # sustained = mean over last third
        out.append(
            dict(
                fd=fd,
                A=A,
                A_sustained=float(A[tail].mean()),
                H_growth=float((H[-1] - H[0]) / (abs(H[0]) + 1e-30)),
                ndir=ndir,
            )
        )
    return out, finite


def response_at(ndir, f_d):
    """Power fraction of the core-director spectrum within ±2 bins of the drive frequency f_d."""
    dt_s = DT * SAMPLE_EVERY
    n = ndir.shape[0]
    win = np.hanning(n)
    freqs = np.fft.rfftfreq(n, d=dt_s)
    k_d = int(np.argmin(np.abs(freqs - f_d)))
    best = 0.0
    for comp in range(3):
        sig = ndir[:, comp] - ndir[:, comp].mean()
        if np.std(sig) < 1e-9:
            continue
        spec = np.abs(np.fft.rfft(sig * win)) ** 2
        spec[0] = 0.0
        if spec.sum() < 1e-30:
            continue
        lo, hi = max(0, k_d - 2), min(len(spec), k_d + 3)
        best = max(best, spec[lo:hi].sum() / spec.sum())  # AC power near the drive freq
    return best


# ====================================================================================
def main():
    print("=" * 82)
    print(
        "M5.7.3 — 9b.1 driven-oscillation preview: can an EM-like lever sustain the defect (A,ω)?"
    )
    print(
        f"  grid {N}³  L={L}  c={C_WAVE}  dt={DT}  steps={N_STEPS}  δ={DELTA}  A_drive={A_DRIVE}  f32"
    )
    print("=" * 82)
    bg = build_biaxial_M()
    # 0.0 = free baseline; sweep around the M5.7.2 natural ~0.1–0.25/t. Env-overridable: M57_FREQS="0.0,0.1,0.2"
    freqs = [float(x) for x in os.environ.get("M57_FREQS", "0.0,0.10,0.20,0.40").split(",")]
    print(
        f"  active {100*bg['active'].mean():.0f}%  shell voxels {bg['shell'].sum()}  "
        f"probe {bg['probe']}  drive freqs {freqs[1:]}/t",
        flush=True,
    )

    res, finite = run_lockstep(bg, freqs, N_STEPS)
    free = res[0]
    driven = res[1:]

    print(f"\n  RESULTS  (A_core = ⟨‖M−M_bg‖⟩ over the shell; sustained = mean over last third):")
    print(
        f"    [FREE   f_d=0.0 ]  A_sustained={free['A_sustained']:.4f}   "
        f"A: {free['A'][0]:.3f}→max {free['A'].max():.3f}→end {free['A'][-1]:.3f}"
    )
    for rr in driven:
        resp = response_at(rr["ndir"], rr["fd"])
        rr["resp"] = resp
        gain = rr["A_sustained"] / (free["A_sustained"] + 1e-30)
        print(
            f"    [DRIVEN f_d={rr['fd']:.2f}]  A_sustained={rr['A_sustained']:.4f}  "
            f"({gain:.1f}× free)  response@f_d={100*resp:.0f}%  H-growth {100*rr['H_growth']:+.0f}%"
        )

    best = max(driven, key=lambda d: d["A_sustained"])
    best_gain = best["A_sustained"] / (free["A_sustained"] + 1e-30)
    print("\n" + "=" * 82)
    print(f"M5.7.3  (N={N}, {N_STEPS} steps, A_drive={A_DRIVE}):  finite={finite}")
    # SUSTAINED: a driven field holds A_core clearly above the free dispersal
    sustained = best_gain >= 1.5 and finite
    # RESONANT: the best frequency is frequency-selective (responds at its f_d, beats neighbors)
    resonant = best.get("resp", 0) >= 0.3
    print(
        f"  SUSTAINED by the lever? {'YES' if sustained else 'NO'} — best f_d={best['fd']:.2f} holds "
        f"A_core {best_gain:.1f}× the free baseline"
    )
    print(
        f"  FREQUENCY-SELECTIVE? {'YES' if resonant else 'WEAK'} — best f_d responds "
        f"{100*best.get('resp',0):.0f}% at its drive frequency (the lever's coupling band)"
    )
    if sustained and resonant:
        print(
            f"  → THE EM LEVER WORKS: a continuous drive at f_d≈{best['fd']:.2f}/t holds the defect in a"
        )
        print(
            "    maintained (A,ω) excess — the field-theoretic basis of heat-as-driven-excess (9b.1)."
        )
    elif sustained:
        print(
            "  → The lever sustains the excitation but the resonance is broad — map the band in 9b.1."
        )
    else:
        print(
            "  → The lever does not clearly sustain the excitation at this A_drive — tune amplitude /"
        )
        print("    frequency, or the sustained-thermal-state needs the 4D clock substrate (M5.8).")
    print("PASS" if finite else "FAIL — blow-up (reduce A_drive or dt)")
    print("=" * 82)
    return 0 if finite else 1


if __name__ == "__main__":
    raise SystemExit(main())
