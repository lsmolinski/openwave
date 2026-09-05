"""M5.32 R14: the two-derivative ENTRANTS of the frozen basis B (ledger 6.3),
as a companion registry of m5_32_lagrangian.py (imported, never modified).

EQUATIONS FIRST
---------------
Conventions (the run's locked table): M(x) real symmetric 4x4 per point,
eta = diag(-1, 1, 1, 1), index 0 = time; raw entries of M and of the jets
A_mu = d_mu M are CONTRAVARIANT internal tensors, d_mu is covariant;
A_0 = omega a0 on a clock.  N := M eta is the (1,1) endomorphism whose
spectrum is Lorentz-invariant; at the certified vacuum M = diag(g, 1, delta, 0)
(code branch s = -1) the spectrum of N is (-g, 1, delta, 0): the TIMELIKE
eigenvalue is lambda_t = -g in this convention.

The three entrants named by the model author (coordination thread, 2026-09-03),
each written with an explicit sign so that a POSITIVE coefficient is the
energy-positive orientation (the same orientation as the certified -4 I1, whose
energy is +4 (U + omega^2 T)):

  K_lambda  E-density = (1/2) sum_a [ sum_i (d_i lambda_a)^2 + omega^2 (d_t lambda_a)^2 ]
            lambda_a the eigenvalues of N;  d_mu lambda_a = (v_a^T eta A_mu eta v_a) / (v_a^T eta v_a)
            (first-order perturbation theory; on a degenerate cluster the sum of the
            squared one-sided derivatives is tr(B^2) of the cluster block B of the
            eta-projected jet).  Inert on every orientation gradient (the spectrum is
            constant on a Lorentz orbit) and on every generator channel (d_t lambda = 0
            for rotations and boosts alike), so it contributes NOTHING to any omega^2
            form: it is a static eigenvalue-channel stiffness only.

  R_G       density = sum_{mu nu} G_cd [ (A_mu)^{nu c} (A_nu)^{mu d} - (A_mu)^{mu c} (A_nu)^{nu d} ]
            the derivative index of one jet contracted with a raw (contravariant)
            internal index of the other by delta (the R0 rule for mixed pairs), the
            remaining internal pair contracted with a COVARIANT (0,2) tensor G:
            G = eta;  G = eta M eta;  G = M^{-1};  G = h_cov = eta + 2 (eta u)(eta u)^T
            (u the timelike unit eigenvector of N, u^T eta u = -1).  No eta^{mu nu}
            appears (both derivative indices are absorbed by internal ones), so the
            signature enters only through G.  The mu = nu = 0 terms cancel identically:
            R_G has NO omega^2 content for any G (at most a B-type omega^1 piece that
            drops from the Hamiltonian).  E-density = c_R R_G (sign free; c_R is the
            scanned coefficient).

  K_P       E-density = (1/2) [ sum_i tr((P A_i eta P)^2) + omega^2 tr((P a0 eta P)^2) ]
            P = (N - p1)(N - p2) with (p1, p2) = (lambda_t, 1) = (-g, 1): the spectral
            polynomial that VANISHES on the timelike and on the unit eigenvalue, so in
            the eigenbasis of N the projected jet keeps only the (2,3) block: boosts
            (timelike-spatial mixing) and tilts (axis-(2,3) mixing) are invisible, the
            (2,3) block carries the phase stiffness [f(lambda_2) f(lambda_3)]^2
            (lambda_2 - lambda_3)^2 |d phi|^2 and the split-modulus gradient, f(x) =
            (x - p1)(x - p2).  In the author's own sign convention the vacuum spectrum
            reads (-g, 1, delta, delta) and the polynomial is written (N - g)(N - 1);
            the LITERAL polynomial (p1, p2) = (+g, 1) in our convention does NOT kill the
            timelike direction (f(-g) = 2 g (g + 1) != 0); R14-0 checks both readings and
            the intended one (roots on the two eigenvalues to be killed) is the entrant.
            In the eta-orthonormal eigenbasis the surviving block is spacelike, so both
            static and omega^2 pieces are >= 0 pointwise ON THE VACUUM ORBIT; off the
            spectrum (eigenvalues deformed) the timelike row leaks in with sigma = -1 and
            the static density is INDEFINITE (selftest PROPERTY line; carried to R14-B).

Imported for the LP basis (definitions and hashes stay with their owners):
  K_T (R7, m5_32_r7_a_kt_form.kt_density_np), the C6 quartics and (F.F)^2 (R8,
  m5_32_r8_a_quartics.d_C6a / d_C6b / d_Fpair-type), R1's thirteen (registry +
  m5_32_terms_ext).

Gradients (exact, complex-step gated in the selftests): K_P through the jets AND
through P(M) (P is a polynomial in N, no eigen-decomposition in the gradient);
R_G for G in {eta, eta M eta, M^{-1}} (h_cov: density only, gradient not
implemented, stated); K_lambda: density only (no rung of R14 relaxes under it).

Selftests (python3 m5_32_r14_terms.py): covariance of every entrant under random
boosts and rotations (<= 1e-10) with the no-eta control FAILING; positivity of
K_P; the K_lambda perturbation formula against finite differences of the sorted
spectrum; the complex-step gradient gate for K_P and R_G; the literal-polynomial
mutant of K_P reddening the boost-blindness line.
"""
from __future__ import annotations
import hashlib
import importlib.util
import json
import os
import sys
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.dirname(HERE)
DATA = os.path.join(RES, "data")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


L0 = _load("m5_32_lagrangian", "m5_32_lagrangian.py")
B3 = L0.B3
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])
ETA_D = np.diag(ETA)
EYE4 = np.eye(4)
T0 = time.time()


def log(m):
    print(f"[{time.time() - T0:8.1f}s] {m}", flush=True)


# ============================================================ eigen helpers
def eig_N(M):
    """eigen-decomposition of N = M eta per cell: lam (...,4) real, V (...,4,4)
    columns v_a with v_a^T eta v_a = sig_a in {-1, +1}."""
    N = M @ ETA
    lam, V = np.linalg.eig(N)
    if np.max(np.abs(lam.imag)) > 1e-8 * max(float(np.max(np.abs(lam.real))), 1.0):
        raise ValueError("complex spectrum of M eta")
    lam, V = lam.real, V.real
    n2 = np.einsum("...ak,a,...ak->...k", V, ETA_D, V)
    V = V / np.sqrt(np.abs(n2))[..., None, :]
    return lam, V, np.sign(n2)


def timelike_lambda(M):
    lam, V, sig = eig_N(M)
    k = np.argmin(sig, axis=-1)
    return np.take_along_axis(lam, k[..., None], axis=-1)[..., 0]


# ============================================================ K_lambda
def dlambda_sq(A, M, tol=1e-6):
    """sum_a (d_mu lambda_a)^2 per cell and mu: (4, ...).  Degenerate clusters
    (|lam_a - lam_b| < tol) use tr(B^2) of the cluster block."""
    lam, V, sig = eig_N(M)
    out = np.zeros(A.shape[:-2])
    # B_mu = S V^T eta A_mu eta V  (S = diag(sig)): the eta-projected jet in the eigenbasis
    VeA = np.einsum("...ak,ab,m...bc->m...kc", V, ETA, A)          # (mu,...,k,c)
    B = np.einsum("m...kc,cd,...dl->m...kl", VeA, ETA, V) * sig[None, ..., :, None]
    diag = np.einsum("m...kk->m...k", B)
    out += np.sum(diag ** 2, axis=-1)
    # degenerate correction: for a near-degenerate pair (a,b) replace B_aa^2 + B_bb^2
    # by tr of the 2x2 block squared = B_aa^2 + B_bb^2 + 2 B_ab B_ba
    for a in range(4):
        for b in range(a + 1, 4):
            close = np.abs(lam[..., a] - lam[..., b]) < tol
            if np.any(close):
                corr = 2.0 * B[..., a, b] * B[..., b, a]
                out += np.where(close[None], corr, 0.0)
    return out


def klam_static(A, M):
    """(1/2) sum_i sum_a (d_i lambda_a)^2 per cell from the PERTURBATION formula on the
    jets A[1:4] (the registry-shaped read; on a lattice orbit field the finite-difference
    jet is not tangent to the orbit, so this read carries an O(h^2) residual there; the
    lattice static energy of K_lambda is taken from klam_static_fd below)."""
    d2 = dlambda_sq(A, M)
    return 0.5 * (d2[1] + d2[2] + d2[3])


def sorted_spectrum(M):
    """eigenvalues of N = M eta per cell, sorted DESCENDING (the vacuum reads (1, delta, 0, -g))."""
    lam = np.linalg.eigvals(M @ ETA)
    return -np.sort(-lam.real, axis=-1)


def klam_static_fd(M, cfg):
    """per-cell h^3-weighted static E of K_lambda from the certified finite differences of
    the sorted eigenvalue fields: (1/2) sum_br wt sum_i sum_a (d_i lambda_a)^2.  Exactly zero
    on any lattice field with a cell-wise constant spectrum (every Lorentz-orbit field)."""
    lam = sorted_spectrum(M)
    h3 = cfg["h"] ** 3
    dens = np.zeros(M.shape[:-2])
    for br, wt in B3.branches(cfg["stencil"]):
        for ax in range(3):
            d = B3.d1(lam, ax, cfg["h"], br)
            dens += wt * 0.5 * np.sum(d * d, axis=-1)
    return h3 * dens


def klam_kin(a0, M):
    """(1/2) sum_a (d_t lambda_a)^2 per cell for A_0 = a0 (the omega^2 coefficient)."""
    A = np.zeros((4,) + M.shape)
    A[0] = a0
    return 0.5 * dlambda_sq(A, M)[0]


# ============================================================ R_G
def G_of(kind, M):
    if kind == "eta":
        return np.broadcast_to(ETA, M.shape)
    if kind == "etaMeta":
        return ETA @ M @ ETA
    if kind == "Minv":
        return np.linalg.inv(M)
    if kind == "hcov":
        lam, V, sig = eig_N(M)
        k = np.argmin(sig, axis=-1)
        u = np.take_along_axis(V, k[..., None, None], axis=-1)[..., 0]
        w = u @ ETA
        return ETA + 2.0 * w[..., :, None] * w[..., None, :]
    raise ValueError(kind)


def rg_tensor(A):
    """T^{cd} = sum_{mu nu} [ (A_mu)^{nu c} (A_nu)^{mu d} - (A_mu)^{mu c} (A_nu)^{nu d} ]
    per cell (...,4,4); R_G = G_cd T^{cd}.  A: (4, ..., 4, 4) jets."""
    # X[mu, nu, c] = (A_mu)^{nu c}
    X = np.stack([A[mu][..., nu, :] for mu in range(4) for nu in range(4)], axis=0)
    X = X.reshape((4, 4) + A.shape[1:-2] + (4,))
    first = np.einsum("mn...c,nm...d->...cd", X, X)
    div = sum(X[mu, mu] for mu in range(4))                       # (A_mu)^{mu c} summed
    second = div[..., :, None] * div[..., None, :]
    return first - second


def rg_density(A, M, kind):
    return np.einsum("...cd,...cd->...", G_of(kind, M), rg_tensor(A))


def rg_grad(M, cfg, kind, weight=1.0):
    """gradient of weight * h^3 sum_cells R_G (static, spatial jets) wrt symmetric M,
    for kind in {eta, etaMeta, Minv}."""
    h = cfg["h"]
    Gm = np.zeros_like(M)
    for br, wt in B3.branches(cfg["stencil"]):
        A = np.zeros((4,) + M.shape)
        for ax in range(3):
            A[1 + ax] = B3.d1(M, ax, h, br)
        G = G_of(kind, M)
        # d R_G / d (A_mu)^{nu c}: from the first term 2 G_cd (A_nu)^{mu d} (symmetric in the
        # pair), from the second -2 delta_{mu nu} G_cd div^d
        X = np.stack([A[mu][..., nu, :] for mu in range(4) for nu in range(4)], axis=0)
        X = X.reshape((4, 4) + M.shape[:-2] + (4,))
        div = sum(X[mu, mu] for mu in range(4))
        Gdiv = np.einsum("...cd,...d->...c", G, div)
        dA = np.zeros_like(A)
        for mu in range(4):
            for nu in range(4):
                t = 2.0 * np.einsum("...cd,...d->...c", G, X[nu, mu])
                if mu == nu:
                    t = t - 2.0 * Gdiv
                dA[mu][..., nu, :] += t
        for ax in range(3):
            Gm += wt * B3.d1_adj(dA[1 + ax], ax, h, br)
        # through G(M)
        T = rg_tensor(A)
        if kind == "etaMeta":
            Gm += wt * (ETA @ T @ ETA)
        elif kind == "Minv":
            Mi = np.linalg.inv(M)
            Gm += wt * (-(Mi @ T @ Mi))
        elif kind != "eta":
            raise NotImplementedError(kind)
    return weight * h ** 3 * B3.sym4(Gm)


# ============================================================ K_P
def P_of(M, roots):
    N = M @ ETA
    p1, p2 = roots
    return (N - p1 * EYE4) @ (N - p2 * EYE4)


def kp_block(A_mu, M, roots):
    """X_mu = P A_mu eta P per cell."""
    P = P_of(M, roots)
    return P @ A_mu @ ETA @ P


def kp_static(A, M, roots):
    """(1/2) sum_i tr((P A_i eta P)^2) per cell."""
    P = P_of(M, roots)
    tot = 0.0
    for i in (1, 2, 3):
        X = P @ A[i] @ ETA @ P
        tot = tot + np.einsum("...ab,...ba->...", X, X)
    return 0.5 * tot


def kp_kin(a0, M, roots):
    P = P_of(M, roots)
    X = P @ a0 @ ETA @ P
    return 0.5 * np.einsum("...ab,...ba->...", X, X)


def kp_energy_grad(M, cfg, roots, c=1.0, a0=None, J=None):
    """E = c h^3 sum_br wt sum_cells kp_static, and its exact gradient wrt symmetric
    M (jets and P(M) both chained; complex-step safe).  Returns (E, grad)."""
    h = cfg["h"]
    N = M @ ETA
    p1, p2 = roots
    s = p1 + p2
    P = (N - p1 * EYE4) @ (N - p2 * EYE4)
    E = 0.0
    Gm = np.zeros_like(M)
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, h, br) for ax in range(3)]
        Zsum = np.zeros_like(M)
        for i in range(3):
            Bm = A[i] @ ETA
            X = P @ Bm @ P
            E = E + wt * 0.5 * np.sum(np.einsum("...ab,...ba->...", X, X))
            # d/dA_i of (1/2) tr(X^2): tr(X P dA eta P) -> (eta P X P)^T
            dAi = (ETA @ P @ X @ P).swapaxes(-1, -2)
            Gm += wt * B3.d1_adj(dAi, i, h, br)
            # through P: d (1/2) tr(X^2) = tr(Y dP), Y = B P X + X P B
            Y = Bm @ P @ X + X @ P @ Bm
            Zsum += wt * (N @ Y + Y @ N - s * Y)
        # tr(Z dN) = tr(Z dM eta) = tr(eta Z dM) -> d/dM = (eta Z)^T
        Gm += (ETA @ Zsum).swapaxes(-1, -2)
    return c * h ** 3 * E, c * h ** 3 * B3.sym4(Gm)


# ============================================================ K_P in H-adjoint form (2026-09-05 correction)
def kp_h_blocks(A, M, roots):
    """B_mu in the eta-orthonormal eigenbasis of N (B = S V^T eta A_mu eta V), weighted by
    f_a f_b: Omega'_mu = f_a f_b B_ab.  The H-adjoint form K_P^h = (1/2) eta^{mu nu}
    tr(Omega_mu^T H Omega_nu H^-1), H = eta + 2 (eta u)(eta u)^T (the (0,2) tensor), H^-1 =
    eta + 2 u u^T, is the FROBENIUS norm of Omega' in that basis (H = 1 there): (1/2) sum_ab
    f_a^2 f_b^2 B_ab^2 >= 0 pointwise.  [R14-0 audit, 2026-09-05: the order matters; the
    transposed placement tr(Omega H Omega^T H^-1), as first written here and on the thread, is
    NOT invariant (a factor 20 under a boost); this function evaluates the invariant one.],
    on and off the vacuum spectrum.  The plain form is (1/2) sum_ab f_a^2 f_b^2 sig_a sig_b
    B_ab^2 (B_ba = sig_a sig_b B_ab), which differs exactly on the timelike-spacelike pairs
    that P kills only on-shell."""
    lam, V, sig = eig_N(M)
    p1, p2 = roots
    f = (lam - p1) * (lam - p2)                                      # (..., 4)
    VeA = np.einsum("...ak,ab,m...bc->m...kc", V, ETA, A)
    B = np.einsum("m...kc,cd,...dl->m...kl", VeA, ETA, V) * sig[None, ..., :, None]
    W = f[None, ..., :, None] * f[None, ..., None, :] * B
    return W, sig


def kp_h_static(A, M, roots):
    W, _ = kp_h_blocks(A, M, roots)
    return 0.5 * np.sum(W[1:4] ** 2, axis=(0, -1, -2))


def kp_h_kin(a0, M, roots):
    A = np.zeros((4,) + M.shape)
    A[0] = a0
    W, _ = kp_h_blocks(A, M, roots)
    return 0.5 * np.sum(W[0] ** 2, axis=(-1, -2))


def _kp_h_cell_grad(A_list, M, roots):
    """E = (1/2) sum_A sum_ab f_a^2 f_b^2 X_ab^2 per cell, X_ab = v_a^T eta A eta v_b (v the
    eta-orthonormal eigenvectors of N = M eta, f_a = f(lambda_a)), summed over the matrices in
    A_list (each (..., 4, 4), held FIXED); returns (E_cells, dE/dM_cells at fixed A (through
    lambda and v, first-order perturbation theory), [dE/dA for each A]).  Degenerate pairs
    (|lambda_a - lambda_m| < 1e-8) are skipped in the eigenvector rotation (guard)."""
    lam, V, sig = eig_N(M)
    p1, p2 = roots
    f = (lam - p1) * (lam - p2)
    fp = 2.0 * lam - (p1 + p2)
    w = f[..., :, None] ** 2 * f[..., None, :] ** 2                      # (..., 4, 4)
    eV = ETA @ V                                                        # columns eta v_a
    E = np.zeros(M.shape[:-2])
    dEdf = np.zeros(lam.shape)
    dEdv = np.zeros(V.shape)                                            # (..., 4, a): dE/dv_a as column a
    dA_out = []
    for A in A_list:
        X = np.einsum("...ca,...cd,...db->...ab", eV, A, eV)            # X_ab = (eta v_a)^T A (eta v_b)
        E += 0.5 * np.sum(w * X * X, axis=(-1, -2))
        WX = w * X
        # dE/dA: d(1/2 sum w X^2) = sum w X dX, dX_ab = (eta v_a)^T dA (eta v_b) -> dE/dA = eV (WX) eV^T
        dA_out.append(np.einsum("...ca,...ab,...db->...cd", eV, WX, eV))
        # dE/df_a = 2 f_a sum_b f_b^2 (X_ab^2 + X_ba^2)
        dEdf += f * np.einsum("...ab,...b->...a", X * X + (X * X).swapaxes(-1, -2), f * f)
        # dE/dv_a = sum_b w_ab X_ab (eta A eta v_b) + sum_b w_ba X_ba (eta A^T eta v_b) = 2 sum_b w_ab X_ab eta A eta v_b (A symmetric, w symmetric)
        EAeV = np.einsum("ab,...bc,...cd->...ad", ETA, A, eV)           # columns eta A eta v_b
        dEdv += 2.0 * np.einsum("...ab,...cb->...ca", WX, EAeV)
    # chain to M: dlambda_a = sig_a (eta v_a)^T dM (eta v_a); dv_a = sum_{m != a} v_m sig_m (eta v_m)^T dM (eta v_a) / (lambda_a - lambda_m)
    G = np.zeros(M.shape)
    coef_l = dEdf * fp * sig                                            # (..., 4)
    G += np.einsum("...a,...ca,...da->...cd", coef_l, eV, eV)
    g_am = np.einsum("...ca,...cm->...am", dEdv, V)                     # (dE/dv_a) . v_m
    den = lam[..., :, None] - lam[..., None, :]                         # lambda_a - lambda_m
    safe = np.abs(den) > 1e-8
    coef = np.where(safe, g_am * sig[..., None, :] / np.where(safe, den, 1.0), 0.0)
    np.einsum("...aa->...a", coef)[...] = 0.0
    G += np.einsum("...am,...cm,...da->...cd", coef, eV, eV)
    G = 0.5 * (G + G.swapaxes(-1, -2))
    return E, G, dA_out


def kp_h_energy_grad(M, cfg, roots, a0=None):
    """K_P^h on the certified stencil: returns (E_stat, G_stat, kin, G_kin) where
    E_stat = h^3 sum_br wt sum_cells (1/2) sum_i ||Om_i||_F^2 (eigenbasis) and G_stat its exact
    gradient wrt symmetric M (jets AND eigenbasis chained), kin = h^3 sum (1/2) ||Om_0||^2 for
    A_0 = a0 (the omega^2 coefficient) and G_kin its gradient at FROZEN a0."""
    h = cfg["h"]
    h3 = h ** 3
    E_stat, G_stat = 0.0, np.zeros_like(M)
    for br, wt in B3.branches(cfg["stencil"]):
        A = [B3.d1(M, ax, h, br) for ax in range(3)]
        Ec, Gc, dA = _kp_h_cell_grad(A, M, roots)
        E_stat += wt * float(np.sum(Ec))
        G_stat += wt * Gc
        for ax in range(3):
            G_stat += wt * B3.d1_adj(dA[ax], ax, h, br)
    kin, G_kin = None, None
    if a0 is not None:
        Ek, Gk, _ = _kp_h_cell_grad([a0], M, roots)
        kin, G_kin = h3 * float(np.sum(Ek)), h3 * Gk
    return h3 * E_stat, h3 * G_stat, kin, G_kin


def klam_energy_grad(M, cfg):
    """K_lambda static: E = h^3 sum_br wt sum_i sum_cells (1/2) sum_a (D_i lambda_a)^2 on the sorted
    spectrum fields (klam_static_fd), and its exact gradient: dE/dlambda_a = h^3 sum_br wt sum_i
    D_i^T (D_i lambda_a) (the stencil adjoint), dlambda_a/dM = sig_a (eta v_a)(eta v_a)^T with v_a
    the eigenvector of the a-th SORTED eigenvalue (first-order perturbation theory)."""
    h = cfg["h"]
    lam_s = sorted_spectrum(M)                                    # (..., 4) descending
    E = 0.0
    dEdl = np.zeros(lam_s.shape)
    for br, wt in B3.branches(cfg["stencil"]):
        for ax in range(3):
            d = B3.d1(lam_s, ax, h, br)
            E += wt * 0.5 * float(np.sum(d * d))
            dEdl += wt * B3.d1_adj(d, ax, h, br)
    # eigenvectors matched to the sorted order
    lam, V, sig = eig_N(M)
    order = np.argsort(-lam, axis=-1)
    Vs = np.take_along_axis(V, order[..., None, :], axis=-1)
    sigs = np.take_along_axis(sig, order, axis=-1)
    eV = ETA @ Vs
    G = np.einsum("...a,...a,...ca,...da->...cd", dEdl, sigs, eV, eV)
    G = 0.5 * (G + G.swapaxes(-1, -2))
    return h ** 3 * E, h ** 3 * G


def rg_hcov_energy_grad(M, cfg):
    """R_hcov static: E = h^3 sum_br wt sum_cells G_cd T^cd with G = eta + 2 w w^T, w = eta u (u the
    timelike unit eigenvector, u^T eta u = -1), and its exact gradient through the jets (as rg_grad)
    and through u(M): dE/du = 4 eta T_sym w, du = sum_{m != t} v_m sig_m (eta v_m)^T dM (eta u) /
    (lambda_t - lambda_m)."""
    h = cfg["h"]
    lam, V, sig = eig_N(M)
    kt = np.argmin(sig, axis=-1)
    u = np.take_along_axis(V, kt[..., None, None], axis=-1)[..., 0]
    w = u @ ETA
    G_ten = ETA + 2.0 * w[..., :, None] * w[..., None, :]
    E = 0.0
    Gm = np.zeros_like(M)
    Tsum = np.zeros_like(M)
    for br, wt in B3.branches(cfg["stencil"]):
        A = np.zeros((4,) + M.shape)
        for ax in range(3):
            A[1 + ax] = B3.d1(M, ax, h, br)
        T = rg_tensor(A)
        E += wt * float(np.sum(np.einsum("...cd,...cd->...", G_ten, T)))
        Tsum += wt * T
        X = np.stack([A[mu][..., nu, :] for mu in range(4) for nu in range(4)], axis=0).reshape((4, 4) + M.shape[:-2] + (4,))
        div = sum(X[mu, mu] for mu in range(4))
        Gdiv = np.einsum("...cd,...d->...c", G_ten, div)
        dA = np.zeros_like(A)
        for mu in range(4):
            for nu in range(4):
                t = 2.0 * np.einsum("...cd,...d->...c", G_ten, X[nu, mu])
                if mu == nu:
                    t = t - 2.0 * Gdiv
                dA[mu][..., nu, :] += t
        for ax in range(3):
            Gm += wt * B3.d1_adj(dA[1 + ax], ax, h, br)
    # through u: dE/dw = 4 T_sym w, w = eta u -> dE/du = eta (4 T_sym w)
    Tsym = 0.5 * (Tsum + Tsum.swapaxes(-1, -2))
    dEdu = np.einsum("ab,...bc,...c->...a", ETA, 4.0 * Tsym, w)
    g_m = np.einsum("...a,...am->...m", dEdu, V)                     # dE/du . v_m
    lam_t = np.take_along_axis(lam, kt[..., None], axis=-1)
    den = lam_t - lam                                               # lambda_t - lambda_m
    safe = np.abs(den) > 1e-8
    coef = np.where(safe, g_m * sig / np.where(safe, den, 1.0), 0.0)
    np.put_along_axis(coef, kt[..., None], 0.0, axis=-1)
    eV = ETA @ V
    eu = u @ ETA
    Gm += np.einsum("...m,...cm,...d->...cd", coef, eV, eu)
    return h ** 3 * E, h ** 3 * B3.sym4(Gm)


def kp_plain_from_blocks(A, M, roots):
    """the plain-trace K_P recomputed from the eigenbasis blocks (a cross-check of kp_static)."""
    W, sig = kp_h_blocks(A, M, roots)
    ss = sig[..., :, None] * sig[..., None, :]
    return 0.5 * np.sum(W[1:4] ** 2 * ss[None], axis=(0, -1, -2))


def kp_plain_kin_grad(M, a0, cfg, roots):
    """the plain K_P omega^2 coefficient kin = h^3 sum (1/2) tr((P a0 eta P)^2) and its gradient wrt
    M at FROZEN a0 through P(M) (polynomial route, complex-step safe)."""
    h3 = cfg["h"] ** 3
    N = M @ ETA
    p1, p2 = roots
    s = p1 + p2
    P = (N - p1 * EYE4) @ (N - p2 * EYE4)
    Bm = a0 @ ETA
    X = P @ Bm @ P
    kin = 0.5 * float(np.sum(np.einsum("...ab,...ba->...", X, X)))
    Y = Bm @ P @ X + X @ P @ Bm
    Z = N @ Y + Y @ N - s * Y
    G = (ETA @ Z).swapaxes(-1, -2)
    return h3 * kin, h3 * B3.sym4(G)


# ============================================================ registry shape
class Entrant:
    def __init__(self, name, definition, static_fn, kin_fn, grad_fn=None):
        self.name = name
        self.definition = definition
        self.hash = hashlib.sha256(definition.encode()).hexdigest()[:12]
        self.static_fn = static_fn        # (A, M) -> per-cell static E-density
        self.kin_fn = kin_fn              # (a0, M) -> per-cell omega^2 E-coefficient
        self.grad_fn = grad_fn            # (M, cfg) -> (E, grad) or None


ROOTS_INTENDED = (-8.0, 1.0)       # (lambda_t, 1) at g = 8 in our convention
ROOTS_LITERAL = (8.0, 1.0)         # the author's formula read literally in our convention


def roots_for(g):
    return (-float(g), 1.0)


ENTRANTS = {}


def _reg(e):
    ENTRANTS[e.name] = e
    return e


_reg(Entrant("K_lambda",
             "K_lambda: E = (1/2) sum_a [sum_i (d_i lambda_a)^2 + omega^2 (d_t lambda_a)^2], "
             "lambda_a eigenvalues of N = M eta, d_mu lambda_a = (v_a^T eta A_mu eta v_a)/(v_a^T eta v_a)",
             klam_static, klam_kin))
for _k in ("eta", "etaMeta", "Minv", "hcov"):
    _reg(Entrant(f"R_{_k}",
                 f"R_G with G = {_k}: sum_mu nu G_cd [(A_mu)^(nu c)(A_nu)^(mu d) - (A_mu)^(mu c)(A_nu)^(nu d)], "
                 "mixed pairs by delta, G covariant (0,2); E-density = c_R R_G",
                 (lambda A, M, k=_k: rg_density(A, M, k)),
                 (lambda a0, M, k=_k: np.zeros(M.shape[:-2])),
                 (lambda M, cfg, k=_k: (None, rg_grad(M, cfg, k))) if _k != "hcov" else None))
_reg(Entrant("K_P_h",
             "K_P^h (H-adjoint form, the author's 2026-09-05 correction, invariant order): E = (1/2)[sum_i tr(Om_i^T H Om_i H^-1) "
             "+ omega^2 tr(Om_0^T H Om_0 H^-1)], Om_mu = P A_mu eta P, H = eta + 2 (eta u)(eta u)^T, H^-1 = eta + 2 u u^T; "
             "= the Frobenius norm of the projected jet in the eta-orthonormal eigenbasis (PSD everywhere)",
             (lambda A, M: kp_h_static(A, M, roots_for(-timelike_lambda_vac(M)))),
             (lambda a0, M: kp_h_kin(a0, M, roots_for(-timelike_lambda_vac(M))))))
_reg(Entrant("K_P",
             "K_P: E = (1/2)[sum_i tr((P A_i eta P)^2) + omega^2 tr((P a0 eta P)^2)], "
             "N = M eta, P = (N - lambda_t)(N - 1), lambda_t the timelike eigenvalue of the vacuum (-g)",
             (lambda A, M: kp_static(A, M, roots_for(-timelike_lambda_vac(M)))),
             (lambda a0, M: kp_kin(a0, M, roots_for(-timelike_lambda_vac(M)))),
             (lambda M, cfg: kp_energy_grad(M, cfg, roots_for(cfg["g"])))))


def timelike_lambda_vac(M):
    """the timelike eigenvalue is taken from the field's own corner cell (the pinned
    vacuum), so the roots follow g without a cfg."""
    Mc = M.reshape((-1, 4, 4))[0]
    return float(timelike_lambda(Mc[None])[0])


# ============================================================ lattice reads
def jets_static(M, cfg, br):
    A = np.zeros((4,) + M.shape)
    for ax in range(3):
        A[1 + ax] = B3.d1(M, ax, cfg["h"], br)
    return A


def static_energy(name, M, cfg):
    """h^3-weighted static E of an entrant on the certified sym stencil."""
    if name == "K_lambda":
        return float(np.sum(klam_static_fd(M, cfg)))
    e = ENTRANTS[name]
    tot = 0.0
    for br, wt in B3.branches(cfg["stencil"]):
        tot = tot + wt * np.sum(e.static_fn(jets_static(M, cfg, br), M))
    return cfg["h"] ** 3 * float(tot)


def static_density(name, M, cfg):
    if name == "K_lambda":
        return klam_static_fd(M, cfg)
    e = ENTRANTS[name]
    d = 0.0
    for br, wt in B3.branches(cfg["stencil"]):
        d = d + wt * e.static_fn(jets_static(M, cfg, br), M)
    return cfg["h"] ** 3 * d


def kin_energy(name, M, a0, cfg):
    """h^3-weighted omega^2 coefficient (E_kin = kin omega^2)."""
    e = ENTRANTS[name]
    return cfg["h"] ** 3 * float(np.sum(e.kin_fn(a0, M)))


def kin_density(name, M, a0, cfg):
    return cfg["h"] ** 3 * ENTRANTS[name].kin_fn(a0, M)


# ============================================================ selftests
def _transform(L, A, M):
    Linv_T = np.linalg.inv(L).T
    Mp = np.einsum("ab,...bc,dc->...ad", L, M, L)
    Ap = np.einsum("mn,n...ab->m...ab", Linv_T,
                   np.einsum("ab,n...bc,dc->n...ad", L, A, L))
    return Ap, Mp


def _random_field(rng, npts, g=8.0, amp=0.3):
    """random symmetric M near the vacuum (real spectrum) and random symmetric jets."""
    M = np.zeros((npts, 4, 4))
    M[:] = np.diag([g, 1.0, 0.3, 0.0])
    dM = rng.normal(size=(npts, 4, 4)) * amp
    M += 0.5 * (dM + dM.swapaxes(-1, -2))
    A = rng.normal(size=(4, npts, 4, 4))
    A = 0.5 * (A + A.swapaxes(-1, -2))
    return A, M


def st_covariance(res, lines, rng):
    A, M = _random_field(rng, 40)
    worst = {}
    for kind in ("boost", "rotation"):
        for k in range(3):
            Lm = L0._lorentz(rng, kind)
            Ap, Mp = _transform(Lm, A, M)
            for nm, e in ENTRANTS.items():
                # only the FULL eta-contracted density is a scalar: L = static(A_i) - kin(A_0)
                # (a boost mixes A_0 into the spatial jets); R_G carries all four jets already
                d0 = e.static_fn(A, M) - e.kin_fn(A[0], M)
                d1 = e.static_fn(Ap, Mp) - e.kin_fn(Ap[0], Mp)
                drift = float(np.max(np.abs(d1 - d0)) / max(float(np.max(np.abs(d0))), 1e-300))
                worst[nm] = max(worst.get(nm, 0.0), drift)
    res["covariance"] = worst
    w = max(worst.values())
    lines.append(("covariance: every entrant (static + omega^2 densities) under 3 boosts + 3 rotations",
                  w, 1e-10, w <= 1e-10))
    # the no-eta control: K_P with N = M (Euclidean) must FAIL under a boost
    Lm = L0._lorentz(rng, "boost")
    Ap, Mp = _transform(Lm, A, M)

    def kp_eucl(A_, M_):
        P = (M_ - ROOTS_INTENDED[0] * EYE4) @ (M_ - ROOTS_INTENDED[1] * EYE4)
        return sum(np.einsum("...ab,...ba->...", P @ A_[i] @ P, P @ A_[i] @ P) for i in (1, 2, 3))
    d0, d1 = kp_eucl(A, M), kp_eucl(Ap, Mp)
    ctrl = float(np.max(np.abs(d1 - d0)) / max(float(np.max(np.abs(d0))), 1e-300))
    res["covariance_control_noeta"] = ctrl
    lines.append(("no-eta control (K_P built on M instead of M eta) FAILS covariance under a boost",
                  ctrl, 1e-3, ctrl > 1e-3))


def st_positivity(res, lines, rng):
    # (a) ON the vacuum orbit (spectrum exactly (-g, 1, delta, 0), random Lorentz frames,
    #     random symmetric jets): the projected block is spacelike, so K_P >= 0 pointwise
    npts = 200
    M0 = np.diag([8.0, 1.0, 0.3, 0.0])
    M = np.zeros((npts, 4, 4)); A = np.zeros((4, npts, 4, 4))
    for k in range(npts):
        Lm = L0._lorentz(rng, "boost" if k % 2 else "rotation", scale=0.8)
        M[k] = Lm @ M0 @ Lm.T
        a = rng.normal(size=(4, 4, 4)); A[:, k] = 0.5 * (a + a.swapaxes(-1, -2))
    ks = kp_static(A, M, ROOTS_INTENDED)
    kk = kp_kin(A[0], M, ROOTS_INTENDED)
    kl = klam_static(A, M)
    mn = float(min(np.min(ks), np.min(kk), np.min(kl)))
    res["positivity_min_on_orbit"] = mn
    lines.append(("positivity ON the vacuum orbit: K_P static, K_P omega^2, K_lambda static >= 0 (200 cells)",
                  mn, -1e-12, mn >= -1e-12))
    # (b) OFF the spectrum (eigenvalues deformed, amp 0.3): the timelike row is no longer
    #     killed and enters with sigma = -1: K_P is INDEFINITE off-shell (a property carried
    #     to R14-B's boundedness probes), K_lambda stays >= 0 (a sum of squares)
    A2, M2 = _random_field(rng, 200)
    ks2 = kp_static(A2, M2, ROOTS_INTENDED)
    kl2 = klam_static(A2, M2)
    res["positivity_off_shell"] = {"K_P_min": float(np.min(ks2)), "K_P_frac_negative": float(np.mean(ks2 < 0)),
                                   "K_lambda_min": float(np.min(kl2))}
    lines.append(("PROPERTY: K_P static density is INDEFINITE off the vacuum spectrum (min < 0 on random cells)",
                  float(np.min(ks2)), 0.0, float(np.min(ks2)) < 0.0))
    lines.append(("K_lambda static density >= 0 off-shell as well (a sum of squares)",
                  float(np.min(kl2)), -1e-12, float(np.min(kl2)) >= -1e-12))
    # (c) the H-adjoint form: PSD on 2000 random off-shell probes; equal to the plain form on
    #     the orbit; the plain form recomputed from the eigenbasis blocks equals kp_static
    A3, M3 = _random_field(rng, 2000)
    kh = kp_h_static(A3, M3, ROOTS_INTENDED)
    res["K_P_h_off_shell"] = {"min": float(np.min(kh)), "frac_negative": float(np.mean(kh < 0)),
                              "plain_frac_negative_same_probes": float(np.mean(kp_static(A3, M3, ROOTS_INTENDED) < 0))}
    lines.append(("K_P^h (H-adjoint) static density >= 0 on 2000 random OFF-shell probes (the plain form's negative fraction on the same probes reported)",
                  float(np.min(kh)), -1e-12, float(np.min(kh)) >= -1e-12))
    d_orbit = float(np.max(np.abs(kp_h_static(A, M, ROOTS_INTENDED) - ks)) / max(float(np.max(np.abs(ks))), 1e-300))
    lines.append(("K_P^h == K_P on the vacuum orbit (rel)", d_orbit, 1e-10, d_orbit <= 1e-10))
    d_plain = float(np.max(np.abs(kp_plain_from_blocks(A2, M2, ROOTS_INTENDED) - ks2)) / max(float(np.max(np.abs(ks2))), 1e-300))
    lines.append(("plain K_P recomputed from the eigenbasis blocks (sig_a sig_b weights) == kp_static off-shell (rel)", d_plain, 1e-10, d_plain <= 1e-10))


def st_klambda_fd(res, lines, rng):
    A, M = _random_field(rng, 30, amp=0.1)
    worst = 0.0
    for mu in range(4):
        eps = 1e-6
        lp = np.sort(np.linalg.eigvals((M + eps * A[mu]) @ ETA).real, axis=-1)
        lm = np.sort(np.linalg.eigvals((M - eps * A[mu]) @ ETA).real, axis=-1)
        fd = np.sum(((lp - lm) / (2 * eps)) ** 2, axis=-1)
        an = dlambda_sq(A, M)[mu]
        worst = max(worst, float(np.max(np.abs(fd - an) / np.maximum(np.abs(fd), 1e-8))))
    res["klambda_fd_rel"] = worst
    lines.append(("K_lambda: perturbation formula vs finite differences of the sorted spectrum (rel)",
                  worst, 1e-5, worst <= 1e-5))


def st_gradients(res, lines, rng):
    cfg = B3.base_cfg(s=-1.0, g=8.0, n=8, L=12.0)
    M = np.zeros((8, 8, 8, 4, 4)) + np.diag([8.0, 1.0, 0.3, 0.0])
    dM = rng.normal(size=M.shape) * 0.2
    M += 0.5 * (dM + dM.swapaxes(-1, -2))
    out = {}
    for nm in ("K_P", "R_etaMeta", "R_Minv", "R_eta"):
        e = ENTRANTS[nm]
        if nm == "K_P":
            E0, G = kp_energy_grad(M, cfg, roots_for(8.0))
            def Ef(Mx):
                return kp_energy_grad(Mx, cfg, roots_for(8.0))[0]
        else:
            kind = nm[2:]
            G = rg_grad(M, cfg, kind)
            def Ef(Mx, kind=kind):
                tot = 0.0
                for br, wt in B3.branches(cfg["stencil"]):
                    Ax = np.zeros((4,) + Mx.shape, dtype=Mx.dtype)
                    for ax in range(3):
                        Ax[1 + ax] = B3.d1(Mx, ax, cfg["h"], br)
                    tot = tot + wt * np.sum(rg_density_cs(Ax, Mx, kind))
                return cfg["h"] ** 3 * tot
        worst = 0.0
        for _ in range(4):
            D = rng.normal(size=M.shape)
            D = 0.5 * (D + D.swapaxes(-1, -2))
            hcs = 1e-20
            Ec = Ef(M + 1j * hcs * D)
            dcs = float(np.imag(Ec) / hcs)
            dan = float(np.sum(G * D))
            worst = max(worst, abs(dcs - dan) / max(abs(dcs), 1e-300))
        out[nm] = worst
        lines.append((f"complex-step gradient gate: {nm} (rel, 4 random directions, n = 8)",
                      worst, 1e-9, worst <= 1e-9))
    res["gradient_gates"] = out


def rg_density_cs(A, M, kind):
    """complex-safe R_G density (no eig): kinds eta / etaMeta / Minv."""
    if kind == "eta":
        G = np.broadcast_to(ETA, M.shape)
    elif kind == "etaMeta":
        G = ETA @ M @ ETA
    elif kind == "Minv":
        G = np.linalg.inv(M)
    else:
        raise ValueError(kind)
    return np.einsum("...cd,...cd->...", G, rg_tensor(A))


def st_channel_blindness(res, lines):
    """K_P at the vacuum: nonzero on the (2,3) clock only; the literal roots see the boosts."""
    M = np.diag([8.0, 1.0, 0.3, 0.0])[None]
    gens = {}
    Jx = np.zeros((4, 4)); Jx[2, 3], Jx[3, 2] = -1.0, 1.0          # rot in (2,3): the clock
    Jy = np.zeros((4, 4)); Jy[1, 3], Jy[3, 1] = 1.0, -1.0           # rot in (1,3): tilt
    Jz = np.zeros((4, 4)); Jz[1, 2], Jz[2, 1] = -1.0, 1.0           # rot in (1,2): tilt
    for i, nm in enumerate(("boost_1", "boost_2", "boost_3")):
        K = np.zeros((4, 4)); K[0, 1 + i] = K[1 + i, 0] = 1.0
        gens[nm] = K
    gens.update({"rot_23_clock": Jx, "rot_13_tilt": Jy, "rot_12_tilt": Jz})
    tab = {}
    for nm, X in gens.items():
        a0 = (X @ M + M @ X.T)                                        # the tangent X M + M X^T
        tab[nm] = {"intended": float(kp_kin(a0, M, ROOTS_INTENDED)[0]),
                   "literal": float(kp_kin(a0, M, ROOTS_LITERAL)[0]),
                   "K_lambda": float(klam_kin(a0, M)[0])}
    res["vacuum_channels"] = tab
    blind = all(abs(tab[k]["intended"]) < 1e-12 for k in tab if k != "rot_23_clock")
    clock = tab["rot_23_clock"]["intended"]
    f2, f3 = (0.3 + 8.0) * (0.3 - 1.0), (0.0 + 8.0) * (0.0 - 1.0)
    pred = 0.5 * 2.0 * (f2 * f3) ** 2 * (0.3 - 0.0) ** 2             # (1/2) tr(X^2), X 2x2 antisym-like block
    lines.append(("K_P (intended roots) at the vacuum: zero on 3 boosts + 2 tilts, nonzero on the (2,3) clock",
                  clock, 1e-12, blind and clock > 1e-12))
    lines.append(("K_P clock stiffness at the vacuum = [f(delta) f(0)]^2 (delta - 0)^2 (analytic)",
                  abs(clock - pred) / pred, 1e-12, abs(clock - pred) / pred <= 1e-12))
    lit_boost = max(abs(tab[k]["literal"]) for k in ("boost_1", "boost_2", "boost_3"))
    lines.append(("MUTANT: the literal polynomial (N - g)(N - 1) in our convention is NOT blind to boosts",
                  lit_boost, 1e-12, lit_boost > 1e-12))
    kl = max(abs(v["K_lambda"]) for v in tab.values())
    lines.append(("K_lambda omega^2 coefficient is exactly 0 on all six vacuum generator channels",
                  kl, 1e-14, kl <= 1e-14))


def selftest(write=True):
    rng = np.random.default_rng(14140)
    res, lines = {"entrants": {k: {"definition": e.definition, "hash": e.hash} for k, e in ENTRANTS.items()}}, []
    st_covariance(res, lines, rng)
    st_positivity(res, lines, rng)
    st_klambda_fd(res, lines, rng)
    st_gradients(res, lines, rng)
    st_channel_blindness(res, lines)
    npass = sum(1 for l in lines if l[3])
    for name, val, thr, ok in lines:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}: {val:.3e} (thr {thr:.1e})")
    res["lines"] = [{"check": n, "value": v, "threshold": t, "pass": bool(ok)} for n, v, t, ok in lines]
    res["summary"] = {"pass": npass, "total": len(lines)}
    log(f"selftest {npass}/{len(lines)}")
    if write:
        json.dump(res, open(os.path.join(DATA, "m5_32_r14_terms_selftest.json"), "w"), indent=1)
    return res


if __name__ == "__main__":
    selftest()
