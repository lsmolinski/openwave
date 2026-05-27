# M5.5 + M5.6 — Paper Lagrangian, KG Emergence & Faber Regularization (math reference)

**Purpose:** the confirmed mathematical foundation for **M5.5** (the Eq.18 action) and **M5.6** (KG-from-twist emergence). §1–4: Duda's Eq.18 action, the building-block operators, the Eq.35 Euler–Lagrange evolution of the matrix field `M`, the matrix Hamiltonian, the `V(M)` options, and the transcription of Duda's Mathematica source (Fig.9) reducing the twist equation to the hedgehog Klein–Gordon — prototyped in `sandbox_v5`. §5 + §5a–§5f: the **M5.6 findings** — the KG mass is *geometric* (minimal coupling to the hedgehog connection `Â`, M5.6.1), the biaxial hedgehog's curvature `C_μν~1/r²` sources it dynamically (M5.6.2), Faber's `Λ=q₀⁶/r₀⁴` regularization pins the mass scale `E₀∝1/r₀` (M5.6.3), the EM/tilt sector reproduces Maxwell by both routes (M5.6.4), the biaxial seeder is ported to production behind an analytic eigensolver fix (M5.6.5a, §5e), turning V on confines the amplitude via a `b=0` well — the 3-term Eq.13 has no biaxial minimum (M5.6.5c, §5f), and the faithful Eq.18 kinetic differs from the shipped `½‖Ṁ‖²` only in physical-mode inertia (the twist/clock frequency, for M5.8) — not gauge slosh (M5.6.5d, §5g).

**Source:** Duda, *Framework for liquid crystal based particle models* (arxiv:2108.07896 v7), §II–IV + Fig.9 (math reading **confirmed by Rodrigo 2026-05-26**); Faber & Golubich, *Universe* 11/2025/113 (regularization, §5c).

**Sister docs:** [4a_convo_2026.05.12.md](4a_convo_2026.05.12.md) (paper digest, eigenvalue map), [0b_M5_roadmap.md § M5.5–M5.6](0b_M5_roadmap.md), [1a_lagrangian_framework.md](1a_lagrangian_framework.md), [[reference_faber_regularization]].

---

## 1. The action (Eq.18) — a Maxwell-analog `E² − B² − V` on the matrix field

```text
ℒ = Σ_{μ=1..3} ‖F_{μ0}‖²_F  −  Σ_{1≤μ<ν≤3} ‖F_{μν}‖²_F  −  V(M)
        time/"electric" curvatures (+)    spatial/"magnetic" curvatures (−)
```

`‖·‖_F` is the Frobenius norm, `A•B = Tr(AB^T)`, `‖A‖_F = √(A•A)` (Eq.17). In 4D the
Frobenius norm picks up the signature `ξ = diag(−1,1,1,1)`: `‖X‖²_ξ = Tr(XξX^Tξ)` (Eq.41) —
that is the M5.8 extension, NOT M5.5 (M5.5 is the 3D case).

## 2. Building blocks (Eq.19–20) — why the EL equation is higher-order

```text
M_μ  := ∂_μ M                                              matrix derivative (have it)
A_μ  := [M, M_μ] = M·∂_μM − ∂_μM·M                          Eq.19 — the 4-potential analogue
F_μν := ∂_μ A_ν − ∂_ν A_μ = 2(∂_μM·∂_νM − ∂_νM·∂_μM)        Eq.20 — Maurer–Cartan, =: 2F_μν
```

Because `A_μ` is itself a commutator of `M` with `∂_μM`, the curvature `F` carries `(∂M)²`
terms and the action carries `(∂M)⁴` — the EL equation is correspondingly higher-order.
This is the core reason we prototype in sympy/numpy before Taichi.

`F_{μ0}` (time index) are the "electric" curvatures `E_i`; `F_{μν}` (spatial) are the
"magnetic" curvatures `B_i`. In the rotation-generator basis the curvatures are 3-vectors
`R = F*` (dual tensor), split by energy scale into `B^1/E^1` (high-energy EM tilts),
`B^2/E^2`, `B^3/E^3` (low-energy QM twist) — see [4a §8](4a_convo_2026.05.12.md).

## 3. Potential `V(M)` (Eq.12 → Eq.13) — the Q7-open part, simple→graduate

```text
Eq.12 (start):  V(M) = Σ_i (λ_i − Λ_i)²            eigenvalue-preference; Λ = (1, δ, 0)
Eq.13 (LdG):    V_LG(M) = a·Tr(M²) − b·Tr(M³) + c·(Tr(M²))²
```

`V` defines `D = diag(Λ_i)` as the vacuum shape AND regularizes the defect-core singularity
(lets the field deviate from `D` near the core → finite energy, the running-coupling effect).
Exact form is Duda's own **open research question (Q7)** — we start from Eq.12, graduate to
Eq.13, then port Faber's regularization ([reference_faber_regularization](../../../../..)).

## 4. Euler–Lagrange evolution (Eq.35) — the equation M5.5 implements

Varying `M = ODO^T` against the SO(3) generators `G` (zeroing `δℒ`) gives, in the vacuum
case (`V=0`, fixed `D`), with `Γ_μ = O^T∂_μO`, `M̄_μ=[Γ_μ,D]`, `F̄_μν=[M̄_μ,M̄_ν]`,
`G' = [G,D]`:

```text
0 = Σ_{μν} d_{μν} Tr( F̄_μν·([Γ_ν,[M̄_μ,G']] − [Γ_μ,[M̄_ν,G']]) )
        + Tr( F̄_{μν,μ}·[M̄_ν,G'] − F̄_{μν,ν}·[M̄_μ,G'] )           (Eq.35)

    d_{μν} = +1 for μν ∈ {10, 20, 30}   (electric/time)
             −1 for μν ∈ {23, 31, 12}   (magnetic/spatial)
              0 otherwise
```

For the 3 rotation generators this collapses (Eq.37–38) to:

```text
twist  (generator Gx):   Γ²·Γ³ = Γ³·Γ²        →  Klein–Gordon for the phase ψ
tilt1  (Gy):             X¹·Γ³ = 0             →  Maxwell
tilt2  (Gz):             X¹·Γ² = 0             →  Maxwell
   with  X^i = (−∇·B^i,  ∂_0 B^i + ∇×E^i)     (Eq.37)
```

The general case adds a `G'_μ = ∂_μG' = GD_μ + D_μG^T` term (Eq.28, from `∂_μD = D_μ`) and
the `V` variation `∂L/∂(λ_i)` (Higgs-like). M5.5 starts in the vacuum/fixed-D case, then adds V.

## 5. Hedgehog reduction → Klein–Gordon (Fig.9, the M5.6 anchor + our M5.5.1 validation)

For the hedgehog ansatz `Q0 = exp(θ Gz)·exp(φ Gy)·exp(ψ Gx)` with `φ=atan2(y,x)`,
`θ=atan2(√(x²+y²), z)`, and phase `ψ = ψ(t,x,y,z)`, the twist equation reduces to:

```text
2 ∂_tt ψ = [ (∇ − Â^hedg)²  +  (Â^hedg·∇ / ‖Â^hedg‖)² ] ψ          Klein–Gordon-like

    Â^hedg(x,y,z) = (x, y, z) / r²,     r = √(x²+y²+z²)
    dual: Ψ = exp(iψ),  Klein–Gordon clock at E = mc² → ω = mc²/ℏ (Zitterbewegung)
```

The other two (tilt) equations are satisfied identically (`= 0`) → they are the Maxwell
sector. **This is the M5.5.1 validation target** (reproduce it from Eq.35) and the **M5.6
implementation target** (KG emerges from twist, not added by hand).

### 5a. M5.6.1 findings (`sandbox_v6`, 2026-05-27) — the KG mass is GEOMETRIC, not a potential

Anatomy of the §5 operator, verified symbolically (`m5_6_1_kg_operator_check.py`) and
numerically (`m5_6_1b_twist_evolution.py`):

| Finding | Detail |
| --- | --- |
| **Hedgehog connection identities** | `∇·Â = 1/r²` and `‖Â‖² = 1/r²` (sympy). |
| **Explicit mass term CANCELS** | `(∇−Â)²ψ = ∇²ψ − (∇·Â)ψ − 2Â·∇ψ + ‖Â‖²ψ`; the zeroth-order coeff `−(∇·Â)+‖Â‖² = 0` **exactly**. The full operator reduces to `L = 2∂_rr + (1/r²)Δ_Ω` (checked on `r⁴`, `cos r`, `e^{−r}`). **Bare phase ψ is MASSLESS** (numeric: residual `→0` as dx→0). |
| **Mass is geometric, lives in `Ψ=e^{iψ}`** | The KG mass is minimal coupling to the connection `Â` (the `k→k−Â` shift), **NOT** an added `V_ψ`. `‖Â‖²=1/r²` survives as a position-dependent mass² only in the dual complex field via `(∇−iÂ)²Ψ`. This is the literal statement of *"KG from twist, not from V"* — and pinpoints the **M5.2 error** (we added mass to a potential; the mass was always in the geometry). |
| **Regularization GENERATES the finite mass** | With a core-regularized `Â = x/(r²+r_c²)`, the cancellation no longer holds: `−(∇·Â)+‖Â‖² = −3r_c²/(r²+r_c²)²` (numeric matches analytic to 5%). ⇒ emergent **mass²(r) = 3r_c²/(2(r²+r_c²)²)** — scale `~1/r_c²` at the core, `→0` far out. The particle mass is *set by the core regularization* (why Faber/M5.6.3 is load-bearing for lepton masses/M5.9). A uniform ψ then oscillates (bounded, real `+`mass²). |
| **Natural conserved measure is `1/r²`-weighted** | `L = 2∂_rr + (1/r²)Δ_Ω` is self-adjoint w.r.t. `dμ = d³x/r²`, not flat `d³x`. Numeric: flat-measure energy drifts 190% on a twist packet, the `1/r²`-weighted energy drifts 6%. Carry this measure into the M5.6.2 core treatment. |
| **Physical field = covariant gradient `(∇ψ−Â)`** | A phase that winds with the connection (`∇ψ ≈ Â`, e.g. `ψ=½ln(r²+r_c²)`) has `‖∇ψ−Â‖/‖∇ψ‖ ≈ 1.5e-3` — the gauge-invariant `(∇ψ−Â)` is the physical observable (the massless vacuum), not bare `∇ψ`. |

**Consequence for M5.6.2/.3:** the core treatment isn't just numerical hygiene — it *creates*
the mass. M5.6.2 handles the disclination + core on the `1/r²` measure; M5.6.3 (Faber) replaces
the ad-hoc `r_c` with the physical running-coupling scale that pins the lepton mass.

### 5b. M5.6.2a findings (`sandbox_v6/m5_6_2a_biaxial_hedgehog.py`, 2026-05-27) — `C_μν` is the matrix-level mass source

The scalar result (§5a) has a matrix-field counterpart on the **biaxial hedgehog** frame
`O = [r̂ | e_Θ | e_Φ]`, `D = diag(1, δ, 0)` (eigenvalue-1 axis radial; δ-twist + null axes
in the tangent plane). Verified numerically:

| Finding | Detail |
| --- | --- |
| **Frame structure** | `O` orthonormal in the bulk (`‖O^TO−I‖=4e-13`, det=+1); `M_bg=ODO^T` recovers eigenvalues `(1, δ, 0)`; principal director · r̂ = 1.0000. |
| **`C_μν = [M_μ^bg, M_ν^bg] ≠ 0`** | the background curvature is **nonzero** — the gauge source of the mass. **Contrast M5.5.2**: the single-generator bump `O=Ry(γ)` has `C_μν ≡ 0` (one generator commutes with its own derivatives) → no mass, which is exactly why M5.5.2 was massless. Biaxiality (multiple generators rotating across space) is what turns the mass on. |
| **`‖C‖ ∝ r^(−1.96)`** | the matrix curvature scales as **`1/r²`** — the same profile as §5a's scalar geometric mass `‖Â‖²=1/r²`. Two independent routes (scalar twist operator + matrix background curvature) give the same position-dependent mass. `C_μν` IS the matrix-level realization of the geometric KG mass. |
| **z-axis disclination** | biaxiality forces a hairy-ball line singularity on the z-axis (`e_Φ` winds). Regularized by a clamped smoothstep: the secondary `(δ, 0)` axes are full-length for `ρ ≥ ρ_c` (exact hedgehog) and **melt smoothly to 0 inside `ρ < ρ_c`** (biaxiality melts in the disclination core, like a nematic). Frame stays orthonormal outside the core; `‖∂O‖²` peak is **capped `∝ 1/ρ_c`** (sweep: `24→9.3→4.9` for `ρ_c = 0.4/0.8/1.2`). |

**M5.6.2b — dynamical confirmation** (`sandbox_v6/m5_6_2b_biaxial_evolution.py`): running the
validated M5.5.2 leapfrog (`2K ψ_tt = Σ_μ ∂_μ J_μ`, `J_μ=−32Σ_ν F̃_μν•P_ν`) on the biaxial
hedgehog, disclination-masked:

| Result | Detail |
| --- | --- |
| **`C_μν` SOURCES the twist** | The `C_μν` piece of `J_μ` is a ψ-independent source `S_μ=−32Σ_ν C_μν•P_ν`. At ψ=0: `max\|force\|` = **0.74 with `C_μν`, exactly 0.000 with `C=0`**; from ψ=0 the twist grows to 0.127 with `C` and stays **0.000** without. The biaxial hedgehog cannot sit static at ψ=0 — it drives its own twist. |
| **restoring / mass** | a seeded twist grows then oscillates around a balance (bounded, `0.12→0.22→0.17`); the massless M5.5.2 bump (`C=0`) had no restoring force. The mass scale tracks `‖C‖ ~ 1/r²`. |
| **conservation** | the FULL Hamiltonian `H=Σ(Kψ_t²+U)` drifts **0.76%** over 1500 steps (conservative). The ψ-sector energy grows +131% — NOT a drift: it's the `C`-drive pumping energy from the background curvature into the twist. (0.76% is non-tiny because the system is driven + disclination-masked; some leakage at the active-region Dirichlet boundary.) |

**Interpretation (flagged, not claimed):** `C_μν≠0` ⇒ ψ=0 is not a static solution ⇒ the defect
intrinsically oscillates. This is the *no-static-soliton / time-periodic* principle made dynamical
([[feedback_no_static_solitons]]) and is the plausible **M5.8 Zitterbewegung-clock seed** — the
clock *frequency* (`ω=2mc²/ℏ`) is M5.8's measurement, after M5.6.3 (Faber) pins the mass scale.

### 5c. M5.6.3a findings (`sandbox_v6/m5_6_3a_faber_regularization.py`, 2026-05-27) — Faber's regularization pins the mass scale

Faber's *Model of Topological Fermions* (Faber & Golubich, *Universe* 11/2025/113) gives the
**physical** regularization that sets the mass scale M5.6.1/.2 left as an ad-hoc `r_c`. Ported
natively (per [[reference_faber_regularization]], "port don't reinvent"):

| Faber (his Eq.) | Mapping / result |
| --- | --- |
| Field `Q = q₀ − iσ·q⃗`, unit 4-vector on S³ (`q₀=cos α`, `q⃗=n̂ sin α`) | the SU(2)/SO(3) twist field; `α` = rotation angle, `n̂` = orientation |
| **Regularization potential `Λ = q₀⁶/r₀⁴`** (Eq.4) | depends on `q₀=cos α`, i.e. the order-parameter **amplitude** `‖q⃗‖=sin α` — **not** the rotation direction. Consistent with the M5.5.3 finding that `V(M)` is rotation-invariant (acts on the eigenvalue/amplitude sector only). |
| Soliton `n̂=x̂`, `α=arctan(r/r₀)` ⇒ `q⃗ = x⃗/√(r₀²+r²)` | `‖q⃗‖→0` at the core (radial vector shrinks) — **identical to the M5.6.2a disclination/core melt**. Faber independently validates that core handling. |
| **Rest energy `E₀ = α_f ℏc · π/(4r₀)`** (Eq.8) | reproduced numerically (units `α_f ℏc=1`): `E₀ → π/4` as N↑ (**4.6% at N=181**, finite-difference + box-truncation limited); **`E₀·r₀` exactly constant (CV 0.0%) ⇒ `E₀ ∝ 1/r₀`**. |
| Physical anchor: `r₀ = 2.2132 fm` (= `π/4 ×` classical e⁻ radius) | `⇒ E₀ = 0.511 MeV` (electron). The regularization radius `r₀` IS the mass knob. |

**The deliverable:** the mass scale is no longer an ad-hoc `r_c` — it's `r₀`, fixed by the `Λ=q₀⁶/r₀⁴`
potential and tied to `α_f` + the mass via `E₀ ∝ 1/r₀`. This is the M5.9 lepton-mass-calibration handle.

**M5.6.3b — Faber's `Λ` mapped onto Duda's M** (`sandbox_v6/m5_6_3b_faber_on_M.py`): the
amplitude mapping is realized by spatially-melting eigenvalues. With `s(r)=sin α=r/√(r₀²+r²)`,
`D(s) = D_iso + s·(D_full − D_iso)` (`D_full=diag(1,δ,0)`, `D_iso=(1+δ)/3·I`), `M = O·D(s)·O^T`:

| Result | Detail |
| --- | --- |
| amplitude = eigenvalue spread | `mean\|spread − s(r)\| = 8e-16` — the melt makes M's eigenvalue spread *exactly* `s(r)`, →0 at the core. |
| both singularities melt | `M→scalar·I` at the core ⇒ point core **and** z-axis disclination regularized together; curvature energy `∫‖[M_μ,M_ν]‖²` drops **67%** vs the rigid `D_full` and is finite. |
| mass pinned `E ∝ 1/r₀` | `E·r₀` constant across `r₀` — Faber's `E₀∝1/r₀` reproduced on Duda's matrix substrate. |

**Caveat (honest):** `E·r₀` is *exactly* constant because the construction imposes Faber's
scale-covariant profile (`s=r/√(r₀²+r²)`, box `∝r₀`) rather than dynamically minimizing — so
`E∝1/r₀` is analytic for that profile, exactly as in Faber's own setup (he imposes `α=arctan(r/r₀)`).
This confirms the **M-framework reproduces Faber's scaling**; independently re-deriving the profile
via energy minimization is an M5.9-calibration BVP, not done here.

### 5d. M5.6.4 findings — the EM/Maxwell sector (the "1"-axis tilts)

The QM/twist (δ) sector gave Klein–Gordon (§5/§5a). The **EM sector is the "1"-axis tilts** (Duda's
eigenvalue map, `4a §8`): voxel-to-voxel tilts around the unity-eigenvalue axis. M5.6.4 verifies
these obey Maxwell, by two routes (`4a §11b`). **Structural note (the abelian/non-abelian split):**
Faber's curvature `R_μν=Γ_μ×Γ_ν` is *intrinsically non-abelian* (quadratic in the connection),
with `*F_μν = −(e₀/4πε₀c₀)R_μν` (Faber Eq.9) — the curvature is the **dual** field strength. Ordinary
Maxwell is the **long-range abelian limit**; the short-range non-abelian corrections are exactly
Faber's **running coupling** (`α_sol(d)`, reproduced in §5c/3a). The soliton is a *dual monopole* —
the topological winding sources Gauss's law.

**M5.6.4a — hydrodynamics ↔ EM dictionary** (`sandbox_v6/m5_6_4a_hydro_em.py`, the clean abelian
route, `4a §11b.1`): an incompressible tilt-flow `u=∇×A` with `ω=∇×u` (↔`B`), Lamb `l=ω×u` (↔`E`),
verified spectrally (periodic box) to reproduce Maxwell's structure to machine precision:

| Correspondence | Result |
| --- | --- |
| `∇·ω = 0` (↔ `∇·B=0`, no monopoles) | `3e-13` (kinematic) |
| `∇·l = u·(∇×ω) − ‖ω‖²` (↔ `∇·E=ρ`) | rel err `7.7e-14` — the turbulent-charge (Gauss) identity |
| `∂_tω = −∇×l` (↔ Faraday) via Euler step | rel err `2e-15` — vorticity transport = Faraday (curl of Lamb-form Euler; `∇×∇φ=0`) |
| Coriolis `−2(v×Ω)` ↔ Lorentz `q(v×B)` | `F∝v×ω`, magnitudes equal — the `v×field` force law (`Ω=ω/2 ⇒ B↔ω`) |

⇒ the hydro reading of "tilt = flow, vorticity = B" reproduces sourceless Maxwell + the Gauss
charge + the Lorentz force, all abelian and exact.

**M5.6.4b — Faber matrix-curvature route** (`sandbox_v6/m5_6_4b_faber_curvature_em.py`, the primary
route): on the regularized Faber hedgehog (3a), `Γ⃗_i = q0∂_iq⃗ − (∂_iq0)q⃗ + q⃗×∂_iq⃗` (Eq.6),
`R⃗_ij = Γ⃗_i×Γ⃗_j` (Eq.5; `*F_μν∝R_μν`, Eq.9):

| Check | Result |
| --- | --- |
| Maurer–Cartan / Bianchi | `dΓ_ij = ∂_iΓ_j−∂_jΓ_i = (1.995)·R_ij` (resid 1.1%) — `R` closed (the ≈2 is the su(2) factor) ⇒ `dR=0` ⇒ homogeneous Maxwell holds; the tilt curvature is a genuine field strength |
| abelian Coulomb far field | `‖R‖ ∝ r^(−1.99)` ⇒ `\|E\|~1/r²` at long range |
| running-coupling onset | `‖R‖·r²` = `0.45→0.73→0.94→0.99→1.00` (`r=0.5→7`): plateaus (abelian) far, rolls off within `r₀` (regularized non-abelian core) = effective coupling runs at scale `r₀`, matching Faber's `α_sol(d)` (§5c/3a) |

**M5.6.4 conclusion:** the EM/tilt sector reproduces Maxwell by both routes — abelian-exact via the
hydro dictionary (4a), and as Faber's closed non-abelian curvature that is abelian-Coulomb at long
range with a `r₀`-scale running coupling (4b). Together with the QM/twist KG sector (§5a–5c), the
eigenvalue map's two main axes (δ=QM twist, 1=EM tilt) are both verified emergent from the matrix field.

### §5e — M5.6.5a production port: biaxial seeder + the `ti.sym_eig` → analytic-eigensolver fix

The sandbox biaxial hedgehog (§5b) is now in the production engine. Two pieces landed, the second a
latent-bug fix that hardens the whole render+tracker pipeline.

**(1) `seed_biaxial_hedgehog_M`** (`engine1_seeds.py`) — builds, per voxel,

```text
M(x) = O(x)·D(s(r))·O(x)ᵀ,   O = [r̂ | e_Θ | e_Φ],   D = diag(1, δ, 0)
s(r) = r/√(r²+r₀²)          (radial eigenvalue melt → isotropic core, the §5c Faber profile)
e_Φ  *= smoothstep(ρ/ρc)    (clamped: biaxiality melts inside the z-axis disclination ρ<ρc)
```

Wired into `_launcher.py` as `TOPOLOGY_SEED["MODE"]="biaxial_hedgehog"` (config `xparameters/_biaxial1.py`,
knobs `R0_FRACTION`, `RHOC_VOXELS`, `BIAXIAL_DELTA`). **No auto-relax** for this mode — the M5.1
`relax_director_step` rebuilds `M` *uniaxially* from the principal director and would destroy the
biaxial structure; the biaxial `M` is constructed directly and is its own seed.

Headless verification (`/tmp/m5_6_5a_check.py`, N=47³, δ=0.3): `M` symmetric+finite; far-field
eigenvalues `(0.995, 0.301, 0.004) ≈ (1, δ, 0)`; principal director `· r̂ = 1.0000`; core melts to
isotropic (spread 0.598 near core vs 0.991 far); `C_μν=[M_μ,M_ν] ≠ 0` (`Σ‖C‖²>0` ⇒ the §5b mass
source is present in the production field).

**(2) ⚠️ Critical fix — `ti.sym_eig` is wrong for biaxial `M` on Metal/f32.** The first headless run
gave director recovery only **0.976**, not 1.0. Diagnosis (`/tmp/symeig_diag.py`): Taichi's `ti.sym_eig`
is accurate for **uniaxial/degenerate** `M` (`(1,δ,δ)`: err ~6e-8) but **catastrophically wrong for
biaxial/non-degenerate** `M` (`(3,2,1)`: eigenvalue err ~0.48). This is why the M5.4 feasibility spike
"passed" — it only tested the degenerate case. `f64` is not an escape: the `f64` `sym_eig` kernel
SPIRV-fails to compile on Metal.

| | uniaxial `(1,δ,δ)` | biaxial `(1,δ,0)` |
| --- | --- | --- |
| `ti.sym_eig` eigenvalue err | ~6e-8 ✅ | ~0.48 ❌ |
| director recovery | 1.0000 ✅ | 0.976 ❌ |
| analytic Cardano (the fix) | 1.0000 ✅ | 1.0000 ✅ |

Fix: replaced `principal_director` in `engine2_pde.py` with an **analytic symmetric-3×3 eigensolver**
(Cardano closed form — trace-shift `q`, deviatoric scale `p`, `φ=⅓·acos(det(B)/2)`, three eigenvalues
`q+2p·cos(φ + 2πk/3)`; principal eigenvector from the largest cross-product of `(M−λ₁I)` rows).
Validated against numpy `eigh` over 20 000 random symmetric matrices (f32): max eigenvalue err **6e-6**,
max director err **2e-7**. No regression — the uniaxial path is still 1.0000; the biaxial path goes
0.976 → **1.0000**. Since `eigen_decompose` is the lynchpin every render/tracker reads from, this is a
prerequisite for the M5.6.5b biaxial-ellipsoid glyph and the M5.8 clock (both genuinely biaxial states).

### §5f — M5.6.5c: turning V on — amplitude confinement, and why Eq.13 can't pin biaxiality

Rodrigo's M5.5.4 on-screen observation: with V off, Evolve-PDE makes the hedgehog *slosh*
and its energy **dilutes over a growing radius** — bounded and energy-conserving, but not
localized (no restoring force against amplitude spread). M5.6.5c turns the production `V_M`
(Eq.13 LdG, off by default) on to supply that force.

**The structural finding.** `V(M) = a·Tr(M²) − b·Tr(M³) + c·(Tr(M²))²` is rotation-invariant
(acts on eigenvalues only, §5/M5.5.3). Its eigenvalue-gradient is

```text
∂V/∂λ_i = λ_i·(2a − 3b·λ_i + 4c·s₂)  ,  s₂ = Tr(M²) = Σλ²
```

At a critical point each `λ_i` is either **0** or the single root `λ* = (2a+4c·s₂)/(3b)` —
one linear equation, shared `s₂` — so *all nonzero eigenvalues equal `λ*`*. The anisotropic
critical points are therefore **uniaxial** `(λ*,λ*,0)`, `(λ*,0,0)`. **The canonical 3-term
Eq.13 LdG cannot have a biaxial `(1,δ,0)` minimum** (three distinct eigenvalues). Verified
numerically for four `(a,b,c)` sets (`m5_6_5c_potential_confinement.py` Stage 1: max
nonzero-eigenvalue spread at the minimum = 0 over 120 random starts each). Consequence: a
`b≠0` term confines the amplitude but **pulls δ toward a uniaxial value — eroding the very
biaxiality** the C_μν mass source needs (§5b).

**The clean confinement (the production choice).** Set `b=0`: `V = a·Tr(M²) + c·(Tr(M²))²`
depends only on `s₂`, with minimum at `s₂* = −a/(2c)`. Choose `s₂* = Tr(diag(1,δ,0)²) = 1+δ²`.
This pins the amplitude (confines) and is **exactly flat in the biaxiality direction**
(`V` constant on the `s₂` sphere → `|V(biaxial) − V(uniaxial)| = 0` at equal `s₂`, Stage 2).

| Metric (full-M leapfrog, biaxial hedgehog) | V OFF | V ON (b=0 well) |
| --- | --- | --- |
| amplitude dev `⟨\|Tr(M²)−s₂*\|⟩` start→max→end | 0.025 → **0.158** → 0.158 (wanders 6.4×) | 0.025 → 0.025 → **0.022** (pinned) |
| energy RMS radius start→end | 3.33 → 3.92 (+18%) | 3.31 → 3.65 (+10%) |

**Production calibration (the units bridge).** The sandbox coefficients are dimensionless;
production `evolve_M` uses physical `dx_am`, `c_amrs`, and `dt_rs² ≈ 3.34·dx²`. The matrix
LdG force balances the **F²-curvature** force `c²·div(G) ~ c²·M³/dx⁴` (cubic, 4 gradient
orders) — NOT a Laplacian — so the coefficient scale is **`c²/dx⁴`**, not the scalar φ⁴'s
`(c/dx)²`. A sweep with the real kernel (`m5_6_5c_prod_scale.py`) confirms `K ∈ [0.5, 25]·c²/dx⁴`
all confine ~3.3× vs V-off with no blow-up (dt²-stable). The launcher computes

```text
ldg_c = K · c_amrs²/dx_am⁴ ,  ldg_a = −2·ldg_c·(1+δ²) ,  ldg_b = 0
```

from the xparameter `LDG_STIFFNESS_K` (off = 0). Configs: `_topo_biaxial1.py` (K=0, seed
smoke test) and `_topo_biaxial_v1.py` (K=1, the V-on confinement A/B demo).

**Energy-display fix (vacuum shift).** The b=0 well bottom is **negative**: `V(vacuum) =
−c·s₂*² ≈ −1.8e-6`. The production curvature density is tiny (`~1e-11` at `dx≈15`), so with
V on the constant well-bottom **swamps** the Hamiltonian — `compute_energyH_density_M` returns
a uniform `≈ −1.8e-6` and the flux mesh renders a featureless floor (looked like "energyH = 0"
on screen). Fix: subtract the vacuum potential `v0 = V_M(D_vacuum)` in the **display only**
(`compute_energyH_density_M` gained a `v0` arg; the launcher computes it from the vacuum
eigenvalues). A constant shift does not touch the force `−dV_M`, so dynamics/conservation are
identical. Shifted, the field is `≥ 0` with structure `[1.6e-11, 1.8e-6]` — vacuum at 0,
brightest at the core/disclination where M collapses (most deviated from vacuum). This is also
what makes the **confinement visible**: under Evolve-PDE the V-pinned core energy stays gathered
(it is the amplitude V pins) while the directors slosh (frame rotation — see below).

**What V does NOT do (the directors still slosh).** V is rotation-invariant (§5/M5.5.3): it
pins the eigenvalue **amplitude**, not the frame **orientation**. So under Evolve-PDE the
director glyphs keep sloshing even with V on — that motion is the dynamical twist sector (the
would-be clock), driven by the F²-curvature force, which V cannot and should not freeze.
Orientation containment is the job of the **gauge-correct `O(x)∈SO(3)` kinetic (5d)**, not V.
The disclination line also carries some energy outward regardless of V. So "fully contained"
is not achievable from V alone — V confines the amplitude component only.

**Q7 flag for Duda.** A fully biaxial-STABLE vacuum needs an *extra invariant* in V (the
3-term Eq.13 has only uniaxial minima). The `b=0` amplitude well is the interim — it confines
without uniaxializing, but leaves δ as a flat (un-pinned) direction. The biaxiality-stabilizing
term is an open question for Duda.

### §5g — M5.6.5d: the faithful (gauge-correct) kinetic vs the simple `½‖Ṁ‖²` we ship

Production `evolve_M` uses the simple kinetic `½‖Ṁ‖²` ⇒ `M̈ = c²·div(G) − dV` (every one of M's
6 symmetric components gets inertia ½). Duda's Eq.18 time-curvature gives the faithful kinetic
(§9): `T = 4Σ_μ‖[M_μ,Ṁ]‖² = ⟨Ṁ, G Ṁ⟩`, with the per-voxel metric `G = 4Σ_μ(−ad²_{M_μ})`.
`m5_6_5d_faithful_kinetic.py` characterizes the difference on the biaxial hedgehog:

| Finding | Result |
| --- | --- |
| **G is degenerate** | exactly **1 null mode per voxel** (median eigenvalues `[0, 0.08, 0.13, 0.27, 0.31, 0.90]`); the simple `½·I₆` has none |
| **The null mode is the TRACE** | null eigenvector · `(I/√3)` = **1.000** — `[M_μ, I]=0`, so the isotropic/dilation direction is non-dynamical under the faithful kinetic |
| **Simple kinetic stays physical** | the curvature force `div(G)` is **traceless**, so it never sources the trace mode: the simple scheme's motion is **0%** in the null space. ⇒ `½‖Ṁ‖²` does **NOT** generate spurious gauge slosh — it is a well-behaved approximation |
| **The real gap is physical-mode inertia (dispersion)** | the 5 physical eigenvalues span (5–95%) `[0.05, 1.45]` vs the simple uniform `0.5` ⇒ the twist/clock frequency is mis-set by **×[0.6, 3.0]** under `½‖Ṁ‖²` (`ω ∝ 1/√inertia`) |

**Consequences (corrects the earlier framing).** The on-screen director slosh under Evolve-PDE is
**physical twist** (the dynamical clock), not a gauge artifact — `½‖Ṁ‖²` does not animate spurious
modes. What the simple kinetic gets wrong is the **frequency** of that physical twist (an O(1)
inertia-weighting error), which matters only for the **M5.8 quantitative clock** `ω = 2mc²/ℏ`. The
faithful kinetic is the `O(x)∈SO(3)` metric already validated on the twist ψ DoF in `m5_6_1b`/`m5_6_2b`
(`2Kψ_tt = Σ∂_μ J_μ`, `K = 4Σ‖[M_μ,M_ψ]‖²`).

**Production recommendation:** do NOT rewrite `evolve_M` to the faithful kinetic — its degenerate
metric makes a full-M leapfrog implicit (per-voxel project onto `range(G)` + invert), a large change
that would not alter the qualitative GUI behaviour. Keep `½‖Ṁ‖²` for qualitative production runs;
measure the M5.8 clock frequency with the faithful **ψ-evolution** (`m5_6_2b` path). This closes
M5.6.5d as a *diagnosis* and routes the faithful kinetic to where it is actually needed (M5.8).

## 6. Matrix Hamiltonian (Eq.23) — the M5.4-carry-over `compute_energyH_density`

```text
ℋ = Σ_{0≤μ<ν} ‖F_μν‖²_F  +  V(M)  +  Σ_μ F_{0μ}•∂_μA_0
```

The last term vanishes in vacuum (integration by parts, as in EM). This replaces the M5.4
placeholder `compute_energyH_density` (which reads the dormant ψ buffer → uniform `¼λ`) and
gets the deferred physical-energy-scaling factor. Lands with the matrix leapfrog in the
launcher's `compute_propagation` no-op.

---

## 7. Duda's Mathematica source (the reference derivation)

The authoritative symbolic derivation lives in Duda's GitHub
`github.com/JarekDuda/liquid-crystals-particle-models`, file **`liquid crystal particles -
3D equations and hedgehog.nb`** (Wolfram, 705 KB / 12.7k lines — open in Mathematica for the
full version). Its operative code is **printed in the paper as Fig.9** (page 10); transcription:

### 7a. 3D evolution derivation (Fig.9 left)

```mathematica
d  = DiagonalMatrix[{1, δ, 0}];              (* ellipsoid shape Λ=(1,δ,0); δ→ℏ *)
Gx = {{0,0,0},{0,0,-1},{0,1,0}};             (* twist generator  *)
Gy = {{0,0,1},{0,0,0},{-1,0,0}};             (* tilt1 generator  *)
Gz = {{0,-1,0},{1,0,0},{0,0,0}};             (* tilt2 generator  *)
Ga = {{1,0,0},{0,0,0},{0,0,0}};              (* 3 elongation generators *)
Gb = {{0,0,0},{0,1,0},{0,0,0}};
Gc = {{0,0,0},{0,0,0},{0,0,1}};
Gpt = Join[Table[G.d + d.Transpose[G], {G, {Gx,Gy,Gz}}], {Ga,Gb,Gc}];  (* G' = GD+DG^T *)
com[A_, B_] := A.B - B.A;                    (* commutator *)
cd   = {{3,2},{1,3},{2,1}};
vect[m_] := Table[m[[ cd[[i,1]], cd[[i,2]] ]], {i,3}];   (* antisym matrix → rotation vector *)
(* Γ_μ affine connection (matrix form), M_μ = com[Γ_μ, …], F_μν = Simplify[com[M_μ, M_ν]] *)
(* vrip = Tr[ F_μν·( com[Γ_μ, M_ν.Gp] − com[Γ_ν, M_μ.Gp] ) + F_{μν,μ}·M_ν.Gp − F_{μν,ν}·M_μ.Gp ]  *)
(*        — integrate-by-parts form of Eq.35, looped over μ,ν, generators Gp *)
vr   = Simplify[Series[vrip, {δ, 0, 0}]];    (* low-order in δ *)
fin  = Table[Sum[…,{μ},{ν}] − Sum[…], {i,3}, {v, vr}];   (* + Lagrangian → evolution terms *)
(* rename R-curvatures → B/E fields → Column shows: ~Klein-Gordon, ~Maxwell1, ~Maxwell2 *)
```

### 7b. Hedgehog application (Fig.9 right)

> **Transcription caveat (M5.5.1, 2026-05-26):** the exact Euler-angle assignment below
> was read off the small Fig.9 image and is **ambiguous** — using `θ→ArcTan[√(x²+y²),z]`
> on `Gz` and `φ→ArcTan[x,y]` on `Gy` as written gave a NON-radial director in the sympy
> mirror (`m5_5_1_evolution_symbolic.py`, `mean|n̂·r̂|=0.40`). The hedgehog is therefore
> built **physics-first** in the sandbox — `O = [r̂ | e_Θ | e_Φ]·Rx(ψ)` with the radial
> director `r̂` as the first column by construction (verified `|n̂·r̂|=1.000`, `det O=+1`,
> connection `Γ_i=O^T∂_iO` antisymmetric + `~1/r`). The precise Euler parameterization is a
> convention detail; any `O` with radial director + ψ-twist about it is the hedgehog. Pull
> the exact angles from the actual `.nb` only if a line-by-line Mathematica match is needed.

```mathematica
sph = {x -> r Cos[θ] Cos[φ], y -> r Cos[θ] Sin[φ], z -> r Sin[θ]};
Q0  = FullSimplify[ MatrixExp[θ Gz].MatrixExp[φ Gy].MatrixExp[ψ Gx] ]
        /. {φ -> ArcTan[x, y], θ -> ArcTan[Sqrt[x^2+y^2], z]};   (* hedgehog *)
Q   = Q0 /. ψ -> ψ[t, x, y, z];                                  (* assume phase dependence *)
BE  = Simplify[Table[ vect[Transpose[Q].D[Q, v]], {v, {t,x,y,z}} ], r > 0];  (* B,E fields *)
(* substitute BE + derivatives into fin[[1;;3]]; first eq → Klein-Gordon-like, 2nd/3rd → 0 *)
A   = {x, y, z}/r^2;                                              (* Â^hedg *)
gmA[f_] := Grad[f, {x,y,z}] - A f;
Adg[f_] := (A + r) . Grad[f, {x,y,z}];
(* verify:  fne[[1]]/r^2  ==  Σ_i gmA[gmA[ψ]][i] + Adg[Adg[ψ]] + 2 ∂_tt ψ   *)
```

---

## 8. sandbox_v5 plan (mirror the Mathematica in sympy/numpy, then Taichi)

| Sub-step | Deliverable | Validates against |
| --- | --- | --- |
| 5.5.0 ✅ | this doc — action + Eq.35 + Fig.9 source transcribed; math confirmed | — |
| 5.5.1 ✅ | sympy (`m5_5_1_evolution_symbolic.py`): operator identity `F_μν=2[M_μ,M_ν]` (Eq.20) ✓; radial hedgehog (director=r̂, det O=+1) ✓; connection `Γ_i=O^T∂_iO` antisymmetric + `~1/r` (the singular `Â^hedg`) ✓ | the symbolic FOUNDATION is verified. Full action→KG EL reduction is NOT done symbolically (a `(∂M)⁴` action) — **moved to numerical validation** (decision 2026-05-26): dispersion `ω(k)` in M5.5.2/M5.6. Also caught + corrected a Fig.9 angle-transcription error (see §7b caveat). |
| 5.5.2 ✅ | numpy (`m5_5_2_twist_evolution.py`): **first numerical evolution of the Eq.18 action** (twist sector, V=0) on a smooth non-singular background. `K(x)=4Σ‖[M_μ^bg,M_ψ]‖²>0` (twist dynamical where the bg varies) ✓; **energy drift 0.21%** over 1500 leapfrog steps ✓; stable. Derived from `ℒ=T−U` via `F_μ0=2[M_μ,Ṁ]`. | **MACHINERY + energy conservation validated.** Static Coulomb already validated (M5.4 page-18 = `¼Σ‖F_μν‖²`). **KG mass gap MOVED to M5.6** — it needs the biaxial HEDGEHOG, whose `δ≠0` perpendicular frame carries a z-axis **DISCLINATION singularity** (`e_Φ~1/ρ`, hairy-ball — unavoidable for biaxial) + the core; disclination handling is M5.6-level. Minor: amplitude grows at the `K→0` active-region edge (metric-degeneracy artifact) → damp in M5.6. |
| 5.5.3 ✅ | `V(M)` Eq.12/13 + regularization (`m5_5_3_potential_regularization.py`). **Finding: V(M) is rotation-invariant** (`V(ODOᵀ)=V(D)`, eigenvalues fixed) → it does NOTHING to the twist/rotation sector (why M5.5.2 needed V=0); V acts ONLY on the eigenvalue-deformation sector = **the regularization**. V defines vacuum Λ (perturbed λ relax to Λ) ✓. **Eq.18 `F²`+V are Derrick-stable**: rigid core `∫‖F‖²` diverges, deformed core finite; `E(L)=E_F/L+E_V·L³` → finite minimum (core size + mass), no extra Skyrme term (F² IS the Skyrme-family kinetic). | exact `Λ=(1,δ,0)` LdG coeffs (Eq.13 a,b,c) + Faber's exact scheme remain **Q7/Q8** (Duda open); baseline shown. |
| 5.5.4 ✅ | **Taichi port** done. `engine2_pde`: `compute_curvature_flux` (G_α=8Σ[[M_α,M_ν],M_ν]) + `evolve_M` (leapfrog) + `V_M`/`dV_M`; `medium`: `curv_flux_*` + `swap_matrix_buffers`; `engine3`: `compute_energyH_density_M`; launcher: `compute_propagation` → matrix leapfrog + `eigen_decompose` (Evolve PDE LIVE), `compute_energyH_density_M` (WAVE_MENU=4 = matrix ℋ). **Energy-conserving** (`m5_5_4_matrix_evolution_check.py`: secular drift 2.15%→1.13%→0.03% as dt→0 = symplectic; the flat 6% osc is the velocity-measurement convention). Integrated path verified (M+director evolve). **Resolves both M5.4 carry-overs.** | e_scale=1 (bare); physical-energy calibration tied to a reference mass → M5.9. Production-scale dynamics speed via SIM_SPEED/c_amrs. **On-screen VERIFIED (Rodrigo 2026-05-26) — "bounded-not-bug"**: a hedgehog under "Evolve PDE" *sloshes* (not the old explode-and-propagate wave); a headless mirror of the GUI scenario (63³, dx≈15.2, V off, 1200 steps) holds **H conserved to 5 digits** with `max‖M−D‖_F` bounded `0.7→2.2→1.6`, finite — correct nonlinear curvature dynamics, NOT a blow-up. The wave-propagation look belonged to the retired *linear* ψ leapfrog. |
| 5.5.5 → **M5.6** | EM-from-tilts cross-check + Faber EM Lagrangian (4a §11b) — **folded into M5.6** (overlaps its Maxwell-sector verification, Eq.37–38) | superfluid-vorticity ↔ Maxwell dictionary |

**Exit criterion (M5.5):** the full Eq.18 action runs; defect dynamics are governed by the
real proposed action (not the M5.1 Frank-energy approximation); KG emerges from the twist mode.

---

## 9. M5.5.2 evolution — EOM derivation + the kinetic-degeneracy finding

**Kinetic term (derived 2026-05-26, confirmed by Rodrigo).** The "electric"/time curvature
parallels Eq.20:

```text
A_0 = [M, Ṁ]                                       (Ṁ = ∂_0 M)
F_μ0 = ∂_μA_0 − ∂_0A_μ = 2[M_μ, Ṁ]                 (the [M,·] terms cancel: ∂_μṀ = ∂_0∂_μM)
```

So the action splits into a proper kinetic + potential:

```text
ℒ = T − U
T = Σ_{μ=1..3} ‖F_μ0‖²_F = 4 Σ_μ ‖[M_μ, Ṁ]‖²_F            kinetic (quadratic in Ṁ)
U = Σ_{μ<ν} ‖F_μν‖²_F + V(M) = 4 Σ_{μ<ν} ‖[M_μ,M_ν]‖²_F + V(M)
EOM:  ∂_0(∂T/∂Ṁ) = ∂T/∂M − ∂U/∂M
```

`M`, `M_μ`, `Ṁ` are all real-symmetric → `[M_μ,Ṁ]` is antisymmetric; `T ≥ 0`.

**THE FINDING — the kinetic metric is DEGENERATE (gauge structure):** `T = 4Σ_μ‖[M_μ,Ṁ]‖²`
vanishes whenever `Ṁ` commutes with every spatial gradient `M_μ`. Consequences that shape
M5.5.2:

| Implication | Detail |
| --- | --- |
| Evolve the rotation DoF, not M's 6 components | The dynamical variable is `O(x) ∈ SO(3)` (3 angles/voxel; `D` frozen). `M`'s 6 free components include gauge/non-dynamical directions — a naive leapfrog on `M` hits the degenerate (non-invertible) metric. Parameterize by `O` (or angular velocity `ω = OᵀȮ`). |
| The twist is dynamical ONLY on a non-uniform director background | A pure single-axis twist `M=Rx(ψ)D Rx(ψ)ᵀ` has `M_x ∝ M'` and `Ṁ ∝ M'` → `[M_x,Ṁ]=0` → **T=0**. So a uniform-axis 1D twist is in the kinetic null space (NOT dynamical). The KG-for-twist emerges only when the twist couples to a position-dependent director — i.e. the **hedgehog**. |
| Minimal KG test is inherently 3D hedgehog + twist | This is exactly the Fig.9 case (twist phase `ψ(t,x,y,z)` on the hedgehog → KG with `Â^hedg`). So M5.5.2's KG-dispersion validation converges with M5.6's headline. |

**Refined M5.5.2 plan:** evolve `O(x) ∈ SO(3)` (gauge-clean DoF) on a small 3D grid, with the
hedgehog as the background + a small twist-phase perturbation; the field-dependent kinetic
metric `T = 4Σ‖[M_μ,Ṁ]‖²` is non-degenerate on the rotation DoF. Validate Eq.23 energy
conservation + KG dispersion `ω(k)` of the twist. (Static Coulomb already validated via the
M5.4 page-18 energy = `¼Σ‖F_μν‖²`.) This is a substantial numerical build — the M5.5↔M5.6 core.
