# M5.5 — Paper Lagrangian + Evolution Equation (math reference)

**Purpose:** the confirmed mathematical foundation for M5.5 — Duda's Eq.18 action, the building-block operators, the Eq.35 Euler–Lagrange evolution of the matrix field `M`, the matrix Hamiltonian, and the `V(M)` options. Plus the transcription of Duda's Mathematica source (paper Fig.9) that derives the evolution equation symbolically and reduces it to the hedgehog Klein–Gordon — the reference we mirror in the `sandbox_v4` sympy prototype.

**Source:** Duda, *Framework for liquid crystal based particle models* (arxiv:2108.07896 v7), §II–IV + Fig.9. Math reading **confirmed by Rodrigo 2026-05-26**.

**Sister docs:** [4a_convo_2026.05.12.md](4a_convo_2026.05.12.md) (paper digest, eigenvalue map), [0b_M5_roadmap.md § M5.5](0b_M5_roadmap.md), [1a_lagrangian_framework.md](1a_lagrangian_framework.md).

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

## 8. sandbox_v4 plan (mirror the Mathematica in sympy/numpy, then Taichi)

| Sub-step | Deliverable | Validates against |
| --- | --- | --- |
| 5.5.0 ✅ | this doc — action + Eq.35 + Fig.9 source transcribed; math confirmed | — |
| 5.5.1 ✅ | sympy (`m5_5_1_evolution_symbolic.py`): operator identity `F_μν=2[M_μ,M_ν]` (Eq.20) ✓; radial hedgehog (director=r̂, det O=+1) ✓; connection `Γ_i=O^T∂_iO` antisymmetric + `~1/r` (the singular `Â^hedg`) ✓ | the symbolic FOUNDATION is verified. Full action→KG EL reduction is NOT done symbolically (a `(∂M)⁴` action) — **moved to numerical validation** (decision 2026-05-26): dispersion `ω(k)` in M5.5.2/M5.6. Also caught + corrected a Fig.9 angle-transcription error (see §7b caveat). |
| 5.5.2 | numpy: evolve `M` on a small grid with Eq.12 `V(M)`; energy conservation; KG dispersion `ω(k)` | **NOTE:** the STATIC Coulomb under the real action is ALREADY validated — the M5.4 `m5_4_coulomb_page18.py` energy `Σ‖[∂_iM,∂_jM]‖²` IS `¼Σ‖F_μν‖²_F` (Eq.20), R²=0.9959. So 5.5.2's genuinely-new piece is the DYNAMICAL evolution + energy (Eq.23) conservation. |
| 5.5.3 | graduate `V` → Eq.13 LdG; port Faber regularization (running coupling) | Faber's running-coupling deformation |
| 5.5.4 | **Taichi port** — matrix leapfrog → `compute_propagation`; matrix ℋ → `compute_energyH_density` + physical scaling (M5.4 carry-overs) | parity with the sandbox |
| 5.5.5 | EM-from-tilts cross-check + Faber EM Lagrangian (4a §11b) | superfluid-vorticity ↔ Maxwell dictionary |

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
