# LAGRANGIAN FRAMEWORK

A research sub-project evaluating whether a Lagrangian formulation could replace OpenWave's empirical wave-equation search with a first-principles derivation.

---

## WHY WE'RE DOING THIS

### The Problem

Even if OpenWave's M2 / Laplacian-Propagation Method already used a Lagrangian — implicitly. M2's PDE solver evolves `∂²ψ/∂t² = c²∇²ψ`, which IS the Euler-Lagrange equation of the simplest possible field Lagrangian: `L = ½(∂ψ/∂t)² - ½c²(∇ψ)²` (the free wave Lagrangian, V=0). So we've already done Lagrangian field theory — just with the most basic version, and on a scalar field. M2 produced a self-stabilized standing-wave background, but the wave centers (modeled as Dirichlet ψ=0 boundary points) were invisible to the field — the wave passes through without reorganizing. Three structural limitations explain this: (1) the linear Lagrangian has no soliton solutions, (2) the scalar field has no topological structure (no winding numbers), (3) WCs as boundary conditions are passive constraints rather than active features built into the field. The new approach fixes all three: non-linear potentials V(ψ) create solitons, vector fields enable topology, and topological defects (hedgehogs) ARE part of the field configuration rather than imposed externally. See [`0_WAVE_EQUATION.md`](0_WAVE_EQUATION.md) for the full M2/Lagrangian analysis.

OpenWave's M3 scalar method tested 5 wave equation candidates empirically (Wolff, LaFreniere-Marcotte, Phase-warped, Combined W-L, Weighted PSW). The best candidate (Combined Wolff-LaFreniere) produces particle lock-in and K=10 tetrahedron stability at perfect placement, but fails under perturbation and cannot produce far-field Coulomb (sinc barriers flip the force direction every λ/2).

We've been searching for the perfect wave equation as the **goal** of our research. A Lagrangian formulation inverts this: the wave equation becomes the **consequence** of a deeper choice (the Lagrangian itself, chosen from symmetries and physical principles).

### What a Lagrangian Would Give Us

If a valid Lagrangian for OpenWave can be found, it would tell us:

- **The exact wave equation** — derived from Euler-Lagrange, not guessed empirically
- **The exact energy functional** — from the Hamiltonian `H = T + V`, possibly different from our current `E = ρV(fA)²`
- **Conservation laws for free** — via Noether's theorem (energy, momentum, angular momentum automatically conserved)
- **Why charge is quantized** — as a topological invariant of the Lagrangian's symmetry group
- **What conservation laws hold** — and which don't, under which conditions
- **Connection to established physics** — QFT, GR, Standard Model formalisms

This is exactly the gap Duda identified when he asked "show me your Lagrangian?"

### The Two Challenges from Duda

1. **"Show me your Lagrangian"** — we test wave equations empirically; we need to derive them from a variational principle
1. **"Without charge quantization your electron explodes"** — our `cos(source_offset)` is imposed as ±1; topology forces integer charges naturally via Gauss-Bonnet
1. **"Particles as topological defects, not standing waves?"** — defects (hedgehogs) protected by topology vs. our standing waves protected by interference

### The Unification Insight

Duda explicitly confirmed that **topology and standing waves are NOT mutually exclusive** — both are needed:

- **Topology** → charge and spin (integer-valued conserved quantities)
- **Standing waves** → orbit quantization (validated by Couder bouncing droplet experiments)

So OpenWave's M3 results (lock-in, K=10 tetrahedron, annihilation) stay valid as the dynamical content — topology adds the geometric skeleton that handles charge. This is additive, not replacement.

### The Approach

We will implement 8 numpy research scripts in `research_hub/sandbox_phase2_lagrangian/` (same pattern as Phase 1 vector/scalar scripts) to test the core ideas before committing to any architecture change. None of these tests require refactoring M3 or M4 — they're standalone exploration scripts that can validate (or rule out) the Lagrangian framework quickly.

See the [What We Can Test in OpenWave](#what-we-can-test-in-openwave) section below for the full test plan.

---

## EMAIL THREAD

On 4/7/2026 6:50 PM, jeffsyee wrote:

Hello all — Jeff here. I’ve been following discussions in this group and wanted to share a project I’m involved in that may be of interest.

An open source simulation framework called OpenWave (developed by Rodrigo Griesi) models particle emergence and interactions using classical wave physics in an aether-like medium.

It already supports wave propagation, standing wave formation, and early particle-like stability, along with interactions such as attraction/repulsion and annihilation behavior.

The goal is to provide a platform where different wave-based or aether-based models can be tested and explored through simulation.

If anyone here is working on similar ideas, it may be useful for running experiments or validating assumptions.

Overview video:
<https://www.youtube.com/watch?v=m51-OQ4mJ_Q>

Repository:
<https://github.com/openwave-labs/openwave>

Would be very interested in feedback or thoughts from this group.

--
You received this message because you are subscribed to the Google Groups "Models of particles" group.
To unsubscribe from this group and stop receiving emails from it, send an email to <models-of-particles+unsubscribe@googlegroups.com>.
To view this discussion visit <https://groups.google.com/d/msgid/models-of-particles/718f73e8-4127-48a6-8d61-1c732015ab54n%40googlegroups.com>.
For more options, visit <https://groups.google.com/d/optout>.

---
On Apr 7, 2026, at 11:55 AM, Jarek Duda <dudajar@gmail.com> wrote:

Hi Jeff,

Thanks, looks very nice!

Many people here have own models, simulations ... but without collaborations these are mostly single person models, dying with the author ...

Maybe let's try to finally find some common language - could you show here you main Lagrangian/PDE so we can discuss its choice?

How do you resolve some basic questions like charge quantization?

Seems most people here agree e.g. that nuclei are knots of topological vortices - do you agree/disagree? Why?

With best regards,
Jarek

---
On 4/7/2026 9:41 PM, Jeff Yee wrote:

Jarek,

Really appreciate the thoughtful questions.

On the core formulation:
OpenWave, designed by Rodrigo Griesi, already operates using PDE-based wave-field simulations (both scalar and vector field approaches, similar in spirit to lattice methods). The system evolves wave dynamics over a discretized grid and explores how structure emerges from these interactions.

However, it is not yet built around a single canonical Lagrangian. The current approach is more exploratory - testing different wave-field formulations and observing which configurations produce stable, physically meaningful behavior. Converging toward a more formal underlying Lagrangian is definitely a goal, but still in progress.

On charge quantization:
Charge is not imposed as a discrete input in the model. The working hypothesis is that it emerges from stable standing wave configurations - particularly phase structure and symmetry.

At a basic level, we are observing that standing wave systems tend to organize into phase-opposed configurations, which may correspond to what we interpret as positive and negative charge. In that sense, the ± symmetry could arise naturally from the allowed phase states of stable structures.

Quantization would then follow from stability constraints - only certain configurations persist. That said, this is still early, and a rigorous derivation of charge quantization has not yet been completed.

On nuclei as topological vortices:
That’s a very interesting direction. OpenWave is currently modeling particles as standing wave structures rather than explicitly as topological knots or vortices. However, there may be overlap - especially if stable configurations end up having topological invariants. This is something that could be explored further within the framework.

More broadly, I completely agree with your point about fragmentation. One of the motivations behind OpenWave being open source is to provide a shared simulation platform where different approaches (including vortex-based models) can be tested in a common environment.

Jeff

---

On Apr 7, 2026, at 7:51 PM, Jarek Duda <dudajar@gmail.com> wrote:

Jeff,

Thanks for the answer, but finally you need Lagrangian - we can discuss if have any candidates?

Regarding charge quantization, standing wave is perfect for orbit quantization, but here you need quantized charges as vacuum excitations, with Coulomb interaction - how would you like to get it with standing waves?

As for vortices, they rather need topological quantization - both obtained experimentally in liquid crystals together with Coulomb: <https://en.wikipedia.org/wiki/Draft:Liquid_crystal_particle_analogs>

Best,
Jarek

![alt text](images/triangulo.png)

---
On 4/8/2026 6:12 AM, Jeff Yee wrote:
Jarek,

Thanks, these are great points.

On the Lagrangian:
I agree that ultimately a well-defined Lagrangian is needed to properly frame any theory and connect it to established formalisms. OpenWave currently operates at the PDE / field evolution level, exploring candidate wave-field dynamics and their emergent behavior.

Since Rodrigo Griesi is leading the implementation and formulation work, I’ll defer to him on current or candidate Lagrangian directions -  he’ll be able to speak to that more precisely.

On charge quantization:
Yes - this is a key challenge. Standing waves naturally give quantized modes, but that alone does not explain quantized charge or Coulomb coupling.

The direction being explored is that stable configurations are constrained by phase relationships at nodes, leading to phase-opposed states that could correspond to positive and negative charge. Interaction effects would then arise from interference in the surrounding field.

That said, deriving discrete charge magnitude and a Coulomb-like interaction rigorously is still an open problem.

On vortices / topology:
Thanks for sharing that reference - very interesting. OpenWave currently models particles as standing wave structures, but we have not yet framed them explicitly in terms of topological invariants or vortex solutions.

That said, if stable configurations correspond to topologically protected structures, there could be a natural connection here. This is definitely an area worth exploring further as the model develops.

More broadly, one of the goals of OpenWave is to provide a shared computational platform where different approaches - wave-based, topological, or otherwise - can be explored and compared in a common environment.

Would be very interested to hear how others here are approaching charge and topology in their models.

Jeff

---
On Apr 7, 2026, at 9:51 PM, Jarek Duda <dudajar@gmail.com> wrote:

Jeff,

Without charge quantization you don't have particles, your charged e.g. electron would split into parts - literally exploding.

The only charge quantization mechanisms I have seen e.g. in this mailing list are topological - used in liquid crystal analogs (e.g. in <https://en.wikipedia.org/wiki/Draft:Liquid_crystal_particle_analogs> ), and easy to formalize Faber's way:

define curvature of deeper field as electric field (generally dual F tensor), this way Gauss law counts topological charge of this deeper field - becoming <https://en.wikipedia.org/wiki/Gauss%E2%80%93Bonnet_theorem>

If you don't have any, maybe consider this mechanism - a few persons here would gladly discuss, collaborate.

Best,
Jarek

![alt text](images/electromagnetism.png)

---
On 4/8/2026 6:51 PM, Jeff Yee wrote:
Thanks, Jarek.

Really appreciate you sharing this, it’s a very interesting direction and will definitely consider it. This is exactly why we built OpenWave as an open source framework - to make it easier for ideas like this to be tested and explored collaboratively, especially across wave-based models.

Jeff

---
On Apr 7, 2026, at 10:00 PM, Jarek Duda <dudajar@gmail.com> wrote:
Great, Jeff,

Let us know if something, would gladly discuss - especially the choice of Lagrangian.

You can start with recreating EM for quantized topological charges like in liquid crystals - e.g. integrating Hamiltonian for two charges in various distances leads to Coulomb potential (fig. 2 in <https://arxiv.org/pdf/2108.07896> ), you could also simulate dynamically.

While many of us agree here, we disagree what to do next - e.g. I believe we should distinguish twists of these vectors, as in basic liquid crystal model: <https://en.wikipedia.org/wiki/Landau%E2%80%93de_Gennes_theory>  - these low energy twists behave like quantum phase, e.g. getting Klein-Gordon-like equation, unifying EM+QM.

Then it is natural to extend to 4D liquid crystal as in Einstein's teleparallelism ( <https://en.wikipedia.org/wiki/Teleparallelism> ) - adding boosts getting 2nd set of Maxwell equations for gravity as in GEM ( <https://en.wikipedia.org/wiki/Gravitoelectromagnetism> ).

Beside point-like topological charges corresponding to elementary electric, there are also topological vortices - in particle physics corresponding to quark string/gluon flux tube, e.g. decaying into particles in LHC string hadronization ( <http://www.scholarpedia.org/article/Parton_shower_Monte_Carlo_event_generators#String_model> ), so we just need to find particle-defect correspondence, e.g. various size knots - various size nuclei ...

Let me know if you agree/disagree here?

Best,
Jarek

Robert, have you maybe finally got charge quantization - explained why your electron doesn't explode?

---

On Apr 8, 2026, Robert Close <robert.close@classicalmatter.org> wrote:

Jeff,

Thank you for sharing this tool. I think it could be quite useful. I have proposed a simple nonlinear vector wave equation (and equivalent first-order Dirac equation) for the evolution of spin density, consistent with the dynamical operators of relativistic quantum mechanics. Spin density, in classical and quantum physics, is the vector field whose curl is equal to twice the intrinsic (i.e. aether) momentum density. I don't know when I'll find time to work on solutions, so I invite you and others to try. This paper describes the theory:
Plane Wave Solutions to a Proposed "Equation of Everything" - Foundations of Physics
<https://doi.org/10.1007/s10701-025-00839-0>

I suggest starting with a spherical harmonic linear wave solution and see what evolves.

Best regards,

Robert Close

---

On 4/9/2026, Rodrigo Griesi replied with OpenWave status, open challenges, 5 proposed tests, and compatibility assessment (see INITIAL ANALYSIS below for full content).

---

On Apr 9, 2026, Jarek Duda <dudajar@gmail.com> wrote:

Hi Rodrigo,

Thanks, sure we need both topological type quantization (e.g. charge, spin), but also different standing wave type - especially for orbit quantization, see e.g. <http://www.pnas.org/content/107/41/17515.short> , <https://www.nature.com/articles/ncomms4219>

While you assume oscillations, finally we need to derive them from deeper theory, explaining why e.g resting electron/neutrinos oscillate - as in original definition of time crystals, we automatically get such mechanism - toy model: <https://arxiv.org/pdf/2501.04036>

5 concrete experiments we can set up:
Hedgehog energy vs distance [...]

Indeed getting Coulomb for quantized charges should be the first step of particle models ... but this calculation is simplified, e.g. missing regularization leading to running coupling effect: <https://www.mdpi.com/2076-3417/16/2/1030>

Klein-Gordon is effective, e.g. having probability distribution of particles - we should search for deeper model averaging e.g. to Klein-Gordon, e.g. Fig. 9 in <https://arxiv.org/pdf/2108.07896> has such derivation around electron.

Lagrangian should be assumed, comparing its consequences with known physics - I assume Landau-de Gennes field with Skyrme-like kinetic term for Lorentz covariance, but there are open questions e.g. how to choose its Higgs-like potential.

Indeed in 3D with 3 distinguishable axes we have 3 families of topological defects, e.g. elementary charge as hedgehog of one of 3 axes - same electric field, different mass, can decay to the lightest (by field rotation) - like 3 leptons.

With best regards,
Jarek

---

On Apr 9, 2026, Rodrigo Griesi replied acknowledging the new territory, committing to study the theory and work through simulations.

---

On Apr 9, 2026, Jarek Duda <dudajar@gmail.com> wrote:

Hi Rodrigo,

Sure, for intuitions about QM study walking droplet hydrodynamical quantum analogs ... which need to be combined with particles: stable localized field configurations, e.g. topological like in these liquid crystal particle analogs - e.g. preventing electron from exploding.

<https://en.wikipedia.org/wiki/Sine-Gordon_equation> is a good starting point for such stable massive particles, with pair creation/annihilation and entire special relativity ... having mechanical models as lattice of coupled pendula: <https://www.youtube.com/watch?v=nl5Qq5kUbEE>

To combine them, we need to get propulsion of its oscillations - for electron: <https://en.wikipedia.org/wiki/Zitterbewegung> , direct experimental confirmation: <https://link.springer.com/article/10.1007/s10701-008-9225-1> ... this Faber's approach in 4D automatically gives such propulsion: <https://arxiv.org/pdf/2501.04036>

Best,
Jarek

---

## DUDA'S FIRST REPLY — KEY INSIGHTS

### 1. Both Topological AND Standing Wave Quantization Are Needed

Duda explicitly agrees: topological quantization handles charge and spin, but standing wave quantization handles **orbit quantization** — different mechanisms for different physics. This validates our dual approach.

He cites two landmark experiments showing orbit quantization in classical wave-particle systems (bouncing droplets on vibrating fluids):

- **Fort, Couder et al. (2010)**: bouncing droplets on a vibrating fluid surface generate quantized orbital states when the wave-field "memory" spans the full trajectory. Macroscopic angular momentum quantization from path memory.
- **Perrard, Couder et al. (2014)**: walking droplets confined in a harmonic potential self-organize into discrete stable orbits with quantized spatial extent and angular momentum — quantum-like eigenstates from classical wave-particle dynamics.

**Connection to OpenWave**: these experiments are the closest physical analog to what we simulate — wave centers (droplets) interacting with their own wave field (standing waves on fluid surface) and exhibiting quantized behavior. Our sinc lock-in wells ARE orbit quantization. The topological charge (Duda) and orbit quantization (our standing waves) are complementary, not competing.

### 2. Time Crystals — Why Resting Particles Oscillate

Duda points out that we *assume* oscillations (our wave centers oscillate at frequency f₀) but we should *derive* them from a deeper theory. His time crystal mechanism provides this: a scalar field with curvature coupling naturally produces mass-driven periodic oscillations (kinks that oscillate without external driving).

- **Duda (2025)**: "Time crystal phi-4 kinks by curvature coupling as toy model for mechanism of oscillations propelled by mass" — a phi-4 field `V(ψ) = k·(ψ²-1)²` with curvature coupling creates kinks that oscillate with frequency proportional to mass. This is the *origin* of the oscillation, not an assumption.

**Connection to OpenWave**: this directly addresses our Phase 3 (Time Dynamics) question: "why does the medium oscillate?" In our framework, we take f₀ = 10²⁵ Hz as given. Duda's time crystal mechanism would *derive* this frequency from the field's potential. The phi-4 potential `(ψ²-1)²` is closely related to Smolinski's quartic potential `k·ψ⁴` — both are non-linear potentials that come from valid Lagrangians.

### 3. Coulomb Calculation Needs Regularization

Duda confirms hedgehog energy vs distance is "the first step" but warns the simple calculation is missing **regularization** leading to the **running coupling effect** — the effective charge strength changes with distance (energy scale). This is a known QFT phenomenon (the fine-structure constant α is not truly constant — it "runs" with energy).

**Connection to OpenWave**: our energy `E = ρV(fA)²` with variable λ(r) from Yee & Hauger shells already has a form of regularization (the wavelength changes near the core, modifying the effective coupling). The running coupling may be the topological equivalent of our variable-λ(r) profile.

### 4. Klein-Gordon Is Effective, Not Fundamental

Duda clarifies: Klein-Gordon describes the *probability distribution* of particles (effective/averaged behavior), not the deeper field dynamics. The deeper model should *average to* Klein-Gordon — see his Fig. 9 in arXiv:2108.07896 for the derivation around an electron.

**Implication for Test 4**: we should not expect the twist dynamics to exactly reproduce Klein-Gordon. Instead, we should check whether the time-averaged twist behavior statistically converges to Klein-Gordon — a weaker but more physically correct test.

### 5. Lagrangian: LdG + Skyrme Kinetic Term, Open Potential

Duda's specific Lagrangian choice: **Landau-de Gennes field with Skyrme-like kinetic term** (for Lorentz covariance). The open question is the Higgs-like potential V(M) — its specific form is not yet determined. This is where our simulation experiments could contribute: test different potentials and see which produces the correct physics.

### 6. Three Lepton Families — Confirmed

Duda confirms the mechanism: 3 distinguishable axes in 3D → 3 families of topological defects. Each is a hedgehog of one axis. Same electric field (same topology), different mass (different axis energy). The heavier ones can decay to the lightest by field rotation — exactly electron/muon/tau behavior.

---

## DUDA'S SECOND REPLY — KEY INSIGHTS

### 7. Sine-Gordon Equation as Starting Point

Duda recommends the **Sine-Gordon equation** as the entry point for stable massive particles:

```text
∂²φ/∂t² - c²∂²φ/∂x² + (m²c²/ℏ²)·sin(φ) = 0
```

This equation is remarkable because it produces:

- **Soliton (kink) solutions**: stable, localized, particle-like structures that maintain their shape during propagation — exactly the "stable localized field configurations" Duda says we need
- **Pair creation/annihilation**: kink + anti-kink collisions can annihilate or create new pairs — directly analogous to our electron-positron annihilation simulations
- **Entire special relativity**: solitons experience Lorentz contraction and time dilation naturally — relativistic behavior emerges from the wave equation, not imposed
- **Mechanical model**: a lattice of coupled pendula reproduces Sine-Gordon physics visually (Duda's YouTube link) — the closest macroscopic analog to what we simulate

**Connection to OpenWave**: the Sine-Gordon equation is a non-linear wave equation with a sinusoidal potential `V(φ) = 1 - cos(φ)`. This sits between:

- Our free wave equation `V = 0` (Combined W-L, linear, sinc nodes)
- Smolinski's quartic `V = k·ψ⁴/4` (cubic non-linearity)
- Close's nonlinear vector equation (Dirac from elastic solid)

The Sine-Gordon potential creates TOPOLOGICALLY PROTECTED solitons — the kink connects two different vacuum states (φ = 0 and φ = 2π), and you cannot continuously deform it away. This is **exactly Duda's topological charge quantization** realized in 1D.

**Key insight**: the coupled pendulum mechanical model could be our simplest Test 1 — implement Sine-Gordon in 1D, verify kink stability, pair creation/annihilation, and Lorentz contraction. This would give us intuition before tackling the full 3D hedgehog tests.

### 8. Zitterbewegung — Electron Trembling Motion

Duda points to **Zitterbewegung** ("trembling motion") as the physical mechanism that needs to emerge: the electron oscillates at frequency `2mc²/ℏ ≈ 1.6 × 10²¹ Hz` even at rest. This is the origin of the electron's spin and magnetic moment.

- **Wikipedia**: <https://en.wikipedia.org/wiki/Zitterbewegung>
- **Experimental confirmation**: Gerritsma et al. (2010) observed Zitterbewegung in a trapped-ion simulation of the Dirac equation (<https://link.springer.com/article/10.1007/s10701-008-9225-1>)

**Connection to OpenWave**: our wave centers oscillate at f₀ = 10²⁵ Hz — but we ASSUME this frequency. Zitterbewegung says the electron's oscillation frequency is `2mc²/ℏ` — DERIVED from its mass. This connects to Duda's time crystal mechanism: the oscillation is propelled by mass, not assumed.

Our f₀ is the medium's fundamental frequency. The electron's Zitterbewegung frequency is different (10²¹ vs 10²⁵ Hz). The relationship between the two frequencies (ratio ~ 10⁴) may relate to the electron's K=10 structure or the fine-structure constant.

### 9. Faber's 4D Approach → Automatic Propulsion

Duda says Faber's approach (his own 4D extension to teleparallelism, arXiv:2501.04036) "automatically gives such propulsion" — meaning the Zitterbewegung oscillation emerges naturally from the 4D field dynamics without being imposed. The kink in the phi-4 field oscillates because of the curvature coupling, and this oscillation IS the electron's trembling motion.

### Duda's Recommended Learning Path

Duda is giving us a clear pedagogical sequence:

1. **Walking droplets** → intuition for QM from classical waves (Couder experiments)
1. **Sine-Gordon equation** → stable massive particles, pair creation, special relativity (1D first)
1. **Liquid crystal particle analogs** → topological charge quantization (3D hedgehogs)
1. **Zitterbewegung** → why particles oscillate (mass-driven trembling)
1. **Faber's 4D approach** → combines all of the above (LdG + Skyrme + teleparallelism)

This maps to a test sequence for OpenWave:

- **Test 1**: Sine-Gordon 1D — kink solitons, pair creation, Lorentz contraction
- **Tests 2-3**: Hedgehog energy + topological charge (3D director field)
- **Test 4**: Klein-Gordon from twist dynamics
- **Test 5**: Lagrangian derivation
- **Test 6**: Three lepton families
- **Test 7**: Close's nonlinear vector wave equation
- **Test 8**: Smolinski's non-linear Ψ³ — direct K-selectivity test on familiar scalar setup

---

## INITIAL ANALYSIS

### Who is Jaroslaw Duda?

Dr. Jaroslaw Duda is a computer scientist and mathematician at Jagiellonian University, Krakow. He is the inventor of **Asymmetric Numeral Systems (ANS)** — a family of entropy coding methods now foundational in data compression (used in zstd, Apple LZFSE, JPEG XL, Google Draco, Linux kernel). Beyond information theory, he works on particle physics models based on liquid crystal analogs and topological field theory.

- Wikipedia: <https://en.wikipedia.org/wiki/Jaros%C5%82aw_Duda_(computer_scientist)>
- Key physics paper: <https://arxiv.org/pdf/2108.07896> (Coulomb from topological charges in liquid crystal framework)
- Liquid crystal particle analogs: <https://en.wikipedia.org/wiki/Draft:Liquid_crystal_particle_analogs>

OpenWave does NOT use ANS (it's data compression, unrelated to physics simulation). However, if we ever need to compress large voxel/time-series datasets, ANS-based compressors (zstd) would be ideal.

---

## Duda's Three Challenges to OpenWave

### Challenge 1: "Show me your Lagrangian"

**The problem**: OpenWave currently operates at the PDE / field evolution level, testing candidate wave equations empirically (5 wave equations tested in M3, selected by best K=10 stability). We do NOT have a canonical Lagrangian from which our wave equations are derived.

**Why it matters**:

- A Lagrangian DERIVES the equations of motion via Euler-Lagrange equations — instead of testing 5+ equations empirically, the correct one would emerge from first principles
- Noether's theorem automatically provides conservation laws (energy, momentum, angular momentum) — we currently check these manually
- It constrains which PDEs are physically valid — not all PDEs come from a valid Lagrangian
- It connects to established physics formalisms (QFT, GR, Standard Model) — making OpenWave's results interpretable by the broader physics community
- It enables variational methods and stability analysis — deriving which configurations are stable minima of the action

**Our current state**: the Combined Wolff-LaFreniere equation `ψ = 2A·sin(kr/2)·cos(kr/2-(ωt+φ))/r` was selected empirically for K=10 tetrahedron stability. Its phase component (`A·sin(kr)/r`) *is* a free-wave solution — it's equivalent to a pure standing wave `2A·sin(kr)·cos(ωt+φ)/(kr)`, which decomposes into in-wave + out-wave spherical solutions of `∂²ψ/∂t² = c²∇²ψ`. **But its quadrature component `A·(1-cos(kr))/r·sin(ωt+φ)` is NOT a free-wave solution** — verified by sympy in Experiment 5. The full Combined W-L produces a residual `-A·c²k²·sin(ωt+φ)/r ≠ 0` when substituted into the d'Alembertian, meaning the formula implicitly includes an external source term or comes from a non-trivial Lagrangian we have yet to identify. This is an important empirical fact: the LaFreniere quadrature piece is a modeling choice, not a free-field solution. Boundary conditions, source terms, and nonlinear extensions are also modeling choices, not derived from a Lagrangian.

### Challenge 2: "Without charge quantization your electron explodes"

**The problem**: in our model, `charge_sign = cos(source_offset)` is imposed as ±1. Nothing in the wave physics FORCES charge to be integer-valued. A standing wave configuration could in principle have any phase offset (not just 0 or π), which would give fractional charge. There's no topological protection preventing the particle from dissipating.

**Why it matters**:

- Quantized modes from standing waves explain orbit quantization (electron shells), but NOT charge quantization
- Without a mechanism that FORCES charge to be exactly ±1 (or 0), the model can't explain why all electrons have exactly the same charge
- A charged particle held together only by standing wave interference has no protection against perturbations that would spread the charge out — it would "explode"
- This is distinct from K-selectivity (which K values form stable particles) — even if K=10 is stable, the charge on that particle needs to be quantized

**Our current state**: charge is imposed via source_offset (0 = positron, π = electron). Phase 1 confirmed charge is NOT emergent from wave interference in the scalar model. Phase 1c found L→T spin conversion distinguishes charges but doesn't quantize them. This is an unsolved fundamental problem.

### Challenge 3: "Particles as topological defects, not standing waves?"

**The problem**: Duda argues particles should be understood as topological defects (hedgehogs, vortices) in a field, similar to defects in liquid crystals — not as standing wave configurations.

**Key distinction**:

- **Standing waves** (our model): particles are localized interference patterns. Stability comes from constructive interference at nodes. Can be disrupted by perturbation (as we demonstrated — K=10 breaks under perturbation)
- **Topological defects** (Duda's model): particles are field configurations with integer winding numbers. Stability comes from topology — you CANNOT continuously deform a hedgehog into a smooth field. Protected against ALL perturbations, not just small ones

**Our current state**: we model particles as standing wave structures. But several of our observations hint at topological structure:

- The dual vs non-dual geometry (K=1,8,20 neutral vs K=10,28,50 charged) has a topological flavor — it's about the SYMMETRY of the configuration, not its amplitude
- Smolinski's winding number classification of particles (sphere S² vs torus T²) is explicitly topological
- Butto's vortex electron model (toroidal flow, spin-1/2 from differential rotation) is a topological vortex

---

## Duda's Framework — The Mathematics

### The Field: Director Field n(x)

Instead of a scalar displacement ψ(r,t), the fundamental object is a **unit vector field** `n(x)` at every point in space — a "director" that specifies the local orientation of the medium. This is the order parameter of a nematic liquid crystal.

More generally, Duda uses a **symmetric matrix field** `M(x) = O·D·O^T` where O is orthogonal (rotation) and `D = diag(λ₁, λ₂, λ₃)` encodes the ellipsoid shape — directly analogous to our 6-phasor ellipse model in M4.

### The Lagrangian (LdGS model)

```text
L = -Σ F_μν_αβ · F^μν_αβ  -  V(M)

where:
  F_μν_αβ = [∂_μ M, ∂_ν M]_αβ     (matrix commutator of spacetime derivatives)
  V(M) = Landau-de Gennes potential  (prefers specific eigenvalue set)
```

The potential in Landau-de Gennes form:

```text
V_LdG = a·Tr(M²) - b·Tr(M³) + c·(Tr(M²))²
```

This has the same structure as the Higgs potential — it regularizes singularities to finite energy. The minima of V define the "vacuum state" (undisturbed medium).

### Topological Charge = Winding Number

A **hedgehog** defect is a field configuration where `n(x) = x̂` (director points radially outward from a center). The topological charge is the winding number:

```text
Q = (1/4π) ∮_S (∂_u n × ∂_v n) · n  du dv
```

This integral over any closed surface S surrounding the defect returns an **integer** — guaranteed by the Gauss-Bonnet theorem. Q = +1 for hedgehog (positive charge), Q = -1 for anti-hedgehog (negative charge). No fractional charges possible.

This is Duda's answer to charge quantization: define the electric field as the curvature of the director field, then Gauss's law `∮ E·dA = Q/ε₀` becomes the Gauss-Bonnet theorem counting topological winding numbers.

### Coulomb Potential from Topology (no sinc!)

From Duda's paper (arXiv:2108.07896, Figure 2): integrating the total field energy (Hamiltonian) for two hedgehog defects at various distances gives:

```text
E(d) ≈ 1590 + 25.6/d
```

This is a clean **1/r Coulomb potential** — no sinc oscillation, no cos(k·Δr) barriers. The Coulomb interaction emerges from the **global topological structure** of the field, not from pairwise wave interference.

**Why no sinc**: the sinc oscillation in our model comes from coherent monochromatic wave interference. The topological approach doesn't use wave superposition for the Coulomb mechanism — it uses the field curvature around defects. Different mathematical structure, different result.

### Three Forces from Three Degrees of Freedom

The vacuum state `D = diag(g, 1, δ, 0)` with `g >> 1 >> δ > 0` separates curvature into three energy scales:

| Degree of freedom | Energy scale | Physics | Governing equation |
| --- | --- | --- | --- |
| Tilts (field direction) | High | Electromagnetism | Maxwell equations |
| Twists (field rotation) | Low | Quantum phase | Klein-Gordon equation |
| Boosts (0th time axis, 4D) | Very low | Gravity | GEM (2nd Maxwell set) |

**EM from tilts**: the two tilt degrees of freedom of the director field satisfy Maxwell's equations. Electric field = tilt curvature. Magnetic field = tilt curl. Charge = topological winding number of tilts.

**QM from twists**: the twist degree of freedom (rotation of the director around its own axis) satisfies the **Klein-Gordon equation**. This is the quantum mechanical wave equation for a massive scalar field. The twist phase IS the quantum phase ψ = e^(iφ).

**Gravity from boosts**: extending from 3D SO(3) to 4D SO(1,3) by adding a time axis (Einstein's teleparallelism), perturbations of the 0th axis satisfy gravitoelectromagnetic (GEM) equations — a second copy of Maxwell's equations with G replacing 1/ε₀.

---

## Connection to OpenWave's Existing Framework

| OpenWave concept | Duda's equivalent | Notes |
| --- | --- | --- |
| M4 vector displacement (ψ_x, ψ_y, ψ_z) | Director field n(x) | Both are 3D vector fields on a grid |
| 6-phasor ellipse (R, Φ per axis) | Order parameter tensor Q_ij | Both encode ellipsoidal shape + orientation |
| L→T spin conversion | Tilt → twist conversion | Both are mode coupling in the vector field |
| Standing wave particle | Topological defect (hedgehog) | Different stability mechanism (interference vs topology) |
| Phase offset (0 vs π) | Winding number (+1 vs -1) | Imposed ±1 vs topologically quantized ±1 |
| Smolinski's Degraded EMC Wall | Defect core boundary | Both separate internal dynamics from external field |
| K=10, K=28, K=50 hierarchy | Lepton families from biaxial hedgehog | 3 distinguishable axes → 3 mass scales |
| Sinc nodes → lock-in | Not present in topological Coulomb | Fundamentally different force mechanism |
| Energy `E = ρV(fA)²` | LdG free energy `F = elastic + V(Q)` | Both compute energy from field configuration |
| `F = -∇E` force | Same: `F = -∇H` from Hamiltonian | Same variational principle |

**Key compatibility**: Duda's framework does NOT contradict standing waves — it adds a **topological layer** on top. The standing wave structure could be the internal dynamics of a topological defect. The defect's core has complex wave dynamics (our near-field), while the far-field is governed by topology (giving clean Coulomb). This maps to Smolinski's two-domain model: Energy Domain (internal, toroidal, non-linear) vs EMC Domain (external, spherical, topological).

---

## What We Can Test in OpenWave

### Test 1: Sine-Gordon 1D Solitons (intuition foundation)

Added after Duda's second reply. Simplest entry point into the topological soliton framework — 1D Sine-Gordon equation with kink/anti-kink solutions. Build intuition before tackling 3D hedgehogs.

```text
Setup:
  - 1D grid with field φ(x,t)
  - Sine-Gordon equation: ∂²φ/∂t² - c²∂²φ/∂x² + (m²c²/ℏ²)·sin(φ) = 0
  - Initial condition: kink solution φ(x,0) = 4·arctan(exp((x-x₀)/L))
  - Time evolution: leapfrog or RK4
  - Tests:
    1. Verify kink stability (no dispersion over many periods)
    2. Kink + anti-kink collision: pair annihilation or pass-through
    3. Two kinks moving toward each other: Lorentz contraction visible?
    4. Compute kink rest energy E = 8·m·c² and verify
    5. Mechanical analog: implement coupled pendula chain, verify same physics

Reference: coupled pendula visualization https://www.youtube.com/watch?v=nl5Qq5kUbEE
```

**Success criteria**:

- Kink propagates without dispersing (topologically stable)
- Kink + anti-kink can annihilate (pair annihilation analog)
- Moving kinks show Lorentz-like contraction (relativistic kinematics from wave equation)

**Why this matters**: Sine-Gordon is the simplest equation that produces all of: stable particles (kinks), pair creation/annihilation, and special relativity. It's the closest 1D analog to what we want in 3D. If we can reproduce this, we have a solid foundation for the 3D hedgehog tests.

**Infrastructure**: standalone numpy script in `sandbox_phase2_lagrangian/`, ~200 lines. 1D grid, simple PDE evolution. No M3/M4 refactor.

### Test 2: Hedgehog Energy vs Distance (Coulomb verification)

Implement two hedgehog defects in a 3D director field on M4's grid. Compute total field energy as a function of separation. Verify E(d) ~ const + C/d (1/r Coulomb) without sinc oscillation.

```text
Setup:
  - 3D grid (existing M4 infrastructure)
  - Director field n(x) = unit vector at each voxel
  - Hedgehog 1 at position r₁: n = (x-r₁)/|x-r₁|
  - Hedgehog 2 at position r₂: n = (x-r₂)/|x-r₂|
    (or anti-hedgehog: n = -(x-r₂)/|x-r₂|)
  - Relax field to minimize elastic energy (gradient descent or conjugate gradient)
  - Measure total energy at separations from 2λ to 20λ
  - Compare: E(d) vs 1/r reference

Hamiltonian: H = Σ_voxels [ K₁(∇·n)² + K₃(n×(∇×n))² ]  (Frank elastic energy)
  where K₁ = splay, K₃ = bend elastic constants
  (one-constant approximation: K₁ = K₃ = K → H = K·Σ|∇n|²)
```

**Success criterion**: E(d) fits const + C/d with R² > 0.99 across all separations, no oscillation.

### Test 3: Topological Charge Quantization

Compute the winding number integral on closed surfaces surrounding defects. Verify it returns integers only.

```text
Q = (1/4π) ∮_S (∂_u n × ∂_v n) · n  du dv

Discretized on a spherical mesh surrounding the defect center.
Should return Q = +1 for hedgehog, Q = -1 for anti-hedgehog, Q = 0 for no defect.
Test: perturb the field (noise, deformation) and verify Q remains integer.
```

**Success criterion**: Q returns integer ±1 regardless of surface shape or field perturbation.

### Test 4: Klein-Gordon from Twist Dynamics

In the uniaxial limit (one distinguished axis), evolve the twist degree of freedom dynamically. Verify it satisfies the Klein-Gordon equation `(∂²/∂t² - c²∇² + m²)φ = 0` where φ is the twist angle.

```text
Setup:
  - Director field with small twist perturbation
  - Evolve using Euler-Lagrange equations from LdG Lagrangian
  - Measure dispersion relation ω²(k)
  - Klein-Gordon predicts: ω² = c²k² + m²  (massive dispersion)
  - Mass m related to the potential curvature (eigenvalue splitting δ)

Compare with our standing wave: ω = ck (massless dispersion).
The mass term would come from the Landau-de Gennes potential.
```

**Success criterion**: dispersion relation fits ω² = c²k² + m² with measurable mass gap.

### Test 5: Can Our Combined W-L Equation Be Derived from a Lagrangian?

Check whether our empirically-selected wave equation is the Euler-Lagrange equation of a known Lagrangian. If yes, that Lagrangian is our candidate. If no, our equation may not conserve energy properly.

```text
Standard wave equation Lagrangian:
  L = ½(∂ψ/∂t)² - ½c²(∇ψ)²
  → Euler-Lagrange: ∂²ψ/∂t² = c²∇²ψ  ✓  (our Combined W-L is a solution)

With Landau-de Gennes potential V(ψ):
  L = ½(∂ψ/∂t)² - ½c²(∇ψ)² - V(ψ)
  → Euler-Lagrange: ∂²ψ/∂t² = c²∇²ψ - V'(ψ)

With Smolinski's cubic non-linearity:
  V(ψ) = k·ψ⁴/4  →  V'(ψ) = k·ψ³
  → Euler-Lagrange: ∂²ψ/∂t² = c²∇²ψ - k·ψ³  (Smolinski's soliton equation!)
```

**Key insight**: Smolinski's non-linear wave equation `(∂²/∂t² - c²∇²)Ψ + k·Ψ³ = 0` DOES come from a valid Lagrangian with quartic potential. This connects our existing non-linear research (Phase 1d) directly to the Lagrangian framework Duda is asking for.

### Test 6: Three Lepton Families from Biaxial Hedgehog

A biaxial nematic hedgehog (3 distinguishable axes instead of 1) should produce 3 types of hedgehog configurations with the same topological charge but different energy/mass. Compute their energy ratios and compare to electron/muon/tau mass ratios.

```text
Setup:
  - Biaxial order parameter: D = diag(λ₁, λ₂, λ₃) with λ₁ ≠ λ₂ ≠ λ₃
  - Three hedgehog types: defect aligned with axis 1, 2, or 3
  - Each has Q = ±1 (same charge) but different energy (different mass)
  - Compute energy of each type
  - Compare ratios to: m_e = 0.511 MeV, m_μ = 105.7 MeV, m_τ = 1776.8 MeV

EWT comparison: K=10 (electron), K=28 (muon), K=50 (tau).
Duda's approach: same topology, different axis orientation → different mass.
```

---

## Implementation Feasibility

All 8 tests are **doable without refactoring M4 or M3**. The director field / liquid crystal / non-linear scalar physics is fundamentally different from our wave propagation — they live as standalone numpy research scripts (same pattern as `sandbox_phase1_vector/`, `sandbox_phase1_scalar/`). Scripts go in `research_hub/sandbox_phase2_lagrangian/`.

| Test | Approach | Effort | M4 refactor? |
| --- | --- | --- | --- |
| 1. Sine-Gordon 1D solitons | Standalone numpy. 1D Sine-Gordon PDE, kink/anti-kink creation, pair annihilation, Lorentz contraction | Low-Medium (~200 lines) | No |
| 2. Hedgehog energy vs distance | Standalone numpy. 3D grid with unit vector field `n(x)`, Frank elastic energy `H = K·Σ\|∇n\|²`, gradient descent relaxation, energy vs separation sweep | Medium (~300 lines) | No |
| 3. Topological charge quantization | Function inside Test 2 script. Discretize sphere, compute winding number integral, verify integer under perturbation | Low (~50-100 lines on top of Test 2) | No |
| 4. Klein-Gordon from twist | Standalone numpy. Initialize twist perturbation, evolve with LdG Euler-Lagrange (leapfrog), FFT for dispersion ω(k) | Medium (~200 lines) | No |
| 5. Lagrangian derivation | Math + optional sympy. Verify Combined W-L is Euler-Lagrange of standard wave Lagrangian. Check Smolinski Ψ³ connection | Low (pen-and-paper + ~50 lines sympy) | No |
| 6. Three lepton families | Extension of Test 2 with 3x3 matrix field. Most complex — 6 independent components per voxel, more complex energy functional | Medium-High (~400 lines) | No |
| 7. Close's nonlinear vector wave eq | Standalone numpy. 3D vector spin density field, seed with spherical harmonic, time-domain PDE evolution | Medium (~300-400 lines) | No |
| 8. Smolinski's non-linear Ψ³ | Standalone numpy. 3D scalar grid, time-domain PDE with `-k·Ψ³` term, multiple WCs, direct K-selectivity test | Medium (~300 lines) | No |

**Decision strategy**: Run all 8 numpy script tests first. After reviewing the results, select the winning equation/approach and implement it in M4 (or a new M5 method if the architecture change is significant).

**Recommended order**: Test 1 (Sine-Gordon 1D, build intuition) → Tests 2+3 together (hedgehog energy + topological charge — Coulomb from topology vs our sinc Coulomb) → Test 5 (math, informs everything) → Test 8 (Smolinski Ψ³, direct K-selectivity test on familiar scalar M3-like setup) → Test 4 (Klein-Gordon dynamics) → Test 6 (three lepton families, most complex) → Test 7 (Close's equation).

**Future scaling**: if a test validates and we want GPU-accelerated 3D visualization, that would mean new Taichi kernels (possibly a new method M5 alongside M3/M4) — but that's a future decision based on results, not a prerequisite.

---

## Compatibility Assessment: Standing Waves vs Topological Defects

**Not mutually exclusive.** The two frameworks may describe different aspects of the same physics:

| Aspect | Standing waves (EWT/OpenWave) | Topological defects (Duda) |
| --- | --- | --- |
| What IS the particle? | Interference pattern of in/out waves | Winding number of director field |
| What stabilizes it? | Constructive interference at nodes | Topological protection (integer winding) |
| What is charge? | Phase offset (imposed ±1) | Winding number (quantized by topology) |
| What is Coulomb force? | Sinc-modulated interference (problematic) | Field curvature around defect (clean 1/r) |
| What is spin? | L→T wave conversion | Handedness of the defect |
| What is mass? | Energy in standing waves (E=mc²) | Total field energy of the defect |
| Near-field behavior | Sinc wells → lock-in (strong force) | Core structure of defect |
| Far-field behavior | Sinc barriers (problematic) | Clean 1/r Coulomb |

**Possible unification**: the standing wave structure IS the internal dynamics of a topological defect. The near-field (where sinc wells create lock-in) is the defect core. The far-field (where we need clean Coulomb) is governed by the topological structure of the surrounding field. This would explain:

- Why sinc is correct for near-field (lock-in, annihilation) but fails for far-field (Coulomb)
- Why charge needs to be integer (topology) but phase offset determines matter/antimatter
- Why K-selectivity exists (only certain winding configurations are topologically stable)

---

## DUDA EXPLICITLY CONFIRMS: Both Topology AND Standing Waves Needed

A critical clarification from Duda's first reply (line 194):

> *"Thanks, sure we need both topological type quantization (e.g. charge, spin), but also different standing wave type - especially for orbit quantization, see e.g. [Fort/Couder 2010, Perrard 2014]"*

This is not "topology REPLACES standing waves" — it's "topology PLUS standing waves." Two different physical phenomena requiring two different mechanisms:

| Quantization type | Mechanism | What it explains |
| --- | --- | --- |
| **Charge / spin** | Topology (winding numbers) | Why charge is exactly ±1, why spin is ±½ (integer-valued conserved quantities) |
| **Orbits / energy levels** | Standing waves (interference) | Why electron orbitals have discrete radii, Couder droplet eigenstates (discrete bound states) |

Duda points to the Couder droplet experiments (Fort 2010, Perrard 2014) as the proof: those experiments show **classical wave-particle systems exhibiting orbit quantization from standing wave dynamics**, not from topology. So orbit quantization doesn't need quantum mechanics OR topology — it emerges from wave interference alone.

**Why this matters for OpenWave**: this validates our entire approach so far. The standing wave physics we've been demonstrating (lock-in wells, sinc nodes, K=10 tetrahedron) IS doing useful work — it's the right mechanism for one of the two quantization types. We just need to ADD the topological layer to handle charge.

> Duda's message (paraphrased): "Your standing waves are correct for what they explain. You just need to add topology for what they CAN'T explain."

This is additive, not replacement. OpenWave's current M3 results stay valid. The Lagrangian sub-project is about extending — not abandoning — the standing wave foundation.

---

## PHYSICAL MAPPING: WL Standing Waves ↔ Hedgehog Topology

How exactly do Wolff-LaFreniere standing waves and Duda's hedgehog defects relate? They describe **different layers of the same particle**.

### Two Layers of One Particle

| Layer | Quantity | Behavior | Time-dependence | What it explains |
| --- | --- | --- | --- | --- |
| **Topology** (Duda) | Vector direction `n(x)` | Field winds around center | Static | Charge, winding number, integer quantization, topological protection |
| **Wave** (Wolff-LaFreniere) | Scalar amplitude `ψ(r,t)` | In-wave + out-wave interfere | Oscillates at ω | Lock-in wells, energy storage, far-field traveling waves, orbit quantization |

```text
TOPOLOGY LAYER (Duda)              WAVE LAYER (Wolff-LaFreniere)
─────────────────────             ─────────────────────────────
n(x) = r̂  (hedgehog)              ψ(r,t) = sin(kr)·cos(ωt)/r

Static structure                   Dynamic oscillation
WHERE the defect is                HOW it vibrates
WHAT kind of charge (winding)      ENERGY content (amplitude)
Topologically protected            Lock-in wells from interference
Cannot be deformed away            Can interact via wave interference
```

Both apply to the same particle simultaneously. The hedgehog provides the topological skeleton, and the standing waves are the dynamical content living on top of it.

### Visualizing the Combination

A hedgehog defect has director arrows pointing radially outward. Now imagine those arrows OSCILLATING in length around their equilibrium radial directions. The oscillation pattern around the hedgehog is what produces the WL standing wave structure:

```text
Hedgehog (static):       Hedgehog + WL oscillation (dynamic):

    ↑                         ↑       (long arrow at antinode)
  ↖ ↑ ↗                     ↖ ↑ ↗     (medium)
  ←•→                       ← • →     (short arrow at node)
  ↙ ↓ ↘                     ↙ ↓ ↘     (medium)
    ↓                         ↓       (long arrow)

Topology only             Topology + standing wave amplitude
```

The arrows still point radially (hedgehog), but their MAGNITUDE varies with r following the WL sinc envelope — antinodes where the arrows are long, nodes where the arrows are short.

### The In+Out Wave Connection

In WL, the standing wave is in-wave + out-wave:

- **In-wave**: spherical wave traveling toward the center (k vector points inward)
- **Out-wave**: spherical wave traveling outward (k vector points outward)

For a longitudinal wave, the **displacement** is along the propagation direction. So:

- In-wave granules oscillate **inward** (negative radial direction)
- Out-wave granules oscillate **outward** (positive radial direction)

When you superpose them, the time-averaged DIRECTION pattern of granule displacement IS radial — that's a hedgehog-like structure. The standing wave's spatial envelope (sinc nodes) is overlaid on a fundamentally radial vector field.

**Key insight**: WL's in+out spherical waves naturally create a hedgehog topology in their displacement field, and the standing wave nodes are the secondary structure on top. We just never noticed because M3's scalar field collapses the directional info into magnitude.

### Why M3 Can't See This

Our M3 stores a **scalar** ψ at each voxel — there's no direction, so there's no hedgehog to have. M3 has the wave layer but is missing the topology layer entirely. M4 stores a **vector** (ψ_x, ψ_y, ψ_z) at each voxel — this CAN represent both the direction (hedgehog) AND the oscillation magnitude (WL). That's why Duda's framework requires M4 architecture.

### Practical Implication for Test 2

When we implement Test 2 (hedgehog energy vs distance), we have a choice:

1. **Pure static hedgehog**: seed `n(x) = r̂` everywhere, relax field, measure energy vs separation
2. **Hedgehog + WL dynamics**: seed `n(x) = r̂` with WL in+out spherical wave structure superposed (oscillating amplitude on top of radial direction), then measure

**Option 2 may be the unification**: not "topology vs waves" but "topology IS the geometry, waves ARE the dynamics" — combining both gives us:

- **Static hedgehog topology** → integer winding, charge quantization, Coulomb 1/r
- **WL standing wave dynamics** → lock-in wells, near-field interactions, energy

If this works, M3's K=10 tetrahedron lock-in physics would carry over directly to the topological framework — not as a replacement but as the dynamical content of the topological defect.

---

## Open Questions

1. Can standing wave configurations be reinterpreted as topological defects? What is the winding number of a K=10 tetrahedron?
1. Does the Landau-de Gennes potential reproduce the energy well structure we observe in M3?
1. Can the Klein-Gordon equation from twist dynamics reproduce our standing wave particle?
1. How does the topological framework handle K-selectivity? Do only certain winding numbers correspond to stable defects?
1. Is the sinc oscillation (near-field lock-in) compatible with topological Coulomb (far-field)? Can both coexist in one framework?
1. How does Duda's 4D extension (teleparallelism → GEM) connect to our gravity-from-spin-deficit model?
1. Can we implement the LdGS Lagrangian on our existing M4 Taichi grid infrastructure?

---

## References

Duda's papers and references:

- Duda, J. (2021). "Four-dimensional understanding of quantum mechanics and Bell violation." <https://arxiv.org/pdf/2108.07896> — core paper: Coulomb from topological charges (Fig. 2), Klein-Gordon derivation around electron (Fig. 9), LdGS Lagrangian
- Duda, J. (2025). "Time crystal phi-4 kinks by curvature coupling as toy model for mechanism of oscillations propelled by mass." <https://arxiv.org/pdf/2501.04036> — derives WHY resting particles oscillate (time crystal mechanism), phi-4 potential with curvature coupling
- Duda, J. (2026?). Running coupling / regularization for Coulomb. <https://www.mdpi.com/2076-3417/16/2/1030> — (URL may be pre-publication or incorrect, verify with Duda)

Orbit quantization in classical wave-particle systems (cited by Duda):

- Fort, E., Eddi, A., Boudaoud, A., Moukhtar, J., Couder, Y. (2010). "Path-memory induced quantization of classical orbits." PNAS 107(41). <http://www.pnas.org/content/107/41/17515.short> — bouncing droplets on vibrating fluid generate quantized orbital states from wave-field memory
- Perrard, S., Labousse, M., Miskin, M., Fort, E., Couder, Y. (2014). "Self-organization into quantized eigenstates of a classical wave-driven particle." Nature Communications 5:3219. <https://www.nature.com/articles/ncomms4219> — walking droplets in harmonic potential self-organize into discrete stable orbits with quantized angular momentum

Sine-Gordon, Zitterbewegung, and soliton physics (cited by Duda, second reply):

- Sine-Gordon equation: <https://en.wikipedia.org/wiki/Sine-Gordon_equation> — soliton kinks, pair creation/annihilation, special relativity from wave equation
- Coupled pendula mechanical model (Sine-Gordon visualization): <https://www.youtube.com/watch?v=nl5Qq5kUbEE>
- Zitterbewegung (electron trembling motion): <https://en.wikipedia.org/wiki/Zitterbewegung> — electron oscillates at 2mc²/h even at rest
- Zitterbewegung experimental confirmation: <https://link.springer.com/article/10.1007/s10701-008-9225-1>

Robert Close:

- Close, R.A. (2025). "Plane Wave Solutions to a Proposed 'Equation of Everything'." Foundations of Physics 55, 27. <https://doi.org/10.1007/s10701-025-00839-0> — Dirac equation from classical elastic solid waves, Lagrangian with classical interpretation, nonlinear vector wave equation for spin density

Wikipedia / background:

- Liquid crystal particle analogs: <https://en.wikipedia.org/wiki/Draft:Liquid_crystal_particle_analogs>
- Landau-de Gennes theory: <https://en.wikipedia.org/wiki/Landau%E2%80%93de_Gennes_theory>
- Teleparallelism: <https://en.wikipedia.org/wiki/Teleparallelism>
- Gravitoelectromagnetism: <https://en.wikipedia.org/wiki/Gravitoelectromagnetism>
- String hadronization: <http://www.scholarpedia.org/article/Parton_shower_Monte_Carlo_event_generators#String_model>

## INTEGRATING THE LAGRANGIAN CONCEPT IN OPENWAVE

## Impact on Force & Motion (xforce_motion.py) — None

Lagrangian mechanics does NOT replace `F = -∇E` and `F = ma`. It's the layer **above** them that *derives* them.

```text
LAGRANGIAN (top — defines the physics)
    L = T - V = ½m·v² - V(x)
           │
           ▼
EULER-LAGRANGE EQUATION (variational calculus)
    d/dt(∂L/∂v) = ∂L/∂x
           │
           ▼
NEWTON'S LAWS (what falls out)
    m·a = -∇V  →  F = -∇E  →  F = ma
           │
           ▼
LEAPFROG INTEGRATOR (numerical solution)  ← xforce_motion.py lives HERE
    v(t+½dt) = v(t-½dt) + a·dt
    x(t+dt)  = x(t) + v(t+½dt)·dt
```

What we have now (xforce_motion.py): we compute `E = ρV(fA)²` from the wave field, then `F = -∇E` via finite differences, then integrate with leapfrog. This is correct Newtonian mechanics. The force-motion script stays.

What a Lagrangian adds: it tells us *what E should be*. Right now we chose `E = ρV(fA)²` from EWT reasoning and our wave equation empirically. A Lagrangian would derive both:

1. The wave equation itself (what ψ does in space and time)
1. The energy functional (what E looks like)
1. Conservation laws for free (energy, momentum — via Noether's theorem)

| What | Current (no Lagrangian) | With Lagrangian |
| --- | --- | --- |
| Wave equation | Chosen empirically (Combined W-L) | Derived from Euler-Lagrange |
| Energy formula | `E = ρV(fA)²` (EWT postulate) | Derived from Hamiltonian `H = T + V` |
| Force | `F = -∇E` (correct) | Same — `F = -∇E` falls out automatically |
| Integration | Leapfrog (correct) | Same — leapfrog is still the numerical method |
| Conservation | Checked manually | Guaranteed by Noether's theorem |

So `xforce_motion.py` doesn't change. What changes is the **justification** for why we compute energy the way we do, and the **wave equation** that feeds into it. The Lagrangian sits upstream of everything in that script.

The real payoff: if Duda's LdG Lagrangian is correct, it would tell us:

- The exact wave equation (no more testing 5 candidates)
- The exact energy functional (possibly different from `ρV(fA)²`)
- Why charge is quantized (topological invariant of the Lagrangian's symmetry group)
- What conservation laws hold (and which don't under which conditions)

The force-motion code is the engine. The Lagrangian is the fuel specification.

---

## Impact on Wave Engine (wave_engine.py) — This Is What Changes

wave_engine.py is exactly what a Lagrangian would affect. This is the script that lives at the Euler-Lagrange level:

```text
LAGRANGIAN
    L = ½(∂ψ/∂t)² - ½c²(∇ψ)² - V(ψ)
           │
           ▼
EULER-LAGRANGE  →  wave_engine.py lives HERE
    ∂²ψ/∂t² = c²∇²ψ - V'(ψ)
           │
           ▼
FORCE & MOTION  →  xforce_motion.py lives HERE
    F = -∇E,  leapfrog
```

Right now, wave_engine.py has 5 wave equations chosen empirically — and most of the research effort has been finding which one "works best." A Lagrangian would tell us which one is *correct*:

| Lagrangian | V(ψ) | Euler-Lagrange → wave equation | What it gives |
| --- | --- | --- | --- |
| Free wave | 0 | `∂²ψ/∂t² = c²∇²ψ` | Our Combined W-L is a solution. But no mass, no soliton stability |
| Quartic (Smolinski) | `k·ψ⁴/4` | `∂²ψ/∂t² = c²∇²ψ - k·ψ³` | Soliton stability, non-linear. Changes spatial structure from pure sinc |
| LdG (Duda) | `a·Tr(M²) - b·Tr(M³) + c·(Tr(M²))²` | Maxwell + Klein-Gordon + GEM | Topological charges, Coulomb, mass gap |

**What would change in wave_engine.py**:

- The oscillator equation itself — instead of `2A·sin(kr/2)·cos(kr/2-(ωt+φ))/r`, the spatial function could be different (non-sinc from the V(ψ) potential term)
- The energy computation — instead of `E = ρV(fA)²`, the Hamiltonian `H = ½(∂ψ/∂t)² + ½c²(∇ψ)² + V(ψ)` would define energy
- The phasor superposition might not apply — if V(ψ) is non-linear (like ψ³), superposition breaks and we'd need actual time-stepping (like M2's Laplacian mode)
- The signed envelope logic — topology could replace it entirely (charge from winding number, not from imposed sign)

**What would NOT change**:

- Grid infrastructure (Taichi fields, voxel layout)
- Rendering / visualization
- xforce_motion.py (`F = -∇E`, leapfrog)
- xparameters / experiment design
- The launcher

**The key tension**: our current wave_engine.py uses **analytical phasor precomputation** (fast, exact, GPU-friendly). A non-linear Lagrangian (ψ³ or LdG) would require **time-domain PDE evolution** (like M2's Laplacian mode) because superposition doesn't hold for non-linear equations. That's not a refactor — it's a different computation strategy. We already have both patterns in the codebase (M2 = PDE evolution, M3 = analytical phasor).

So the Lagrangian tests in `sandbox_phase2_lagrangian/` would first validate *which* Lagrangian is right, and only then would we port the winning equation into wave_engine.py — potentially as a new method (M5) rather than modifying M3.

---

## Impact on Base Wave Architecture — M2 vs M3 Philosophy

The Lagrangian / topological framework has a strong preference for HOW the background field is modeled. This maps differently to our existing methods:

### The Vacuum State IS the Base Wave

In Duda's liquid crystal model, the **vacuum state** (all directors aligned, minimum of V(M)) is the base wave equivalent. It's not empty space — it's a fully ordered field. A topological defect (hedgehog = particle) is a deformation OF this field. **The hedgehog cannot exist without the background field** — a winding number only makes sense relative to a vacuum state.

```text
Our base wave (M1/M2):  ψ_base(x,t) = A₀·cos(kx)·cos(ωt)  ← oscillating wave
Duda's vacuum:          n(x) = ẑ  (everywhere aligned)       ← static ground state

Perturbations of vacuum → propagating waves (Klein-Gordon)
Topological defects     → particles (hedgehogs)
Time crystal mechanism  → WHY defects oscillate (mass-driven, not assumed)
```

The oscillations we assume in M3 (f₀ = 10²⁵ Hz) would need to be *derived* from the time crystal phi-4 kink mechanism, not assumed. The base wave isn't waves from all matter in the universe (EWT concept) — it's the ground state of the field's potential.

### How Each OpenWave Method Maps

| Concept | M2 (Laplacian) | M3 (Wolff-LF) | Duda (LdG) |
| --- | --- | --- | --- |
| Background field | Standing wave from reflections | None — WCs emit into void | Ordered vacuum state (LdG minimum) |
| Particle | WC disturbance in field | Sum of in+out analytical waves | Topological defect (hedgehog) |
| Oscillations | Emerge from PDE evolution | Assumed at f₀ (phasor) | Derived from time crystal / phi-4 |
| Charge | Not emergent | Imposed ±1 from phase | Winding number (topological) |
| Stability | Energy conservation | Well depth (fragile) | Topological protection (robust) |
| Force | F = -∇E from Laplacian field | F = -∇E from phasor amplitude | F = -∇H from field curvature |

### M2 Is More Compatible Than M3

**M2's philosophy**: "the base wave exists, let it self-stabilize, WCs are disturbances" — defects in an existing field. This aligns with Duda's vacuum + hedgehog model.

**M3's philosophy**: "WCs ARE the field, no background" — there is no vacuum field to have topology in. A hedgehog requires a field to deform. M3 emits waves into nothing — there's no ordered state to measure winding numbers against.

But M2 has its own limitation (validated in 1D): PDE/Laplacian solvers with Dirichlet BC cannot create standing waves around absorber WCs. The wave passes through without reorganizing.

### The Resolution: Duda's Vacuum Is Static, Not Oscillating

The key difference from our M1/M2 base wave: Duda's vacuum is a **static ordered state**, not an oscillating wave field. The hierarchy is:

1. **Vacuum** (ordered director field, static) — the ground state
1. **Waves** (perturbations of vacuum) — propagating disturbances, satisfy Klein-Gordon
1. **Particles** (topological defects) — non-perturbative, large deformations, protected by topology
1. **Particle oscillations** (time crystal) — mass-driven kink oscillation, derived from phi-4 potential

Our M2 base wave oscillates at f₀. Duda's vacuum doesn't oscillate — the oscillations are a *consequence* of the defect's mass interacting with the field potential (time crystal mechanism). This is more fundamental: it explains WHY particles oscillate rather than assuming it.

### Practical Implication

All 5 Lagrangian tests (sandbox_phase2_lagrangian/) require a background vacuum field. They're inherently M2-like (field exists, defects are disturbances), not M3-like (no background). The scripts will initialize a uniform director field `n(x) = ẑ` as the vacuum, then create hedgehog defects as deformations.

This does NOT invalidate M3's results — lock-in, annihilation, K=10 stability are real near-field wave physics demonstrated in M3. But it suggests the far-field (Coulomb, charge quantization) may need a fundamentally different architecture: one where particles are topological features of a background field, not self-contained wave emitters.

### Vector Field Required — M4, Not M3

Duda's model requires a **vector field** (director `n(x)` = unit vector at every point). A scalar field has no direction, no winding number, no topology.

| Requirement | M3 (scalar) | M4 (vector) | Duda needs |
| --- | --- | --- | --- |
| Displacement per voxel | 1 f32 (scalar ψ) | 3 f32 (ψ_x, ψ_y, ψ_z) | 3 f32 minimum (n_x, n_y, n_z) |
| Can represent director n(x)? | No — no direction | Yes — 3 components = unit vector | Director field is fundamental |
| Can compute winding number? | No — scalar has no topology | Yes — vector field has winding | Topological charge = winding of n(x) |
| Can distinguish L/T? | No — scalar magnitude only | Yes — radial vs perpendicular | Tilts vs twists = EM vs QM |
| Full biaxial Q_ij tensor? | No | Partial — 3 components, needs 6 | 6 independent (3x3 symmetric matrix) |

**For initial tests** (hedgehog, uniaxial nematic): M4's 3-component vector field is sufficient. The director `n(x) = (nx, ny, nz)` maps directly to M4's `displacement_am` vector field.

**For full biaxial model** (Test 6, three lepton families): needs 6 independent components per voxel (symmetric 3x3 matrix). M4 stores 3 — but our 6-phasor model (R_x, R_y, R_z, Φ_x, Φ_y, Φ_z) already stores 6 numbers per voxel encoding an ellipsoidal shape. That IS the order parameter tensor Q_ij.

**For numpy research scripts** (sandbox_phase2_lagrangian/): infrastructure doesn't matter — we allocate whatever arrays we need. But for eventual GPU implementation (M5 on Taichi), it would extend M4's vector grid infrastructure.

#### Future M5 Method, or an upgrade to M4

If the Lagrangian tests validate, a future M4/M5 method would combine:

- **M2's philosophy**: background field exists, particles are disturbances
- **M3's near-field results**: standing wave lock-in, energy well structure
- **M4's vector infrastructure**: 3-component displacement → director field, extendable to 6-component for biaxial
- **Duda's topology**: charge quantization, clean Coulomb, topological stability
- **Time crystal dynamics**: derived oscillations, not assumed
- **Close's Lagrangian**: Dirac equation from elastic solid waves, classical interpretation of every term

---

## RELATED WORK: CLOSE (2025) — DIRAC EQUATION FROM ELASTIC SOLID

### Paper

Close, R.A. (2025). "Plane Wave Solutions to a Proposed 'Equation of Everything'." Foundations of Physics 55, 27. <https://doi.org/10.1007/s10701-025-00839-0> (open access)

**Note**: the "Robert" in Duda's thread (*"Robert, have you maybe finally got charge quantization — explained why your electron doesn't explode?"*) is almost certainly **Robert Close**, the author of this paper. He is in the same "Models of particles" mailing group. Duda is challenging him on the same charge quantization problem he challenged us on.

### Core Thesis

Close derives the **Dirac equation from classical wave mechanics in an ideal elastic solid** (an aether). Spin angular momentum density in this elastic medium obeys a nonlinear vector wave equation. When factored into first-order form, it produces the Dirac equation — with every term having a clear classical physics interpretation:

- **Wave propagation** term
- **Convection** term (medium flow)
- **Rotation of the medium** → rotational kinetic energy operator
- **Rotation of wave velocity relative to the medium** → conventional potential energy operator
- **Potential energy** = half the mass term of the free electron Dirac equation
- **Electron rest energy** = twice the conventional potential energy in the elastic solid model

### Key Results

#### The Medium Is Our Medium

Close models an "ideal elastic solid" — a vacuum/aether with density, elasticity, and angular momentum. His spin density is the angular momentum of this medium. This maps directly to our spacetime medium (ρ = 3.86 × 10²² kg/m³, wave speed c). He quotes Robert Laughlin: *"Relativity actually says nothing about the existence or nonexistence of matter pervading the universe, only that any such matter must have relativistic symmetry. It turns out that such matter exists."*

#### Lagrangian and Hamiltonian with Classical Interpretation

Close constructs a Lagrangian and Hamiltonian density where each term corresponds to a classical physical quantity:

- **Hamiltonian** = total energy = rotational kinetic energy + potential energy
- **Rotational kinetic energy**: associated with rotation of the wave function WITH the medium
- **Potential energy**: associated with wave propagation and rotation of wave velocity RELATIVE to the medium
- **Intrinsic momentum**: from the Belinfante-Rosenfeld stress-energy tensor, generator of translations

This is a **candidate Lagrangian** for OpenWave — directly answering Duda's "show me your Lagrangian" challenge, but from a different angle (elastic solid rather than liquid crystal).

#### Nonlinearity → Quantized Amplitudes

Close's wave equation is **nonlinear**. He argues this is WHY amplitudes are quantized — because multiplying a solution by a constant factor does not generally yield another solution (unlike linear equations where any scalar multiple is also a solution). Soliton/breather solutions exhibit particle-like behavior.

This provides a **second mechanism for quantization** alongside Duda's topology:

- **Duda**: charge is quantized because topology forces integer winding numbers (Gauss-Bonnet)
- **Close**: amplitudes are quantized because nonlinearity prevents continuous scaling

Both may be needed: topology quantizes the charge TYPE (integer), nonlinearity quantizes the charge MAGNITUDE (fixed amplitude).

#### Vector Wave Equation for Spin Density

Close uses a **second-order vector wave equation** for spin density — scalar fields cannot represent spin. This independently confirms our conclusion that M4 (vector field) is required, not M3 (scalar).

#### Phase Shifts as Interaction Potentials

Close proposes that interaction potentials arise from **phase shifts** between wave functions. For a phase shift `δ = (mφ - ωt)` with integer m, the derivatives correspond to a magnetic vector potential. This is a different approach to our Coulomb problem:

- **Our current approach**: forces from amplitude gradients `F = -∇E`
- **Close's approach**: forces from phase shift derivatives → vector potentials

This connects to LaFreniere's core phase shift concept: the electron core creates a λ/2 phase shift from medium compression. Close formalizes this as an interaction potential.

#### Pauli Exclusion from Wave Interference

Two classical wave functions superposed produce interference cross-terms:

```text
(ψ_A + ψ_B)†(ψ_A + ψ_B) = |ψ_A|² + |ψ_B|² + ψ_A†·ψ_B + ψ_B†·ψ_A
```

Phase shifts can cancel the interference terms without changing individual magnitudes. This cancellation is mathematically equivalent to anticommutation of wave functions — the **Pauli exclusion principle emerges from classical wave mechanics**. No quantum postulate needed.

#### Couder Droplet Foundation

Close cites the same Fort/Couder bouncing droplet papers that Duda cited (orbit quantization, pilot-wave hydrodynamics, tunneling). This research community is converging on the same experimental foundation: classical wave-particle systems exhibiting quantum-like behavior.

### Connection to OpenWave

| Close's approach | Our approach | Overlap |
| --- | --- | --- |
| Ideal elastic solid medium | Spacetime medium (EWT) | Same physical model |
| Vector spin density field | M4 vector displacement | Same field type needed |
| Nonlinear wave eq → Dirac | Nonlinear ψ³ (Smolinski) | Both need nonlinearity |
| Phase shifts → potentials | Sinc wells → lock-in | Different force mechanisms (complementary) |
| Lagrangian constructed | Lagrangian needed (Duda challenge) | Close has one we can evaluate |
| Dirac equation derived | Wave equations empirical | Close further along on formalism |
| Soliton/breather particles | Standing wave particles | Similar particle concept |
| Pauli exclusion from phases | Not yet addressed | New avenue for OpenWave |

### Two Lagrangian Candidates Now

| Source | Lagrangian type | Field | Quantization mechanism | Strength |
| --- | --- | --- | --- | --- |
| **Duda** | Landau-de Gennes + Skyrme kinetic | Director field M(x) | Topological (winding numbers) | Clean Coulomb, charge quantization |
| **Close** | Elastic solid spin density | Vector spin density | Nonlinear (fixed amplitudes) | Dirac equation, classical interpretation |

Both are valid candidates. They may address different aspects: Duda's topology handles charge quantization and far-field Coulomb; Close's nonlinearity handles amplitude quantization and the Dirac formalism. A complete theory may need both.

### Test 7: Close's Nonlinear Vector Wave Equation

Close explicitly invites others to solve his equation: *"I don't know when I'll find time to work on solutions, so I invite you and others to try."* His suggestion: *"starting with a spherical harmonic linear wave solution and see what evolves."*

```text
Setup:
  - Implement Close's nonlinear vector wave equation for spin density evolution
  - 3D grid with vector spin density field s(x) = (sx, sy, sz)
  - Seed with a spherical harmonic linear wave solution (Y_l^m)
  - Evolve using time-domain PDE (leapfrog or RK4)
  - Observe: do particle-like soliton/breather structures emerge?
  - Measure: does the evolved field match Dirac bispinor plane wave solutions?

Key equation (from Close's paper):
  Second-order vector wave equation for spin density where temporal changes
  are attributable to convection, rotation, and torque density.
  Factors into first-order Dirac equation for bispinor fields.
```

**Success criterion**: spherical harmonic initial condition evolves into a localized, stable, particle-like structure (soliton/breather) rather than dispersing. If so, this is particle formation from wave dynamics — exactly what OpenWave demonstrates with standing waves, but now derived from a formal Lagrangian.

**Infrastructure**: requires vector field (M4-like), time-domain PDE evolution (M2-like), nonlinear terms. Standalone numpy script in `sandbox_phase2_lagrangian/`, ~300-400 lines. No M4 refactor needed.

### Test 8: Smolinski's Non-linear Soliton Wave Equation (Direct K-Selectivity Test)

Directly address K-selectivity avenue #2 by numerically evolving Smolinski's non-linear Soliton Wave Equation with multiple wave centers on a 3D grid. Test whether the cubic non-linearity creates K-dependent stability that pure M3 (linear) cannot.

**The equation** — from Smolinski's *MagnetismGravity v2*, Section 6.1, Eq. 18-19 (page 18):

```text
General form (Eq. 18):
  (∂²/∂t² - c²∇²) Ψ(r,t) + F(Ψ, ε_G, |ε_M|, N_ν) = 0

Example NLS stabilizing term (Eq. 19):
  F = k(|ε_M|) · Ψ³

Variables:
  Ψ(r,t)  — Soliton wave function (the WC / particle)
  c²∇²     — standard d'Alembertian (linear wave operator)
  F        — nonlinear stabilizing functional (prevents dispersion)
  ε_G      — geometric parameter (Soliton internal geometry)
  |ε_M|    — Dynamic Magnetic Deficit ≡ 1/(N·π³)  (Eq. 22)
  N_ν      — effective volume deficit (≈ 4.66 × 10⁵⁴)
  k(|ε_M|) — nonlinear elasticity coefficient (depends on magnetic deficit)
```

The nonlinearity is explicitly called out as **NLS (Nonlinear Schrödinger) form** — the canonical soliton-supporting cubic term. Smolinski states: "Stability arises from geometric self-stabilization via internal magnetic/spin coherence."

**Boundary condition**: the topology parameter K acts as the boundary condition for the wave equation (Sec 8.1.2, p.21). The Soliton radius and packing volume are defined by K — this is potentially how K-selectivity could emerge naturally from the equation.

**Radial density profile** (Sec 5.6, p.15): `ρ(r) = ρ₀ · (1 - (r/r_ν)^k)^p · Θ(r_ν - r)` — with Heaviside truncation at the Soliton boundary.

**Test setup**:

```text
Setup:
  - 3D scalar grid with field Ψ(x,y,z,t)
  - Non-linear wave equation: ∂²Ψ/∂t² = c²∇²Ψ - k·Ψ³
  - Time-domain PDE evolution (leapfrog or RK4) — NOT phasor superposition
    (superposition breaks for non-linear equations)
  - Initial condition: K wave centers placed at geometric positions
    (K=2 line, K=3 triangle, K=10 tetrahedron, etc.)
  - Start with coefficient k as a free parameter (sweep its value)
  - Observe long-time stability under perturbation

Tests to run:
  1. Single WC: does it self-stabilize into a soliton? (breather solution)
  2. Two WCs (K=2): stable or decays?
  3. Tetrahedron (K=10): stable AND more stable than K=2..9?
  4. Perturbation test: does K=10 resist deformation?
  5. Sweep k values: find the coefficient range where K=10 is uniquely stable
  6. Compare with K=3, K=5 at same k values
```

**Success criterion**: K=10 tetrahedron is stable under perturbation AT SOME k value, while K=2..9 are not. If YES, this is **K-selectivity from nonlinearity alone** — no topology needed. If NO, it confirms that topology is required (Duda's position).

**Why this test is critical**:

- Directly tests K-selectivity avenue **#2** from our open avenues list
- Minimal extension of M3 — same scalar field, same grid, just add `-k·Ψ³` term
- Computation strategy changes: phasor superposition (M3) → time-domain PDE (M2-like) because superposition breaks for non-linear equations
- Can be compared head-to-head with M3's Combined W-L results
- Connects directly to the Lagrangian framework: Smolinski's equation comes from a quartic potential Lagrangian `V(ψ) = k·ψ⁴/4`, so this is also a partial Test 5 result

**Open questions**:

- What value of `k` should we use? It depends on `|ε_M|` (the Dynamic Magnetic Deficit) — we may need to sweep
- Should `Ψ` be real scalar, or do we need complex (to match NLS form)? Start with real, test complex if real fails
- Does the K boundary condition (from Sec 8.1.2) need to be explicitly imposed, or does it emerge from the dynamics?

**Infrastructure**: standalone numpy script in `sandbox_phase2_lagrangian/`, ~300 lines. Time-domain PDE evolution (like M2 Laplacian mode) on a 3D scalar grid. No M3/M4 refactor.

#### Smolinski's F Functional is a Placeholder

An important nuance about Smolinski's equation: the `F(Ψ, ε_G, |ε_M|, N_ν)` functional in Eq. 18 is a **placeholder**. Smolinski is being honest — he doesn't fully know the form, only that SOME nonlinear stabilizing term must exist.

```text
General form (Smolinski doesn't know the exact form):
  (∂²/∂t² - c²∇²) Ψ + F(Ψ, ε_G, |ε_M|, N_ν) = 0
                      ─────────────────────
                      "some nonlinear functional of these 4 things"

Specific example (Eq. 19, NLS form):
  F = k(|ε_M|) · Ψ³
      ─────────
      simplest non-trivial case borrowed from Nonlinear Schrödinger
```

Smolinski explicitly says the `Ψ³` form "is known from NLS equations" — meaning it's the canonical soliton-supporting term. It's the **simplest non-trivial nonlinear choice** that respects the expected ±Ψ symmetry (quadratic terms are usually excluded because the equation should be invariant under Ψ → −Ψ).

**What the four variables represent**:

| Variable | Physical meaning |
| --- | --- |
| **Ψ** | The wave function itself (self-interaction) |
| **ε_G** | Soliton internal geometry parameter (shape: spherical, tetrahedral, etc.) |
| **\|ε_M\|** | Dynamic Magnetic Deficit ≡ 1/(N·π³) — vacuum medium stiffness |
| **N_ν** | Effective volume deficit ≈ 4.66 × 10⁵⁴ — how much volume the soliton displaces |

**What the simplification hides**: when Smolinski writes `F = k(|ε_M|)·Ψ³`, he's making three assumptions:

1. **Only the cubic term matters** → ignoring Ψ⁵, Ψ⁷, and higher powers
2. **ε_G and N_ν are hidden inside k** → geometry and volume deficit might be absorbed into the coefficient
3. **No mixed terms** → F doesn't have things like `ε_G · Ψ²` that couple geometry to the wave function directly

The full F could be much richer — something like:

```text
F_full = k₁(|ε_M|)·Ψ³ + k₂(ε_G)·Ψ⁵ + k₃(N_ν)·(∇Ψ)²·Ψ + ...
         ───────────   ───────────   ──────────────────
         cubic NLS     quintic       gradient-coupled
```

**Think of F as a Taylor expansion** of an unknown nonlinear function:

```text
F(Ψ, ...) ≈ a·Ψ + b·Ψ² + c·Ψ³ + d·Ψ⁴ + e·Ψ⁵ + ...
            ─     ─     ─
            linear      ← first non-trivial nonlinear term that
            (already    respects ±Ψ symmetry — Smolinski's choice
            in wave eq)
            ─
            quadratic (usually excluded by symmetry)
```

So `k·Ψ³` is the **leading nonlinear term** in a perturbative expansion. It's what you try FIRST, and only add higher-order terms if the cubic alone doesn't capture the physics.

#### What Test 8 Should Look For — Progressive Complexity

Given that Smolinski's F is a placeholder, Test 8 should proceed through a **ladder of nonlinear forms**, starting simple and adding complexity only if simpler forms fail:

**Level 1 — Simplest form** (start here):

```text
F = k · Ψ³     (constant k, no dependence on other parameters)
```

- Sweep k values across multiple orders of magnitude
- Test K=2..10 tetrahedron stability under perturbation
- **Success criterion**: K=10 is stable AT SOME k value while K=2..9 are not

**Level 2 — Spatially varying k** (if Level 1 fails):

```text
F = k(r) · Ψ³    (k depends on distance from WC, encoding ε_M or ε_G implicitly)
```

- k could depend on local energy density, local wavelength, or distance to nearest WC
- Tests whether the K-selectivity emerges when k tracks the Yee & Hauger shells
- Connects to Phase 1d variable-λ research

**Level 3 — Higher-order terms** (if Level 2 fails):

```text
F = k₁·Ψ³ + k₂·Ψ⁵
```

- Add quintic term (next odd power preserving ±Ψ symmetry)
- Tests whether the cubic is just too weak to distinguish K values

**Level 4 — Mixed geometric coupling** (if Level 3 fails):

```text
F = f(ε_G) · g(Ψ)
```

- Explicitly encode the geometric parameter ε_G (could be the K-count itself)
- Test whether the K value needs to appear directly in the equation, not just emerge

**Level 5 — Fallback: topology required**:

If no scalar nonlinear form produces K-selectivity, the conclusion is that **topology (Tests 1-2) is the right answer, not nonlinearity alone**. Smolinski's cubic term might stabilize a soliton against dispersion, but it cannot discriminate K values without topological structure.

**The value of Test 8**: regardless of which level succeeds (or if none does), we learn something important:

- **Level 1 works** → K-selectivity comes from simple cubic nonlinearity (huge win, minimal complexity)
- **Level 2-3 works** → nonlinearity is the right mechanism but needs more structure
- **Level 4 works** → geometric coupling is required, hints at vector field structure
- **No level works** → topology is required, validates Duda's topological approach over scalar nonlinearity

Each negative result narrows the search space. Each positive result gives us a concrete equation to port to M4/M5.
