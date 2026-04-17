# PHASE 3 — CONCEPT REVIEW

A conceptual walk-through of the Lagrangian / topological framework from Rodrigo's intuitive-build-up perspective. This document captures where the EWT / M1-M2-M3 physical picture has to bend to accommodate Duda's vacuum + defect framing, and where it stays the same.

**Status** (2026-04-17): the 8 sandbox experiments of Phase 3 are now complete. Every conceptual claim below that had a testable prediction has been checked numerically. The consolidated "empirical validation" scorecard is in the **[EMPIRICAL VALIDATION](#empirical-validation--what-the-sandbox-experiments-confirmed)** section near the end of this document. Core findings in one line: **topology is load-bearing, Klein-Gordon-like wave dynamics validate the mass-gap mechanism, pure scalar/vector nonlinearity without topology is insufficient**.

Related reading:

- [3_LAGRANGIAN_FRAMEWORK.md](3_LAGRANGIAN_FRAMEWORK.md) — full Lagrangian framework evaluation, email thread, Duda/Close context
- [3b_lagrangian_experiments.md](3b_lagrangian_experiments.md) — numerical experiment results
- [3c_path_to_m5.md](3c_path_to_m5.md) — M5 / Lagrangian-Wave Method implementation plan
- [0_WAVE_EQUATION.md](0_WAVE_EQUATION.md) — M2/M3/M4 vs. Lagrangian comparisons

---

## GLOSSARY — FIELD, VACUUM, SPACETIME, AETHER, GRANULES

These five terms get used interchangeably in physics writing but mean different things. This is the disambiguation for OpenWave's vocabulary.

### The ontological stack (bottom to top)

```text
╔══════════════════════════════════════════════════════════════╗
║  FIELD  ψ(x,t), n(x,t), φ(x,t)                               ║
║         ← MATHEMATICAL description / observable              ║
║         (what we compute with in code)                       ║
╠══════════════════════════════════════════════════════════════╣
║  VACUUM                                                      ║
║         ← STATE of the aether (field at its ground state)    ║
║         not a substance — a condition                        ║
╠══════════════════════════════════════════════════════════════╣
║  SPACETIME  ≈  AETHER                                        ║
║         ← the SUBSTANCE, viewed at large scales              ║
║         (OpenWave rejects empty spacetime — spacetime IS     ║
║          the aether)                                         ║
╠══════════════════════════════════════════════════════════════╣
║  GRANULES (Planck-scale)                                     ║
║         ← PHYSICAL microphysics — discrete particles         ║
║         Planck mass, oscillate, actually exist               ║
║         Bottom of the ontology                               ║
╚══════════════════════════════════════════════════════════════╝
```

### One-line definitions

| Term | One-liner | Type |
| --- | --- | --- |
| **Granule** | A discrete Planck-scale particle. The actual "stuff" of the universe. Has mass (Planck mass), oscillates, can move. | Physical object |
| **Aether** | The ensemble / sea of granules viewed as a substance. "The medium". | Substance |
| **Spacetime** | The aether viewed at large scales — where geometry lives. In OpenWave, spacetime IS the aether, not a backdrop for it. | Substance (= aether) |
| **Medium** | Synonym for aether. Used when emphasizing wave propagation. | Substance (= aether) |
| **Vacuum** | The *state* of the aether when no particles are present. Ground state of the field's potential V(ψ). Not "empty space" — a condition of the aether. | State |
| **Field** | The mathematical object we use to describe the aether's observable state. ψ(x,t), n(x,t), φ(x,t). A coarse-grained view of granule dynamics. | Mathematical construct |

### The four key conceptual moves

**1. Spacetime = aether.** Standard physics treats spacetime as an empty manifold; EWT / OpenWave treats it as the physical aether. Same thing at different vocabulary layers: physicists say "spacetime" for the geometry, we say "aether" when emphasizing the physical stuff, but they're the same substance.

**2. Vacuum ≠ aether.** "Vacuum" is *not another substance* — it's a *state*. Specifically, the state where the aether is at the minimum of its potential V(ψ), with no particles / defects / excitations. The aether (granule sea) still exists in vacuum; it's just undisturbed.

**3. Field ≠ aether.** The field ψ(x,t) is a *mathematical description* of the aether's state, not the aether itself. Different Lagrangians use different fields (scalar displacement, director, spin density, angle) — all observables of the same underlying granule ensemble. The aether is the *territory*; the field is the *map*.

**4. Granules are the microphysics; aether is the macro.** Granules are discrete; the aether (viewed at length scales >> Planck) appears continuous. Same relationship as molecules → fluid, atoms → crystal. The coarse-graining gives you a field; the underlying reality is granular.

### Where each term lives in OpenWave

| Term | Where you see it in code / docs |
| --- | --- |
| **Granule** | M1 engine (`m1_granule_motion/`), `constants.py` (Planck mass, granule density), scientific source papers |
| **Aether / Medium** | EWT papers (aether density ρ = 3.86×10²² kg/m³), README.md ("aether-like medium") |
| **Spacetime** | README.md ("spacetime subatomic simulator"), docstrings ("@spacetime module"). Often used synonymously with aether |
| **Vacuum** | Phase 3 / Lagrangian framework (`3_LAGRANGIAN_FRAMEWORK.md`, this doc) — the Duda / Close language |
| **Field** | M2 engine (scalar ψ), M4 engine (vector ψ), sandbox experiments (director n, kink φ) |

### Common confusions to watch out for

- **"Vacuum" does NOT mean "nothing there".** The granules are still present; they're just in their ground-state configuration. Vacuum has density, stiffness, supports waves. Classical physics calls this "full aether"; modern QFT calls it "vacuum with nonzero expectation values".
- **"Spacetime is curved" (GR) ↔ "the aether has a deformation" (OpenWave).** Same physics, different vocabulary. Einstein's "fields" are curvature of spacetime; EWT's fields are configurations of the aether.
- **"Field ψ" vs "aether"** — ψ is not the aether. ψ is our *mathematical bookkeeping* about what the aether is doing. Change Lagrangian → ψ's meaning changes. The aether doesn't change.
- **Granules vs. particles (electron, etc.)** — granules are the *microscopic* building blocks (Planck scale, ~10⁻³⁵ m). Particles like electrons (~10⁻¹⁵ m) are *composite topological structures* built from many granules. Don't confuse "granule" with "electron" — they're 20 orders of magnitude apart.

### Short cheat sheet

```text
PHYSICAL SUBSTANCE:  granules → aether → spacetime (same thing, different scales)
STATE OF SUBSTANCE:  vacuum (ground state) vs excited (contains particles/waves)
MATHEMATICAL MAP:    field ψ, n, φ — varies by Lagrangian
PARTICLES:           topological defects / excitations of the aether
WAVES:               perturbations of the aether around the vacuum state
```

---

## THE INTUITIVE LAYERING — USER'S FIRST PASS

Rodrigo's initial mental model, after working through the sandbox experiments:

1. **Medium / field** — spacetime / vacuum is an elastic medium at rest, static, but carries mass/density so it can transfer energy over distance
2. **WC** — wave centers are defects in the medium, a topological property of the field (not wave reflectors as in EWT, not external objects imposed on the field)
3. **Energy / waves** — waves are generated by the interaction of the medium with WCs, a disturbance in the field, kinks / twists … but where does the wave energy come from? If there is no base wave flowing through spacetime anymore, and waves are consequences of the defects, how exactly? And where do standing waves come in?

Layers 1 and 2 are correct. Layer 3 is the one the Lagrangian framework reframes.

---

## CORRECTED CONCEPTUAL LAYERING

```text
Layer 0 — Vacuum (static)                                [✅ Exps 2, 3, 4, 7]
    The elastic medium in its ground state. Not moving. Not flowing.
    What sits everywhere when "nothing is happening".
    In Duda's language: the MINIMUM of V(ψ), where all directors
    align peacefully.
    Validated: all experiments that seeded a vacuum (n = ẑ or Q = 0)
    confirmed it as a stable background.

Layer 1 — Topological defects (static structure)         [✅ Exps 1, 2, 3]
    Places where the vacuum's orientation CANNOT be smoothly returned
    to the ground state. These are the WCs. Charge, spin, and identity
    come from here. Still no motion.
    Validated: Exp 2 gave clean 1/d Coulomb between hedgehog + anti-
    hedgehog (R²=0.993, no sinc); Exp 3 confirmed winding number
    Q = ±1 is an integer, surface-independent, robust to 50% noise;
    Exp 1 showed 1D kink analogs carry stored mass energy (8·mc²).

Layer 2 — Perturbations around the vacuum (small, linear) [✅ Exp 4, partial Exp 7]
    Disturb the vacuum a little and it oscillates around the minimum.
    These oscillations propagate — that's a wave. Small perturbations
    satisfy Klein-Gordon:  (∂² − c²∇² + m²)ψ = 0.
    This is where waves come from in the Lagrangian picture.
    Validated: Exp 4 confirmed ω² = c²k² + m² to R² = 0.999982 across
    9 modes; Exp 7 (v2) confirmed Close's Eq. 19 gives transverse
    elastic-solid wave propagation — both are consistent with Layer 2.

Layer 3 — Non-perturbative dynamics (defect self-motion)  [⚠️ hinted, not directly tested]
    A defect is a LARGE, topologically protected perturbation. It has
    mass (stored energy, E = 8·m·c² for Sine-Gordon kinks — validated
    numerically in Experiment 1). Via the curvature of V(ψ) at the
    defect, it "trembles" at ω = 2mc²/ℏ on its own.
    This is the time-crystal mechanism — the defect oscillates without
    anyone pushing it.
    Validated: Exp 1 measured exactly E = 8·mc² for the static kink
    (0.06% off), and the mass-gap mechanism is confirmed at linear
    order by Exp 4. The full time-crystal Zitterbewegung oscillation
    requires a driven defect simulation, not done in sandbox.

Layer 4 — Radiation, standing waves, bound states         [⚠️ partial]
    Oscillating defects shake the surrounding vacuum → emit waves
    outward. When two defects are near each other, their emitted waves
    interfere → standing-wave patterns in the space between them.
    This is orbit quantization (hydrogen, Couder droplets).
    Validated: Exp 2 showed two defects interact at a distance via the
    intervening field (the 1/d Coulomb is a manifestation of this).
    Standing-wave lock-in between multiple defects wasn't directly
    tested — reserved for M5 multi-defect runs.
```

---

## WHERE WAVE ENERGY COMES FROM

Three sources, in order of subtlety:

### (a) From perturbing the vacuum

The vacuum is at the minimum of V(ψ). Push it slightly off the minimum and it oscillates — same reason a pendulum oscillates after you pull it sideways. The energy of the wave is the energy you put in when you created the disturbance.

**Physical analogy**: a drumhead at rest. Flat, no motion. Hit it → waves spread. The drumhead wasn't flowing before; *hitting it created the wave*. The waves carry the energy you deposited.

### (b) From defects themselves

A defect is a localized package of stored twist-energy (Sine-Gordon kink carries 8·m·c² of stored elastic energy — measured in Experiment 1). If the defect moves or accelerates, it shakes the surrounding vacuum, and that shaking propagates outward as waves.

**A static defect does not radiate** — same as a static electron doesn't emit light. Only when it accelerates, oscillates, or gets disturbed does it generate outgoing waves.

### (c) From the time-crystal mechanism (Duda's insight)

This is the deeper answer. A defect isn't *quite* static. Its internal structure interacts with the vacuum potential V(ψ) in a way that forces it to oscillate at its own rest frequency ω = 2·m·c²/ℏ (Zitterbewegung). Nobody pushes it; the vacuum curvature at the defect's location makes trembling intrinsic to having mass.

**This is what makes electrons oscillate at ~10²⁵ Hz in OpenWave's EWT picture — it's not an assumption, it's a consequence of the potential.**

---

## WHERE STANDING WAVES COME IN

Exactly where they come in for drumheads, musical instruments, and Couder walking droplets:

1. An oscillating defect emits waves outward
2. Those waves hit something — another defect, or the boundary of a bound state
3. Reflection creates inward waves
4. Inward + outward waves of the same frequency interfere
5. **Interference pattern = standing wave**

A hydrogen atom is an electron (defect) + proton (defect) whose mutually-emitted waves interfere. The standing-wave pattern that emerges has discrete modes — those are the allowed orbits. This is exactly what Couder's walking droplets reproduce: bouncing droplets on a vibrating fluid create standing-wave patterns that quantize their orbits.

**M3's standing-wave physics isn't discarded — it's reframed.** M3 said "particles ARE standing waves". Duda / Lagrangian says "particles are defects; standing waves form BETWEEN multiple defects, quantizing their mutual orbits". Same standing waves, different ontological role.

---

## WHAT EWT's "BASE WAVE" BECOMES

In the old framing (EWT, M1, M2):

```text
Base wave everywhere  ← assumed "always there" universal oscillation
WC = reflector / emitter
Particle = WC + its reflected waves
```

In the Lagrangian / topological framing:

```text
Vacuum everywhere     ← static ground state of V(ψ)
WC = topological defect (part of the field itself)
Wave = perturbation of vacuum (emerges from defect motion)
Standing wave = interference between multi-defect emissions
Time-crystal trembling  ← why defects oscillate → radiate at all
```

**EWT's "waves flowing everywhere" is replaced by "vacuum at rest; defects oscillate and make waves".** The universe isn't pre-filled with flowing energy — it's filled with a *potential* for energy, sitting quietly at the V(ψ) minimum, waiting for defects (matter) to stir it up.

---

## WHY THIS IS MORE FUNDAMENTAL, NOT LESS

The EWT framing had to *assume* the base wave exists ("all matter in the universe creates it") — fine as a model, but **unexplained**. The Lagrangian framing *derives* waves from two simple ingredients:

1. A potential V(ψ) with a minimum (vacuum)
2. Defects (particles) that disturb the vacuum

Everything else — wave frequencies, masses, forces, oscillations, charge quantization — follows as consequences of these two. No "why is there a base wave?" question; only "what is the potential?" (which is what Duda's challenge amounts to).

---

## THE CORRECTED MENTAL MODEL (side-by-side)

| Layer | Rodrigo's first pass | Corrected version |
| --- | --- | --- |
| Medium | Elastic medium at rest, has mass/density | ✅ right — *static*, not pre-oscillating |
| WC | Defects in the medium, topological | ✅ right — *intrinsic to the field*, not external |
| Waves | Disturbance in field, kinks/twists | Waves are **perturbations around vacuum** (small oscillations); **kinks/twists are the defects themselves**, not waves |
| Wave energy | ??? | Comes from: (a) external disturbance, (b) moving/accelerating defects, (c) time-crystal trembling of the defect itself |
| Standing waves | ??? | Interference between defect-emitted outgoing waves and their reflections / other defects' emissions — exactly like a drumhead |

---

## IMPLICATIONS FOR OPENWAVE'S DOCUMENTATION

This reframing prompts a revisit of the **Energy Layers** hierarchy (see [2_ENERGY_LAYERS.md](2_ENERGY_LAYERS.md), [0_OVERVIEW.md](0_OVERVIEW.md), [README.md](../README.md)). The current hierarchy starts at "Layer 1: Fundamental energy wave" and assumes a pre-existing base wave. A more complete picture inserts a **Layer 0 — Vacuum at rest** (static ground state of V(ψ)) before the wave physics begins.

The proposed updated hierarchy:

| Layer | Old (EWT / M1-M2-M3) | New (Lagrangian-aware) |
| --- | --- | --- |
| **0 — Vacuum** | — (implicit) | Static ordered vacuum (LdG minimum) |
| **1 — Fundamental wave** | Base wave from all matter | Perturbations around vacuum (Klein-Gordon) |
| **2 — WC disturbance** | Wave reflector creates standing wave | Topological defect creates surrounding field distortion |
| **3 — Particle formation** | Standing-wave standalone | Defect + interference from emissions / other defects |
| **4-7** | Unchanged | Unchanged — force emergence, gravity, composite particles |

Follow-up TODO tracked in [0_ROADMAP.md](0_ROADMAP.md) PHASE 3 — rewrite the layers in `2_ENERGY_LAYERS.md`, `0_OVERVIEW.md`, `README.md` with a Layer 0 prepended and the Layer 1–3 semantics refined.

---

## WHAT IS ψ, REALLY? (DISPLACEMENT OF WHAT?)

A natural follow-up question: **ψ is still displacement from rest, right? Displacement of what — field, vacuum, spacetime, aether, Planck-scale granules?**

This matters because the Lagrangian framework uses ψ (or φ, or n) as a central object, and its physical meaning is subtly different from what it was in M1-M4.

### Short answer

**In OpenWave's M1-M4: yes, ψ = granule displacement from equilibrium.** That intuition is correct there.

**In the Lagrangian framework: not always.** ψ is whichever field quantity the chosen Lagrangian operates on — often an *orientation* (unit vector, angle) rather than a spatial displacement. Every field theory is a *lens onto the granule dynamics*, and the Lagrangian determines which aspect of granule motion we track.

### What each method's ψ actually is

| Framework | Field symbol | What it physically represents | Units |
| --- | --- | --- | --- |
| **M1 — Granule Motion** | r(t) per granule | Position vector of each granule (actual particle) | length (am) |
| **M2 — Laplace-Propagation** | ψ(x,t) | Scalar displacement of the medium (coarse-grained) | length |
| **M3 — Wolff-LaFreniere** | ψ(x,t) | Scalar wave amplitude (analytical) | length |
| **M4 — Vector-Wave** | ψ(x,t) = (ψx, ψy, ψz) | Vector displacement (3 components per voxel) | length |
| **Sine-Gordon** (Exp 1) | φ(x,t) | **Angle / orientation** of a local rotational degree of freedom | radians |
| **LdG / hedgehog** (Exp 2, 3) | n(x,t) | **Unit vector / orientation** (director — which way a local axis points) | dimensionless |
| **Close's spin density** (Exp 7) | s(x,t) | **Local axis of rotational motion** of the elastic medium | angular momentum / volume |
| **Smolinski Ψ³** (Exp 8) | Ψ(x,t) | **Scalar amplitude** (closest to classical displacement) | length or derived |
| **Klein-Gordon** (Exp 4) | ψ(x,t) | Generic scalar field — interpretation depends on context | varies |

### What ψ actually "displaces" in each case

- **In M1/M2/M4**: granules physically move; ψ measures how far they've moved from their rest positions. Pure translational displacement.
- **In Sine-Gordon**: the field *rotates* locally. φ=0 means "at rest in the vacuum orientation"; φ=π means "rotated 180°". Nothing physically translates — the "displacement" is an angular shift of some local degree of freedom.
- **In LdG (our Exp 2/3)**: the field *aims* in a direction. n = ẑ means "pointing up"; a hedgehog means "pointing radially outward from a center". Nothing translates — the "displacement" is from the ground-state orientation.
- **In Close / Duda combined**: a mixture — some displacement-like + rotation-like degrees of freedom together.

### Connection to the granule picture: we never abandoned granules

This is the key reconciliation. The Lagrangian field theories **don't replace granules — they are different coarse-grained observables of the same underlying granule dynamics**:

```text
GRANULE LEVEL (always true, underneath everything):
  Each granule at position r(t) traces an elliptical trajectory
  Ellipse has:
    - shape (amplitude / eccentricity)
    - orientation (major axis direction, orbital-plane normal)
    - phase (where on orbit at t=0)
    - handedness (CW / CCW)

FIELD LEVEL (depends on which Lagrangian we pick):
  M2/M4 track:   granule DISPLACEMENT (position vector)
  LdG tracks:    granule ELLIPSE AXIS  (major axis direction — orientation only)
  Close tracks:  granule ROTATION AXIS (normal to ellipse plane)
  Sine-Gordon:   a scalar ANGLE         (aspect of rotation — setup-dependent)
```

This is the same "M4 is a superset" insight from [0_WAVE_EQUATION.md § M4's Elliptical Motion and the Director Field](0_WAVE_EQUATION.md#m4s-elliptical-granule-motion-and-the-director-field): the director field Duda operates on is a *projection* of M4's 6-phasor ellipse data. If M4 captures the full elliptical trajectory, we can extract any Lagrangian-picture observable (director, spin density, scalar amplitude) by reading off the appropriate aspect.

### What "at rest" means at each level

Subtle but important:

```text
GRANULE LEVEL:    granules NEVER truly rest — they oscillate at ~10²⁵ Hz
                  (Zitterbewegung / time crystal; ω = 2mc²/ℏ from the mass)

FIELD LEVEL:      "at rest" = all granule ellipses aligned the same way
                  (director field uniform, all pointing ẑ)
                  → the MACROSCOPIC orientation is at rest,
                  → even though microscopically every granule is still bouncing
```

**Analogy**: consider a gas in a room. Individual molecules are zipping around at ~500 m/s. But if their bulk motion averages to zero, we say "the air is at rest." The molecules aren't at rest — their *collective observable* (bulk velocity) is. Same here: granules always oscillate, but their *collective observable* (the director orientation pattern, the stress field, whatever ψ tracks) can be uniform = vacuum.

**"Vacuum at rest" is a statement at the ψ level, not the granule level.** This is a key reframing: it lets us say "the field is at rest" while the underlying Planck-scale microphysics is still perpetually bouncing. Both are true simultaneously, at different scales.

### Why this matters for M5 design

When M5 is built ([3c_path_to_m5.md](3c_path_to_m5.md)), we must pick: *what does ψ physically represent in our production engine?*

- **Option A (closest to M4)**: ψ = granule displacement (3D vector). The director is derived as the major axis of each voxel's elliptical trajectory (computed from phase/amplitude trackers).
- **Option B (pure LdG / Duda)**: ψ = director n(x) directly (3D unit vector). Model only the orientation, not the underlying ellipse.
- **Option C (Close / elastic solid)**: ψ = spin density (local rotation axis). Evolves under Close's nonlinear vector wave equation.

Experiments 2 and 3 used **Option B** (director only) — which is why they worked cleanly for topology. Exps 4–8 have now run: Exp 4 confirmed Klein-Gordon dispersion on a scalar ψ (Option A/B-neutral); Exp 7 v2 confirmed Close's vector-Q approach works as transverse elastic-solid dynamics (Option C). **The resolved M5 recipe** (per [3b § Winning Approach for M5](3b_lagrangian_experiments.md#winning-approach-for-m5) and [3c](3c_path_to_m5.md)) uses **Option B (director) for topology** (from Exps 2, 3) combined with **Option C (Close's Eq. 19) for wave dynamics** on top — both layers validated in the sandbox.

### The fundamental reframe

The single most important conceptual shift the Lagrangian perspective introduces:

> **The field quantity ψ doesn't have to be a displacement.** It just has to be (1) *something that has a well-defined minimum-energy configuration* (the vacuum) and (2) *something that can wind around a defect* (topology).

Displacement is one valid choice. Orientation is another. Angle is another. What unifies them all is that they're observables of the granule ensemble, and the one we pick becomes the ψ of our Lagrangian.

---

## EMPIRICAL VALIDATION — WHAT THE SANDBOX EXPERIMENTS CONFIRMED

With all 8 Phase 3 sandbox experiments complete (see [3b_lagrangian_experiments.md](3b_lagrangian_experiments.md) for full write-ups), we can now mark each conceptual claim in this document with its experimental status.

### Claim-by-claim scorecard

| Conceptual claim | Status | Validating experiment(s) |
| --- | --- | --- |
| Vacuum is a stable ground state (Layer 0) | ✅ Confirmed | All experiments — `n = ẑ` / `Q = 0` remains stable |
| Topological defects are real, with integer charge (Layer 1) | ✅ Strongly confirmed | Exp 3 (Q = ±1, robust to 50% noise) |
| Defects produce far-field Coulomb 1/d interaction | ✅ Strongly confirmed | Exp 2 (R² = 0.993, no sinc) |
| Defects carry stored mass-energy | ✅ Confirmed | Exp 1 (kink energy = 8·mc², 0.06% off) |
| Defects exhibit relativistic kinematics (Lorentz contraction) | ✅ Confirmed | Exp 1 (measured v = 0.4997c at input 0.5c; width = L/γ) |
| Perturbations of vacuum give massive Klein-Gordon waves (Layer 2) | ✅ Confirmed | Exp 4 (ω² = c²k² + m² to R² = 0.999982) |
| Mass emerges from the curvature of V(ψ) | ✅ Confirmed at linear order | Exp 4 (quadratic potential → mass gap observed) |
| Close's vector wave equation = transverse elastic-solid wave | ✅ Confirmed | Exp 7 v2 (Eq. 19 implementation gives dispersing transverse waves) |
| Smolinski's Ψ³ is derived from a quartic Lagrangian | ✅ Confirmed (sympy) | Exp 5 Test 2 |
| Noether gives energy conservation `H = T + V` | ✅ Confirmed (sympy) | Exp 5 Test 3 |
| Combined W-L product form (docs) is a free-wave solution | ❌ Falsified | Exp 5 Test 1d — residual `−A·c²k²·sin(ωt+φ)/r ≠ 0` |
| Nonlinearity alone can produce K-selectivity | ❌ Falsified | Exp 8 (K=1, 2, 4, 6, 8 all behave the same under Ψ³) |
| Biaxial geometry can produce three distinct energy scales | ⚠️ Mechanism confirmed, specific ratios deferred | Exp 6 (E ∝ K linear; lepton ratios reproduced by construction, not derived) |
| Close's equation supports soliton emergence from harmonic seeds | ⚠️ Not confirmed — but Close's framework doesn't actually require it | Exp 7 v2 (particles in Close's picture = plane-wave bispinors, not solitons) |
| Time-crystal Zitterbewegung emerges naturally | ⚠️ Hinted, not directly tested | Mass-gap mechanism (Layer 3) validated by Exp 4 at linear order |
| Standing-wave lock-in between multiple defects (Layer 4) | 🚧 Reserved for M5 | Not in sandbox scope |

### What the pattern of results tells us

1. **Topology (Layer 1) is the most strongly validated ingredient.** Exps 1, 2, 3 each decisively confirmed different aspects (stability, Coulomb, quantization). Topology is the *load-bearing structural element* of the whole framework
2. **Klein-Gordon wave dynamics (Layer 2) are also cleanly validated** — but this is mostly textbook physics we reproduced numerically to make sure our leapfrog PDE solver is correct. The interesting step is that the *mass gap from a potential* (not assumed, not imposed) emerges automatically
3. **Nonlinearity alone is NOT sufficient** — Exp 8 killed the "Smolinski Ψ³ gives K-selectivity" hypothesis cleanly, and Exp 7 v2 showed that a nonlinear vector equation without topology also doesn't spontaneously localize from harmonic seeds. This constrains M5: it must have topology, not just nonlinearity
4. **Documentation integrity check**: Exp 5 caught a concrete error — the product form of Combined W-L in our docs does not actually satisfy the free-wave equation. The code's sum form is correct; the docs needed updating

### Updated conceptual layering (after the sandbox)

The corrected layering in [CORRECTED CONCEPTUAL LAYERING](#corrected-conceptual-layering) above now has inline `[✅ Exp N]` markers showing which experiments validated each layer. Layers 0, 1, 2 are empirically solid; Layer 3 (time-crystal) is hinted at by Layer 2 + Exp 1; Layer 4 (standing waves between defects) is reserved for M5 multi-defect runs on the production engine.

### Refinements to the mental model

The one conceptual shift that came from the experiments (vs. just validating what we expected):

- **"Combined W-L is an exact free-wave solution"** — what the M3/M4 docs claimed — was wrong. The product form `2A·sin(kr/2)·cos(kr/2−(ωt+φ))/r` contains a quadrature piece `(1−cos(kr))/r·sin(ωt+φ)` that is NOT a free-wave solution. The M4 code's equivalent sum form IS a solution; they're not algebraically identical. For M5, we keep the sum form and discard the claim that the product form is Lagrangian-derived. See [0_WAVE_EQUATION.md § M3: The Quest for the Wave Equation](0_WAVE_EQUATION.md#m3-the-quest-for-the-wave-equation-vs-lagrangian-mechanics) and [3b § Experiment 5](3b_lagrangian_experiments.md#experiment-5-lagrangian-derivation-combined-w-l-from-a-lagrangian) for the full reasoning

This is exactly the kind of finding the Lagrangian framework was designed to surface: when you formalize the equation's origin, latent inconsistencies in empirical modeling choices show up.

---

## KEY TAKEAWAYS

1. **Vacuum first, waves second.** The elastic medium at rest is the starting point; waves are *oscillations of* the vacuum, not something that *flows through* it. **[✅ Exp 4: Klein-Gordon dispersion emerges from perturbing a quadratic-potential vacuum]**
2. **Defects are features of the field, not external objects.** A WC is a topologically protected configuration of the vacuum, not a particle-with-stuff-attached. **[✅ Exps 1, 2, 3: kinks, hedgehogs, winding numbers all validated as field configurations]**
3. **Waves emerge from disturbing the vacuum** — via external perturbation, defect motion, or intrinsic time-crystal trembling. **[✅ Layer 2 confirmed in Exp 4, Layer 3 hinted at by Exp 1 mass measurement]**
4. **Standing waves are interference patterns between defect emissions** — the same mechanism as musical instruments or Couder droplets. M3's lock-in physics is preserved in this role. **[🚧 Multi-defect lock-in reserved for M5 — not directly tested in sandbox]**
5. **Energy sources are explainable**, not assumed. The time-crystal mechanism tells us *why* electrons oscillate at ω = 2mc²/ℏ instead of assuming a universal base wave. **[✅ Mass-gap mechanism validated at linear order by Exp 4 (`ω² = c²k² + m²`)]**
6. **ψ is not always a displacement.** In M1-M4 it is; in the Lagrangian framework it can be a displacement, an orientation, an angle, or a rotation axis. The granules are always underneath; ψ is a coarse-grained observable chosen by the Lagrangian. "Vacuum at rest" is a statement at the ψ level, not the granule level — granules never truly rest. **[Framework-level claim — confirmed by which observables each experiment used: Exp 1 used angle φ, Exps 2/3 used director n̂, Exp 4 used scalar ψ, Exp 7 used vector Q]**
7. **Topology is the load-bearing ingredient — nonlinearity alone is insufficient.** This was the emergent lesson from the full Phase 3 program. **[Pattern evident from Exps 2, 3 (✅ topology) + Exp 8 (❌ pure nonlinearity) + Exp 7 v2 (⚠️ nonlinear vector equation without topology)]**
