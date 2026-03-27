# LEVEL 3: PARTICLE EMERGENCE

Near-field standing wave physics — demonstrating that particles emerge from wave interference.

Task checklists in [00_ROADMAP.md](00_ROADMAP.md). Building blocks hierarchy in [00_OVERVIEW.md](00_OVERVIEW.md).

---

## WAVE EQUATION: Weighted Partial Standing Wave

Superposition of two counter-propagating spherical waves with a spatial blending function:

```text
ψ(r,t) = A · [w(r)·sin(kr + ωt + φ) + sin(kr - ωt - φ)] / kr
              \______ in-wave ______/   \_____ out-wave _____/
```

### What IS grounded in wave physics

- **`sin(kr ± ωt) / kr`** — exact analytical solution to the 3D wave equation in spherical coordinates. Textbook.
- **Superposition principle** — adding two solutions gives another solution. Textbook.
- **In-wave + out-wave** — Wolff's Space Resonance model. The in-wave is the Huygens reconstruction from all other matter in the universe (every out-wave from every other particle converges back). Established wave physics (Huygens-Fresnel principle).
- **Standing wave at core** — when w=1, the superposition gives `2·sin(kr)·cos(ωt+φ)/kr`, a pure standing wave. The particle IS the standing wave. Central claim of Wolff/EWT.
- **Traveling wave at distance** — when w=0, only the out-wave remains: `sin(kr-ωt-φ)/kr`. The force field IS the traveling wave.
- **Sinc lock-in** — energy wells at `kr = nπ` from `sin(kr)/kr` nodes. Direct consequence of the wave equation, not imposed.
- **Annihilation** — opposite-phase superposition cancels. Direct consequence of superposition.

### What's a modeling choice (not derived from first principles)

1. **The weight function** `w(r) = 1 / (1 + (r / 1.25λ)^8)` — the shape, the power (8), the transition distance (1.25λ). Physically motivated (in-wave convergence weakens with distance) but the specific function is chosen, not predicted.
1. **Monochromatic** (single frequency) — real particles may be broadband. This is the root cause of the sinc barrier issue at far-field.
1. **The 1.25λ transition distance** — chosen to include the first quadrature lobe in the standing wave zone. Not derived from EWT constants.

### Amplitude profile: sinc steepness and energy conservation

The standing wave spatial function `sin(kr)/kr` (sinc) creates a steep amplitude rise toward the center. This is correct physics for 3D spherical geometry — not a numerical artifact.

**Energy in a thin shell at radius r:**

```text
E_shell = ρ · (f·A(r))² · 4πr² · dr

With A(r) = 2A₀ · sin(kr)/kr:

E_shell ∝ r² × (sin(kr)/kr)²  =  sin²(kr) / k²
```

The `r²` from shell surface area **exactly cancels** the `1/r²` from the sinc envelope. Each half-wavelength shell contains roughly the same total energy. Amplitude concentrates toward the center because the same energy packs into a smaller volume (like a whirlpool: speed ∝ 1/r, cross-section ∝ r², total flow constant).

**Amplitude values across the standing wave:**

| Distance | kr | sinc | Envelope (×2A₀) | Shell volume ∝ r² |
| --- | --- | --- | --- | --- |
| center | 0 | 1.000 | 2.00 | tiny |
| λ/4 | π/2 | 0.637 | 1.27 | 16× larger |
| λ/2 | π | 0 | 0 (node) | 64× larger |
| 3λ/4 | 3π/2 | 0.212 | 0.42 | 144× larger |

**Why it doesn't diverge at r=0:** The sinc function is the spherical Bessel function j₀(kr), the regular solution that stays finite at the origin. Unlike bare `1/r`, the `sin(kr)/kr → 1` at center (L'Hôpital), giving a finite peak amplitude of 2A₀.

**Variable λ(r) — future refinement:** The current sinc assumes constant λ everywhere. The Yee & Hauger shells and the energy steepness concept (`A/λ ≈ const`) suggest λ grows near the center too: `A(r)↑, λ(r)↑`. This would make the profile less steep than pure sinc while conserving energy steepness. Deferred to Block 1b (variable λ for tetrahedron stability).

### Sensitivity to investigate

- Does changing the transition distance (0.5λ, 1.25λ, 2λ) qualitatively change lock-in physics, or just shift barrier distances?
- Does changing the power (2, 4, 8) affect stability?

---

## SIGNED ENVELOPE — Energy & Force

The force is computed from `F = -∇E` where `E = ρ·V·(f·A)²` and A is the signed envelope.

The envelope is the **single-WC phasor magnitude** from the same wave equation, multiplied by a charge sign:

```text
envelope = charge_sign × A_eff × |phasor|

|phasor| = √( (w+1)²·sinc² + (w-1)²·cosc² )

where:  sinc = sin(kr)/kr,  cosc = cos(kr)/kr
```

### Behavior by regime (automatic from the weight function)

| Region | Weight | \|phasor\| | Physics |
| --- | --- | --- | --- |
| Center (r=0) | 1.0 | 2.0 | Peak standing wave |
| Near-field (r < 1.25λ) | ~1.0 | `2·\|sin(kr)\|/kr` | Sinc oscillation → lock-in wells |
| Transition | 0→1 | blended | Smooth rolloff |
| Far-field (r >> 1.25λ) | ~0.0 | `1/kr` | Smooth 1/r (no sinc) |

### Why far-field is smooth

In the far-field (w=0), sin and cos are in **quadrature** (90° apart):

```text
|phasor| = √( sin²(kr)/kr² + cos²(kr)/kr² )
         = √( 1/kr² )       ← sin² + cos² = 1
         = 1/kr              ← smooth, no oscillation
```

The standing wave (w=1) has **no quadrature component** (`S_n = 0`), so the sin zeros survive as sinc oscillation. The traveling wave has both components → oscillation cancels.

### Charge sign — imposed, not emergent

`charge_sign = cos(source_offset)` gives ±1 based on phase (0° = positron, 180° = electron). This is NOT emergent from wave interference — it's an imposed label that gives correct Coulomb direction as a placeholder.

- Same sign → constructive superposition → high E → repulsion gradient
- Opposite sign → destructive superposition → low E → attraction gradient

True emergent charge from spin-based L→T conversion is deferred to **Block 2 / M4 vector field**.

### The sinc barrier problem (Phase 1 conclusion)

Coherent monochromatic wave interference ALWAYS produces `cos(k·Δr)` oscillation in the phasor cross-term. No single-frequency wave equation avoids this (all 10 models tested in Phase 1). The signed envelope sidesteps this for far-field by using the single-WC phasor magnitude (smooth 1/r) instead of the multi-WC phasor (oscillating cross-terms). Near-field sinc structure is preserved for lock-in physics.

Carry-over approaches for emergent far-field Coulomb (Block 2): 3D flux, variable λ(r), non-linear Ψ³, K=10 scale averaging.

---

## BLOCK 1 CLAIMS — What we're demonstrating

All claims emerge from wave interference (not imposed):

1. **Lock-in** — same-phase WCs lock into sinc energy wells at λ/2 intervals
1. **Annihilation** — opposite-phase WCs cancel via superposition (deepest well at r=0)
1. **K=10 stability** — 1-3-6 tetrahedral geometry is the first configuration where all WCs sit near standing wave nodes
1. **K=2..9 instability** — intermediate geometries can't accommodate all WCs at nodes → decay

What's acknowledged as placeholder:

- Far-field force direction from imposed charge sign (Block 2)
- Weight function parameters (sensitivity testing needed)

---

## M3 SIMULATION PIPELINE

Scripts at `openwave/xperiments/m3_wolff_lafreniere/`. Launched from `_launcher.py`.

### Per-frame computation

1. **`wave_engine.propagate_wave()`** — computes displacement (weighted partial standing wave), phasor RMS, signed envelope, EMA-RMS, energy field
1. **`xforce_motion.compute_force_vector()`** — `F = -∇E` with weighted multi-shell gradient (radius=3, 1/d² falloff). Energy from signed envelope
1. **`xforce_motion.integrate_motion_leapfrog()`** — symplectic integrator, 0.995 damping
1. **`xforce_motion.detect_annihilation()`** — opposite-phase pairs within 5 grid units

### Xperiments

| Xperiment | Setup | Goal |
| --- | --- | --- |
| annihilation1 | 2 WC, opposite phase, head-on | Demonstrate wave cancellation |
| annihilation2 | 2 WC, opposite phase, diagonal | Demonstrate off-axis annihilation |
| lock_in | 2 WC, same phase, close | Demonstrate standing wave lock-in |
| tetrahedron | 10 WC, 1-3-6 geometry | Stabilize K=10 electron |

### Known issues

- **K=10 tetrahedron unstable** — 15/45 WC pairs at non-node distances (√3×λ/2, √2×λ/2). Needs variable λ(r) and/or wave equation tuning.
- **Damping (0.995)** — drains KE over ~500 frames. May need tuning per xperiment.
- **Annihilation threshold (5 grid units)** — may need calibration for different grid resolutions.
