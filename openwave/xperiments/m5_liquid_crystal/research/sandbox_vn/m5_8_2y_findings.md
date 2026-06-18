# M5.8.2y — The 28× absolute-ω gap: lever budget (issue #217)

> Does turning V(M) on lift the clock ω toward the electron ZBW, closing the 28× gap (#208)? **Answer: no — the gap is largely structural.** Driver: [`m5_8_2y_gap_levers.py`](m5_8_2y_gap_levers.py); the faithful-kinetic input re-run via [`../sandbox_v6/m5_6_5d_faithful_kinetic.py`](../sandbox_v6/m5_6_5d_faithful_kinetic.py). Date 2026-06-17.
>
> Built on VALIDATED M5.6.5 diagnoses (5c rotation-invariance + 5d faithful-kinetic metric), not a fresh fragile V-on evolution. The V-on null is **analytic** (below), which is stronger than a single noisy run.

## The two named levers, settled

### Lever 1 — V(M) on (the LdG amplitude potential): **NULL** on frequency

The LdG potential `V(M) = a·Tr(M²) + c·(Tr(M²))²` (the `b=0` confining form, M5.6.5c) is a function of the **rotation-invariants** `Tr(M²)`, `Tr(M³)` only. The de Broglie clock is a **pure rotation** of the director frame at fixed eigenvalues (the orientation/twist slosh). Therefore:

```
∂V/∂(clock coordinate) = 0   exactly   (V is rotation-invariant; the clock is a rotation)
```

so V-on adds **zero DIRECT force** to the clock motion → zero direct frequency shift. This is analytic, and it matches the M5.6.5c empirical observation: *"V pins the amplitude, NOT the frame orientation, so the directors still slosh with V on — that is the dynamical twist sector."* The clock lives precisely in the sector V cannot touch.

(Honest caveat: the *direct* force is exactly zero; an *indirect* effect is possible via the amplitude profile setting the effective core size that the twist happens in. But 5c shows V-on pins `Tr(M²)` near `s₂*` with only a modest profile change, so the indirect ω shift is small — and even a generous indirect effect cannot supply a 28× lift. V-on is null-to-small on the clock frequency either way.)

⇒ **The #208 guess that the V-on core geometry lifts ω is refuted for the V lever.** V confines amplitude; it does not move the clock.

### Lever 2 — the faithful F-commutator kinetic: the **one real lever, bounded ×3.04**

The shipped kinetic `½‖Ṁ‖²` gives all modes the uniform inertia 0.5. The faithful Eq.18 metric `G = 4Σ_μ(−ad²_{M_μ})` weights them correctly. Re-run of M5.6.5d (2026-06-17, δ=0.3):

| Quantity | Value |
| --- | --- |
| faithful physical-mode inertia (5 modes, 5–95% span) | **[0.054, 1.447]** |
| simple kinetic inertia (uniform) | 0.5 |
| ⇒ clock frequency mis-set by `×√(0.5/inertia)` | **×[0.59, 3.04]** |
| best case (clock = softest mode 0.054) | ω rises **×3.04** |
| mid case (clock = median mode 0.265) | ω ×1.37 |

The sign is clock-mode-dependent, but the **bound is ×3.04**: the faithful kinetic can close **at most ~3×** of the 28×.

## The lever budget vs the 28× gap

| Lever | Factor on ω | Closes |
| --- | --- | --- |
| V-on (LdG amplitude) | **1.00 (null)** | nothing (rotation-invariant) |
| faithful kinetic | ×[0.59, **3.04**] | at most ~3× |
| **Residual (best case)** | — | **≥ 9.3×** still unaccounted |
| Residual (mid case) | — | ~20.5× |

## Verdict — the gap is largely STRUCTURAL

| Component | Size | Nature |
| --- | --- | --- |
| faithful-kinetic correction | ~×3 (max) | a genuine **dynamics** lever (the one physical piece) |
| residual | **~9–20×** | **structural** — the missing LENGTH ANCHOR at V=0 |

The dominant residual is **not** a dynamics error. At V=0 there is no independent length scale fixing the absolute frequency — the action↔ℏ postulate carries it, and that is exactly where the ~9–20× sits. Closing it needs the **Faber `r₀` length chain** (`E·r₀=const`, `r₀=2.2132 fm`), an **independent length calibration**, not the V-on or kinetic physics.

So #208's "structural" read is **confirmed and now quantified**:

> lever budget ≈ **3× physical** (faithful kinetic) + **9–20× structural** (length anchor). The sim-units-ratio fallback stands; absolute Hz needs the `r₀` length fix, not a dynamics change.

## DoD status

- [x] The two named levers' effect on ω settled: V-on **null** (analytic + 5c), faithful kinetic **×[0.59, 3.04]** (re-run 5d).
- [x] Residual quantified per piece: ≥9× (best) to ~20× (mid) unaccounted by the physics levers.
- [x] Verdict: **largely structural** — ~3× physical lever budget, ~9–20× = the missing length anchor (the Faber `r₀` chain, an independent length calibration).

## Honest scope / follow-ups

| Item | Note |
| --- | --- |
| V-on null is analytic | exact (V = f(rotation-invariants); clock = pure rotation). An empirical V-on clock evolution would only confirm it; not needed for the verdict |
| The faithful-kinetic factor is a **bound** | pinning the exact clock-mode factor within [0.59, 3.04] needs mapping `Ṁ_clock` onto G's eigenmodes — a refinement; it cannot exceed ×3.04, so it cannot change the verdict |
| The real closure path | the **Faber `r₀` independent length calibration** (the residual is a length-anchor deficit, not a dynamics one). That is the next absolute-scale step, not another V-on/kinetic run |

Model: M5 Liquid Crystal. Physics-only, headless. Follows #208 (the calibration framework).
