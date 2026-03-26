# ✅ Phase 1a: Signed Disturbance (forced charge sign)

Modeled WCs as signed amplitude modulators (not actual wave disturbances): `A_total = A₀ + Σ q·A_peak·δ(r)` with `q = cos(phase) = ±1`.

**Test results (equation #6 in wave_engine_1D_v2.py)**:

- **Same charge repulsion**: 9/9 correct direction, near-constant Coulomb ratio ✓
- **Opposite charge attraction**: 0/9 correct — asymmetric energy landscape, Newton's 3rd law violated (forces differ by ~40x)

**Root cause — the linear cross-term gives wrong charge dependence**:

```text
E = ρVf²(A₀ + q₁δ₁ + q₂δ₂)²

Expanding:
  A₀²              → constant, no gradient, no force
  2A₀·q₁δ₁         → force ∝ q₁ (individual charge, gravity-like)  ← DOMINANT
  2A₀·q₂δ₂         → force ∝ q₂ (individual charge, gravity-like)  ← DOMINANT
  q₁²δ₁²            → self-energy, zero gradient at WC center
  2q₁q₂·δ₁δ₂       → force ∝ q₁q₂ (charge product, Coulomb-like)  ← CORRECT but small
  q₂²δ₂²            → self-energy, zero gradient at WC center
```

The dominant force term `2A₀·q·∇δ` depends on the **individual** charge sign q, not the **product** q₁·q₂. This is equivalent to the previously ruled-out "smooth envelope with imposed charge sign" — the charge enters as a ±1 label, not through wave interference. It doesn't actually simulate a base wave — it's a signed potential, not a wave disturbance.

**The emergence test**: if you can replace `q = cos(phase)` with `q = +1` or `q = -1` as a manual input and get the same result, the charge sign is not emergent.

**Implementation**: equation #6 "Signed Disturbance" in `wave_engine_1D_v2.py` with `BASE_AMPLITUDE_RATIO` parameter. Kept for experimentation.
