# PHASE 1d: Vector Wave Force (M4 displacement direction)

Divergence/curl/flux from M4 vector displacement; recovers charge sign from rotation direction

**Problem**: F = -∇(|ψ|²) uses scalar magnitude, which discards vector direction information. On-axis, vector reduces to scalar — no help for the standard test case.

**Opportunity**: vector displacement carries information beyond magnitude — ellipse rotation direction (handedness), divergence, curl, energy flux direction. These are **signed quantities** that could recover charge-phase information.

Force must be computed from a **different quantity** than |ψ|²:

- **Divergence** (∇·ψ): compression/rarefaction — scalar but signed
- **Curl** (∇×ψ): rotational displacement — vector, related to magnetic field
- **Energy flux** (ψ × ∂ψ/∂t or similar): directional energy flow
- **Per-component amplitude** (A_x, A_y, A_z separately): preserves directional structure

**One force, different directions**: F = -∇E is one force — electric (longitudinal), magnetic (transverse), gravitational (density deficit) are projections onto different components. Scalar collapses all directional information into magnitude — correct scaling (1/r²) but wrong direction (charge sign). See [04_magnetic_vector.md](04_magnetic_vector.md#why-scalar-is-insufficient-monopole--longitudinal-only) for full analysis.

**Connection to non-linear equations (Phase 1c)**: non-linear Ψ³ soliton, toroidal wave flows (r⁵), spin-as-vortex all require vector displacement. Phases 1c and 1d may converge.

Maybe this path needs 3D simulation, definitely 1D is not enough, but possibly 2D is not enough either. On the force unification concept, there is only ONE force, but at human scale (inertial frame / mass scale / frequencies that this scale of mass can experience), this single force appears at defined conditions that makes us perceive them as separate forces, so we named and describe them as so, but if they are a single 3D elliptical behavior that can be decomposed into 2 major amplitudes (90 degrees apart) and this elliptical form can be oriented in multiple orders in 3D space, this opens up the possibility. 2D simulation won't capture that.

## Other Possible Solutions

### Energy Flux (Radiation Pressure)

Force can also arise from **energy flux** — the directional flow of energy through the medium:

- **Energy density** (current: `F = -∇E`): energy *stored* per voxel. Force from the energy landscape shape
- **Energy flux** (`S = E · v_group`): energy *flowing* through a surface per unit time (W/m²). Has direction
- **Radiation pressure** (`P_rad = S/c`): force per unit area from wave momentum transfer (LaFreniere's mechanism)

For constant c: `∇P_rad = ∇E` — both approaches give the same force. They diverge when c varies spatially or waves are directional.

Energy flux could **naturally separate standing wave (near-field) from traveling wave (far-field) contributions** — standing waves have zero net flux, traveling waves have nonzero flux.

---

## Findings from Phase 1b Step 2c (hints for Phase 1d implementation)

Phase 1b tested L→T spin conversion (Option G) using the quadrature base wave's two phasor channels (P, Q) as a proxy for longitudinal and transverse components. The WC converts a fraction of L (P) into T (Q) with charge-dependent direction: `q = +1 → L→+T` (CW), `q = -1 → L→-T` (CCW).

**Result — the ONLY model that distinguishes charges**:

- Opposite charge: 12/24 attract, 12/24 repel (oscillates — sinc partially disrupted)
- Same charge: 24/24 unclear (both forces same direction — distinct from opposite)
- No other model (7 tested: 4 passive, 3 elastic) produced different behavior for opposite vs same charge

**Why the quadrature proxy is limited**: the two phasor channels P and Q are not truly independent — they combine into magnitude via `RMS = √(P² + Q²)/√2`. Converting energy between P and Q changes their individual values but the magnitude (and therefore energy) is still dominated by the sinc-modulated interference pattern. The L→T conversion shifts energy between channels, but `P² + Q²` doesn't fully decouple them.

**What Phase 1d needs**: true two-component displacement — independent L and T fields that contribute to energy **separately**:

```text
E_total = E_L + E_T = ρV(f·A_L)² + ρV(f·A_T)²
```

With independent components:

- L→T conversion at the WC reduces A_L and increases A_T locally
- The T component is NEW (not in the incoming field) → breaks isotropic cancellation (M2's core finding)
- The T component doesn't have the same spatial sinc structure as L → its contribution to the energy gradient is different
- Force direction depends on HOW the L/T balance changes with charge sign → charge-dependent force from wave physics, not imposed labels

**Implementation path**: extend the 1D engine with a second displacement field (ψ_T alongside ψ_L). Each has its own phasor (P_L, Q_L, P_T, Q_T). At WC positions, apply L→T conversion with charge-dependent direction. Energy computed from both: `E = E_L + E_T`. This is a minimal vector extension — 1D spatial domain but 2D displacement (longitudinal + transverse), sufficient to test the spin hypothesis without full 3D.

**Connection to M2 spin theory**: M2's `14_spin_theory.md` proposed exactly this mechanism and attempted implementation (`interact_wc_spinUP/DOWN`) but it "never worked correctly" in the 3D isotropic Laplacian field. The Phase 1b quadrature proxy test confirms the concept works (charge discrimination achieved) but the implementation needs true independent components. Phase 1d should revisit the M2 spin code with the improved understanding from Phase 1b testing.

**Convergence with Phase 1c**: Phases 1c and 1d may converge into one solution — variable λ(r) in the energy equation (1c) + vector displacement with L→T conversion (1d). The L→T spin creates charge-dependent asymmetry, while λ variation creates force from wavelength gradients. Both are needed: L→T alone still oscillates in scalar RMS, λ variation alone has no charge sensitivity. Together they could produce charge-dependent force from wavelength + mode gradients.
