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
