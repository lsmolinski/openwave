"""
XPERIMENT PARAMETERS — Boost-DRESSED Hedgehog + 4D Minkowski evolve (M5.8.2c)

THE TIME AXIS GOES LIVE. Seeds the M5.8.2b-1 ground state — the biaxial hedgehog
dressed with a core-localized boost exp(b*·w(r)·B₁) mixing the e_Θ eigen-axis with
the time axis (b* = 0.13, the GEM dip: the dressed defect sits BELOW the bare one)
— and evolves it under the 4D SIGNED curvature flux (F → ηFη, the 5a §10d sign
rule), replacing the M5.8.1 time-freeze clamp with the soft coherent-drift guard.

GHOST GUARD (the 2b-2/2c-1 lesson chain): the signed kernel is ill-posed on the
"fuel shell" where the gradient-stiffness form Q(∇ψ) loses positive-definiteness
(the linear-propulsion runaway, λ≈15.6/t in the sandbox). The launcher computes a
per-voxel stable_mask ONCE post-seed (analytic Cardano, not ti.sym_eig): the
stable region (~77% at b*) evolves Minkowski-signed; the fuel shell falls back to
the always-stable Euclidean flux. Its genuine propulsion dynamics needs the
constrained faithful kernel (sandbox m5_8_2c1) — future production iteration.

CLOCK_KICK starts the de Broglie clock (the 2c-1 finding: the kicked clock HOLDS
and concentrates into the core; spontaneous start is under diagnosis). Set 0.0
for the pure static-survival run (the G-2c-1 dispersal test).

V(M): ON with the standard b=0 amplitude well (2c-0: the GEM dip SURVIVES the
whole κ scan; the clock rides a known washboard ~20% of the dip per unit κ).
Set LDG_STIFFNESS_K: 0.0 for the V-off baseline.

DT_SCALE_4D derates dt (the fastest 4D collective mode runs above the clock rate
— the 2b-1 timescale-hierarchy lesson). Boots PAUSED on the Hamiltonian view:
press Evolve PDE and watch whether the dressed defect HOLDS (G-2c-1).
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # ~262k voxels — small for fast smoke-test

TOPOLOGY_SEED = {
    "MODE": "dressed_hedgehog",
    "CENTER": [0.50, 0.50, 0.50],
    "R0_FRACTION": 0.06,
    "RHOC_VOXELS": 3.0,
    "BIAXIAL_DELTA": 0.30,
    "B_STAR": 0.13,  # the 2b-1 GEM-dip dressing amplitude
    "RW_FRACTION": 0.29,  # dressing width / box edge (~r_w=3.5 of L=12 sandbox)
    "CLOCK_KICK": 0.0,  # rad; >0 starts the clock (e.g. 0.05)
    "DT_SCALE_4D": 0.5,  # dt derate for the 4D fast modes
    "KM_INERTIA_4D": 30.0,  # diagonal faithful-lite inertia m=1+km·dx²Σ‖∂M‖² (2c-2)
    # SIGNED_FLUX_4D: the Minkowski-signed force. OFF (default) = the SAFE v1:
    # Euclidean flux with the time axis live (dressing renders + holds). ON is
    # UNSTABLE pending the constrained spectral-projection kernel (2c-1) port —
    # the cheap-inertia signed system grows on ~700-step timescales (2c-2).
    "SIGNED_FLUX_4D": False,
    "AUTO_RELAX_STEPS": 0,
    "LDG_STIFFNESS_K": 1.0,  # V ON — the 2c-0-validated amplitude well
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "Hedgehog-Dressed 4D",
        "DESCRIPTION": "Boost-dressed hedgehog + LIVE time axis — 4D signed flux, stable-masked",
    },
    "camera": {
        "INITIAL_POSITION": [1.10, 1.46, 0.81],
    },
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],
        "TARGET_VOXELS": TARGET_VOXELS,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "VIZ_STRIDE": 1,
        "SHOW_GLYPHS": 2,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 2,
        "WARP_MESH": 0,
        "SHOW_GRANULES": False,
        "SIM_SPEED": 1.0,
        "PAUSED": True,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 4,  # Hamiltonian energy — watch the dressed defect HOLD (G-2c-1)
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
