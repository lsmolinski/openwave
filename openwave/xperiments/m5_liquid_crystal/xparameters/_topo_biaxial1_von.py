"""
XPERIMENT PARAMETERS — Biaxial Hedgehog + V(M) confinement (M5.6.5c)

Same biaxial hedgehog seed as _topo_biaxial1.py, but with the Eq.13 LdG potential
V(M) turned ON via LDG_STIFFNESS_K. This is the A/B partner to _topo_biaxial1 (V off):
press Evolve PDE on both and compare how the energy behaves.

WHAT V DOES HERE (M5.6.5c finding).  With V off, Evolve PDE makes the field slosh and
its energy DILUTES over a growing radius — there is no restoring force against amplitude
spread (bounded, energy-conserving, but not localized). V supplies that restoring force.

We use the b=0 amplitude well

    V(M) = a·Tr(M²) + c·(Tr(M²))²  ,  minimum at Tr(M²) = s₂* = 1 + δ²

which pins the field amplitude (confines the energy) while leaving the biaxiality δ
EXACTLY flat. This matters: the canonical 3-term Eq.13 LdG (a·Tr(M²) − b·Tr(M³) +
c·(Tr(M²))²) has NO biaxial minimum — for any (a,b,c) the nonzero eigenvalues collapse
to a single λ* (5a §5f) — so a b≠0 term would confine but erode δ toward uniaxial. b=0
confines without uniaxializing. (A fully biaxial-STABLE vacuum needs an extra invariant
in V — Duda's open Q7.)

SCALING.  The launcher computes the production coefficients from LDG_STIFFNESS_K:
    ldg_c = K · c_amrs²/dx_am⁴ ,  ldg_a = −2·ldg_c·(1+δ²) ,  ldg_b = 0
The c²/dx⁴ unit is the cubic-balance scale matching the F²-curvature force. A production
sweep (m5_6_5c_prod_scale) confirmed K ∈ [0.5, 25] all confine ~3.3× vs V-off with no
blow-up; K=1.0 is the clean default.

Boots PAUSED with the Hamiltonian energy view (WAVE_MENU=4) so the containment is visible
when you press Evolve PDE: the energy stays gathered near the core instead of diluting.
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # ~262k voxels — small for fast smoke-test

TOPOLOGY_SEED = {
    "MODE": "biaxial_hedgehog",
    "CENTER": [0.50, 0.50, 0.50],
    "R0_FRACTION": 0.06,
    "RHOC_VOXELS": 3.0,
    "BIAXIAL_DELTA": 0.30,
    "AUTO_RELAX_STEPS": 0,
    "LDG_STIFFNESS_K": 1.0,  # V ON — amplitude well at s₂*=1+δ², in units of c²/dx⁴
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "Hedgehog-Biaxial +V (1)",
        "DESCRIPTION": "Biaxial M = O·diag(1,δ,0)·Oᵀ + V(M) amplitude well — energy confinement",
    },
    "camera": {
        "INITIAL_POSITION": [1.10, 1.46, 0.81],
    },
    "universe": {
        "SIZE": [UNIVERSE_EDGE, UNIVERSE_EDGE, UNIVERSE_EDGE],
        "TARGET_VOXELS": TARGET_VOXELS,
    },
    "wave_centers": {
        "COUNT": 1,
        "POSITION": [[0.50, 0.50, 0.50]],
        "PHASE_OFFSETS_DEG": [0],
        "INIT_VELOCITY": [[0.0, 0.0, 0.0]],
        "APPLY_MOTION": False,
    },
    "ui_defaults": {
        "SHOW_AXIS": False,
        "TICK_SPACING": 0.25,
        "SHOW_GRID": False,
        "SHOW_EDGES": False,
        "FLUX_MESH_PLANES": [0.5, 0.5, 0.5],
        "SHOW_FLUX_MESH": 2,
        "WARP_MESH": 0,
        "SHOW_DIRECTORS": 2,
        "VIZ_STRIDE": 1,
        "PARTICLE_SHELL": False,
        "SHOW_GRANULES": False,
        "SIM_SPEED": 1.0,
        "PAUSED": True,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 4,  # Hamiltonian energy — watch it stay gathered under Evolve PDE
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
