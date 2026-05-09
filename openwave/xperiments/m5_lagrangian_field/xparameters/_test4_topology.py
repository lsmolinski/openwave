"""
XPERIMENT PARAMETERS — Topology Smoke Test (M5.1)

Seeds either a uniform vacuum or one/two hedgehog defects via M5.1's
`seed_vacuum` / `seed_hedgehog` kernels (validated port from Exp 2).
The simulator runs PAUSED — this xperiment exists to *visually verify*
the seeded director field before the M5.1 director-glyph renderer (task 4)
or the gradient-descent relaxation (task 6) land.

Three configurations selectable via TOPOLOGY_MODE below:
- "vacuum"        : uniform n = ẑ everywhere (sanity check)
- "hedgehog_1"  : single +1 hedgehog at the domain center
- "hedgehog_2" : +1 / −1 pair offset along x (M5.4 prelude)

Visual check (without director-glyph renderer): flux_mesh with WAVE_MENU=1
shows |ψ| ≈ 1 everywhere (since seeds are unit-length); director structure
is invisible to magnitude rendering. This xperiment is mainly a
correctness gate — it confirms the kernels run, fields are populated
without NaN/inf, and the simulator boots cleanly with the new seed paths.

The director structure becomes visible once M5.1 task 4 (director-glyph
renderer per `research_hub/3e_director_glyph_rendering.md`) lands.
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # ~262k voxels — small for fast smoke-test

# Switch the mode here:
# TOPOLOGY_MODE = "vacuum"
TOPOLOGY_MODE = "hedgehog_1"
# TOPOLOGY_MODE = "hedgehog_2"

if TOPOLOGY_MODE == "vacuum":
    TOPOLOGY_SEED = {
        "MODE": "vacuum",
    }
elif TOPOLOGY_MODE == "hedgehog_1":
    TOPOLOGY_SEED = {
        "MODE": "hedgehog",
        "DEFECTS": [
            {"CENTER": [0.50, 0.50, 0.50], "SIGN": +1},
        ],
        "DOMAIN_QUARTER_FRACTION": 0.25,  # w_vac falloff at ~D/4
    }
elif TOPOLOGY_MODE == "hedgehog_2":
    TOPOLOGY_SEED = {
        "MODE": "hedgehog",
        "DEFECTS": [
            {"CENTER": [0.30, 0.50, 0.50], "SIGN": +1},  # left, outward
            {"CENTER": [0.70, 0.50, 0.50], "SIGN": -1},  # right, inward
        ],
        "DOMAIN_QUARTER_FRACTION": 0.20,  # tighter blend so pair stays distinct
    }
else:
    raise ValueError(f"Unknown TOPOLOGY_MODE: {TOPOLOGY_MODE}")


XPARAMETERS = {
    "meta": {
        "X_NAME": f"Topo Test ({TOPOLOGY_MODE})",
        "DESCRIPTION": "Seeded director field — visual verification of seed_vacuum / seed_hedgehog",
    },
    "camera": {
        "INITIAL_POSITION": [1.20, 1.20, 1.20],
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
        "SHOW_FLUX_MESH": 0,  # off by default — directors are the primary view in M5.1
        "WARP_MESH": 0,  # no Z-warp — director field has |n|=1 uniformly, warp would be misleading
        "SHOW_DIRECTORS": 3,  # M5.1: 0=off, 1=XY, 2=+XZ, 3=all three planes
        "VIZ_STRIDE": 2,  # shared every-Nth-voxel sampling stride for directors AND granules
        "PARTICLE_SHELL": False,
        "SHOW_GRANULES": False,
        "SIM_SPEED": 1.0,
        "PAUSED": True,  # static field, no time evolution in M5.1
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 1,
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
