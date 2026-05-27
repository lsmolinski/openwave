"""
XPERIMENT PARAMETERS — Biaxial Hedgehog Smoke Test (M5.6.5a)

Seeds a single BIAXIAL hedgehog matrix field via M5.6.5a's
`seed_biaxial_hedgehog_M` kernel. The biaxial hedgehog is

    M(x) = O(x) · D(s(r)) · O(x)ᵀ ,  O = [r̂ | e_Θ | e_Φ] ,  D = diag(1, δ, 0)

i.e. three DISTINCT eigenvalues carried by an orthonormal hedgehog frame.
This is the M5.6.2 sandbox structure ported to production. Contrast with
the uniaxial `seed_hedgehog` (M5.1), where D = diag(1, ½, ½) is degenerate
and the field reduces to a single director n̂ = r̂.

Why biaxial matters (M5.6.2 finding):
- The matrix connection C_μν = [M_μ, M_ν] is NON-zero for a biaxial frame
  (it is identically zero for the single-generator uniaxial hedgehog).
  C is the matrix-level mass source — ‖C‖ ∝ 1/r² matches the scalar
  geometric KG mass of M5.6.1. So the biaxial hedgehog is the seed that
  dynamically SOURCES its own twist (ψ = 0 is NOT static) → the M5.8
  Zitterbewegung-clock seed.

Two melt scales (both regularize the otherwise-singular seed):
- r0  (radial)        : eigenvalues melt to isotropic D→(1+δ)/3·I as r→0.
                         Set via R0_FRACTION × max(grid). This is the
                         M5.6.3 Faber core-size knob (the physical mass
                         handle for M5.9 calibration).
- ρc  (disclination)  : the e_Φ azimuthal frame field winds like 1/ρ on
                         the z-axis (hairy-ball). A clamped smoothstep
                         melts biaxiality inside ρ < ρc so the frame stays
                         orthonormal in the bulk but is regular on-axis.

NO auto-relax: the biaxial M is built directly from the frame + melt.
The M5.1 relax_director_step rebuilds M *uniaxially* from the director,
which would destroy the biaxial structure — so AUTO_RELAX_STEPS is forced
to 0 for this mode (the launcher also skips relax for biaxial_hedgehog).

Runs PAUSED — visual correctness gate. The biaxial-ellipsoid glyph
renderer (M5.6.5b) makes the three distinct eigenvalues visible; until
then SHOW_DIRECTORS shows the principal (largest-eigenvalue) axis, which
the analytic Cardano eigensolver (M5.6.5a fix) now recovers correctly
on Metal/f32 even for biaxial M (ti.sym_eig was catastrophically wrong
for the non-degenerate case).
"""

UNIVERSE_EDGE = 1e-15  # m
TARGET_VOXELS = 64**3  # ~262k voxels — small for fast smoke-test

TOPOLOGY_SEED = {
    "MODE": "biaxial_hedgehog",
    "CENTER": [0.50, 0.50, 0.50],  # domain center (fractional)
    "R0_FRACTION": 0.06,  # radial eigenvalue-melt core, × max(grid) → voxels
    "RHOC_VOXELS": 3.0,  # z-axis disclination melt radius, voxels
    "BIAXIAL_DELTA": 0.30,  # middle eigenvalue δ → D = diag(1, 0.30, 0)
    "AUTO_RELAX_STEPS": 0,  # forced 0 — relax would re-uniaxialize M
}


XPARAMETERS = {
    "meta": {
        "X_NAME": f"Hedgehog-Biaxial (1)",
        "DESCRIPTION": "Seeded biaxial matrix field M = O·diag(1,δ,0)·Oᵀ — visual verification",
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
        "SHOW_FLUX_MESH": 1,
        "WARP_MESH": 0,
        "SHOW_DIRECTORS": 3,  # all three planes — principal axis of the biaxial frame
        "VIZ_STRIDE": 1,
        "PARTICLE_SHELL": False,
        "SHOW_GRANULES": False,
        "SIM_SPEED": 1.0,
        "PAUSED": True,
    },
    "color_defaults": {
        "COLOR_THEME": "OCEAN",
        "WAVE_MENU": 5,
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
