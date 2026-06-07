"""
XPERIMENT PARAMETERS — Dressed Hedgehog + the SIGNED constrained integrator
(M5.8.2c Option B — the real Minkowski dynamics in production)

THE MINKOWSKI-SIGNED DYNAMICS GOES LIVE. Same boost-dressed biaxial hedgehog
seed as _topo_dressed4d, but evolved under INTEGRATOR_4D "constrained" — the
spectral-projection kernel (engine2_pde solve_constrained_4d): per voxel the
faithful inertia operator A(Ṁ) = 4Σ[η[Ṁ,M_i]η,M_i] is eigendecomposed (own
10×10 Jacobi), only positive-inertia directions evolve, and the momentum is
projected onto the kept subspace every step. This is the ONLY integrator that
is stable under the signed flux (the 2c-2 structural verdict: every cheap
inertia merely delays the blow-up). Validated machine-exact against the f64
sandbox reference (sandbox_v8/m5_8_2cb_taichi_constrained.py, B-1 gates).

WHAT TO WATCH (vs the safe v1 _topo_dressed4d): the v1's wandering director
glyphs / dissipating energyF should now COHERE — orientation cohesion is the
signed dynamics' physics (2b-2: the signature coheres the clock; 2c-1: best
core retention under Minkowski). energyH stays the containment monitor.

CLOCK_KICK 0.05 starts the de Broglie clock (the demonstrated regime — the
kicked clock holds and concentrates into the core; spontaneous start is under
diagnosis, G-2c-3). Set 0.0 for the pure static-survival run.

DT_SCALE_4D 0.007 ≈ the spike-validated dt ratio (c·dt/dx = 0.0038 — the
constrained system's stiff local modes run far above the wave CFL; the 2b-1
timescale-hierarchy lesson). The kernel costs ~15 ms/step at 64³ Metal
(69 steps/s — the B-4 measurement), so the GUI stays interactive.

KM_INERTIA_4D and SIGNED_FLUX_4D are leapfrog-path devices — the constrained
integrator ignores both (it is intrinsically signed; the spectral projection
IS the exact ghost treatment, no stable-mask blend).

Boots PAUSED on the Hamiltonian view: press Evolve PDE and watch the dressed
defect HOLD with a coherent core clock.
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
    "CLOCK_KICK": 0.05,  # rad; the kicked-clock regime (0.0 = static survival)
    "DT_SCALE_4D": 0.007,  # ≈ the B-1-validated spike ratio c·dt/dx ≈ 0.0038
    "KM_INERTIA_4D": 0.0,  # UNUSED by the constrained path (faithful A instead)
    "SIGNED_FLUX_4D": True,  # informational: the constrained path is intrinsically signed
    "INTEGRATOR_4D": "constrained",  # OPTION B — the spectral-projection kernel
    "AUTO_RELAX_STEPS": 0,
    "LDG_STIFFNESS_K": 1.0,  # dressed-t* well (2c-0-validated; scaled to the 4× convention)
}


XPARAMETERS = {
    "meta": {
        "X_NAME": "Hedgehog 4D SIGNED",
        "DESCRIPTION": "Dressed hedgehog + Minkowski-SIGNED constrained integrator (Option B)",
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
        "WAVE_MENU": 4,  # Hamiltonian energy — watch containment + core clock
    },
    "analytics": {
        "INSTRUMENTATION": False,
        "EXPORT_VIDEO": False,
        "VIDEO_FRAMES": 24,
    },
    "topology_seed": TOPOLOGY_SEED,
}
