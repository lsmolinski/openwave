# ENVELOPE FIELD for force calculation (smooth 1/r, no oscillations)
# Tracks signed amplitude envelope: sum of (charge_sign * A₀/r) from each source
# This gives smooth 1/r² force law matching predictions
# self.amp_local_envelope_am = ti.field(dtype=ti.f32, shape=grid_size)  # am, signed


# ================================================================
# OLD - ANALYTICAL SIGNED AMPLITUDE ENVELOPE
# Particles don't respond to 10²⁵ Hz oscillation frequencies.
# Particle's mass (inertia) acts as a low-pass filter, averaging out the rapid
# oscillations and responding only to the time-averaged energy-density (envelope).
# This envelope drives the force calculations, computed directly from wave functions.
# Applies superposition principle for multiple wave-centers, with signed charge sign.
# Avoids computationally expensive real-time tracking methods (RMS, zero-crossing).
# Also avoids instability from real-time EMS calculations of moving wave-centers.
# ================================================================
# Charge sign: cos(0)=+1 (eg: positron), cos(π)=-1 (eg: electron)
# charge_sign = ti.cos(source_offset)

# TODO: Investigate why these constants work well for envelope
# golden_ratio = (1 + ti.sqrt(5)) / 2  # ~1.6180339887
# weight_factor = 2.0 * ti.math.pi**2  # ~19.7392088, decay, damping
# offset_factor = 1 / (wavelength_grid * golden_ratio)

# # SPIKED 1/r ==================================
# if r_grid < 0.5:  # CENTER VOXEL only, avoids singularity
#     trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#         base_amplitude_am * wave_field.scale_factor * k_grid * 2
#     )  # finite value at center, k_grid / constant + offset
# else:  # FAR-FIELD: smooth 1/r decay
#     trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#         base_amplitude_am * wave_field.scale_factor * 1.0 / r_grid
#     )  # spiked 1/r decay

# # SMOOTHED 1/r ==================================
# trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#     base_amplitude_am
#     * wave_field.scale_factor
#     * k_grid
#     / ti.sqrt((k_grid * r_grid) ** 2 + 1)
# )  # smoothed 1/r decay

# # DAMPED SMOOTHED 1/r ==================================
# trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#     base_amplitude_am
#     * wave_field.scale_factor
#     * k_grid
#     / ti.sqrt((k_grid * r_grid) ** 2 + (2 * ti.math.pi) ** 2)
# )  # smoothed 1/r decay

# # WOLFF-ORIGINAL ENVELOPE ==================================
# if r_grid < 0.5:  # CENTER VOXEL only, avoids singularity
#     trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#         base_amplitude_am * wave_field.scale_factor * k_grid
#     )  # finite value at center, k_grid
# else:
#     trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#         base_amplitude_am * wave_field.scale_factor * ti.sin(k_grid * r_grid) / r_grid
#     )  # standing-smooth sin(kr)/r decay

# # WOLFF ONLY AT NEAR-FIELD ==================================
# if r_grid < 0.5:  # CENTER VOXEL only, avoids singularity
#     trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#         base_amplitude_am * wave_field.scale_factor * k_grid
#     )  # finite value at center, k_grid
# else:
#     if r_grid <= (2.5 * ti.math.pi / k_grid):  # NEAR-FIELD: time-dilated-1.25λ
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am
#             * wave_field.scale_factor
#             * ti.sin(k_grid * r_grid)
#             / r_grid
#         )  # standing-smooth sin(kr)/r decay
#     else:  # FAR-FIELD: smooth 1/r decay
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am * wave_field.scale_factor * 1.0 / r_grid
#         )  # smooth 1/r decay

# # ABS WOLFF ONLY NEAR-FIELD ==================================
# if r_grid < 0.5:  # CENTER VOXEL only, avoids singularity
#     trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#         base_amplitude_am * wave_field.scale_factor * k_grid
#     )  # finite value at center, k_grid
# else:
#     if r_grid <= (2.5 * ti.math.pi / k_grid):  # NEAR-FIELD: time-dilated-1.25λ
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am
#             * wave_field.scale_factor
#             * ti.abs(ti.sin(k_grid * r_grid))
#             / r_grid
#         )  # standing-smooth sin(kr)/r decay
#     else:  # FAR-FIELD: smooth 1/r decay
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am * wave_field.scale_factor * 1.0 / r_grid
#         )  # smooth 1/r decay

# # DAMPED+OFFSET WOLFF NEAR-FIELD ==================================
# if r_grid < 0.5:  # CENTER VOXEL only, avoids singularity
#     trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#         base_amplitude_am * wave_field.scale_factor * k_grid / (2 * ti.math.pi)
#         + 1 / (wavelength_grid * golden_ratio)
#     )  # finite value at center, k_grid / constant + offset
# else:
#     if r_grid <= (2.5 * ti.math.pi / k_grid):  # NEAR-FIELD: time-dilated-1.25λ
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am
#             * wave_field.scale_factor
#             * ti.sin(k_grid * r_grid)
#             / (r_grid * 2 * ti.math.pi)
#             + 1 / (wavelength_grid * golden_ratio)
#         )  # standing-smooth sin(kr)/(r·constant) + offset decay
#     else:  # FAR-FIELD: smooth 1/r decay
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am * wave_field.scale_factor * 1.0 / r_grid
#         )  # smooth 1/r decay

# # ABS DAMPED+OFFSET WOLFF NEAR-FIELD ==================================
# if r_grid < 0.5:  # CENTER VOXEL only, avoids singularity
#     trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#         base_amplitude_am * wave_field.scale_factor * k_grid / (2)
#         + 1 / (wavelength_grid * golden_ratio)
#     )  # finite value at center, k_grid / constant + offset
# else:
#     if r_grid <= (2.5 * ti.math.pi / k_grid):  # NEAR-FIELD: time-dilated-1.5λ
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am
#             * wave_field.scale_factor
#             * ti.abs(ti.sin(k_grid * r_grid))
#             / (r_grid * 2)
#             + 1 / (wavelength_grid * golden_ratio)
#         )  # standing-smooth sin(kr)/(r·constant) + offset decay
#     else:  # FAR-FIELD: smooth 1/r decay
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am * wave_field.scale_factor * 1.0 / r_grid
#         )  # smooth 1/r decay

# # DAMPED + WOLFF ==================================
# if r_grid < 0.5:  # CENTER VOXEL only, avoids singularity
#     trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#         base_amplitude_am * wave_field.scale_factor * (k_grid / (2 * ti.math.pi))
#     )  # finite value at center, k_grid / constant
# else:
#     if r_grid <= (1.25 * 2 * ti.math.pi / k_grid):  # NEAR-FIELD: time-dilated-1.25λ
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am
#             * wave_field.scale_factor
#             * (
#                 k_grid / ti.sqrt((k_grid * r_grid) ** 2 + (48))
#                 + ti.sin(k_grid * r_grid) / (r_grid * 4)
#             )
#         )  # smoothed 1/r decay
#     else:  # FAR-FIELD: smooth 1/r decay
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am * wave_field.scale_factor * 1.0 / r_grid
#         )  # smooth 1/r decay

# # FLAT NEAR-FIELD ==================================
# if r_grid < 0.5:  # CENTER VOXEL only, avoids singularity
#     trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#         base_amplitude_am * wave_field.scale_factor * k_grid / (2 * ti.math.pi * 1.25)
#     )  # finite value at center, k_grid
# else:
#     if r_grid <= (2.5 * ti.math.pi / k_grid):  # NEAR-FIELD: time-dilated-1.25λ
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am
#             * wave_field.scale_factor
#             * k_grid
#             / (2 * ti.math.pi * 1.25)
#         )  # standing-smooth sin(kr)/r decay
#     else:  # FAR-FIELD: smooth 1/r decay
#         trackers.amp_local_envelope_am[i, j, k] += charge_sign * (
#             base_amplitude_am * wave_field.scale_factor * 1.0 / r_grid
#         )  # smooth 1/r decay
