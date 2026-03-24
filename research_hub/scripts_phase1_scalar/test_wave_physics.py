"""
Physics invariant tests for 1D wave sandbox — Phase 0 validation.

Tests cover:
- Phasor RMS boundary limits (center, far-field)
- Single WC sinc envelope shape
- Same/opposite charge interference patterns
- Symmetry properties
- Energy density dimensional consistency
- Force field consistency
- EMA-RMS equivalence with phasor RMS
- Weight function behavior

Run from project root:
    pytest openwave/xperiments/m3_wolff_lafreniere/research/FORCE_UNIFICATION/test_wave_physics.py -v
"""

import numpy as np
import pytest

import wave_engine_1D_v2 as we

# ================================================================
# Aliases
# ================================================================

A0 = we.A0_am
lam = we.lam_am
k = we.k
omega = we.omega
f_rHz = we.f_rHz
rho = we.rho_qgam


# ================================================================
# Helpers
# ================================================================


@pytest.fixture(autouse=True)
def _reset_wave_centers():
    """Save and restore module-level wave_centers after each test."""
    original = we.wave_centers.copy()
    yield
    we.wave_centers.clear()
    we.wave_centers.extend(original)


def set_wcs(centers):
    """Replace module-level wave_centers with given list."""
    we.wave_centers.clear()
    we.wave_centers.extend(centers)


def domain(half_lam=10, npoints=8000):
    """Create a 1D spatial domain centered at 0."""
    half = half_lam * lam
    return np.linspace(-half, half, npoints)


# ================================================================
# Phasor RMS: Boundary Limits
# ================================================================


class TestPhasorBoundaryLimits:
    """Phasor RMS at limiting cases (center, far-field)."""

    def test_single_wc_center_amplitude(self):
        """At WC center (kr=0), RMS = A0 * 2 / sqrt(2).

        w(0)=1 → C_n = (1+1)*A*sin(kr)/kr → 2A, S_n → 0.
        Peak = 2A, RMS = 2A/sqrt(2).
        """
        x = domain(5, 4001)  # odd → x=0 exactly on grid
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0, amplitude=A0)])
        rms = we.compute_phasor_rms(x)

        center_idx = np.argmin(np.abs(x))
        expected = A0 * 2.0 / np.sqrt(2.0)
        np.testing.assert_allclose(rms[center_idx], expected, rtol=1e-6)

    def test_center_rms_phase_independent(self):
        """Center RMS magnitude is the same regardless of source_offset."""
        x = domain(5, 4001)  # odd → x=0 exactly on grid
        center_idx = np.argmin(np.abs(x))
        expected = A0 * 2.0 / np.sqrt(2.0)

        for phase in [0.0, np.pi, np.pi / 2, np.pi / 4]:
            set_wcs([we.WaveCenter(x_am=0.0, phase=phase, amplitude=A0)])
            rms = we.compute_phasor_rms(x)
            np.testing.assert_allclose(
                rms[center_idx],
                expected,
                rtol=1e-6,
                err_msg=f"phase={phase:.3f}",
            )

    def test_farfield_envelope(self):
        """Far from WC (w→0), RMS → A/(kr*sqrt(2)).

        C_n ≈ A*sin(kr)/kr, S_n ≈ -A*cos(kr)/kr
        Peak = sqrt(C² + S²) = A/kr, RMS = A/(kr*sqrt(2)).
        """
        x = domain(20, 16000)
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0, amplitude=A0)])
        rms = we.compute_phasor_rms(x)

        # Sample well beyond transition (TRANSITION_LAM = 1.25λ)
        for r_lam in [8, 10, 12, 15]:
            r = r_lam * lam
            idx = np.argmin(np.abs(x - r))
            kr = k * r
            expected_rms = A0 / (kr * np.sqrt(2.0))
            np.testing.assert_allclose(
                rms[idx],
                expected_rms,
                rtol=0.05,
                err_msg=f"r={r_lam}λ",
            )

    def test_rms_nonnegative_single_wc(self):
        """Phasor RMS >= 0 everywhere (single WC)."""
        x = domain(10, 8000)
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0)])
        rms = we.compute_phasor_rms(x)
        assert np.all(rms >= 0), "RMS contains negative values"

    def test_rms_nonnegative_two_wc(self):
        """Phasor RMS >= 0 everywhere (two opposite-charge WCs)."""
        x = domain(10, 8000)
        set_wcs([
            we.WaveCenter(x_am=-2 * lam, phase=0.0),
            we.WaveCenter(x_am=+2 * lam, phase=np.pi),
        ])
        rms = we.compute_phasor_rms(x)
        assert np.all(rms >= 0), "RMS contains negative values"


# ================================================================
# Sinc Envelope Shape
# ================================================================


class TestSincEnvelope:
    """Single WC produces expected sinc-like envelope."""

    def test_near_field_nodes(self):
        """Near WC center (w≈1), standing wave nodes at kr = nπ (r = nλ/2).

        At these nodes sin(kr)=0, so RMS should dip sharply.
        """
        x = domain(5, 20000)
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0, amplitude=A0)])
        rms = we.compute_phasor_rms(x)

        center_rms = rms[np.argmin(np.abs(x))]
        for n in [1, 2]:
            r_node = n * lam / 2
            idx = np.argmin(np.abs(x - r_node))
            ratio = rms[idx] / center_rms
            assert ratio < 0.15, (
                f"Node at r={n}λ/2: ratio to center = {ratio:.4f}, expected < 0.15"
            )

    def test_farfield_peaks_decrease(self):
        """Far-field sinc lobe peaks decrease with distance (1/kr envelope)."""
        x = domain(20, 16000)
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0, amplitude=A0)])
        rms = we.compute_phasor_rms(x)

        # Between sinc zeros: peaks at roughly r = (n+0.75)λ
        peak_positions = [2.75, 4.75, 6.75, 8.75]
        peaks = []
        for r_lam in peak_positions:
            idx = np.argmin(np.abs(x - r_lam * lam))
            peaks.append(rms[idx])

        for i in range(1, len(peaks)):
            assert peaks[i] < peaks[i - 1], (
                f"Peak at {peak_positions[i]}λ ({peaks[i]:.6e}) "
                f">= peak at {peak_positions[i-1]}λ ({peaks[i-1]:.6e})"
            )


# ================================================================
# Charge Interference Patterns
# ================================================================


class TestChargeInterference:
    """Constructive/destructive interference from charge phase."""

    def _midpoint_rms(self, phase1, phase2, sep_lam=6.0):
        """RMS at midpoint between two WCs separated by sep_lam."""
        x = domain(15, 12000)
        sep = sep_lam * lam
        set_wcs([
            we.WaveCenter(x_am=-sep / 2, phase=phase1),
            we.WaveCenter(x_am=+sep / 2, phase=phase2),
        ])
        rms = we.compute_phasor_rms(x)
        return rms[np.argmin(np.abs(x))]

    def _single_wc_rms_at(self, distance_am):
        """RMS from a single WC at given distance."""
        x = domain(15, 12000)
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0)])
        rms = we.compute_phasor_rms(x)
        return rms[np.argmin(np.abs(x - distance_am))]

    def test_same_charge_midpoint_higher(self):
        """Same-charge: P adds constructively → midpoint > single WC."""
        mid = self._midpoint_rms(0.0, 0.0, sep_lam=6.0)
        single = self._single_wc_rms_at(3.0 * lam)
        assert mid > single, (
            f"Same-charge midpoint ({mid:.6e}) should be > single WC ({single:.6e})"
        )

    def test_opposite_charge_midpoint_lower(self):
        """Opposite-charge: P cancels → midpoint < single WC."""
        mid = self._midpoint_rms(0.0, np.pi, sep_lam=6.0)
        single = self._single_wc_rms_at(3.0 * lam)
        assert mid < single, (
            f"Opposite-charge midpoint ({mid:.6e}) should be < single WC ({single:.6e})"
        )

    def test_same_higher_than_opposite(self):
        """At midpoint, same-charge RMS > opposite-charge RMS."""
        same = self._midpoint_rms(0.0, 0.0, sep_lam=6.0)
        opp = self._midpoint_rms(0.0, np.pi, sep_lam=6.0)
        assert same > opp, (
            f"Same ({same:.6e}) should be > opposite ({opp:.6e})"
        )


# ================================================================
# Symmetry
# ================================================================


class TestSymmetry:
    """Symmetry properties of the wave field."""

    def test_same_charge_symmetric_rms(self):
        """Two same-charge WCs at ±d/2 → RMS symmetric about x=0."""
        x = domain(10, 8000)
        sep = 5 * lam
        set_wcs([
            we.WaveCenter(x_am=-sep / 2, phase=0.0),
            we.WaveCenter(x_am=+sep / 2, phase=0.0),
        ])
        rms = we.compute_phasor_rms(x)
        np.testing.assert_allclose(rms, rms[::-1], atol=1e-12)

    def test_opposite_charge_symmetric_rms(self):
        """Two opposite-charge WCs at ±d/2 → RMS still symmetric (magnitude)."""
        x = domain(10, 8000)
        sep = 5 * lam
        set_wcs([
            we.WaveCenter(x_am=-sep / 2, phase=0.0),
            we.WaveCenter(x_am=+sep / 2, phase=np.pi),
        ])
        rms = we.compute_phasor_rms(x)
        np.testing.assert_allclose(rms, rms[::-1], atol=1e-12)

    def test_single_wc_symmetric(self):
        """Single WC at origin → symmetric RMS."""
        x = domain(10, 8000)
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0)])
        rms = we.compute_phasor_rms(x)
        np.testing.assert_allclose(rms, rms[::-1], atol=1e-12)


# ================================================================
# Energy Density
# ================================================================


class TestEnergyDensity:
    """Energy density E = ρ·V·(f·A)²."""

    def test_nonnegative(self):
        """E >= 0 everywhere."""
        x = domain(10, 8000)
        set_wcs([
            we.WaveCenter(x_am=-2 * lam, phase=0.0),
            we.WaveCenter(x_am=+2 * lam, phase=np.pi),
        ])
        rms = we.compute_phasor_rms(x)
        energy = we.compute_energy_density(rms)
        assert np.all(energy >= 0)

    def test_zero_with_no_wc(self):
        """No wave centers → zero energy."""
        x = domain(5, 4000)
        set_wcs([])
        rms = we.compute_phasor_rms(x)
        energy = we.compute_energy_density(rms)
        np.testing.assert_allclose(energy, 0.0, atol=1e-30)

    def test_peak_at_wc_center(self):
        """Energy peaks at WC center (highest amplitude)."""
        x = domain(5, 4000)
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0)])
        rms = we.compute_phasor_rms(x)
        energy = we.compute_energy_density(rms)

        center_idx = np.argmin(np.abs(x))
        assert energy[center_idx] == np.max(energy)

    def test_scales_with_amplitude_squared(self):
        """E ∝ A² — doubling amplitude quadruples energy."""
        x = domain(5, 4000)
        center_idx = np.argmin(np.abs(x))

        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0, amplitude=A0)])
        e1 = we.compute_energy_density(we.compute_phasor_rms(x))[center_idx]

        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0, amplitude=2 * A0)])
        e2 = we.compute_energy_density(we.compute_phasor_rms(x))[center_idx]

        np.testing.assert_allclose(e2 / e1, 4.0, rtol=1e-6)


# ================================================================
# Force Field
# ================================================================


class TestForceField:
    """Force field F = -∇E where E = ρ·V·(f·A)²."""

    def test_zero_at_symmetric_center(self):
        """Single WC at origin → ~zero force at center (gradient of symmetric peak)."""
        x = domain(10, 8000)
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0)])
        rms = we.compute_phasor_rms(x)
        force = we.compute_force_field(rms)

        center_idx = np.argmin(np.abs(x))
        f_max = np.max(np.abs(force))
        assert abs(force[center_idx]) < f_max * 0.01, (
            f"Force at center = {force[center_idx]:.3e}, max = {f_max:.3e}"
        )

    def test_antisymmetric_single_wc(self):
        """Single WC at origin: F(x) = -F(-x) (gradient of symmetric fn)."""
        x = domain(10, 8000)
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0)])
        rms = we.compute_phasor_rms(x)
        force = we.compute_force_field(rms)

        np.testing.assert_allclose(
            force, -force[::-1], atol=1e-10 * np.max(np.abs(force))
        )

    def test_antisymmetric_two_wc(self):
        """Two WCs at ±d/2: force antisymmetric about midpoint."""
        x = domain(10, 8000)
        sep = 5 * lam
        set_wcs([
            we.WaveCenter(x_am=-sep / 2, phase=0.0),
            we.WaveCenter(x_am=+sep / 2, phase=np.pi),
        ])
        rms = we.compute_phasor_rms(x)
        force = we.compute_force_field(rms)

        np.testing.assert_allclose(
            force, -force[::-1], atol=1e-10 * np.max(np.abs(force))
        )


# ================================================================
# EMA-RMS Equivalence
# ================================================================


class TestEmaRmsEquivalence:
    """Phasor RMS matches time-averaged sqrt(<ψ²>) over one full period."""

    @staticmethod
    def _time_averaged_rms(x, n_samples=200):
        """Compute RMS by sampling displacement over one full period."""
        psi_sq_sum = np.zeros_like(x)
        for i in range(n_samples):
            t = (i / n_samples) * we.period_rs
            psi = we.compute_displacement(x, t)
            psi_sq_sum += psi**2
        return np.sqrt(psi_sq_sum / n_samples)

    def test_single_wc(self):
        """Single WC: phasor RMS ≈ time-averaged RMS."""
        x = domain(8, 4000)
        set_wcs([we.WaveCenter(x_am=0.0, phase=0.0)])

        rms_phasor = we.compute_phasor_rms(x)
        rms_time = self._time_averaged_rms(x)

        mask = rms_phasor > 0.01 * np.max(rms_phasor)
        np.testing.assert_allclose(rms_time[mask], rms_phasor[mask], rtol=0.02)

    def test_two_wc_opposite(self):
        """Two opposite-charge WCs: phasor RMS ≈ time-averaged RMS."""
        x = domain(10, 4000)
        sep = 5 * lam
        set_wcs([
            we.WaveCenter(x_am=-sep / 2, phase=0.0),
            we.WaveCenter(x_am=+sep / 2, phase=np.pi),
        ])

        rms_phasor = we.compute_phasor_rms(x)
        rms_time = self._time_averaged_rms(x)

        mask = rms_phasor > 0.01 * np.max(rms_phasor)
        np.testing.assert_allclose(rms_time[mask], rms_phasor[mask], rtol=0.02)

    def test_two_wc_same(self):
        """Two same-charge WCs: phasor RMS ≈ time-averaged RMS."""
        x = domain(10, 4000)
        sep = 5 * lam
        set_wcs([
            we.WaveCenter(x_am=-sep / 2, phase=0.0),
            we.WaveCenter(x_am=+sep / 2, phase=0.0),
        ])

        rms_phasor = we.compute_phasor_rms(x)
        rms_time = self._time_averaged_rms(x)

        mask = rms_phasor > 0.01 * np.max(rms_phasor)
        np.testing.assert_allclose(rms_time[mask], rms_phasor[mask], rtol=0.02)


# ================================================================
# Weight Function
# ================================================================


class TestWeightFunction:
    """Weight function w(r) = 1/(1 + (r/(transition·λ))^8)."""

    def test_center_is_one(self):
        """w(0) = 1."""
        w = we.compute_weight(np.array([0.0]))
        np.testing.assert_allclose(w[0], 1.0, rtol=1e-10)

    def test_transition_is_half(self):
        """w(transition·λ) = 1/(1+1^8) = 0.5."""
        r = np.array([we.TRANSITION_LAM * lam])
        w = we.compute_weight(r)
        np.testing.assert_allclose(w[0], 0.5, rtol=1e-10)

    def test_far_approaches_zero(self):
        """w(100λ) ≈ 0."""
        w = we.compute_weight(np.array([100 * lam]))
        assert w[0] < 1e-10

    def test_monotonically_decreasing(self):
        """w(r) decreases monotonically with r."""
        r = np.linspace(0, 20 * lam, 1000)
        w = we.compute_weight(r)
        assert np.all(np.diff(w) <= 0)
