"""
UnitSystem contract and its implementations.

A UnitSystem carries every dimensional and geometric constant the engine
needs. Processors read from `ctx.data.require(UnitSystem)` and must not
contain a dimensional literal themselves; the rule is enforced by a linter.

See M4_PIPELINE_PLAN.md, Section 3, for the design rationale.

Three implementations are provided:

    NaturalUnitSystem    lambda_nu = 1, c = 1. Default for research.
    OpenWaveUnitSystem   attometres and rontoseconds, as in the legacy
                         xparameters.
    SIUnitSystem         metres and seconds. For output conversion only;
                         f32 field precision is not sufficient for
                         simulation at the soliton scale.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from openwave.common import constants

# =============================================================================
# Physical anchors, used by the conversion methods only.
# =============================================================================

_LAMBDA_PHYSICAL_M = constants.EWAVE_LENGTH   # fundamental wavelength, m
_C_PHYSICAL_MS = constants.WAVE_SPEED         # wave speed, m/s
_ATTOMETER = constants.ATTOMETER              # m per am
_RONTOSECOND = constants.RONTOSECOND          # s per rs

# Electron rest mass, CODATA 2022, kg. Local anchor: the openwave
# constants module does not expose it, and it is used only by the energy
# conversion in NaturalUnitSystem.
_M_E = 9.1093837015e-31

# Precomputed derived anchors.
_C_AMRS = _C_PHYSICAL_MS * _RONTOSECOND / _ATTOMETER   # ~0.2998 am/rs
_LAMBDA_AM = _LAMBDA_PHYSICAL_M / _ATTOMETER            # ~28.5 am
_ELECTRON_REST_ENERGY_J = _M_E * _C_PHYSICAL_MS**2      # ~8.19e-14 J

# =============================================================================
# Geometric constants.
#
# Pure numbers, independent of the unit system. All three implementations
# return the same values.
#
# Reference: manuscript v5.0.x, Section "Geometric Equation of the
# Fine-Structure Constant and the Deficit Terms"; M4.7 emergence engine.
# =============================================================================

_PI = math.pi
_SQRT2 = math.sqrt(2.0)
_SQRT3 = math.sqrt(3.0)
_EULER = math.e

# BCC packing fraction and ideal stiffness.
_ETA_BCC = _SQRT3 * _PI / 8.0
_N_IDEAL = 8.0 * _PI**4

# Lattice impedance and effective stiffness.
_ZETA = (1.0 - _ETA_BCC) / (_ETA_BCC * _N_IDEAL)
_N_GEOM = _N_IDEAL * (1.0 - _ZETA)

# Geometric core of the soliton.
_A_PI = 4.0 * _PI**3 + _PI**2 + _PI

# Magnetic deficit and its inverse.
_EPS_M = 1.0 / (_N_GEOM * _PI**3)
_GAMMA = 1.0 / _EPS_M

# Geometric fine-structure constant.
_ALPHA_INV = _A_PI - _EPS_M
_ALPHA = 1.0 / _ALPHA_INV

# Unified coupling operator and dilution factor.
_L_P_GEOM = 2.0 / _SQRT3
_K_WC = 10
_C_UNIF = 1.0 / _K_WC + 1.0 + _ALPHA / (_PI * _L_P_GEOM)
_X_EFF = _A_PI * 3.0 * _K_WC * _SQRT2 / _C_UNIF

# Neutrino radius and Planck length, from the M4.7 self-consistent chain.
# Hardcoded here as the current source of truth; a runtime provider will
# replace them once the emergence engine is packaged (plan item 1.11).
_R_NU = 2.8179354e-17      # m
_LAMBDA_L = 1.6166464e-35  # m

# Statutory and effective EMC densities (dimensionless counts).
_N_NU_STAT = (_R_NU / (2.0 * _LAMBDA_L * _EULER))**3
_N_NU_EFF = _N_NU_STAT / _X_EFF


# =============================================================================
# Contract
# =============================================================================


@runtime_checkable
class UnitSystem(Protocol):
    """
    A unit system for M4 simulations.

    Every field is a dimensional or geometric constant the engine needs.
    Every method converts an engine quantity to SI for output. Processors
    must read all dimensional values through this interface; a processor
    that writes a dimensional literal directly is a bug.
    """

    # --- Core constants ---------------------------------------------------
    c: float
    """Wave speed in these units."""

    wavelength: float
    """Fundamental wavelength (lambda_nu) in these units."""

    dx: float
    """Grid step in these units."""

    dt: float
    """Time step in these units. Must respect the 3D CFL bound."""

    rho_0: float
    """Statutory background density (N_nu,stat) in these units."""

    # --- Geometric constants, unit-independent ----------------------------
    A_pi: float
    """Geometric core of the soliton: 4*pi^3 + pi^2 + pi."""

    eps_M: float
    """Magnetic deficit: 1 / (N_geom * pi^3)."""

    N_geom: float
    """Effective BCC stiffness: 8*pi^4 * (1 - zeta)."""

    gamma: float
    """Nonlinear coupling: 1 / eps_M."""

    X_eff: float
    """Geometric dilution factor."""

    N_nu_eff: float
    """Effective volume deficit inside a soliton."""

    # --- Conversion to SI -------------------------------------------------
    def to_physical_length(self, x: float) -> float:
        """Convert a length in engine units to metres."""
        ...

    def to_physical_time(self, t: float) -> float:
        """Convert a time in engine units to seconds."""
        ...

    def to_physical_energy(self, E: float) -> float:
        """Convert an energy in engine units to joules."""
        ...

    def to_physical_density(self, rho: float) -> float:
        """Convert a number density in engine units to 1/m^3."""
        ...


# =============================================================================
# Shared geometric constants
# =============================================================================


class _GeometricMixin:
    """
    Geometric constants shared by every unit system.

    These are pure numbers derived from the BCC lattice geometry and the
    manuscript's fine-structure derivation. They do not depend on the
    choice of unit system.
    """

    @property
    def A_pi(self) -> float:
        return _A_PI

    @property
    def eps_M(self) -> float:
        return _EPS_M

    @property
    def N_geom(self) -> float:
        return _N_GEOM

    @property
    def gamma(self) -> float:
        return _GAMMA

    @property
    def X_eff(self) -> float:
        return _X_EFF

    @property
    def N_nu_eff(self) -> float:
        return _N_NU_EFF


# =============================================================================
# Natural units
# =============================================================================


@dataclass(frozen=True)
class NaturalUnitSystem(_GeometricMixin):
    """
    Natural units: lambda_nu = 1, c = 1.

    Every engine quantity is a dimensionless multiple of the fundamental
    wavelength and the wave speed. The grid step is chosen to give a
    configurable number of voxels per wavelength.

    Parameters
    ----------
    grid_voxels_per_lambda : int, default 20
        Grid resolution. 20 voxels per lambda is a reasonable default for
        a second-order-accurate leapfrog scheme. The minimum enforced is 2.
    cfl_safety : float, default 0.9
        Fraction of the 3D CFL bound to use. The 3D leapfrog stability
        bound for a scalar wave equation is dt <= dx / (c * sqrt(3)); the
        safety factor keeps the run strictly below it.
    """

    grid_voxels_per_lambda: int = 20
    cfl_safety: float = 0.9

    def __post_init__(self) -> None:
        if self.grid_voxels_per_lambda < 2:
            raise ValueError(
                "grid_voxels_per_lambda must be >= 2, got "
                f"{self.grid_voxels_per_lambda}"
            )
        if not 0.0 < self.cfl_safety <= 1.0:
            raise ValueError(
                f"cfl_safety must be in (0, 1], got {self.cfl_safety}"
            )

    @property
    def c(self) -> float:
        return 1.0

    @property
    def wavelength(self) -> float:
        return 1.0

    @property
    def dx(self) -> float:
        return 1.0 / self.grid_voxels_per_lambda

    @property
    def dt(self) -> float:
        return self.cfl_safety * self.dx / (self.c * _SQRT3)

    @property
    def rho_0(self) -> float:
        return _N_NU_STAT

    def to_physical_length(self, x: float) -> float:
        return x * _LAMBDA_PHYSICAL_M

    def to_physical_time(self, t: float) -> float:
        return t * _LAMBDA_PHYSICAL_M / _C_PHYSICAL_MS

    def to_physical_energy(self, E: float) -> float:
        # Convention: one engine energy unit is the electron rest energy.
        # The physical interpretation of the engine's energy variable
        # depends on field normalization, which is set elsewhere.
        return E * _ELECTRON_REST_ENERGY_J

    def to_physical_density(self, rho: float) -> float:
        # One engine density unit is one per cubic wavelength.
        return rho / _LAMBDA_PHYSICAL_M**3


# =============================================================================
# OpenWave units
# =============================================================================


@dataclass(frozen=True)
class OpenWaveUnitSystem(_GeometricMixin):
    """
    Attometres and rontoseconds, as in the legacy OpenWave xparameters.

    Length in attometres (am), time in rontoseconds (rs). The wave speed
    is c ~ 0.3 am/rs, which equals 3e8 m/s in SI.

    Parameters
    ----------
    grid_voxels_per_lambda : int, default 12
        Grid resolution. 12 voxels per lambda is the legacy default.
    cfl_safety : float, default 0.9
        Fraction of the 3D CFL bound to use.
    """

    grid_voxels_per_lambda: int = 12
    cfl_safety: float = 0.9

    def __post_init__(self) -> None:
        if self.grid_voxels_per_lambda < 2:
            raise ValueError(
                "grid_voxels_per_lambda must be >= 2, got "
                f"{self.grid_voxels_per_lambda}"
            )
        if not 0.0 < self.cfl_safety <= 1.0:
            raise ValueError(
                f"cfl_safety must be in (0, 1], got {self.cfl_safety}"
            )

    @property
    def c(self) -> float:
        return _C_AMRS

    @property
    def wavelength(self) -> float:
        return _LAMBDA_AM

    @property
    def dx(self) -> float:
        return self.wavelength / self.grid_voxels_per_lambda

    @property
    def dt(self) -> float:
        return self.cfl_safety * self.dx / (self.c * _SQRT3)

    @property
    def rho_0(self) -> float:
        return _N_NU_STAT

    def to_physical_length(self, x: float) -> float:
        return x * _ATTOMETER

    def to_physical_time(self, t: float) -> float:
        return t * _RONTOSECOND

    def to_physical_energy(self, E: float) -> float:
        # One engine energy unit is one attojoule.
        return E * 1e-18

    def to_physical_density(self, rho: float) -> float:
        # One engine density unit is one per cubic attometre.
        return rho / _ATTOMETER**3


# =============================================================================
# SI units
# =============================================================================


@dataclass(frozen=True)
class SIUnitSystem(_GeometricMixin):
    """
    Metres and seconds.

    For output conversion and cross-checking against CODATA. Not suitable
    for simulation: f32 fields lose precision at ~1e-17 m, the natural
    scale of the soliton core.

    Parameters
    ----------
    grid_voxels_per_lambda : int, default 20
        Grid resolution. 20 voxels per lambda is a reasonable default.
    cfl_safety : float, default 0.9
        Fraction of the 3D CFL bound to use.
    """

    grid_voxels_per_lambda: int = 20
    cfl_safety: float = 0.9

    def __post_init__(self) -> None:
        if self.grid_voxels_per_lambda < 2:
            raise ValueError(
                "grid_voxels_per_lambda must be >= 2, got "
                f"{self.grid_voxels_per_lambda}"
            )
        if not 0.0 < self.cfl_safety <= 1.0:
            raise ValueError(
                f"cfl_safety must be in (0, 1], got {self.cfl_safety}"
            )

    @property
    def c(self) -> float:
        return _C_PHYSICAL_MS

    @property
    def wavelength(self) -> float:
        return _LAMBDA_PHYSICAL_M

    @property
    def dx(self) -> float:
        return self.wavelength / self.grid_voxels_per_lambda

    @property
    def dt(self) -> float:
        return self.cfl_safety * self.dx / (self.c * _SQRT3)

    @property
    def rho_0(self) -> float:
        return _N_NU_STAT

    def to_physical_length(self, x: float) -> float:
        return x

    def to_physical_time(self, t: float) -> float:
        return t

    def to_physical_energy(self, E: float) -> float:
        return E

    def to_physical_density(self, rho: float) -> float:
        return rho