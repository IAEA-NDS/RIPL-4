# -*- coding: utf-8 -*-
"""Core dataclasses for optical model potential components.

This module defines the fundamental data structures used to represent
optical model potentials, including coefficient arrays, energy ranges,
potential components, and Coulomb parameters.
"""

from dataclasses import dataclass as _dataclass, field as _field
from typing import List, Optional, ClassVar as _ClassVar

from . import config


@_dataclass
class RadiusCoefficients:
    """Coefficients for calculating the radius parameter R.

    The radius R(i,j) is calculated as:
        R = {rco[0] + rco[2]*eta + rco[3]/A + rco[4]/sqrt(A) + ...} * A^(1/3)

    See om-parameter-u.readme for full formula details.

    Attributes:
        values: Array of 13 coefficients (rco[k], k=0-12)
    """
    values: List[float] = _field(default_factory=lambda: [0.0] * config.NUM_RADIUS_COEFFS)

    _field_info: _ClassVar[dict] = {
        'values': 'Radius coefficients rco(k), k=1-13 [various units]',
    }

    def __post_init__(self):
        if len(self.values) != config.NUM_RADIUS_COEFFS:
            raise ValueError(f"Expected {config.NUM_RADIUS_COEFFS} radius coefficients, got {len(self.values)}")

    @property
    def r0(self) -> float:
        """Base radius coefficient r0 [fm]."""
        return abs(self.values[0])

    @property
    def energy_dep(self) -> float:
        """Energy dependence coefficient."""
        return self.values[1]

    @property
    def eta_dep(self) -> float:
        """Asymmetry (eta) dependence coefficient."""
        return self.values[2]


@_dataclass
class DiffusenessCoefficients:
    """Coefficients for calculating the diffuseness parameter a.

    The diffuseness a(i,j) is calculated as:
        a = aco[0] + aco[1]*E + aco[2]*eta + aco[3]/A + ...

    See om-parameter-u.readme for full formula details.

    Attributes:
        values: Array of 13 coefficients (aco[k], k=0-12)
    """
    values: List[float] = _field(default_factory=lambda: [0.0] * config.NUM_DIFFUSENESS_COEFFS)

    _field_info: _ClassVar[dict] = {
        'values': 'Diffuseness coefficients aco(k), k=1-13 [various units]',
    }

    def __post_init__(self):
        if len(self.values) != config.NUM_DIFFUSENESS_COEFFS:
            raise ValueError(f"Expected {config.NUM_DIFFUSENESS_COEFFS} diffuseness coefficients, got {len(self.values)}")

    @property
    def a0(self) -> float:
        """Base diffuseness coefficient a0 [fm]."""
        return abs(self.values[0])

    @property
    def energy_dep(self) -> float:
        """Energy dependence coefficient."""
        return self.values[1]

    @property
    def eta_dep(self) -> float:
        """Asymmetry (eta) dependence coefficient."""
        return self.values[2]


@_dataclass
class PotentialStrength:
    """Coefficients for calculating the potential strength V.

    The potential strength V(i,j) can be calculated using several forms
    (standard, Smith, Varner, extended Koning) depending on pot[21-24].

    See om-parameter-u.readme for full formula details.

    Attributes:
        values: Array of 25 coefficients (pot[k], k=0-24)
    """
    values: List[float] = _field(default_factory=lambda: [0.0] * config.NUM_POTENTIAL_COEFFS)

    _field_info: _ClassVar[dict] = {
        'values': 'Potential strength coefficients pot(k), k=1-25 [various units]',
    }

    def __post_init__(self):
        if len(self.values) != config.NUM_POTENTIAL_COEFFS:
            raise ValueError(f"Expected {config.NUM_POTENTIAL_COEFFS} potential coefficients, got {len(self.values)}")

    @property
    def V0(self) -> float:
        """Base potential depth [MeV]."""
        return self.values[0]

    @property
    def energy_dep(self) -> float:
        """Linear energy dependence coefficient."""
        return self.values[1]

    @property
    def is_smith_form(self) -> bool:
        """Check if using Smith parameterization."""
        return self.values[21] != 0.0

    @property
    def is_varner_form(self) -> bool:
        """Check if using Varner parameterization."""
        return self.values[22] != 0.0

    @property
    def is_koning_form(self) -> bool:
        """Check if using extended Koning parameterization."""
        return self.values[23] != 0.0

    @property
    def fermi_energy(self) -> Optional[float]:
        """Fermi energy if explicitly specified [MeV]."""
        ef = self.values[17]
        return ef if ef != 0.0 else None


@_dataclass
class EnergyRange:
    """Potential parameters for a specific energy range.

    Each potential component can have multiple energy ranges, each with
    its own set of radius, diffuseness, and strength coefficients.

    Attributes:
        epot: Upper energy limit for this range [MeV]
        radius: Radius coefficients (13 values)
        diffuseness: Diffuseness coefficients (13 values)
        strength: Potential strength coefficients (25 values)
    """
    epot: float = 0.0
    radius: RadiusCoefficients = _field(default_factory=RadiusCoefficients)
    diffuseness: DiffusenessCoefficients = _field(default_factory=DiffusenessCoefficients)
    strength: PotentialStrength = _field(default_factory=PotentialStrength)

    _field_info: _ClassVar[dict] = {
        'epot': 'Upper energy limit for this range [MeV]',
        'radius': 'Radius coefficients',
        'diffuseness': 'Diffuseness coefficients',
        'strength': 'Potential strength coefficients',
    }

    def __repr__(self) -> str:
        return f"EnergyRange(epot={self.epot:.1f}, r0={self.radius.r0:.4f}, a0={self.diffuseness.a0:.4f}, V0={self.strength.V0:.2f})"


@_dataclass
class PotentialComponent:
    """A single component of the optical model potential.

    The OMP has 6 components:
        1. Real volume (Woods-Saxon)
        2. Imaginary volume (Woods-Saxon)
        3. Real surface derivative
        4. Imaginary surface derivative
        5. Real spin-orbit
        6. Imaginary spin-orbit

    Each component can have multiple energy ranges or be absent (jrange=0).

    Attributes:
        component_type: Component index (0-5)
        jrange: Number of energy ranges (0 if component not used)
        use_volume_integrals: If True, pot values are volume integrals not depths
        energy_ranges: List of energy range specifications
    """
    component_type: int = 0
    jrange: int = 0
    use_volume_integrals: bool = False
    energy_ranges: List[EnergyRange] = _field(default_factory=list)

    _field_info: _ClassVar[dict] = {
        'component_type': 'Component index (0=real vol, 1=imag vol, 2=real surf, 3=imag surf, 4=real SO, 5=imag SO)',
        'jrange': 'Number of energy ranges (0 if unused, negative for volume integrals)',
        'use_volume_integrals': 'If True, strength values are volume integrals',
        'energy_ranges': 'List of energy range specifications',
    }

    @property
    def name(self) -> str:
        """Get component name."""
        return config.COMPONENT_NAMES.get(self.component_type, f'component_{self.component_type}')

    @property
    def description(self) -> str:
        """Get component description."""
        return config.COMPONENT_DESCRIPTIONS.get(self.component_type, '')

    @property
    def is_used(self) -> bool:
        """Check if this component is used in the potential."""
        return self.jrange != 0

    def __repr__(self) -> str:
        if not self.is_used:
            return f"PotentialComponent({self.name}, unused)"
        return f"PotentialComponent({self.name}, {abs(self.jrange)} ranges)"


@_dataclass
class CoulombParameters:
    """Coulomb potential parameters for a specific energy range.

    Attributes:
        ecoul: Maximum energy for this Coulomb range [MeV]
        rcoul0: Coulomb radius coefficient for A^(-1/3)
        rcoul: Base Coulomb radius coefficient
        rcoul1: Coulomb radius coefficient for A^(-2/3)
        rcoul2: Coulomb radius coefficient for A^(-5/3)
        rcoul3: Coulomb radius coefficient for A
        beta: Nonlocality range [fm]
        acoul: Diffuseness of Woods-Saxon charge distribution [fm]
    """
    ecoul: float = 0.0
    rcoul0: float = 0.0
    rcoul: float = 0.0
    rcoul1: float = 0.0
    rcoul2: float = 0.0
    rcoul3: float = 0.0
    beta: float = 0.0
    acoul: float = 0.0

    _field_info: _ClassVar[dict] = {
        'ecoul': 'Maximum energy for Coulomb range [MeV]',
        'rcoul0': 'Coulomb radius coefficient for A^(-1/3)',
        'rcoul': 'Base Coulomb radius coefficient',
        'rcoul1': 'Coulomb radius coefficient for A^(-2/3)',
        'rcoul2': 'Coulomb radius coefficient for A^(-5/3)',
        'rcoul3': 'Coulomb radius coefficient for A',
        'beta': 'Nonlocality range [fm]',
        'acoul': 'Diffuseness of charge distribution [fm]',
    }

    def coulomb_radius(self, A: int) -> float:
        """Calculate Coulomb radius for given mass number.

        Args:
            A: Target mass number

        Returns:
            Coulomb radius RC [fm]
        """
        A13 = A ** (1.0 / 3.0)
        return (self.rcoul0 * A ** (-1.0 / 3.0) + self.rcoul +
                self.rcoul1 * A ** (-2.0 / 3.0) + self.rcoul2 * A ** (-5.0 / 3.0) +
                self.rcoul3 * A) * A13

    def __repr__(self) -> str:
        return f"CoulombParameters(ecoul={self.ecoul:.1f}, rcoul={self.rcoul:.4f}, beta={self.beta:.4f})"


@_dataclass
class RotationalLevel:
    """A rotational level in the coupled-channel rotational model.

    Attributes:
        ex: Excitation energy [MeV]
        spin: Level spin
        parity: Level parity (+1 or -1)
    """
    ex: float = 0.0
    spin: float = 0.0
    parity: int = 1

    _field_info: _ClassVar[dict] = {
        'ex': 'Excitation energy [MeV]',
        'spin': 'Level spin',
        'parity': 'Level parity (+1 or -1)',
    }

    def __repr__(self) -> str:
        par_str = '+' if self.parity > 0 else '-'
        return f"RotationalLevel(Ex={self.ex:.4f}, J={self.spin}{par_str})"


@_dataclass
class VibrationalLevel:
    """A vibrational level in the coupled-channel vibrational model.

    Attributes:
        ex: Excitation energy [MeV]
        spin: Level spin
        parity: Level parity (+1 or -1)
        nph: Phonon number (1=one-phonon, 2=two-phonon, 3=mixed)
        defv: Vibrational deformation parameter
        thetm: Mixing parameter [degrees] (for nph=3)
    """
    ex: float = 0.0
    spin: float = 0.0
    parity: int = 1
    nph: int = 1
    defv: float = 0.0
    thetm: float = 0.0

    _field_info: _ClassVar[dict] = {
        'ex': 'Excitation energy [MeV]',
        'spin': 'Level spin',
        'parity': 'Level parity (+1 or -1)',
        'nph': 'Phonon number (1=1-phonon, 2=2-phonon, 3=mixed)',
        'defv': 'Vibrational deformation parameter',
        'thetm': 'Mixing parameter [degrees]',
    }

    def __repr__(self) -> str:
        par_str = '+' if self.parity > 0 else '-'
        return f"VibrationalLevel(Ex={self.ex:.4f}, J={self.spin}{par_str}, nph={self.nph})"


@_dataclass
class IsotopeData:
    """Target isotope data for coupled-channel calculations.

    This holds the deformation parameters and collective levels for
    a specific target isotope used in coupled-channel calculations.

    Attributes:
        Z: Atomic number
        A: Mass number
        ncoll: Number of collective states
        lmax: Maximum multipole order
        idef: Largest order of deformation
        bandk: K quantum number for rotational band
        deformations: Deformation parameters (beta2, beta4, beta6, ...)
        levels: List of collective levels
    """
    Z: int = 0
    A: int = 0
    ncoll: int = 0
    lmax: int = 0
    idef: int = 0
    bandk: float = 0.0
    deformations: List[float] = _field(default_factory=list)
    levels: List = _field(default_factory=list)  # Can be RotationalLevel or VibrationalLevel

    _field_info: _ClassVar[dict] = {
        'Z': 'Atomic number',
        'A': 'Mass number',
        'ncoll': 'Number of collective states',
        'lmax': 'Maximum multipole order',
        'idef': 'Largest order of deformation',
        'bandk': 'K quantum number for rotational band',
        'deformations': 'Deformation parameters [beta2, beta4, ...]',
        'levels': 'List of collective levels',
    }

    def __repr__(self) -> str:
        return f"IsotopeData(Z={self.Z}, A={self.A}, ncoll={self.ncoll}, levels={len(self.levels)})"
