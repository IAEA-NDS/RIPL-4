# -*- coding: utf-8 -*-
"""Parser for spherical optical model potentials (imodel=0).

Spherical potentials are the most common type in the RIPL library (~461 of 581).
This module handles parsing of the full potential specification including
header, validity ranges, 6 potential components, and Coulomb parameters.
"""

from dataclasses import dataclass as _dataclass, field as _field
from typing import List, Optional, Iterator, ClassVar as _ClassVar

from . import config
from .core import (
    RadiusCoefficients, DiffusenessCoefficients, PotentialStrength,
    EnergyRange, PotentialComponent, CoulombParameters
)
from .reader import (
    parse_fortran_float, read_floats_from_line, read_ints_from_line,
    read_text_line, read_text_lines, read_potential_block
)


@_dataclass
class OMPHeader:
    """Header metadata for an optical model potential.

    Attributes:
        iref: Unique library reference number
        author: Author names (up to 80 characters)
        reference: Publication reference (up to 80 characters)
        summary: Description of the potential (4 lines, 320 characters total)
    """
    iref: int = 0
    author: str = ''
    reference: str = ''
    summary: str = ''

    _field_info: _ClassVar[dict] = {
        'iref': 'Unique library reference number',
        'author': 'Author names',
        'reference': 'Publication reference',
        'summary': 'Description of the potential',
    }

    def __repr__(self) -> str:
        return f"OMPHeader(iref={self.iref}, author='{self.author[:30]}...')"


@_dataclass
class ValidityRanges:
    """Validity ranges for an optical model potential.

    Attributes:
        E_min: Minimum projectile energy [MeV]
        E_max: Maximum projectile energy [MeV]
        Z_min: Minimum target atomic number
        Z_max: Maximum target atomic number
        A_min: Minimum target mass number
        A_max: Maximum target mass number
    """
    E_min: float = 0.0
    E_max: float = 0.0
    Z_min: int = 0
    Z_max: int = 0
    A_min: int = 0
    A_max: int = 0

    _field_info: _ClassVar[dict] = {
        'E_min': 'Minimum projectile energy [MeV]',
        'E_max': 'Maximum projectile energy [MeV]',
        'Z_min': 'Minimum target atomic number',
        'Z_max': 'Maximum target atomic number',
        'A_min': 'Minimum target mass number',
        'A_max': 'Maximum target mass number',
    }

    def is_valid_for(self, E: float, Z: int, A: int) -> bool:
        """Check if potential is valid for given conditions.

        Args:
            E: Projectile energy [MeV]
            Z: Target atomic number
            A: Target mass number

        Returns:
            True if potential is valid
        """
        return (self.E_min <= E <= self.E_max and
                self.Z_min <= Z <= self.Z_max and
                self.A_min <= A <= self.A_max)


@_dataclass
class ModelFlags:
    """Model and physics flags for an optical model potential.

    Attributes:
        imodel: Model type (0=spherical, 1-5=coupled-channel)
        iz_proj: Projectile atomic number (0=neutron, 1=proton/d/t, 2=He3/alpha)
        ia_proj: Projectile mass number
        irel: Relativistic flag (0=non-rel, 1=rel kinematics, 2=rel+scaled)
        idr: Dispersion relation flag (0=none, 1-3=various types)
    """
    imodel: int = 0
    iz_proj: int = 0
    ia_proj: int = 1
    irel: int = 0
    idr: int = 0

    _field_info: _ClassVar[dict] = {
        'imodel': 'Model type (0=spherical, 1-5=coupled-channel)',
        'iz_proj': 'Projectile atomic number',
        'ia_proj': 'Projectile mass number',
        'irel': 'Relativistic flag (0=non-rel, 1=rel, 2=rel+scaled)',
        'idr': 'Dispersion relation flag',
    }

    @property
    def projectile_name(self) -> str:
        """Get projectile name."""
        return config.projectile_name(self.iz_proj, self.ia_proj)

    @property
    def model_name(self) -> str:
        """Get model type name."""
        return config.model_name(self.imodel)

    @property
    def is_spherical(self) -> bool:
        """Check if spherical model."""
        return self.imodel == config.MODEL_SPHERICAL

    @property
    def is_relativistic(self) -> bool:
        """Check if relativistic."""
        return config.is_relativistic(self.irel)

    @property
    def uses_dispersion(self) -> bool:
        """Check if uses dispersion relations."""
        return config.uses_dispersion(self.idr)


@_dataclass
class SphericalOMP:
    """A complete spherical optical model potential.

    This is the main data structure for spherical potentials (imodel=0).
    It contains all the information needed to calculate the optical potential
    for a given target nucleus and projectile energy.

    Attributes:
        header: Metadata (iref, author, reference, summary)
        validity: Energy, Z, and A validity ranges
        flags: Model and physics flags
        components: List of 6 potential components (real/imag × vol/surf/SO)
        coulomb: List of Coulomb parameter sets (one per energy range)
    """
    header: OMPHeader = _field(default_factory=OMPHeader)
    validity: ValidityRanges = _field(default_factory=ValidityRanges)
    flags: ModelFlags = _field(default_factory=ModelFlags)
    components: List[PotentialComponent] = _field(default_factory=list)
    coulomb: List[CoulombParameters] = _field(default_factory=list)

    _field_info: _ClassVar[dict] = {
        'header': 'Metadata (iref, author, reference, summary)',
        'validity': 'Energy, Z, and A validity ranges',
        'flags': 'Model and physics flags',
        'components': 'List of 6 potential components',
        'coulomb': 'List of Coulomb parameter sets',
    }

    @property
    def iref(self) -> int:
        """Get library reference number."""
        return self.header.iref

    @property
    def projectile(self) -> str:
        """Get projectile name."""
        return self.flags.projectile_name

    def get_component(self, comp_type: int) -> Optional[PotentialComponent]:
        """Get a specific potential component.

        Args:
            comp_type: Component type (0-5, see config.COMP_*)

        Returns:
            PotentialComponent or None if not found
        """
        for comp in self.components:
            if comp.component_type == comp_type:
                return comp
        return None

    @property
    def real_volume(self) -> Optional[PotentialComponent]:
        """Get real volume component."""
        return self.get_component(config.COMP_REAL_VOLUME)

    @property
    def imag_volume(self) -> Optional[PotentialComponent]:
        """Get imaginary volume component."""
        return self.get_component(config.COMP_IMAG_VOLUME)

    @property
    def real_surface(self) -> Optional[PotentialComponent]:
        """Get real surface component."""
        return self.get_component(config.COMP_REAL_SURFACE)

    @property
    def imag_surface(self) -> Optional[PotentialComponent]:
        """Get imaginary surface component."""
        return self.get_component(config.COMP_IMAG_SURFACE)

    @property
    def real_spinorbit(self) -> Optional[PotentialComponent]:
        """Get real spin-orbit component."""
        return self.get_component(config.COMP_REAL_SPINORBIT)

    @property
    def imag_spinorbit(self) -> Optional[PotentialComponent]:
        """Get imaginary spin-orbit component."""
        return self.get_component(config.COMP_IMAG_SPINORBIT)

    def __repr__(self) -> str:
        n_active = sum(1 for c in self.components if c.is_used)
        return (f"SphericalOMP(iref={self.iref}, {self.projectile}, "
                f"Z={self.validity.Z_min}-{self.validity.Z_max}, "
                f"E={self.validity.E_min:.1f}-{self.validity.E_max:.1f}, "
                f"{n_active} active components)")


def parse_header(lines_iter: Iterator[str]) -> OMPHeader:
    """Parse the header section of an OMP entry.

    Args:
        lines_iter: Iterator over entry lines

    Returns:
        OMPHeader with parsed data
    """
    # Line 1: iref
    iref_line = read_text_line(lines_iter)
    iref = int(iref_line.strip())

    # Line 2: author
    author = read_text_line(lines_iter)

    # Line 3: reference
    reference = read_text_line(lines_iter)

    # Lines 4-7: summary (4 lines)
    summary_lines = read_text_lines(lines_iter, 4)
    summary = '\n'.join(summary_lines)

    return OMPHeader(iref=iref, author=author, reference=reference, summary=summary)


def parse_validity_ranges(lines_iter: Iterator[str]) -> ValidityRanges:
    """Parse the validity range section.

    Args:
        lines_iter: Iterator over entry lines

    Returns:
        ValidityRanges with parsed data
    """
    # Line 1: emin, emax
    e_values = read_floats_from_line(next(lines_iter))
    E_min, E_max = e_values[0], e_values[1]

    # Line 2: izmin, izmax
    z_values = read_ints_from_line(next(lines_iter))
    Z_min, Z_max = z_values[0], z_values[1]

    # Line 3: iamin, iamax
    a_values = read_ints_from_line(next(lines_iter))
    A_min, A_max = a_values[0], a_values[1]

    return ValidityRanges(
        E_min=E_min, E_max=E_max,
        Z_min=Z_min, Z_max=Z_max,
        A_min=A_min, A_max=A_max
    )


def parse_model_flags(lines_iter: Iterator[str]) -> ModelFlags:
    """Parse the model flags line.

    Args:
        lines_iter: Iterator over entry lines

    Returns:
        ModelFlags with parsed data
    """
    # Line: imodel, izproj, iaproj, irel, idr
    values = read_ints_from_line(next(lines_iter))

    return ModelFlags(
        imodel=values[0],
        iz_proj=values[1],
        ia_proj=values[2],
        irel=values[3],
        idr=values[4]
    )


def parse_potential_component(lines_iter: Iterator[str], comp_type: int) -> PotentialComponent:
    """Parse a single potential component.

    Args:
        lines_iter: Iterator over entry lines
        comp_type: Component type index (0-5)

    Returns:
        PotentialComponent with parsed data
    """
    # First line: jrange (number of energy ranges)
    jrange_line = read_text_line(lines_iter)
    jrange = int(jrange_line.strip())

    # Check for volume integrals (negative jrange)
    use_volume_integrals = jrange < 0
    n_ranges = abs(jrange)

    # Read energy ranges
    energy_ranges = []
    for _ in range(n_ranges):
        epot, rco, aco, pot = read_potential_block(lines_iter)

        energy_range = EnergyRange(
            epot=epot,
            radius=RadiusCoefficients(values=rco),
            diffuseness=DiffusenessCoefficients(values=aco),
            strength=PotentialStrength(values=pot)
        )
        energy_ranges.append(energy_range)

    return PotentialComponent(
        component_type=comp_type,
        jrange=jrange,
        use_volume_integrals=use_volume_integrals,
        energy_ranges=energy_ranges
    )


def parse_coulomb_parameters(lines_iter: Iterator[str]) -> List[CoulombParameters]:
    """Parse the Coulomb parameters section.

    Args:
        lines_iter: Iterator over entry lines

    Returns:
        List of CoulombParameters (one per energy range)
    """
    # First line: jcoul (number of Coulomb energy ranges)
    jcoul_line = read_text_line(lines_iter)
    jcoul = int(jcoul_line.strip())

    coulomb_list = []
    for _ in range(jcoul):
        # Each line: ecoul, rcoul0, rcoul, rcoul1, rcoul2, beta, acoul, rcoul3
        values = read_floats_from_line(next(lines_iter))

        # Ensure we have enough values (pad with zeros if needed)
        while len(values) < 8:
            values.append(0.0)

        coulomb = CoulombParameters(
            ecoul=values[0],
            rcoul0=values[1],
            rcoul=values[2],
            rcoul1=values[3],
            rcoul2=values[4],
            beta=values[5],
            acoul=values[6],
            rcoul3=values[7] if len(values) > 7 else 0.0
        )
        coulomb_list.append(coulomb)

    return coulomb_list


def parse_spherical_omp(lines: List[str]) -> SphericalOMP:
    """Parse a complete spherical OMP entry.

    Args:
        lines: List of lines for one OMP entry (between separators)

    Returns:
        SphericalOMP with all parsed data

    Raises:
        ValueError: If entry is not a spherical potential
    """
    lines_iter = iter(lines)

    # Parse header
    header = parse_header(lines_iter)

    # Parse validity ranges
    validity = parse_validity_ranges(lines_iter)

    # Parse model flags
    flags = parse_model_flags(lines_iter)

    # Check that this is a spherical potential
    if not flags.is_spherical:
        raise ValueError(f"Entry {header.iref} is not spherical (imodel={flags.imodel})")

    # Parse 6 potential components
    components = []
    for comp_type in range(config.NUM_COMPONENTS):
        comp = parse_potential_component(lines_iter, comp_type)
        components.append(comp)

    # Parse Coulomb parameters
    coulomb = parse_coulomb_parameters(lines_iter)

    return SphericalOMP(
        header=header,
        validity=validity,
        flags=flags,
        components=components,
        coulomb=coulomb
    )


def read_spherical_omp_by_iref(filepath: str, target_iref: int) -> Optional[SphericalOMP]:
    """Read a specific spherical OMP from the file by reference number.

    Args:
        filepath: Path to om-parameter-u.dat
        target_iref: Library reference number to find

    Returns:
        SphericalOMP if found and is spherical, None otherwise
    """
    from .reader import OMPFileReader

    with OMPFileReader(filepath) as reader:
        for entry_lines in reader.entries():
            if not entry_lines:
                continue

            # Quick check: first line should be the iref
            try:
                iref = int(entry_lines[0].strip())
                if iref == target_iref:
                    try:
                        return parse_spherical_omp(entry_lines)
                    except ValueError:
                        # Not a spherical potential
                        return None
            except ValueError:
                continue

    return None
