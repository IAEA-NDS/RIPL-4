# -*- coding: utf-8 -*-
"""Parser for coupled-channel optical model potentials (imodel=1-5).

Coupled-channel potentials extend spherical potentials with additional
target isotope data including deformation parameters and collective levels.

Model types:
    - imodel=1: Rigid rotor (rotational model)
    - imodel=2: Vibrational (harmonic oscillator)
    - imodel=3: Soft rotor (vibrational-rotational)
    - imodel=4: Rigid-soft rotor (rigid multiband)
    - imodel=5: Soft deformed rotor (soft multiband)
"""

from dataclasses import dataclass as _dataclass, field as _field
from typing import List, Optional, Iterator, Union, ClassVar as _ClassVar

from . import config
from .core import (
    RadiusCoefficients, DiffusenessCoefficients, PotentialStrength,
    EnergyRange, PotentialComponent, CoulombParameters,
    RotationalLevel, VibrationalLevel, IsotopeData
)
from .spherical import (
    OMPHeader, ValidityRanges, ModelFlags,
    parse_header, parse_validity_ranges, parse_model_flags,
    parse_potential_component, parse_coulomb_parameters
)
from .reader import (
    parse_fortran_float, read_floats_from_line, read_ints_from_line,
    read_text_line
)


# =============================================================================
# Soft Rotor Hamiltonian Parameters
# =============================================================================

@_dataclass
class SoftRotorHamiltonian:
    """Hamiltonian parameters for soft rotor models (imodel=3, 5).

    These parameters define the nuclear collective Hamiltonian used in
    soft rotor coupled-channel calculations.

    Attributes:
        hw: Phonon energy [MeV]
        amb0: Mass parameter for beta vibration
        amg0: Mass parameter for gamma vibration
        gam0: Ground state gamma deformation [degrees]
        bet0: Ground state beta deformation
        bet4: Hexadecapole deformation
        bb42: Beta4-beta2 coupling
        gamg: Gamma band gamma deformation
        delg: Gamma band energy shift
        bet3: Octupole deformation
        et0: Reference energy
        amu0: Reduced mass parameter
        hw0: Zero-point energy
        bb32: Beta3-beta2 coupling
        gamde: Gamma deformation energy
        dpar: Parity splitting parameter
        gshape: Shape parameter
    """
    hw: float = 0.0
    amb0: float = 0.0
    amg0: float = 0.0
    gam0: float = 0.0
    bet0: float = 0.0
    bet4: float = 0.0
    bb42: float = 0.0
    gamg: float = 0.0
    delg: float = 0.0
    bet3: float = 0.0
    et0: float = 0.0
    amu0: float = 0.0
    hw0: float = 0.0
    bb32: float = 0.0
    gamde: float = 0.0
    dpar: float = 0.0
    gshape: float = 0.0

    _field_info: _ClassVar[dict] = {
        'hw': 'Phonon energy [MeV]',
        'amb0': 'Mass parameter for beta vibration',
        'amg0': 'Mass parameter for gamma vibration',
        'gam0': 'Ground state gamma deformation [degrees]',
        'bet0': 'Ground state beta deformation',
        'bet4': 'Hexadecapole deformation',
    }


@_dataclass
class SoftRotorLevel:
    """A level in the soft rotor model with quantum numbers.

    Attributes:
        ex: Excitation energy [MeV]
        spin: Level spin
        parity: Level parity (+1 or -1)
        ntu: Phonon quantum number (ntu)
        nnb: Beta phonon quantum number
        nng: Gamma phonon quantum number
        nno: Octupole phonon quantum number
    """
    ex: float = 0.0
    spin: float = 0.0
    parity: int = 1
    ntu: int = 0
    nnb: int = 0
    nng: int = 0
    nno: int = 0

    _field_info: _ClassVar[dict] = {
        'ex': 'Excitation energy [MeV]',
        'spin': 'Level spin',
        'parity': 'Level parity (+1 or -1)',
        'ntu': 'Phonon quantum number',
        'nnb': 'Beta phonon quantum number',
        'nng': 'Gamma phonon quantum number',
        'nno': 'Octupole phonon quantum number',
    }

    def __repr__(self) -> str:
        par_str = '+' if self.parity > 0 else '-'
        return f"SoftRotorLevel(Ex={self.ex:.4f}, J={self.spin}{par_str})"


@_dataclass
class RigidSoftLevel:
    """A level in the rigid-soft model (imodel=4, 5).

    Attributes:
        ex: Excitation energy [MeV]
        spin: Level spin
        parity: Level parity (+1 or -1)
        ko: K of the band (multiplied by 2)
        nca: Flag for isobar-analogue states
        numb: Ordinal number of the band
        defv: Effective deformation
        defr: Single-particle states overlapping
    """
    ex: float = 0.0
    spin: float = 0.0
    parity: int = 1
    ko: int = 0
    nca: int = 0
    numb: int = 0
    defv: float = 0.0
    defr: float = 0.0
    # For imodel=5, also include soft rotor quantum numbers
    ntu: int = 0
    nnb: int = 0
    nng: int = 0
    nno: int = 0

    _field_info: _ClassVar[dict] = {
        'ex': 'Excitation energy [MeV]',
        'spin': 'Level spin',
        'parity': 'Level parity (+1 or -1)',
        'ko': 'K of the band (multiplied by 2)',
        'nca': 'Flag for isobar-analogue states',
        'numb': 'Ordinal number of the band',
        'defv': 'Effective deformation',
        'defr': 'Single-particle states overlapping',
    }

    def __repr__(self) -> str:
        par_str = '+' if self.parity > 0 else '-'
        return f"RigidSoftLevel(Ex={self.ex:.4f}, J={self.spin}{par_str}, K={self.ko/2})"


# =============================================================================
# Isotope Data Extensions
# =============================================================================

@_dataclass
class RotationalIsotope(IsotopeData):
    """Isotope data for rigid rotor model (imodel=1).

    Extends IsotopeData with rotational band information.
    """
    levels: List[RotationalLevel] = _field(default_factory=list)


@_dataclass
class VibrationalIsotope(IsotopeData):
    """Isotope data for vibrational model (imodel=2)."""
    levels: List[VibrationalLevel] = _field(default_factory=list)


@_dataclass
class SoftRotorIsotope(IsotopeData):
    """Isotope data for soft rotor model (imodel=3).

    Includes Hamiltonian parameters and soft rotor levels.
    """
    hamiltonian: SoftRotorHamiltonian = _field(default_factory=SoftRotorHamiltonian)
    levels: List[SoftRotorLevel] = _field(default_factory=list)


@_dataclass
class RigidSoftIsotope(IsotopeData):
    """Isotope data for rigid-soft model (imodel=4).

    Includes deformation parameters and rigid-soft levels.
    """
    levels: List[RigidSoftLevel] = _field(default_factory=list)


@_dataclass
class SoftDeformedIsotope(IsotopeData):
    """Isotope data for soft deformed model (imodel=5).

    Includes Hamiltonian parameters, deformations, and extended levels.
    """
    hamiltonian: SoftRotorHamiltonian = _field(default_factory=SoftRotorHamiltonian)
    levels: List[RigidSoftLevel] = _field(default_factory=list)


# =============================================================================
# Coupled-Channel OMP
# =============================================================================

@_dataclass
class CoupledChannelOMP:
    """A complete coupled-channel optical model potential.

    This extends the spherical OMP with target isotope data for
    coupled-channel calculations.

    Attributes:
        header: Metadata (iref, author, reference, summary)
        validity: Energy, Z, and A validity ranges
        flags: Model and physics flags (imodel > 0)
        components: List of 6 potential components
        coulomb: List of Coulomb parameter sets
        isotopes: List of target isotope data with collective structure
    """
    header: OMPHeader = _field(default_factory=OMPHeader)
    validity: ValidityRanges = _field(default_factory=ValidityRanges)
    flags: ModelFlags = _field(default_factory=ModelFlags)
    components: List[PotentialComponent] = _field(default_factory=list)
    coulomb: List[CoulombParameters] = _field(default_factory=list)
    isotopes: List[IsotopeData] = _field(default_factory=list)

    _field_info: _ClassVar[dict] = {
        'header': 'Metadata (iref, author, reference, summary)',
        'validity': 'Energy, Z, and A validity ranges',
        'flags': 'Model and physics flags',
        'components': 'List of 6 potential components',
        'coulomb': 'List of Coulomb parameter sets',
        'isotopes': 'List of target isotope data',
    }

    @property
    def iref(self) -> int:
        """Get library reference number."""
        return self.header.iref

    @property
    def projectile(self) -> str:
        """Get projectile name."""
        return self.flags.projectile_name

    @property
    def model_name(self) -> str:
        """Get model type name."""
        return self.flags.model_name

    @property
    def n_isotopes(self) -> int:
        """Get number of target isotopes."""
        return len(self.isotopes)

    def get_component(self, comp_type: int) -> Optional[PotentialComponent]:
        """Get a specific potential component."""
        for comp in self.components:
            if comp.component_type == comp_type:
                return comp
        return None

    def get_isotope(self, Z: int, A: int) -> Optional[IsotopeData]:
        """Get isotope data for specific target.

        Args:
            Z: Target atomic number
            A: Target mass number

        Returns:
            IsotopeData if found, None otherwise
        """
        for iso in self.isotopes:
            if iso.Z == Z and iso.A == A:
                return iso
        return None

    def __repr__(self) -> str:
        n_active = sum(1 for c in self.components if c.is_used)
        return (f"CoupledChannelOMP(iref={self.iref}, {self.projectile}, "
                f"{self.model_name}, "
                f"Z={self.validity.Z_min}-{self.validity.Z_max}, "
                f"E={self.validity.E_min:.1f}-{self.validity.E_max:.1f}, "
                f"{n_active} components, {self.n_isotopes} isotopes)")


# =============================================================================
# Parsing Functions
# =============================================================================

def parse_rigid_rotor_isotopes(lines_iter: Iterator[str]) -> List[RotationalIsotope]:
    """Parse isotope data for rigid rotor model (imodel=1).

    Format:
        nisotopes
        iz(n), ia(n), ncoll(n), lmax(n), idef(n), bandk(n), [def(j), j=2,idef,2]
        ex(k), spin(k), ipar(k)  [for k=1,ncoll]
    """
    # Read number of isotopes
    nisotopes_line = read_text_line(lines_iter)
    nisotopes = int(nisotopes_line.strip())

    isotopes = []
    for _ in range(nisotopes):
        # Read isotope header line
        header_values = read_floats_from_line(next(lines_iter))

        # Parse fixed values
        Z = int(header_values[0])
        A = int(header_values[1])
        ncoll = int(header_values[2])
        lmax = int(header_values[3])
        idef = int(header_values[4])
        bandk = header_values[5]

        # Parse deformation parameters (beta2, beta4, beta6, ...)
        # These come after bandk, starting at index 6
        deformations = []
        for i in range(6, len(header_values)):
            deformations.append(header_values[i])

        # Read collective levels
        levels = []
        for _ in range(ncoll):
            level_values = read_floats_from_line(next(lines_iter))
            level = RotationalLevel(
                ex=level_values[0],
                spin=level_values[1],
                parity=int(level_values[2])
            )
            levels.append(level)

        isotope = RotationalIsotope(
            Z=Z, A=A, ncoll=ncoll, lmax=lmax, idef=idef,
            bandk=bandk, deformations=deformations, levels=levels
        )
        isotopes.append(isotope)

    return isotopes


def parse_vibrational_isotopes(lines_iter: Iterator[str]) -> List[VibrationalIsotope]:
    """Parse isotope data for vibrational model (imodel=2).

    Format:
        nisotopes
        iz(n), ia(n), nvib(n)
        exv(k), spinv(k), iparv(k), nph(k), defv(k), thetm(k)  [for k=1,nvib]
    """
    nisotopes_line = read_text_line(lines_iter)
    nisotopes = int(nisotopes_line.strip())

    isotopes = []
    for _ in range(nisotopes):
        # Read isotope header
        header_values = read_ints_from_line(next(lines_iter))
        Z = header_values[0]
        A = header_values[1]
        nvib = header_values[2]

        # Read vibrational levels
        levels = []
        for _ in range(nvib):
            level_values = read_floats_from_line(next(lines_iter))
            level = VibrationalLevel(
                ex=level_values[0],
                spin=level_values[1],
                parity=int(level_values[2]),
                nph=int(level_values[3]),
                defv=level_values[4],
                thetm=level_values[5] if len(level_values) > 5 else 0.0
            )
            levels.append(level)

        isotope = VibrationalIsotope(
            Z=Z, A=A, ncoll=nvib, levels=levels
        )
        isotopes.append(isotope)

    return isotopes


def parse_soft_rotor_isotopes(lines_iter: Iterator[str]) -> List[SoftRotorIsotope]:
    """Parse isotope data for soft rotor model (imodel=3).

    Format:
        nisotopes
        iz(n), ia(n), ncoll(n)
        SR_hw, SR_amb0, SR_amg0, SR_gam0, SR_bet0, SR_bet4  (Hamiltonian line 1)
        SR_bb42, SR_gamg, SR_delg, SR_bet3, SR_et0, SR_amu0 (Hamiltonian line 2)
        SR_hw0, SR_bb32, SR_gamde, SR_dpar, SR_gshape       (Hamiltonian line 3)
        exv(k), spinv(k), iparv(k), SR_ntu, SR_nnb, SR_nng, SR_nno  [for k=1,ncoll]
    """
    nisotopes_line = read_text_line(lines_iter)
    nisotopes = int(nisotopes_line.strip())

    isotopes = []
    for _ in range(nisotopes):
        # Read isotope header
        header_values = read_ints_from_line(next(lines_iter))
        Z = header_values[0]
        A = header_values[1]
        ncoll = header_values[2]

        # Read Hamiltonian parameters (3 lines)
        ham1 = read_floats_from_line(next(lines_iter))
        ham2 = read_floats_from_line(next(lines_iter))
        ham3 = read_floats_from_line(next(lines_iter))

        hamiltonian = SoftRotorHamiltonian(
            hw=ham1[0], amb0=ham1[1], amg0=ham1[2],
            gam0=ham1[3], bet0=ham1[4], bet4=ham1[5] if len(ham1) > 5 else 0.0,
            bb42=ham2[0], gamg=ham2[1], delg=ham2[2],
            bet3=ham2[3], et0=ham2[4], amu0=ham2[5] if len(ham2) > 5 else 0.0,
            hw0=ham3[0], bb32=ham3[1], gamde=ham3[2],
            dpar=ham3[3] if len(ham3) > 3 else 0.0,
            gshape=ham3[4] if len(ham3) > 4 else 0.0
        )

        # Read soft rotor levels
        levels = []
        for _ in range(ncoll):
            level_values = read_floats_from_line(next(lines_iter))
            level = SoftRotorLevel(
                ex=level_values[0],
                spin=level_values[1],
                parity=int(level_values[2]),
                ntu=int(level_values[3]),
                nnb=int(level_values[4]),
                nng=int(level_values[5]),
                nno=int(level_values[6]) if len(level_values) > 6 else 0
            )
            levels.append(level)

        isotope = SoftRotorIsotope(
            Z=Z, A=A, ncoll=ncoll,
            hamiltonian=hamiltonian, levels=levels
        )
        isotopes.append(isotope)

    return isotopes


def parse_rigid_soft_isotopes(lines_iter: Iterator[str]) -> List[RigidSoftIsotope]:
    """Parse isotope data for rigid-soft model (imodel=4).

    Format:
        nisotopes
        iz(n), ia(n), ncoll(n), lmax(n), idef(n), bandk(n), [def(j), j=2,idef,2]
        exv(k), spinv(k), iparv(k), OPT_ko, OPT_nca, OPT_numb, defv, defr  [for k=1,ncoll]
    """
    nisotopes_line = read_text_line(lines_iter)
    nisotopes = int(nisotopes_line.strip())

    isotopes = []
    for _ in range(nisotopes):
        # Read isotope header (same as rigid rotor)
        header_values = read_floats_from_line(next(lines_iter))
        Z = int(header_values[0])
        A = int(header_values[1])
        ncoll = int(header_values[2])
        lmax = int(header_values[3])
        idef = int(header_values[4])
        bandk = header_values[5]

        # Parse deformations
        deformations = []
        for i in range(6, len(header_values)):
            deformations.append(header_values[i])

        # Read rigid-soft levels
        levels = []
        for _ in range(ncoll):
            level_values = read_floats_from_line(next(lines_iter))
            level = RigidSoftLevel(
                ex=level_values[0],
                spin=level_values[1],
                parity=int(level_values[2]),
                ko=int(level_values[3]),
                nca=int(level_values[4]),
                numb=int(level_values[5]),
                defv=level_values[6] if len(level_values) > 6 else 0.0,
                defr=level_values[7] if len(level_values) > 7 else 0.0
            )
            levels.append(level)

        isotope = RigidSoftIsotope(
            Z=Z, A=A, ncoll=ncoll, lmax=lmax, idef=idef,
            bandk=bandk, deformations=deformations, levels=levels
        )
        isotopes.append(isotope)

    return isotopes


def parse_soft_deformed_isotopes(lines_iter: Iterator[str]) -> List[SoftDeformedIsotope]:
    """Parse isotope data for soft deformed model (imodel=5).

    Format:
        nisotopes
        iz(n), ia(n), ncoll(n), lmax(n), idef(n), bandk(n), [def(j), j=2,idef,2]
        SR_hw, SR_amb0, SR_amg0, SR_gam0, SR_bet0, SR_bet4  (Hamiltonian line 1)
        SR_bb42, SR_gamg, SR_delg, SR_bet3, SR_et0, SR_amu0 (Hamiltonian line 2)
        SR_hw0, SR_bb32, SR_gamde, SR_dpar, SR_gshape       (Hamiltonian line 3)
        exv, spinv, iparv, OPT_ko, OPT_nca, OPT_numb, defv, defr, SR_ntu, SR_nnb, SR_nng, SR_nno

    Note:
        Some entries in RIPL have truncated data (nisotopes > actual isotopes).
        This function handles such cases gracefully by returning partial data.
    """
    nisotopes_line = read_text_line(lines_iter)
    nisotopes = int(nisotopes_line.strip())

    isotopes = []
    for _ in range(nisotopes):
        try:
            # Read isotope header (same as rigid rotor)
            header_values = read_floats_from_line(next(lines_iter))
            Z = int(header_values[0])
            A = int(header_values[1])
            ncoll = int(header_values[2])
            lmax = int(header_values[3])
            idef = int(header_values[4])
            bandk = header_values[5]

            # Parse deformations
            deformations = []
            for i in range(6, len(header_values)):
                deformations.append(header_values[i])

            # Read Hamiltonian parameters (3 lines)
            ham1 = read_floats_from_line(next(lines_iter))
            ham2 = read_floats_from_line(next(lines_iter))
            ham3 = read_floats_from_line(next(lines_iter))

            hamiltonian = SoftRotorHamiltonian(
                hw=ham1[0], amb0=ham1[1], amg0=ham1[2],
                gam0=ham1[3], bet0=ham1[4], bet4=ham1[5] if len(ham1) > 5 else 0.0,
                bb42=ham2[0], gamg=ham2[1], delg=ham2[2],
                bet3=ham2[3], et0=ham2[4], amu0=ham2[5] if len(ham2) > 5 else 0.0,
                hw0=ham3[0], bb32=ham3[1], gamde=ham3[2],
                dpar=ham3[3] if len(ham3) > 3 else 0.0,
                gshape=ham3[4] if len(ham3) > 4 else 0.0
            )

            # Read extended levels
            levels = []
            for _ in range(ncoll):
                level_values = read_floats_from_line(next(lines_iter))
                level = RigidSoftLevel(
                    ex=level_values[0],
                    spin=level_values[1],
                    parity=int(level_values[2]),
                    ko=int(level_values[3]),
                    nca=int(level_values[4]),
                    numb=int(level_values[5]),
                    defv=level_values[6] if len(level_values) > 6 else 0.0,
                    defr=level_values[7] if len(level_values) > 7 else 0.0,
                    ntu=int(level_values[8]) if len(level_values) > 8 else 0,
                    nnb=int(level_values[9]) if len(level_values) > 9 else 0,
                    nng=int(level_values[10]) if len(level_values) > 10 else 0,
                    nno=int(level_values[11]) if len(level_values) > 11 else 0
                )
                levels.append(level)

            isotope = SoftDeformedIsotope(
                Z=Z, A=A, ncoll=ncoll, lmax=lmax, idef=idef,
                bandk=bandk, deformations=deformations,
                hamiltonian=hamiltonian, levels=levels
            )
            isotopes.append(isotope)
        except StopIteration:
            # Truncated data - return what we have
            break

    return isotopes


def parse_coupled_channel_omp(lines: List[str]) -> CoupledChannelOMP:
    """Parse a complete coupled-channel OMP entry.

    Args:
        lines: List of lines for one OMP entry (between separators)

    Returns:
        CoupledChannelOMP with all parsed data

    Raises:
        ValueError: If entry is a spherical potential (imodel=0)
    """
    lines_iter = iter(lines)

    # Parse header (same as spherical)
    header = parse_header(lines_iter)

    # Parse validity ranges (same as spherical)
    validity = parse_validity_ranges(lines_iter)

    # Parse model flags
    flags = parse_model_flags(lines_iter)

    # Check that this is NOT a spherical potential
    if flags.is_spherical:
        raise ValueError(f"Entry {header.iref} is spherical (imodel=0), not coupled-channel")

    # Parse 6 potential components (same as spherical)
    components = []
    for comp_type in range(config.NUM_COMPONENTS):
        comp = parse_potential_component(lines_iter, comp_type)
        components.append(comp)

    # Parse Coulomb parameters (same as spherical)
    coulomb = parse_coulomb_parameters(lines_iter)

    # Parse isotope data based on model type
    imodel = flags.imodel
    if imodel == config.MODEL_RIGID_ROTOR:
        isotopes = parse_rigid_rotor_isotopes(lines_iter)
    elif imodel == config.MODEL_VIBRATIONAL:
        isotopes = parse_vibrational_isotopes(lines_iter)
    elif imodel == config.MODEL_SOFT_ROTOR:
        isotopes = parse_soft_rotor_isotopes(lines_iter)
    elif imodel == config.MODEL_RIGID_SOFT:
        isotopes = parse_rigid_soft_isotopes(lines_iter)
    elif imodel == config.MODEL_SOFT_DEFORMED:
        isotopes = parse_soft_deformed_isotopes(lines_iter)
    else:
        raise ValueError(f"Unknown model type: imodel={imodel}")

    return CoupledChannelOMP(
        header=header,
        validity=validity,
        flags=flags,
        components=components,
        coulomb=coulomb,
        isotopes=isotopes
    )


def read_coupled_channel_omp_by_iref(filepath: str, target_iref: int) -> Optional[CoupledChannelOMP]:
    """Read a specific coupled-channel OMP from the file by reference number.

    Args:
        filepath: Path to om-parameter-u.dat
        target_iref: Library reference number to find

    Returns:
        CoupledChannelOMP if found and is coupled-channel, None otherwise
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
                        return parse_coupled_channel_omp(entry_lines)
                    except ValueError:
                        # Not a coupled-channel potential
                        return None
            except ValueError:
                continue

    return None
