# -*- coding: utf-8 -*-
"""Main optical model potential database.

This module provides the primary interface for loading and accessing
optical model potentials from the RIPL library.
"""

import os as _os
from typing import Dict, List, Optional, Iterator, Union

from riplpy.db import Database
from . import config
from .spherical import SphericalOMP, parse_spherical_omp
from .coupled_channel import CoupledChannelOMP, parse_coupled_channel_omp
from .reader import OMPFileReader

# Type alias for any OMP type
OMPType = Union[SphericalOMP, CoupledChannelOMP]


class OMPDatabase(Database):
    """Database of optical model potentials.

    This class provides access to the full library of optical model potentials,
    supporting lookup by reference number and filtering by various criteria.

    Supports all model types:
        - Spherical potentials (imodel=0): ~461 potentials
        - Rigid rotor (imodel=1): ~99 potentials
        - Vibrational (imodel=2): ~10 potentials
        - Soft rotor (imodel=3): ~7 potentials
        - Rigid-soft (imodel=4): ~5 potentials
        - Soft deformed (imodel=5): ~2 potentials

    Example:
        >>> db = OMPDatabase()
        >>> db.load('/path/to/om-parameter-u.dat')
        >>> pot = db.get(2405)
        >>> print(pot.projectile)
        'n'
    """

    entry = SphericalOMP  # Default entry type

    def __init__(self, data: Dict[int, OMPType] = None):
        """Initialize the OMP database.

        Args:
            data: Optional dictionary mapping iref -> potential
        """
        super().__init__(data)
        self._n_spherical = 0
        self._n_coupled_channel = 0

    def __repr__(self) -> str:
        """Return string representation of the database."""
        return f"<OMPDatabase: {len(self.data)} potentials ({self._n_spherical} spherical, {self._n_coupled_channel} CC)>"

    def load(self, fpath: str) -> None:
        """Load optical model potentials from the parameter file.

        Args:
            fpath: Path to om-parameter-u.dat file

        Loads all potential types (spherical and coupled-channel).
        """
        self.data = {}
        self._n_spherical = 0
        self._n_coupled_channel = 0

        with OMPFileReader(fpath) as reader:
            for entry_lines in reader.entries():
                if not entry_lines:
                    continue

                try:
                    # Get imodel from flags line to determine parser
                    flags_line = entry_lines[10].split()
                    imodel = int(flags_line[0])

                    if imodel == 0:
                        omp = parse_spherical_omp(entry_lines)
                        self._n_spherical += 1
                    else:
                        omp = parse_coupled_channel_omp(entry_lines)
                        self._n_coupled_channel += 1

                    self.data[omp.iref] = omp
                except Exception:
                    # Skip malformed entries
                    pass

    def get(self, iref: int) -> OMPType:
        """Get a potential by its reference number.

        Args:
            iref: Library reference number

        Returns:
            SphericalOMP or CoupledChannelOMP for the specified potential

        Raises:
            EntryNotFoundError: If iref not found in database (a KeyError
                subclass, so ``except KeyError`` still works).
        """
        if iref not in self.data:
            from riplpy.exceptions import EntryNotFoundError
            raise EntryNotFoundError(f"OMP reference {iref} not found in database")
        return self.data[iref]

    def get_by_iref(self, iref: int) -> Optional[OMPType]:
        """Get a potential by reference number, returning None if not found.

        Args:
            iref: Library reference number

        Returns:
            OMP if found, None otherwise
        """
        return self.data.get(iref)

    @property
    def irefs(self) -> List[int]:
        """Get list of all reference numbers in the database."""
        return sorted(self.data.keys())

    @property
    def n_spherical(self) -> int:
        """Get count of spherical potentials."""
        return self._n_spherical

    @property
    def n_coupled_channel(self) -> int:
        """Get count of coupled-channel potentials."""
        return self._n_coupled_channel

    def filter_by_projectile(self, projectile: str) -> "OMPDatabase":
        """Filter potentials by projectile type.

        Args:
            projectile: Projectile name ('n', 'p', 'd', 't', 'He3', 'a')

        Returns:
            New OMPDatabase with filtered potentials
        """
        # Map projectile names to (Z, A) pairs
        proj_map = {
            'n': (0, 1),
            'p': (1, 1),
            'd': (1, 2),
            't': (1, 3),
            'He3': (2, 3),
            '3He': (2, 3),
            'a': (2, 4),
            'alpha': (2, 4),
        }

        if projectile not in proj_map:
            raise ValueError(f"Unknown projectile: {projectile}. "
                           f"Valid options: {list(proj_map.keys())}")

        iz, ia = proj_map[projectile]
        filtered = {
            k: v for k, v in self.data.items()
            if v.flags.iz_proj == iz and v.flags.ia_proj == ia
        }
        result = OMPDatabase(filtered)
        return result

    def filter_by_target(self, Z: int = None, A: int = None) -> "OMPDatabase":
        """Filter potentials by target nucleus.

        Args:
            Z: Target atomic number (filters Z_min <= Z <= Z_max)
            A: Target mass number (filters A_min <= A <= A_max)

        Returns:
            New OMPDatabase with filtered potentials
        """
        filtered = {}
        for k, v in self.data.items():
            if Z is not None:
                if not (v.validity.Z_min <= Z <= v.validity.Z_max):
                    continue
            if A is not None:
                if not (v.validity.A_min <= A <= v.validity.A_max):
                    continue
            filtered[k] = v
        return OMPDatabase(filtered)

    def filter_by_energy(self, E: float) -> "OMPDatabase":
        """Filter potentials by projectile energy.

        Args:
            E: Projectile energy [MeV]

        Returns:
            New OMPDatabase with filtered potentials
        """
        filtered = {}
        for k, v in self.data.items():
            E_min = v.validity.E_min
            E_max = v.validity.E_max
            # Handle None values (treat as unbounded)
            if E_min is None:
                E_min = 0.0
            if E_max is None:
                E_max = float('inf')
            if E_min <= E <= E_max:
                filtered[k] = v
        return OMPDatabase(filtered)

    def filter_relativistic(self) -> "OMPDatabase":
        """Filter to relativistic potentials only.

        Returns:
            New OMPDatabase with only relativistic potentials
        """
        filtered = {
            k: v for k, v in self.data.items()
            if v.flags.is_relativistic
        }
        return OMPDatabase(filtered)

    def filter_dispersive(self) -> "OMPDatabase":
        """Filter to dispersive potentials only.

        Returns:
            New OMPDatabase with only dispersive potentials
        """
        filtered = {
            k: v for k, v in self.data.items()
            if v.flags.uses_dispersion
        }
        return OMPDatabase(filtered)

    def filter_spherical(self) -> "OMPDatabase":
        """Filter to spherical potentials only.

        Returns:
            New OMPDatabase with only spherical potentials
        """
        filtered = {
            k: v for k, v in self.data.items()
            if isinstance(v, SphericalOMP)
        }
        return OMPDatabase(filtered)

    def filter_coupled_channel(self) -> "OMPDatabase":
        """Filter to coupled-channel potentials only.

        Returns:
            New OMPDatabase with only coupled-channel potentials
        """
        filtered = {
            k: v for k, v in self.data.items()
            if isinstance(v, CoupledChannelOMP)
        }
        return OMPDatabase(filtered)

    def filter_by_model(self, imodel: int) -> "OMPDatabase":
        """Filter potentials by model type.

        Args:
            imodel: Model type (0=spherical, 1=rigid rotor, 2=vibrational,
                    3=soft rotor, 4=rigid-soft, 5=soft deformed)

        Returns:
            New OMPDatabase with filtered potentials
        """
        filtered = {
            k: v for k, v in self.data.items()
            if v.flags.imodel == imodel
        }
        return OMPDatabase(filtered)

    def find_for_reaction(self, projectile: str, Z: int, A: int, E: float = None) -> List[OMPType]:
        """Find all potentials valid for a specific reaction.

        Args:
            projectile: Projectile type ('n', 'p', etc.)
            Z: Target atomic number
            A: Target mass number
            E: Projectile energy [MeV]. If None, no energy filtering.

        Returns:
            List of OMP objects valid for this reaction
        """
        results = self.filter_by_projectile(projectile)
        results = results.filter_by_target(Z=Z, A=A)
        if E is not None:
            results = results.filter_by_energy(E)
        return list(results.data.values())

    def info(self) -> str:
        """Return detailed information about this database."""
        lines = ["OMP Database"]
        lines.append("-" * 40)
        lines.append(f"Total potentials: {len(self.data)}")
        lines.append(f"Spherical: {self._n_spherical}")
        lines.append(f"Coupled-channel: {self._n_coupled_channel}")

        # Count by model type
        model_counts = {}
        for v in self.data.values():
            model = v.flags.model_name
            model_counts[model] = model_counts.get(model, 0) + 1
        lines.append("\nBy model type:")
        for model, count in sorted(model_counts.items()):
            lines.append(f"  {model}: {count}")

        # Count by projectile
        proj_counts = {}
        for v in self.data.values():
            proj = v.projectile
            proj_counts[proj] = proj_counts.get(proj, 0) + 1
        lines.append("\nBy projectile:")
        for proj, count in sorted(proj_counts.items()):
            lines.append(f"  {proj}: {count}")

        # Energy range
        if self.data:
            e_min = min(v.validity.E_min for v in self.data.values())
            e_max = max(v.validity.E_max for v in self.data.values())
            lines.append(f"\nEnergy range: {e_min:.1f} - {e_max:.1f} MeV")

        return "\n".join(lines)

    def __iter__(self) -> Iterator[int]:
        """Iterate over reference numbers."""
        return iter(self.data.keys())


def load(directory: str = None, file_path: str = None) -> OMPDatabase:
    """Load the optical model potential database.

    Args:
        directory: RIPL directory path. If provided, constructs full path.
        file_path: Direct path to parameter file. Takes precedence over directory.

    Returns:
        OMPDatabase with loaded potentials

    Raises:
        FileNotFoundError: If the parameter file cannot be found
    """
    if file_path:
        path = file_path
    elif directory:
        path = _os.path.join(directory, config.get_data_file_path('parameters'))
    else:
        # Try default RIPL path
        from riplpy import config as riplpy_config
        ripl_path = riplpy_config.get_path()
        if ripl_path:
            path = _os.path.join(ripl_path, config.get_data_file_path('parameters'))
        else:
            raise FileNotFoundError(
                "Cannot locate OMP parameter file. Please provide directory or file_path."
            )

    if not _os.path.exists(path):
        raise FileNotFoundError(f"OMP parameter file not found: {path}")

    db = OMPDatabase()
    db.load(path)
    return db


def _parse_single_omp_file(filepath: str) -> Optional[OMPType]:
    """Parse a single OMP file (modified potential format).

    Args:
        filepath: Path to the OMP file (omp-XXXXXX.dat format)

    Returns:
        SphericalOMP or CoupledChannelOMP, or None if parsing fails
    """
    with open(filepath, 'r', encoding='utf-8', errors='replace') as fp:
        lines = fp.readlines()

    # Clean up lines
    entry_lines = []
    for line in lines:
        line = line.replace('\r', '').rstrip()
        # Stop at separator if present
        if line.startswith('+' * 20):
            break
        entry_lines.append(line)

    if len(entry_lines) < 11:
        return None

    try:
        # Get imodel from flags line to determine parser
        flags_line = entry_lines[10].split()
        imodel = int(flags_line[0])

        if imodel == 0:
            return parse_spherical_omp(entry_lines)
        else:
            return parse_coupled_channel_omp(entry_lines)
    except Exception:
        return None


def load_modified_potentials(directory: str = None, mod_dir: str = None) -> OMPDatabase:
    """Load modified optical model potentials from the mod-potentials directory.

    The mod-potentials directory contains individual OMP files (omp-XXXXXX.dat)
    that are updated or corrected versions of potentials in the main database.

    Args:
        directory: RIPL directory path. If provided, constructs full path.
        mod_dir: Direct path to mod-potentials directory. Takes precedence.

    Returns:
        OMPDatabase with modified potentials

    Raises:
        FileNotFoundError: If the mod-potentials directory cannot be found

    Example:
        >>> mod_db = load_modified_potentials()
        >>> print(mod_db.irefs)
        [1480, 1481, 1482, 2412, 2415, 4609, 4610]
    """
    if mod_dir:
        path = mod_dir
    elif directory:
        path = _os.path.join(directory, config.get_data_file_path('mod_potentials'))
    else:
        # Try default RIPL path
        from riplpy import config as riplpy_config
        ripl_path = riplpy_config.get_path()
        if ripl_path:
            path = _os.path.join(ripl_path, config.get_data_file_path('mod_potentials'))
        else:
            raise FileNotFoundError(
                "Cannot locate mod-potentials directory. Please provide directory or mod_dir."
            )

    if not _os.path.exists(path):
        raise FileNotFoundError(f"mod-potentials directory not found: {path}")

    db = OMPDatabase()
    db.data = {}

    # Find all .dat files in the directory
    for filename in _os.listdir(path):
        if filename.startswith('omp-') and filename.endswith('.dat'):
            filepath = _os.path.join(path, filename)
            omp = _parse_single_omp_file(filepath)
            if omp is not None:
                db.data[omp.iref] = omp
                if isinstance(omp, SphericalOMP):
                    db._n_spherical += 1
                else:
                    db._n_coupled_channel += 1

    return db


def load_with_modifications(directory: str = None, apply_modifications: bool = True) -> OMPDatabase:
    """Load the OMP database with optional modifications applied.

    This is a convenience function that loads the main database and optionally
    replaces entries with their modified versions from the mod-potentials directory.

    Args:
        directory: RIPL directory path.
        apply_modifications: If True, replace potentials with modified versions
                            where available. Default True.

    Returns:
        OMPDatabase with potentials (modified versions applied if requested)

    Example:
        >>> db = load_with_modifications()
        >>> # Potential 1480 will use the modified version if available
        >>> pot = db.get(1480)
    """
    # Load main database
    db = load(directory=directory)

    if apply_modifications:
        try:
            mod_db = load_modified_potentials(directory=directory)
            # Replace potentials with modified versions
            for iref, omp in mod_db.data.items():
                if iref in db.data:
                    # Update counts
                    old_omp = db.data[iref]
                    if isinstance(old_omp, SphericalOMP):
                        db._n_spherical -= 1
                    else:
                        db._n_coupled_channel -= 1
                    if isinstance(omp, SphericalOMP):
                        db._n_spherical += 1
                    else:
                        db._n_coupled_channel += 1
                db.data[iref] = omp
        except FileNotFoundError:
            # mod-potentials directory doesn't exist, skip
            pass

    return db
