# -*- coding: utf-8 -*-
"""Parser for optical model deformation parameters (om-deformations.dat).

This module provides access to excited-level deformation parameters (beta-2, beta-3)
derived from experiments and nuclear data evaluations.
"""

import os as _os
from dataclasses import dataclass as _dataclass
from typing import Dict, List, Optional, Tuple, ClassVar as _ClassVar

from riplpy.db import Database
from . import config


@_dataclass
class DeformationEntry:
    """An entry in the deformation parameters database.

    Attributes:
        Z: Atomic number
        A: Mass number
        element: Element symbol
        Ex: Excitation energy [MeV]
        spin: Spin of the excited level (None if unknown)
        parity: Parity of the excited level (+1, -1, or 0 if unknown)
        L: Order of the deformation parameter (2 for beta-2, 3 for beta-3)
        beta: Deformation parameter value
        reference: Source of the data (e.g., 'JENDL-3.2', 'ENSDF(BE2)', 'Raman')
    """
    Z: int = 0
    A: int = 0
    element: str = ''
    Ex: float = 0.0
    spin: Optional[float] = None
    parity: int = 0
    L: int = 2
    beta: float = 0.0
    reference: str = ''

    _field_info: _ClassVar[dict] = {
        'Z': 'Atomic number',
        'A': 'Mass number',
        'element': 'Element symbol',
        'Ex': 'Excitation energy [MeV]',
        'spin': 'Spin of the excited level',
        'parity': 'Parity (+1, -1, or 0 if unknown)',
        'L': 'Order of deformation (2=quadrupole, 3=octupole)',
        'beta': 'Deformation parameter value',
        'reference': 'Data source',
    }

    @property
    def is_quadrupole(self) -> bool:
        """Check if this is a quadrupole (L=2) deformation."""
        return self.L == 2

    @property
    def is_octupole(self) -> bool:
        """Check if this is an octupole (L=3) deformation."""
        return self.L == 3

    @property
    def is_ground_state(self) -> bool:
        """Check if this is for the ground state."""
        return self.Ex < 0.001  # Essentially zero

    def __repr__(self) -> str:
        par_str = '+' if self.parity > 0 else ('-' if self.parity < 0 else '?')
        spin_str = f"{self.spin}" if self.spin is not None else "?"
        return (f"DeformationEntry({self.element}-{self.A}, Ex={self.Ex:.3f}, "
                f"J={spin_str}{par_str}, L={self.L}, beta={self.beta:.4f})")


class DeformationDatabase(Database):
    """Database of excited-level deformation parameters.

    This database contains beta-2 and beta-3 deformation parameters for
    collective levels, derived from experiments and nuclear data evaluations.
    """

    entry = DeformationEntry

    def __init__(self, data: Dict[Tuple[int, int, float, int], DeformationEntry] = None):
        """Initialize the deformation database.

        Args:
            data: Optional dictionary mapping (Z, A, Ex, L) -> DeformationEntry
        """
        super().__init__(data)

    def __repr__(self) -> str:
        return f"<DeformationDatabase: {len(self.data)} entries>"

    def get(self, Z: int, A: int, Ex: float = 0.0, L: int = 2) -> DeformationEntry:
        """Get a deformation entry.

        Args:
            Z: Atomic number
            A: Mass number
            Ex: Excitation energy [MeV] (default 0.0 for ground state)
            L: Deformation order (2 or 3, default 2)

        Returns:
            DeformationEntry for the specified level

        Raises:
            KeyError: If entry not found
        """
        # Try exact match first
        key = (Z, A, Ex, L)
        if key in self.data:
            return self.data[key]

        # Try finding closest Ex match
        best_key = None
        best_diff = float('inf')
        for k in self.data.keys():
            if k[0] == Z and k[1] == A and k[3] == L:
                diff = abs(k[2] - Ex)
                if diff < best_diff:
                    best_diff = diff
                    best_key = k

        if best_key is not None and best_diff < 0.01:  # Within 10 keV
            return self.data[best_key]

        from riplpy.exceptions import NucleusNotFoundError
        raise NucleusNotFoundError(
            f"Deformation not found for Z={Z}, A={A}, Ex={Ex}, L={L}"
        )

    def get_for_nucleus(self, Z: int, A: int) -> List[DeformationEntry]:
        """Get all deformation entries for a nucleus.

        Args:
            Z: Atomic number
            A: Mass number

        Returns:
            List of DeformationEntry objects for the nucleus
        """
        return [v for k, v in self.data.items() if k[0] == Z and k[1] == A]

    def get_ground_state(self, Z: int, A: int, L: int = 2) -> Optional[DeformationEntry]:
        """Get ground state deformation for a nucleus.

        Args:
            Z: Atomic number
            A: Mass number
            L: Deformation order (default 2)

        Returns:
            DeformationEntry for ground state, or None if not found
        """
        entries = self.get_for_nucleus(Z, A)
        for e in entries:
            if e.is_ground_state and e.L == L:
                return e
        return None

    def filter_by_element(self, Z: int) -> "DeformationDatabase":
        """Filter to entries for a specific element.

        Args:
            Z: Atomic number

        Returns:
            New DeformationDatabase with filtered entries
        """
        filtered = {k: v for k, v in self.data.items() if k[0] == Z}
        return DeformationDatabase(filtered)

    def filter_quadrupole(self) -> "DeformationDatabase":
        """Filter to quadrupole (L=2) deformations only.

        Returns:
            New DeformationDatabase with only L=2 entries
        """
        filtered = {k: v for k, v in self.data.items() if k[3] == 2}
        return DeformationDatabase(filtered)

    def filter_octupole(self) -> "DeformationDatabase":
        """Filter to octupole (L=3) deformations only.

        Returns:
            New DeformationDatabase with only L=3 entries
        """
        filtered = {k: v for k, v in self.data.items() if k[3] == 3}
        return DeformationDatabase(filtered)

    def info(self) -> str:
        """Return detailed information about this database."""
        lines = ["Deformation Database"]
        lines.append("-" * 40)
        lines.append(f"Total entries: {len(self.data)}")

        # Count by L
        l2_count = sum(1 for k in self.data.keys() if k[3] == 2)
        l3_count = sum(1 for k in self.data.keys() if k[3] == 3)
        lines.append(f"Quadrupole (L=2): {l2_count}")
        lines.append(f"Octupole (L=3): {l3_count}")

        # Count unique nuclei
        nuclei = set((k[0], k[1]) for k in self.data.keys())
        lines.append(f"Unique nuclei: {len(nuclei)}")

        # Z range
        if self.data:
            z_min = min(k[0] for k in self.data.keys())
            z_max = max(k[0] for k in self.data.keys())
            lines.append(f"Z range: {z_min} - {z_max}")

        return "\n".join(lines)


def _parse_deformation_line(line: str) -> Optional[DeformationEntry]:
    """Parse a single line from om-deformations.dat.

    Format: (2i4,1x,a2,1x,f10.6,1x,f4.1,i3,i2,1x,f10.6,2x,a13)

    Args:
        line: A line from the deformation file

    Returns:
        DeformationEntry if valid data line, None otherwise
    """
    # Skip comment and header lines
    line = line.replace('\r', '')
    if line.startswith('#') or not line.strip():
        return None

    try:
        # Fixed-width format parsing
        # Columns: Z(1-4), A(5-8), space, El(10-11), space, Ex(13-22), space,
        #          J(24-27), P(28-30), L(31-32), space, beta(34-43), spaces, Ref(46+)

        Z = int(line[0:4].strip())
        A = int(line[4:8].strip())
        element = line[9:11].strip()
        Ex = float(line[12:22].strip())

        # Spin may be blank
        spin_str = line[23:27].strip()
        spin = float(spin_str) if spin_str else None

        parity = int(line[27:30].strip()) if line[27:30].strip() else 0
        L = int(line[30:32].strip())
        beta = float(line[33:43].strip())
        reference = line[45:].strip()

        return DeformationEntry(
            Z=Z, A=A, element=element, Ex=Ex,
            spin=spin, parity=parity, L=L, beta=beta,
            reference=reference
        )
    except (ValueError, IndexError):
        return None


def read_deformations(filepath: str) -> DeformationDatabase:
    """Read the optical model deformations file.

    Args:
        filepath: Path to om-deformations.dat

    Returns:
        DeformationDatabase with all entries

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    if not _os.path.exists(filepath):
        raise FileNotFoundError(f"Deformations file not found: {filepath}")

    data = {}
    with open(filepath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            entry = _parse_deformation_line(line)
            if entry is not None:
                # Key: (Z, A, Ex, L)
                key = (entry.Z, entry.A, entry.Ex, entry.L)
                # Keep first entry if duplicate (some have multiple sources)
                if key not in data:
                    data[key] = entry

    return DeformationDatabase(data)


def load(directory: str = None, file_path: str = None) -> DeformationDatabase:
    """Load the optical model deformations database.

    Args:
        directory: RIPL directory path. If provided, constructs full path.
        file_path: Direct path to deformations file. Takes precedence over directory.

    Returns:
        DeformationDatabase with loaded entries

    Raises:
        FileNotFoundError: If the deformations file cannot be found
    """
    if file_path:
        return read_deformations(file_path)

    if directory:
        path = _os.path.join(directory, config.get_data_file_path('deformations'))
        return read_deformations(path)

    # Try default RIPL path
    from riplpy import config as riplpy_config
    ripl_path = riplpy_config.get_path()
    if ripl_path:
        path = _os.path.join(ripl_path, config.get_data_file_path('deformations'))
        return read_deformations(path)

    raise FileNotFoundError(
        "Cannot locate deformations file. Please provide directory or file_path."
    )
