# -*- coding: utf-8 -*-
"""Parser for the optical model index file (om-index.txt).

The index file provides a summary of all optical model potentials in the library,
including metadata about projectile type, model type, validity ranges, and references.
"""

import os as _os
from dataclasses import dataclass as _dataclass, field as _field
from typing import List, Dict, Optional, ClassVar as _ClassVar

import riplpy.db as _db
from . import config


@_dataclass
class IndexEntry:
    """An entry in the optical model index.

    Represents metadata for a single optical model potential.

    Attributes:
        iref: Library reference number (unique ID)
        projectile: Incident particle ('n', 'p', 'd', 't', 'He3', 'a')
        model_type: Model type string ('spher.', 'CC rot.', 'CC vib.', etc.)
        dispersive: Whether dispersion relations are used
        relativistic: Whether relativistic kinematics are used
        Z_min: Minimum target atomic number
        Z_max: Maximum target atomic number
        A_min: Minimum target mass number
        A_max: Maximum target mass number
        E_min: Minimum energy [MeV]
        E_max: Maximum energy [MeV]
        ref_num: Reference number
        first_author: First author name
    """
    iref: int = 0
    projectile: str = ''
    model_type: str = ''
    dispersive: bool = False
    relativistic: bool = False
    Z_min: int = 0
    Z_max: int = 0
    A_min: int = 0
    A_max: int = 0
    E_min: float = 0.0
    E_max: float = 0.0
    ref_num: int = 0
    first_author: str = ''

    _field_info: _ClassVar[dict] = {
        'iref': 'Library reference number (unique ID)',
        'projectile': "Incident particle ('n', 'p', 'd', 't', 'He3', 'a')",
        'model_type': "Model type ('spher.', 'CC rot.', 'CC vib.', 'CC soft')",
        'dispersive': 'Uses dispersion relations',
        'relativistic': 'Uses relativistic kinematics',
        'Z_min': 'Minimum target atomic number',
        'Z_max': 'Maximum target atomic number',
        'A_min': 'Minimum target mass number',
        'A_max': 'Maximum target mass number',
        'E_min': 'Minimum valid energy [MeV]',
        'E_max': 'Maximum valid energy [MeV]',
        'ref_num': 'Reference number',
        'first_author': 'First author name',
    }

    @property
    def is_spherical(self) -> bool:
        """Check if this is a spherical potential."""
        return 'spher' in self.model_type.lower()

    @property
    def is_coupled_channel(self) -> bool:
        """Check if this is a coupled-channel potential."""
        return 'cc' in self.model_type.lower()

    @property
    def is_neutron(self) -> bool:
        """Check if projectile is neutron."""
        return self.projectile == 'n'

    @property
    def is_proton(self) -> bool:
        """Check if projectile is proton."""
        return self.projectile == 'p'

    def valid_for(self, Z: int, A: int, E: float) -> bool:
        """Check if potential is valid for given target and energy.

        Args:
            Z: Target atomic number
            A: Target mass number
            E: Projectile energy [MeV]

        Returns:
            True if potential is valid for this target/energy
        """
        return (self.Z_min <= Z <= self.Z_max and
                self.A_min <= A <= self.A_max and
                self.E_min <= E <= self.E_max)

    def __repr__(self) -> str:
        return (f"IndexEntry(iref={self.iref}, {self.projectile}, {self.model_type}, "
                f"Z={self.Z_min}-{self.Z_max}, A={self.A_min}-{self.A_max}, "
                f"E={self.E_min:.1f}-{self.E_max:.1f})")


class OMPIndex(_db.Database):
    """Database of optical model potential index entries.

    This provides fast lookup and filtering of OMP metadata without
    loading the full potential parameters.
    """

    entry = IndexEntry

    def __init__(self, data: Dict[int, IndexEntry] = None):
        """Initialize the index database.

        Args:
            data: Dictionary mapping iref -> IndexEntry
        """
        super().__init__(data)

    def get(self, iref: int) -> IndexEntry:
        """Get index entry by reference number.

        Args:
            iref: Library reference number

        Returns:
            IndexEntry for the specified potential

        Raises:
            KeyError: If iref not found
        """
        if iref not in self.data:
            raise KeyError(f"OMP reference {iref} not found in index")
        return self.data[iref]

    def filter_by_projectile(self, projectile: str) -> "OMPIndex":
        """Filter index by projectile type.

        Args:
            projectile: Projectile name ('n', 'p', 'd', 't', 'He3', 'a')

        Returns:
            New OMPIndex with filtered entries
        """
        filtered = {k: v for k, v in self.data.items() if v.projectile == projectile}
        return OMPIndex(filtered)

    def filter_by_target(self, Z: int = None, A: int = None) -> "OMPIndex":
        """Filter index by target nucleus.

        Args:
            Z: Target atomic number (filters Z_min <= Z <= Z_max)
            A: Target mass number (filters A_min <= A <= A_max)

        Returns:
            New OMPIndex with filtered entries
        """
        filtered = {}
        for k, v in self.data.items():
            if Z is not None and not (v.Z_min <= Z <= v.Z_max):
                continue
            if A is not None and not (v.A_min <= A <= v.A_max):
                continue
            filtered[k] = v
        return OMPIndex(filtered)

    def filter_by_energy(self, E: float) -> "OMPIndex":
        """Filter index by projectile energy.

        Args:
            E: Projectile energy [MeV]

        Returns:
            New OMPIndex with filtered entries
        """
        filtered = {k: v for k, v in self.data.items() if v.E_min <= E <= v.E_max}
        return OMPIndex(filtered)

    def filter_spherical(self) -> "OMPIndex":
        """Filter to spherical potentials only.

        Returns:
            New OMPIndex with only spherical potentials
        """
        filtered = {k: v for k, v in self.data.items() if v.is_spherical}
        return OMPIndex(filtered)

    def filter_coupled_channel(self) -> "OMPIndex":
        """Filter to coupled-channel potentials only.

        Returns:
            New OMPIndex with only coupled-channel potentials
        """
        filtered = {k: v for k, v in self.data.items() if v.is_coupled_channel}
        return OMPIndex(filtered)

    def filter_dispersive(self) -> "OMPIndex":
        """Filter to dispersive potentials only.

        Returns:
            New OMPIndex with only dispersive potentials
        """
        filtered = {k: v for k, v in self.data.items() if v.dispersive}
        return OMPIndex(filtered)

    def find_for_target(self, projectile: str, Z: int, A: int, E: float) -> List[IndexEntry]:
        """Find all potentials valid for a specific reaction.

        Args:
            projectile: Projectile type ('n', 'p', etc.)
            Z: Target atomic number
            A: Target mass number
            E: Projectile energy [MeV]

        Returns:
            List of IndexEntry objects valid for this reaction
        """
        results = []
        for entry in self.data.values():
            if entry.projectile == projectile and entry.valid_for(Z, A, E):
                results.append(entry)
        return results

    def list_projectiles(self) -> List[str]:
        """Get list of unique projectile types.

        Returns:
            List of projectile names
        """
        return sorted(set(v.projectile for v in self.data.values()))

    def list_authors(self) -> List[str]:
        """Get list of unique first authors.

        Returns:
            List of author names
        """
        return sorted(set(v.first_author for v in self.data.values()))

    def summary(self) -> str:
        """Get a summary of the index contents.

        Returns:
            Summary string
        """
        lines = [f"OMP Index: {len(self.data)} potentials"]
        lines.append("-" * 40)

        # Count by projectile
        proj_counts = {}
        for v in self.data.values():
            proj_counts[v.projectile] = proj_counts.get(v.projectile, 0) + 1
        lines.append("By projectile:")
        for proj, count in sorted(proj_counts.items()):
            lines.append(f"  {proj}: {count}")

        # Count by model type
        model_counts = {}
        for v in self.data.values():
            model_counts[v.model_type] = model_counts.get(v.model_type, 0) + 1
        lines.append("By model type:")
        for model, count in sorted(model_counts.items()):
            lines.append(f"  {model}: {count}")

        # Count dispersive/relativistic
        n_disp = sum(1 for v in self.data.values() if v.dispersive)
        n_rel = sum(1 for v in self.data.values() if v.relativistic)
        lines.append(f"Dispersive: {n_disp}")
        lines.append(f"Relativistic: {n_rel}")

        return "\n".join(lines)


def _parse_index_line(line: str) -> Optional[IndexEntry]:
    """Parse a single line from the index file.

    Args:
        line: A line from om-index.txt

    Returns:
        IndexEntry if line contains valid data, None otherwise
    """
    import re

    # Remove DOS line endings
    line = line.replace('\r', '')

    # Skip header lines and empty lines
    if not line.strip() or line.strip().startswith('Lib') or line.strip().startswith('No.'):
        return None

    # The format is semi-fixed width. Use regex to parse.
    # Example lines:
    #   15   n    spher.   no   no    6- 6   12- 12    .0- 65.0     1  M.B.Chadwick
    # 2602   n    CC soft  no   yes  12-50   24-124    .0-200.0     7  E.Soukhovitskii
    # 1484   n    CC rig.  yes  yes  25-25   50- 56    .0-200.0    26  E.Soukhovitskii

    # Pattern breakdown:
    # (\d+)\s+           - iref (library number)
    # (\S+)\s+           - projectile (n, p, d, t, 3He, 4He)
    # ([\w\s.]+?)\s+     - model type (spher., CC soft, CC rig., CC vib.)
    # (yes|no)\s+        - dispersive
    # (yes|no)\s+        - relativistic
    # (\d+)-\s*(\d+)\s+  - Z range
    # (\d+)-\s*(\d+)\s+  - A range
    # ([\d.]+)-([\d.]+)\s+ - E range
    # (\d+)\s+           - reference number
    # (.+)               - first author

    pattern = (
        r'^\s*(\d+)\s+'           # iref
        r'(\S+)\s+'               # projectile
        r'(spher\.|CC\s+[\w/]+\.?)\s+'  # model type (includes CC r/vi)
        r'(yes|no)\s+'            # dispersive
        r'(yes|no)\s+'            # relativistic
        r'(\d+)\s*-\s*(\d+)\s+'   # Z range
        r'(\d+)\s*-\s*(\d+)\s+'   # A range
        r'([\d.]+)\s*-\s*([\d.]+)\s+'  # E range
        r'(\d+)\s+'               # ref number
        r'(.+?)\s*$'              # author
    )

    match = re.match(pattern, line)
    if not match:
        return None

    try:
        return IndexEntry(
            iref=int(match.group(1)),
            projectile=match.group(2),
            model_type=match.group(3).strip(),
            dispersive=match.group(4).lower() == 'yes',
            relativistic=match.group(5).lower() == 'yes',
            Z_min=int(match.group(6)),
            Z_max=int(match.group(7)),
            A_min=int(match.group(8)),
            A_max=int(match.group(9)),
            E_min=float(match.group(10)),
            E_max=float(match.group(11)),
            ref_num=int(match.group(12)),
            first_author=match.group(13).strip(),
        )
    except (ValueError, IndexError):
        return None


def read_index(filepath: str) -> OMPIndex:
    """Read the optical model index file.

    Args:
        filepath: Path to om-index.txt

    Returns:
        OMPIndex database with all entries

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    if not _os.path.exists(filepath):
        raise FileNotFoundError(f"Index file not found: {filepath}")

    data = {}
    with open(filepath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            entry = _parse_index_line(line)
            if entry is not None:
                data[entry.iref] = entry

    return OMPIndex(data)


def load(directory: str = None, file_path: str = None) -> OMPIndex:
    """Load the optical model index.

    Args:
        directory: RIPL directory path. If provided, constructs full path.
        file_path: Direct path to index file. Takes precedence over directory.

    Returns:
        OMPIndex database

    Raises:
        FileNotFoundError: If the index file cannot be found
    """
    if file_path:
        return read_index(file_path)

    if directory:
        path = _os.path.join(directory, config.get_data_file_path('index'))
        return read_index(path)

    # Try default RIPL path
    from riplpy import config as riplpy_config
    ripl_path = riplpy_config.get_path()
    if ripl_path:
        path = _os.path.join(ripl_path, config.get_data_file_path('index'))
        return read_index(path)

    raise FileNotFoundError("Cannot locate OMP index file. Please provide directory or file_path.")
