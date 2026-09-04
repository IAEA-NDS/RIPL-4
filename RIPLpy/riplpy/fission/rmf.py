# -*- coding: utf-8 -*-
"""Python objects which provide access to the RMF fission-path databases of RIPL-4.

Two layouts share the same per-nucleus ASCII format:

* ``fission/RMF/Path_Axial``    -- axial RMF fission paths.
* ``fission/RMF/Path_Triaxial`` -- triaxial RMF fission paths.

There is no RMF readme; the format is derived from the files. Each file is
named after a single nucleus (element symbol + mass number, e.g. ``Am239.dat``)
and contains a ``#``-banner header::

    #--------------------------------------------------------------------
    # Z = 95  A = 239 Am  nbeta= 131  Egs= 0.000 MeV  E(0+)= 0.0000 MeV
    #     beta20      beta22      beta30       E-EGS          Mu
    #--------------------------------------------------------------------

followed by ``nbeta`` data rows of ``beta20 beta22 beta30 E-EGS[MeV] Mu``.
The header is parsed by regex; the filename gives sym+A as a cross-check
while Z is trusted from the header.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load the axial RMF fission paths:
        $ import riplpy.fission.rmf as rmf
        $ db = rmf.load_axial(directory='/path/to/RIPL-4')
        $ print(db.data)

"""

# ========================

# OS
import os as _os

# Regular expressions
import re as _re

# Dataclasses
from dataclasses import dataclass as _dataclass
from dataclasses import field as _field

# Typing
from typing import ClassVar as _ClassVar

# Logging
import logging as _logging

# RIPLpy
import riplpy.config as _config
import riplpy.collections as _c
import riplpy.db as _db

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = (
    'CONFIG_KEYS', 'parse_header', 'read_ascii_file',
    'Entry', 'Database',
    'load_axial', 'load_triaxial', 'load',
)

# ========================

# Map of layout name -> config key for the per-nucleus directory
CONFIG_KEYS = {
    'axial':    'fission_rmf_axial_dir',
    'triaxial': 'fission_rmf_triaxial_dir',
}

# Header: '# Z = 95  A = 239 Am  nbeta= 131  Egs= 0.000 MeV  E(0+)= 0.0000 MeV'
_HEADER_RE = _re.compile(
    r"Z\s*=\s*(\d+).*?A\s*=\s*(\d+)\s*([A-Za-z]{1,2})?.*?"
    r"nbeta\s*=\s*(\d+).*?Egs\s*=\s*(-?\d+(?:\.\d+)?).*?"
    r"E\(0\+\)\s*=\s*(-?\d+(?:\.\d+)?)",
    _re.IGNORECASE,
)

# ========================


def parse_header(line: str) -> tuple | None:
    """Parse an RMF fission-path header line.

    Returns ``(Z, A, sym, nbeta, Egs, E0)`` or ``None`` if not a header.
    """
    m = _HEADER_RE.search(line)
    if not m:
        return None
    Z = int(m.group(1))
    A = int(m.group(2))
    sym = (m.group(3) or "").strip()
    nbeta = int(m.group(4))
    Egs = float(m.group(5))
    E0 = float(m.group(6))
    return Z, A, sym, nbeta, Egs, E0


def read_ascii_file(fpath: str) -> dict:
    """Read a per-nucleus RMF fission-path file and return a dictionary."""
    d = {}

    cur_n = None
    cur = None

    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for raw in fp:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith('#'):
                hdr = parse_header(stripped)
                if hdr is not None:
                    Z, A, sym, nbeta, Egs, E0 = hdr
                    cur_n = _c.Nuclide(Z=Z, A=A)
                    cur = {
                        'n': cur_n,
                        'Z': Z, 'A': A, 'sym': sym,
                        'nbeta': nbeta, 'Egs': Egs, 'E0': E0,
                        'beta20': [], 'beta22': [], 'beta30': [],
                        'dE': [], 'Mu': [],
                    }
                    d[cur_n] = cur
                # Banner / column-header comment lines are skipped.
                continue

            if cur is None:
                continue

            tokens = stripped.split()
            if len(tokens) < 5:
                continue
            try:
                beta20 = float(tokens[0])
                beta22 = float(tokens[1])
                beta30 = float(tokens[2])
                dE = float(tokens[3])
                Mu = float(tokens[4])
            except (ValueError, IndexError):
                continue

            cur['beta20'].append(beta20)
            cur['beta22'].append(beta22)
            cur['beta30'].append(beta30)
            cur['dE'].append(dE)
            cur['Mu'].append(Mu)

    return d


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in an RMF fission-path database.

    Attributes:
        n: Target nucleus
        Z: Charge number (trusted from the header)
        A: Mass number
        sym: Element symbol from the header
        nbeta: Number of beta deformation points declared in the header
        Egs: Ground-state energy [MeV]
        E0: E(0+) band-head energy [MeV]
        beta20: Quadrupole deformation parameter (parallel array)
        beta22: Quadrupole (triaxial) deformation parameter (parallel array)
        beta30: Octupole deformation parameter (parallel array)
        dE: Energy relative to ground state E-EGS [MeV] (parallel array)
        Mu: Collective inertia Mu (parallel array)
    """

    n     : _c.Nuclide = None
    Z     : int = 0
    A     : int = 0
    sym   : str = ""
    nbeta : int = 0
    Egs   : float = 0.0
    E0    : float = 0.0
    beta20: list = _field(default_factory=list)
    beta22: list = _field(default_factory=list)
    beta30: list = _field(default_factory=list)
    dE    : list = _field(default_factory=list)
    Mu    : list = _field(default_factory=list)

    _field_info: _ClassVar[dict] = {
        'n':      'Target nucleus',
        'Z':      'Charge number',
        'A':      'Mass number',
        'sym':    'Element symbol from the header',
        'nbeta':  'Number of beta deformation points',
        'Egs':    'Ground-state energy [MeV]',
        'E0':     'E(0+) band-head energy [MeV]',
        'beta20': 'Quadrupole deformation parameter',
        'beta22': 'Quadrupole (triaxial) deformation parameter',
        'beta30': 'Octupole deformation parameter',
        'dE':     'Energy relative to ground state E-EGS [MeV]',
        'Mu':     'Collective inertia Mu',
    }


class Database(_db.Database):
    """An RMF fission-path database (per-nucleus .dat files)."""

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = None

    def load_all(self, directory: str) -> None:
        """Load every ``*.dat`` file found in ``directory``."""
        if not _os.path.isdir(directory):
            _logger.warning(f"RMF fission-path directory not found: {directory}")
            return
        for fn in sorted(_os.listdir(directory)):
            if fn.endswith('.dat'):
                self.load(_os.path.join(directory, fn))

    def load(self, fpath: str) -> None:
        """Read ``fpath`` and merge its nucleus into the database."""
        data = type(self).reader(fpath)
        for key, payload in data.items():
            self.data[key] = self.entry(**payload)


# ========================


def _layout_directory(layout: str, directory: str = None) -> str:
    """Resolve the absolute per-nucleus directory for ``layout``."""
    layout = layout.lower()
    if layout not in CONFIG_KEYS:
        raise ValueError(f"Unknown RMF layout: {layout!r} (use 'axial' or 'triaxial')")
    base = _config.resolve_directory(directory)
    rel = _config.get_data_file_path(CONFIG_KEYS[layout])
    return _os.path.join(base, rel)


def load_axial(directory: str = None) -> "Database":
    """Load all axial RMF fission-path files."""
    db = Database()
    db.load_all(_layout_directory('axial', directory))
    return db


def load_triaxial(directory: str = None) -> "Database":
    """Load all triaxial RMF fission-path files."""
    db = Database()
    db.load_all(_layout_directory('triaxial', directory))
    return db


def load(layout: str = 'axial', directory: str = None) -> "Database":
    """Load all RMF fission-path files for ``layout`` ('axial' or 'triaxial')."""
    if layout.lower() == 'triaxial':
        return load_triaxial(directory)
    return load_axial(directory)
