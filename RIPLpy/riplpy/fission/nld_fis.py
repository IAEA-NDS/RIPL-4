# -*- coding: utf-8 -*-
"""Python objects which provide access to the BSkG3 fission saddle/well level
densities of RIPL-4 (``nld-fis-bskg3``).

The directory ``fission/nld-fis-bskg3`` contains five sub-directories of
per-Z files (``z090.dat`` ...):

* ``Max1``, ``Max2``, ``Max3`` -- inner / first-outer / second-outer saddle points.
* ``Min1``, ``Min2``           -- first / second superdeformed wells.

Each per-Z file concatenates several isotopes; per isotope there are TWO
blocks, positive parity then negative parity (S. Goriely, 19 August 2025).
Each block has a 3-line ``*****`` banner whose middle line carries the
metadata, e.g.::

    *  Z= 90 A=200: Positive-Parity Spin-dependent Level Density [MeV-1] for TH200 b2 = 0.66 b2 = 0.66 b30= 0.05 b40= 0.41 Icr=  25.48  *

followed by a ``U[MeV] T[MeV] NCUMUL RHOOBS RHOTOT J=00 ...`` column header
and data rows ``U T NCUMUL RHOOBS RHOTOT rho_J(50)``. The banner is parsed
by regex (Z, A, parity, b2, b30, b40, Icr). The number of spin columns
follows the file (50 for even-A integer spins, 50 half-integer for odd-A).

This reuses the parsing approach of :mod:`riplpy.densities.hfb`.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load the inner-saddle (Max1) level densities:
        $ import riplpy.fission.nld_fis as nld_fis
        $ db = nld_fis.load('Max1', directory='/path/to/RIPL-4')
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
    'CONFIG_KEY', 'SADDLES', 'parse_header', 'read_ascii_file',
    'ParityData', 'Entry', 'Database',
    'load_element', 'load',
)

# ========================

# Config key for the parent directory (subdirs Max1/Max2/Max3/Min1/Min2)
CONFIG_KEY = 'fission_nld_bskg3_dir'

# Recognised saddle/well sub-directories
SADDLES = ('Max1', 'Max2', 'Max3', 'Min1', 'Min2')

# Banner middle line:
# '*  Z= 90 A=200: Positive-Parity ... b2 = 0.66 b2 = 0.66 b30= 0.05 b40= 0.41 Icr=  25.48  *'
_HEADER_RE = _re.compile(
    r"Z\s*=\s*(\d+)\s+A\s*=\s*(\d+)\s*:\s*(Positive|Negative)-Parity"
    r".*?b2\s*=\s*(-?\d+(?:\.\d+)?)\s+b2\s*=\s*(-?\d+(?:\.\d+)?)"
    r"\s+b30\s*=\s*(-?\d+(?:\.\d+)?)\s+b40\s*=\s*(-?\d+(?:\.\d+)?)"
    r"\s+Icr\s*=\s*(-?\d+(?:\.\d+)?)",
    _re.IGNORECASE,
)

# ========================


def parse_header(line: str) -> tuple | None:
    """Parse a banner middle line.

    Returns ``(Z, A, parity, b2, b30, b40, Icr)`` where parity is +1/-1,
    or ``None`` if the line is not a recognisable banner.
    """
    m = _HEADER_RE.search(line)
    if not m:
        return None
    Z = int(m.group(1))
    A = int(m.group(2))
    parity = 1 if m.group(3).lower() == 'positive' else -1
    b2 = float(m.group(4))
    b30 = float(m.group(6))
    b40 = float(m.group(7))
    Icr = float(m.group(8))
    return Z, A, parity, b2, b30, b40, Icr


def read_ascii_file(fpath: str, label: str) -> dict:
    """Read a per-Z saddle/well level-density file and return a dictionary.

    Each isotope has a positive-parity block followed by a negative-parity
    block. The returned payloads carry the saddle/well ``label`` plus
    ``positive_parity`` and ``negative_parity`` tables.
    """
    d = {}

    cur_n = None
    cur_parity = None
    cur_meta = None
    cur_data = None
    in_data = False

    def _flush():
        if cur_n is None or cur_data is None:
            return
        if not cur_data['U']:
            return
        key = 'positive_parity' if cur_parity == 1 else 'negative_parity'
        d[cur_n][key] = {
            'b2': cur_meta[0], 'b30': cur_meta[1], 'b40': cur_meta[2],
            'Icr': cur_meta[3],
            'U': cur_data['U'], 'T': cur_data['T'],
            'NCUMUL': cur_data['NCUMUL'],
            'RHOOBS': cur_data['RHOOBS'], 'RHOTOT': cur_data['RHOTOT'],
            'rho_J': cur_data['rho_J'],
        }

    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for raw in fp:
            stripped = raw.strip()
            if not stripped:
                continue

            if stripped.startswith('*') and 'Z=' in stripped and 'A=' in stripped:
                hdr = parse_header(stripped)
                if hdr is None:
                    continue
                # Close any in-progress block before opening the next one.
                _flush()

                Z, A, parity, b2, b30, b40, Icr = hdr
                cur_n = _c.Nuclide(Z=Z, A=A)
                cur_parity = parity
                cur_meta = (b2, b30, b40, Icr)
                if cur_n not in d:
                    d[cur_n] = {
                        'n': cur_n, 'Z': Z, 'A': A,
                        'label': label,
                        'positive_parity': None,
                        'negative_parity': None,
                    }
                cur_data = {
                    'U': [], 'T': [], 'NCUMUL': [],
                    'RHOOBS': [], 'RHOTOT': [], 'rho_J': [],
                }
                in_data = False
                continue

            # Banner rule lines (only asterisks).
            if stripped.startswith('*'):
                continue

            # Column-header line.
            if 'U[MeV]' in stripped or 'RHOTOT' in stripped:
                in_data = True
                continue

            if in_data and cur_data is not None:
                tokens = stripped.split()
                if len(tokens) < 5:
                    continue
                try:
                    cur_data['U'].append(float(tokens[0]))
                    cur_data['T'].append(float(tokens[1]))
                    cur_data['NCUMUL'].append(float(tokens[2]))
                    cur_data['RHOOBS'].append(float(tokens[3]))
                    cur_data['RHOTOT'].append(float(tokens[4]))
                    cur_data['rho_J'].append([float(v) for v in tokens[5:]])
                except (ValueError, IndexError):
                    _logger.debug(f"Skipping malformed line: {stripped}")
                    continue

    # Finalise the last block.
    _flush()

    return d


@_dataclass
class ParityData:
    """Saddle/well level-density data for a single parity.

    Attributes:
        b2: Quadrupole axial deformation parameter
        b30: Octupole axial deformation parameter
        b40: Hexadecapole deformation parameter
        Icr: Cranking moment of inertia
        U: Excitation energy grid [MeV]
        T: Nuclear temperature [MeV]
        NCUMUL: Cumulative number of levels
        RHOOBS: Total (observed) level density [MeV^-1]
        RHOTOT: Total state density [MeV^-1]
        rho_J: Spin-dependent level densities (list of lists) [MeV^-1]
    """

    b2    : float = 0.0
    b30   : float = 0.0
    b40   : float = 0.0
    Icr   : float = 0.0
    U     : list = _field(default_factory=list)
    T     : list = _field(default_factory=list)
    NCUMUL: list = _field(default_factory=list)
    RHOOBS: list = _field(default_factory=list)
    RHOTOT: list = _field(default_factory=list)
    rho_J : list = _field(default_factory=list)

    _field_info: _ClassVar[dict] = {
        'b2':     'Quadrupole axial deformation parameter',
        'b30':    'Octupole axial deformation parameter',
        'b40':    'Hexadecapole deformation parameter',
        'Icr':    'Cranking moment of inertia',
        'U':      'Excitation energy grid [MeV]',
        'T':      'Nuclear temperature [MeV]',
        'NCUMUL': 'Cumulative number of levels',
        'RHOOBS': 'Total (observed) level density [MeV^-1]',
        'RHOTOT': 'Total state density [MeV^-1]',
        'rho_J':  'Spin-dependent level densities [MeV^-1]',
    }


def _to_parity(block: dict | None) -> ParityData | None:
    if block is None:
        return None
    return ParityData(
        b2=block['b2'], b30=block['b30'], b40=block['b40'], Icr=block['Icr'],
        U=block['U'], T=block['T'], NCUMUL=block['NCUMUL'],
        RHOOBS=block['RHOOBS'], RHOTOT=block['RHOTOT'], rho_J=block['rho_J'],
    )


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in a BSkG3 fission saddle/well level-density database.

    Attributes:
        n: Target nucleus
        Z: Charge number
        A: Mass number
        label: Saddle/well label (Max1, Max2, Max3, Min1, Min2)
        positive_parity: ParityData for positive parity (or None)
        negative_parity: ParityData for negative parity (or None)
    """

    n              : _c.Nuclide = None
    Z              : int = 0
    A              : int = 0
    label          : str = ""
    positive_parity: ParityData = None
    negative_parity: ParityData = None

    _field_info: _ClassVar[dict] = {
        'n':               'Target nucleus',
        'Z':               'Charge number',
        'A':               'Mass number',
        'label':           'Saddle/well label (Max1/Max2/Max3/Min1/Min2)',
        'positive_parity': 'Positive-parity level density (ParityData)',
        'negative_parity': 'Negative-parity level density (ParityData)',
    }


class Database(_db.Database):
    """A BSkG3 fission saddle/well level-density database (per-Z .dat files)."""

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = None

    def load_file(self, fpath: str, label: str) -> None:
        """Read ``fpath`` for saddle/well ``label`` and merge isotopes in."""
        data = type(self).reader(fpath, label)
        for key, payload in data.items():
            self.data[key] = self.entry(
                n=payload['n'], Z=payload['Z'], A=payload['A'],
                label=payload['label'],
                positive_parity=_to_parity(payload['positive_parity']),
                negative_parity=_to_parity(payload['negative_parity']),
            )

    def load_all(self, directory: str, label: str) -> None:
        """Load every ``z*.dat`` file in the ``label`` sub-directory."""
        sub = _os.path.join(directory, label)
        if not _os.path.isdir(sub):
            _logger.warning(f"Fission NLD directory not found: {sub}")
            return
        for fn in sorted(_os.listdir(sub)):
            if fn.endswith('.dat') and fn.startswith('z'):
                self.load_file(_os.path.join(sub, fn), label)


# ========================


def _base_directory(directory: str = None) -> str:
    """Resolve the absolute ``nld-fis-bskg3`` parent directory."""
    base = _config.resolve_directory(directory)
    rel = _config.get_data_file_path(CONFIG_KEY)
    return _os.path.join(base, rel)


def _check_saddle(saddle: str) -> str:
    if saddle not in SADDLES:
        raise ValueError(
            f"Unknown saddle/well {saddle!r}; expected one of {SADDLES}"
        )
    return saddle


def load_element(Z: int, saddle: str, directory: str = None) -> "Database":
    """Load a single per-Z file from the ``saddle`` sub-directory."""
    _check_saddle(saddle)
    sub = _os.path.join(_base_directory(directory), saddle)
    fpath = _os.path.join(sub, f"z{Z:03d}.dat")
    db = Database()
    if _os.path.exists(fpath):
        db.load_file(fpath, saddle)
    else:
        _logger.warning(f"Fission NLD file not found: {fpath}")
    return db


def load(saddle: str, directory: str = None) -> "Database":
    """Load all per-Z files for ``saddle`` (Max1/Max2/Max3/Min1/Min2)."""
    _check_saddle(saddle)
    db = Database()
    db.load_all(_base_directory(directory), saddle)
    return db
