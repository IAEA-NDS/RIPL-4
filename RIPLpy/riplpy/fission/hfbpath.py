# -*- coding: utf-8 -*-
"""Python objects which provide access to the HFB fission-path databases of RIPL-4.

Two microscopic models are supported, both sharing the same ASCII layout:

* ``hfbpath-bskg3`` -- BSkG3 fission paths (S. Goriely, 21 August 2025),
  ~2500 nuclei with 90 <= Z <= 118 (per-Z files ``z090.dat`` ... ``z118.dat``).
* ``hfbpath-d1m``   -- D1M fission paths (J.F. Lemaitre, March 2026),
  per-Z files ``z090.dat`` ... ``z096.dat``.

Each per-Z file concatenates several isotopes. Per isotope there is a title
line followed by ``nbeta`` data lines. The readme states the title Fortran
format is ``(4x,i4,4x,i4,1x,a2,9x,i4,7x,f10.3)`` (Z, A, symbol, nbeta, Egs)
and the data Fortran format is ``(5f11.3,6x,i2)`` (beta20, beta22, beta30,
E-Egs [MeV], MI_ATD [hbar^2/MeV], idx). Because the bskg3 files wrap the
title in a ``#`` comment banner and add a column-header comment line while
the d1m files use a single ``# Z= .. A= .. nbeta= .. Egs= ..`` line, the
header is parsed robustly by regex rather than rigid columns.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load the D1M fission paths:
        $ import riplpy.fission.hfbpath as hfbpath
        $ db = hfbpath.load_d1m(directory='/path/to/RIPL-4')
        $ print(db.data)

    (2) Load a single element from the BSkG3 model:
        $ db = hfbpath.load_element(Z=92, model='bskg3')

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
    'load_element', 'load_bskg3', 'load_d1m', 'load',
)

# ========================

# Map of model name -> config key for the per-Z directory
CONFIG_KEYS = {
    'bskg3': 'fission_path_bskg3_dir',
    'd1m':   'fission_path_d1m_dir',
}

# Data Fortran format is (5f11.3,6x,i2): five 11-char floats, 6 spaces, i2.
_DATA_WIDTHS = (11, 11, 11, 11, 11)


def _parse_data_line(stripped: str, raw: str) -> tuple | None:
    """Parse a data line, falling back to fixed Fortran columns.

    The file is written with Fortran format ``(5f11.3,6x,i2)``. Adjacent
    fixed-width fields can run together when a value overflows its 11-char
    slot (no separating space), so whitespace splitting may yield too few
    tokens. In that case slice by column widths instead.
    """
    tokens = stripped.split()
    if len(tokens) >= 6:
        try:
            return (
                float(tokens[0]), float(tokens[1]), float(tokens[2]),
                float(tokens[3]), float(tokens[4]), int(tokens[5]),
            )
        except (ValueError, IndexError):
            pass

    # Fixed-width fallback on the original (non-stripped) line.
    line = raw.rstrip("\n")
    try:
        vals = []
        pos = 0
        for w in _DATA_WIDTHS:
            vals.append(float(line[pos:pos + w]))
            pos += w
        idx = int(line[pos + 6:pos + 8])
        return (vals[0], vals[1], vals[2], vals[3], vals[4], idx)
    except (ValueError, IndexError):
        return None


# Header lines contain Z=, A=, nbeta and Egs; parse the numbers robustly.
_HEADER_RE = _re.compile(
    r"Z\s*=\s*(\d+).*?A\s*=\s*(\d+)\s*([A-Za-z]{1,2})?.*?"
    r"nbeta\s*=\s*(\d+).*?Egs\s*=\s*(-?\d+(?:\.\d+)?)",
    _re.IGNORECASE,
)

# ========================


def parse_header(line: str) -> tuple | None:
    """Parse a fission-path title line.

    Returns ``(Z, A, sym, nbeta, Egs)`` or ``None`` if the line is not a
    recognisable header (e.g. a banner or column-header comment).
    """
    m = _HEADER_RE.search(line)
    if not m:
        return None
    Z = int(m.group(1))
    A = int(m.group(2))
    sym = (m.group(3) or "").strip()
    nbeta = int(m.group(4))
    Egs = float(m.group(5))
    return Z, A, sym, nbeta, Egs


def read_ascii_file(fpath: str) -> dict:
    """Read a per-Z fission-path file and return a dictionary keyed by Nuclide.

    Each file concatenates several isotopes. A header line (containing Z=,
    A=, nbeta and Egs) is followed by ``nbeta`` data lines. Comment lines
    starting with ``#`` that are not headers (banner rules, column labels)
    are skipped.
    """
    d = {}

    cur_n = None
    cur = None
    nbeta_expected = 0

    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for raw in fp:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if not stripped:
                continue

            hdr = parse_header(stripped)
            if hdr is not None:
                Z, A, sym, nbeta, Egs = hdr
                cur_n = _c.Nuclide(Z=Z, A=A)
                cur = {
                    'n': cur_n,
                    'Z': Z, 'A': A, 'sym': sym,
                    'nbeta': nbeta, 'Egs': Egs,
                    'beta20': [], 'beta22': [], 'beta30': [],
                    'dE': [], 'MI': [], 'idx': [],
                }
                nbeta_expected = nbeta
                d[cur_n] = cur
                continue

            # Skip comment / banner lines that are not headers.
            if stripped.startswith('#') or stripped.startswith('*'):
                continue

            if cur is None:
                continue

            parsed = _parse_data_line(stripped, raw)
            if parsed is None:
                continue
            beta20, beta22, beta30, dE, MI, idx = parsed

            cur['beta20'].append(beta20)
            cur['beta22'].append(beta22)
            cur['beta30'].append(beta30)
            cur['dE'].append(dE)
            cur['MI'].append(MI)
            cur['idx'].append(idx)

            if len(cur['beta20']) >= nbeta_expected:
                # Done with this isotope; the next header opens a new one.
                cur = None
                cur_n = None

    return d


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in an HFB fission-path database.

    Attributes:
        n: Target nucleus
        Z: Charge number
        A: Mass number
        sym: Element symbol (from the title line, may be empty)
        nbeta: Number of beta deformation points declared in the header
        Egs: Ground-state binding energy [MeV]
        beta20: Quadrupole deformation parameter (parallel array)
        beta22: Quadrupole (triaxial) deformation parameter (parallel array)
        beta30: Octupole deformation parameter (parallel array)
        dE: Energy above the ground state E-Egs [MeV] (parallel array)
        MI: Collective inertia MI_ATD [hbar^2/MeV] wrt line number (parallel array)
        idx: Index flag: -1 ground state, 1 barrier, 2 well (parallel array)
    """

    n     : _c.Nuclide = None
    Z     : int = 0
    A     : int = 0
    sym   : str = ""
    nbeta : int = 0
    Egs   : float = 0.0
    beta20: list = _field(default_factory=list)
    beta22: list = _field(default_factory=list)
    beta30: list = _field(default_factory=list)
    dE    : list = _field(default_factory=list)
    MI    : list = _field(default_factory=list)
    idx   : list = _field(default_factory=list)

    _field_info: _ClassVar[dict] = {
        'n':      'Target nucleus',
        'Z':      'Charge number',
        'A':      'Mass number',
        'sym':    'Element symbol from the title line',
        'nbeta':  'Number of beta deformation points',
        'Egs':    'Ground-state binding energy [MeV]',
        'beta20': 'Quadrupole deformation parameter',
        'beta22': 'Quadrupole (triaxial) deformation parameter',
        'beta30': 'Octupole deformation parameter',
        'dE':     'Energy above ground state E-Egs [MeV]',
        'MI':     'Collective inertia MI_ATD [hbar^2/MeV] (wrt line number)',
        'idx':    'Index flag: -1 ground state, 1 barrier, 2 well',
    }


class Database(_db.Database):
    """An HFB fission-path database (per-Z .dat files, multi-isotope)."""

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = None

    def load_all(self, directory: str) -> None:
        """Load every ``z*.dat`` file found in ``directory``."""
        if not _os.path.isdir(directory):
            _logger.warning(f"HFB fission-path directory not found: {directory}")
            return
        for fn in sorted(_os.listdir(directory)):
            if fn.endswith('.dat') and fn.startswith('z'):
                self.load(_os.path.join(directory, fn))

    def load(self, fpath: str) -> None:
        """Read ``fpath`` and merge its isotopes into the database."""
        data = type(self).reader(fpath)
        for key, payload in data.items():
            self.data[key] = self.entry(**payload)


# ========================


def _model_directory(model: str, directory: str = None) -> str:
    """Resolve the absolute per-Z directory for ``model`` ('bskg3'|'d1m')."""
    model = model.lower()
    if model not in CONFIG_KEYS:
        raise ValueError(f"Unknown HFB fission-path model: {model!r} (use 'bskg3' or 'd1m')")
    base = _config.resolve_directory(directory)
    rel = _config.get_data_file_path(CONFIG_KEYS[model])
    return _os.path.join(base, rel)


def load_element(Z: int, model: str, directory: str = None) -> "Database":
    """Load a single per-Z file for ``model`` ('bskg3' or 'd1m')."""
    db_dir = _model_directory(model, directory)
    fpath = _os.path.join(db_dir, f"z{Z:03d}.dat")
    db = Database()
    if _os.path.exists(fpath):
        db.load(fpath)
    else:
        _logger.warning(f"HFB fission-path file not found: {fpath}")
    return db


def load_bskg3(directory: str = None) -> "Database":
    """Load all BSkG3 per-Z fission-path files."""
    db = Database()
    db.load_all(_model_directory('bskg3', directory))
    return db


def load_d1m(directory: str = None) -> "Database":
    """Load all D1M per-Z fission-path files."""
    db = Database()
    db.load_all(_model_directory('d1m', directory))
    return db


def load(model: str = 'bskg3', directory: str = None) -> "Database":
    """Load all per-Z fission-path files for ``model`` ('bskg3' or 'd1m')."""
    if model.lower() == 'd1m':
        return load_d1m(directory)
    return load_bskg3(directory)
