# -*- coding: utf-8 -*-
"""Python objects which provide access to the theoretical BSkG3 fission barrier database of the Reference Input Parameter Library (RIPL).

The BSkG3 fission barriers are microscopic HFB-based predictions (inner,
outer, and isomer sections) shipped with both the full RIPL distribution and
the RIPL-4 GitHub release (``fission/barriers-bskg3.dat``).

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load the BSkG3 fission barrier database and read an actinide::
       $ import riplpy.fission as fission
       $ from riplpy.collections import Nuclide
       $ bskg3 = fission.bskg3.load()
       $ entry = bskg3.get(Nuclide(Z=92, A=238))
       $ print(entry.inner, entry.outer1)

"""

# ========================

# Dataclasses
from dataclasses import dataclass as _dataclass

# Typing
from typing import ClassVar as _ClassVar

# Functools
from functools import partial as _partial

# Fortran format
import fortranformat as _ff

# RIPLpy
import riplpy.config as _config
import riplpy.collections as _c
import riplpy.db as _db

# ========================

__all__ = ('FILE_PATH_KEY', 'LOCAL_FILE_PATH', 'FILE_HEADER', 'FORTRAN_FORMAT', 'reader', 'writer', 'read_ascii_file', 'write_ascii_file', 'Entry', 'Database', 'load')

# ========================

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'fission_barriers_bskg3'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = (
    "#             |        Inner            |            Outer 1         |         Outer 2           |     Isomer 1             |       Isomer 2\n"
    "# Z   N   A   |  E[MeV] B20   B22   B30 |  E[MeV] B20   B22   B30    |  E[MeV] B20   B22   B30   | E[MeV] B20   B22   B30   |   E[MeV] B20   B22   B30\n"
)
# Two file formats are shipped under RIPL-4:
#   - GitHub release: Z, N, A and five sections (Inner, Outer 1, Outer 2,
#     Isomer 1, Isomer 2) -- 23 columns total.
#   - Full RIPL-4 release: Z, N and four sections (Inner, Outer 1, Outer 2,
#     Isomer) -- 18 columns total, per ``barriers-bskg3.readme``.
# The reader autodetects which layout it is reading.
FORTRAN_FORMAT = '(i3,i4,i5,5(f10.3,3f6.2))'
reader = None  # parsing is whitespace-tolerant; fortranformat is not used
writer = _ff.FortranRecordWriter(FORTRAN_FORMAT)

# ========================


def _parse_section(tokens: list[str], offset: int) -> dict:
    return {
        'E[MeV]': float(tokens[offset]),
        'B20':    float(tokens[offset + 1]),
        'B22':    float(tokens[offset + 2]),
        'B30':    float(tokens[offset + 3]),
    }


def read_ascii_file(fpath: str) -> dict:
    """Read data from the BSkG3 ASCII file and return a dictionary.

    Two file layouts ship under different RIPL-4 distributions:

    * GitHub release (``barriers-bskg3.dat``, ~2449 nuclei): 23 columns laid
      out as ``Z N A`` plus five sections (Inner, Outer 1, Outer 2,
      Isomer 1, Isomer 2). Each section is ``E[MeV] B20 B22 B30``.
    * Full RIPL-4 release (``barriers-bskg3.dat``, 45 actinides): 18 columns
      laid out as ``Z N`` plus four sections (Inner, Outer 1, Outer 2,
      Isomer).

    The number of whitespace-separated tokens on each data row determines
    which layout is in effect. Lines that match neither are silently skipped.
    """
    d = {}
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            if not line.strip():
                continue
            if line.lstrip().startswith('#'):
                continue

            tokens = line.split()

            try:
                if len(tokens) >= 23:
                    # GitHub-release layout: Z, N, A + 5 sections
                    Z = int(tokens[0])
                    N = int(tokens[1])
                    A = int(tokens[2])
                    inner   = _parse_section(tokens, 3)
                    outer1  = _parse_section(tokens, 7)
                    outer2  = _parse_section(tokens, 11)
                    isomer  = _parse_section(tokens, 15)
                    isomer2 = _parse_section(tokens, 19)
                elif len(tokens) >= 18:
                    # Full RIPL-4 layout: Z, N + 4 sections (no second isomer)
                    Z = int(tokens[0])
                    N = int(tokens[1])
                    A = Z + N
                    inner   = _parse_section(tokens, 2)
                    outer1  = _parse_section(tokens, 6)
                    outer2  = _parse_section(tokens, 10)
                    isomer  = _parse_section(tokens, 14)
                    isomer2 = None
                else:
                    continue
            except (ValueError, IndexError):
                continue

            n = _c.Nuclide(Z=Z, A=A)
            d[n] = {
                'n': n,
                'inner': inner,
                'outer1': outer1,
                'outer2': outer2,
                'isomer': isomer,
                'isomer2': isomer2,
            }
    return d


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write the BSkG3 database to a whitespace-delimited ASCII file.

    Always emits the GitHub-release layout (Z, N, A plus five sections).
    The ``isomer2`` section falls back to all-zeros if the source entry
    came from the four-section legacy layout.
    """
    zero_section = {'E[MeV]': 0.0, 'B20': 0.0, 'B22': 0.0, 'B30': 0.0}
    with open(fpath, 'w', encoding='utf-8') as fp:
        fp.write(FILE_HEADER)
        for n in data.keys():
            e = data[n]
            isomer2 = e.isomer2 if getattr(e, 'isomer2', None) is not None else zero_section
            parts = [f"{e.n.Z:3d}{e.n.N:4d}{e.n.A:5d}"]
            for sec in (e.inner, e.outer1, e.outer2, e.isomer, isomer2):
                parts.append(f"{sec['E[MeV]']:10.3f}{sec['B20']:6.2f}{sec['B22']:6.2f}{sec['B30']:6.2f}")
            fp.write(''.join(parts) + "\n")


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the BSkG3 fission barrier database.

    Contains microscopic HFB-based fission barrier parameters from BSkG3,
    including inner, outer, and isomer saddle points with deformation.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        inner: Inner barrier {E[MeV], B20, B22, B30}
        outer1: First outer barrier {E[MeV], B20, B22, B30}
        outer2: Second outer barrier {E[MeV], B20, B22, B30}
        isomer: First shape isomer {E[MeV], B20, B22, B30}
        isomer2: Second shape isomer {E[MeV], B20, B22, B30}
            (None when sourced from the 4-section legacy file)
    """

    n      : _c.Nuclide = None
    inner  : dict = None
    outer1 : dict = None
    outer2 : dict = None
    isomer : dict = None
    isomer2: dict = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'inner': 'Inner barrier {E[MeV], B20, B22, B30}',
        'outer1': 'First outer barrier {E[MeV], B20, B22, B30}',
        'outer2': 'Second outer barrier {E[MeV], B20, B22, B30}',
        'isomer': 'First shape isomer {E[MeV], B20, B22, B30}',
        'isomer2': 'Second shape isomer {E[MeV], B20, B22, B30} (None on the legacy layout)',
    }

    @property
    def as_tuple(self) -> tuple:
        """Return a tuple representation matching the ASCII file format."""
        zero = {'E[MeV]': 0.0, 'B20': 0.0, 'B22': 0.0, 'B30': 0.0}
        iso2 = self.isomer2 if self.isomer2 is not None else zero
        return (self.n.Z, self.n.N, self.n.A,
                self.inner['E[MeV]'], self.inner['B20'], self.inner['B22'], self.inner['B30'],
                self.outer1['E[MeV]'], self.outer1['B20'], self.outer1['B22'], self.outer1['B30'],
                self.outer2['E[MeV]'], self.outer2['B20'], self.outer2['B22'], self.outer2['B30'],
                self.isomer['E[MeV]'], self.isomer['B20'], self.isomer['B22'], self.isomer['B30'],
                iso2['E[MeV]'], iso2['B20'], iso2['B22'], iso2['B30'])


class Database(_db.Database):

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
