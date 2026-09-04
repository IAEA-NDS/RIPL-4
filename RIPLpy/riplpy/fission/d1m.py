# -*- coding: utf-8 -*-
"""Python objects which provide access to the D1M HFB fission barrier database of RIPL-4.

This module reads ``barriers-d1m_lep.dat`` (Lemaitre, 2026), which contains
HFB fission barrier heights and deformation parameters obtained within the
Energy Density Functional theory with the D1M effective Gogny interaction,
extracted from the 2-dimensional least-energy path (LEP).

Python 3.10+ is expected to run this code properly.

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
FILE_PATH_KEY = 'fission_barriers_d1m'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = (
    "             |        Inner            |            Outer\n"
    " Z   N   A   | Binner  B20   B22   B30 |   Bouter  B20   B22   B30\n"
)
FORTRAN_FORMAT = '(i3,2i4,2(f10.3,3f6.2))'
reader = _ff.FortranRecordReader(FORTRAN_FORMAT)
writer = _ff.FortranRecordWriter(FORTRAN_FORMAT)

# ========================


def read_ascii_file(fpath: str) -> dict:
    """Read data from the D1M LEP ASCII file and return a dictionary."""
    d = {}
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            stripped = line.rstrip("\n")
            if not stripped.strip():
                continue
            # Skip header lines (start with whitespace + non-digit) or contain
            # the column labels.
            first = stripped.lstrip()[:1]
            if first in ('|',) or 'Inner' in stripped or 'Binner' in stripped:
                continue
            if not (first.isdigit() or first == '-'):
                continue

            # Tokenize the line and rely on whitespace separation; the file
            # is space-delimited with fixed-column data.
            tokens = stripped.split()
            if len(tokens) < 9:
                continue
            try:
                Z = int(tokens[0])
                N = int(tokens[1])
                A = int(tokens[2])
                inner_e = float(tokens[3])
                inner_b20 = float(tokens[4])
                inner_b22 = 0.0  # D1M LEP file does not include B22
                inner_b30 = float(tokens[5])
                outer_e = float(tokens[6])
                outer_b20 = float(tokens[7])
                outer_b22 = 0.0
                outer_b30 = float(tokens[8])
            except (ValueError, IndexError):
                continue

            # The readme labels columns Binner, B20, B30, Bouter, B20, B30
            # (no B22 in the D1M LEP file).
            inner = {'E[MeV]': inner_e, 'B20': inner_b20, 'B22': inner_b22, 'B30': inner_b30}
            outer = {'E[MeV]': outer_e, 'B20': outer_b20, 'B22': outer_b22, 'B30': outer_b30}

            n = _c.Nuclide(Z=Z, A=A)
            d[n] = {'n': n, 'inner': inner, 'outer': outer}
    return d


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write the database to an ASCII file in the D1M LEP layout."""
    with open(fpath, 'w', encoding='utf-8') as fp:
        fp.write(FILE_HEADER)
        for n in data.keys():
            e = data[n]
            row = (
                e.n.Z, e.n.N, e.n.A,
                e.inner['E[MeV]'], e.inner['B20'], e.inner['B22'], e.inner['B30'],
                e.outer['E[MeV]'], e.outer['B20'], e.outer['B22'], e.outer['B30'],
            )
            fp.write(writer.write(row) + "\n")


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the D1M HFB fission barrier database.

    Attributes:
        n: Target nucleus
        inner: Inner barrier {E[MeV], B20, B22, B30}
        outer: Outer barrier {E[MeV], B20, B22, B30}
    """

    n    : _c.Nuclide = None
    inner: dict = None
    outer: dict = None

    _field_info: _ClassVar[dict] = {
        'n':     'Target nucleus',
        'inner': 'Inner barrier {E[MeV], B20, B22, B30}',
        'outer': 'Outer barrier {E[MeV], B20, B22, B30}',
    }

    @property
    def as_tuple(self) -> tuple:
        return (self.n.Z, self.n.N, self.n.A,
                self.inner['E[MeV]'], self.inner['B20'], self.inner['B22'], self.inner['B30'],
                self.outer['E[MeV]'], self.outer['B20'], self.outer['B22'], self.outer['B30'])


class Database(_db.Database):
    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
