# -*- coding: utf-8 -*-
"""Python objects which provide access to the RIPL-4 empirical fission barrier database.

This module reads ``empirical-barriers-ripl4.dat`` (Capote & Sin, 2025), which
contains empirical double-humped fission barrier parameters: inner and outer
barrier heights and curvatures with estimated uncertainties and saddle-point
symmetry designators.

Python 3.10+ is expected to run this code properly.

"""

# ========================

# Dataclasses
from dataclasses import dataclass as _dataclass

# Typing
from typing import ClassVar as _ClassVar

# Functools
from functools import partial as _partial

# RIPLpy
import riplpy.config as _config
import riplpy.collections as _c
import riplpy.db as _db

# ========================

__all__ = ('FILE_PATH_KEY', 'LOCAL_FILE_PATH', 'FILE_HEADER', 'read_ascii_file', 'write_ascii_file', 'Entry', 'Database', 'load')

# ========================

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'fission_barriers_empirical'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = (
    "#-----------------------------------------------------------------------\n"
    "#  Inner and outer barrier heights, widths with uncertainties are in MeV\n"
    "#  Z   A  s   Syma   Va  dVa  hwa dwa  Symb   Vb  dVb hwb dwb\n"
    "#-----------------------------------------------------------------------\n"
)

# ========================


def _parse_float(token: str | None):
    """Parse a token to float, returning None for blanks."""
    if token is None:
        return None
    token = token.strip()
    if not token:
        return None
    try:
        return float(token)
    except ValueError:
        return None


def read_ascii_file(fpath: str) -> dict:
    """Read data from the RIPL-4 empirical-barriers ASCII file and return a dictionary.

    The file may include records with only inner-barrier information, or with
    both inner and outer barrier data. Trailing fields are treated as optional.
    """
    d = {}
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            stripped = line.rstrip("\n")
            if not stripped.strip():
                continue
            if stripped.lstrip().startswith('#'):
                continue
            # Tokenize on whitespace; the format is whitespace separated in
            # the data records and tolerant of optional trailing fields.
            tokens = stripped.split()
            # Require at least: Z A sym Syma Va
            if len(tokens) < 5:
                continue
            try:
                Z = int(tokens[0])
                A = int(tokens[1])
            except ValueError:
                continue
            s = tokens[2]
            syma = tokens[3]
            Va = _parse_float(tokens[4])
            dVa = _parse_float(tokens[5]) if len(tokens) > 5 else None
            hwa = _parse_float(tokens[6]) if len(tokens) > 6 else None
            dwa = _parse_float(tokens[7]) if len(tokens) > 7 else None
            symb = tokens[8] if len(tokens) > 8 else None
            Vb = _parse_float(tokens[9]) if len(tokens) > 9 else None
            dVb = _parse_float(tokens[10]) if len(tokens) > 10 else None
            hwb = _parse_float(tokens[11]) if len(tokens) > 11 else None
            dwb = _parse_float(tokens[12]) if len(tokens) > 12 else None

            n = _c.Nuclide(Z=Z, A=A)
            d[n] = {
                'n': n, 's': s,
                'syma': syma, 'Va': Va, 'dVa': dVa, 'hwa': hwa, 'dwa': dwa,
                'symb': symb, 'Vb': Vb, 'dVb': dVb, 'hwb': hwb, 'dwb': dwb,
            }
    return d


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write the database to an ASCII file using a whitespace-separated layout."""
    with open(fpath, 'w', encoding='utf-8') as fp:
        fp.write(FILE_HEADER)
        for n in data.keys():
            e = data[n]
            line = f"{e.n.Z:4d}{e.n.A:4d} {e.s or e.n.element_symbol:<2s}  "
            line += f" {(e.syma or '  '):<2s}{(e.Va if e.Va is not None else 0.0):8.2f}{(e.dVa if e.dVa is not None else 0.0):4.1f}"
            if e.hwa is not None:
                line += f"{e.hwa:5.2f}{(e.dwa if e.dwa is not None else 0.0):5.2f}"
            if e.Vb is not None or e.symb:
                line += f"  {(e.symb or '  '):<2s}{(e.Vb if e.Vb is not None else 0.0):8.2f}{(e.dVb if e.dVb is not None else 0.0):4.1f}"
                if e.hwb is not None:
                    line += f"{e.hwb:5.2f}{(e.dwb if e.dwb is not None else 0.0):5.2f}"
            fp.write(line + "\n")


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the RIPL-4 empirical fission barrier database.

    Attributes:
        n: Target nucleus
        s: Element symbol
        syma: Inner saddle-point symmetry (S, GA, MA)
        Va: Inner barrier height [MeV]
        dVa: Inner barrier height uncertainty [MeV]
        hwa: Inner barrier curvature [MeV]
        dwa: Inner barrier curvature uncertainty [MeV]
        symb: Outer saddle-point symmetry (S, GA, MA)
        Vb: Outer barrier height [MeV]
        dVb: Outer barrier height uncertainty [MeV]
        hwb: Outer barrier curvature [MeV]
        dwb: Outer barrier curvature uncertainty [MeV]
    """

    n   : _c.Nuclide = None
    s   : str   = None
    syma: str   = None
    Va  : float = None
    dVa : float = None
    hwa : float = None
    dwa : float = None
    symb: str   = None
    Vb  : float = None
    dVb : float = None
    hwb : float = None
    dwb : float = None

    # Backwards-compatibility aliases (older API used Ea/Eb).
    @property
    def Ea(self) -> float:
        return self.Va

    @property
    def Eb(self) -> float:
        return self.Vb

    _field_info: _ClassVar[dict] = {
        'n':    'Target nucleus',
        's':    'Element symbol',
        'syma': 'Inner saddle-point symmetry (S/GA/MA)',
        'Va':   'Inner barrier height [MeV]',
        'dVa':  'Uncertainty on Va [MeV]',
        'hwa':  'Inner barrier curvature [MeV]',
        'dwa':  'Uncertainty on hwa [MeV]',
        'symb': 'Outer saddle-point symmetry (S/GA/MA)',
        'Vb':   'Outer barrier height [MeV]',
        'dVb':  'Uncertainty on Vb [MeV]',
        'hwb':  'Outer barrier curvature [MeV]',
        'dwb':  'Uncertainty on hwb [MeV]',
    }

    @property
    def as_tuple(self) -> tuple:
        return (self.n.Z, self.n.A, self.s or self.n.element_symbol,
                self.syma, self.Va, self.dVa, self.hwa, self.dwa,
                self.symb, self.Vb, self.dVb, self.hwb, self.dwb)


class Database(_db.Database):
    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file

    @property
    def heights(self) -> list:
        """Return the inner/outer barrier heights as a list of (Z, A, Va, Vb)."""
        ret = []
        for n in self.data.keys():
            ret.append((n.Z, n.A, self.data[n].Va, self.data[n].Vb))
        return ret

    @property
    def curvatures(self) -> list:
        """Return the inner/outer barrier curvatures as a list of (Z, A, hwa, hwb)."""
        ret = []
        for n in self.data.keys():
            ret.append((n.Z, n.A, self.data[n].hwa, self.data[n].hwb))
        return ret

    def first_barrier(self, n: _c.Nuclide) -> float:
        return self.data[n].Va

    def second_barrier(self, n: _c.Nuclide) -> float:
        return self.data[n].Vb


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
