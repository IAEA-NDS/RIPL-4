# -*- coding: utf-8 -*-
"""Python objects that provide access to the D1M mass model database.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load and print the D1M database::
       $ import riplpy.masses as masses
       $ d1m = masses.frdm.load() # The D1M database
       $ print(d1m.data)
    (2) Access a particular entry of the mass database using the get method::
       $ from riplpy.collections import Nuclide
       $ n = Nuclide(z=82, a=208)
       $ print(d1m.get(n))
    (3) 
"""

# ========================

# Dataclasses
from dataclasses import dataclass as _dataclass

# Typing
from typing import ClassVar as _ClassVar

# Functools
from functools import partial as _partial

# Logging
import logging as _logging

# Fortran format
import fortranformat as _ff

# RIPLpy
import riplpy.config as _config
import riplpy.collections as _c
import riplpy.db as _db
from . import core as _core

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = ('FILE_PATH_KEY', 'LOCAL_FILE_PATH', 'FILE_HEADER', 'FORTRAN_FORMAT', 'reader', 'writer', 'read_ascii_file', 'write_ascii_file', 'Entry', 'Database', 'load')

# ========================

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'mass_d1m'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = "#\n#  Z   A s fl     Mexp      Err       Mth    beta20   Rch\n#                [MeV]     [MeV]     [MeV]            [fm]\n#----------------------------------------------------------------------\n"
FORTRAN_FORMAT = '(2i4,1x,a2,1x,i1,3f10.3,f8.3,f8.3)'
reader = _ff.FortranRecordReader(FORTRAN_FORMAT)
writer = _ff.FortranRecordWriter(FORTRAN_FORMAT)

# ========================

def read_ascii_file(fpath: str) -> dict:
    """Read data from an ASCII file and return a dictionary."""
    data_dict = {}

    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp.readlines():
            # Skip header or empty lines
            if line.startswith("#") or not line.strip():
                continue
            try:
                data = reader.read(line)
                Z, A = data[0], data[1]
                n = _c.Nuclide(Z=Z, A=A)
                package = {
                    'n': n,
                    'flag': data[3],
                    'Mexp': data[4],
                    'Err': data[5],
                    'Mth': data[6],
                    'beta20': data[7],
                    'Rch': data[8],
                }
                data_dict[n] = package
            except Exception as e:
                _logger.warning("Failed to parse line: %s. Error: %s", line.strip(), e)
    return data_dict


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write ASCII file given a filename and a data dictionary."""
    with open(fpath, 'w', encoding='utf-8') as fp:
        # Write the header line
        fp.write(FILE_HEADER)
        for n, entry in data.items():
            try:
                line = (
                    n.Z,
                    n.A,
                    n.symbol,
                    entry.flag,
                    entry.Mexp,
                    entry.Err,
                    entry.Mth,
                    entry.beta20,
                    entry.Rch,
                )
                fp.write(writer.write(line) + "\n")
            except Exception as e:
                _logger.warning("Failed to write entry: %s. Error: %s", entry, e)


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the D1M mass database.

    Contains theoretical atomic masses from the D1M (Gogny) model along with
    experimental reference values and nuclear structure parameters.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        flag: Data quality flag
        Mexp: Experimental mass excess [MeV]
        Err: Uncertainty on mass excess [MeV]
        Mth: Calculated D1M mass excess [MeV]
        beta20: Quadrupole deformation
        Rch: Charge radius [fm]
    """
    n      : _c.Nuclide = None
    flag   : int = None
    Mexp   : float = None
    Err    : float = None
    Mth    : float = None
    beta20 : float = None
    Rch    : float = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'flag': 'Data quality flag',
        'Mexp': 'Experimental mass excess [MeV]',
        'Err': 'Uncertainty on mass excess [MeV]',
        'Mth': 'Calculated D1M mass excess [MeV]',
        'beta20': 'Quadrupole deformation',
        'Rch': 'Charge radius [fm]',
    }

    @classmethod
    def from_dict(cls, data: dict) -> "Entry":
        """Create an entry from a dictionary."""
        return cls(
            n=data['n'],
            flag=data['flag'],
            Mexp=data['Mexp'],
            Err=data['Err'],
            Mth=data['Mth'],
            beta20=data['beta20'],
            Rch=data['Rch'],
        )


class Database(_core.MassDatabase):
    """The D1M ground state masses database. """

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
