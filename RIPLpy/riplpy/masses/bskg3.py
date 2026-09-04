# -*- coding: utf-8 -*-
"""Python objects that provide access to the BSkG3 mass model database.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load and print the BSkG3 database::
       $ import riplpy.masses as masses
       $ bskg3 = masses.bskg3.load() # The BSkG3 database
       $ print(bskg3.data)
    (2) Access a particular entry of the mass database using the get method::
       $ from riplpy.collections import Nuclide
       $ n = Nuclide(z=82, a=208)
       $ print(bskg3.get(n))
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
FILE_PATH_KEY = 'mass_bskg3'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = "#\n#  Z   A s fl     Mexp      Err       Mth    beta20  beta22  Gamma  beta30  beta32  beta4    Rch\n#                [MeV]     [MeV]     [MeV]                                                   [fm]\n#----------------------------------------------------------------------------------------------------\n"
FORTRAN_FORMAT = '(2i4,1x,a2,1x,i1,3f10.3,7f8.3)'
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
                    'beta22': data[8],
                    'Gamma': data[9],
                    'beta30': data[10],
                    'beta32': data[11],
                    'beta4': data[12],
                    'Rch': data[13],
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
                    entry.beta22,
                    entry.Gamma,
                    entry.beta30,
                    entry.beta32,
                    entry.beta4,
                    entry.Rch,
                )
                fp.write(writer.write(line) + "\n")
            except Exception as e:
                _logger.warning("Failed to write entry: %s. Error: %s", entry, e)


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the BSkG3 mass database.

    Contains theoretical atomic masses from the BSkG3 model along with
    experimental reference values and nuclear structure parameters.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        flag: Data quality flag
        Mexp: Experimental mass excess [MeV]
        Err: Uncertainty on mass excess [MeV]
        Mth: Calculated BSkG3 mass excess [MeV]
        beta20: Quadrupole deformation (axial)
        beta22: Quadrupole deformation (triaxial)
        Gamma: Triaxial deformation angle [deg]
        beta30: Octupole deformation (axial)
        beta32: Octupole deformation (triaxial)
        beta4: Hexadecapole deformation
        Rch: Charge radius [fm]
    """
    n       : _c.Nuclide = None
    flag    : int = None
    Mexp    : float = None
    Err     : float = None
    Mth     : float = None
    beta20  : float = None
    beta22  : float = None
    Gamma   : float = None
    beta30  : float = None
    beta32  : float = None
    beta4   : float = None
    Rch     : float = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'flag': 'Data quality flag',
        'Mexp': 'Experimental mass excess [MeV]',
        'Err': 'Uncertainty on mass excess [MeV]',
        'Mth': 'Calculated BSkG3 mass excess [MeV]',
        'beta20': 'Quadrupole deformation (axial)',
        'beta22': 'Quadrupole deformation (triaxial)',
        'Gamma': 'Triaxial deformation angle [deg]',
        'beta30': 'Octupole deformation (axial)',
        'beta32': 'Octupole deformation (triaxial)',
        'beta4': 'Hexadecapole deformation',
        'Rch': 'Charge radius [fm]',
    }

    @classmethod
    def from_dict(cls, data: dict) -> "Entry":
        """Create an entry from a dictionary."""
        return cls(**data)


class Database(_core.MassDatabase):
    """The BSkG3 ground state masses database. """

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
