# -*- coding: utf-8 -*-
"""Python objects that provide access to the natural abundances. 

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load and print the natural abundance database::
       $ import riplpy.masses as masses
       $ natab = masses.ab.load() # The natural abundances database
       $ print(natab.data)
    (2) Access a particular entry of the natural abundances database using the get method::
       $ from riplpy.collections import Nuclide
       $ n = Nuclide(z=50, a=120)
       $ print(natab.get(n))
    (3) 
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
FILE_PATH_KEY = 'abundances'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER    = "#      Natural abundances\n#  Z  A  El   abundance  uncert.\n#               [%]       [%]\n#-------------------------------\n"
FORTRAN_FORMAT = '(2i4,1x,a2,1x,2f10.6)'
reader   = _ff.FortranRecordReader(FORTRAN_FORMAT)
writer   = _ff.FortranRecordWriter(FORTRAN_FORMAT)

# ========================

def read_ascii_file(fpath: str) -> dict:
    """Read data from an ASCII file and return a dictionary. """
    # Placeholder dictionary
    d = {}

    # Read data from file
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for i, line in enumerate(fp.readlines()):
            # Skip the header line(s)
            if "#" in line[0]:
                continue
            data  = reader.read(line)
            Z     = data[0]
            A     = data[1]
            n       = _c.Nuclide(Z=Z, A=A)
            package = {'n': n, 'abundance': data[3], 'uncertainty': data[4]}
            d[n] = package
    return d


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write ASCII file given a filename and a data dictionary. """

    # Write data to file
    with open(fpath, 'w', encoding='utf-8') as fp:
        # Write header line
        fp.write(FILE_HEADER)
        # Loop over data
        for n in data.keys():
            line = data[n].as_tuple
            fp.write(writer.write(line)+"\n")


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the natural abundances database.

    Contains natural isotopic abundances for stable and long-lived nuclides.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        abundance: Natural isotopic abundance [%]
        uncertainty: Uncertainty on abundance [%]
    """

    n          : _c.Nuclide = None
    abundance  : float      = None
    uncertainty: float      = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'abundance': 'Natural isotopic abundance [%]',
        'uncertainty': 'Uncertainty on abundance [%]',
    }


class Database(_db.NuclideDatabase):
    """The natural abundances database. """

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
