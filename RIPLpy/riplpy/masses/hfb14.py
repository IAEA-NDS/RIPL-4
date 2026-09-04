# -*- coding: utf-8 -*-
"""Python objects that provide access to the HFB-14 mass model.

RIPL-3 LEGACY: The HFB-14 mass table (``masses/mass-hfb14.dat``) is a RIPL-3
legacy product. It ships with the full RIPL distribution but is NOT part of
the RIPL-4 GitHub release. The reader/writer are retained for users with the
full distribution; ``riplpy.masses.load()`` skips this model with a warning
when the data file is absent. For RIPL-4 use ``masses.hfb27`` (HFB-27)
instead.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load and print the HFB14 database (requires the full RIPL release)::
       $ import riplpy.masses as masses
       $ hfb14 = masses.hfb14.load() # The HFB-14 database
       $ print(hfb14.data)
    (2) Access a particular entry of the mass database using the get method::
       $ from riplpy.collections import Nuclide
       $ n = Nuclide(z=82, a=208)
       $ print(hfb14.get(n))
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
import riplpy.exceptions as _exceptions
import riplpy.collections as _c
import riplpy.db as _db
from . import core as _core

# ========================

__all__ = ('FILE_PATH_KEY', 'LOCAL_FILE_PATH', 'FILE_HEADER', 'FORTRAN_FORMAT', 'reader', 'writer', 'read_ascii_file', 'write_ascii_file', 'Entry', 'Database', 'load')

# ========================

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'mass_hfb14'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = "#\n#  Z   A s fl     Mexp      Err       Mth    beta2   beta4     rhon     rn       an       rhop     rp       ap   \n#                [MeV]     [MeV]     [MeV]                    [fm-3]   [fm]     [fm]     [fm-3]   [fm]     [fm]  \n#-----------------------------------------------------------------------------------------------------------------\n"
FORTRAN_FORMAT   = '(2i4,1x,a2,1x,i1,3f10.3,2f8.3,6f9.4)'
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
            package = {'n': n, 'flag': data[3], 'Mexp': data[4], 'Err': data[5], 'Mth': data[6], 
                       'beta2': data[7], 'beta4': data[8], 'rhon': data[9], 'rn': data[10],
                       'an': data[11], 'rhop': data[12], 'rp': data[13], 'ap': data[14],
                      }
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
    """An entry in the HFB-14 mass database.

    Contains theoretical atomic masses from the HFB-14 model along with
    experimental reference values and nuclear structure parameters.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        flag: Data quality flag (0=no exp, 1=recommended, 2=measured)
        Mexp: Experimental mass excess [MeV]
        Err: Uncertainty on mass excess [MeV]
        Mth: Calculated HFB-14 mass excess [MeV]
        beta2: Quadrupole deformation
        beta4: Hexadecapole deformation
        rhon: Neutron density amplitude [fm^-3]
        rn: Neutron density radius [fm]
        an: Neutron density diffuseness [fm]
        rhop: Proton density amplitude [fm^-3]
        rp: Proton density radius [fm]
        ap: Proton density diffuseness [fm]
    """

    n     : _c.Nuclide = None
    flag  : int        = None
    Mexp  : float      = None
    Err   : float      = None
    Mth   : float      = None
    beta2 : float      = None
    beta4 : float      = None
    rhon  : float      = None
    rn    : float      = None
    an    : float      = None
    rhop  : float      = None
    rp    : float      = None
    ap    : float      = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'flag': 'Data quality flag (0=no exp, 1=recommended, 2=measured)',
        'Mexp': 'Experimental mass excess [MeV]',
        'Err': 'Uncertainty on mass excess [MeV]',
        'Mth': 'Calculated HFB-14 mass excess [MeV]',
        'beta2': 'Quadrupole deformation',
        'beta4': 'Hexadecapole deformation',
        'rhon': 'Neutron density amplitude [fm^-3]',
        'rn': 'Neutron density radius [fm]',
        'an': 'Neutron density diffuseness [fm]',
        'rhop': 'Proton density amplitude [fm^-3]',
        'rp': 'Proton density radius [fm]',
        'ap': 'Proton density diffuseness [fm]',
    }


class Database(_core.MassDatabase):

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
