# -*- coding: utf-8 -*-
"""Python objects that provide access to the Constant Temperature (CT) level density model.

RIPL-3 LEGACY: The CT level density parameters
(``densities/level-densities-ctmeff.dat``) are a RIPL-3 legacy product. They
ship with the full RIPL distribution but are NOT part of the RIPL-4 GitHub
release. The reader/writer are retained for users with the full distribution;
``riplpy.densities.load()`` skips this model with a warning when the data file
is absent. For RIPL-4 use the EGSM model (``densities.egsm``) or the
microscopic combinatorial tables.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load and print the CT level density database::
       $ import riplpy.densities as dens
       $ ct = dens.ct.load() # The CT level density database
       $ print(ct.data)
    (2) Access a particular entry of the CT level density database using the get method::
       $ from riplpy.collections import Nuclide
       $ n = Nuclide(z=17, a=36)
       $ print(ct.get(n))
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
import riplpy.db as _db
import riplpy.collections as _c

# Local Namespace
from . import bsfg as _bsfg

# ========================

__all__ = ('FILE_PATH_KEY', 'LOCAL_FILE_PATH', 'FILE_HEADER', 'FORTRAN_FORMAT', 'reader', 'writer', 'read_ascii_file', 'write_ascii_file', 'Entry', 'Database', 'load')

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'densities_ct'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = "   Z  A  El   I0    Bn       D0        Derr    Nlow  Ulow   Ntop  Utop      dW      gamma      ainf     aerr    pairing   Ematch      E0         T\n"
FORTRAN_FORMAT   = '(2i4,1x,a2,1x,f4.1,2x,f6.3,1x,1pe10.3,1x,1pe10.3,0p,1x,i3,2x,f6.3,2x,i3,2x,f6.3,3f10.5,f8.3,4f10.5)'
reader   = _ff.FortranRecordReader(FORTRAN_FORMAT)
writer   = _ff.FortranRecordWriter(FORTRAN_FORMAT)

# ========================

def read_ascii_file(fpath: str) -> dict:
    """Read Constant Temperature (CT) data from an ASCII file and return a dictionary. """
    # Placeholder dictionary
    d = {}

    # Skip this many header lines
    skiprows = 1

    # Read data from file
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for i, line in enumerate(fp.readlines()):
            # Skip the header line(s)
            if i < skiprows:
                continue
            data  = reader.read(line)
            Z       = data[0]
            A       = data[1]
            n       = _c.Nuclide(Z=Z, A=A)
            package = {'n': n, 'Io': data[3], 'Bn': data[4], 'Do': data[5], 'Derr': data[6], 'Nlow': data[7], 'Ulow': data[8],
                       'Ntop': data[9], 'Utop': data[10], 'dW': data[11], 'gamma': data[12], 'ainf': data[13], 'aerr': data[14],
                       'pairing': data[15], 'Ematch': data[16], 'E0': data[17], 'T': data[18],
                      }
            d[n] = package
    return d


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write a Constant Temperature (CT) ASCII file given a filename and a data dictionary. """

    # Write data to file
    with open(fpath, 'w', encoding='utf-8') as fp:
        # Write header line
        fp.write(FILE_HEADER)
        # Loop over data
        for n in data.keys():
            #line = (n.Z, n.A, n.sym, data[n].Io, data[n].Bn, data[n].Do, data[n].Derr, data[n].Nlow, data[n].Ulow, data[n].Ntop, data[n].Utop, data[n].dW, data[n].gamma, data[n].ainf, data[n].aerr, data[n].pairing, data[n].Ematch, data[n].E0, data[n].T)
            line = data[n].as_tuple
            fp.write(writer.write(line)+"\n")


@_dataclass
class Entry(_bsfg.Entry):
    """An entry in the Constant Temperature (CT) database.

    Extends BSFG Entry with additional parameters for the constant
    temperature model at low excitation energies.

    Attributes:
        Ematch: Matching energy between CT and FG formulas [MeV]
        E0: Energy shift for constant temperature formula [MeV]
        T: Nuclear temperature [MeV]
    """

    Ematch: float = None
    E0    : float = None
    T     : float = None

    # Field descriptions for help/info display (extends parent)
    _field_info: _ClassVar[dict] = {
        **_bsfg.Entry._field_info,
        'Ematch': 'Matching energy between CT and FG formulas [MeV]',
        'E0': 'Energy shift for constant temperature formula [MeV]',
        'T': 'Nuclear temperature [MeV]',
    }


class Database(_bsfg.Database):
    """The database object for the Constant Temperature (CT) level density model. """

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
