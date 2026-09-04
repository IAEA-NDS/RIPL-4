# -*- coding: utf-8 -*-
"""Python objects that provide access to the Back-Shifted Fermi-Gas (BSFG) level density model.

   RIPL-3 LEGACY: The BSFG level density parameters
   (``densities/level-densities-bfmeff.dat``) are a RIPL-3 legacy product.
   They ship with the full RIPL distribution but are NOT part of the RIPL-4
   GitHub release. The reader/writer are retained for users with the full
   distribution; ``riplpy.densities.load()`` skips this model with a warning
   when the data file is absent. For RIPL-4 use the EGSM model
   (``densities.egsm``) or the microscopic combinatorial tables.

   Level density parameters are provided for the BSFG model without explicit treatment of collective effects.
   The parameters are obtained by fitting the Fermi-gas model formula both to the RIPL II recommended spacings of s-wave neutron resonances D0 and to the cumulative number of low-lying levels evaluated from the analysis of nuclear levels. 
   Since the goal is to reproduce as well as possible both the discrete levels and the D0 values, the adopted parameters do not generally provide a perfect fit of D0. 
   Instead, we obtain a theoretical value Dth. 
   Therefore, uncertainty aerr is deduced by fitting Dth+Derr or Dth-Derr (See TECHDOC for more explanations). 

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load and print the BSFG level density database::
       $ import riplpy.densities as dens
       $ bsfg = dens.bsfg.load() # The BSFG level density database
       $ print(bsfg.data)
    (2) Access a particular entry of the BSFG level density database using the get method::
       $ from riplpy.collections import Nuclide
       $ n = Nuclide(z=50, a=125)
       $ print(bsfg.get(n))
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

# ========================

__all__ = ('FILE_PATH_KEY', 'LOCAL_FILE_PATH', 'FILE_HEADER', 'FORTRAN_FORMAT', 'reader', 'writer', 'read_ascii_file', 'write_ascii_file', 'Entry', 'Database', 'load')

# ========================

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'densities_bsfg'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = "   Z  A  El   I0    Bn       D0        Derr    Nlow  Ulow   Ntop  Utop      dW      gamma      ainf     aerr    pairing\n"
FORTRAN_FORMAT = '(2i4,1x,a2,1x,f4.1,2x,f6.3,1x,1pe10.3,1x,1pe10.3,0p,1x,i3,2x,f6.3,2x,i3,2x,f6.3,3f10.5,f8.3,1f10.5)'
reader = _ff.FortranRecordReader(FORTRAN_FORMAT)
writer = _ff.FortranRecordWriter(FORTRAN_FORMAT)

# ========================

def read_ascii_file(fpath: str) -> dict:
    """Read BSFG data from an ASCII formatted file and return as a dictionary. """

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
                       'pairing': data[15],
                      }
            d[n] = package
    return d


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write a Back-Shifted Fermi-Gas (BSFG) ASCII file given a filename and a data dictionary. """

    # Write data to file
    with open(fpath, 'w', encoding='utf-8') as fp:
        # Write header line
        fp.write(FILE_HEADER)
        # Loop over data
        for n in data.keys():
            #line = (n.Z, n.A, n.sym, data[n].Io, data[n].Bn, data[n].Do, data[n].Derr, data[n].Nlow, data[n].Ulow, data[n].Ntop, data[n].Utop, data[n].dW, data[n].gamma, data[n].ainf, data[n].aerr, data[n].pairing)
            line = data[n].as_tuple
            fp.write(writer.write(line)+"\n")


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the Back-Shifted Fermi-Gas (BSFG) database.

    Contains level density parameters for the BSFG model fitted to
    s-wave neutron resonance spacings and discrete level data.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        Io: Ground state spin of target nucleus
        Bn: Neutron binding energy [MeV]
        Do: Average s-wave resonance spacing [eV]
        Derr: Uncertainty on resonance spacing [eV]
        Nlow: Lowest discrete level used in fit
        Ulow: Excitation energy of Nlow [MeV]
        Ntop: Highest discrete level used in fit
        Utop: Excitation energy of Ntop [MeV]
        dW: Shell correction energy [MeV]
        gamma: Ignatyuk damping parameter
        ainf: Asymptotic level density parameter [MeV^-1]
        aerr: Uncertainty on ainf [MeV^-1]
        pairing: Pairing energy shift [MeV]
    """

    n      : _c.Nuclide = None
    Io     : float      = None
    Bn     : float      = None
    Do     : float      = None
    Derr   : float      = None
    Nlow   : int        = None
    Ulow   : float      = None
    Ntop   : int        = None
    Utop   : float      = None
    dW     : float      = None
    gamma  : float      = None
    ainf   : float      = None
    aerr   : float      = None
    pairing: float      = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'Io': 'Ground state spin of target nucleus',
        'Bn': 'Neutron binding energy [MeV]',
        'Do': 'Average s-wave resonance spacing [eV]',
        'Derr': 'Uncertainty on resonance spacing [eV]',
        'Nlow': 'Lowest discrete level used in fit',
        'Ulow': 'Excitation energy of Nlow [MeV]',
        'Ntop': 'Highest discrete level used in fit',
        'Utop': 'Excitation energy of Ntop [MeV]',
        'dW': 'Shell correction energy [MeV]',
        'gamma': 'Ignatyuk damping parameter',
        'ainf': 'Asymptotic level density parameter [MeV^-1]',
        'aerr': 'Uncertainty on ainf [MeV^-1]',
        'pairing': 'Pairing energy shift [MeV]',
    }


class Database(_db.NuclideDatabase):
    """The database object for the Back-Shifted Fermi-Gas (BSFG) level density model. """

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
