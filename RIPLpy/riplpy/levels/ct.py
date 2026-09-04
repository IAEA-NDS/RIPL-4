

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
FILE_PATH_KEY = 'levels_param'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER    = "#Constant Temperature fit results; RIPL III discrete level library; December 2021\n#  Z   A El         T        dT        U0       dU0NlevNmax  N0  Nc     Umax       Uc     Chi Fit Flag NoX  Xm  EX   sigma\n#               [MeV]     [MeV]     [MeV]     [MeV]                    [MeV]    [MeV]                         [MeV]\n#-------------------------------------------------------------------------------------------------------------------------\n"
FORTRAN_FORMAT = '(2i4,1x,a2,4(1x,f9.5),4i4,1x,f8.5,1x,f8.5,(1x,e10.3),1x,A1,1x,A1,2i4,f7.4,1x,f6.3)'
reader = _ff.FortranRecordReader(FORTRAN_FORMAT)
writer = _ff.FortranRecordWriter(FORTRAN_FORMAT)

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
            package = {'n': n,            # Nucleus
                       'T': data[3],      # Temperature T in the CT model
                       'dT': data[4],     # Uncertainty of T
                       'U0': data[5],     # Back-shift in CT model
                       'dU0': data[6],    # Uncertainty of U0
                       'Nlev': data[7],   # Number of levels in the ENSDF data set
                       'Nmax': data[8],   # Maximum level number up to which the level scheme is complete
                       'N0': data[9],     # Minimum level number from which the fit is considered
                       'Nc': data[10],    # Number of level at which unique spin sequence ends
                       'Umax': data[11],  # Energy corresponds to Nmax
                       'Uc': data[12],    # Energy corresponds to Nc
                       'Chi': data[13],   # Measure of the quality of fit (blank if there was no fit)
                       'Fit': data[14],   # If the record comes from the fit of pre-selected nuclei, which provided the T(A) function = blank if there was no fit
                       'Flag': data[15],  # 'F' if Chi>0.05 (bad fit); blank otherwise
                       'NoX': data[16],   # Number of levels with +X, +Y, +Z ... notation (X,Y,Z... are unknown energy values)
                       'Xm': data[17],    # Level number where the first +X, +Y, +Z ... notation appears
                       'Ex': data[18],    # Level energy where the first +X, +Y, +Z ... notation appears
                       'sigma': data[19]  # Spin cut-off values extracted from the discrete level library
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
    """An entry in the Constant Temperature (CT) fit database.

    Contains CT model parameters fitted to discrete level data for
    extrapolating level densities to higher excitation energies.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        T: Nuclear temperature [MeV]
        dT: Uncertainty on T [MeV]
        U0: Back-shift energy [MeV]
        dU0: Uncertainty on U0 [MeV]
        Nlev: Total number of ENSDF levels
        Nmax: Completeness limit (level number)
        N0: First level included in fit
        Nc: Level where unique spin sequence ends
        Umax: Energy at completeness limit [MeV]
        Uc: Energy at Nc [MeV]
        Chi: Fit quality chi-squared
        Fit: Fit flag (blank if no fit)
        Flag: 'F' if poor fit (Chi>0.05)
        NoX: Number of levels with +X notation
        Xm: First level with +X notation
        Ex: Energy of first +X level [MeV]
        sigma: Spin cutoff parameter
    """

    n    : _c.Nuclide = None
    T    : float      = None
    dT   : float      = None
    U0   : float      = None
    dU0  : float      = None
    Nlev : int        = None
    Nmax : int        = None
    N0   : int        = None
    Nc   : int        = None
    Umax : float      = None
    Uc   : float      = None
    Chi  : float      = None
    Fit  : str        = None
    Flag : str        = None
    NoX  : int        = None
    Xm   : int        = None
    Ex   : float      = None
    sigma: float      = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'T': 'Nuclear temperature [MeV]',
        'dT': 'Uncertainty on T [MeV]',
        'U0': 'Back-shift energy [MeV]',
        'dU0': 'Uncertainty on U0 [MeV]',
        'Nlev': 'Total number of ENSDF levels',
        'Nmax': 'Completeness limit (level number)',
        'N0': 'First level included in fit',
        'Nc': 'Level where unique spin sequence ends',
        'Umax': 'Energy at completeness limit [MeV]',
        'Uc': 'Energy at Nc [MeV]',
        'Chi': 'Fit quality chi-squared',
        'Fit': 'Fit flag (blank if no fit)',
        'Flag': 'F if poor fit (Chi>0.05)',
        'NoX': 'Number of levels with +X notation',
        'Xm': 'First level with +X notation',
        'Ex': 'Energy of first +X level [MeV]',
        'sigma': 'Spin cutoff parameter',
    }


class Database(_db.NuclideDatabase):
    """Database of Constant Temperature (CT) fit of nuclear level schemes. """

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
