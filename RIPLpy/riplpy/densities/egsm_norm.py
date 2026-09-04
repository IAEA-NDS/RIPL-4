
"""
Normalization factors that can be applied to the level density parameter 'a' 
to improve precision of the global systematics. These factors were obtained
for 83 individual elements (8<Z<99) by calculating ratio between experimental 
and systematics' value of 'a' at neutron binding and averaging these ratios 
over the whole isotopic chain. The factors are intended to be used for nuclei 
for which there is no experimental value for S-wave resonance spacing (otherwise 
EGSM code normalizes systematics to reproduce experimental 'a' as given in the 
file 'level-densities-egsm.dat'). 

The average (over Z) value of the normalization factor is 1.00 +- 0.06 (one
standard deviation).
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
FILE_PATH_KEY = 'densities_egsm_norm'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = "# Normalization factors for the level density parameter in the EGSM \n# a_exp/a_sys ratio averaged over isotopes of an element\n#  Z   factor\n#------------------------------------------------------------------------\n"
FORTRAN_FORMAT = '(i5,f8.4)'
reader = _ff.FortranRecordReader(FORTRAN_FORMAT)
writer = _ff.FortranRecordWriter(FORTRAN_FORMAT)

# ========================

def read_ascii_file(fpath: str) -> dict:
    """Read Enhanced Generalized Superfluid Model (EGSM) normalization data from an ASCII formatted file and return as a dictionary. """

    # Placeholder dictionary
    d = {}

    # Read data from file
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            # Normalise line endings (handle CRLF) and skip blanks / comments
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            # Whitespace tokenisation is more robust than the historical
            # fixed-width fortran format (the github release has CRLF endings).
            tokens = stripped.split()
            if len(tokens) < 2:
                continue
            Z = int(tokens[0])
            factor = float(tokens[1])
            ele = _c.Element(Z=Z)
            d[ele] = {'element': ele, 'factor': factor}
    return d


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write a Enhanced Generalized Superfluid Model (EGSM) normalization ASCII file given a filename and a data dictionary. """

    # Write data to file
    with open(fpath, 'w', encoding='utf-8') as fp:
        # Write header line
        fp.write(FILE_HEADER)
        # Loop over data
        for n in data.keys():
            line = data[n].as_tuple
            fp.write(writer.write(line)+"\n")


@_dataclass
class Entry(_db.DatabaseEntry):
    """An entry in the EGSM normalization database.

    Contains normalization factors for the level density parameter 'a'
    to improve global systematics predictions for nuclei without
    experimental resonance data.

    Attributes:
        element: Target element (Element object with Z property)
        factor: Normalization factor (a_exp/a_sys averaged over isotopes)
    """

    element: _c.Element = None
    factor : float      = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'element': 'Target element',
        'factor': 'Normalization factor (a_exp/a_sys averaged over isotopes)',
    }

    @property
    def Z(self, ) -> int:
        return self.element.Z

    @property
    def sym(self, ) -> str:
        return self.element.sym

    @property
    def symbol(self, ) -> str:
        return self.element.symbol

    @property
    def name(self, ) -> str:
        return self.element.name

    @property
    def as_tuple(self, ) -> tuple:
        return tuple([self.Z, self.factor])


class Database(_db.NuclideDatabase):
    """The database object for the Enhanced Generalized Superfluid Model (EGSM) level density model. """

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
