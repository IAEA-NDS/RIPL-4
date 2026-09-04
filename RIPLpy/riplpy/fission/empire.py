# -*- coding: utf-8 -*-
"""Python objects which provide access to the EMPIRE empirical fission barrier database of the Reference Input Parameter Library (RIPL). 

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load the EMPIRE fission barrier compilation and read an entry::
       $ import riplpy.fission as fission
       $ from riplpy.collections import Nuclide
       $ empire = fission.empire.load()
       $ entry = empire.get(Nuclide(Z=92, A=235))
       $ print(entry.Va, entry.Vb)  # inner / outer barrier heights [MeV]

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
FILE_PATH_KEY = 'fission_barriers_empire'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = "# The tabulation is given for fissioning nuclei (those marked * are 3H)\n#                    All barrier heights, widths and uncertainties are in MeV \n#  Z   A  s   Syma    Va     hwa   Symb    Vb     hwb    Vc  	   Deltaf\n#-------------------------------------------------------------------------------  \n"
FORTRAN_FORMAT = '(A1,I3,1X,I3,1X,A3,2X,A3,3X,F6.2,F6.2,4X,A2,2X,F6.2,2X,F6.2,2X,F6.3,2X,F6.3)'
reader = _ff.FortranRecordReader(FORTRAN_FORMAT)
writer = _ff.FortranRecordWriter(FORTRAN_FORMAT)

# ========================

def _parse_float(token: str | None):
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
    """Read data from the EMPIRE ASCII file and return a dictionary.

    The file contains a mix of whitespace and tab characters used as
    column-alignment fillers. We tokenize on whitespace and treat all
    fields beyond the inner barrier as optional.
    """
    d = {}

    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for raw in fp:
            line = raw.replace('\t', ' ')
            if not line.strip():
                continue
            if line.lstrip().startswith('#'):
                continue

            tokens = line.split()
            # A valid data line begins with Z A symbol (e.g. "80 198 Hg ...").
            if len(tokens) < 5:
                continue
            try:
                Z = int(tokens[0])
                A = int(tokens[1])
            except ValueError:
                continue

            symbol = tokens[2]
            symmetry = tokens[3]
            Va = _parse_float(tokens[4])
            hwa = _parse_float(tokens[5]) if len(tokens) > 5 else None
            Symb = tokens[6] if len(tokens) > 6 else None
            Vb = _parse_float(tokens[7]) if len(tokens) > 7 else None
            hwb = _parse_float(tokens[8]) if len(tokens) > 8 else None
            Vc = _parse_float(tokens[9]) if len(tokens) > 9 else None
            Deltaf = _parse_float(tokens[10]) if len(tokens) > 10 else None

            n = _c.Nuclide(Z=Z, A=A)
            d[n] = {'n': n, 'symmetry': symmetry, 'Va': Va, 'hwa': hwa,
                    'Symb': Symb, 'Vb': Vb, 'hwb': hwb, 'Vc': Vc, 'Deltaf': Deltaf}

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
    """An entry in the EMPIRE empirical fission barrier database.

    Contains fission barrier parameters formatted for use with the EMPIRE code.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        symmetry: Saddle point symmetry designation
        Va: First barrier height [MeV]
        hwa: First barrier curvature [MeV]
        Symb: Second barrier symmetry designation
        Vb: Second barrier height [MeV]
        hwb: Second barrier curvature [MeV]
        Vc: Third barrier height [MeV] (if present)
        Deltaf: Pairing correction [MeV]
    """

    n       : _c.Nuclide = None
    symmetry: str   = None
    Va      : float = None
    hwa     : float = None
    Symb    : str   = None
    Vb      : float = None
    hwb     : float = None
    Vc      : float = None
    Deltaf  : float = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'symmetry': 'Saddle point symmetry designation',
        'Va': 'First barrier height [MeV]',
        'hwa': 'First barrier curvature [MeV]',
        'Symb': 'Second barrier symmetry designation',
        'Vb': 'Second barrier height [MeV]',
        'hwb': 'Second barrier curvature [MeV]',
        'Vc': 'Third barrier height [MeV]',
        'Deltaf': 'Pairing correction [MeV]',
    }

    @property
    def as_tuple(self) -> tuple:
        """Return a tuple representation matching the ASCII file format."""
        marker = '*' if self.Vc is not None else ' '
        return (marker, self.n.Z, self.n.A, self.n.element_symbol, self.symmetry, self.Va, self.hwa,
                self.Symb, self.Vb, self.hwb, self.Vc, self.Deltaf)


class Database(_db.Database):

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
