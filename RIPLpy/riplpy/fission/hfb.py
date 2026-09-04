# -*- coding: utf-8 -*-
"""Python objects which provide access to the HFB fission barrier database of the Reference Input Parameter Library (RIPL).

RIPL-3 LEGACY: The single-file HFB fission barrier table
(``fission/empirical-hfb-barriers.dat``) is a RIPL-3 legacy product. It ships
with the full RIPL distribution but is NOT part of the RIPL-4 GitHub release,
which provides the BSkG3 and D1M barriers instead (see ``fission.bskg3`` and
``fission.d1m``). ``riplpy.fission.load()`` returns an empty database with a
warning when this legacy file is absent.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load the HFB fission barrier database (requires the full RIPL release)::
       $ import riplpy.fission as fission
       $ from riplpy.collections import Nuclide
       $ hfb = fission.hfb.load()
       $ entry = hfb.get(Nuclide(Z=92, A=238))
       $ print(entry.Bin, entry.Bout)

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
FILE_PATH_KEY = 'fission_barriers_hfb'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = "# Z   N   A  s     Bin     hwin   alphain deltain  Bout   hwout  alphaout deltaout Bout2  hwout2 alphaout2 deltaout2\n#                 [MeV]   [MeV] [MeV^-1/2] [MeV]  [MeV]   [MeV] [MeV^-1/2] [MeV]  [MeV]   [MeV] [MeV^-1/2] [MeV]\n#-------------------------------------------------------------------------------------------------------------------\n"
FORTRAN_FORMAT = '(3i4,1x,a2,12f8.2)'
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
            data      = reader.read(line)
            Z         = data[0]
            N         = data[1]
            A         = data[2]
            n         = _c.Nuclide(Z=Z, A=A)
            Bin       = data[4]
            hwin      = data[5]
            alphain   = data[6]
            deltain   = data[7]
            Bout      = data[8]
            hwout     = data[9]
            alphaout  = data[10]
            deltaout  = data[11]
            Bout2     = data[12]
            hwout2    = data[13]
            alphaout2 = data[14]
            deltaout2 = data[15]

            package = {'n': n, 'Bin': Bin, 'hwin': hwin, 'alphain': alphain, 'deltain': deltain,
                       'Bout': Bout, 'hwout': hwout, 'alphaout': alphaout, 'deltaout': deltaout,
                       'Bout2': Bout2, 'hwout2': hwout2, 'alphaout2': alphaout2, 'deltaout2': deltaout2,
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
    """An entry in the HFB fission barrier database.

    Contains microscopic HFB-based fission barrier parameters including
    inner and outer barriers with asymmetry parameters.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        Bin: Inner barrier height [MeV]
        hwin: Inner barrier curvature [MeV]
        alphain: Inner barrier asymmetry [MeV^-1/2]
        deltain: Inner barrier pairing correction [MeV]
        Bout: First outer barrier height [MeV]
        hwout: First outer barrier curvature [MeV]
        alphaout: First outer barrier asymmetry [MeV^-1/2]
        deltaout: First outer barrier pairing correction [MeV]
        Bout2: Second outer barrier height [MeV]
        hwout2: Second outer barrier curvature [MeV]
        alphaout2: Second outer barrier asymmetry [MeV^-1/2]
        deltaout2: Second outer barrier pairing correction [MeV]
    """

    n        : _c.Nuclide = None
    Bin      : float = None
    hwin     : float = None
    alphain  : float = None
    deltain  : float = None
    Bout     : float = None
    hwout    : float = None
    alphaout : float = None
    deltaout : float = None
    Bout2    : float = None
    hwout2   : float = None
    alphaout2: float = None
    deltaout2: float = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'Bin': 'Inner barrier height [MeV]',
        'hwin': 'Inner barrier curvature [MeV]',
        'alphain': 'Inner barrier asymmetry [MeV^-1/2]',
        'deltain': 'Inner barrier pairing correction [MeV]',
        'Bout': 'First outer barrier height [MeV]',
        'hwout': 'First outer barrier curvature [MeV]',
        'alphaout': 'First outer barrier asymmetry [MeV^-1/2]',
        'deltaout': 'First outer barrier pairing correction [MeV]',
        'Bout2': 'Second outer barrier height [MeV]',
        'hwout2': 'Second outer barrier curvature [MeV]',
        'alphaout2': 'Second outer barrier asymmetry [MeV^-1/2]',
        'deltaout2': 'Second outer barrier pairing correction [MeV]',
    }

    @property
    def as_tuple(self) -> tuple:
        """Return a tuple representation matching the ASCII file format."""
        return (self.n.Z, self.n.N, self.n.A, self.n.element_symbol,
                self.Bin, self.hwin, self.alphain, self.deltain,
                self.Bout, self.hwout, self.alphaout, self.deltaout,
                self.Bout2, self.hwout2, self.alphaout2, self.deltaout2)


class Database(_db.Database):

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file

    # @property
    # def heights(self, ) -> list:
    #     """Return the barrier heights in the database as a list of tuples (Z,A,F_A,F_B,F_C). """
    #     ret = []
    #     for nuc in self.data.keys():
    #         Z, A = nuc.Z, nuc.A
    #         ret.append((Z,A,self.data[nuc].A.height,self.data[nuc].B.height,self.data[nuc].C.height))
    #     return ret

    # @property
    # def curvatures(self, ) -> list:
    #     """Return the barrier curvatures in the database as list of tuples (Z,A,C_A,C_B,C_C). """
    #     ret = []
    #     for nuc in self.data.keys():
    #         Z, A = nuc.Z, nuc.A
    #         ret.append((Z,A,self.data[nuc].A.curvature,self.data[nuc].B.curvature,self.data[nuc].C.curvature))
    #     return ret

    # @property
    # def as_nuclides(self, ) -> list:
    #     """Return the database as a list of Nuclides. Each nuclide contains the relevant barrier properties. """
    #     ret = []
    #     for _n, _v in self.data.items():
    #         # Copy nucleus to a new memory location
    #         n = _copy.deepcopy(_n)
    #         # First barrier
    #         n.Bin, n.hwin, n.alphain, n.deltain = _v.A.height, _v.A.curvature, _v.A.alpha, _v.A.delta
    #         # Second barrier
    #         n.Bout, n.hwout, n.alphaout, n.deltaout = _v.B.height, _v.B.curvature, _v.B.alpha, _v.B.delta
    #         # Third barrier
    #         n.Bout2, n.hwout2, n.alphaout2, n.deltaout2 = _v.C.height, _v.C.curvature, _v.C.alpha, _v.C.delta
    #         ret.append(n)
    #     return ret


# Create the local 'load' function with config, database object, and file_path_key pre-filled.
#
# The legacy ``empirical-hfb-barriers.dat`` file is no longer included in the
# RIPL-4 release layout. Wrap the loader so that callers receive an empty
# Database rather than a hard exception when the file is absent.
def _safe_load(directory: str | None = None, file_path: str | None = None) -> 'Database':
    """Load the HFB fission barrier database, returning an empty database if the file is missing."""
    try:
        return _db.loader(
            directory=directory, file_path=file_path,
            config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY,
        )
    except Exception as exc:
        # Importing inside the function avoids any startup-time circular import.
        from riplpy.exceptions import RiplFileNotFoundError
        if isinstance(exc, (FileNotFoundError, RiplFileNotFoundError)):
            try:
                import riplpy as _riplpy
                _riplpy.logger.warning(
                    "fission.hfb: data file not found (%s); using empty database.", exc,
                )
            except Exception:
                pass
            return Database()
        raise


load = _safe_load
