# -*- coding: utf-8 -*-
"""Python objects which provide access to the theoretical giant dipole resonance (GDR) parameters of the Reference Input Parameter Library (RIPL).

RIPL-3 LEGACY: The single-file theoretical GDR parameter table
(``gamma/gdr-parameters-theor.dat``) is a RIPL-3 legacy product. It ships
with the full RIPL distribution but is NOT part of the RIPL-4 GitHub release,
which replaces it with the per-nucleus D1M+QRPA predictions under
``gamma/d1m/`` (see ``riplpy.gamma.d1m``). ``riplpy.gamma.load()`` returns an
empty database with a warning when this legacy file is absent.

Predictions of the Giant Dipole Resonance (GDR) energies and widths for about 6000 nuclei with 14<=Z<=110 lying between the proton and the neutron driplines.
The GDR is represented in the Goldhaber-Teller model [1] where the neutron and proton densities perform an out-of-phase vibration around their centre of mass. 
The dynamics of the oscillation is assumed to be dominated by the np-interaction [2]. 
The present table gives the GDR energies predicted by [2] with a renormalized np-interaction of strength derived from a least-square fit to the experimental GDR energies [3]. 
The nucleon density distribution and ground-state deformation are taken from the Extended Thomas-Fermi plus Strutinsky Integral (ETFSI) compilation [4,5]. 
The expression for the shell-dependent GDR width is taken from [6] using the newly-determined GDR energies and the ETFSI shell corrections. 
Such predictions include the shell-dependent GDR broadening due to the coupling between the dipole oscillations and the quadrupole surface vibrations. 
Comparison between predicted and experimental GDR energies and widths can be found in [3]. 
In case of deformed nuclei, the GDR splits into two peaks for oscillations parallel to the axis of rotational symmetry and perpendicular to it.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load the legacy theoretical GDR database (full RIPL release only;
        returns an empty database on the GitHub release)::
       $ import riplpy.gamma as gamma
       $ from riplpy.collections import Nuclide
       $ gdr = gamma.gdr.load()
       $ entry = gdr.get(Nuclide(Z=82, A=208))
       $ print(entry.E1, entry.W1)  # GDR energy / width [MeV]

"""

# ========================

# Functools
from functools import partial as _partial

# OS
import os as _os

# Logging
import logging as _logging

# Fortran format
import fortranformat as _ff

# RIPLpy
import riplpy.config as _config
import riplpy.collections as _c
import riplpy.db as _db
from riplpy.config import resolve_directory as _resolve_directory
from riplpy.exceptions import RiplFileNotFoundError as _RiplFileNotFoundError
from . import core as _core

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = ('FILE_PATH_KEY', 'LOCAL_FILE_PATH', 'FILE_HEADER', 'FORTRAN_FORMAT', 'reader', 'writer', 'read_ascii_file', 'write_ascii_file', 'Entry', 'Database', 'load')

# ========================

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'gdr_parameters_theor'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = "#  Theoretical Giant Dipole Resonance parameters\n#  Z   A El    Eta    E1     W1     E2     W2\n#                   [MeV]  [MeV]  [MeV]  [MeV] \n#--------------------------------------------------\n"
FORTRAN_FORMAT = '(2i4,1x,a2,f7.3,4f7.2)'
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
            package = {'n': n, 'eta': data[3], 
                       'E1': data[4], 'W1': data[5], 
                       'E2': data[6], 'W2': data[7],
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


class Entry(_core.Theory_GDR_Parameter, _db.NuclideDatabaseEntry):
    """An object which represents an entry in the theoretical GDR parameter database. """

    @property
    def as_tuple(self) -> tuple:
        """Return a tuple representation matching the ASCII file format: Z, A, El, eta, E1, W1, E2, W2"""
        return (self.n.Z, self.n.A, self.n.element_symbol, self.eta, self.E1, self.W1, self.E2, self.W2)


class Database(_db.NuclideDatabase):

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


def load(directory: str = None, file_path: str = None) -> Database:
    """Load the theoretical GDR parameter database.

    The legacy ``gdr-parameters-theor.dat`` file is NOT shipped in the
    RIPL-4 github release. When the file is absent this loader returns
    an empty :class:`Database` (and logs a warning) rather than raising,
    so that the section ``load()`` driver keeps working.
    """
    # Resolve target path
    db_path = file_path
    if db_path is None:
        # Fall back to the configured RIPL path when no directory is given.
        # _resolve_directory raises only when nothing is configured at all;
        # a configured-but-missing legacy file still yields an empty database
        # with a warning (handled below).
        directory = _resolve_directory(directory)
        rel = _config.DATA_FILES.get(FILE_PATH_KEY)
        db_path = _os.path.join(directory, rel) if rel else None

    if db_path is None or not _os.path.exists(db_path):
        _logger.warning(
            "Theoretical GDR file 'gdr-parameters-theor.dat' is not "
            "available; returning empty database. (RIPL-4 github release "
            "replaces this product with the per-nucleus D1M+QRPA tables under "
            "'gamma/d1m/'.)"
        )
        return Database()

    return _db.load(Database, db_path)
