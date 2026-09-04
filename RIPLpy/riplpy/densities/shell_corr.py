# -*- coding: utf-8 -*-
"""Shell correction databases from the Reference Input Parameter Library (RIPL).

Two compilations are provided:

* Myers-Swiatecki (MS) shell corrections — ``densities/shellcor-ms.dat``.
  Available in both the full RIPL distribution and the RIPL-4 GitHub release.
  Loaded via ``load_ms()``.
* Mengoni-Nakajima (MK) shell corrections — ``densities/shellmengoninakajima.dat``.
  RIPL-3 LEGACY: ships with the full RIPL distribution but is NOT part of the
  RIPL-4 GitHub release. Loaded via ``load_mk()``; ``riplpy.densities.load()``
  skips it with a warning when the data file is absent.

Python 3.10+ is expected to run this code properly.

"""

# ========================

# OS
import os as _os

# Dataclasses
from dataclasses import dataclass as _dataclass

# Typing
from typing import ClassVar as _ClassVar

# Fortran format
import fortranformat as _ff

# RIPLpy
import riplpy.db as _db
import riplpy.collections as _c

# ========================

__all__ = ('db_ms_file_path', 'db_mk_file_path', 'FILE_MS_HEADER', 'FORTRAN_READ_MS_FORMAT', 'FORTRAN_WRITE_MS_FORMAT', 'FORTRAN_READ_MK_FORMAT', 'FORTRAN_WRITE_MK_FORMAT', 'reader_ms', 'reader_mk', 'writer_ms', 'writer_mk', 'read_mk_ascii_file', 'read_ms_ascii_file', 'write_ms_ascii_file', 'write_mk_ascii_file', 'Entry', 'ExtendedEntry', 'MS_Database', 'MK_Database', 'load_ms', 'load_mk')

# ========================

# Relative path to database files
db_ms_file_path = _os.path.join("densities", "shellcor-ms.dat")
db_mk_file_path = _os.path.join("densities", "shellmengoninakajima.dat")

# ========================

FILE_MS_HEADER = "#     Myers-Swiatecki shell corrections \n#  Z   A El     Shell   Defcor  bet2  bet4\n#               [MeV]    [MeV]\n#========================================== \n"
FILE_MK_HEADER = ""
FORTRAN_READ_MS_FORMAT  = '(2i4,1x,a2,2(1x,f8.3),2(1x,f7.3))'
FORTRAN_WRITE_MS_FORMAT = '(2i4,1x,a2,2(1x,f8.3),2(1x,f7.3))'
FORTRAN_READ_MK_FORMAT  = '(9x,i3,9x,i3,1x,f8.3)'
FORTRAN_WRITE_MK_FORMAT = '(9x,i3,9x,i3,1x,f8.3)'
reader_ms = _ff.FortranRecordReader(FORTRAN_READ_MS_FORMAT)
writer_ms = _ff.FortranRecordWriter(FORTRAN_WRITE_MS_FORMAT)
reader_mk = _ff.FortranRecordReader(FORTRAN_READ_MK_FORMAT)
writer_mk = _ff.FortranRecordWriter(FORTRAN_WRITE_MK_FORMAT)

# ========================

def read_mk_ascii_file(fpath: str) -> dict:
    """Read Mengonina-Kajima shell correction data from an ASCII formatted file and return as a dictionary. """

    # Placeholder dictionary
    d = {}

    # Read data from file
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for i, line in enumerate(fp.readlines()):
            # Skip the header line(s)
            if "#" in line[0]:
                continue
            data  = reader_mk.read(line)
            Z       = data[0]
            A       = data[1]
            n       = _c.Nuclide(Z=Z, A=A)
            package = {'n': n, 'shell': float(data[2])}
            d[n] = package
    return d


def read_ms_ascii_file(fpath: str) -> dict:
    """Read Myers-Swiatecki shell correction data from an ASCII formatted file and return as a dictionary. """

    # Placeholder dictionary
    d = {}

    # Read data from file
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for i, line in enumerate(fp.readlines()):
            # Skip the header line(s)
            if "#" in line[0]:
                continue
            data  = reader_ms.read(line)
            Z       = data[0]
            A       = data[1]
            n       = _c.Nuclide(Z=Z, A=A)
            package = {'n': n, 'shell': data[3], 'corr': data[4], 'beta2': data[5], 'beta4': data[6]}
            d[n] = package
    return d


def write_ms_ascii_file(fpath: str, data: dict) -> None:
    """Write a shell correction ASCII file given a filename and a data dictionary. """

    # Write data to file
    with open(fpath, 'w', encoding='utf-8') as fp:
        # Write header line
        fp.write(FILE_MS_HEADER)
        # Loop over data
        for n in data.keys():
            line = data[n].as_tuple
            fp.write(writer_ms.write(line)+"\n")


def write_mk_ascii_file(fpath: str, data: dict) -> None:
    """Write a shell correction ASCII file given a filename and a data dictionary. """

    # Write data to file
    with open(fpath, 'w', encoding='utf-8') as fp:
        # Write header line (empty for MK format)
        fp.write(FILE_MK_HEADER)
        # Loop over data - use format matching the reader: (9x,i3,9x,i3,1x,f8.3)
        for n in data.keys():
            entry = data[n]
            # Format: 9 spaces, Z (3 chars), 9 spaces, A (3 chars), 1 space, shell (8.3 float)
            line = f"         {entry.n.Z:3d}         {entry.n.A:3d} {entry.shell:8.3f}"
            fp.write(line+"\n")


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the shell correction database.

    Contains shell correction energies for nuclear level density calculations.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        shell: Shell correction energy [MeV]
    """

    n      : _c.Nuclide = None
    shell  : float      = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'shell': 'Shell correction energy [MeV]',
    }


@_dataclass
class ExtendedEntry(Entry):
    """An extended entry in the Myers-Swiatecki shell correction database.

    Contains shell corrections with deformation parameters.

    Attributes:
        corr: Deformation correction energy [MeV]
        beta2: Quadrupole deformation parameter
        beta4: Hexadecapole deformation parameter
    """

    corr   : float      = None
    beta2  : float      = None
    beta4  : float      = None

    # Field descriptions for help/info display (extends parent)
    _field_info: _ClassVar[dict] = {
        **Entry._field_info,
        'corr': 'Deformation correction energy [MeV]',
        'beta2': 'Quadrupole deformation parameter',
        'beta4': 'Hexadecapole deformation parameter',
    }


class MS_Database(_db.NuclideDatabase):

    reader: object = read_ms_ascii_file
    entry : object = ExtendedEntry
    writer: object = write_ms_ascii_file


class MK_Database(_db.NuclideDatabase):

    reader: object = read_mk_ascii_file
    entry : object = Entry
    writer: object = write_mk_ascii_file


def load_ms(directory: str = None, file_path: str = None) -> MS_Database:
    """Load the database into memory and return it. """
    if file_path is not None:
        return _db.load(MS_Database, file_path)
    elif directory is not None:
        file_path = _os.path.join(directory, db_ms_file_path)
        return _db.load(MS_Database, file_path)


def load_mk(directory: str = None, file_path: str = None) -> MK_Database:
    """Load the database into memory and return it. """
    if file_path is not None:
        return _db.load(MK_Database, file_path)
    elif directory is not None:
        file_path = _os.path.join(directory, db_mk_file_path)
        return _db.load(MK_Database, file_path)
