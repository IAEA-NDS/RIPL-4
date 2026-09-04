# -*- coding: utf-8 -*-
"""TLO (Triple-Lorentzian Open) E1 photon strength function tables (RIPL-4).

The github release contains per-Z files under ``gamma/tlo/`` named
``FE1_z<NN>.dat``. Each file packs many nuclei separated by headers of the
form::

     Z= 34 A= 64 bet= 0.192 gam= 29.000   EFF
      U[MeV]  fE1[MeV**-3]
       0.100  0.000  0.264E-09  0.000E+00
       ...

For each nucleus we store ``U``, ``fE1`` (the second strength column, the
unconditional fE1 value), the deformation parameters ``beta`` / ``gamma``,
the calculation ``mode`` flag, and the full per-nucleus ``fE1_T`` array.

"""

# OS
import os as _os

# Logging
import logging as _logging

# Re
import re as _re

# RIPLpy
import riplpy.collections as _c
import riplpy.db as _db
from riplpy.config import resolve_directory as _resolve_directory

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = ('local_data_dir', 'read_ascii_file', 'Database', 'load',
           'load_all', 'load_element')

local_data_dir = _os.path.join('gamma', 'tlo')

_FNAME_RE = _re.compile(r"FE1_z(\d+)", _re.IGNORECASE)


def _parse_header(line: str) -> dict | None:
    """Parse a header line such as ``Z= 34 A= 64 bet= 0.192 gam= 29.000   EFF``."""
    s = line.replace('=', ' = ')
    tokens = s.split()
    info: dict = {}
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in ('Z', 'A') and i + 2 < len(tokens) and tokens[i + 1] == '=':
            try:
                info[tok] = int(tokens[i + 2])
            except ValueError:
                pass
            i += 3
            continue
        if tok in ('bet', 'gam') and i + 2 < len(tokens) and tokens[i + 1] == '=':
            try:
                info[tok] = float(tokens[i + 2])
            except ValueError:
                pass
            i += 3
            continue
        if tok in ('EFF', 'INTER'):
            info['mode'] = tok
        i += 1
    if 'Z' not in info or 'A' not in info:
        return None
    return info


def read_ascii_file(fpath: str) -> dict:
    """Read a per-Z TLO file containing multiple nuclei."""
    d: dict = {}
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        # The TLO files use CRLF line endings; readlines handles that.
        lines = fp.readlines()

    n = None
    pkt = None
    awaiting_columns = False
    for raw in lines:
        line = raw.rstrip('\r\n')
        if not line.strip():
            continue
        header = _parse_header(line) if line.lstrip().startswith('Z=') or 'Z=' in line[:6] else None
        if header is not None:
            if n is not None and pkt is not None:
                d[n] = pkt
            Z = header['Z']
            A = header['A']
            n = _c.Nuclide(Z=Z, A=A)
            pkt = {
                'n': n, 'U': [], 'fE1': [], 'fE1_T': [],
                'beta':  header.get('bet'),
                'gamma': header.get('gam'),
                'mode':  header.get('mode'),
            }
            awaiting_columns = True
            continue
        if awaiting_columns:
            # The line right after the Z/A header is the column header text
            awaiting_columns = False
            continue
        if pkt is None:
            continue
        parts = line.split()
        try:
            energy = float(parts[0])
        except (ValueError, IndexError):
            continue
        try:
            row = [float(x) for x in parts[1:]]
        except ValueError:
            continue
        pkt['U'].append(energy)
        pkt['fE1_T'].append(row)
        # The 2nd column (parts[2]) is the "fE1" strength; first column is a
        # zero/auxiliary value
        if len(row) >= 2:
            pkt['fE1'].append(row[1])
        elif row:
            pkt['fE1'].append(row[0])
        else:
            pkt['fE1'].append(None)

    if n is not None and pkt is not None:
        d[n] = pkt
    return d


class Entry(_db.PacketEntry):
    """A TLO (Triple-Lorentzian Open) E1 photon strength function entry."""

    _field_info = {
        'n': 'Target nucleus',
        'U': 'Photon-energy grid [MeV]',
        'fE1': 'E1 dipole strength function [MeV^-3]',
        'fE1_T': 'Full per-nucleus strength rows (rows = U) [MeV^-3]',
        'beta': 'Ground-state beta deformation parameter (Delaroche et al.)',
        'gamma': 'Ground-state gamma deformation parameter [deg]',
        'mode': 'Deformation mode flag (EFF = reduced near shells, INTER = '
                'interpolated odd nuclei, 5DCH = Delaroche et al.)',
    }


class Database(_db.NuclideDatabase):

    reader: object = read_ascii_file
    entry : object = Entry
    local_data_dir: str = local_data_dir

    def load(self, fpath: str) -> None:
        data = type(self).reader(fpath)
        entry_cls = type(self).entry
        for n, pkt in data.items():
            self.data[n] = entry_cls(pkt)

    def load_all(self, directory: str) -> None:
        db_loc = _os.path.join(directory, self.local_data_dir)
        if not _os.path.isdir(db_loc):
            _logger.warning(f"TLO data directory not found: {db_loc}")
            return
        for fn in sorted(_os.listdir(db_loc)):
            if not _FNAME_RE.search(fn):
                continue
            fpath = _os.path.join(db_loc, fn)
            if _os.path.isfile(fpath):
                try:
                    self.load(fpath)
                except Exception as exc:  # noqa: BLE001
                    _logger.debug(f"Failed to parse TLO file {fpath}: {exc}")


def load(directory: str = None, file_path: str = None) -> Database:
    db = Database()
    if file_path is not None:
        db.load(file_path)
    else:
        directory = _resolve_directory(directory)
        db.load_all(directory)
    return db


def load_all(directory: str = None) -> Database:
    return load(directory=directory)


def load_element(Z: int, directory: str = None) -> Database:
    directory = _resolve_directory(directory)
    db = Database()
    fpath = _os.path.join(directory, local_data_dir, f"FE1_z{Z:03d}.dat")
    if _os.path.exists(fpath):
        db.load(fpath)
    else:
        _logger.warning(f"TLO file not found: {fpath}")
    return db
