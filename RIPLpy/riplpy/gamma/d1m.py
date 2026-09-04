# -*- coding: utf-8 -*-
"""D1M+QRPA microscopic E1 / M1 photon strength function tables (RIPL-4).

The RIPL-4 github release ships per-Z files under ``gamma/d1m/`` named
``z<NNN>_e1`` and ``z<NNN>_m1``. Each file contains many nuclei concatenated,
each separated by a header like::

     Z=  8 A= 12 E1 PSF in MeV-3 from D1M+QRPA: 01/01/2019
       E[MeV]   U=000MeV    U=002MeV   ...  U=100MeV
       0.100   3.849E-12   1.944E-10   ...
       ...

For each nucleus we store:

* ``U``     -- the photon-energy grid (column 1)
* ``fE1``   -- the cold-nucleus strength (the ``U=000`` column)
* ``T``     -- the temperature/excitation labels from the per-nucleus header
* ``fE1_T`` -- the full 2-D array of strengths (rows = E, cols = T-labels)

The ``U`` / ``fE1`` schema matches the legacy ``gamma-strength-micro`` GSF
reader so existing user code continues to work.

"""

# OS
import os as _os

# Logging
import logging as _logging

# RIPLpy
import riplpy.collections as _c
import riplpy.db as _db
from riplpy.config import resolve_directory as _resolve_directory

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = ('local_data_dir', 'read_ascii_file', 'Database', 'load',
           'load_all', 'load_element')

local_data_dir = _os.path.join('gamma', 'd1m')


def _parse_zheader(line: str) -> tuple[int, int] | None:
    """Extract (Z, A) from a header beginning with ``Z=`` ``A=``."""
    tokens = line.replace('=', ' = ').split()
    Z = A = None
    i = 0
    while i < len(tokens):
        if tokens[i] == 'Z' and i + 2 < len(tokens):
            try:
                Z = int(tokens[i + 2])
            except ValueError:
                pass
        elif tokens[i] == 'A' and i + 2 < len(tokens):
            try:
                A = int(tokens[i + 2])
            except ValueError:
                pass
        i += 1
    if Z is None or A is None:
        return None
    return Z, A


def _parse_column_labels(line: str) -> list[str]:
    """Parse a column-header like ``E[MeV] U=000MeV U=002MeV ...``."""
    labels = []
    for token in line.split():
        if token.startswith('U=') or token.startswith('T=') or token.startswith('T>'):
            labels.append(token)
    return labels


def read_ascii_file(fpath: str) -> dict:
    """Read a per-Z D1M file containing multiple nuclei.

    Returns ``{Nuclide: package}`` where each ``package`` carries ``U``,
    ``fE1``, ``T``, and ``fE1_T``.
    """
    d: dict = {}
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        lines = fp.readlines()

    n = None
    pkt = None
    labels: list[str] = []
    awaiting_columns = False

    for raw in lines:
        line = raw.rstrip('\n')
        if not line.strip():
            continue
        zaheader = _parse_zheader(line)
        if zaheader is not None:
            # Flush the previous nucleus first
            if n is not None and pkt is not None:
                d[n] = pkt
            Z, A = zaheader
            n = _c.Nuclide(Z=Z, A=A)
            pkt = {'n': n, 'U': [], 'fE1': [], 'T': [], 'fE1_T': []}
            awaiting_columns = True
            continue
        if awaiting_columns:
            # The line immediately after a Z/A header is the column header
            labels = _parse_column_labels(line)
            if pkt is not None:
                pkt['T'] = labels
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
        pkt['fE1'].append(row[0] if row else None)

    # Flush the final nucleus
    if n is not None and pkt is not None:
        d[n] = pkt
    return d


class Entry(_db.PacketEntry):
    """A D1M+QRPA microscopic E1 photon strength function entry."""

    _field_info = {
        'n': 'Target nucleus',
        'U': 'Photon-energy grid [MeV]',
        'fE1': 'Cold-nucleus (T=0) E1 photon strength function [MeV^-3]',
        'fE1_T': 'Per-temperature E1 strength rows (rows = U, cols = T) [MeV^-3]',
        'T': 'Temperature/excitation column labels [MeV]',
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
            _logger.warning(f"D1M data directory not found: {db_loc}")
            return
        for fn in sorted(_os.listdir(db_loc)):
            fpath = _os.path.join(db_loc, fn)
            if _os.path.isfile(fpath):
                try:
                    self.load(fpath)
                except Exception as exc:  # noqa: BLE001
                    _logger.debug(f"Failed to parse D1M file {fpath}: {exc}")


def load(directory: str = None, file_path: str = None) -> Database:
    """Load and return the D1M+QRPA database (whole directory or one file)."""
    db = Database()
    if file_path is not None:
        db.load(file_path)
    else:
        directory = _resolve_directory(directory)
        db.load_all(directory)
    return db


def load_all(directory: str = None) -> Database:
    """Load the whole D1M+QRPA directory."""
    return load(directory=directory)


def load_element(Z: int, directory: str = None,
                  multipolarity: str = 'e1') -> Database:
    """Load D1M data for a single element (E1 or M1)."""
    directory = _resolve_directory(directory)
    db = Database()
    fpath = _os.path.join(directory, local_data_dir,
                          f"z{Z:03d}_{multipolarity.lower()}")
    if _os.path.exists(fpath):
        db.load(fpath)
    else:
        _logger.warning(f"D1M file not found: {fpath}")
    return db
