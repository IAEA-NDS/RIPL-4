# -*- coding: utf-8 -*-
"""SMLO M1 photon strength function tables (RIPL-4).

Per-Z files under ``gamma/smlo_M1/`` (e.g. ``z026_m1``). Each file packs many
nuclei with the same column layout as the D1M+QRPA tables, so we reuse the
multi-nucleus reader from :mod:`riplpy.gamma.d1m`.

Each per-nucleus block uses two strength columns labelled ``T=0`` and ``T>0``
(cold versus finite-temperature limit), so ``fE1`` (alias of ``fM1``) is just
the ``T=0`` column.

"""

# OS
import os as _os

# Logging
import logging as _logging

# RIPLpy
import riplpy.db as _db
from riplpy.config import resolve_directory as _resolve_directory
from . import d1m as _d1m

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = ('local_data_dir', 'Database', 'load', 'load_all', 'load_element')

local_data_dir = _os.path.join('gamma', 'smlo_M1')


class Entry(_db.PacketEntry):
    """A SMLO M1 photon strength function entry.

    Reuses the D1M+QRPA multi-nucleus column layout; the two strength
    columns are the cold (``T=0``) and finite-temperature (``T>0``) limits.
    """

    _field_info = {
        'n': 'Target nucleus',
        'U': 'Photon-energy grid [MeV]',
        'fE1': 'Cold-nucleus (T=0) M1 photon strength function [MeV^-3]',
        'fE1_T': 'M1 strength rows (rows = U, cols = T=0 / T>0) [MeV^-3]',
        'T': 'Strength column labels (T=0, T>0)',
    }


class Database(_d1m.Database):

    reader: object = _d1m.read_ascii_file
    entry : object = Entry
    local_data_dir: str = local_data_dir

    def load_all(self, directory: str) -> None:
        db_loc = _os.path.join(directory, self.local_data_dir)
        if not _os.path.isdir(db_loc):
            _logger.warning(f"SMLO M1 data directory not found: {db_loc}")
            return
        for fn in sorted(_os.listdir(db_loc)):
            fpath = _os.path.join(db_loc, fn)
            if _os.path.isfile(fpath):
                try:
                    self.load(fpath)
                except Exception as exc:  # noqa: BLE001
                    _logger.debug(f"Failed to parse SMLO M1 file {fpath}: {exc}")


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
    fpath = _os.path.join(directory, local_data_dir, f"z{Z:03d}_m1")
    if _os.path.exists(fpath):
        db.load(fpath)
    else:
        _logger.warning(f"SMLO M1 file not found: {fpath}")
    return db
