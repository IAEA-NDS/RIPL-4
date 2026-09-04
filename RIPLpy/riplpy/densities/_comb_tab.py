# -*- coding: utf-8 -*-
"""Shared reader for the combinatorial spin-dependent level density tables
(BSk14, BSkG3, QRPA-BE, T-HFB) provided by S. Goriely in RIPL-4.

These per-Z files (zXXX.tab) follow the same ASCII layout as the HFB
``level-densities-hfb`` tables, so we reuse :mod:`riplpy.densities.hfb`
parsing primitives. Each model just points at a different directory.

"""

# ========================

# OS
import os as _os

# Logging
import logging as _logging

# RIPLpy
import riplpy.db as _db
from riplpy.config import resolve_directory as _resolve_directory
from . import hfb as _hfb

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = ('Database', 'load')


class Database(_hfb.Database):
    """A combinatorial level-density database (per-Z .tab files).

    Subclasses set ``local_data_dir`` to the relative directory containing
    the ``zXXX.tab`` files for the model.
    """

    # Overridden by per-model subclasses
    local_data_dir: str = ""

    def load_all(self, directory: str) -> None:
        """Load all per-Z .tab files from the model directory."""
        db_loc = _os.path.join(directory, self.local_data_dir)
        if not _os.path.exists(db_loc):
            _logger.warning(f"Combinatorial level density directory not found: {db_loc}")
            return
        for fn in sorted(_os.listdir(db_loc)):
            if fn.endswith('.tab'):
                self.load(_os.path.join(db_loc, fn))

    def load_element(self, Z: int, directory: str) -> None:
        """Load the .tab file for a specific element."""
        fpath = _os.path.join(directory, self.local_data_dir, f"z{Z:03d}.tab")
        if _os.path.exists(fpath):
            self.load(fpath)
        else:
            _logger.warning(f"Combinatorial level density file not found: {fpath}")


def load(db_cls, directory: str = None, file_path: str = None, Z: int = None) -> "Database":
    """Generic loader for combinatorial level density databases.

    Args:
        db_cls: Concrete Database subclass for this model.
        directory: Path to the RIPL directory. If provided without Z, loads all elements.
        file_path: Direct path to a specific zXXX.tab file.
        Z: Element number. If provided with directory, loads only that element.

    Returns:
        An instance of ``db_cls`` populated with the loaded data.
    """
    db = db_cls()
    if file_path is not None:
        db.load(file_path)
        return db
    directory = _resolve_directory(directory)
    if Z is not None:
        db.load_element(Z, directory)
    else:
        db.load_all(directory)
    return db
