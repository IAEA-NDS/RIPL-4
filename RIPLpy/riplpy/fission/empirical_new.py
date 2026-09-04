# -*- coding: utf-8 -*-
"""Python objects which provide access to the RIPL-4 empirical fission barrier database.

This module is an alias of :mod:`riplpy.fission.empirical` retained for
backward-compatibility. In RIPL-3 the empirical barriers were split across two
separate files (``empirical-barriers.dat`` and ``empirical-barriers-new.dat``).
RIPL-4 consolidates both into ``empirical-barriers-ripl4.dat``, so both
``fission.empirical`` and ``fission.empirical_new`` now read the same file.

Python 3.10+ is expected to run this code properly.

"""

# ========================

# Functools
from functools import partial as _partial

# RIPLpy
import riplpy.config as _config
import riplpy.db as _db

# Re-export the new RIPL-4 reader/writer/entry from the empirical module
from .empirical import (  # noqa: F401
    read_ascii_file,
    write_ascii_file,
    Entry,
    Database as _BaseDatabase,
    FILE_HEADER,
)

# ========================

__all__ = ('FILE_PATH_KEY', 'LOCAL_FILE_PATH', 'FILE_HEADER', 'read_ascii_file', 'write_ascii_file', 'Entry', 'Database', 'load')

# ========================

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'fission_barriers_empirical_new'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================


class Database(_BaseDatabase):
    """Alias of :class:`riplpy.fission.empirical.Database` for backward-compatibility."""


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
