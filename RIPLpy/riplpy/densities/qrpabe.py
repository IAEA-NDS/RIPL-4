# -*- coding: utf-8 -*-
"""QRPA-BE spin-dependent nuclear level densities.

Data files (per-Z ``zXXX.tab``) provided by S. Goriely under
``densities/total/qrpabe/`` in the RIPL-4 release. The file format
matches the HFB level density tables, so parsing is reused.

"""

# OS
import os as _os

# RIPLpy
from . import _comb_tab as _comb_tab
from . import hfb as _hfb

# ========================

__all__ = ('local_data_dir', 'Database', 'load', 'load_all', 'load_element')

# Relative path to database folder
local_data_dir = _os.path.join("densities", "total", "qrpabe")


class Database(_comb_tab.Database):
    """QRPA-BE level density database (per-Z .tab files)."""

    reader: object = _hfb.read_ascii_file
    entry: object = _hfb.Entry
    writer: object = _hfb.write_ascii_file
    local_data_dir: str = local_data_dir


def load(directory: str = None, file_path: str = None, Z: int = None) -> Database:
    """Load and return the QRPA-BE level density database."""
    return _comb_tab.load(Database, directory=directory, file_path=file_path, Z=Z)


def load_all(directory: str = None) -> Database:
    """Load all QRPA-BE level densities under ``directory``."""
    return load(directory=directory)


def load_element(Z: int, directory: str = None) -> Database:
    """Load QRPA-BE level densities for a single element."""
    return load(directory=directory, Z=Z)
