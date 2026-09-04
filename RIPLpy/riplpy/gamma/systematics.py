# -*- coding: utf-8 -*-
"""SLO / SMLO experimental + systematics GDR parameter compilations (RIPL-4).

These files (``gamma/gdr_parameters_exp&systematics/``) extend the experimental
SLO and SMLO fits with global systematics for ~8980 isotopes. Each record is a
single line with the Fortran layout ``(2I4, 9F9.3, I5)``::

    Z  A  Er1  Wr1  S1  Er2  Wr2  S2  S  CSp1  CSp2  In

where ``In=1`` flags experimental data and ``In=0`` flags systematics-only
values. Spherical-nucleus rows leave the deformed-component fields blank.

"""

# OS
import os as _os

# Logging
import logging as _logging

# RIPLpy
import riplpy.collections as _c
import riplpy.db as _db
from riplpy.config import resolve_directory as _resolve_directory
from . import core as _core

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = (
    'slo_local_file_path', 'smlo_local_file_path',
    'read_ascii_file', 'Entry', 'Database', 'load_slo', 'load_smlo',
)

# Relative paths
slo_local_file_path  = _os.path.join(
    'gamma', 'gdr_parameters_exp&systematics', 'gdr-parameters_exp&systematics_slo.dat'
)
smlo_local_file_path = _os.path.join(
    'gamma', 'gdr_parameters_exp&systematics', 'gdr-parameters_exp&systematics_smlo.dat'
)


def _safe_float(s: str) -> float | None:
    s = s.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def read_ascii_file(fpath: str) -> dict:
    """Read a ``gdr-parameters_exp&systematics_<slo|smlo>.dat`` file."""
    d = {}
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            try:
                Z = int(line[0:4])
                A = int(line[4:8])
                # Nine 9-wide float fields
                Er1  = _safe_float(line[8:17])
                Wr1  = _safe_float(line[17:26])
                S1   = _safe_float(line[26:35])
                Er2  = _safe_float(line[35:44])
                Wr2  = _safe_float(line[44:53])
                S2   = _safe_float(line[53:62])
                S    = _safe_float(line[62:71])
                CSp1 = _safe_float(line[71:80])
                CSp2 = _safe_float(line[80:89])
                In   = int(line[89:94]) if line[89:94].strip() else None
            except (ValueError, IndexError) as exc:
                _logger.debug(f"Skipping malformed line in {fpath}: {exc}")
                continue
            n = _c.Nuclide(Z=Z, A=A)
            d[n] = {
                'n': n,
                'E1': Er1, 'W1': Wr1,
                'E2': Er2, 'W2': Wr2,
                'Sr1': S1, 'CSp1': CSp1,
                'Sr2': S2, 'CSp2': CSp2,
                'Sr':  S,
                'In':  In,
            }
    return d


class Entry(_core.Systematics_GDR_Parameter, _db.NuclideDatabaseEntry):
    """An experimental-or-systematics GDR parameter entry."""
    pass


class Database(_db.NuclideDatabase):

    reader: object = read_ascii_file
    entry : object = Entry


def _load(rel_path: str, directory: str = None, file_path: str = None) -> Database:
    db = Database()
    db_path = file_path
    if db_path is None:
        directory = _resolve_directory(directory)
        db_path = _os.path.join(directory, rel_path)
    if db_path is None or not _os.path.exists(db_path):
        _logger.warning(
            f"systematics: data file not found at {db_path!r}; "
            "returning empty database."
        )
        return db
    db.load(db_path)
    return db


def load_slo(directory: str = None, file_path: str = None) -> Database:
    """Load the SLO experimental+systematics GDR parameter database."""
    return _load(slo_local_file_path, directory=directory, file_path=file_path)


def load_smlo(directory: str = None, file_path: str = None) -> Database:
    """Load the SMLO experimental+systematics GDR parameter database."""
    return _load(smlo_local_file_path, directory=directory, file_path=file_path)
