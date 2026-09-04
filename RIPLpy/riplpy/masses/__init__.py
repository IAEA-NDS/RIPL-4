
# Logging
import logging as _logging

# RIPLpy
from riplpy.config import resolve_directory as _resolve_directory
from riplpy.db import _DbAccessor
from riplpy.exceptions import RiplFileNotFoundError as _RiplFileNotFoundError

# Module-level logger (parent: riplpy)
_logger = _logging.getLogger('riplpy.masses')

# Local
from .core import *
from . import ame20 as ame20
from . import bskg3 as bskg3
from . import d1m as d1m
from . import frdm12 as frdm12
from . import frdm95 as frdm95
from . import hfb14 as hfb14
from . import hfb27 as hfb27
from . import ws4 as ws4
from . import ab as ab
from . import density_bskg3 as density_bskg3
from . import density_d1m as density_d1m
from . import deformations as deformations

# ========================

db = _DbAccessor()

# Shorthand <-> longhand aliases so db.frdm12 resolves to db.frdm2012, etc.
# (the canonical attribute carries the full year; the abbreviation is the alias)
db.add_alias('frdm12', 'frdm2012')
db.add_alias('frdm95', 'frdm1995')
db.add_alias('ame2020', 'ame20')

# ========================

# Model map: user-friendly name -> (db attribute name, mass excess field name)
# This allows users to use shorthand names like 'frdm12' instead of 'frdm2012'
MASS_MODELS = {
    'ame20': ('ame20', 'Mexp'),
    'frdm12': ('frdm2012', 'Mth'),
    'frdm2012': ('frdm2012', 'Mth'),
    'frdm95': ('frdm1995', 'Mth'),
    'frdm1995': ('frdm1995', 'Mth'),
    'hfb14': ('hfb14', 'Mth'),
    'hfb27': ('hfb27', 'Mth'),
    'bskg3': ('bskg3', 'Mth'),
    'd1m': ('d1m', 'Mth'),
    'ws4': ('ws4', 'Mth'),
}

# ========================

def _safe_load(name: str, loader_callable, directory: str, required: bool = True):
    """Load a mass model database, returning an empty Database if the data file is missing.

    The official RIPL-4 GitHub release does not currently ship every mass model
    data file (notably ``mass-frdm95.dat`` and ``mass-hfb14.dat``). To keep
    ``riplpy.load()`` usable against the GitHub layout while still supporting
    users who have a full release, we catch the missing-file error, log a
    warning, and return an empty database for that model.

    Arguments:
        name: Human-readable model name (used in the warning message).
        loader_callable: The model's ``load`` partial (e.g. ``frdm95.load``).
        directory: Path passed through to the loader.
        required: If True (default), missing files re-raise. Set to False for
            optional models that are known to be absent from some releases.
    """
    try:
        return loader_callable(directory)
    except _RiplFileNotFoundError as exc:
        if required:
            raise
        _logger.warning(
            "Skipping %s mass database: data file not found in this RIPL "
            "release (%s). The reader/writer remain available for users with "
            "the full RIPL distribution.", name, exc
        )
        # Return an empty Database of the appropriate type so attribute access
        # still works (e.g. for in_ripl / list_databases iteration).
        return loader_callable.keywords['db_obj']()


def load_only_masses(directory: str) -> None:
    """Loads only the mass database of the Masses section of RIPL.

       Once this method is called, the databases cannot be unloaded.

       Missing data files for ``frdm95`` and ``hfb14`` are skipped with a
       warning; those models are not part of the current RIPL-4 GitHub release
       but the readers/writers stay available for users with a full release.
    """
    # AME2020 mass evaluation
    db.ame20 = _safe_load('AME2020', ame20.load, directory)

    # BSkG3 mass database
    db.bskg3 = _safe_load('BSkG3', bskg3.load, directory)

    # D1M mass database
    db.d1m = _safe_load('D1M', d1m.load, directory)

    # FRDM1995 mass database (not in current github release)
    db.frdm1995 = _safe_load('FRDM1995', frdm95.load, directory, required=False)

    # FRDM2012 mass database
    db.frdm2012 = _safe_load('FRDM2012', frdm12.load, directory)

    # HFB14 mass database (not in current github release)
    db.hfb14 = _safe_load('HFB14', hfb14.load, directory, required=False)

    # HFB27 mass database
    db.hfb27 = _safe_load('HFB27', hfb27.load, directory)

    # WS4 mass database
    db.ws4 = _safe_load('WS4', ws4.load, directory)


def load(directory: str = None, include_heavy: bool = False) -> None:
    """Load the entire database of the Masses section of RIPL.

    Args:
        directory: Path to the RIPL database. If None, the configured path
            is used (set_path(), RIPL_LOCATION env var, ~/.riplpyrc, or an
            auto-detected location).
        include_heavy: If True, also load the per-Z auxiliary density tables
            ``masses/Density-bskg3/`` and ``masses/Density-d1m/`` (200+ files
            in total). Default is False because these tables can take ~60s to
            load and are rarely needed alongside the mass-model fits. Use
            ``riplpy.masses.density_bskg3.load_all(directory)`` or
            ``density_d1m.load_element(Z, directory)`` to load on demand.

    Once this method is called, the databases cannot be unloaded.
    """
    directory = _resolve_directory(directory)

    # Load masses
    load_only_masses(directory)

    # Natural abundance database
    db.natab = ab.load(directory)

    # Ground state deformations
    db.deformations = deformations.load(directory)

    # Heavy per-Z auxiliary density tables (skip by default).
    if include_heavy:
        db.density_bskg3 = density_bskg3.load_all(directory)
        db.density_d1m   = density_d1m.load_all(directory)
    else:
        db.density_bskg3 = None
        db.density_d1m   = None
