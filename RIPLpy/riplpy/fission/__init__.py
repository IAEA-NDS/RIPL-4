
# RIPLpy
import riplpy as _riplpy
from riplpy.db import _DbAccessor
from riplpy.exceptions import RiplFileNotFoundError as _RiplFileNotFoundError

# Local
from .core import *
from . import bskg3 as bskg3
from . import d1m as d1m
from . import empirical as emp
from . import empirical_new as emp_new
from . import empire as empire
from . import hfb
from . import hfbpath as hfbpath
from . import rmf as rmf
from . import nld_fis as nld_fis

# ========================

saddle_shorthand = {'S' : 'symmetric',
                    'GA': 'axially asymmetric (triaxial)',
                    'MA': 'mass asymmetric',
                    '  ': '',
                   }

db = _DbAccessor()

# ========================

# Model map: user-friendly name -> db attribute name
BARRIER_MODELS = {
    'empirical': 'empirical_barriers',
    'empirical_new': 'empirical_barriers_new',
    'empire': 'empire_barriers',
    'hfb': 'hfb_barriers',
    'bskg3': 'bskg3_barriers',
    'd1m': 'd1m_barriers',
}

# Shorthand aliases so db.bskg3 / db.empirical resolve to the canonical
# '*_barriers' attributes (mirrors the BARRIER_MODELS shorthand names).
for _short, _canonical in BARRIER_MODELS.items():
    db.add_alias(_short, _canonical)

# ========================


def _safe_load(loader, model_name: str, directory: str, db_cls):
    """Load a database, returning an empty instance if the data file is missing.

    Some RIPL-4 release layouts no longer include legacy data files. Rather
    than failing hard, we log a warning and proceed with an empty database so
    that other models remain usable.
    """
    try:
        return loader(directory=directory)
    except (FileNotFoundError, _RiplFileNotFoundError) as exc:
        _riplpy.logger.warning(
            "fission: %s data file not found (%s); using empty database.",
            model_name, exc,
        )
        return db_cls()


def load(directory: str = None, include_paths: bool = False) -> None:
    """Load the entire database of the Fission section of RIPL.

       Args:
           directory: Path to the RIPL database. If None, the configured path
               is used (set_path(), RIPL_LOCATION env var, ~/.riplpyrc, or an
               auto-detected location).
           include_paths: If True, also eagerly load the per-isotope fission
               -path / saddle level-density directories (``hfbpath-bskg3``,
               ``hfbpath-d1m``, ``RMF/Path_Axial``, ``RMF/Path_Triaxial``,
               ``nld-fis-bskg3/Max1``). Default is False because these
               directories contain thousands of nuclei and slow
               ``riplpy.load()`` substantially. The per-element / per-model
               module loaders (e.g. ``fission.hfbpath.load_d1m()``,
               ``fission.rmf.load_axial()``, ``fission.nld_fis.load('Max1')``)
               remain available regardless of this flag.

       The exposed databases (post-load):

       * ``db.hfbpath_bskg3`` -- BSkG3 fission paths (heavy)
       * ``db.hfbpath_d1m``   -- D1M fission paths (heavy)
       * ``db.rmf_axial``     -- axial RMF fission paths (heavy)
       * ``db.rmf_triaxial``  -- triaxial RMF fission paths (heavy)
       * ``db.nld_fis``       -- BSkG3 inner-saddle (Max1) level densities (heavy)

       Once this method is called, the databases cannot be unloaded.
    """
    directory = _riplpy.config.resolve_directory(directory)

    # Empirical fission barrier databases (RIPL-4)
    db.empirical_barriers     = _safe_load(emp.load,     'empirical',     directory, emp.Database)
    db.empirical_barriers_new = _safe_load(emp_new.load, 'empirical_new', directory, emp_new.Database)
    db.empire_barriers        = _safe_load(empire.load,  'empire',        directory, empire.Database)

    # Theoretical BSkG3 fission barrier database
    db.bskg3_barriers = _safe_load(bskg3.load, 'bskg3', directory, bskg3.Database)

    # Theoretical HFB fission barrier database (legacy; file may be absent)
    db.hfb_barriers = _safe_load(hfb.load, 'hfb', directory, hfb.Database)

    # Theoretical D1M HFB fission barrier database (RIPL-4)
    db.d1m_barriers = _safe_load(d1m.load, 'd1m', directory, d1m.Database)

    # Heavy per-isotope fission-path / saddle level-density directories
    # (thousands of nuclei). Default to skipping these — callers opt in via
    # ``include_paths=True`` or use the per-element / per-model loaders.
    if include_paths:
        db.hfbpath_bskg3 = _safe_load(hfbpath.load_bskg3,   'hfbpath_bskg3', directory, hfbpath.Database)
        db.hfbpath_d1m   = _safe_load(hfbpath.load_d1m,     'hfbpath_d1m',   directory, hfbpath.Database)
        db.rmf_axial     = _safe_load(rmf.load_axial,       'rmf_axial',     directory, rmf.Database)
        db.rmf_triaxial  = _safe_load(rmf.load_triaxial,    'rmf_triaxial',  directory, rmf.Database)
        db.nld_fis       = _safe_load(
            lambda directory=None: nld_fis.load('Max1', directory=directory),
            'nld_fis', directory, nld_fis.Database,
        )
    else:
        db.hfbpath_bskg3 = None
        db.hfbpath_d1m   = None
        db.rmf_axial     = None
        db.rmf_triaxial  = None
        db.nld_fis       = None
