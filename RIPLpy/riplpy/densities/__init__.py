
# RIPLpy
from riplpy import logger as _logger
from riplpy.config import resolve_directory as _resolve_directory
from riplpy.db import _DbAccessor
from riplpy.exceptions import RiplFileNotFoundError

# Local
from . import bsfg as bsfg
from . import ct as ct
from . import egsm as egsm
from . import egsm_norm as egsm_norm
from . import shell_corr as shell_corr
from . import hfb as hfb
from . import bsk14_comb as bsk14_comb
from . import bskg3_comb as bskg3_comb
from . import qrpabe as qrpabe
from . import thfb_comb as thfb_comb

# ========================

db = _DbAccessor()

# ========================

# Model map: user-friendly name -> db attribute name
LEVEL_DENSITY_MODELS = {
    'bsfg': 'bsfg',
    'ct': 'ct',
    'egsm': 'egsm',
    'egsm_norm': 'egsm_norm',
    'hfb': 'hfb',
    'bsk14_comb': 'bsk14_comb',
    'bskg3_comb': 'bskg3_comb',
    'qrpabe': 'qrpabe',
    'thfb_comb': 'thfb_comb',
}

# ========================

def _safe_load(name: str, loader, *args, **kwargs):
    """Run a database loader, returning None and warning if data is missing.

    The github RIPL-4 release does not ship every legacy data file (e.g.
    ``level-densities-bfmeff.dat`` for BSFG). This helper lets the rest of the
    section load without raising when a single backing file is absent.
    """
    try:
        return loader(*args, **kwargs)
    except (FileNotFoundError, RiplFileNotFoundError) as exc:
        _logger.warning(f"Densities: skipping {name} (data file not available): {exc}")
        return None


def load(directory: str = None) -> None:
    """Load the entire database of the Densities section of RIPL.

       Args:
           directory: Path to the RIPL database. If None, the configured path
               is used (set_path(), RIPL_LOCATION env var, ~/.riplpyrc, or an
               auto-detected location).

       Once this method is called, the databases cannot be unloaded.

       Loaders for files that are not present in the github RIPL-4 release are
       gracefully skipped with a warning, so a partial installation will still
       populate the available databases.
    """
    directory = _resolve_directory(directory)

    # The BSFG level density database (not in github release; legacy only)
    db.bsfg = _safe_load('bsfg', bsfg.load, directory=directory)

    # The CT level density database (not in github release; legacy only)
    db.ct = _safe_load('ct', ct.load, directory=directory)

    # The EGSM database
    db.egsm = _safe_load('egsm', egsm.load, directory=directory)

    # The EGSM normalization database
    db.egsm_norm = _safe_load('egsm_norm', egsm_norm.load, directory=directory)

    # The shell correction energies databases
    db.shellcorr_ms = _safe_load('shellcor_ms', shell_corr.load_ms, directory=directory)
    db.shellcorr_mk = _safe_load('shellcor_mk', shell_corr.load_mk, directory=directory)

    # The HFB level density database (not in github release; legacy only)
    db.hfb = _safe_load('hfb', hfb.load_all, directory=directory)

    # Microscopic combinatorial level density databases (per-Z .tab files)
    db.bsk14_comb = _safe_load('bsk14_comb', bsk14_comb.load_all, directory=directory)
    db.bskg3_comb = _safe_load('bskg3_comb', bskg3_comb.load_all, directory=directory)
    db.qrpabe = _safe_load('qrpabe', qrpabe.load_all, directory=directory)
    db.thfb_comb = _safe_load('thfb_comb', thfb_comb.load_all, directory=directory)
