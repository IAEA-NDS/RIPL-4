
# RIPLpy
from riplpy import logger as _logger
from riplpy.config import resolve_directory as _resolve_directory
from riplpy.db import _DbAccessor
from riplpy.exceptions import RiplFileNotFoundError as _RiplFileNotFoundError

# Local
from .core import *
from . import gdr
from . import gsf
from . import exp
from . import systematics
from . import d1m
from . import smlo_e1
from . import smlo_m1
from . import tlo
from . import psf

# ========================

db = _DbAccessor()

# ========================

def _safe_load(name: str, loader, *args, **kwargs):
    """Run a sub-loader, returning ``None`` and warning if data is missing.

    The github RIPL-4 release reorganises the gamma section substantially and
    does not ship every legacy file. This helper lets ``load(directory)``
    populate as many databases as are available without raising.
    """
    try:
        return loader(*args, **kwargs)
    except (FileNotFoundError, _RiplFileNotFoundError) as exc:
        _logger.warning(f"Gamma: skipping {name} (data file not available): {exc}")
        return None
    except Exception as exc:  # noqa: BLE001
        _logger.warning(f"Gamma: skipping {name} (loader error): {exc}")
        return None


def load(directory: str = None, include_heavy: bool = False) -> None:
    """Load the entire Gamma section of RIPL from ``directory``.

    Loaders for files that are absent in the github RIPL-4 release are
    skipped with a warning, so a partial installation will still populate
    the databases that are available.

    Args:
        directory: Path to the RIPL database. If None, the configured path
            is used (set_path(), RIPL_LOCATION env var, ~/.riplpyrc, or an
            auto-detected location).
        include_heavy: If True, also eagerly load the per-nucleus tables in
            ``gamma/d1m``, ``gamma/smlo_E1`` (~8980 files), ``gamma/smlo_M1``,
            ``gamma/tlo``, and the PSFDatabase-v2024.1 directory. Default is
            False because these directories contain thousands of files and
            slow ``riplpy.load()`` substantially. Each of these databases is
            still available on demand via ``riplpy.gamma.<module>.load_all()``
            or ``load_element()`` / ``load_nucleus()`` helpers.

    The exposed databases (post-load):

    * ``db.gsf``                       -- legacy gamma-strength-micro (falls back to D1M when heavy)
    * ``db.theory_gdr``                -- legacy theoretical GDR fits
    * ``db.experiment_slo``            -- RIPL-4 recommended SLO experimental fits
    * ``db.experiment_smlo``           -- RIPL-4 recommended SMLO experimental fits
    * ``db.experiment_mlo``            -- alias of ``experiment_smlo`` (back-compat)
    * ``db.experiment_slo_errors``     -- SLO fits with 1-sigma uncertainties
    * ``db.experiment_smlo_errors``    -- SMLO fits with 1-sigma uncertainties
    * ``db.experiment_systematics_slo``  -- broader SLO experiment+systematics fits
    * ``db.experiment_systematics_smlo`` -- broader SMLO experiment+systematics fits
    * ``db.gsf_d1m``                   -- D1M+QRPA per-Z PSF tables (heavy)
    * ``db.smlo_e1``                   -- SMLO E1 photoabsorption tables, per-nucleus (heavy)
    * ``db.smlo_m1``                   -- SMLO M1 strength tables, per-Z (heavy)
    * ``db.tlo``                       -- TLO E1 strength tables, per-Z (heavy)
    * ``db.psf``                       -- PSFDatabase-v2024.1 experimental PSF data (heavy)
    """

    directory = _resolve_directory(directory)

    # Theoretical GDR parameters (legacy single-file product)
    db.theory_gdr = _safe_load('theory_gdr', gdr.load, directory=directory)

    # Experimental SLO/SMLO parameter databases (new RIPL-4 layout)
    db.experiment_slo  = _safe_load('experiment_slo',  exp.load_slo,  directory=directory)
    db.experiment_smlo = _safe_load('experiment_smlo', exp.load_smlo, directory=directory)
    # Backwards-compatible alias (MLO was replaced by SMLO in RIPL-4)
    db.experiment_mlo  = db.experiment_smlo

    # Versions with 1-sigma uncertainties (two-line records)
    db.experiment_slo_errors  = _safe_load(
        'experiment_slo_errors',  exp.load_slo,  directory=directory, errors=True,
    )
    db.experiment_smlo_errors = _safe_load(
        'experiment_smlo_errors', exp.load_smlo, directory=directory, errors=True,
    )

    # Broader systematics-extended fits
    db.experiment_systematics_slo  = _safe_load(
        'experiment_systematics_slo', systematics.load_slo, directory=directory,
    )
    db.experiment_systematics_smlo = _safe_load(
        'experiment_systematics_smlo', systematics.load_smlo, directory=directory,
    )

    # Heavy per-nucleus / per-Z databases (thousands of files). Default to
    # skipping these — callers can opt in via ``include_heavy=True`` or use
    # the per-element / per-nucleus loaders directly.
    if include_heavy:
        db.gsf      = _safe_load('gsf',      gsf.load_all,     directory)
        db.gsf_d1m  = _safe_load('gsf_d1m',  d1m.load_all,     directory)
        db.smlo_e1  = _safe_load('smlo_e1',  smlo_e1.load_all, directory)
        db.smlo_m1  = _safe_load('smlo_m1',  smlo_m1.load_all, directory)
        db.tlo      = _safe_load('tlo',      tlo.load_all,     directory)
    else:
        db.gsf      = None
        db.gsf_d1m  = None
        db.smlo_e1  = None
        db.smlo_m1  = None
        db.tlo      = None

    # Experimental PSF database (partial coverage) -- also heavy
    if include_heavy:
        db.psf  = _safe_load('psf', psf.load_all, directory)
    else:
        db.psf  = None
