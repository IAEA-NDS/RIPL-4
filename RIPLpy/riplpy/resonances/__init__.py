"""Average parameters of s- and p-wave neutron resonances.

The average resonance parameters recommended for RIPL-4 were taken from RIPL-3
and Brookhaven National Laboratory (BNL, Mughabghab 2018) and compiled by
S. Goriely (August 2025). Each record provides the recommended values from
both evaluations side-by-side: the average resonance spacing, the average
total radiative width, and the neutron strength function, each with its
quoted uncertainty.

The recommended parameters are given separately for s- and p-wave resonances
in files ``resonances_L0.dat`` and ``resonances_L1.dat`` respectively. The
s-wave file covers 324 nuclides; the p-wave file covers 248 nuclides.

These RIPL-4 filenames differ from the legacy RIPL-3 distribution (which ships
the data in an older, incompatible format as ``resonances0.dat`` /
``resonances1.dat``). When the RIPL-4 files are absent, the section is skipped
with a warning rather than aborting the load, mirroring the densities, fission
and gamma sections, so a partial installation still loads cleanly.

References
----------
S.F. Mughabghab, *Atlas of Neutron Resonances*, Elsevier (2018).
R. Capote *et al.*, *Nuclear Data Sheets* **110**, 3107 (2009).
"""

# ========================

# Logging
from riplpy import logger as _logger

# RIPLpy
from riplpy.config import resolve_directory as _resolve_directory
from riplpy.db import _DbAccessor
from riplpy.exceptions import RiplFileNotFoundError as _RiplFileNotFoundError

# Local
from . import swave as swave
from . import pwave as pwave

# ========================

db = _DbAccessor()

# ========================

def _safe_load(name: str, loader, *args, **kwargs):
    """Run a sub-loader, returning ``None`` and warning if data is missing.

    The s-/p-wave RIPL-4 files (``resonances_L0.dat`` / ``resonances_L1.dat``)
    are absent under the legacy RIPL-3 layout, and the older ``resonances0.dat``
    / ``resonances1.dat`` files use a different, unsupported format. This helper
    lets ``load(directory)`` populate whatever is available without raising.
    """
    try:
        return loader(*args, **kwargs)
    except (FileNotFoundError, _RiplFileNotFoundError) as exc:
        _logger.warning(f"Resonances: skipping {name} (data file not available): {exc}")
        return None
    except Exception as exc:  # noqa: BLE001
        _logger.warning(f"Resonances: skipping {name} (loader error): {exc}")
        return None


def load(directory: str = None) -> None:
    """Load the entire database of the Resonances section of RIPL.

       Loaders for files that are absent (or in the incompatible legacy RIPL-3
       format) are skipped with a warning, so a partial installation still
       loads without raising.

       Args:
           directory: Path to the RIPL database. If None, the configured path
               is used (set_path(), RIPL_LOCATION env var, ~/.riplpyrc, or an
               auto-detected location).

       Once this method is called, the databases cannot be unloaded.
    """
    directory = _resolve_directory(directory)

    # The s-/p-wave resonance databases
    db.swave = _safe_load('swave', swave.load, directory)
    db.pwave = _safe_load('pwave', pwave.load, directory)
