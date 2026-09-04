
# RIPLpy
from riplpy.config import resolve_directory as _resolve_directory
from riplpy.db import _DbAccessor

# Local
from .core import *
from . import ct as ct
from . import discrete as discrete

# ========================

db = _DbAccessor()

# Shorthand aliases for the verbose canonical names.
db.add_alias('ct', 'constant_temperature')
db.add_alias('discrete', 'discrete_levels')

# ========================

def load(directory: str = None) -> None:
    """Load the entire database of the Levels section of RIPL.

       Args:
           directory: Path to the RIPL database. If None, the configured path
               is used (set_path(), RIPL_LOCATION env var, ~/.riplpyrc, or an
               auto-detected location).

       Once this method is called, the databases cannot be unloaded.
       Warning: this method can take a long time to complete.
    """
    directory = _resolve_directory(directory)

    # The entire discrete level database
    db.discrete_levels = discrete.load_all(directory)

    # The constant temperature fit database
    db.constant_temperature = ct.load(directory)

