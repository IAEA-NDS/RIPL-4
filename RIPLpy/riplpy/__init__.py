

# Code version
__version__ = "0.6.0"

# Logging
import logging as _logging

# Set up library-level logger with NullHandler (users can configure their own handlers)
logger = _logging.getLogger(__name__)
logger.addHandler(_logging.NullHandler())

# Local
from . import config as config
from .collections import *
from .models import Models
from . import densities
from . import fission
from . import gamma
from . import levels
from . import masses
from . import optical
from . import resonances


# Supported RIPL sections
sections = ('densities', 'fission', 'gamma', 'levels', 'masses', 'optical', 'resonances')

# Public API surface (explicit for tooling, autodoc, and ML/AI introspection)
__all__ = [
    # Core types
    'Nuclide', 'Nucleus', 'Element', 'Elements', 'Models',
    # Sections (submodules)
    'densities', 'fission', 'gamma', 'levels', 'masses', 'optical',
    'resonances', 'config', 'sections',
    # Configuration
    'set_path', 'get_path', 'load',
    # Discovery / introspection
    'in_ripl', 'in_sections', 'in_dbs', 'list_sections', 'list_databases',
    'list_nuclei', 'get_database',
    # Uniform ML/AI ingestion entry points
    'to_dataframe', 'to_records', 'to_numpy',
    # Scalar convenience getters
    'get_mass', 'get_mass_entry', 'get_level_density', 'get_gdr',
    'get_resonance', 'get_fission_barrier',
    # Batch getters
    'get_masses', 'get_mass_entries', 'get_level_densities', 'get_gdrs',
    'get_resonances', 'get_fission_barriers',
    # Optical model
    'get_omp', 'list_omps', 'find_omp', 'get_deformation',
    'get_omp_reference',
]


def set_path(path: str) -> None:
    """Set the RIPL database path.

    This is the recommended way to configure RIPLpy if the RIPL_LOCATION
    environment variable is not set.

    Args:
        path: Absolute path to the RIPL database directory

    Raises:
        ValueError: If the path does not exist or is not a directory

    Examples:
        >>> import riplpy
        >>> riplpy.set_path('/path/to/RIPL-4')
        >>> riplpy.load()
    """
    config.set_path(path)


def get_path() -> str:
    """Get the current RIPL database path.

    Returns:
        The resolved RIPL path, or empty string if not configured

    Examples:
        >>> import riplpy
        >>> print(riplpy.get_path())
        '/usr/local/share/RIPL-4'
    """
    return config.get_path()


def _parallel_load(section: str, directory: str = None) -> None:
   """Load a single section by name, used by the threaded quick loader.

   Runs inside a worker thread (see :func:`load` with ``quick=True``). Threads
   share the interpreter's memory, so the section's module-level ``db`` is
   populated in the same process the caller is using — unlike a process pool,
   whose workers would discard their results.
   """
   if section in sections:
      try:
         module = globals()[section]
         module.load(directory=directory)
      except Exception as exc:  # noqa: BLE001
         logger.warning(
            f"load(): skipping section '{section}' (loader error): {exc}"
         )


def load(directory: str = None, quick: bool = False) -> None:
   """Load the entire RIPL database from a specified directory on the local machine.

   Args:
       directory: Path to the RIPL database. If not provided, uses the
                  configured path (from set_path(), RIPL_LOCATION env var,
                  ~/.riplpyrc, or auto-detected location).
       quick: If True, load the sections concurrently with a thread pool.
              The readers are I/O-bound, so this overlaps their file parsing
              while still populating the databases in this process.

   Raises:
       ValueError: If no directory is provided and no path is configured.

   Examples:
       >>> import riplpy
       >>> riplpy.load('/path/to/RIPL')  # Explicit path
       >>> # Or configure first, then load
       >>> riplpy.set_path('/path/to/RIPL')
       >>> riplpy.load()

   Note:
       Once this method is called, the databases cannot be unloaded.

       A section whose backing data is absent (e.g. a legacy RIPL-3 layout, or
       a partial install) is skipped with a warning rather than aborting the
       whole load, so the databases that *are* available still populate.
   """
   if directory is None:
      directory = config.get_path()
      if not directory:
         raise ValueError(
            "No RIPL directory specified and no path configured. "
            "Either pass a directory to load(), or configure the path first:\n"
            "  - riplpy.set_path('/path/to/RIPL')\n"
            "  - Set RIPL_LOCATION environment variable\n"
            "  - Create ~/.riplpyrc config file"
         )

   if quick:
      # Load sections concurrently with threads. The readers are I/O-bound
      # (parsing data files), so they release the GIL during file reads and
      # genuinely overlap. Threads share memory, so each section's db is
      # populated in this process -- a process pool would load into workers
      # and discard the results. Each section writes to its own module-level
      # db accessor, so there is no cross-thread contention.
      from concurrent.futures import ThreadPoolExecutor as _ThreadPoolExecutor
      with _ThreadPoolExecutor(max_workers=len(sections)) as _executor:
         _futures = [
            _executor.submit(_parallel_load, section, directory)
            for section in sections
         ]
         for _future in _futures:
            _future.result()  # surface unexpected (non-loader) exceptions
   else:
      for section in sections:
         module = globals()[section]
         try:
            module.load(directory=directory)
         except Exception as exc:  # noqa: BLE001
            logger.warning(
               f"load(): skipping section '{section}' (loader error): {exc}"
            )


def in_ripl(x: Nuclide | Element) -> bool:
   """Check if the given Nuclide (or Element) appears anywhere in the RIPL database. Returns True if so, False otherwise. 
      This method may only be used once riplpy.load() has been issued.

      Input:
         x [Nuclide | Element]: a Nuclide or Element object that may be in the RIPL database

      Output:
         [bool]: True if x in RIPL, False otherwise

      Examples:
         (1) Check if input is in the RIPL database::
            $ import riplpy
            $ riplpy.load()
            $ n = riplpy.Nuclide(z=50, a=130)
            $ riplpy.in_ripl(n) # True
   """
   # Loop over each section of RIPL
   for section in sections:
      # Check the databases in this section
      for db in getattr(globals()[section], 'db')():
          db_obj = getattr(getattr(globals()[section], 'db'), db)
          if db_obj is None:
              # Section loader skipped this database (data file missing).
              continue
          if x in db_obj.data:
              return True
   return False


def in_sections(x: Nuclide | Element) -> list:
   """Check if the given Nuclide (or Element) is in various sections of the RIPL database. Returns a list of sections where True. 
      This method may only be used once riplpy.load() has been issued.

      Input:
         x [Nuclide | Element]: a Nuclide or Element object that may be in the RIPL database

      Output:
         [list]: A list of section names (as a string) x is in; empty otherwise

      Examples:
         (1) Return the section names that the input is in::
            $ import riplpy
            $ riplpy.load() # This will take some time
            $ n = riplpy.Nuclide(Z=46, a=119)
            $ riplpy.in_sections(n) # ['densities', 'gamma', 'levels', 'masses']
   """
   
   _sections = []
   # Loop over each section of RIPL
   for section in sections:
      # Check the databases in this section
      for db in getattr(globals()[section], 'db')():
          db_obj = getattr(getattr(globals()[section], 'db'), db)
          if db_obj is None:
              continue
          if x in db_obj.data:
              _sections.append(section)
              break
   return _sections


def in_dbs(x: Nuclide | Element) -> list:
   """Check if the given Nuclide (or Element) is in various databases of RIPL. Returns a list of database names where the object can be found.
      This method may only be used once riplpy.load() has been issued.

      Input:
         x [Nuclide | Element]: a Nuclide or Element object that may be in the RIPL database

      Output:
         [list]: A list of tuples that includes the section name (as a string) and the database name (as a string) that x is in; empty otherwise
   """
   dbs = []
   # Loop over each section of RIPL
   for section in sections:
      # Check the databases in this section
      for db in getattr(globals()[section], 'db')():
          db_obj = getattr(getattr(globals()[section], 'db'), db)
          if db_obj is None:
              continue
          if x in db_obj.data:
              dbs.append((section, db))
   return dbs


def list_sections() -> list:
    """List all available RIPL sections.

    Returns:
        List of section names (strings)

    Examples:
        >>> riplpy.list_sections()
        ['densities', 'fission', 'gamma', 'levels', 'masses', 'optical', 'resonances']
    """
    return list(sections)


def list_databases(section: str = None) -> list | dict:
    """List available databases, optionally filtered by section.

    Args:
        section: Section name to filter by (e.g., 'masses', 'densities').
                 If None, returns a dict of all sections and their databases.

    Returns:
        If section is provided: List of database names in that section
        If section is None: Dict mapping section names to lists of database names

    Examples:
        >>> riplpy.list_databases('masses')
        ['ame20', 'bskg3', 'd1m', 'frdm1995', 'frdm2012', ...]
        >>> riplpy.list_databases()
        {'densities': ['bsfg', 'ct', ...], 'masses': ['ame20', ...], ...}

    Raises:
        ValueError: If the section name is not recognized

    Note:
        Databases must be loaded first via riplpy.load() for this to return results.
    """
    if section is not None:
        if section not in sections:
            raise ValueError(f"Unknown section '{section}'. Available: {list(sections)}")
        module = globals()[section]
        return list(getattr(module, 'db')())

    # Return all sections
    result = {}
    for sec in sections:
        module = globals()[sec]
        result[sec] = list(getattr(module, 'db')())
    return result


def list_nuclei(database: str, section: str = None) -> list:
    """List all nuclei available in a specific database.

    Args:
        database: Database name (e.g., 'ame20', 'bsfg'). Can also use
                  'section.database' format (e.g., 'masses.ame20').
        section: Section name. Required if database name is ambiguous
                 across sections, optional if using 'section.database' format.

    Returns:
        List of Nuclide objects available in the database

    Examples:
        >>> riplpy.list_nuclei('ame20')  # Unambiguous
        [Nuclide(Z=1, A=1), Nuclide(Z=1, A=2), ...]
        >>> riplpy.list_nuclei('masses.ame20')  # Explicit section
        [Nuclide(Z=1, A=1), Nuclide(Z=1, A=2), ...]
        >>> riplpy.list_nuclei('ct', section='densities')  # Disambiguate
        [Nuclide(Z=11, A=24), ...]

    Raises:
        ValueError: If the database is not found or is ambiguous

    Note:
        Databases must be loaded first via riplpy.load() for this to return results.
    """
    # Handle 'section.database' format
    if '.' in database:
        parts = database.split('.', 1)
        section = parts[0]
        database = parts[1]

    # If section is specified, use it directly
    if section is not None:
        if section not in sections:
            raise ValueError(f"Unknown section '{section}'. Available: {list(sections)}")
        module = globals()[section]
        db_accessor = getattr(module, 'db')
        if not db_accessor.has(database):
            available = list(db_accessor())
            raise ValueError(f"Database '{database}' not found in section '{section}'. Available: {available}")
        db = getattr(db_accessor, db_accessor.resolve(database))
        return list(db.data.keys())

    # Search all sections for the database
    found_in = []
    for sec in sections:
        module = globals()[sec]
        db_accessor = getattr(module, 'db')
        if db_accessor.has(database):
            found_in.append(sec)

    if len(found_in) == 0:
        raise ValueError(f"Database '{database}' not found in any section")
    elif len(found_in) > 1:
        raise ValueError(f"Database '{database}' is ambiguous, found in sections: {found_in}. "
                         f"Please specify section, e.g., list_nuclei('{found_in[0]}.{database}')")

    # Unambiguous - return the nuclei
    module = globals()[found_in[0]]
    db_accessor = getattr(module, 'db')
    db = getattr(db_accessor, db_accessor.resolve(database))
    return list(db.data.keys())


def get_database(database: str, section: str = None) -> object:
    """Get a database object by name.

    Args:
        database: Database name (e.g., 'ame20', 'bsfg'). Can also use
                  'section.database' format (e.g., 'masses.ame20').
        section: Section name. Required if database name is ambiguous
                 across sections, optional if using 'section.database' format.

    Returns:
        The database object

    Examples:
        >>> db = riplpy.get_database('ame20')
        >>> db = riplpy.get_database('masses.ame20')
        >>> db = riplpy.get_database('ct', section='levels')

    Raises:
        ValueError: If the database is not found or is ambiguous

    Note:
        Databases must be loaded first via riplpy.load() for this to return results.
    """
    # Handle 'section.database' format
    if '.' in database:
        parts = database.split('.', 1)
        section = parts[0]
        database = parts[1]

    # If section is specified, use it directly
    if section is not None:
        if section not in sections:
            raise ValueError(f"Unknown section '{section}'. Available: {list(sections)}")
        module = globals()[section]
        db_accessor = getattr(module, 'db')
        if not db_accessor.has(database):
            available = list(db_accessor())
            raise ValueError(f"Database '{database}' not found in section '{section}'. Available: {available}")
        return getattr(db_accessor, db_accessor.resolve(database))

    # Search all sections for the database
    found_in = []
    for sec in sections:
        module = globals()[sec]
        db_accessor = getattr(module, 'db')
        if db_accessor.has(database):
            found_in.append(sec)

    if len(found_in) == 0:
        raise ValueError(f"Database '{database}' not found in any section")
    elif len(found_in) > 1:
        raise ValueError(f"Database '{database}' is ambiguous, found in sections: {found_in}. "
                         f"Please specify section, e.g., get_database('{found_in[0]}.{database}')")

    # Unambiguous - return the database
    module = globals()[found_in[0]]
    db_accessor = getattr(module, 'db')
    return getattr(db_accessor, db_accessor.resolve(database))


def to_dataframe(database: str, section: str = None) -> "pandas.DataFrame":
    """Return a loaded database as a pandas DataFrame.

    A single uniform ML-ingestion entry point: every section's databases —
    scalar or array/spectral — flatten to a tidy DataFrame with the nuclide
    Z/A/symbol columns first. Array fields are kept as native Python lists so
    downstream code can consume them directly.

    Args:
        database: Database name or 'section.database' (e.g. 'masses.ame20',
                  'gamma.experiment_slo', 'fission.bskg3').
        section: Optional section name if the database name is ambiguous.

    Returns:
        pandas.DataFrame of the flattened entries.

    Raises:
        ValueError: If the database is unknown/ambiguous, or is not loaded
            (call riplpy.load() or the section/per-element loader first).

    Examples:
        >>> riplpy.load()
        >>> df = riplpy.to_dataframe('masses.ame20')
        >>> df = riplpy.to_dataframe('fission.bskg3')

    Note:
        Databases must be loaded first (riplpy.load(), the section load(), or
        a per-element/per-nucleus helper for the heavy array datasets).
    """
    db = get_database(database, section)
    if db is None:
        raise ValueError(
            f"Database '{database}' is not loaded (it may be a RIPL-3 legacy "
            f"product absent from this release, or a heavy dataset not yet "
            f"loaded). Load it first via riplpy.load() or the per-element "
            f"loader."
        )
    return db.to_dataframe()


def to_records(database: str, section: str = None) -> list:
    """Return a loaded database as a list of serialized entry dictionaries.

    Like :func:`to_dataframe` but returns plain JSON-serializable records
    (Nuclide objects become ``{'Z','A','symbol'}``; arrays stay lists) — a
    uniform, dependency-free ML/AI ingestion format.

    Args:
        database: Database name or 'section.database'.
        section: Optional section name if ambiguous.

    Returns:
        list[dict]: One serialized record per entry (per dataset for the
        list-valued PSF database).

    Raises:
        ValueError: If the database is unknown/ambiguous or not loaded.
    """
    db = get_database(database, section)
    if db is None:
        raise ValueError(
            f"Database '{database}' is not loaded (RIPL-3 legacy or a heavy "
            f"dataset not yet loaded). Load it first."
        )
    return db.to_list()


def to_numpy(database: str, section: str = None, structured: bool = True):
    """Return a loaded database as a numpy array for ML/AI pipelines.

    A uniform numpy ingestion entry point parallel to :func:`to_dataframe`.

    Args:
        database: Database name or 'section.database' (e.g. 'masses.ame20',
                  'gamma.experiment_slo', 'fission.bskg3_barriers').
        section: Optional section name if the database name is ambiguous.
        structured: If True (default), return a numpy structured (record)
            array with every column preserved (Z/A/symbol first). If False,
            return ``(X, columns)`` where ``X`` is a dense float64 feature
            matrix of the numeric columns only.

    Returns:
        numpy.ndarray when ``structured`` is True; otherwise a
        ``(numpy.ndarray, list[str])`` tuple.

    Raises:
        ValueError: If the database is unknown/ambiguous or not loaded.
        ImportError: If numpy is not installed (pip install 'riplpy[numpy]').

    Examples:
        >>> arr = riplpy.to_numpy('masses.ame20')
        >>> X, cols = riplpy.to_numpy('masses.ame20', structured=False)
    """
    db = get_database(database, section)
    if db is None:
        raise ValueError(
            f"Database '{database}' is not loaded (RIPL-3 legacy or a heavy "
            f"dataset not yet loaded). Load it first."
        )
    return db.to_numpy(structured=structured)


def _resolve_nucleus(Z_or_nuclide, A=None) -> Nuclide:
    """Resolve input to a Nuclide object.

    Accepts either:
        - A Nuclide/Nucleus object directly
        - Z and A as separate arguments
    """
    if isinstance(Z_or_nuclide, (Nuclide, Nucleus)):
        return Z_or_nuclide
    if A is None:
        raise ValueError("Must provide either a Nuclide or both Z and A")
    return Nuclide(Z=Z_or_nuclide, A=A)


def get_mass(Z_or_nuclide, A: int = None, model: str = 'ame20') -> float:
    """Get the mass excess for a nucleus from a specified mass model.

    Args:
        Z_or_nuclide: Either a Nuclide/Nucleus object, or the atomic number Z
        A: Mass number (required if Z_or_nuclide is an integer)
        model: Mass model to use. Options: 'ame20', 'frdm12', 'frdm95',
               'hfb14', 'hfb27', 'bskg3', 'd1m'. Default: 'ame20'

    Returns:
        Mass excess in MeV

    Examples:
        >>> riplpy.get_mass(82, 208)  # Using Z, A
        -21.749
        >>> riplpy.get_mass(Nuclide(Z=82, A=208), model='frdm12')
        -21.352

    Raises:
        ValueError: If the model is not recognized
        NucleusNotFoundError: If the nucleus is not in the database
    """
    n = _resolve_nucleus(Z_or_nuclide, A)

    if model not in masses.MASS_MODELS:
        raise ValueError(f"Unknown mass model '{model}'. Available: {list(masses.MASS_MODELS.keys())}")

    db_name, field = masses.MASS_MODELS[model]
    db = getattr(masses.db, db_name)
    entry = db.get(n)
    return getattr(entry, field)


def get_mass_entry(Z_or_nuclide, A: int = None, model: str = 'ame20') -> object:
    """Get the full mass entry for a nucleus from a specified mass model.

    Args:
        Z_or_nuclide: Either a Nuclide/Nucleus object, or the atomic number Z
        A: Mass number (required if Z_or_nuclide is an integer)
        model: Mass model to use. Options: 'ame20', 'frdm12', 'frdm95',
               'hfb14', 'hfb27', 'bskg3', 'd1m'. Default: 'ame20'

    Returns:
        Entry object with all available mass data for that model

    Raises:
        ValueError: If the model is not recognized
        NucleusNotFoundError: If the nucleus is not in the database
    """
    n = _resolve_nucleus(Z_or_nuclide, A)

    if model not in masses.MASS_MODELS:
        raise ValueError(f"Unknown mass model '{model}'. Available: {list(masses.MASS_MODELS.keys())}")

    db_name, _ = masses.MASS_MODELS[model]
    db = getattr(masses.db, db_name)
    return db.get(n)


def get_level_density(Z_or_nuclide, A: int = None, model: str = 'egsm') -> object:
    """Get level density parameters for a nucleus from a specified model.

    Args:
        Z_or_nuclide: Either a Nuclide/Nucleus object, or the atomic number Z
        A: Mass number (required if Z_or_nuclide is an integer)
        model: Level density model. RIPL-4 options: 'egsm', 'egsm_norm',
               'bsk14_comb', 'bskg3_comb', 'qrpabe', 'thfb_comb'. RIPL-3
               legacy options (full RIPL distribution only): 'bsfg', 'ct',
               'hfb'. Default: 'egsm'

    Returns:
        Entry object containing level density parameters (model-dependent fields)

    Raises:
        ValueError: If the model is not recognized
        NucleusNotFoundError: If the nucleus is not in the database
    """
    n = _resolve_nucleus(Z_or_nuclide, A)

    if model not in densities.LEVEL_DENSITY_MODELS:
        raise ValueError(f"Unknown level density model '{model}'. Available: {list(densities.LEVEL_DENSITY_MODELS.keys())}")

    from riplpy.exceptions import NucleusNotFoundError
    db = getattr(densities.db, densities.LEVEL_DENSITY_MODELS[model])
    if db is None:
        raise NucleusNotFoundError(
            f"Level density model '{model}' database not loaded (data file missing)"
        )
    return db.get(n)


def get_gdr(Z_or_nuclide, A: int = None) -> object:
    """Get Giant Dipole Resonance (GDR) parameters for a nucleus.

    Args:
        Z_or_nuclide: Either a Nuclide/Nucleus object, or the atomic number Z
        A: Mass number (required if Z_or_nuclide is an integer)

    Returns:
        Entry object containing GDR parameters

    Raises:
        NucleusNotFoundError: If the nucleus is not in the database
    """
    from riplpy.exceptions import NucleusNotFoundError
    n = _resolve_nucleus(Z_or_nuclide, A)
    db = gamma.db.theory_gdr
    if db is None:
        raise NucleusNotFoundError(
            "Theoretical GDR database not loaded (gdr-parameters-theor.dat missing)"
        )
    return db.get(n)


def get_resonance(Z_or_nuclide, A: int = None, wave: str = 's') -> object:
    """Get resonance parameters for a nucleus.

    Args:
        Z_or_nuclide: Either a Nuclide/Nucleus object, or the atomic number Z
        A: Mass number (required if Z_or_nuclide is an integer)
        wave: Resonance type - 's' (s-wave) or 'p' (p-wave). Default: 's'

    Returns:
        Entry object containing resonance parameters

    Raises:
        ValueError: If the wave type is not recognized
        RiplFileNotFoundError: If the resonance database was skipped at load
            time (e.g. legacy RIPL-3 layout or a partial install)
        NucleusNotFoundError: If the nucleus is not in the database
    """
    from riplpy.exceptions import RiplFileNotFoundError

    n = _resolve_nucleus(Z_or_nuclide, A)

    if wave == 's':
        wave_db = resonances.db.swave
    elif wave == 'p':
        wave_db = resonances.db.pwave
    else:
        raise ValueError(f"Unknown wave type '{wave}'. Available: 's', 'p'")

    if wave_db is None:
        raise RiplFileNotFoundError(
            f"The {wave}-wave resonance database is not loaded. The RIPL-4 "
            f"file (resonances_L{0 if wave == 's' else 1}.dat) was not found "
            f"at the configured path; the legacy RIPL-3 "
            f"resonances{0 if wave == 's' else 1}.dat format is not supported. "
            f"Point RIPL_LOCATION at a RIPL-4 (github-layout) tree."
        )
    return wave_db.get(n)


def get_fission_barrier(Z_or_nuclide, A: int = None, model: str = 'empirical') -> object:
    """Get fission barrier parameters for a nucleus.

    Args:
        Z_or_nuclide: Either a Nuclide/Nucleus object, or the atomic number Z
        A: Mass number (required if Z_or_nuclide is an integer)
        model: Barrier model. Options: 'empirical', 'empirical_new', 'empire',
               'hfb', 'bskg3'. Default: 'empirical'

    Returns:
        Entry object containing fission barrier parameters

    Raises:
        ValueError: If the model is not recognized
        NucleusNotFoundError: If the nucleus is not in the database
    """
    n = _resolve_nucleus(Z_or_nuclide, A)

    if model not in fission.BARRIER_MODELS:
        raise ValueError(f"Unknown fission barrier model '{model}'. Available: {list(fission.BARRIER_MODELS.keys())}")

    from riplpy.exceptions import NucleusNotFoundError
    db = getattr(fission.db, fission.BARRIER_MODELS[model])
    if db is None:
        raise NucleusNotFoundError(
            f"Fission barrier model '{model}' database not loaded (data file missing)"
        )
    return db.get(n)


def _normalize_nuclei_list(nuclei: list) -> list:
    """Convert a list of nuclei specifications to Nuclide objects.

    Accepts:
        - Nuclide/Nucleus objects
        - Tuples of (Z, A)
        - Lists of [Z, A]

    Args:
        nuclei: List of nuclei in any supported format

    Returns:
        List of Nuclide objects
    """
    result = []
    for item in nuclei:
        if isinstance(item, (Nuclide, Nucleus)):
            result.append(item)
        elif isinstance(item, (tuple, list)) and len(item) >= 2:
            result.append(Nuclide(Z=item[0], A=item[1]))
        else:
            raise ValueError(f"Cannot convert {item} to Nuclide. Expected Nuclide, (Z, A) tuple, or [Z, A] list.")
    return result


def get_masses(nuclei: list, model: str = 'ame20', skip_missing: bool = False) -> dict:
    """Get mass excess values for multiple nuclei.

    Args:
        nuclei: List of nuclei. Each element can be:
                - A Nuclide/Nucleus object
                - A tuple (Z, A)
                - A list [Z, A]
        model: Mass model to use. Options: 'ame20', 'frdm12', 'frdm95',
               'hfb14', 'hfb27', 'bskg3', 'd1m'. Default: 'ame20'
        skip_missing: If True, skip nuclei not in the database.
                      If False (default), raise NucleusNotFoundError.

    Returns:
        Dictionary mapping Nuclide -> mass excess (MeV)

    Examples:
        >>> nuclei = [Nuclide(82, 208), Nuclide(82, 206), (82, 204)]
        >>> masses = riplpy.get_masses(nuclei, model='ame20')
        >>> masses[Nuclide(82, 208)]
        -21.749

        >>> # Using tuples for convenience
        >>> masses = riplpy.get_masses([(26, 56), (28, 58), (50, 120)])

    Raises:
        ValueError: If the model is not recognized
        NucleusNotFoundError: If a nucleus is not found and skip_missing=False
    """
    from riplpy.exceptions import NucleusNotFoundError

    if model not in masses.MASS_MODELS:
        raise ValueError(f"Unknown mass model '{model}'. Available: {list(masses.MASS_MODELS.keys())}")

    db_name, field = masses.MASS_MODELS[model]
    db = getattr(masses.db, db_name)

    nuclide_list = _normalize_nuclei_list(nuclei)
    result = {}

    for n in nuclide_list:
        try:
            entry = db.get(n)
            result[n] = getattr(entry, field)
        except KeyError:
            # NucleusNotFoundError subclasses KeyError; a few Database
            # subclasses still raise the plain KeyError.
            if not skip_missing:
                raise
            # skip_missing=True: just don't include in results

    return result


def get_mass_entries(nuclei: list, model: str = 'ame20', skip_missing: bool = False) -> dict:
    """Get full mass entries for multiple nuclei.

    Args:
        nuclei: List of nuclei. Each element can be:
                - A Nuclide/Nucleus object
                - A tuple (Z, A)
                - A list [Z, A]
        model: Mass model to use. Options: 'ame20', 'frdm12', 'frdm95',
               'hfb14', 'hfb27', 'bskg3', 'd1m'. Default: 'ame20'
        skip_missing: If True, skip nuclei not in the database.
                      If False (default), raise NucleusNotFoundError.

    Returns:
        Dictionary mapping Nuclide -> Entry object

    Examples:
        >>> nuclei = [(82, 208), (82, 206), (82, 204)]
        >>> entries = riplpy.get_mass_entries(nuclei)
        >>> entries[Nuclide(82, 208)].Mexp
        -21.749

    Raises:
        ValueError: If the model is not recognized
        NucleusNotFoundError: If a nucleus is not found and skip_missing=False
    """
    from riplpy.exceptions import NucleusNotFoundError

    if model not in masses.MASS_MODELS:
        raise ValueError(f"Unknown mass model '{model}'. Available: {list(masses.MASS_MODELS.keys())}")

    db_name, _ = masses.MASS_MODELS[model]
    db = getattr(masses.db, db_name)

    nuclide_list = _normalize_nuclei_list(nuclei)
    result = {}

    for n in nuclide_list:
        try:
            result[n] = db.get(n)
        except KeyError:
            # NucleusNotFoundError subclasses KeyError; a few Database
            # subclasses still raise the plain KeyError.
            if not skip_missing:
                raise

    return result


def get_level_densities(nuclei: list, model: str = 'egsm', skip_missing: bool = False) -> dict:
    """Get level density parameters for multiple nuclei.

    Args:
        nuclei: List of nuclei. Each element can be:
                - A Nuclide/Nucleus object
                - A tuple (Z, A)
                - A list [Z, A]
        model: Level density model. RIPL-4 options: 'egsm', 'egsm_norm',
               'bsk14_comb', 'bskg3_comb', 'qrpabe', 'thfb_comb'. RIPL-3
               legacy options (full RIPL distribution only): 'bsfg', 'ct',
               'hfb'. Default: 'egsm'
        skip_missing: If True, skip nuclei not in the database.
                      If False (default), raise NucleusNotFoundError.

    Returns:
        Dictionary mapping Nuclide -> Entry object

    Examples:
        >>> nuclei = [(11, 24), (12, 26), (13, 28)]
        >>> densities = riplpy.get_level_densities(nuclei, model='egsm')

    Raises:
        ValueError: If the model is not recognized
        NucleusNotFoundError: If a nucleus is not found and skip_missing=False
    """
    from riplpy.exceptions import NucleusNotFoundError

    if model not in densities.LEVEL_DENSITY_MODELS:
        raise ValueError(f"Unknown level density model '{model}'. Available: {list(densities.LEVEL_DENSITY_MODELS.keys())}")

    db = getattr(densities.db, densities.LEVEL_DENSITY_MODELS[model])

    nuclide_list = _normalize_nuclei_list(nuclei)
    result = {}

    # The underlying database may be ``None`` if the backing data file is
    # absent from this RIPL layout. Honour ``skip_missing`` in that case.
    if db is None:
        if skip_missing:
            return result
        raise NucleusNotFoundError(
            f"Level density model '{model}' database not loaded (data file missing)"
        )

    for n in nuclide_list:
        try:
            result[n] = db.get(n)
        except KeyError:
            # NucleusNotFoundError subclasses KeyError; a few Database
            # subclasses still raise the plain KeyError.
            if not skip_missing:
                raise

    return result


def get_gdrs(nuclei: list, skip_missing: bool = False) -> dict:
    """Get Giant Dipole Resonance (GDR) parameters for multiple nuclei.

    Args:
        nuclei: List of nuclei. Each element can be:
                - A Nuclide/Nucleus object
                - A tuple (Z, A)
                - A list [Z, A]
        skip_missing: If True, skip nuclei not in the database.
                      If False (default), raise NucleusNotFoundError.

    Returns:
        Dictionary mapping Nuclide -> Entry object

    Examples:
        >>> nuclei = [(82, 208), (50, 120), (28, 58)]
        >>> gdrs = riplpy.get_gdrs(nuclei)

    Raises:
        NucleusNotFoundError: If a nucleus is not found and skip_missing=False
    """
    from riplpy.exceptions import NucleusNotFoundError

    nuclide_list = _normalize_nuclei_list(nuclei)
    result = {}

    db = gamma.db.theory_gdr
    if db is None:
        if skip_missing:
            return result
        raise NucleusNotFoundError(
            "Theoretical GDR database not loaded (gdr-parameters-theor.dat missing)"
        )

    for n in nuclide_list:
        try:
            result[n] = db.get(n)
        except KeyError:
            # NucleusNotFoundError subclasses KeyError; a few Database
            # subclasses still raise the plain KeyError.
            if not skip_missing:
                raise

    return result


def get_resonances(nuclei: list, wave: str = 's', skip_missing: bool = False) -> dict:
    """Get resonance parameters for multiple nuclei.

    Args:
        nuclei: List of nuclei. Each element can be:
                - A Nuclide/Nucleus object
                - A tuple (Z, A)
                - A list [Z, A]
        wave: Resonance type - 's' (s-wave) or 'p' (p-wave). Default: 's'
        skip_missing: If True, skip nuclei not in the database.
                      If False (default), raise NucleusNotFoundError.

    Returns:
        Dictionary mapping Nuclide -> Entry object

    Examples:
        >>> nuclei = [(82, 208), (50, 120)]
        >>> resonances = riplpy.get_resonances(nuclei, wave='s')

    Raises:
        ValueError: If the wave type is not recognized
        NucleusNotFoundError: If a nucleus is not found and skip_missing=False
    """
    from riplpy.exceptions import NucleusNotFoundError

    if wave == 's':
        db = resonances.db.swave
    elif wave == 'p':
        db = resonances.db.pwave
    else:
        raise ValueError(f"Unknown wave type '{wave}'. Available: 's', 'p'")

    nuclide_list = _normalize_nuclei_list(nuclei)
    result = {}

    for n in nuclide_list:
        try:
            result[n] = db.get(n)
        except KeyError:
            # NucleusNotFoundError subclasses KeyError; a few Database
            # subclasses still raise the plain KeyError.
            if not skip_missing:
                raise

    return result


def get_fission_barriers(nuclei: list, model: str = 'empirical', skip_missing: bool = False) -> dict:
    """Get fission barrier parameters for multiple nuclei.

    Args:
        nuclei: List of nuclei. Each element can be:
                - A Nuclide/Nucleus object
                - A tuple (Z, A)
                - A list [Z, A]
        model: Barrier model. Options: 'empirical', 'empirical_new', 'empire',
               'hfb', 'bskg3'. Default: 'empirical'
        skip_missing: If True, skip nuclei not in the database.
                      If False (default), raise NucleusNotFoundError.

    Returns:
        Dictionary mapping Nuclide -> Entry object

    Examples:
        >>> nuclei = [(92, 235), (92, 238), (94, 239)]
        >>> barriers = riplpy.get_fission_barriers(nuclei, model='hfb')

    Raises:
        ValueError: If the model is not recognized
        NucleusNotFoundError: If a nucleus is not found and skip_missing=False
    """
    from riplpy.exceptions import NucleusNotFoundError

    if model not in fission.BARRIER_MODELS:
        raise ValueError(f"Unknown fission barrier model '{model}'. Available: {list(fission.BARRIER_MODELS.keys())}")

    db = getattr(fission.db, fission.BARRIER_MODELS[model])

    nuclide_list = _normalize_nuclei_list(nuclei)
    result = {}

    if db is None:
        if skip_missing:
            return result
        raise NucleusNotFoundError(
            f"Fission barrier model '{model}' database not loaded (data file missing)"
        )

    for n in nuclide_list:
        try:
            result[n] = db.get(n)
        except KeyError:
            # NucleusNotFoundError subclasses KeyError; a few Database
            # subclasses still raise the plain KeyError.
            if not skip_missing:
                raise

    return result


def get_omp(iref: int) -> object:
    """Get an optical model potential by its reference number.

    Args:
        iref: Reference number of the potential (1-2500+)

    Returns:
        SphericalOMP or CoupledChannelOMP object

    Examples:
        >>> pot = riplpy.get_omp(2405)
        >>> print(pot.projectile)
        'n'
        >>> print(pot.header.author)
        'Koning'

    Raises:
        KeyError: If the potential is not found

    Note:
        The optical database must be loaded first via riplpy.load()
    """
    return optical.db.potentials.get(iref)


def list_omps(projectile: str = None) -> list:
    """List all optical model potential reference numbers.

    Args:
        projectile: Optional projectile filter ('n', 'p', 'd', 't', 'h', 'a').
                    If None, returns all potentials.

    Returns:
        List of iref values (integers)

    Examples:
        >>> all_irefs = riplpy.list_omps()
        >>> neutron_irefs = riplpy.list_omps(projectile='n')

    Note:
        The optical database must be loaded first via riplpy.load()
    """
    db = optical.db.potentials
    if projectile is not None:
        db = db.filter_by_projectile(projectile)
    return list(db.data.keys())


def find_omp(projectile: str, Z: int, A: int, E: float = None) -> list:
    """Find optical model potentials suitable for a specific reaction.

    Args:
        projectile: Projectile type ('n', 'p', 'd', 't', 'h', 'a')
        Z: Target atomic number
        A: Target mass number
        E: Optional projectile energy in MeV. If provided, filters to
           potentials valid at this energy.

    Returns:
        List of matching OMP objects, sorted by relevance

    Examples:
        >>> # Find neutron potentials for Pb-208
        >>> pots = riplpy.find_omp('n', 82, 208)
        >>> print(len(pots))
        15

        >>> # Find potentials valid at 14 MeV
        >>> pots = riplpy.find_omp('n', 82, 208, E=14.0)

    Note:
        The optical database must be loaded first via riplpy.load()
    """
    return optical.db.potentials.find_for_reaction(projectile, Z, A, E)


def get_deformation(Z: int, A: int, Ex: float = 0.0, L: int = 2) -> object:
    """Get nuclear deformation parameter.

    Args:
        Z: Atomic number
        A: Mass number
        Ex: Excitation energy in MeV (default 0.0 for low-lying states)
        L: Deformation order (2 for quadrupole, 3 for octupole). Default: 2

    Returns:
        DeformationEntry object with beta value and metadata

    Examples:
        >>> # Get quadrupole deformation for Pu-239
        >>> deform = riplpy.get_deformation(94, 239)
        >>> print(deform.beta)
        0.3699

        >>> # Get excited state deformation for U-238
        >>> deform = riplpy.get_deformation(92, 238, Ex=0.045)

    Raises:
        KeyError: If the deformation is not found

    Note:
        The optical database must be loaded first via riplpy.load().
        The database contains excited-level deformations derived from
        experiments. Use get_for_nucleus() for all available entries.
    """
    return optical.db.deformations.get(Z, A, Ex, L)


def get_omp_reference(ref_num: int) -> object:
    """Get bibliographic reference for an optical model potential.

    Args:
        ref_num: Reference number

    Returns:
        Reference object with citation text

    Examples:
        >>> ref = riplpy.get_omp_reference(100)
        >>> print(ref.citation)
        'A.J.Koning and J.P.Delaroche, Nucl.Phys. A713 (2003) 231'
        >>> print(ref.year)
        2003

    Raises:
        KeyError: If the reference is not found

    Note:
        The optical database must be loaded first via riplpy.load()
    """
    return optical.db.references.get(ref_num)
