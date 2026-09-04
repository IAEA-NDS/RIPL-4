# -*- coding: utf-8 -*-
"""Optical Model Potential (OMP) section of RIPL.

This module provides access to the RIPL optical model parameter library,
containing ~584 optical model potentials for various projectiles and targets.

Submodules:
    config: Configuration and file paths
    core: Base dataclasses for OMP components
    index: Index file parser
    spherical: Spherical potential parser
    coupled_channel: Coupled-channel potential parser
    omp: Main OMP database
    deformations: Excited-level deformation parameters
    references: Bibliographic references
    rop2013: Avrigeanu revised alpha ROP table (ROP2013za.dat)
    atomki: Lazy reader for the per-nucleus ATOMKI alpha-OMP gnu files

Example:
    >>> import riplpy.optical as optical
    >>> optical.load(directory='/path/to/RIPL')
    >>> # Get a specific potential by reference number
    >>> pot = optical.db.potentials.get(2405)
    >>> print(pot.projectile)
    'n'
    >>> # Get deformation for Pu-239
    >>> deform = optical.db.deformations.get(94, 239)

    >>> # Load modified potentials
    >>> mod_db = optical.load_modified_potentials()
    >>> print(mod_db.irefs)  # [1480, 1481, 1482, 2412, 2415, 4609, 4610]

    >>> # Lazy ATOMKI per-nucleus reader (no preload)
    >>> entry = optical.atomki.load_nucleus(Z=26, A=42)
"""

# RIPLpy
from riplpy.db import _DbAccessor

# Local imports
from . import config
from . import core
from . import index
from . import omp
from . import spherical
from . import coupled_channel
from . import deformations
from . import references
from . import rop2013
from . import atomki

# Database accessor
db = _DbAccessor()


def load(directory: str = None) -> None:
    """Load the optical model database from the RIPL directory.

    Each underlying loader is wrapped in a try/except so that a missing or
    malformed data file emits a warning via :data:`riplpy.logger` instead of
    aborting the whole section. Subsequent attributes on :data:`optical.db`
    will simply be absent.

    Args:
        directory: Path to the RIPL database directory. If not provided,
                   uses the configured path from riplpy.config.

    Note:
        Once loaded, databases are accessible via optical.db:
            - optical.db.index: Index of all potentials (metadata)
            - optical.db.potentials: Full potential data (spherical + CC)
            - optical.db.deformations: Excited-level deformation parameters
            - optical.db.references: Bibliographic references
            - optical.db.rop2013: ROP2013za.dat (Avrigeanu alpha ROP table)

        The ATOMKI alpha-OMP gnu files are not preloaded; use
        :func:`riplpy.optical.atomki.load_nucleus` for on-demand access.
    """
    from riplpy import config as riplpy_config
    from riplpy import logger

    if directory is None:
        directory = riplpy_config.get_path()

    # Load the index
    try:
        db.index = index.load(directory=directory)
    except FileNotFoundError as exc:
        logger.warning("optical: index file missing, skipping: %s", exc)

    # Load the main OMP library
    try:
        db.potentials = omp.load(directory=directory)
    except FileNotFoundError as exc:
        logger.warning("optical: potential parameter file missing, skipping: %s", exc)

    # Load deformations
    try:
        db.deformations = deformations.load(directory=directory)
    except FileNotFoundError as exc:
        logger.warning("optical: deformations file missing, skipping: %s", exc)

    # Load references
    try:
        db.references = references.load(directory=directory)
    except FileNotFoundError as exc:
        logger.warning("optical: references file missing, skipping: %s", exc)

    # Load the ROP2013 (Avrigeanu revised alpha ROP) table
    try:
        db.rop2013 = rop2013.load(directory=directory)
    except FileNotFoundError as exc:
        logger.warning("optical: ROP2013za.dat missing, skipping: %s", exc)


def load_modified_potentials(directory: str = None):
    """Load modified optical model potentials from the mod-potentials directory.

    The mod-potentials directory contains individual OMP files (omp-XXXXXX.dat)
    that are updated or corrected versions of potentials in the main database.

    Args:
        directory: RIPL directory path. If None, uses configured path.

    Returns:
        OMPDatabase with modified potentials

    Example:
        >>> mod_db = optical.load_modified_potentials()
        >>> print(mod_db.irefs)
        [1480, 1481, 1482, 2412, 2415, 4609, 4610]
    """
    return omp.load_modified_potentials(directory=directory)


def load_with_modifications(directory: str = None, apply_modifications: bool = True):
    """Load the OMP database with optional modifications applied.

    This is a convenience function that loads the main database and optionally
    replaces entries with their modified versions from the mod-potentials directory.

    Args:
        directory: RIPL directory path. If None, uses configured path.
        apply_modifications: If True, replace potentials with modified versions
                            where available. Default True.

    Returns:
        OMPDatabase with potentials (modified versions applied if requested)

    Example:
        >>> db = optical.load_with_modifications()
        >>> # Potential 1480 will use the modified version if available
        >>> pot = db.get(1480)
    """
    return omp.load_with_modifications(directory=directory, apply_modifications=apply_modifications)


def find_for_reaction(projectile: str, Z: int, A: int, E: float = None):
    """Find OMPs valid for a given projectile + target (and optional energy).

    Convenience wrapper over :meth:`omp.OMPDatabase.find_for_reaction` that
    operates on the main library loaded into :data:`optical.db.potentials`.
    Use this to clearly select which incident particle the optical potential
    is for, e.g. ``find_for_reaction('alpha', 26, 56)`` for alpha + Fe-56.

    Args:
        projectile: Incident particle name. One of ``'n'``, ``'p'``, ``'d'``,
            ``'t'``, ``'He3'`` (or ``'3He'``), ``'a'`` (or ``'alpha'``).
        Z: Target atomic number.
        A: Target mass number.
        E: Projectile energy [MeV]. If ``None``, no energy filtering.

    Returns:
        List of OMP objects (``SphericalOMP`` / ``CoupledChannelOMP``) valid
        for this reaction.

    Raises:
        RuntimeError: If the main OMP library has not been loaded yet.

    Example:
        >>> import riplpy.optical as optical
        >>> optical.load()
        >>> omps = optical.find_for_reaction('alpha', Z=26, A=56, E=20.0)
        >>> [p.projectile for p in omps]
        ['a', 'a', ...]
    """
    if not hasattr(db, 'potentials'):
        raise RuntimeError(
            "Main OMP library not loaded. Call optical.load() first "
            "(or load a database with optical.omp.load())."
        )
    return db.potentials.find_for_reaction(projectile, Z=Z, A=A, E=E)
