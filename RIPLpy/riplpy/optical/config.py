# -*- coding: utf-8 -*-
"""Configuration and constants for the optical model section.

This module defines file paths, projectile types, model types, and other
constants used throughout the optical model parsing code.
"""

import os as _os

from riplpy import config as _riplpy_config

# =============================================================================
# File Paths (relative to RIPL directory)
# =============================================================================

# Main data directory
DATA_DIR = _os.path.join('optical', 'om-data')

# Data files
FILE_PATHS = {
    'parameters': _os.path.join(DATA_DIR, 'om-parameter-u.dat'),
    'index': _os.path.join(DATA_DIR, 'om-index.txt'),
    'index_by_z': _os.path.join(DATA_DIR, 'om-index-by-Z.txt'),
    'references': _os.path.join(DATA_DIR, 'om-references.txt'),
    'deformations': _os.path.join(DATA_DIR, 'om-deformations.dat'),
    'mod_potentials': _os.path.join(DATA_DIR, 'mod-potentials'),
    'rop2013': _os.path.join(DATA_DIR, 'ROP2013za.dat'),
    'atomki': _os.path.join('optical', 'atomki'),
}

# Entry separator in the main parameter file
ENTRY_SEPARATOR = '+' * 80


def get_data_file_path(key: str) -> str:
    """Get the relative path to a data file.

    Args:
        key: File identifier (e.g., 'parameters', 'index', 'references')

    Returns:
        Relative file path from RIPL directory

    Raises:
        KeyError: If the key is not recognized
    """
    if key not in FILE_PATHS:
        raise KeyError(f"Unknown file key '{key}'. Available: {list(FILE_PATHS.keys())}")
    return FILE_PATHS[key]


# =============================================================================
# Projectile Types
# =============================================================================

# Map of (Z, A) to projectile name
PROJECTILES = {
    (0, 1): 'n',      # neutron
    (1, 1): 'p',      # proton
    (1, 2): 'd',      # deuteron
    (1, 3): 't',      # triton
    (2, 3): 'He3',    # helion
    (2, 4): 'a',      # alpha
}

# Reverse map: name to (Z, A)
PROJECTILE_ZA = {v: k for k, v in PROJECTILES.items()}


def projectile_name(z: int, a: int) -> str:
    """Get the projectile name from Z and A.

    Args:
        z: Projectile atomic number
        a: Projectile mass number

    Returns:
        Projectile name (e.g., 'n', 'p', 'd', 't', 'He3', 'a')
    """
    return PROJECTILES.get((z, a), f'Z{z}A{a}')


# =============================================================================
# Model Types
# =============================================================================

# Model type codes
MODEL_SPHERICAL = 0
MODEL_RIGID_ROTOR = 1
MODEL_VIBRATIONAL = 2
MODEL_SOFT_ROTOR = 3
MODEL_RIGID_SOFT = 4
MODEL_SOFT_DEFORMED = 5

# Model type names
MODEL_NAMES = {
    MODEL_SPHERICAL: 'spherical',
    MODEL_RIGID_ROTOR: 'rigid_rotor',
    MODEL_VIBRATIONAL: 'vibrational',
    MODEL_SOFT_ROTOR: 'soft_rotor',
    MODEL_RIGID_SOFT: 'rigid_soft',
    MODEL_SOFT_DEFORMED: 'soft_deformed',
}

# Model type descriptions
MODEL_DESCRIPTIONS = {
    MODEL_SPHERICAL: 'Spherical optical model potential',
    MODEL_RIGID_ROTOR: 'Coupled-channel rotational model (rigid rotor)',
    MODEL_VIBRATIONAL: 'Coupled-channel vibrational model (harmonic oscillator)',
    MODEL_SOFT_ROTOR: 'Coupled-channel soft rotor (vibrational-rotational)',
    MODEL_RIGID_SOFT: 'Coupled-channel rigid multiband model',
    MODEL_SOFT_DEFORMED: 'Coupled-channel soft multiband model',
}


def model_name(imodel: int) -> str:
    """Get the model name from model code.

    Args:
        imodel: Model type code (0-5)

    Returns:
        Model name string
    """
    return MODEL_NAMES.get(imodel, f'unknown_{imodel}')


def is_coupled_channel(imodel: int) -> bool:
    """Check if model type is coupled-channel.

    Args:
        imodel: Model type code

    Returns:
        True if coupled-channel model
    """
    return imodel in (MODEL_RIGID_ROTOR, MODEL_VIBRATIONAL, MODEL_SOFT_ROTOR,
                      MODEL_RIGID_SOFT, MODEL_SOFT_DEFORMED)


# =============================================================================
# Potential Component Types
# =============================================================================

# Component indices (1-6 in Fortran, 0-5 in Python)
COMP_REAL_VOLUME = 0       # Real volume (Woods-Saxon)
COMP_IMAG_VOLUME = 1       # Imaginary volume (Woods-Saxon)
COMP_REAL_SURFACE = 2      # Real surface derivative
COMP_IMAG_SURFACE = 3      # Imaginary surface derivative
COMP_REAL_SPINORBIT = 4    # Real spin-orbit
COMP_IMAG_SPINORBIT = 5    # Imaginary spin-orbit

NUM_COMPONENTS = 6

COMPONENT_NAMES = {
    COMP_REAL_VOLUME: 'real_volume',
    COMP_IMAG_VOLUME: 'imag_volume',
    COMP_REAL_SURFACE: 'real_surface',
    COMP_IMAG_SURFACE: 'imag_surface',
    COMP_REAL_SPINORBIT: 'real_spinorbit',
    COMP_IMAG_SPINORBIT: 'imag_spinorbit',
}

COMPONENT_DESCRIPTIONS = {
    COMP_REAL_VOLUME: 'Real central volume potential (Woods-Saxon)',
    COMP_IMAG_VOLUME: 'Imaginary central volume potential (Woods-Saxon)',
    COMP_REAL_SURFACE: 'Real surface derivative potential',
    COMP_IMAG_SURFACE: 'Imaginary surface derivative potential',
    COMP_REAL_SPINORBIT: 'Real spin-orbit potential',
    COMP_IMAG_SPINORBIT: 'Imaginary spin-orbit potential',
}


# =============================================================================
# Coefficient Array Sizes
# =============================================================================

NUM_RADIUS_COEFFS = 13      # rco(i,j,k), k=1,13
NUM_DIFFUSENESS_COEFFS = 13  # aco(i,j,k), k=1,13
NUM_POTENTIAL_COEFFS = 25    # pot(i,j,k), k=1,25
NUM_COULOMB_PARAMS = 8       # Parameters per Coulomb energy range


# =============================================================================
# Dispersion Relation Flags
# =============================================================================

# idr values
DR_NONE = 0                  # No dispersion relations
DR_EQUIVALENT_VOLUME = 1     # Equivalent volume real potential
DR_EXACT_ANALYTIC = 2        # Exact analytic dispersion
DR_EXACT_ANALYTIC_SO = 3     # Exact analytic + spin-orbit
DR_EXACT_NUMERIC = -2        # Exact numeric dispersion
DR_EXACT_NUMERIC_SO = -3     # Exact numeric + spin-orbit


def uses_dispersion(idr: int) -> bool:
    """Check if potential uses dispersion relations.

    Args:
        idr: Dispersion relation flag

    Returns:
        True if dispersion relations are used
    """
    return idr != DR_NONE


# =============================================================================
# Relativistic Flags
# =============================================================================

# irel values
REL_NONRELATIVISTIC = 0      # Non-relativistic
REL_RELATIVISTIC = 1         # Relativistic kinematics
REL_RELATIVISTIC_SCALED = 2  # Relativistic + scaled potentials


def is_relativistic(irel: int) -> bool:
    """Check if potential uses relativistic kinematics.

    Args:
        irel: Relativistic flag

    Returns:
        True if relativistic
    """
    return irel != REL_NONRELATIVISTIC
