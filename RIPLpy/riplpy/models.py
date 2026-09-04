# -*- coding: utf-8 -*-
"""Database model constants and enums for RIPLpy.

This module provides type-safe constants for referencing database models,
enabling IDE autocompletion and catching typos at development time.

Examples:
    >>> from riplpy import Models
    >>> riplpy.get_mass(82, 208, model=Models.Mass.AME20)
    >>> riplpy.get_level_density(82, 208, model=Models.Density.BSFG)
"""

from enum import Enum as _Enum


class MassModel(str, _Enum):
    """Mass model identifiers."""
    AME20 = 'ame20'
    FRDM12 = 'frdm12'
    FRDM2012 = 'frdm2012'
    FRDM95 = 'frdm95'
    FRDM1995 = 'frdm1995'
    HFB14 = 'hfb14'
    HFB27 = 'hfb27'
    BSKG3 = 'bskg3'
    D1M = 'd1m'


class DensityModel(str, _Enum):
    """Level density model identifiers."""
    BSFG = 'bsfg'
    CT = 'ct'
    EGSM = 'egsm'
    EGSM_NORM = 'egsm_norm'
    HFB = 'hfb'
    BSK14_COMB = 'bsk14_comb'
    BSKG3_COMB = 'bskg3_comb'
    QRPABE = 'qrpabe'
    THFB_COMB = 'thfb_comb'


class FissionModel(str, _Enum):
    """Fission barrier model identifiers."""
    EMPIRICAL = 'empirical'
    EMPIRICAL_NEW = 'empirical_new'
    EMPIRE = 'empire'
    HFB = 'hfb'
    BSKG3 = 'bskg3'
    D1M = 'd1m'


class ResonanceWave(str, _Enum):
    """Resonance wave types."""
    S = 's'
    P = 'p'


class Models:
    """Namespace for all model constants.

    Provides type-safe access to database model identifiers with
    IDE autocompletion support.

    Examples:
        >>> from riplpy import Models
        >>> riplpy.get_mass(82, 208, model=Models.Mass.AME20)
        >>> riplpy.get_level_density(26, 56, model=Models.Density.BSFG)
        >>> riplpy.get_fission_barrier(92, 238, model=Models.Fission.HFB)
        >>> riplpy.get_resonance(82, 208, wave=Models.Wave.S)
    """
    Mass = MassModel
    Density = DensityModel
    Fission = FissionModel
    Wave = ResonanceWave


__all__ = ['Models', 'MassModel', 'DensityModel', 'FissionModel', 'ResonanceWave']
