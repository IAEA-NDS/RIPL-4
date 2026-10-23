# -*- coding: utf-8 -*-
"""Tests that the Models enums match the model names the getters accept,
and that the package source compiles without escape-sequence warnings."""

import pathlib
import warnings

import pytest

import riplpy
from riplpy import Models
import riplpy.densities as densities
import riplpy.fission as fission
import riplpy.gamma as gamma
import riplpy.masses as masses


@pytest.mark.parametrize('enum, names', [
    (Models.Mass, masses.MASS_MODELS),
    (Models.Density, densities.LEVEL_DENSITY_MODELS),
    (Models.Fission, fission.BARRIER_MODELS),
    (Models.Gdr, ['auto'] + list(gamma.GDR_MODELS)),
])
def test_enum_matches_accepted_models(enum, names):
    """Each enum lists exactly the model names its getter accepts.

    Models.Mass used to lack WS4, although get_mass(model='ws4') works.
    """
    assert {m.value for m in enum} == set(names)


def test_ws4_member(ripl_path):
    masses.load(directory=ripl_path)
    assert riplpy.get_mass(82, 208, model=Models.Mass.WS4) == \
        riplpy.get_mass(82, 208, model='ws4')


def test_docstring_examples_use_ripl4_models():
    """The Models examples must run on RIPL-4, which lacks BSFG and the
    legacy HFB fission barriers."""
    import riplpy.models as models
    text = models.__doc__ + models.Models.__doc__
    assert 'Models.Density.BSFG' not in text
    assert 'Models.Fission.HFB)' not in text


def test_no_invalid_escape_sequences():
    """No module emits SyntaxWarning (an error in future Python versions)."""
    root = pathlib.Path(riplpy.__file__).parent
    problems = []
    for path in sorted(root.rglob('*.py')):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            compile(path.read_text(encoding='utf-8'), str(path), 'exec')
        problems += [f"{path.name}:{w.lineno}: {w.message}" for w in caught]
    assert problems == []
