# -*- coding: utf-8 -*-
"""Tests that databases left empty by missing data files are reported.

A database whose data file is absent from the RIPL tree is loaded as empty
(or None) so that the rest of RIPL still loads. These tests check that
riplpy.load() says so with a warning shown by default, and that
riplpy.missing_databases() lists them, without listing the heavy tables
that load() skips on purpose.
"""

import warnings

import pytest

import riplpy
from riplpy.exceptions import MissingDataWarning


@pytest.fixture(scope="module")
def missing(ripl_path):
    """Load RIPL once, capturing the warning, and return (warning, list)."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        riplpy.load(directory=ripl_path)
    record = [w for w in caught if issubclass(w.category, MissingDataWarning)]
    return record, riplpy.missing_databases()


def test_load_warns_with_missing_names(missing):
    """The warning text names every database that missing_databases() lists."""
    record, names = missing
    if not names:
        assert record == []
        pytest.skip("This RIPL tree has every data file; nothing to report")
    assert len(record) == 1
    message = str(record[-1].message)
    for name in names:
        assert name in message


def test_missing_databases_are_empty(missing):
    """Every listed database is really None or empty."""
    _, names = missing
    for name in names:
        section, db_name = name.split('.')
        value = getattr(getattr(riplpy, section).db, db_name)
        assert value is None or len(value) == 0


def test_skipped_heavy_tables_not_listed(missing):
    """Tables that load() skips by default are not reported as missing."""
    _, names = missing
    for name in ('gamma.gsf', 'gamma.tlo', 'gamma.psf',
                 'fission.hfbpath_bskg3', 'fission.nld_fis',
                 'masses.density_bskg3'):
        assert name not in names


def test_loaded_databases_not_listed(missing):
    """Databases that loaded entries are not reported."""
    _, names = missing
    assert len(riplpy.masses.db.ame20) > 0
    assert 'masses.ame20' not in names
