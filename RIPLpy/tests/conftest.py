# -*- coding: utf-8 -*-
"""Pytest configuration and fixtures for RIPLpy tests.

The RIPL-4 data directory is resolved from the environment rather than a
hardcoded path, so the suite runs anywhere the data is available:

    RIPL_LOCATION=/path/to/RIPL-4/github python -m pytest

Resolution order: ``RIPLPY_TEST_PATH`` (test-only override), then
``RIPL_LOCATION`` (the library's standard variable), then whatever
``riplpy.config`` can auto-detect (e.g. ``~/.riplpyrc``). If none resolves to
an existing directory, data-backed tests skip with an explanatory message.
"""

import os
import pytest


def _resolve_ripl_path():
    """Return the RIPL-4 data directory for the session, or None."""
    for env in ("RIPLPY_TEST_PATH", "RIPL_LOCATION"):
        p = os.environ.get(env)
        if p and os.path.isdir(p):
            return p
    try:
        import riplpy.config as _config
        p = _config.get_path()
        if p and os.path.isdir(p):
            return p
    except Exception:
        pass
    return None


RIPL_PATH = _resolve_ripl_path()


@pytest.fixture(scope="session")
def ripl_path():
    """Return the resolved RIPL-4 data directory, or skip if unavailable."""
    if not RIPL_PATH:
        pytest.skip(
            "No RIPL-4 data directory found. Set RIPL_LOCATION (or "
            "RIPLPY_TEST_PATH) to a RIPL-4 tree to run data-backed tests."
        )
    return RIPL_PATH


@pytest.fixture(scope="session")
def temp_output_dir(tmp_path_factory):
    """Create a temporary directory for test output files."""
    return tmp_path_factory.mktemp("riplpy_test_output")


@pytest.fixture(scope="session", autouse=True)
def set_ripl_env():
    """Point riplpy.config at the resolved RIPL path for the test session."""
    if not RIPL_PATH:
        yield
        return

    old_value = os.environ.get("RIPL_LOCATION")
    os.environ["RIPL_LOCATION"] = RIPL_PATH
    try:
        import riplpy.config as _config
        _config.set_path(RIPL_PATH)
    except Exception:
        pass

    yield

    if old_value is not None:
        os.environ["RIPL_LOCATION"] = old_value
    else:
        os.environ.pop("RIPL_LOCATION", None)
