# -*- coding: utf-8 -*-
"""Tests for the RIPLpy core functionality.

Tests the core collections (Nucleus, Nuclide, Element) and overall API.
"""

import os
import pytest

import riplpy
from riplpy.collections import Nucleus, Nuclide, Element
from riplpy.elements import Elements, ZtoSym, SymtoZ


class TestNucleus:
    """Tests for the Nucleus class."""

    def test_create_by_z_n(self):
        """Test creating nucleus by Z and N."""
        n = Nucleus(Z=82, N=126)
        assert n.Z == 82
        assert n.N == 126
        assert n.A == 208

    def test_create_default(self):
        """Test creating default placeholder nucleus."""
        n = Nucleus()
        assert n.Z is None
        assert n.N is None


class TestNuclide:
    """Tests for the Nuclide class."""

    def test_create_by_z_a(self):
        """Test creating nuclide by Z and A."""
        n = Nuclide(Z=82, A=208)
        assert n.Z == 82
        assert n.A == 208
        assert n.N == 126

    def test_create_by_z_n(self):
        """Test creating nuclide by Z and N."""
        n = Nuclide(Z=82, N=126)
        assert n.Z == 82
        assert n.N == 126
        assert n.A == 208

    def test_element_symbol(self):
        """Test getting element symbol."""
        n = Nuclide(Z=82, A=208)
        assert n.element_symbol == "Pb"

    def test_equality(self):
        """Test nuclide equality."""
        n1 = Nuclide(Z=82, A=208)
        n2 = Nuclide(Z=82, A=208)
        assert n1 == n2

    def test_hash(self):
        """Test nuclide hashing (for use as dict key)."""
        n1 = Nuclide(Z=82, A=208)
        n2 = Nuclide(Z=82, A=208)
        d = {n1: "test"}
        assert d[n2] == "test"

    def test_out_of_range_z_does_not_crash_formatting(self):
        """str()/repr of a Nuclide with no element symbol must not crash.

        Regression: formatting used self.element_symbol unguarded, so any
        Z > 118 raised a confusing bare KeyError from inside repr -- which in
        turn corrupted error messages of any caller that f-string-formatted
        the nuclide (e.g. NucleusNotFoundError construction in db.get).
        """
        from riplpy.exceptions import ElementNotFoundError
        n = Nuclide(Z=200, A=400)
        # str/repr must not raise
        s = str(n); r = repr(n)
        assert '200' in s and 'symbol' not in s
        assert '200' in r
        # element_symbol/name now raise ElementNotFoundError, not bare KeyError
        with pytest.raises(ElementNotFoundError):
            n.element_symbol
        with pytest.raises(ElementNotFoundError):
            n.element_name


class TestElement:
    """Tests for the Element class."""

    def test_create_by_z(self):
        """Test creating element by atomic number."""
        e = Element(Z=82)
        assert e.Z == 82
        assert e.symbol == "Pb"

    def test_create_by_symbol(self):
        """Test creating element by symbol."""
        e = Element(symbol="Pb")
        assert e.Z == 82
        assert e.symbol == "Pb"


class TestElementMappings:
    """Tests for element symbol/Z mappings."""

    def test_z_to_sym(self):
        """Test Z to symbol mapping."""
        assert ZtoSym[82] == "Pb"
        assert ZtoSym[26] == "Fe"
        assert ZtoSym[1] == "H"

    def test_sym_to_z(self):
        """Test symbol to Z mapping."""
        assert SymtoZ["Pb"] == 82
        assert SymtoZ["Fe"] == 26
        assert SymtoZ["H"] == 1

    def test_elements_class(self):
        """Test Elements class methods."""
        assert Elements.symbol(82) == "Pb"
        assert Elements.name(82) == "Lead"


class TestRIPLpyLoad:
    """Tests for the main RIPLpy load functionality."""

    def test_load_all_sections(self, ripl_path):
        """Test loading all RIPL sections at once."""
        riplpy.load(directory=ripl_path)

        # Check that databases are accessible
        assert riplpy.masses.db is not None
        assert riplpy.densities.db is not None
        assert riplpy.fission.db is not None
        assert riplpy.gamma.db is not None
        assert riplpy.levels.db is not None
        assert riplpy.resonances.db is not None

    def test_in_ripl_function(self, ripl_path):
        """Test the in_ripl helper function."""
        riplpy.load(directory=ripl_path)

        # Pb-208 should be in RIPL
        n = Nuclide(Z=82, A=208)
        assert riplpy.in_ripl(n) == True

    def test_in_sections_function(self, ripl_path):
        """Test the in_sections helper function."""
        riplpy.load(directory=ripl_path)

        # Common nucleus should be in multiple sections
        n = Nuclide(Z=26, A=56)  # Fe-56
        sections = riplpy.in_sections(n)
        assert len(sections) > 0

    def test_in_dbs_function(self, ripl_path):
        """Test the in_dbs helper function."""
        riplpy.load(directory=ripl_path)

        # Common nucleus should be in multiple databases
        n = Nuclide(Z=26, A=56)  # Fe-56
        dbs = riplpy.in_dbs(n)
        assert len(dbs) > 0
        # Each entry should be a tuple (section, db_name)
        for section, db_name in dbs:
            assert section in riplpy.sections


class TestVersion:
    """Tests for version information."""

    def test_version_defined(self):
        """Test that version is defined."""
        assert hasattr(riplpy, "__version__")
        assert riplpy.__version__ is not None
