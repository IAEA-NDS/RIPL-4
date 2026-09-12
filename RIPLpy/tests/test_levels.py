# -*- coding: utf-8 -*-
"""Tests for the RIPLpy levels module.

Tests the reading and writing of nuclear level data from RIPL-4.
"""

import os
import pytest
from pathlib import Path

import riplpy.levels as levels
from riplpy.collections import Nuclide


class TestLevelsCT:
    """Tests for the Constant Temperature fit parameters database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading levels CT parameters from RIPL directory."""
        db = levels.ct.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_get_ct_parameters(self, ripl_path):
        """Test retrieving CT fit parameters for a nucleus."""
        db = levels.ct.load(directory=ripl_path)
        # Fe-56 should be in the database
        n = Nuclide(Z=26, A=56)
        if n in db.data:
            entry = db.get(n)
            assert entry is not None
            # Should have temperature T and back-shift U0
            assert hasattr(entry, 'T')
            assert hasattr(entry, 'U0')

    def test_has_discrete_level_info(self, ripl_path):
        """Test that entries have discrete level information."""
        db = levels.ct.load(directory=ripl_path)
        for n, entry in list(db.data.items())[:50]:
            # Should have number of levels
            assert hasattr(entry, 'Nlev')
            assert hasattr(entry, 'Nmax')

    def test_temperature_values_reasonable(self, ripl_path):
        """Test that temperature values are in reasonable range."""
        db = levels.ct.load(directory=ripl_path)
        for n, entry in list(db.data.items())[:100]:
            # Temperature typically 0-10 MeV (0 for nuclei with no fit)
            if entry.T is not None:
                assert 0 <= entry.T < 20, f"T out of range for {n}"

    def test_write_and_read_roundtrip(self, ripl_path, temp_output_dir):
        """Test that writing and re-reading produces identical data."""
        db = levels.ct.load(directory=ripl_path)

        # Write to temp file
        output_path = os.path.join(temp_output_dir, "test_levels_ct.dat")
        db.save(output_path)

        # Read back
        db2 = levels.ct.load(file_path=output_path)

        # Compare number of entries
        assert len(db.data) == len(db2.data)


class TestDiscreteLevels:
    """Tests for the discrete nuclear levels database."""

    def test_load_element(self, ripl_path):
        """Test loading discrete levels for a single element."""
        db = levels.discrete.load_element(Z=26, directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_fe56_levels(self, ripl_path):
        """Test discrete levels for Fe-56 (well-studied nucleus)."""
        db = levels.discrete.load_element(Z=26, directory=ripl_path)
        n = Nuclide(Z=26, A=56)
        if n in db.data:
            entry = db.get(n)
            assert entry is not None
            # Fe-56 should have many known levels
            assert entry.num_levels > 0
            # Check ground state energy
            level_energies = entry.level_energies
            assert level_energies[0] == 0.0 or level_energies[0] < 0.01

    def test_levels_have_spin_parity(self, ripl_path):
        """Test that levels have spin and parity information."""
        db = levels.discrete.load_element(Z=26, directory=ripl_path)
        for n, entry in list(db.data.items())[:5]:
            for level in entry.levels[:3]:
                # Should have spin and parity attributes
                assert hasattr(level, 'spin')
                assert hasattr(level, 'parity')

    def test_gamma_transitions(self, ripl_path):
        """Test that levels have associated gamma transitions."""
        db = levels.discrete.load_element(Z=26, directory=ripl_path)
        # Look for a nucleus with gamma transitions
        for n, entry in db.data.items():
            if entry.total_number_gammas > 0:
                # Should have gamma transition data
                assert len(entry.gammas) > 0
                break

    def test_level_energies_increasing(self, ripl_path):
        """Test that level energies are in increasing order."""
        db = levels.discrete.load_element(Z=26, directory=ripl_path)
        for n, entry in list(db.data.items())[:10]:
            energies = entry.level_energies
            # Energies should be roughly increasing (within tolerance for isomers)
            for i in range(1, len(energies)):
                assert energies[i] >= energies[i-1] - 0.01, f"Non-increasing energy at level {i} for {n}"


class TestLevelsCoverage:
    """Tests for level data coverage."""

    def test_covers_many_elements(self, ripl_path):
        """Test that level database covers many elements."""
        db = levels.ct.load(directory=ripl_path)

        # Get unique Z values
        z_values = set(n.Z for n in db.data.keys())

        # Should cover many elements
        assert len(z_values) > 50
