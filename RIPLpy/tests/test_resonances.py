# -*- coding: utf-8 -*-
"""Tests for the RIPLpy resonances module.

Tests the reading and writing of neutron resonance data from RIPL-4.
"""

import os
import pytest
from pathlib import Path

import riplpy.resonances as resonances
from riplpy.collections import Nuclide


class TestSwaveResonances:
    """Tests for the S-wave neutron resonance database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading S-wave resonances from RIPL directory."""
        db = resonances.swave.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_get_resonance_spacing(self, ripl_path):
        """Test retrieving resonance spacing D0."""
        db = resonances.swave.load(directory=ripl_path)
        # Na-23 should be in the database
        n = Nuclide(Z=11, A=23)
        if n in db.data:
            entry = db.get(n)
            assert entry is not None
            # Should have resonance spacing D
            assert entry.D is not None

    def test_resonance_values_reasonable(self, ripl_path):
        """Test that resonance parameters are in reasonable range."""
        db = resonances.swave.load(directory=ripl_path)
        for n, entry in list(db.data.items())[:50]:
            # Binding energy typically 5-12 MeV
            if entry.Bn is not None:
                assert 0 < entry.Bn < 20, f"Bn out of range for {n}"
            # Resonance spacing positive
            if entry.D is not None:
                assert entry.D > 0, f"D should be positive for {n}"

    def test_write_and_read_roundtrip(self, ripl_path, temp_output_dir):
        """Test that writing and re-reading produces identical data."""
        db = resonances.swave.load(directory=ripl_path)

        # Write to temp file
        output_path = os.path.join(temp_output_dir, "test_swave.dat")
        db.save(output_path)

        # Read back
        db2 = resonances.swave.load(file_path=output_path)

        # Compare number of entries
        assert len(db.data) == len(db2.data)

        # Spot-check a representative nuclide for value-level fidelity
        n = Nuclide(Z=11, A=23)
        if n in db.data and n in db2.data:
            e1 = db.get(n)
            e2 = db2.get(n)
            assert e1.D == e2.D
            assert e1.Bn == e2.Bn
            assert e1.D_BNL == e2.D_BNL


class TestPwaveResonances:
    """Tests for the P-wave neutron resonance database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading P-wave resonances from RIPL directory."""
        db = resonances.pwave.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_pwave_has_fewer_entries(self, ripl_path):
        """Test that P-wave has fewer entries than S-wave (as expected)."""
        swave_db = resonances.swave.load(directory=ripl_path)
        pwave_db = resonances.pwave.load(directory=ripl_path)

        # P-wave data is typically more limited
        assert len(pwave_db.data) < len(swave_db.data)


class TestResonanceComparison:
    """Tests comparing S-wave and P-wave resonances."""

    def test_swave_and_pwave_have_common_nuclei(self, ripl_path):
        """Test that S-wave and P-wave share some common nuclei."""
        swave = resonances.swave.load(directory=ripl_path)
        pwave = resonances.pwave.load(directory=ripl_path)

        # Find common nuclei
        common = set(swave.data.keys()) & set(pwave.data.keys())
        # Should have some overlap
        assert len(common) > 0
