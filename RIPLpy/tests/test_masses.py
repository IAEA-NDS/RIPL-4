# -*- coding: utf-8 -*-
"""Tests for the RIPLpy masses module.

Tests the reading and writing of nuclear mass data from RIPL-4.
"""

import os
import pytest
from pathlib import Path

import riplpy.masses as masses
from riplpy.collections import Nuclide


class TestAME20:
    """Tests for the AME2020 atomic mass evaluation database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading AME20 database from RIPL directory."""
        db = masses.ame20.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_load_from_file_path(self, ripl_path):
        """Test loading AME20 database from explicit file path."""
        file_path = os.path.join(ripl_path, "masses", "mass-ame20.dat")
        db = masses.ame20.load(file_path=file_path)
        assert db is not None
        assert len(db.data) > 0

    def test_get_entry(self, ripl_path):
        """Test retrieving a specific entry from the database."""
        db = masses.ame20.load(directory=ripl_path)
        # Pb-208 is a well-known stable nucleus
        n = Nuclide(Z=82, A=208)
        entry = db.get(n)
        assert entry is not None
        assert entry.n == n
        # Pb-208 should have experimental mass
        assert entry.Mexp is not None

    def test_nucleus_not_found(self, ripl_path):
        """Test that missing nucleus raises appropriate error."""
        db = masses.ame20.load(directory=ripl_path)
        # Non-existent nucleus (Z=200, A=500)
        n = Nuclide(Z=200, A=500)
        with pytest.raises(Exception):  # NucleusNotFoundError or KeyError
            db.get(n)

    def test_write_and_read_roundtrip(self, ripl_path, temp_output_dir):
        """Test that writing and re-reading produces identical data."""
        db = masses.ame20.load(directory=ripl_path)

        # Write to temp file
        output_path = os.path.join(temp_output_dir, "test_ame20.dat")
        db.save(output_path)

        # Read back
        db2 = masses.ame20.load(file_path=output_path)

        # Compare a few entries
        for n in list(db.data.keys())[:10]:
            assert db.get(n).Mexp == pytest.approx(db2.get(n).Mexp, rel=1e-6)


class TestBSkG3:
    """Tests for the BSkG3 mass model database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading BSkG3 database from RIPL directory."""
        db = masses.bskg3.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_get_entry_with_deformation(self, ripl_path):
        """Test retrieving entry with deformation parameters."""
        db = masses.bskg3.load(directory=ripl_path)
        # U-238 is a deformed nucleus
        n = Nuclide(Z=92, A=238)
        if n in db.data:
            entry = db.get(n)
            assert entry is not None
            # BSkG3 should have deformation parameters
            assert hasattr(entry, 'beta20') or 'beta20' in entry.__dict__


class TestHFB27:
    """Tests for the HFB-27 mass model database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading HFB27 database from RIPL directory."""
        db = masses.hfb27.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0


class TestFRDM12:
    """Tests for the FRDM2012 mass model database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading FRDM12 database from RIPL directory."""
        db = masses.frdm12.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0


class TestFRDM95:
    """Tests for the FRDM1995 mass model database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading FRDM95 database from RIPL directory.

        ``mass-frdm95.dat`` is not part of the current RIPL-4 GitHub release.
        We skip when the file is absent so the reader/writer code remains
        validated for users that still have the full release.
        """
        file_path = os.path.join(ripl_path, "masses", "mass-frdm95.dat")
        if not os.path.exists(file_path):
            pytest.skip("mass-frdm95.dat not present in this RIPL release")
        db = masses.frdm95.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0


class TestD1M:
    """Tests for the D1M mass model database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading D1M database from RIPL directory."""
        db = masses.d1m.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0


class TestWS4:
    """Tests for the WS4 mass model database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading WS4 database from RIPL directory."""
        db = masses.ws4.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0


class TestHFB14:
    """Tests for the HFB-14 mass model database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading HFB14 database from RIPL directory.

        ``mass-hfb14.dat`` is not part of the current RIPL-4 GitHub release.
        We skip when the file is absent so the reader/writer code remains
        validated for users that still have the full release.
        """
        file_path = os.path.join(ripl_path, "masses", "mass-hfb14.dat")
        if not os.path.exists(file_path):
            pytest.skip("mass-hfb14.dat not present in this RIPL release")
        db = masses.hfb14.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0


class TestDeformations:
    """Tests for the ground-state deformations database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading deformations database from RIPL directory."""
        db = masses.deformations.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_get_deformation(self, ripl_path):
        """Test retrieving deformation for a known deformed nucleus."""
        db = masses.deformations.load(directory=ripl_path)
        # U-238 is a well-known deformed nucleus
        n = Nuclide(Z=92, A=238)
        if n in db.data:
            entry = db.get(n)
            assert entry.beta2 is not None
            # U-238 should have a positive quadrupole deformation
            assert entry.beta2 > 0


class TestAbundances:
    """Tests for the natural abundances database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading abundances database from RIPL directory."""
        db = masses.ab.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_get_hydrogen_abundance(self, ripl_path):
        """Test retrieving H-1 abundance (most abundant isotope)."""
        db = masses.ab.load(directory=ripl_path)
        n = Nuclide(Z=1, A=1)
        entry = db.get(n)
        assert entry is not None
        # H-1 is about 99.98% abundant
        assert entry.abundance > 99.0

    def test_abundance_sum_for_element(self, ripl_path):
        """Test that abundances for an element sum to approximately 100%."""
        db = masses.ab.load(directory=ripl_path)
        # Get all carbon isotopes (Z=6)
        carbon_isotopes = [n for n in db.data.keys() if n.Z == 6]
        total = sum(db.get(n).abundance for n in carbon_isotopes)
        # Should sum to 100% within small tolerance
        assert 99.9 < total < 100.1

    def test_all_abundances_positive(self, ripl_path):
        """Test that all abundances are non-negative."""
        db = masses.ab.load(directory=ripl_path)
        for n, entry in db.data.items():
            assert entry.abundance >= 0, f"Negative abundance for {n}"

    def test_write_and_read_roundtrip(self, ripl_path, temp_output_dir):
        """Test that writing and re-reading produces identical data."""
        db = masses.ab.load(directory=ripl_path)

        # Write to temp file
        output_path = os.path.join(temp_output_dir, "test_abundance.dat")
        db.save(output_path)

        # Read back
        db2 = masses.ab.load(file_path=output_path)

        # Compare entries
        for n in list(db.data.keys())[:20]:
            orig = db.get(n)
            new = db2.get(n)
            assert orig.abundance == pytest.approx(new.abundance, rel=1e-4)


class TestMassComparison:
    """Tests comparing different mass models."""

    def test_mass_models_have_common_nuclei(self, ripl_path):
        """Test that different mass models share common nuclei."""
        ame = masses.ame20.load(directory=ripl_path)
        frdm = masses.frdm12.load(directory=ripl_path)

        # Find common nuclei
        common = set(ame.data.keys()) & set(frdm.data.keys())
        assert len(common) > 100  # Should have many common nuclei

    def test_mass_values_reasonable(self, ripl_path):
        """Test that mass excess values are in reasonable range."""
        db = masses.ame20.load(directory=ripl_path)

        for n, entry in list(db.data.items())[:100]:
            # Mass excess should generally be between -100 and +100 MeV
            # (with some exceptions for very light/heavy nuclei)
            if entry.Mexp is not None:
                assert -200 < entry.Mexp < 200, f"Mass excess out of range for {n}"
