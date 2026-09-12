# -*- coding: utf-8 -*-
"""Tests for the RIPLpy densities module.

Tests the reading and writing of nuclear level density data from RIPL-4.
"""

import os
import pytest
from pathlib import Path

import riplpy.densities as densities
from riplpy.collections import Nuclide
from riplpy.exceptions import RiplFileNotFoundError


def _bsfg_available(ripl_path):
    return os.path.exists(os.path.join(ripl_path, 'densities', 'level-densities-bfmeff.dat'))


def _ct_available(ripl_path):
    return os.path.exists(os.path.join(ripl_path, 'densities', 'level-densities-ctmeff.dat'))


def _hfb_available(ripl_path):
    return os.path.isdir(os.path.join(ripl_path, 'densities', 'total', 'level-densities-hfb'))


def _mk_available(ripl_path):
    return os.path.exists(os.path.join(ripl_path, 'densities', 'shellmengoninakajima.dat'))


class TestBSFG:
    """Tests for the Back-Shifted Fermi Gas level density database.

    BSFG data is not shipped in the public RIPL-4 github release, so these
    tests skip cleanly when the backing file is absent.
    """

    def test_load_from_directory(self, ripl_path):
        if not _bsfg_available(ripl_path):
            pytest.skip("BSFG data file (level-densities-bfmeff.dat) not present")
        db = densities.bsfg.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_get_entry(self, ripl_path):
        if not _bsfg_available(ripl_path):
            pytest.skip("BSFG data file (level-densities-bfmeff.dat) not present")
        db = densities.bsfg.load(directory=ripl_path)
        n = Nuclide(Z=26, A=56)
        if n in db.data:
            entry = db.get(n)
            assert entry is not None
            assert entry.n == n
            assert entry.ainf is not None

    def test_write_and_read_roundtrip(self, ripl_path, temp_output_dir):
        if not _bsfg_available(ripl_path):
            pytest.skip("BSFG data file (level-densities-bfmeff.dat) not present")
        db = densities.bsfg.load(directory=ripl_path)
        output_path = os.path.join(temp_output_dir, "test_bsfg.dat")
        db.save(output_path)
        db2 = densities.bsfg.load(file_path=output_path)
        for n in list(db.data.keys())[:10]:
            orig = db.get(n)
            new = db2.get(n)
            if orig.ainf is not None and new.ainf is not None:
                assert orig.ainf == pytest.approx(new.ainf, rel=1e-4)


class TestCT:
    """Tests for the Constant Temperature level density database.

    CT data is not shipped in the public RIPL-4 github release; these tests
    skip cleanly when the backing file is absent.
    """

    def test_load_from_directory(self, ripl_path):
        if not _ct_available(ripl_path):
            pytest.skip("CT data file (level-densities-ctmeff.dat) not present")
        db = densities.ct.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_has_temperature_parameters(self, ripl_path):
        if not _ct_available(ripl_path):
            pytest.skip("CT data file (level-densities-ctmeff.dat) not present")
        db = densities.ct.load(directory=ripl_path)
        for n, entry in list(db.data.items())[:10]:
            assert hasattr(entry, 'T')
            assert hasattr(entry, 'E0')


class TestEGSM:
    """Tests for the Enhanced Generalized Superfluid Model database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading EGSM database from RIPL directory."""
        db = densities.egsm.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_get_level_density_parameter(self, ripl_path):
        """Test retrieving level density parameter 'a'."""
        db = densities.egsm.load(directory=ripl_path)
        for n, entry in list(db.data.items())[:10]:
            # EGSM should have level density parameter 'a'
            assert hasattr(entry, 'a')


class TestEGSMNorm:
    """Tests for the EGSM normalization factors database."""

    def test_load_from_directory(self, ripl_path):
        """Test loading EGSM normalization database from RIPL directory."""
        db = densities.egsm_norm.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_normalization_factors_reasonable(self, ripl_path):
        """Test that normalization factors are close to 1.0."""
        db = densities.egsm_norm.load(directory=ripl_path)
        for ele, entry in db.data.items():
            # Normalization factors should be close to 1.0 (within ~20%)
            assert 0.5 < entry.factor < 1.5, f"Factor out of range for element {ele}"


class TestShellCorrections:
    """Tests for the shell correction databases."""

    def test_load_ms_from_directory(self, ripl_path):
        """Test loading Myers-Swiatecki shell corrections from RIPL directory."""
        db = densities.shell_corr.load_ms(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_load_mk_from_directory(self, ripl_path):
        """Test loading Mengoni-Nakajima shell corrections from RIPL directory."""
        if not _mk_available(ripl_path):
            pytest.skip("MK shell corrections (shellmengoninakajima.dat) not present")
        db = densities.shell_corr.load_mk(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_ms_shell_corrections_have_deformations(self, ripl_path):
        """Test that MS shell corrections include deformation parameters."""
        db = densities.shell_corr.load_ms(directory=ripl_path)
        for n, entry in list(db.data.items())[:10]:
            assert hasattr(entry, 'shell')
            assert hasattr(entry, 'beta2')
            assert hasattr(entry, 'beta4')

    def test_shell_corrections_reasonable_values(self, ripl_path):
        """Test that shell correction values are in reasonable range."""
        db = densities.shell_corr.load_ms(directory=ripl_path)
        for n, entry in list(db.data.items())[:50]:
            # Shell corrections typically -10 to +10 MeV
            assert -20 < entry.shell < 20, f"Shell correction out of range for {n}"


class TestHFBLevelDensities:
    """Tests for the HFB spin-dependent level density database.

    The level-densities-hfb directory is not in the public github release;
    skip when absent.
    """

    def test_load_element(self, ripl_path):
        if not _hfb_available(ripl_path):
            pytest.skip("HFB level-densities-hfb directory not present")
        db = densities.hfb.load_element(Z=26, directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_has_positive_and_negative_parity(self, ripl_path):
        if not _hfb_available(ripl_path):
            pytest.skip("HFB level-densities-hfb directory not present")
        db = densities.hfb.load_element(Z=26, directory=ripl_path)
        for n, entry in list(db.data.items())[:5]:
            assert entry.positive_parity is not None or entry.negative_parity is not None
            if entry.positive_parity:
                assert len(entry.positive_parity.U) > 0
                assert len(entry.positive_parity.rho_tot) > 0

    def test_get_total_density(self, ripl_path):
        if not _hfb_available(ripl_path):
            pytest.skip("HFB level-densities-hfb directory not present")
        db = densities.hfb.load_element(Z=26, directory=ripl_path)
        n = Nuclide(Z=26, A=56)
        if n in db.data:
            entry = db.get(n)
            density = entry.get_total_density(5.0)
            assert density > 0

    def test_spin_dependent_densities(self, ripl_path):
        if not _hfb_available(ripl_path):
            pytest.skip("HFB level-densities-hfb directory not present")
        db = densities.hfb.load_element(Z=26, directory=ripl_path)
        for n, entry in list(db.data.items())[:3]:
            if entry.positive_parity:
                assert len(entry.positive_parity.rho_J) > 0
                assert len(entry.positive_parity.rho_J[0]) > 0


class TestCombinatorialLevelDensities:
    """Tests for the microscopic combinatorial level density databases
    (BSk14, BSkG3, QRPA-BE, T-HFB) shipped under densities/total/."""

    @pytest.mark.parametrize("module_name,Z", [
        ('bsk14_comb', 8),
        ('bskg3_comb', 8),
        ('qrpabe', 10),
        ('thfb_comb', 8),
    ])
    def test_load_element(self, ripl_path, module_name, Z):
        """Each model should load a single element directory of .tab files."""
        module = getattr(densities, module_name)
        db = module.load_element(Z=Z, directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0
        # Confirm we got entries for the requested Z
        assert any(n.Z == Z for n in db.data)
        # Verify spin-dependent data is parsed
        entry = next(iter(db.data.values()))
        assert entry.positive_parity is not None or entry.negative_parity is not None

    @pytest.mark.parametrize("module_name", [
        'bsk14_comb', 'bskg3_comb', 'qrpabe', 'thfb_comb',
    ])
    def test_load_all(self, ripl_path, module_name):
        """Each model should load thousands of nuclei across all elements."""
        module = getattr(densities, module_name)
        db = module.load_all(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 100, f"{module_name} loaded only {len(db.data)} nuclei"


class TestDensityComparison:
    """Tests comparing different level density models."""

    def test_bsfg_and_ct_have_common_nuclei(self, ripl_path):
        """Test that BSFG and CT share common nuclei."""
        if not (_bsfg_available(ripl_path) and _ct_available(ripl_path)):
            pytest.skip("BSFG/CT data files not present in this RIPL release")
        bsfg = densities.bsfg.load(directory=ripl_path)
        ct = densities.ct.load(directory=ripl_path)
        common = set(bsfg.data.keys()) & set(ct.data.keys())
        assert len(common) > 50


class TestDensitiesSectionLoad:
    """Smoke test for ``riplpy.densities.load``."""

    def test_section_load_does_not_raise(self, ripl_path):
        """The section-level loader should tolerate missing legacy files."""
        densities.load(ripl_path)
        # EGSM and the combinatorial databases should always be present
        assert densities.db.egsm is not None
        assert densities.db.bsk14_comb is not None
        assert densities.db.bskg3_comb is not None
        assert densities.db.qrpabe is not None
        assert densities.db.thfb_comb is not None
