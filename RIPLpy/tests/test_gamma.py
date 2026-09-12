# -*- coding: utf-8 -*-
"""Tests for the RIPLpy gamma module.

Tests the reading and writing of gamma-ray strength function data from RIPL-4.
"""

import os
import pytest
from pathlib import Path

import riplpy.gamma as gamma
from riplpy.collections import Nuclide


class TestGDRTheoretical:
    """Tests for the theoretical Giant Dipole Resonance parameters.

    The legacy ``gdr-parameters-theor.dat`` file is not shipped in the
    public github RIPL-4 release; tests are skipped when absent.
    """

    def test_load_from_directory(self, ripl_path):
        db = gamma.gdr.load(directory=ripl_path)
        assert db is not None
        if len(db.data) == 0:
            pytest.skip(
                "gdr-parameters-theor.dat not present in github RIPL-4 layout"
            )
        assert len(db.data) > 0

    def test_get_gdr_parameters(self, ripl_path):
        db = gamma.gdr.load(directory=ripl_path)
        if len(db.data) == 0:
            pytest.skip("gdr-parameters-theor.dat not present")
        n = Nuclide(Z=14, A=28)
        if n in db.data:
            entry = db.get(n)
            assert entry is not None
            assert entry.E1 is not None
            assert entry.W1 is not None

    def test_gdr_values_reasonable(self, ripl_path):
        db = gamma.gdr.load(directory=ripl_path)
        if len(db.data) == 0:
            pytest.skip("gdr-parameters-theor.dat not present")
        for n, entry in list(db.data.items())[:50]:
            if entry.E1 is not None:
                assert 5 < entry.E1 < 35, f"E1 out of range for {n}"
            if entry.E2 is not None:
                assert 5 < entry.E2 < 35, f"E2 out of range for {n}"
            if entry.W1 is not None:
                assert 0 < entry.W1 < 20, f"W1 out of range for {n}"

    def test_write_and_read_roundtrip(self, ripl_path, temp_output_dir):
        db = gamma.gdr.load(directory=ripl_path)
        if len(db.data) == 0:
            pytest.skip("gdr-parameters-theor.dat not present")
        output_path = os.path.join(temp_output_dir, "test_gdr.dat")
        db.save(output_path)
        # gdr.load() does not accept file_path; build a fresh DB and reload
        db2 = gamma.gdr.Database()
        db2.load(output_path)
        for n in list(db.data.keys())[:10]:
            orig = db.get(n)
            new = db2.get(n)
            if orig.E1 is not None and new.E1 is not None:
                assert orig.E1 == pytest.approx(new.E1, rel=1e-3)


class TestGammaStrengthFunction:
    """Tests for the gamma-ray strength function database.

    The legacy ``gamma/gamma-strength-micro`` directory is replaced in the
    github release by per-Z D1M+QRPA tables under ``gamma/d1m/``. The GSF
    loader automatically falls back to those tables when the legacy directory
    is absent, preserving the ``data[n]['U']`` / ``data[n]['fE1']`` schema.
    """

    def test_load_all(self, ripl_path):
        db = gamma.gsf.load_all(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_gsf_energy_grid(self, ripl_path):
        db = gamma.gsf.load_all(directory=ripl_path)
        for n in list(db.data.keys())[:5]:
            data = db.data[n]
            assert 'U' in data
            assert 'fE1' in data
            assert len(data['U']) > 100
            assert min(data['U']) < 1.0
            # D1M tables go up to 30 MeV; older legacy went past 20 MeV
            assert max(data['U']) > 20.0

    def test_gsf_values_positive(self, ripl_path):
        db = gamma.gsf.load_all(directory=ripl_path)
        for n in list(db.data.keys())[:10]:
            data = db.data[n]
            for fE1 in data['fE1']:
                if fE1 is not None:
                    assert fE1 >= 0, f"Negative GSF value for {n}"


class TestExperimentalGDR:
    """Tests for the RIPL-4 SLO / SMLO experimental GDR fits."""

    def test_load_slo_from_directory(self, ripl_path):
        db = gamma.exp.load_slo(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_load_smlo_from_directory(self, ripl_path):
        db = gamma.exp.load_smlo(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_load_mlo_alias(self, ripl_path):
        """Loading 'MLO' should redirect to SMLO data."""
        smlo = gamma.exp.load_smlo(directory=ripl_path)
        mlo  = gamma.exp.load_mlo(directory=ripl_path)
        assert len(mlo.data) == len(smlo.data)

    def test_load_slo_errors(self, ripl_path):
        """Two-line errors file should populate dEr1/dWr1 etc."""
        db = gamma.exp.load_slo(directory=ripl_path, errors=True)
        assert db is not None
        # Should contain uncertainty values
        for entry in list(db.data.values())[:5]:
            if entry.Er1 is not None:
                assert hasattr(entry, 'dEr1')

    def test_experimental_vs_systematics_overlap(self, ripl_path):
        exp_db = gamma.exp.load_slo(directory=ripl_path)
        sys_db = gamma.systematics.load_slo(directory=ripl_path)
        common = set(exp_db.data.keys()) & set(sys_db.data.keys())
        assert len(common) > 0


class TestSystematicsGDR:
    """Tests for the SLO/SMLO experimental+systematics compilation."""

    def test_load_slo(self, ripl_path):
        db = gamma.systematics.load_slo(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 1000  # ~8980 entries in the github release

    def test_load_smlo(self, ripl_path):
        db = gamma.systematics.load_smlo(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 1000

    def test_In_flag_values(self, ripl_path):
        db = gamma.systematics.load_slo(directory=ripl_path)
        seen = set()
        for entry in list(db.data.values())[:200]:
            if entry.In is not None:
                seen.add(entry.In)
        # Should see both 0 (systematics) and 1 (experimental)
        assert seen <= {0, 1}


class TestD1M:
    """Tests for the D1M+QRPA microscopic PSF tables."""

    def test_load_all(self, ripl_path):
        db = gamma.d1m.load_all(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 0

    def test_entry_shape(self, ripl_path):
        db = gamma.d1m.load_all(directory=ripl_path)
        # Check the first few entries
        for n in list(db.data.keys())[:3]:
            pkt = db.data[n]
            assert 'U' in pkt and 'fE1' in pkt
            assert len(pkt['U']) > 100
            # 10 temperature columns expected
            assert isinstance(pkt['fE1_T'][0], list)

    def test_load_element(self, ripl_path):
        db = gamma.d1m.load_element(Z=8, directory=ripl_path)
        assert len(db.data) > 0
        for n in db.data:
            assert n.Z == 8


class TestSMLO_E1:
    """Tests for the per-nucleus SMLO E1 photoabsorption tables."""

    def test_load_single_nucleus(self, ripl_path):
        # Just load one nucleus rather than the whole 8980-file directory
        db = gamma.smlo_e1.load_nucleus(Nuclide(Z=8, A=16), directory=ripl_path)
        if len(db.data) == 0:
            pytest.skip("smlo_E1 not available in this RIPL layout")
        assert Nuclide(Z=8, A=16) in db.data
        pkt = db.data[Nuclide(Z=8, A=16)]
        assert 'U' in pkt and 'fE1' in pkt
        assert len(pkt['U']) > 100


class TestSMLO_M1:
    """Tests for the per-Z SMLO M1 strength tables."""

    def test_load_element(self, ripl_path):
        db = gamma.smlo_m1.load_element(Z=26, directory=ripl_path)
        if len(db.data) == 0:
            pytest.skip("smlo_M1 not available")
        for n in db.data:
            assert n.Z == 26


class TestTLO:
    """Tests for the per-Z TLO E1 strength tables."""

    def test_load_element(self, ripl_path):
        db = gamma.tlo.load_element(Z=34, directory=ripl_path)
        if len(db.data) == 0:
            pytest.skip("TLO data not available")
        for n in db.data:
            assert n.Z == 34
            pkt = db.data[n]
            assert 'U' in pkt and 'fE1' in pkt
            assert 'beta' in pkt


class TestPSFDatabase:
    """Smoke tests for the partial PSF experimental database loader."""

    def test_load_category_arcdrc(self, ripl_path):
        db = gamma.psf.load_category('arcdrc', directory=ripl_path)
        if len(db.data) == 0:
            pytest.skip("PSF arcdrc data not available")
        assert len(db.data) > 0
        # Entries may be a list of packages (multiple datasets per nucleus)
        for n, entries in list(db.data.items())[:3]:
            assert n.Z > 0 and n.A > 0


class TestGammaSectionLoad:
    """Verify ``gamma.load(directory)`` populates the expected attributes."""

    def test_section_load_does_not_raise(self, ripl_path):
        gamma.load(directory=ripl_path)
        # Attributes that should always be present
        for attr in (
            'gsf', 'theory_gdr', 'experiment_slo', 'experiment_smlo',
            'experiment_mlo', 'experiment_systematics_slo',
            'experiment_systematics_smlo', 'gsf_d1m', 'smlo_m1', 'tlo', 'psf',
        ):
            assert hasattr(gamma.db, attr), f"db.{attr} missing"
