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

    def test_natural_rows_do_not_overwrite_isotopes(self, ripl_path):
        """Natural-composition rows (A = 0) are kept apart from isotopes.

        They used to be stored under A = 2Z, which replaced the C-12, O-16,
        Mg-24, Si-28, S-32 and Ca-40 rows. File rows (SMLO, recommended):
            6  12 C      22.90    3.69
            6   0 C      23.14    4.24
        """
        db = gamma.exp.load_smlo(directory=ripl_path)

        c12 = db.data[Nuclide(Z=6, A=12)]
        assert (c12.Er1, c12.Wr1, c12.Id) == (pytest.approx(22.90), pytest.approx(3.69), None)

        nat_c = db.data[Nuclide(Z=6, A=0)]
        assert (nat_c.Er1, nat_c.Wr1, nat_c.Id) == (pytest.approx(23.14), pytest.approx(4.24), 'nat')

        # Every row of the recommended file is loaded: 144 isotopes and
        # 19 natural elements (readme), and no made-up A = 2Z nuclides.
        assert len(db.data) == 163
        assert sum(1 for n in db.data if n.A == 0) == 19
        assert Nuclide(Z=82, A=164) not in db.data

    def test_natural_rows_round_trip(self, ripl_path, tmp_path):
        """Writing the database back keeps A = 0 for natural rows."""
        db = gamma.exp.load_smlo(directory=ripl_path)
        out = tmp_path / 'smlo.dat'
        gamma.exp.write_smlo_ascii_file(str(out), db.data)
        again = gamma.exp.SMLO_Database()
        again.load(str(out))
        assert set(again.data) == set(db.data)

    @pytest.mark.parametrize('loader', ['load_slo', 'load_smlo'])
    def test_errors_file_keeps_every_fit(self, ripl_path, loader):
        """The errors files hold several fits per nuclide; all are kept.

        The readme gives 478 entries (one per nuclide and type of
        experimental data, the Id column). Previously only the last fit per
        nuclide was kept.
        """
        db = getattr(gamma.exp, loader)(directory=ripl_path, errors=True)
        assert len(db.data) == 163
        assert sum(len(f) for f in db.fits.values()) == 478

        # Li-6 rows in file order: Id 2, then Id 4.
        li6 = Nuclide(Z=3, A=6)
        assert [f.Id for f in db.get_fits(li6)] == ['2', '4']
        assert db.data[li6] is db.get_fits(li6)[0]

    @pytest.mark.parametrize('loader', ['load_slo', 'load_smlo'])
    def test_errors_file_first_fit_is_recommended(self, ripl_path, loader):
        """data[n] is the first fit, which the readme calls recommended.

        It must match the recommended file for every nuclide.
        """
        load = getattr(gamma.exp, loader)
        rec = load(directory=ripl_path)
        err = load(directory=ripl_path, errors=True)
        assert set(err.data) == set(rec.data)
        fields = ('Er1', 'Wr1', 'Sr1', 'CSp1', 'Er2', 'Wr2', 'Sr2', 'CSp2', 'Sr')
        for n, r in rec.data.items():
            e = err.data[n]
            for f in fields:
                assert getattr(e, f) == pytest.approx(getattr(r, f)), (n, f)

    def test_filter_keeps_fits(self, ripl_path):
        """A filtered errors database still has every fit of its nuclides."""
        db = gamma.exp.load_smlo(directory=ripl_path, errors=True)
        li = db.filter(Z=3)
        assert set(li.fits) == set(li.data)
        assert [f.Id for f in li.get_fits(Nuclide(Z=3, A=6))] == ['2', '4']

    @pytest.mark.parametrize('loader', ['load_slo', 'load_smlo'])
    def test_errors_table_round_trip(self, ripl_path, tmp_path, loader):
        """Saving an errors table writes the two-line layout with every fit.

        It used to write the one-line recommended layout, which kept only
        the recommended fit of each nuclide and dropped the uncertainties.
        """
        db = getattr(gamma.exp, loader)(directory=ripl_path, errors=True)
        out = tmp_path / 'errors.dat'
        db.save(str(out))

        again = getattr(gamma.exp, loader)(file_path=str(out))
        assert again.errors_layout
        assert set(again.fits) == set(db.fits)
        fields = ('Id', 'Er1', 'Wr1', 'Sr1', 'CSp1', 'Er2', 'Wr2', 'Sr2', 'CSp2', 'Sr',
                  'dEr1', 'dWr1', 'dSr1', 'dCSp1', 'dEr2', 'dWr2', 'dSr2', 'dCSp2', 'dSr',
                  'Emin', 'Emax', 'reference')
        for n, fits in db.fits.items():
            assert len(again.fits[n]) == len(fits), n
            for a, b in zip(fits, again.fits[n]):
                for f in fields:
                    # nan_ok: the SLO file has a literal NaN (V-51, Id 4b, dCSp2).
                    assert getattr(b, f) == pytest.approx(getattr(a, f), nan_ok=True), (n, a.Id, f)

    def test_filtered_errors_table_saves_errors_layout(self, ripl_path, tmp_path):
        db = gamma.exp.load_smlo(directory=ripl_path, errors=True).filter(Z=82)
        out = tmp_path / 'pb.dat'
        db.save(str(out))
        again = gamma.exp.load_smlo(file_path=str(out))
        assert len(again.get_fits(Nuclide(Z=82, A=208))) == 9

    def test_get_fits_missing_nuclide_raises(self, ripl_path):
        from riplpy.exceptions import NucleusNotFoundError
        db = gamma.exp.load_smlo(directory=ripl_path, errors=True)
        with pytest.raises(NucleusNotFoundError):
            db.get_fits(Nuclide(Z=26, A=56))

    def test_recommended_file_has_one_fit_each(self, ripl_path):
        """The single-line recommended file gives one fit per nuclide."""
        db = gamma.exp.load_smlo(directory=ripl_path)
        assert all(len(f) == 1 for f in db.fits.values())
        assert len(db.fits) == len(db.data)

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
