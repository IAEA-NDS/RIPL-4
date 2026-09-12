# -*- coding: utf-8 -*-
"""Tests for the RIPLpy optical module.

Tests the reading and parsing of optical model potential data from RIPL-4.
"""

import os
import pytest

import riplpy.optical as optical
from riplpy.optical.spherical import SphericalOMP, parse_spherical_omp
from riplpy.optical.coupled_channel import CoupledChannelOMP, parse_coupled_channel_omp
from riplpy.optical.reader import parse_fortran_float, OMPFileReader
from riplpy.optical.index import IndexEntry, OMPIndex
from riplpy.optical import config


class TestFortranFloatParser:
    """Tests for Fortran-style float parsing."""

    def test_standard_float(self):
        """Test parsing standard float notation."""
        assert parse_fortran_float("1.23456") == pytest.approx(1.23456)
        assert parse_fortran_float("-1.23456") == pytest.approx(-1.23456)
        assert parse_fortran_float(".5") == pytest.approx(0.5)
        assert parse_fortran_float("-.5") == pytest.approx(-0.5)

    def test_fortran_exponent_notation(self):
        """Test parsing Fortran-style exponent notation (no 'e')."""
        assert parse_fortran_float("1.00000-3") == pytest.approx(1.0e-3)
        assert parse_fortran_float("1.00000+3") == pytest.approx(1.0e+3)
        assert parse_fortran_float("-1.00000-3") == pytest.approx(-1.0e-3)
        assert parse_fortran_float("2.5-2") == pytest.approx(2.5e-2)

    def test_zero_values(self):
        """Test parsing zero and near-zero values."""
        assert parse_fortran_float("0.0") == pytest.approx(0.0)
        assert parse_fortran_float(".00000") == pytest.approx(0.0)
        assert parse_fortran_float(".00000+0") == pytest.approx(0.0)

    def test_empty_string(self):
        """Test that empty string returns 0.0."""
        assert parse_fortran_float("") == pytest.approx(0.0)
        assert parse_fortran_float("   ") == pytest.approx(0.0)


class TestOMPIndex:
    """Tests for the OMP index parser."""

    def test_load_index(self, ripl_path):
        """Test loading the OMP index file."""
        idx = optical.index.load(directory=ripl_path)
        assert idx is not None
        assert len(idx.data) > 500  # Should have ~581 entries

    def test_index_entry_count(self, ripl_path):
        """Test that we parse all index entries."""
        idx = optical.index.load(directory=ripl_path)
        # RIPL-4 has 581 entries in the index
        assert len(idx.data) >= 581

    def test_get_entry_by_iref(self, ripl_path):
        """Test retrieving a specific index entry."""
        idx = optical.index.load(directory=ripl_path)
        # Entry 2405 is a well-known Koning-Delaroche potential
        entry = idx.get(2405)
        assert entry is not None
        assert entry.iref == 2405
        assert entry.projectile == 'n'
        assert entry.is_spherical

    def test_filter_by_projectile(self, ripl_path):
        """Test filtering index by projectile type."""
        idx = optical.index.load(directory=ripl_path)

        neutron_idx = idx.filter_by_projectile('n')
        assert len(neutron_idx.data) > 200  # Most potentials are for neutrons

        proton_idx = idx.filter_by_projectile('p')
        assert len(proton_idx.data) > 50

        # All entries should have matching projectile
        for entry in neutron_idx.data.values():
            assert entry.projectile == 'n'

    def test_filter_by_target(self, ripl_path):
        """Test filtering index by target Z and A."""
        idx = optical.index.load(directory=ripl_path)

        # Filter for Pb (Z=82)
        pb_idx = idx.filter_by_target(Z=82)
        assert len(pb_idx.data) > 10

        # All entries should cover Z=82
        for entry in pb_idx.data.values():
            assert entry.Z_min <= 82 <= entry.Z_max

    def test_filter_spherical(self, ripl_path):
        """Test filtering for spherical potentials only."""
        idx = optical.index.load(directory=ripl_path)

        spherical_idx = idx.filter_spherical()
        assert len(spherical_idx.data) > 400  # ~461 spherical

        for entry in spherical_idx.data.values():
            assert entry.is_spherical

    def test_filter_coupled_channel(self, ripl_path):
        """Test filtering for coupled-channel potentials."""
        idx = optical.index.load(directory=ripl_path)

        cc_idx = idx.filter_coupled_channel()
        assert len(cc_idx.data) > 100  # ~120 CC

        for entry in cc_idx.data.values():
            assert entry.is_coupled_channel

    def test_find_for_target(self, ripl_path):
        """Test finding potentials for a specific reaction."""
        idx = optical.index.load(directory=ripl_path)

        # Find potentials for n + Pb-208 at 14 MeV
        matches = idx.find_for_target('n', Z=82, A=208, E=14.0)
        assert len(matches) > 5  # Should find several potentials

        for entry in matches:
            assert entry.projectile == 'n'
            assert entry.Z_min <= 82 <= entry.Z_max
            assert entry.A_min <= 208 <= entry.A_max
            assert entry.E_min <= 14.0 <= entry.E_max

    def test_index_summary(self, ripl_path):
        """Test the summary method."""
        idx = optical.index.load(directory=ripl_path)
        summary = idx.summary()
        assert "OMP Index" in summary
        assert "potentials" in summary


class TestSphericalOMP:
    """Tests for spherical optical model potential parsing."""

    def test_load_spherical_potential(self, ripl_path):
        """Test loading a specific spherical potential."""
        from riplpy.optical.spherical import read_spherical_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        omp = read_spherical_omp_by_iref(file_path, 2405)

        assert omp is not None
        assert isinstance(omp, SphericalOMP)
        assert omp.iref == 2405
        assert omp.projectile == 'n'

    def test_spherical_components(self, ripl_path):
        """Test that spherical potential has 6 components."""
        from riplpy.optical.spherical import read_spherical_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        omp = read_spherical_omp_by_iref(file_path, 2405)

        assert len(omp.components) == 6
        # Check component accessors
        assert omp.real_volume is not None
        assert omp.imag_volume is not None

    def test_spherical_validity_ranges(self, ripl_path):
        """Test validity range checking."""
        from riplpy.optical.spherical import read_spherical_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        omp = read_spherical_omp_by_iref(file_path, 2405)

        # Koning-Delaroche is valid for wide range
        assert omp.validity.E_min < 1.0
        assert omp.validity.E_max > 100.0
        assert omp.validity.is_valid_for(E=14.0, Z=82, A=208)

    def test_spherical_model_flags(self, ripl_path):
        """Test model flags for spherical potential."""
        from riplpy.optical.spherical import read_spherical_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        omp = read_spherical_omp_by_iref(file_path, 2405)

        assert omp.flags.imodel == 0  # Spherical
        assert omp.flags.is_spherical
        assert omp.flags.iz_proj == 0  # Neutron
        assert omp.flags.ia_proj == 1

    def test_non_spherical_raises_error(self, ripl_path):
        """Test that parsing CC potential as spherical raises error."""
        from riplpy.optical.spherical import read_spherical_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        # Entry 1 is a rigid rotor (imodel=1)
        omp = read_spherical_omp_by_iref(file_path, 1)
        assert omp is None  # Should return None for non-spherical


class TestCoupledChannelOMP:
    """Tests for coupled-channel optical model potential parsing."""

    def test_load_rigid_rotor(self, ripl_path):
        """Test loading a rigid rotor potential (imodel=1)."""
        from riplpy.optical.coupled_channel import read_coupled_channel_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        omp = read_coupled_channel_omp_by_iref(file_path, 1)

        assert omp is not None
        assert isinstance(omp, CoupledChannelOMP)
        assert omp.flags.imodel == 1
        assert omp.model_name == 'rigid_rotor'
        assert omp.n_isotopes >= 1

    def test_load_vibrational(self, ripl_path):
        """Test loading a vibrational potential (imodel=2)."""
        from riplpy.optical.coupled_channel import read_coupled_channel_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        omp = read_coupled_channel_omp_by_iref(file_path, 2)

        assert omp is not None
        assert omp.flags.imodel == 2
        assert omp.model_name == 'vibrational'
        assert omp.n_isotopes >= 1

        # Check vibrational levels
        if omp.isotopes:
            iso = omp.isotopes[0]
            assert len(iso.levels) > 0

    def test_load_soft_rotor(self, ripl_path):
        """Test loading a soft rotor potential (imodel=3)."""
        from riplpy.optical.coupled_channel import read_coupled_channel_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        omp = read_coupled_channel_omp_by_iref(file_path, 609)

        assert omp is not None
        assert omp.flags.imodel == 3
        assert omp.model_name == 'soft_rotor'

        # Soft rotor should have Hamiltonian parameters
        if omp.isotopes:
            iso = omp.isotopes[0]
            assert hasattr(iso, 'hamiltonian')

    def test_cc_isotope_data(self, ripl_path):
        """Test that CC potentials have isotope data."""
        from riplpy.optical.coupled_channel import read_coupled_channel_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        omp = read_coupled_channel_omp_by_iref(file_path, 1)

        assert omp.n_isotopes > 0
        iso = omp.isotopes[0]
        assert iso.Z > 0
        assert iso.A > 0
        assert len(iso.levels) > 0

    def test_cc_get_isotope(self, ripl_path):
        """Test getting isotope data by Z and A."""
        from riplpy.optical.coupled_channel import read_coupled_channel_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        omp = read_coupled_channel_omp_by_iref(file_path, 1)

        # Get the first isotope's Z and A
        iso = omp.isotopes[0]
        found = omp.get_isotope(iso.Z, iso.A)
        assert found is not None
        assert found.Z == iso.Z

    def test_spherical_as_cc_raises_error(self, ripl_path):
        """Test that parsing spherical as CC raises error."""
        from riplpy.optical.coupled_channel import read_coupled_channel_omp_by_iref

        file_path = os.path.join(ripl_path, config.get_data_file_path('parameters'))
        # Entry 2405 is spherical (imodel=0)
        omp = read_coupled_channel_omp_by_iref(file_path, 2405)
        assert omp is None  # Should return None for spherical


class TestOMPDatabase:
    """Tests for the main OMP database."""

    def test_load_database(self, ripl_path):
        """Test loading the full OMP database."""
        db = optical.omp.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 500  # Should have ~584 potentials

    def test_database_counts(self, ripl_path):
        """Test database contains expected counts."""
        db = optical.omp.load(directory=ripl_path)

        assert db.n_spherical > 400  # ~461 spherical
        assert db.n_coupled_channel > 100  # ~123 CC
        assert db.n_spherical + db.n_coupled_channel == len(db.data)

    def test_get_by_iref(self, ripl_path):
        """Test retrieving potential by reference number."""
        db = optical.omp.load(directory=ripl_path)

        # Get spherical
        omp = db.get(2405)
        assert omp is not None
        assert omp.iref == 2405
        assert isinstance(omp, SphericalOMP)

        # Get CC
        omp = db.get(1)
        assert omp is not None
        assert omp.iref == 1
        assert isinstance(omp, CoupledChannelOMP)

    def test_get_nonexistent_raises_error(self, ripl_path):
        """Test that getting nonexistent iref raises error."""
        db = optical.omp.load(directory=ripl_path)

        with pytest.raises(KeyError):
            db.get(999999)

    def test_get_by_iref_returns_none(self, ripl_path):
        """Test get_by_iref returns None for missing entries."""
        db = optical.omp.load(directory=ripl_path)

        omp = db.get_by_iref(999999)
        assert omp is None

    def test_filter_by_projectile(self, ripl_path):
        """Test filtering database by projectile."""
        db = optical.omp.load(directory=ripl_path)

        neutron_db = db.filter_by_projectile('n')
        assert len(neutron_db.data) > 200

        for omp in neutron_db.data.values():
            assert omp.projectile == 'n'

    def test_filter_by_target(self, ripl_path):
        """Test filtering database by target."""
        db = optical.omp.load(directory=ripl_path)

        # Filter for Pb (Z=82)
        pb_db = db.filter_by_target(Z=82)
        assert len(pb_db.data) > 10

        for omp in pb_db.data.values():
            assert omp.validity.Z_min <= 82 <= omp.validity.Z_max

    def test_filter_by_energy(self, ripl_path):
        """Test filtering database by energy."""
        db = optical.omp.load(directory=ripl_path)

        # Filter for 14 MeV
        filtered = db.filter_by_energy(14.0)
        assert len(filtered.data) > 100

        for omp in filtered.data.values():
            assert omp.validity.E_min <= 14.0 <= omp.validity.E_max

    def test_filter_spherical(self, ripl_path):
        """Test filtering for spherical potentials."""
        db = optical.omp.load(directory=ripl_path)

        spherical_db = db.filter_spherical()
        assert len(spherical_db.data) == db.n_spherical

        for omp in spherical_db.data.values():
            assert isinstance(omp, SphericalOMP)

    def test_filter_coupled_channel(self, ripl_path):
        """Test filtering for coupled-channel potentials."""
        db = optical.omp.load(directory=ripl_path)

        cc_db = db.filter_coupled_channel()
        assert len(cc_db.data) == db.n_coupled_channel

        for omp in cc_db.data.values():
            assert isinstance(omp, CoupledChannelOMP)

    def test_filter_by_model(self, ripl_path):
        """Test filtering by model type."""
        db = optical.omp.load(directory=ripl_path)

        # Filter for rigid rotor (imodel=1)
        rr_db = db.filter_by_model(1)
        assert len(rr_db.data) > 90  # ~99 rigid rotor

        for omp in rr_db.data.values():
            assert omp.flags.imodel == 1

    def test_find_for_reaction(self, ripl_path):
        """Test finding potentials for a specific reaction."""
        db = optical.omp.load(directory=ripl_path)

        # Find potentials for n + Pb-208 at 14 MeV
        matches = db.find_for_reaction('n', Z=82, A=208, E=14.0)
        assert len(matches) > 5

        for omp in matches:
            assert omp.projectile == 'n'
            assert omp.validity.Z_min <= 82 <= omp.validity.Z_max
            assert omp.validity.A_min <= 208 <= omp.validity.A_max
            assert omp.validity.E_min <= 14.0 <= omp.validity.E_max

    def test_find_for_actinide_reaction(self, ripl_path):
        """Test finding potentials for actinide reaction."""
        db = optical.omp.load(directory=ripl_path)

        # Find potentials for n + U-238 at 5 MeV (includes CC potentials)
        matches = db.find_for_reaction('n', Z=92, A=238, E=5.0)
        assert len(matches) > 10

        # Should include both spherical and CC potentials
        has_spherical = any(isinstance(m, SphericalOMP) for m in matches)
        has_cc = any(isinstance(m, CoupledChannelOMP) for m in matches)
        assert has_spherical
        assert has_cc

    def test_top_level_find_for_reaction(self, ripl_path):
        """Top-level optical.find_for_reaction selects OMPs by incident particle."""
        optical.load(directory=ripl_path)

        # Alpha + Fe-56 at 20 MeV via the convenience wrapper.
        matches = optical.find_for_reaction('alpha', Z=26, A=56, E=20.0)
        assert len(matches) > 0
        for omp in matches:
            assert omp.projectile == 'a'

        # The 'a' alias selects the same set.
        assert len(optical.find_for_reaction('a', Z=26, A=56, E=20.0)) == len(matches)

    def test_database_info(self, ripl_path):
        """Test the info method."""
        db = optical.omp.load(directory=ripl_path)
        info = db.info()

        assert "OMP Database" in info
        assert "Spherical" in info
        assert "Coupled-channel" in info
        assert "By projectile" in info

    def test_database_iteration(self, ripl_path):
        """Test iterating over database."""
        db = optical.omp.load(directory=ripl_path)

        irefs = list(db)
        assert len(irefs) == len(db.data)

        for iref in irefs[:10]:
            assert iref in db.data


class TestOpticalModuleAPI:
    """Tests for the high-level optical module API."""

    def test_optical_load(self, ripl_path):
        """Test the optical.load() function."""
        optical.load(directory=ripl_path)

        assert optical.db.index is not None
        assert optical.db.potentials is not None
        assert optical.db.deformations is not None
        assert optical.db.references is not None

    def test_optical_db_index(self, ripl_path):
        """Test accessing index via optical.db."""
        optical.load(directory=ripl_path)

        assert len(optical.db.index.data) > 500
        entry = optical.db.index.get(2405)
        assert entry.iref == 2405

    def test_optical_db_potentials(self, ripl_path):
        """Test accessing potentials via optical.db."""
        optical.load(directory=ripl_path)

        assert len(optical.db.potentials.data) > 500
        omp = optical.db.potentials.get(2405)
        assert omp.iref == 2405

    def test_optical_db_deformations(self, ripl_path):
        """Test accessing deformations via optical.db."""
        optical.load(directory=ripl_path)

        assert len(optical.db.deformations.data) > 1000

    def test_optical_db_references(self, ripl_path):
        """Test accessing references via optical.db."""
        optical.load(directory=ripl_path)

        assert len(optical.db.references.data) > 100
        ref = optical.db.references.get(14)
        assert "Becchetti" in ref.citation


class TestDeformations:
    """Tests for the deformations database."""

    def test_load_deformations(self, ripl_path):
        """Test loading the deformations database."""
        db = optical.deformations.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 1000  # Should have ~1400+ entries

    def test_deformation_entry_fields(self, ripl_path):
        """Test that deformation entries have required fields."""
        db = optical.deformations.load(directory=ripl_path)

        # Get any entry
        entry = list(db.data.values())[0]
        assert hasattr(entry, 'Z')
        assert hasattr(entry, 'A')
        assert hasattr(entry, 'Ex')
        assert hasattr(entry, 'L')
        assert hasattr(entry, 'beta')

    def test_get_for_nucleus(self, ripl_path):
        """Test getting all deformations for a nucleus."""
        db = optical.deformations.load(directory=ripl_path)

        # Pb-208 is a doubly-magic nucleus, may have vibrational states
        entries = db.get_for_nucleus(82, 208)
        # Might have entries or might not - just check it returns a list
        assert isinstance(entries, list)

    def test_filter_quadrupole(self, ripl_path):
        """Test filtering for quadrupole (L=2) deformations."""
        db = optical.deformations.load(directory=ripl_path)

        quad_db = db.filter_quadrupole()
        assert len(quad_db.data) > 500  # Most are quadrupole

        for entry in quad_db.data.values():
            assert entry.L == 2

    def test_filter_octupole(self, ripl_path):
        """Test filtering for octupole (L=3) deformations."""
        db = optical.deformations.load(directory=ripl_path)

        oct_db = db.filter_octupole()
        assert len(oct_db.data) > 100  # Some are octupole

        for entry in oct_db.data.values():
            assert entry.L == 3

    def test_deformation_info(self, ripl_path):
        """Test the info method."""
        db = optical.deformations.load(directory=ripl_path)
        info = db.info()

        assert "Deformation Database" in info
        assert "Quadrupole" in info
        assert "Octupole" in info


class TestReferences:
    """Tests for the references database."""

    def test_load_references(self, ripl_path):
        """Test loading the references database."""
        db = optical.references.load(directory=ripl_path)
        assert db is not None
        assert len(db.data) > 100  # Should have ~143 references

    def test_get_reference(self, ripl_path):
        """Test getting a specific reference."""
        db = optical.references.load(directory=ripl_path)

        # Reference 14 is Becchetti-Greenlees
        ref = db.get(14)
        assert ref is not None
        assert ref.ref_num == 14
        assert "Becchetti" in ref.citation

    def test_reference_properties(self, ripl_path):
        """Test reference property extraction."""
        db = optical.references.load(directory=ripl_path)

        ref = db.get(14)
        assert ref.first_author  # Should extract author
        assert ref.year == 1969  # Known year

    def test_search_references(self, ripl_path):
        """Test searching references by text."""
        db = optical.references.load(directory=ripl_path)

        # Search for Koning
        results = db.search("Koning")
        assert len(results.data) > 0

        for ref in results.data.values():
            assert "koning" in ref.citation.lower()

    def test_search_by_author(self, ripl_path):
        """Test searching by author name."""
        db = optical.references.load(directory=ripl_path)

        results = db.search_by_author("Young")
        assert len(results.data) > 0

    def test_reference_info(self, ripl_path):
        """Test the info method."""
        db = optical.references.load(directory=ripl_path)
        info = db.info()

        assert "Reference Database" in info
        assert "references" in info.lower()


class TestDataIntegrity:
    """Tests for data integrity and consistency."""

    def test_all_index_entries_have_potentials(self, ripl_path):
        """Test that index entries correspond to actual potentials."""
        optical.load(directory=ripl_path)

        # Note: Some index entries might not have potentials loaded
        # if parsing fails, so we just check overlap
        index_irefs = set(optical.db.index.data.keys())
        pot_irefs = set(optical.db.potentials.data.keys())

        # Most potentials should be in both
        overlap = index_irefs & pot_irefs
        assert len(overlap) > 500

    def test_potential_projectile_consistency(self, ripl_path):
        """Test that potential projectile matches expected values."""
        db = optical.omp.load(directory=ripl_path)

        valid_projectiles = {'n', 'p', 'd', 't', 'He3', 'a'}
        for omp in db.data.values():
            assert omp.projectile in valid_projectiles

    def test_potential_energy_ranges_valid(self, ripl_path):
        """Test that energy ranges are valid (E_min < E_max)."""
        db = optical.omp.load(directory=ripl_path)

        for omp in db.data.values():
            assert omp.validity.E_min <= omp.validity.E_max
            assert omp.validity.E_min >= 0  # Energies should be non-negative

    def test_potential_za_ranges_valid(self, ripl_path):
        """Test that Z and A ranges are valid."""
        db = optical.omp.load(directory=ripl_path)

        for omp in db.data.values():
            assert omp.validity.Z_min <= omp.validity.Z_max
            assert omp.validity.A_min <= omp.validity.A_max
            assert omp.validity.Z_min >= 0
            assert omp.validity.A_min >= 0

    def test_spherical_has_six_components(self, ripl_path):
        """Test that all spherical potentials have 6 components."""
        db = optical.omp.load(directory=ripl_path)
        spherical_db = db.filter_spherical()

        for omp in spherical_db.data.values():
            assert len(omp.components) == 6

    def test_cc_has_isotopes(self, ripl_path):
        """Test that all CC potentials have at least one isotope."""
        db = optical.omp.load(directory=ripl_path)
        cc_db = db.filter_coupled_channel()

        for omp in cc_db.data.values():
            assert omp.n_isotopes >= 1


class TestModifiedPotentials:
    """Tests for the modified potentials loader."""

    def test_load_modified_potentials(self, ripl_path):
        """Test loading modified potentials from mod-potentials directory."""
        mod_db = optical.omp.load_modified_potentials(directory=ripl_path)

        assert mod_db is not None
        assert len(mod_db.data) == 7  # 7 modified potentials
        assert sorted(mod_db.irefs) == [1480, 1481, 1482, 2412, 2415, 4609, 4610]

    def test_modified_potential_parse_correctly(self, ripl_path):
        """Test that modified potentials parse with correct fields."""
        mod_db = optical.omp.load_modified_potentials(directory=ripl_path)

        # OMP 1480 is a rigid rotor
        pot = mod_db.get(1480)
        assert pot.iref == 1480
        assert pot.projectile == 'n'
        assert "Capote" in pot.header.author
        assert pot.flags.imodel == 1  # Rigid rotor

    def test_load_with_modifications_applies_changes(self, ripl_path):
        """Test that load_with_modifications replaces potentials."""
        # Load without modifications
        db_original = optical.omp.load(directory=ripl_path)
        original_author = db_original.get(1480).header.author

        # Load with modifications
        db_modified = optical.omp.load_with_modifications(directory=ripl_path)
        modified_author = db_modified.get(1480).header.author

        # Total count should be the same
        assert len(db_modified.data) == len(db_original.data)

        # Author should be updated in modified version
        assert "Capote" in modified_author

    def test_load_with_modifications_preserves_other_potentials(self, ripl_path):
        """Test that non-modified potentials remain unchanged."""
        db = optical.omp.load_with_modifications(directory=ripl_path)

        # OMP 2405 is not in modified potentials
        pot = db.get(2405)
        assert pot is not None
        assert pot.iref == 2405

    def test_load_with_modifications_flag_false(self, ripl_path):
        """Test that apply_modifications=False skips modifications."""
        db = optical.omp.load_with_modifications(directory=ripl_path, apply_modifications=False)

        # Should have same count as original
        db_original = optical.omp.load(directory=ripl_path)
        assert len(db.data) == len(db_original.data)

    def test_module_level_functions(self, ripl_path):
        """Test the module-level convenience functions."""
        # Test load_modified_potentials at module level
        mod_db = optical.load_modified_potentials(directory=ripl_path)
        assert len(mod_db.data) == 7

        # Test load_with_modifications at module level
        db = optical.load_with_modifications(directory=ripl_path)
        assert len(db.data) > 500


class TestROP2013:
    """Tests for the Avrigeanu revised alpha ROP table (ROP2013za.dat)."""

    def test_load_rop2013(self, ripl_path):
        """ROP2013 loads and contains the expected number of nuclei."""
        db = optical.rop2013.load(directory=ripl_path)
        assert db is not None
        assert len(db) >= 70  # 79 nuclei in RIPL-4
        assert db.iref == 9999
        assert "Avrigeanu" in db.authors

    def test_rop2013_sc45_entry(self, ripl_path):
        """The Sc-45 entry has the expected tabulated values."""
        db = optical.rop2013.load(directory=ripl_path)
        entry = db.get(21, 45)
        assert entry is not None
        assert entry.Z == 21
        assert entry.A == 45
        assert entry.n_points > 10
        # First row at E = 2 MeV
        assert entry.E[0] == pytest.approx(2.0)
        assert entry.V_R[0] == pytest.approx(164.0)
        assert entry.r_R[0] == pytest.approx(1.204)
        # sigma_R should be positive after ~5 MeV
        assert max(entry.sigma_R) > 100.0

    def test_rop2013_missing_returns_none(self, ripl_path):
        """Unknown (Z, A) returns None."""
        db = optical.rop2013.load(directory=ripl_path)
        assert db.get(999, 999) is None

    def test_rop2013_loaded_via_optical_load(self, ripl_path):
        """optical.load() attaches the ROP2013 database to optical.db."""
        optical.load(directory=ripl_path)
        assert optical.db.rop2013 is not None
        assert len(optical.db.rop2013) > 0


class TestAtomki:
    """Tests for the ATOMKI alpha-OMP gnu lazy reader."""

    def test_atomki_load_nucleus(self, ripl_path):
        """Loading a single nucleus returns a populated profile."""
        atomki_dir = os.path.join(ripl_path, 'optical', 'atomki')
        entry = optical.atomki.load_nucleus(Z=26, A=42, directory=atomki_dir)
        assert entry.Z == 26
        assert entry.A == 42
        assert entry.AP == 4
        assert entry.ZP == 2
        assert entry.JR0 > 0
        assert entry.n_points > 100
        # Mesh values are monotonic and positive
        assert entry.r[0] == pytest.approx(0.025)
        assert entry.r[-1] > entry.r[0]

    def test_atomki_missing_raises(self, ripl_path):
        """A missing (Z, A) raises FileNotFoundError."""
        atomki_dir = os.path.join(ripl_path, 'optical', 'atomki')
        with pytest.raises(FileNotFoundError):
            optical.atomki.load_nucleus(Z=999, A=999, directory=atomki_dir)

    def test_atomki_available_nuclei(self, ripl_path):
        """available_nuclei() returns all 4359 files in the github release."""
        atomki_dir = os.path.join(ripl_path, 'optical', 'atomki')
        nuclei = optical.atomki.available_nuclei(directory=atomki_dir)
        assert len(nuclei) > 4000
        # Every entry is an int tuple
        for z, a in nuclei[:10]:
            assert isinstance(z, int)
            assert isinstance(a, int)
            assert z > 0 and a > 0

    def test_atomki_filename_helper(self):
        """filename_for() produces the canonical 3-digit zero-padded name."""
        assert optical.atomki.filename_for(26, 42) == 'z026a042a_talys_alphaomp9real.gnu'
        assert optical.atomki.filename_for(8, 16) == 'z008a016a_talys_alphaomp9real.gnu'
