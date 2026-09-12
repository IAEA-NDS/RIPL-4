# -*- coding: utf-8 -*-
"""Tests for the RIPLpy API.

Tests cover:
    #1: Convenience functions (get_mass, get_mass_entry, etc.)
    #2: Database discovery/introspection (list_sections, list_databases, list_nuclei)
    #3: Simplified configuration (set_path, get_path)
    #4: Query/filter methods (filter, filter_by_range, find, count)

Note: These tests load only the specific databases needed rather than calling
riplpy.load() to ensure fast test execution and test realistic usage patterns.
"""

import os
import pytest

import riplpy
import riplpy.masses as masses
import riplpy.densities as densities
import riplpy.gamma as gamma
import riplpy.fission as fission
import riplpy.resonances as resonances
from riplpy.collections import Nuclide, Nucleus


# =============================================================================
# API: Convenience Functions (Direct Loading Tests)
# =============================================================================

class TestConvenienceFunctionsDirectLoad:
    """Tests for convenience functions using direct database loading."""

    def test_get_mass_direct_load(self, ripl_path):
        """Test getting mass by loading database directly."""
        # Load only AME20
        db = masses.ame20.load(directory=ripl_path)

        # Get entry directly
        n = Nuclide(Z=82, A=208)
        entry = db.get(n)

        assert entry is not None
        assert entry.Mexp is not None
        # Pb-208 mass excess should be around -21.7 MeV
        assert -22.0 < entry.Mexp < -21.0

    def test_get_mass_frdm12(self, ripl_path):
        """Test getting mass from FRDM12 model."""
        db = masses.frdm12.load(directory=ripl_path)

        n = Nuclide(Z=82, A=208)
        entry = db.get(n)

        assert entry is not None
        assert entry.Mth is not None  # Theoretical mass

    def test_get_mass_different_models_compare(self, ripl_path):
        """Test comparing masses from different models."""
        ame = masses.ame20.load(directory=ripl_path)
        frdm = masses.frdm12.load(directory=ripl_path)

        n = Nuclide(Z=82, A=208)
        ame_mass = ame.get(n).Mexp
        frdm_mass = frdm.get(n).Mth

        # Both should exist and be reasonably close
        assert ame_mass is not None
        assert frdm_mass is not None
        assert abs(ame_mass - frdm_mass) < 3.0  # Within 3 MeV

    def test_get_level_density_direct(self, ripl_path):
        """Test getting level density parameters."""
        from riplpy.exceptions import RiplFileNotFoundError
        try:
            db = densities.bsfg.load(directory=ripl_path)
        except (FileNotFoundError, RiplFileNotFoundError):
            pytest.skip("BSFG data file (level-densities-bfmeff.dat) not present")

        # Na-24 is in the BSFG database
        n = Nuclide(Z=11, A=24)
        entry = db.get(n)

        assert entry is not None
        assert hasattr(entry, 'ainf') or hasattr(entry, 'a')

    def test_get_gdr_direct(self, ripl_path):
        """Test getting GDR parameters."""
        db = gamma.gdr.load(directory=ripl_path)

        # The legacy theoretical GDR file is not shipped in the github layout.
        if len(db.data) == 0:
            pytest.skip("gdr-parameters-theor.dat not present in github layout")

        # Pb-208 should have GDR parameters
        n = Nuclide(Z=82, A=208)
        entry = db.get(n)

        assert entry is not None
        assert hasattr(entry, 'E1')
        assert hasattr(entry, 'W1')

    def test_get_resonance_direct(self, ripl_path):
        """Test getting resonance parameters."""
        db = resonances.swave.load(directory=ripl_path)

        # Check if database has any data
        if len(db.data) > 0:
            # Get first available entry
            n = list(db.data.keys())[0]
            entry = db.get(n)
            assert entry is not None


# =============================================================================
# API: Discovery/Introspection
# =============================================================================

class TestDiscoveryIntrospection:
    """Tests for discovery/introspection functions."""

    def test_list_sections(self):
        """Test list_sections returns all RIPL sections."""
        sections = riplpy.list_sections()

        assert isinstance(sections, list)
        assert 'masses' in sections
        assert 'densities' in sections
        assert 'fission' in sections
        assert 'gamma' in sections
        assert 'levels' in sections
        assert 'resonances' in sections

    def test_list_databases_single_section(self, ripl_path):
        """Test list_databases for a single section after loading."""
        # Load just the masses section
        masses.load(directory=ripl_path)

        mass_dbs = riplpy.list_databases('masses')

        assert isinstance(mass_dbs, list)
        assert 'ame20' in mass_dbs
        assert len(mass_dbs) > 5  # Should have multiple mass models

    def test_list_databases_invalid_section(self):
        """Test list_databases with invalid section raises ValueError."""
        with pytest.raises(ValueError, match="Unknown section"):
            riplpy.list_databases('invalid_section')

    def test_list_nuclei_after_load(self, ripl_path):
        """Test list_nuclei with explicit section after loading database."""
        # Load only AME20
        masses.ame20.load(directory=ripl_path)

        nuclei = riplpy.list_nuclei('ame20', section='masses')

        assert isinstance(nuclei, list)
        assert len(nuclei) > 1000  # AME should have many nuclei
        # All should be Nuclide objects
        assert all(isinstance(n, Nuclide) for n in nuclei)

    def test_list_nuclei_dotted_format(self, ripl_path):
        """Test list_nuclei with 'section.database' format."""
        masses.ame20.load(directory=ripl_path)

        nuclei = riplpy.list_nuclei('masses.ame20')

        assert isinstance(nuclei, list)
        assert len(nuclei) > 1000

    def test_list_nuclei_invalid_database(self, ripl_path):
        """Test list_nuclei with invalid database raises ValueError."""
        masses.ame20.load(directory=ripl_path)

        with pytest.raises(ValueError, match="not found"):
            riplpy.list_nuclei('invalid_database')

    def test_get_database(self, ripl_path):
        """Test get_database returns database object."""
        masses.ame20.load(directory=ripl_path)

        db = riplpy.get_database('ame20')

        assert db is not None
        assert hasattr(db, 'data')
        assert hasattr(db, 'get')
        assert len(db.data) > 0

    def test_get_database_dotted_format(self, ripl_path):
        """Test get_database with 'section.database' format."""
        masses.ame20.load(directory=ripl_path)

        db = riplpy.get_database('masses.ame20')

        assert db is not None
        assert len(db.data) > 0


# =============================================================================
# API: Simplified Configuration
# =============================================================================

class TestConfiguration:
    """Tests for configuration functions."""

    def test_get_path_returns_configured_path(self, ripl_path):
        """Test get_path returns the configured RIPL path."""
        # The conftest sets RIPL_LOCATION
        path = riplpy.get_path()

        assert path is not None
        assert os.path.isdir(path)

    def test_set_path_valid_directory(self, ripl_path):
        """Test set_path with valid directory."""
        # Save original
        original = riplpy.get_path()

        # Set to the same valid path (we know it exists)
        riplpy.set_path(ripl_path)

        assert riplpy.get_path() == ripl_path

        # Restore original
        if original:
            riplpy.set_path(original)

    def test_set_path_invalid_directory(self):
        """Test set_path with invalid directory raises ValueError."""
        with pytest.raises(ValueError, match="does not exist"):
            riplpy.set_path('/nonexistent/path/to/ripl')

    def test_config_module_functions(self):
        """Test config module has required functions."""
        from riplpy import config

        assert hasattr(config, 'set_path')
        assert hasattr(config, 'get_path')
        assert hasattr(config, 'clear_path')
        assert hasattr(config, 'write_config_file')


# =============================================================================
# API: Threaded "quick" loading
# =============================================================================

class TestQuickThreadedLoad:
    """The quick=True path must populate databases in the *calling* process.

    Regression: the old implementation used multiprocessing.Pool, whose worker
    processes loaded the data into their own memory and discarded it, leaving
    the caller's module-level db accessors empty. It also ignored the explicit
    directory argument. The replacement uses a ThreadPoolExecutor (threads
    share memory) and forwards the directory to each section loader.
    """

    def test_parallel_load_populates_in_process(self, ripl_path):
        """_parallel_load (the per-thread worker) populates this process's db."""
        riplpy._parallel_load('resonances', ripl_path)
        assert resonances.db.swave is not None
        assert len(resonances.db.swave.data) > 0

    def test_parallel_load_honours_directory(self, ripl_path):
        """The directory argument is forwarded to the section loader."""
        # A bogus directory must not silently populate from a stale config.
        masses.db.ame20 = None
        riplpy._parallel_load('masses', '/nonexistent/ripl/path')
        # Loader error is swallowed with a warning; db stays unpopulated.
        assert masses.db.ame20 is None
        # The real path repopulates it.
        riplpy._parallel_load('masses', ripl_path)
        assert masses.db.ame20 is not None
        assert len(masses.db.ame20.data) > 0

    def test_quick_load_matches_serial_subset(self, ripl_path):
        """Threaded executor over multiple sections populates each in-process."""
        from concurrent.futures import ThreadPoolExecutor

        subset = ('masses', 'resonances', 'fission')
        with ThreadPoolExecutor(max_workers=len(subset)) as ex:
            for fut in [ex.submit(riplpy._parallel_load, s, ripl_path)
                        for s in subset]:
                fut.result()

        assert len(masses.db.ame20.data) > 0
        assert len(resonances.db.swave.data) > 0
        assert len(fission.db.bskg3_barriers.data) > 0


# =============================================================================
# API: Database shorthand/longhand aliases
# =============================================================================

class TestDatabaseAliases:
    """Shorthand/longhand aliases resolve to the same database, work through
    get_database/list_nuclei, and do not pollute list_databases.
    """

    def test_direct_attribute_aliases(self, ripl_path):
        masses.load(directory=ripl_path)
        fission.load(directory=ripl_path)
        import riplpy.levels as levels
        levels.load(directory=ripl_path)
        assert masses.db.frdm12 is masses.db.frdm2012
        assert masses.db.frdm95 is masses.db.frdm1995
        assert masses.db.ame2020 is masses.db.ame20
        assert fission.db.bskg3 is fission.db.bskg3_barriers
        assert fission.db.empirical is fission.db.empirical_barriers
        assert levels.db.ct is levels.db.constant_temperature

    def test_get_database_accepts_aliases(self, ripl_path):
        riplpy.set_path(ripl_path)
        masses.load(directory=ripl_path)
        import riplpy.levels as levels
        levels.load(directory=ripl_path)
        assert riplpy.get_database('masses.frdm12') is masses.db.frdm2012
        assert riplpy.get_database('levels.ct') is levels.db.constant_temperature

    def test_aliases_not_listed(self, ripl_path):
        masses.load(directory=ripl_path)
        listed = riplpy.list_databases('masses')
        assert 'frdm2012' in listed       # canonical present
        assert 'frdm12' not in listed     # alias not listed

    def test_unknown_database_still_errors(self, ripl_path):
        masses.load(directory=ripl_path)
        with pytest.raises(ValueError):
            riplpy.get_database('masses.does_not_exist')


# =============================================================================
# API: Query/Filter Methods
# =============================================================================

class TestFilterMethods:
    """Tests for filter/query methods."""

    @pytest.fixture
    def ame20_db(self, ripl_path):
        """Fixture that returns a loaded AME20 database."""
        return masses.ame20.load(directory=ripl_path)

    def test_filter_by_z(self, ame20_db):
        """Test filter by atomic number Z."""
        pb_isotopes = ame20_db.filter(Z=82)

        assert len(pb_isotopes.data) > 0
        # All should be Pb (Z=82)
        for n in pb_isotopes.data.keys():
            assert n.Z == 82

    def test_filter_by_a(self, ame20_db):
        """Test filter by mass number A."""
        isobars = ame20_db.filter(A=208)

        assert len(isobars.data) > 0
        # All should have A=208
        for n in isobars.data.keys():
            assert n.A == 208

    def test_filter_by_predicate(self, ame20_db):
        """Test filter with custom predicate function."""
        # Filter for neutron-rich nuclei (N > Z)
        neutron_rich = ame20_db.filter(lambda e: e.n.N > e.n.Z)

        assert len(neutron_rich.data) > 0
        # All should have N > Z
        for n in neutron_rich.data.keys():
            assert n.N > n.Z

    def test_filter_chaining(self, ame20_db):
        """Test chaining multiple filter operations."""
        # Heavy neutron-rich nuclei
        result = ame20_db.filter(lambda e: e.n.A > 100).filter(lambda e: e.n.N > e.n.Z)

        assert len(result.data) > 0
        for n in result.data.keys():
            assert n.A > 100
            assert n.N > n.Z

    def test_filter_by_range(self, ame20_db):
        """Test filter_by_range method."""
        # Light nuclei (A <= 50)
        light = ame20_db.filter_by_range('A', max_val=50)

        assert len(light.data) > 0
        for n in light.data.keys():
            assert n.A <= 50

    def test_filter_by_range_min_max(self, ame20_db):
        """Test filter_by_range with both min and max."""
        # Medium mass nuclei (50 < A < 150)
        medium = ame20_db.filter_by_range('A', min_val=50, max_val=150)

        assert len(medium.data) > 0
        for n in medium.data.keys():
            assert 50 <= n.A <= 150

    def test_find_returns_entry(self, ame20_db):
        """Test find returns matching entry."""
        # Find Fe-56
        entry = ame20_db.find(lambda e: e.n.Z == 26 and e.n.A == 56)

        assert entry is not None
        assert entry.n.Z == 26
        assert entry.n.A == 56

    def test_find_returns_none_for_no_match(self, ame20_db):
        """Test find returns None when no match."""
        # Try to find impossible nucleus
        entry = ame20_db.find(lambda e: e.n.Z == 999)

        assert entry is None

    def test_count_with_predicate(self, ame20_db):
        """Test count with predicate."""
        heavy_count = ame20_db.count(lambda e: e.n.A > 200)

        assert heavy_count > 0
        assert isinstance(heavy_count, int)

    def test_count_with_kwargs(self, ame20_db):
        """Test count with keyword arguments."""
        pb_count = ame20_db.count(Z=82)

        assert pb_count > 0
        assert isinstance(pb_count, int)

    def test_len_support(self, ame20_db):
        """Test __len__ support."""
        length = len(ame20_db)

        assert length > 0
        assert length == len(ame20_db.data)

    def test_contains_support(self, ame20_db):
        """Test __contains__ support."""
        n = Nuclide(Z=82, A=208)

        assert n in ame20_db

        # Non-existent nucleus
        fake = Nuclide(Z=200, A=500)
        assert fake not in ame20_db


# =============================================================================
# Entry Metadata Tests
# =============================================================================

class TestEntryMetadata:
    """Tests for Entry class metadata features."""

    def test_entry_field_info_exists(self):
        """Test that Entry classes have _field_info metadata."""
        from riplpy.masses.ame20 import Entry

        assert hasattr(Entry, '_field_info')
        assert isinstance(Entry._field_info, dict)
        assert 'n' in Entry._field_info
        assert 'Mexp' in Entry._field_info

    def test_entry_field_info_method(self):
        """Test field_info() method."""
        from riplpy.masses.ame20 import Entry

        info = Entry.field_info()
        assert isinstance(info, dict)
        assert 'n' in info
        assert 'description' in info['n']

    def test_entry_summary(self, ripl_path):
        """Test Entry summary() method."""
        db = masses.ame20.load(directory=ripl_path)
        entry = db.get(Nuclide(Z=82, A=208))

        summary = entry.summary()

        assert isinstance(summary, str)
        assert '82' in summary or 'Pb' in summary

    def test_database_repr(self, ripl_path):
        """Test that Database has informative __repr__."""
        db = masses.ame20.load(directory=ripl_path)
        repr_str = repr(db)

        # Should include database info
        assert 'Database' in repr_str or 'database' in repr_str.lower()


# =============================================================================
# Integration Tests
# =============================================================================

class TestAPIIntegration:
    """Integration tests combining multiple API features."""

    def test_filter_and_iterate(self, ripl_path):
        """Test filtering and iterating over results."""
        db = masses.ame20.load(directory=ripl_path)

        # Filter Pb isotopes and check each entry
        pb = db.filter(Z=82)

        for n, entry in pb.data.items():
            assert n.Z == 82
            assert entry.n == n
            assert entry.Mexp is not None or entry.Err is not None

    def test_multiple_databases(self, ripl_path):
        """Test working with multiple databases."""
        from riplpy.exceptions import RiplFileNotFoundError
        ame = masses.ame20.load(directory=ripl_path)
        try:
            other = densities.egsm.load(directory=ripl_path)
        except (FileNotFoundError, RiplFileNotFoundError):
            pytest.skip("EGSM data file not present")

        # Find common nuclei
        ame_nuclei = set(ame.data.keys())
        other_nuclei = set(other.data.keys())
        common = ame_nuclei & other_nuclei

        assert len(common) > 0

    def test_nuclide_as_key(self, ripl_path):
        """Test that Nuclide objects work correctly as dictionary keys."""
        db = masses.ame20.load(directory=ripl_path)

        # Create equivalent Nuclide objects
        n1 = Nuclide(Z=82, A=208)
        n2 = Nuclide(Z=82, A=208)

        # Both should access the same entry
        entry1 = db.get(n1)
        entry2 = db.get(n2)

        assert entry1.Mexp == entry2.Mexp


# =============================================================================
# Models Enum Tests
# =============================================================================

class TestModelsEnum:
    """Tests for the Models enum constants."""

    def test_models_import(self):
        """Test that Models can be imported from riplpy."""
        from riplpy import Models
        assert hasattr(Models, 'Mass')
        assert hasattr(Models, 'Density')
        assert hasattr(Models, 'Fission')
        assert hasattr(Models, 'Wave')

    def test_mass_model_values(self):
        """Test MassModel enum values."""
        from riplpy import Models
        assert Models.Mass.AME20 == 'ame20'
        assert Models.Mass.FRDM12 == 'frdm12'
        assert Models.Mass.HFB27 == 'hfb27'

    def test_density_model_values(self):
        """Test DensityModel enum values."""
        from riplpy import Models
        assert Models.Density.BSFG == 'bsfg'
        assert Models.Density.CT == 'ct'

    def test_fission_model_values(self):
        """Test FissionModel enum values."""
        from riplpy import Models
        assert Models.Fission.EMPIRICAL == 'empirical'
        assert Models.Fission.HFB == 'hfb'

    def test_wave_values(self):
        """Test ResonanceWave enum values."""
        from riplpy import Models
        assert Models.Wave.S == 's'
        assert Models.Wave.P == 'p'

    def test_enum_usable_as_string(self):
        """Test that enum values work as strings in API calls."""
        from riplpy import Models
        # Enum values should equal their string values
        assert Models.Mass.AME20 == 'ame20'
        assert str(Models.Mass.AME20) == 'MassModel.AME20'


# =============================================================================
# Database Contains Method Tests
# =============================================================================

class TestDatabaseContains:
    """Tests for the contains() method on databases."""

    def test_contains_method_exists(self, ripl_path):
        """Test that contains() method exists on Database."""
        db = masses.ame20.load(directory=ripl_path)
        assert hasattr(db, 'contains')
        assert callable(db.contains)

    def test_contains_returns_true_for_existing(self, ripl_path):
        """Test contains() returns True for existing nucleus."""
        db = masses.ame20.load(directory=ripl_path)
        n = Nuclide(Z=82, A=208)
        assert db.contains(n) is True

    def test_contains_returns_false_for_missing(self, ripl_path):
        """Test contains() returns False for non-existent nucleus."""
        db = masses.ame20.load(directory=ripl_path)
        fake = Nuclide(Z=200, A=500)
        assert db.contains(fake) is False

    def test_contains_matches_in_operator(self, ripl_path):
        """Test contains() gives same result as 'in' operator."""
        db = masses.ame20.load(directory=ripl_path)
        n = Nuclide(Z=82, A=208)
        fake = Nuclide(Z=200, A=500)

        assert db.contains(n) == (n in db)
        assert db.contains(fake) == (fake in db)


# =============================================================================
# Database Nuclei Property Tests
# =============================================================================

class TestDatabaseNucleiProperty:
    """Tests for the nuclei property on NuclideDatabase."""

    def test_nuclei_property_exists(self, ripl_path):
        """Test that nuclei property exists on NuclideDatabase."""
        db = masses.ame20.load(directory=ripl_path)
        assert hasattr(db, 'nuclei')

    def test_nuclei_returns_list(self, ripl_path):
        """Test nuclei returns a list."""
        db = masses.ame20.load(directory=ripl_path)
        nuclei = db.nuclei
        assert isinstance(nuclei, list)
        assert len(nuclei) > 0

    def test_nuclei_matches_data_keys(self, ripl_path):
        """Test nuclei matches data.keys()."""
        db = masses.ame20.load(directory=ripl_path)
        assert db.nuclei == list(db.data.keys())

    def test_nuclei_contains_nuclide_objects(self, ripl_path):
        """Test nuclei contains Nuclide objects."""
        db = masses.ame20.load(directory=ripl_path)
        for n in db.nuclei[:10]:  # Check first 10
            assert hasattr(n, 'Z')
            assert hasattr(n, 'A')


# =============================================================================
# Batch Operations Tests
# =============================================================================

class TestBatchOperations:
    """Tests for batch operation functions."""

    def test_get_masses_with_nuclides(self, ripl_path):
        """Test get_masses with Nuclide objects."""
        masses.load(directory=ripl_path)

        nuclei = [Nuclide(Z=82, A=208), Nuclide(Z=82, A=206), Nuclide(Z=82, A=204)]
        result = riplpy.get_masses(nuclei, model='ame20')

        assert len(result) == 3
        assert Nuclide(Z=82, A=208) in result
        assert isinstance(result[Nuclide(Z=82, A=208)], float)

    def test_get_masses_with_tuples(self, ripl_path):
        """Test get_masses with (Z, A) tuples."""
        masses.load(directory=ripl_path)

        nuclei = [(82, 208), (82, 206), (82, 204)]
        result = riplpy.get_masses(nuclei, model='ame20')

        assert len(result) == 3
        # Results should be keyed by Nuclide objects
        assert Nuclide(Z=82, A=208) in result

    def test_get_masses_with_lists(self, ripl_path):
        """Test get_masses with [Z, A] lists."""
        masses.load(directory=ripl_path)

        nuclei = [[82, 208], [82, 206]]
        result = riplpy.get_masses(nuclei, model='ame20')

        assert len(result) == 2

    def test_get_masses_mixed_input(self, ripl_path):
        """Test get_masses with mixed input types."""
        masses.load(directory=ripl_path)

        nuclei = [Nuclide(Z=82, A=208), (82, 206), [82, 204]]
        result = riplpy.get_masses(nuclei, model='ame20')

        assert len(result) == 3

    def test_get_masses_skip_missing(self, ripl_path):
        """Test get_masses with skip_missing=True."""
        masses.load(directory=ripl_path)

        # Include a nucleus that doesn't exist in AME (Z=118 is valid element but A=500 not in database)
        nuclei = [(82, 208), (118, 500)]  # 118, 500 doesn't exist in database
        result = riplpy.get_masses(nuclei, model='ame20', skip_missing=True)

        # Should only have the valid one
        assert len(result) == 1
        assert Nuclide(Z=82, A=208) in result

    def test_get_masses_raises_for_missing(self, ripl_path):
        """Test get_masses raises error for missing nuclei by default."""
        masses.load(directory=ripl_path)
        from riplpy.exceptions import NucleusNotFoundError

        nuclei = [(82, 208), (118, 500)]  # 118, 500 doesn't exist in database

        with pytest.raises(NucleusNotFoundError):
            riplpy.get_masses(nuclei, model='ame20')

    def test_get_masses_invalid_model(self, ripl_path):
        """Test get_masses with invalid model."""
        masses.load(directory=ripl_path)

        with pytest.raises(ValueError, match="Unknown mass model"):
            riplpy.get_masses([(82, 208)], model='invalid')

    def test_get_mass_entries_batch(self, ripl_path):
        """Test get_mass_entries returns full entry objects."""
        masses.load(directory=ripl_path)

        nuclei = [(82, 208), (82, 206)]
        result = riplpy.get_mass_entries(nuclei, model='ame20')

        assert len(result) == 2
        # Check that we get full entry objects
        entry = result[Nuclide(Z=82, A=208)]
        assert hasattr(entry, 'Mexp')
        assert hasattr(entry, 'Err')

    def test_get_level_densities_batch(self, ripl_path):
        """Test get_level_densities batch operation."""
        densities.load(directory=ripl_path)

        # Prefer BSFG when available; otherwise fall back to EGSM (which is
        # shipped in the github RIPL-4 release).
        if densities.db.bsfg is not None:
            model = 'bsfg'
            nuclei = [(11, 24), (12, 26)]
        elif densities.db.egsm is not None:
            model = 'egsm'
            nuclei = [(26, 56), (50, 120)]
        else:
            pytest.skip("No level density database available in this layout")

        result = riplpy.get_level_densities(nuclei, model=model, skip_missing=True)

        # At least one should be found
        assert len(result) >= 1

    def test_get_gdrs_batch(self, ripl_path):
        """Test get_gdrs batch operation."""
        gamma.load(directory=ripl_path)

        nuclei = [(82, 208)]
        result = riplpy.get_gdrs(nuclei, skip_missing=True)

        if len(result) > 0:
            entry = list(result.values())[0]
            assert hasattr(entry, 'E1')

    def test_get_fission_barriers_batch(self, ripl_path):
        """Test get_fission_barriers batch operation."""
        fission.load(directory=ripl_path)

        nuclei = [(92, 235), (92, 238)]
        result = riplpy.get_fission_barriers(nuclei, model='empirical', skip_missing=True)

        # Results depend on what's in database
        assert isinstance(result, dict)

    def test_batch_empty_list(self, ripl_path):
        """Test batch operations with empty list."""
        masses.load(directory=ripl_path)

        result = riplpy.get_masses([], model='ame20')
        assert result == {}

    def test_batch_preserves_order(self, ripl_path):
        """Test that batch results are keyed correctly."""
        masses.load(directory=ripl_path)

        nuclei = [(26, 56), (28, 58), (82, 208)]
        result = riplpy.get_masses(nuclei, model='ame20')

        # All nuclei should be in result
        for z, a in nuclei:
            assert Nuclide(Z=z, A=a) in result


# =============================================================================
# Optical Model Potential Tests
# =============================================================================

class TestOpticalConvenienceFunctions:
    """Tests for optical model convenience functions."""

    def test_optical_in_sections(self):
        """Test that optical is in the sections tuple."""
        assert 'optical' in riplpy.sections

    def test_list_sections_includes_optical(self):
        """Test list_sections includes optical."""
        sections = riplpy.list_sections()
        assert 'optical' in sections

    def test_get_omp(self, ripl_path):
        """Test get_omp returns potential by reference number."""
        riplpy.optical.load(directory=ripl_path)

        pot = riplpy.get_omp(2405)

        assert pot is not None
        assert pot.projectile == 'n'
        assert hasattr(pot, 'header')
        assert hasattr(pot, 'components')

    def test_get_omp_invalid_raises_keyerror(self, ripl_path):
        """Test get_omp raises KeyError for invalid iref."""
        riplpy.optical.load(directory=ripl_path)

        with pytest.raises(KeyError):
            riplpy.get_omp(999999)

    def test_list_omps(self, ripl_path):
        """Test list_omps returns all potential reference numbers."""
        riplpy.optical.load(directory=ripl_path)

        irefs = riplpy.list_omps()

        assert isinstance(irefs, list)
        assert len(irefs) == 584  # Total potentials in database
        assert all(isinstance(i, int) for i in irefs)

    def test_list_omps_by_projectile(self, ripl_path):
        """Test list_omps with projectile filter."""
        riplpy.optical.load(directory=ripl_path)

        neutron_irefs = riplpy.list_omps(projectile='n')
        proton_irefs = riplpy.list_omps(projectile='p')

        assert len(neutron_irefs) > 0
        assert len(proton_irefs) > 0
        assert len(neutron_irefs) > len(proton_irefs)  # More neutron potentials

    def test_find_omp_basic(self, ripl_path):
        """Test find_omp returns potentials for a reaction."""
        riplpy.optical.load(directory=ripl_path)

        pots = riplpy.find_omp('n', 82, 208)

        assert isinstance(pots, list)
        assert len(pots) > 0
        # All should be for neutrons on Pb
        for pot in pots:
            assert pot.projectile == 'n'

    def test_find_omp_with_energy(self, ripl_path):
        """Test find_omp with energy filter."""
        riplpy.optical.load(directory=ripl_path)

        pots_all = riplpy.find_omp('n', 82, 208)
        pots_14MeV = riplpy.find_omp('n', 82, 208, E=14.0)

        # Energy filter should reduce results
        assert len(pots_14MeV) <= len(pots_all)
        assert len(pots_14MeV) > 0

    def test_get_deformation(self, ripl_path):
        """Test get_deformation returns deformation entry."""
        riplpy.optical.load(directory=ripl_path)

        # Pu-239 has a deformation entry
        deform = riplpy.get_deformation(94, 239)

        assert deform is not None
        assert hasattr(deform, 'beta')
        assert hasattr(deform, 'L')
        assert deform.Z == 94
        assert deform.A == 239

    def test_get_deformation_invalid_raises_keyerror(self, ripl_path):
        """Test get_deformation raises KeyError for missing entry."""
        riplpy.optical.load(directory=ripl_path)

        with pytest.raises(KeyError):
            riplpy.get_deformation(1, 1)  # No deformation for H-1

    def test_get_omp_reference(self, ripl_path):
        """Test get_omp_reference returns reference."""
        riplpy.optical.load(directory=ripl_path)

        ref = riplpy.get_omp_reference(100)

        assert ref is not None
        assert hasattr(ref, 'citation')
        assert hasattr(ref, 'ref_num')
        assert ref.ref_num == 100
        assert len(ref.citation) > 0

    def test_get_omp_reference_invalid_raises_keyerror(self, ripl_path):
        """Test get_omp_reference raises KeyError for invalid ref."""
        riplpy.optical.load(directory=ripl_path)

        with pytest.raises(KeyError):
            riplpy.get_omp_reference(999999)

    def test_optical_load_via_riplpy_load(self, ripl_path):
        """Test that riplpy.load() also loads optical."""
        riplpy.load(directory=ripl_path)

        # Optical should be loaded
        assert riplpy.optical.db.potentials is not None
        assert len(riplpy.optical.db.potentials.data) > 0
