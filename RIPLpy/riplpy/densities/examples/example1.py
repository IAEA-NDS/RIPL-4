# -*- coding: utf-8 -*-
"""Load the BSFG level density database from ASCII file and print data.

NOTE: The Back-Shifted Fermi Gas (BSFG) level density parameters are a RIPL-3
legacy product (file ``densities/level-densities-bfmeff.dat``). They ship with
the full RIPL distribution but are NOT part of the RIPL-4 GitHub release. This
example skips automatically when the legacy data file is unavailable.

Usage:
    if __name__ == '__main__':
        # RIPLpy
        import riplpy.densities as dens

        # Load the database (requires the full RIPL distribution)
        bsfg = dens.bsfg.load() # The BSFG level density database
        # Print the entire database
        print(bsfg.data)

"""

# ========================

# Unit testing
import unittest

# RIPLpy
import riplpy.densities as dens
from riplpy.exceptions import RiplFileNotFoundError

# ========================

class TestRIPLpyBSFGDatabase(unittest.TestCase):
    """A simple Unit Test object for the RIPL-3 legacy BSFG database. """

    def test_load_bsfg_database(self) -> None:
        """Test loading the BSFG level density database (legacy; skipped if absent)."""
        # Attempt to load the BSFG level density database. This is a RIPL-3
        # legacy product that is absent from the RIPL-4 GitHub release, so
        # skip the test rather than fail when the data file is missing.
        try:
            bsfg = dens.bsfg.load()  # Load the database
        except (RiplFileNotFoundError, FileNotFoundError):
            self.skipTest(
                "BSFG is RIPL-3 legacy data (level-densities-bfmeff.dat); "
                "not present in this RIPL release"
            )

        # Check that the 'data' attribute is a non-empty dictionary
        self.assertIsNotNone(bsfg.data, "BSFG database 'data' attribute is None")
        self.assertIsInstance(bsfg.data, dict, "BSFG database 'data' is not a dictionary")
        self.assertTrue(len(bsfg.data) > 0, "BSFG database 'data' attribute is empty")


if __name__ == '__main__':
    unittest.main()
