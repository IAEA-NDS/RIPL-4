# -*- coding: utf-8 -*-
"""Test loading all densities databases and accessing the ESGM database.

Usage:
    if __name__ == '__main__':
        # RIPLpy
        import riplpy.densities as dens

        # At this point dens.db is not defined!

        # Define the directory where the RIPL data resides
        # Note that this directory is assumed to contain a standard structure based on RIPL sections
        directory = "<path-to-your-directory>"

        # Load all of the databases
        dens.load(directory=directory)

        # Print the ESGM database
        print(dens.db.egsm.data)

"""

# ========================

# Unit testing
import unittest

# RIPLpy
import riplpy.densities as dens

# ========================

class TestRIPLpyDensities(unittest.TestCase):
    """A simple Unit Test object. """

    def test_load_all_densities_and_access_egsm(self) -> None:
        """Test loading all densities databases and accessing the ESGM database."""
        # Load all densities databases
        try:
            dens.load()  # This loads all the databases
        except Exception as e:
            self.fail(f"Failed to load densities databases: {e}")
        
        # Ensure the dens.db attribute is defined
        self.assertIsNotNone(dens.db, "The 'dens.db' attribute is not defined after loading databases")
        
        # Ensure the ESGM database exists in dens.db and has a data attribute
        self.assertTrue(hasattr(dens.db, 'egsm'), "The 'egsm' database is not available in 'dens.db'")
        self.assertIsNotNone(dens.db.egsm.data, "'egsm.data' is None or not properly loaded")
        self.assertTrue(len(dens.db.egsm.data) > 0, "'egsm.data' is empty")


if __name__ == '__main__':
    unittest.main()
