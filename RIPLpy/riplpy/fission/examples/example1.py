# -*- coding: utf-8 -*-
"""Load the BSkG3 fission barrier database from ASCII file and verify its data.

The BSkG3 fission barriers ship with both the full RIPL distribution and the
RIPL-4 GitHub release.

Usage:
    if __name__ == '__main__':
        # RIPLpy
        import riplpy.fission as fission

        # Load the database
        bskg3 = fission.bskg3.load()  # The BSkG3 fission barrier database
        # Print the entire database
        print(bskg3.data)

"""

# ========================

# Unit testing
import unittest

# RIPLpy
import riplpy.fission as fission

# ========================

class TestRIPLpyBSkG3FissionDatabase(unittest.TestCase):
    """A simple Unit Test object. """

    def test_load_bskg3_database(self) -> None:
        """Test loading the BSkG3 fission barrier database."""
        # Attempt to load the BSkG3 fission barrier database
        try:
            bskg3 = fission.bskg3.load()  # Load the database
        except Exception as e:
            self.fail(f"Failed to load BSkG3 fission barrier database: {e}")

        # Check that the 'data' attribute is a non-empty dictionary
        self.assertIsNotNone(bskg3.data, "BSkG3 database 'data' attribute is None")
        self.assertIsInstance(bskg3.data, dict, "BSkG3 database 'data' is not a dictionary")
        self.assertTrue(len(bskg3.data) > 0, "BSkG3 database 'data' attribute is empty")

        # Each entry should expose inner/outer barrier sections
        n, entry = next(iter(bskg3.data.items()))
        self.assertTrue(hasattr(entry, 'inner'), "BSkG3 entry is missing 'inner' barrier")
        self.assertTrue(hasattr(entry, 'outer1'), "BSkG3 entry is missing 'outer1' barrier")


if __name__ == '__main__':
    unittest.main()
