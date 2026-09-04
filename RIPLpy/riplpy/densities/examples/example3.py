# -*- coding: utf-8 -*-
"""Print the location on disk where the RIPL MS shell correction database resides.

The Myers-Swiatecki shell correction database (``densities/shellcor-ms.dat``)
ships with both the full RIPL distribution and the RIPL-4 GitHub release.

Usage:
    if __name__ == '__main__':
        # RIPLpy
        import riplpy
        import riplpy.densities as dens

        # ``db_ms_file_path`` is relative to the configured RIPL directory
        print(dens.shell_corr.db_ms_file_path)
        # Resolve against the configured RIPL path for an absolute location
        import os
        print(os.path.join(riplpy.get_path(), dens.shell_corr.db_ms_file_path))

"""

# ========================

# OS
import os

# Unit testing
import unittest

# RIPLpy
import riplpy
import riplpy.densities as dens

# ========================

class TestRIPLShellCorrDatabase(unittest.TestCase):
    def test_shell_corr_db_ms_file_path(self):
        """Test the existence and validity of the RIPL MS shell correction database file path."""
        # Retrieve the (RIPL-relative) database file path
        try:
            db_ms_file_path = dens.shell_corr.db_ms_file_path
        except Exception as e:
            self.fail(f"Failed to access shell correction database file path: {e}")

        # Assert that the file path is not None or empty
        self.assertIsNotNone(db_ms_file_path, "The shell correction database file path is None")
        self.assertTrue(len(db_ms_file_path) > 0, "The shell correction database file path is empty")

        # ``db_ms_file_path`` is relative to the configured RIPL directory.
        # Resolve it against the configured RIPL path before checking on disk.
        ripl_path = riplpy.get_path()
        if not ripl_path:
            self.skipTest("RIPL path is not configured (set RIPL_LOCATION)")
        full_path = os.path.join(ripl_path, db_ms_file_path)
        self.assertTrue(
            os.path.exists(full_path),
            f"The shell correction file does not exist at: {full_path}",
        )


if __name__ == '__main__':
    unittest.main()
