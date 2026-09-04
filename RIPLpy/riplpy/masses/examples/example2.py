# Load all of the masses section databases from ASCII files into memory and
# access a mass model.
#
# This example assumes:
#     (1) The RIPL database has been downloaded to the local machine.
#     (2) The RIPL_LOCATION environment variable points to that location.
#
# FRDM2012 (masses.db.frdm2012) ships with both the full RIPL distribution and
# the RIPL-4 GitHub release. The RIPL-3 legacy FRDM1995 model is exposed as
# masses.db.frdm1995 but its data file (mass-frdm95.dat) is only available in
# the full RIPL distribution; with the GitHub release it loads as empty.

if __name__ == '__main__':
    # RIPLpy
    from riplpy.config import RIPL_PATH
    import riplpy.masses as masses

    # At this point masses.db is not defined!

    # Load all of the databases
    masses.load(RIPL_PATH)  # This may take some time!

    # Print the FRDM2012 database (RIPL-4)
    print(masses.db.frdm2012.data)

    # FRDM1995 is RIPL-3 legacy; present only in the full RIPL distribution
    if len(masses.db.frdm1995.data) > 0:
        print(masses.db.frdm1995.data)
    else:
        print("FRDM1995 (RIPL-3 legacy) not available in this RIPL release.")
