# Load the HFB-14 ground state masses database from ASCII file and print data.
#
# NOTE: HFB-14 is a RIPL-3 legacy mass model (file masses/mass-hfb14.dat). It
# ships with the full RIPL distribution but is NOT part of the RIPL-4 GitHub
# release. Use a RIPL-4 model (e.g. masses.hfb27, masses.bskg3, masses.d1m)
# if you only have the GitHub release.
#
# This example assumes:
#     (1) The full RIPL distribution has been downloaded to the local machine.
#     (2) The RIPL_LOCATION environment variable points to that location.

if __name__ == '__main__':
    # RIPLpy
    import riplpy.masses as masses
    from riplpy.exceptions import RiplFileNotFoundError

    # Load the database (RIPL-3 legacy; only in the full RIPL distribution)
    try:
        hfb14 = masses.hfb14.load()  # The HFB-14 ground state masses database
    except (RiplFileNotFoundError, FileNotFoundError):
        print(
            "HFB-14 is RIPL-3 legacy data (mass-hfb14.dat) and is not present "
            "in this RIPL release. Use masses.hfb27 / bskg3 / d1m instead."
        )
    else:
        # Print the entire database
        print(hfb14.data)
