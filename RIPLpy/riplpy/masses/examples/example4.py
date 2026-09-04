# Load a ground state mass database from ASCII file, update entries, and
# save the modified database to a new file.
#
# This example uses the FRDM2012 model, which ships with both the full RIPL
# distribution and the RIPL-4 GitHub release. (For the RIPL-3 legacy FRDM1995
# model use ``masses.frdm95`` instead; that data file is only available in the
# full RIPL distribution.)

if __name__ == '__main__':
    # RIPLpy
    from riplpy.collections import Nuclide
    import riplpy.masses as masses

    # Load the database
    frdm = masses.frdm12.load()  # The FRDM2012 ground state masses database

    # Copy the database
    new_data = frdm.data
    # Update an entry with a new theoretical mass excess [MeV]
    new_data[Nuclide(z=50, a=132)].Mth = -76.112
    # Update another entry
    new_data[Nuclide(z=80, a=190)].Mth = -80.32

    # Save the modified database to a new file
    new_db = masses.frdm12.Database(new_data)
    new_db.save("example4.dat")
