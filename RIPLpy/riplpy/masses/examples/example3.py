# Print the location on disk where the RIPL natural abundance database resides
# This example assumes: 
#     (1) The RIPL database has been downloaded to the user's local machine.
#     (2) The RIPL_LOCATION environmental variable has been set to the above location.

if __name__ == '__main__':
    # RIPLpy
    import riplpy.masses as masses

    # Print the local file path of this particular RIPL database
    print(masses.ab.LOCAL_FILE_PATH)
