# -*- coding: utf-8 -*-
"""Python objects to parse and process the DENSITY-BSKG3 data. 

The data file contains the predictions of the deformed density distribution obtained
within the Hartree-Fock-Bogoliubov method with the BSkG3 Skyrme effective 
interaction [1 and references therein]. The axial densities are expanded in terms 
of spherical harmonics $Y_{\lambda0}$. The 3 first multipoles are provided.

The present BSkG3 compilation of neutron and proton densities includes all nuclei
with N,Z>=8 and Z<=118 between the proton and neutron driplines.
The densities are tabulated in the radius grid in steps of 0.1 fm up to 20 fm.

"""

# OS
import os as _os

# Collections
from collections import defaultdict as _defaultdict

# Dataclasses
from dataclasses import dataclass as _dataclass

# Typing
from typing import ClassVar as _ClassVar

# Logging
import logging as _logging

# Fortran format
import fortranformat as _ff

# RIPLpy
import riplpy.collections as _c
import riplpy.db as _db

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = ('DensityEntry', 'DensityDatabase', 'read_density_file', 'load_density')

# ========================

# Relative path to database folder
local_data_dir = _os.path.join("masses", "Density-bskg3")

# ========================

FORMAT_HEADER  = '(2i4,i5,f7.3)'
FORMAT_DENSITY = '(f7.3,1p,10e13.5)'
header_reader  = _ff.FortranRecordReader(FORMAT_HEADER)
density_reader = _ff.FortranRecordReader(FORMAT_DENSITY)
header_writer  = _ff.FortranRecordWriter(FORMAT_HEADER)
density_writer = _ff.FortranRecordWriter(FORMAT_DENSITY)

# ========================

@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the DENSITY-BSKG3 database.

    Contains radial density distributions from the BSkG3 HFB model,
    expanded in spherical harmonics multipoles.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        dr: Radius grid spacing [fm]
        densities: List of DensityPoint objects at each radius
    """
    n        : _c.Nuclide = None
    dr       : float = None
    densities: list = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'dr': 'Radius grid spacing [fm]',
        'densities': 'List of density multipole values at each radius',
    }

    @property
    def nrho(self, ) -> int:
        """Number of density points points. """
        return len(self.densities)


@_dataclass
class DensityPoint:
    """A density point at a specific radius.

    Contains neutron and proton density multipole expansions.

    Attributes:
        r: Radius [fm]
        rhon0: Neutron monopole density [fm^-3]
        rhon2: Neutron quadrupole density
        rhon4: Neutron hexadecapole density
        rhon6: Neutron 64-pole density
        rhon8: Neutron 256-pole density
        rhop0: Proton monopole density [fm^-3]
        rhop2: Proton quadrupole density
        rhop4: Proton hexadecapole density
        rhop6: Proton 64-pole density
        rhop8: Proton 256-pole density
    """

    r      : float = None
    rhon0  : float = None
    rhon2  : float = None
    rhon4  : float = None
    rhon6  : float = None
    rhon8  : float = None
    rhop0  : float = None
    rhop2  : float = None
    rhop4  : float = None
    rhop6  : float = None
    rhop8  : float = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'r': 'Radius [fm]',
        'rhon0': 'Neutron monopole density [fm^-3]',
        'rhon2': 'Neutron quadrupole density',
        'rhon4': 'Neutron hexadecapole density',
        'rhon6': 'Neutron 64-pole density',
        'rhon8': 'Neutron 256-pole density',
        'rhop0': 'Proton monopole density [fm^-3]',
        'rhop2': 'Proton quadrupole density',
        'rhop4': 'Proton hexadecapole density',
        'rhop6': 'Proton 64-pole density',
        'rhop8': 'Proton 256-pole density',
    }


def read_ascii_file(filepath: str) -> dict:
    """
    Parse the DENSITY-BSKG3 data file and return a dictionary of Entry objects.

    Args:
        filepath (str): Path to the input file.

    Returns:
        dict: A dictionary where keys are Nuclide objects and values are Entry objects containing density data.
    """

    density_data = _defaultdict(list)

    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as file:
            lines = iter(file.readlines())  # Create an iterator for the file lines
            for line in lines:
                if not line.strip():  # Skip empty lines
                    continue

                try:
                    # Parse the header line
                    header_data = header_reader.read(line)
                    Z, A, nrho, dr = header_data
                    n = _c.Nuclide(Z=int(Z), A=int(A))  # Create Nuclide object
                    density_entry = Entry(n=n, dr=dr, densities=[])

                    # Parse the density data lines
                    for _ in range(int(nrho)):
                        density_line = next(lines) # Get the next line
                        density_data_point = density_reader.read(density_line)
                        # Create a DensityPoint object with unpacked data
                        density_point = DensityPoint(*density_data_point)
                        density_entry.densities.append(density_point)

                    # Add the entry to the dictionary
                    density_data[n] = density_entry

                except StopIteration:
                    _logger.warning("Unexpected end of file while reading density points.")
                    break
                except Exception as e:
                    _logger.warning("Failed to parse line: %s. Error: %s", line, e)

    except FileNotFoundError:
        _logger.error("File not found: %s", filepath)
    except IOError as e:
        _logger.error("Error reading file %s: %s", filepath, e)

    return density_data


def write_ascii_file(filepath: str, data: dict) -> None:
    """Write ASCII file given a filename and data dictionary. """
    try:
        with open(filepath, 'w', encoding='utf-8') as fp:
            for nuclide, entry in data.items():
                # Write the header for each nuclide
                Z = nuclide.Z  # Charge number
                A = nuclide.A  # Mass number
                nrho = len(entry.densities)  # Number of radii
                dr = entry.dr  # Radius increment
                header_line = header_writer.write([Z, A, nrho, dr])
                fp.write(header_line + '\n')

                # Write the density points
                for density in entry.densities:
                    density_line = density_writer.write([
                        density.r,  # Radius
                        density.rhon0, density.rhon2, density.rhon4, density.rhon6, density.rhon8,  # Neutron densities
                        density.rhop0, density.rhop2, density.rhop4, density.rhop6, density.rhop8   # Proton densities
                    ])
                    fp.write(density_line + '\n')

    except IOError as e:
        _logger.error("Error writing to file %s: %s", filepath, e)

    
class Database(_db.Database):
    """A database for the density-BSKG3 data."""

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


    def load(self, fpath: str) -> None:
        """Read discrete level information from an ASCII file. 
           Not this can be used in successive calls, each time adding potentially new entries to the database (self.data).
        """
        _data = type(self).reader(fpath)
        for entry in _data.keys():
            self.data[entry] = _data[entry]

    def load_all(self, directory: str) -> None:
        """Load all density files from the directory."""
        # Loop over the element numbers
        for Z in range(1, 119):
            filepath = _os.path.join(directory, local_data_dir, f'z{Z:03d}.dat')
            self.load(filepath)


def load_all(directory: str) -> Database:
    """Load and return the Density-BSKG3 database. """
    d = Database()
    d.load_all(directory)
    return d


def load_nucleus(n: _c.Nuclide, directory: str) -> "Entry":
    d = Database()
    filepath = _os.path.join(directory, local_data_dir, f'z{n.Z:03d}.dat')
    d.load(filepath)
    return d.data[n]
