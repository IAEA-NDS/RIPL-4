
# ========================

# OS
import os as _os

# Collections
from collections import defaultdict as _defaultdict
from dataclasses import dataclass as _dataclass
from dataclasses import field as _field

# Typing
from typing import ClassVar as _ClassVar

# Fortran format
import fortranformat as _ff

# RIPLpy
import riplpy.collections as _c
import riplpy.db as _db
from . import core as _core

# ========================

__all__ = ('local_data_dir', 'IDENT_FORTRAN_FORMAT', 'LEVEL_FORTRAN_FORMAT', 'GAMMA_FORTRAN_FORMAT', 'ident_reader', 'level_reader', 'gamma_reader',
           'read_ascii_file', 'parse_ascii_dict', 'GammaEntry', 'LevelEntry', 'Database', 'load_all', 'load_element')

# ========================

# Relative path to database folder
local_data_dir = "levels"

# ========================

IDENT_FORTRAN_FORMAT = '(a5,6i5,2f12.6)'
LEVEL_FORTRAN_FORMAT = '(i3,1x,f10.6,1x,f5.1,i3,1x,(1pe10.2),i3,1x,a1,1x,a4,1x,a18,i3,10(1x,a2,1x,0pf10.4,1x,a7),f10.6,1x,3(i2))'
GAMMA_FORTRAN_FORMAT = '(39x,i4,1x,f10.3,3(1x,e10.3))'

ident_reader = _ff.FortranRecordReader(IDENT_FORTRAN_FORMAT)
level_reader = _ff.FortranRecordReader(LEVEL_FORTRAN_FORMAT)
gamma_reader = _ff.FortranRecordReader(GAMMA_FORTRAN_FORMAT)

ident_writer = _ff.FortranRecordWriter(IDENT_FORTRAN_FORMAT)
level_writer = _ff.FortranRecordWriter(LEVEL_FORTRAN_FORMAT)
gamma_writer = _ff.FortranRecordWriter(GAMMA_FORTRAN_FORMAT)

# ========================

def read_identification_record(substr: str) -> list:
    """Read the isotope's identification record. """
    return ident_reader.read(substr)


def read_level_record(substr: str) -> list:
    """Read a level record of the particular isotope. """
    return level_reader.read(substr)


def read_ascii_file(fpath: str) -> dict:
    """Read data into the database object from a single ASCII file. """
    # Placeholder for data
    d = {}

    # RIPL discrete level files have no headers
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:    
        # Begin reading data directly from the first line
        for line in fp.readlines():
            # Parse the given line
            try:
                header     = ident_reader.read(line)
                nucleus    = header[0]
                num_levels = int(header[3])
                num1       = header[5]
                num2       = header[6]
                num3       = header[7]
                num4       = header[8]
                extras     = [num1, num2, num3, num4]
                levels     = []
                gammas     = _defaultdict(list)
            except ValueError:
                if "                                      " in line[0:40]:
                    # Find appropriate level to attach to
                    lvl = int(level[0])
                    # Parse data
                    gamma = gamma_reader.read(line)
                    #print(nucleus, lvl, gamma)
                    gammas[lvl].append(gamma)
                else:
                    # Parse data
                    level = level_reader.read(line)
                    levels.append(level)
            # Check if we save data to the dictionary
            if len(levels) == num_levels:
                d[nucleus] = (header, levels, gammas, extras)

        return d


def parse_ascii_dict(data: dict) -> dict:
    """Parse the dictionary to a more Pythonic friendly representation. """
        
    offset = 3 # Offset used for decay mode information
        
    output = {} # Placeholder for return information

    # Loop over the dictionary keys
    for nucleus in data.keys():
        n = _c.Nuclide.from_string(nucleus)
        header, levels, gammas, extras = data[nucleus]
        # Recast levels into objects
        recasted_levels = []
        for lvl in levels:
            e       = float(lvl[1])
            spin    = float(lvl[2])
            parity  = int(lvl[3])
            t12     = float(lvl[4])
            ng      = int(lvl[5])
            jflag   = lvl[6].strip()
            unc     = lvl[7].strip()
            jpistr  = lvl[8].strip()
            n_modes = int(lvl[9])
            modes = []
            for m in range(n_modes):
                dpm     = lvl[10+m*offset].strip()
                percent = float(lvl[11+m*offset])
                mode    = lvl[12+m*offset].strip()
                modes.append(_core.DecayMode(mode, dpm, percent))
            recasted_levels.append(_core.DiscreteLevelInfo(e, spin, parity, t12, ng, jflag, unc, jpistr, modes))

        # Recast the gammas into objects
        recasted_gammas = _defaultdict(list)
        for lvl_id in gammas.keys():
            for gamma in gammas[lvl_id]:
                recasted_gammas[lvl_id].append(_core.DiscreteGamma(*gamma))

        # Save for later
        output[n] = LevelEntry(n, recasted_levels, recasted_gammas, extras)

    return output


def write_ascii_file(fpath: str, data: dict) -> None:
    """"Write the discrete level information to file. This method assumes parse_ascii_dict has been called on the dictionary. """

    # Note: there is no file header
    with open(fpath, 'w', encoding='utf-8') as fp:

        for nucleus, entry in data.items():
            # Extract header information from LevelEntry
            num_levels = entry.num_levels
            num_gammas = entry.total_number_gammas  # Total gammas in all levels

            # Construct the header record for this nucleus
            header_record = [
                str(nucleus.A)+nucleus.symbol, # Nucleus identifier
                nucleus.A,
                nucleus.Z,
                num_levels,                       # Number of levels
                num_gammas,                       # Number of gamma-rays
                entry.extras[0], entry.extras[1], # Nmax, Nc (completeness limits)
                entry.extras[2], entry.extras[3]  # Sn, Sp separation energies [MeV]
            ]
            fp.write(ident_writer.write(header_record) + '\n')

            # Write each nuclear level
            for i,level in enumerate(entry.levels):
                glist = list(entry.gammas.keys())
                level_id = i + 1
                level_data = [
                    level_id,                 # Level index
                    level.energy,             # Energy of the level
                    level.spin,               # Spin value
                    level.parity,             # Parity
                    level.halflife,           # Half-life
                    level.num_gammas,         # Number of gamma transitions
                    level.jflag,              # J-flag
                    level.unc,                # Uncertainty flag
                    level.jpistr,             # Spin-parity string
                    level.number_decay_modes, # Number of decay modes
                ]

                # Append decay modes
                for dm in level.decay_modes:
                    level_data.extend([
                        dm.mode,    # Decay parent mode
                        dm.percent, # Decay probability percentage
                        dm.flag     # Decay mode description
                    ])

                # Write level record to file
                fp.write(level_writer.write(level_data) + '\n')

                if len(entry.gammas) > 0 and level_id in glist:
                    for gamma in entry.gammas[level_id]:
                        gamma_data = [
                            gamma.to_level_number,  # Gamma transition index
                            gamma.energy,          # Gamma energy
                            gamma.prob_g,          # Photon probability
                            gamma.prob_em,         # Probability of the electromagnetic transition (photon, conversion electron, pair creation)
                            gamma.cc               # Internal conversion coefficient of the transition
                        ]
                    
                        fp.write(gamma_writer.write(gamma_data) + '\n')

@_dataclass
class GammaEntry(_core.DiscreteNuclearLevel, _db.NuclideDatabaseEntry):
    """A discrete nuclear level with associated gamma-ray transitions.

    Contains level properties and all gamma-ray de-excitation data.

    Attributes:
        halflife: Level half-life [s]
        jflag: Spin assignment method flag
        unc: Uncertain energy flag
        jpistr: Original ENSDF spin-parity string
        gammas: Dict of gamma transitions keyed by level index
        decay_modes: List of decay modes
        shift: Uncertain energy offset (X value)
        band: Rotational band assignment
    """

    halflife   : float = None
    jflag      : str   = None
    unc        : str   = None
    jpistr     : str   = None
    gammas     : dict  = _field(default_factory=dict)
    decay_modes: list  = _field(default_factory=list)
    shift      : str   = None
    band       : int   = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'halflife': 'Level half-life [s]',
        'jflag': 'Spin assignment method flag',
        'unc': 'Uncertain energy flag',
        'jpistr': 'Original ENSDF spin-parity string',
        'gammas': 'Dict of gamma transitions',
        'decay_modes': 'List of decay modes',
        'shift': 'Uncertain energy offset (X value)',
        'band': 'Rotational band assignment',
    }

    @property
    def number_decay_modes(self, ) -> int:
        """Return the number of decay modes. Returns None for unknown. """
        return len(self.decay_modes)
    
    @property
    def number_gammas(self, ) -> int:
        """Return the number of gamma-rays associated with the de-excitation of this nuclear level. """
        return len(self.gammas)


@_dataclass
class LevelEntry(_db.NuclideDatabaseEntry):
    """An entry in the discrete levels database.

    Contains the complete discrete level scheme for a nucleus including
    all levels, gamma transitions, and metadata.

    Attributes:
        nucleus: Target nucleus (Nuclide object with Z, A, N properties)
        levels: List of DiscreteLevelInfo objects (ordered by energy)
        gammas: Dict of gamma transitions keyed by level index
        extras: Additional metadata from ENSDF header
    """

    nucleus: _c.Nuclide = None
    levels: list = _field(default_factory=list)
    gammas: dict = _field(default_factory=dict)
    extras: list = _field(default_factory=list)

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'nucleus': 'Target nucleus',
        'levels': 'List of discrete levels (ordered by energy)',
        'gammas': 'Dict of gamma transitions by level index',
        'extras': 'Additional ENSDF metadata',
    }

    @property
    def num_levels(self, ) -> int:
        return len(self.levels)

    @property
    def level_energies(self, ) -> list:
        """Return a list of the level energies. """
        return [lvl.energy for lvl in self.levels]

    @property
    def total_number_gammas(self, ) -> int:
        return len((self.gammas.keys()))


class Database(_db.NuclideDatabase):
    """An object which holds the RIPL discrete levels information. """

    reader: object = read_ascii_file
    entry : object = LevelEntry
    writer: object = write_ascii_file
    parser: object = parse_ascii_dict

    def load(self, fpath: str) -> None:
        """Read discrete level information from an ASCII file. 
           Not this can be used in successive calls, each time adding potentially new entries to the database (self.data).
        """
        _data = type(self).reader(fpath)
        _data = type(self).parser(_data)
        for entry in _data.keys():
            self.data[entry] = _data[entry]

    def load_all(self, directory: str) -> None:
        """Read the entire database into memory from many individual ASCII files. """
        # Loop over the element numbers
        for Z in range(1, 119):
            fpath = _os.path.join(directory, local_data_dir, f"z{Z:03d}.dat")
            self.load(fpath)

    def load_element(self, Z: int, directory: str) -> None:
        """Read the database entry for element with proton number (Z). """
        fpath = _os.path.join(directory, local_data_dir, f"z{Z:03d}.dat")
        self.load(fpath)

    def levels(self, n: _c.Nuclide) -> object:
        """Return the levels for the particular nuclide in the database. """
        try:
            return self.data[n].levels
        except KeyError:
            return []

    @property
    def levels_dict(self, ) -> dict:
        """Return the levels in the current database as a dictionary accessible by the nucleus. """
        d = {}
        for n in self.data.keys():
            d[n] = self.data[n].levels
        return d

    @property
    def levels_list(self, ) -> list:
        """Return the levels in the current database as a list (independent of the nuclei). """
        lst = []
        for n in self.data.keys():
            for lvl in self.data[n].levels:
                lst.append(lvl)
        return lst

    def insert_level(self, n: _c.Nuclide, new_level: _core.DiscreteLevelInfo) -> None:
        """Insert a new level into the database for the provided nuclide. """
        # Check if nucleus is already in the database
        if n in self.data.keys():
            # Nucleus is in the database
            # Grab existing energies
            energies = [x.energy for x in self.levels(n)]
            # Check if the level is already in the database via the energy of the level
            new_energy = new_level.energy
            if new_level.energy in energies:
                # Replace the existing level with the new level
                i = energies.index(new_energy)
                self.data[n].levels[i] = new_level
            else:
                # Find the nearest energy level and insert there
                i = min(range(len(energies)), key=lambda x: abs(new_level.energy - energies[x]))
                self.data[n].levels.insert(i, new_level)
        else:
            # We can easily add since the nucleus is not in the database
            self.data[n] = LevelEntry(n, new_level, [])


def load_all(directory: str) -> Database:
    """Load and return the entire discrete level database. """
    d = Database()
    d.load_all(directory)
    return d


def load_element(Z: int, directory: str) -> Database:
    """Load and return the discrete level database for a given element (Z). """
    d = Database()
    d.load_element(Z, directory)
    return d
