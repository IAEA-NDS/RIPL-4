# Enhanced Generalized Superfluid Model (EGSM)

"""
Level density parameters for the Enhanced Generalized Superfluid Model, 
which combines the superfluid model below BCS critical energy and  the
Fermi gas model above. The spin distribution in the Fermi gas part is
calculated by reducing excitation energy by subtracting rotational energy.
Nuclear deformation is used to determine the latter. EGSM
takes into account collective enhancement of nuclear level densities
in addition to shell and superfluid effects. The parameters were
obtained by fitting the corresponding model formulas to the RIPL-3
recommended spacings of s-wave neutron resonances.  Matching of the
level densities to the low-lying levels must be performed by the user. 
To this end, a dedicated option is provided in the EGSM code.

Quantities Sym, Io, Bn, Do, dDo, Esh, +da, and -da are only provided 
for convenience and are not used in the calculations. The EGSM code
reads shell corrections, nuclear deformations, and neutron bindings
directly from the respective RIPL files.
"""

# Example usage:
# d = read_ascii_file("data/level-densities-egsm.dat")
# print(d)

# ========================

# Dataclasses
from dataclasses import dataclass as _dataclass

# Typing
from typing import ClassVar as _ClassVar

# Functools
from functools import partial as _partial

# Fortran format
import fortranformat as _ff

# RIPLpy
import riplpy.config as _config
import riplpy.db as _db
import riplpy.collections as _c

# ========================

__all__ = ('FILE_PATH_KEY', 'LOCAL_FILE_PATH', 'FILE_HEADER', 'FORTRAN_FORMAT', 'reader', 'writer', 'read_ascii_file', 'write_ascii_file', 'Entry', 'Database', 'load')

# ========================

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'densities_egsm'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

FILE_HEADER = (
    "#  Fit of EGSM a-parameter to experimental Do\n"
    "#  Clasic vibr enhan., RIPL-4 Do & Shell corr.\n"
    "#  Z  A  Sym  J   Qn         Dobs          dDobs        Shell         Def          Dcalc        +da_exp  a_exp -da_exp  a_sys    a_exp/a_sys\n"
    "#------------------------------------------------------------------------------------------------------------------------------------------\n"
)
# The RIPL-4 github release reformatted this file (added Def, Dcalc, a_sys
# columns). We use whitespace-tokenised parsing so we remain tolerant of small
# layout tweaks; ``fortranformat`` is kept imported only for symmetry with
# other section readers.
FORTRAN_FORMAT = None
reader = None
writer = None

# ========================


def _parse_line(line: str) -> dict:
    """Parse a single EGSM data line into a field dict."""
    tokens = line.split()
    # Expected 15 columns in the RIPL-4 github release format
    if len(tokens) < 13:
        raise ValueError(f"Unexpected EGSM record: {line!r}")
    Z = int(tokens[0])
    A = int(tokens[1])
    sym = tokens[2]
    J = float(tokens[3])
    Qn = float(tokens[4])
    Dobs = float(tokens[5])
    dDobs = float(tokens[6])
    Shell = float(tokens[7])
    # Older releases had no Def/Dcalc/a_sys columns. Branch on length.
    if len(tokens) >= 15:
        Def = float(tokens[8])
        Dcalc = float(tokens[9])
        dap = float(tokens[10])
        a = float(tokens[11])
        dam = float(tokens[12])
        a_sys = float(tokens[13])
        a_ratio = float(tokens[14])
    else:
        Def = None
        Dcalc = None
        dap = float(tokens[8])
        a = float(tokens[9])
        dam = float(tokens[10])
        a_sys = None
        a_ratio = None
    return {
        'Z': Z, 'A': A, 'sym': sym, 'Io': J, 'Bn': Qn,
        'Do': Dobs, 'Derr': dDobs, 'Esh': Shell,
        'Def': Def, 'Dcalc': Dcalc,
        'dap': dap, 'a': a, 'dam': dam,
        'a_sys': a_sys, 'a_ratio': a_ratio,
    }


def read_ascii_file(fpath: str) -> dict:
    """Read Enhanced Generalized Superfluid Model (EGSM) data from an ASCII formatted file. """

    d = {}
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            if not line.strip():
                continue
            if line.lstrip().startswith('#'):
                continue
            fields = _parse_line(line)
            n = _c.Nuclide(Z=fields['Z'], A=fields['A'])
            package = {
                'n': n,
                'Io': fields['Io'],
                'Bn': fields['Bn'],
                'Do': fields['Do'],
                'Derr': fields['Derr'],
                'Esh': fields['Esh'],
                'Def': fields['Def'],
                'Dcalc': fields['Dcalc'],
                'dap': fields['dap'],
                'a': fields['a'],
                'dam': fields['dam'],
                'a_sys': fields['a_sys'],
                'a_ratio': fields['a_ratio'],
            }
            d[n] = package
    return d


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write an Enhanced Generalized Superfluid Model (EGSM) ASCII file. """

    with open(fpath, 'w', encoding='utf-8') as fp:
        fp.write(FILE_HEADER)
        for n in data.keys():
            e = data[n]
            fp.write(
                f"  {e.n.Z:2d}  {e.n.A:3d} {e.n.element_symbol:<2s}  {e.Io:3.1f}  "
                f"{e.Bn:6.3f}  "
                f"{e.Do:12.5E}  {e.Derr:12.5E}  {e.Esh:12.5E}  "
                f"{(e.Def if e.Def is not None else 0.0):12.5E}  "
                f"{(e.Dcalc if e.Dcalc is not None else 0.0):12.5E}  "
                f"{e.dap:7.4f}  {e.a:6.4f}  {e.dam:6.4f}  "
                f"{(e.a_sys if e.a_sys is not None else 0.0):6.3f}  "
                f"{(e.a_ratio if e.a_ratio is not None else 0.0):6.3f}\n"
            )


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the Enhanced Generalized Superfluid Model (EGSM) database.

    Contains level density parameters for the EGSM model which combines
    the superfluid model below BCS critical energy with the Fermi gas model.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        Io: Ground state spin of target nucleus
        Bn: Neutron binding energy [MeV]
        Do: Experimental s-wave resonance spacing [keV]
        Derr: Uncertainty on resonance spacing [keV]
        Esh: Shell correction energy [MeV]
        dap: Upper uncertainty on level density parameter [MeV^-1]
        a: Level density parameter at Bn [MeV^-1]
        dam: Lower uncertainty on level density parameter [MeV^-1]
    """

    n      : _c.Nuclide = None
    Io     : float      = None
    Bn     : float      = None
    Do     : float      = None
    Derr   : float      = None
    Esh    : float      = None
    Def    : float      = None
    Dcalc  : float      = None
    dap    : float      = None
    a      : float      = None
    dam    : float      = None
    a_sys  : float      = None
    a_ratio: float      = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'Io': 'Ground state spin of target nucleus',
        'Bn': 'Neutron binding energy / Qn [MeV]',
        'Do': 'Experimental s-wave resonance spacing Dobs [keV]',
        'Derr': 'Uncertainty on resonance spacing dDobs [keV]',
        'Esh': 'Shell correction energy [MeV]',
        'Def': 'Ground-state deformation parameter',
        'Dcalc': 'Calculated s-wave resonance spacing [keV]',
        'dap': 'Upper uncertainty on level density parameter [MeV^-1]',
        'a': 'Level density parameter at Bn [MeV^-1]',
        'dam': 'Lower uncertainty on level density parameter [MeV^-1]',
        'a_sys': 'Systematic level density parameter a_sys [MeV^-1]',
        'a_ratio': 'Ratio a_exp / a_sys',
    }


class Database(_db.NuclideDatabase):
    """The database object for the Enhanced Generalized Superfluid Model (EGSM) level density model. """

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
