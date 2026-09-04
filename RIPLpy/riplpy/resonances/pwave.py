
# ========================

# Dataclass
from dataclasses import dataclass as _dataclass

# Typing
from typing import ClassVar as _ClassVar

# Functools
from functools import partial as _partial

# RIPLpy
import riplpy.config as _config
import riplpy.collections as _c
import riplpy.db as _db

# ========================

# Path key of the database file as defined in config.py
FILE_PATH_KEY = 'resonances_pwave'
LOCAL_FILE_PATH = _config.get_data_file_path(FILE_PATH_KEY)

# ========================

# RIPL-4 resonance file format (S. Goriely, 18 August 2025):
#
#   (3i4, 2x, a2, 2x, i2, f5.1, 1x, a1, f10.3, 1p, 12e12.3)
#
# Columns: Z, N, A, Sym, L, J, P, Sn, D(RIPL-3), dD(RIPL-3), D(BNL), dD(BNL),
#          Gg(RIPL-3), dGg(RIPL-3), Gg(BNL), dGg(BNL), S(RIPL-3), dS(RIPL-3),
#          S(BNL), dS(BNL)
#
# Numeric fields may be left blank to indicate "no value".
FORTRAN_FORMAT = '(3i4,2x,a2,2x,i2,f5.1,1x,a1,f10.3,1p,12e12.3)'

FILE_HEADER = (
    "#  Z   N   A Sym   L  J   P     Sn      D(RIPL-3)    Err        D(BNL)       Err"
    "       Gg(RIPL-3)    Err        Gg(BNL)       Err       S(RIPL-3)     Err     "
    "   S(BNL)       Err\n"
    "# []  []  []  []  []  []  []   [MeV]      [eV]       [eV]        [eV]        [eV]"
    "        [eV]        [eV]        [eV]        [eV]        [1E-4]      [1E-4]    "
    "  [1E-4]      [1E-4]\n"
)

# Column slices (0-indexed, end-exclusive) corresponding to the FORTRAN format.
_COL_Z      = (0, 4)
_COL_N      = (4, 8)
_COL_A      = (8, 12)
_COL_SYM    = (14, 16)
_COL_L      = (18, 20)
_COL_J      = (20, 25)
_COL_P      = (26, 27)
_COL_SN     = (27, 37)
# 12 e12.3 numeric fields starting at col 37
_E_FIELDS = [(37 + 12 * i, 37 + 12 * (i + 1)) for i in range(12)]

# ========================

def _parse_float(value):
    """Parse a float value, returning None for empty/invalid values."""
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value) if value == value else None  # NaN check
    s = str(value).strip()
    if s == '' or s == '****' or s == '*':
        return None
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def _parse_int(value):
    """Parse an integer value, returning None for empty/invalid values."""
    if value is None:
        return None
    s = str(value).strip()
    if s == '':
        return None
    try:
        return int(s)
    except (ValueError, TypeError):
        return None


def _slice(line: str, span: tuple) -> str:
    """Safely slice a line within bounds."""
    start, stop = span
    return line[start:stop] if start < len(line) else ''


def read_ascii_file(fpath: str) -> dict:
    """Read resonance data from an ASCII file and return a dictionary.

    The file follows the RIPL-4 (Goriely 2025) format with header lines
    prefixed by ``#`` and fixed-width columns:

        (3i4, 2x, a2, 2x, i2, f5.1, 1x, a1, f10.3, 1p, 12e12.3)
    """
    # Placeholder dictionary
    d = {}

    # Read data from file
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            # Skip header lines and blank lines
            if not line.strip() or line.lstrip().startswith('#'):
                continue

            Z   = _parse_int(_slice(line, _COL_Z))
            N   = _parse_int(_slice(line, _COL_N))
            A   = _parse_int(_slice(line, _COL_A))
            sym = _slice(line, _COL_SYM).strip()
            L   = _parse_int(_slice(line, _COL_L))
            J   = _parse_float(_slice(line, _COL_J))
            P   = _slice(line, _COL_P).strip() or None
            Sn  = _parse_float(_slice(line, _COL_SN))

            # 12 numeric e12.3 fields
            vals = [_parse_float(_slice(line, span)) for span in _E_FIELDS]
            (D_R3,  dD_R3,  D_BNL,  dD_BNL,
             Gg_R3, dGg_R3, Gg_BNL, dGg_BNL,
             S_R3,  dS_R3,  S_BNL,  dS_BNL) = vals

            if Z is None or A is None:
                # Skip malformed rows
                continue

            n = _c.Nuclide(Z=Z, A=A)
            package = {
                'n':       n,
                'sym':     sym,
                'L':       L,
                'Io':      J,        # target ground-state spin
                'parity':  P,
                'Bn':      Sn,       # neutron binding energy [MeV]
                # Primary (RIPL-3) values exposed via the historical field names
                'D':       D_R3,     # average resonance spacing [eV]
                'Derr':    dD_R3,
                'Gam':     Gg_R3,    # average radiative width [eV]
                'Gerr':    dGg_R3,
                'Str':     S_R3,     # neutron strength function [1e-4]
                'Serr':    dS_R3,
                # BNL (Mughabghab 2018) values
                'D_BNL':    D_BNL,
                'D_BNL_err':  dD_BNL,
                'Gam_BNL':    Gg_BNL,
                'Gam_BNL_err':dGg_BNL,
                'Str_BNL':    S_BNL,
                'Str_BNL_err':dS_BNL,
            }
            d[n] = package
    return d


def _format_value(value, fmt_type, width, decimals=None):
    """Format a value for writing, handling None as blank spaces.

    Args:
        value: The value to format
        fmt_type: 'i' for integer, 'f' for float, 'e' for exponential, 'a' for string
        width: Total field width
        decimals: Number of decimal places (for float/exponential)
    """
    if value is None:
        return ' ' * width

    if fmt_type == 'i':
        return f"{int(value):>{width}d}"
    elif fmt_type == 'f':
        return f"{float(value):>{width}.{decimals}f}"
    elif fmt_type == 'e':
        return f"{float(value):>{width}.{decimals}E}"
    elif fmt_type == 'a':
        return f"{str(value):<{width}}"
    return str(value)


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write resonance data to an ASCII file in the RIPL-4 format."""

    with open(fpath, 'w', encoding='utf-8') as fp:
        # Write header
        fp.write(FILE_HEADER)
        # Loop over data
        for n in data.keys():
            entry = data[n]

            sym = getattr(entry, 'sym', None) or entry.n.element_symbol
            parity = getattr(entry, 'parity', None) or ' '

            line_parts = [
                _format_value(entry.n.Z, 'i', 4),       # Z (i4)
                _format_value(entry.n.N, 'i', 4),       # N (i4)
                _format_value(entry.n.A, 'i', 4),       # A (i4)
                '  ',                                    # 2x
                f"{sym:<2}",                            # Sym (a2)
                '  ',                                    # 2x
                _format_value(entry.L, 'i', 2),         # L (i2)
                _format_value(entry.Io, 'f', 5, 1),     # J (f5.1)
                ' ',                                     # 1x
                f"{parity[:1]:<1}",                     # P (a1)
                _format_value(entry.Bn, 'f', 10, 3),    # Sn (f10.3)
                _format_value(entry.D,           'e', 12, 3),
                _format_value(entry.Derr,        'e', 12, 3),
                _format_value(entry.D_BNL,       'e', 12, 3),
                _format_value(entry.D_BNL_err,   'e', 12, 3),
                _format_value(entry.Gam,         'e', 12, 3),
                _format_value(entry.Gerr,        'e', 12, 3),
                _format_value(entry.Gam_BNL,     'e', 12, 3),
                _format_value(entry.Gam_BNL_err, 'e', 12, 3),
                _format_value(entry.Str,         'e', 12, 3),
                _format_value(entry.Serr,        'e', 12, 3),
                _format_value(entry.Str_BNL,     'e', 12, 3),
                _format_value(entry.Str_BNL_err, 'e', 12, 3),
            ]
            fp.write(''.join(line_parts) + "\n")


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the neutron resonance database.

    Contains s-wave or p-wave neutron resonance parameters from
    evaluated experimental data (RIPL-3 and BNL/Mughabghab 2018).

    The primary fields (``D``, ``Derr``, ``Gam``, ``Gerr``, ``Str``, ``Serr``)
    carry the RIPL-3 recommended values; the ``*_BNL`` fields carry the
    Mughabghab (2018) BNL evaluation.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        sym: Element symbol as it appears in the data file
        L: Angular momentum of the incident neutron (0 = s-wave, 1 = p-wave)
        Io: Spin of the target ground state
        parity: Parity of the target ground state ('+' or '-')
        Bn: Neutron binding energy for the compound nucleus [MeV]
        D: Average resonance spacing from RIPL-3 [eV]
        Derr: Uncertainty on D [eV]
        Gam: Average radiative width from RIPL-3 [eV]
        Gerr: Uncertainty on Gam [eV]
        Str: Neutron strength function from RIPL-3 [10^-4]
        Serr: Uncertainty on Str [10^-4]
        D_BNL: Average resonance spacing from BNL [eV]
        D_BNL_err: Uncertainty on D_BNL [eV]
        Gam_BNL: Average radiative width from BNL [eV]
        Gam_BNL_err: Uncertainty on Gam_BNL [eV]
        Str_BNL: Neutron strength function from BNL [10^-4]
        Str_BNL_err: Uncertainty on Str_BNL [10^-4]
    """

    n           : _c.Nuclide = None
    sym         : str        = None
    L           : int        = None
    Io          : float      = None
    parity      : str        = None
    Bn          : float      = None
    D           : float      = None
    Derr        : float      = None
    Gam         : float      = None
    Gerr        : float      = None
    Str         : float      = None
    Serr        : float      = None
    D_BNL       : float      = None
    D_BNL_err   : float      = None
    Gam_BNL     : float      = None
    Gam_BNL_err : float      = None
    Str_BNL     : float      = None
    Str_BNL_err : float      = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n':           'Target nucleus',
        'sym':         'Element symbol of target nucleus',
        'L':           'Angular momentum of incident neutron (0=s-wave, 1=p-wave)',
        'Io':          'Spin of the target ground state',
        'parity':      'Parity of the target ground state',
        'Bn':          'Neutron binding energy [MeV]',
        'D':           'Average resonance spacing (RIPL-3) [eV]',
        'Derr':        'Uncertainty on D (RIPL-3) [eV]',
        'Gam':         'Average radiative width (RIPL-3) [eV]',
        'Gerr':        'Uncertainty on Gam (RIPL-3) [eV]',
        'Str':         'Neutron strength function (RIPL-3) [10^-4]',
        'Serr':        'Uncertainty on Str (RIPL-3) [10^-4]',
        'D_BNL':       'Average resonance spacing (BNL) [eV]',
        'D_BNL_err':   'Uncertainty on D_BNL [eV]',
        'Gam_BNL':     'Average radiative width (BNL) [eV]',
        'Gam_BNL_err': 'Uncertainty on Gam_BNL [eV]',
        'Str_BNL':     'Neutron strength function (BNL) [10^-4]',
        'Str_BNL_err': 'Uncertainty on Str_BNL [10^-4]',
    }

    @property
    def as_tuple(self) -> tuple:
        """Tuple representation matching the on-disk field order."""
        return (
            self.n.Z, self.n.N, self.n.A, self.n.element_symbol,
            self.L, self.Io, self.parity, self.Bn,
            self.D, self.Derr, self.D_BNL, self.D_BNL_err,
            self.Gam, self.Gerr, self.Gam_BNL, self.Gam_BNL_err,
            self.Str, self.Serr, self.Str_BNL, self.Str_BNL_err,
        )


class Database(_db.NuclideDatabase):

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_ascii_file


# Create the local 'load' function with config, database object, and file_path_key pre-filled
load = _partial(_db.loader, config=_config, db_obj=Database, file_path_key=FILE_PATH_KEY)
