# -*- coding: utf-8 -*-
"""Experimental GDR parameter fits within the SLO and SMLO frameworks (RIPL-4).

The RIPL-4 github release replaces the older ``gdr-parameters-exp-SLO.dat`` /
``gdr-parameters-exp-MLO.dat`` files with a redesigned compilation by
V.A. Plujko, O.M. Gorbachenko, K.M. Solodovnyk (TSNU), and R. Capote &
P. Dimitriou (IAEA). The new directory ``gamma/gdr_parameters_exp_new/``
contains four data files:

* ``gdr_parameters_recommended_exp_slo.dat``  -- single-line records, recommended values
* ``gdr_parameters_recommended_exp_smlo.dat`` -- single-line records, recommended values
* ``gdr_parameters&errors_exp_slo.dat``       -- two-line records (values + 1-sigma errors)
* ``gdr_parameters&errors_exp_smlo.dat``      -- two-line records (values + 1-sigma errors)

Each entry exposes the GDR Lorentzian peak energies ``Er1, Er2`` [MeV], widths
``Wr1, Wr2`` [MeV], strengths ``S1, S2`` (in units of the TRK sum rule), peak
cross sections ``CSp1, CSp2`` [mb], total strength ``Sr=S1+S2`` [TRK units],
and the experimental energy range plus reference.

The ``MLO`` API is retained as an alias of ``SMLO`` for backwards
compatibility, since RIPL-4 no longer distributes MLO-only fits.

"""

# ========================

# OS
import os as _os

# Logging
import logging as _logging

# RIPLpy
import riplpy.collections as _c
import riplpy.config as _config
import riplpy.db as _db
from riplpy.exceptions import NucleusNotFoundError as _NucleusNotFoundError
from riplpy.config import resolve_directory as _resolve_directory
from riplpy.gamma import _layout as _layout
from . import core as _core

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = (
    'slo_local_file_path', 'smlo_local_file_path', 'mlo_local_file_path',
    'slo_errors_local_file_path', 'smlo_errors_local_file_path',
    'SLO_FILE_HEADER', 'SMLO_FILE_HEADER', 'MLO_FILE_HEADER',
    'read_recommended_file', 'read_errors_file', 'read_errors_records',
    'write_slo_ascii_file', 'write_smlo_ascii_file',
    'write_slo_errors_ascii_file', 'write_smlo_errors_ascii_file',
    'SLO_ERRORS_FILE_HEADER', 'SMLO_ERRORS_FILE_HEADER',
    'Entry', 'SLO_Database', 'SMLO_Database', 'MLO_Database',
    'load_slo', 'load_smlo', 'load_mlo',
)

# ========================

# Recommended (single-line) files
# RIPL-4 keeps these flat at gamma/; reorganised trees nest them under
# gamma/gdr_parameters_exp_new/. Both are accepted, flat first.
_NESTED = 'gdr_parameters_exp_new'


def _candidates(fname: str) -> tuple[str, str]:
    """Return the flat RIPL-4 path and the nested legacy path for ``fname``."""
    return (_os.path.join('gamma', fname),
            _os.path.join('gamma', _NESTED, fname))


slo_local_file_path  = _os.path.join('gamma', 'gdr_parameters_recommended_exp_slo.dat')
smlo_local_file_path = _os.path.join('gamma', 'gdr_parameters_recommended_exp_smlo.dat')
# Backwards-compat alias: MLO == SMLO in RIPL-4
mlo_local_file_path  = smlo_local_file_path

# Errors files (two-line records)
slo_errors_local_file_path  = _os.path.join('gamma', 'gdr_parameters&errors_exp_slo.dat')
smlo_errors_local_file_path = _os.path.join('gamma', 'gdr_parameters&errors_exp_smlo.dat')

_SLO_CANDIDATES         = _candidates('gdr_parameters_recommended_exp_slo.dat')
_SMLO_CANDIDATES        = _candidates('gdr_parameters_recommended_exp_smlo.dat')
_MLO_CANDIDATES         = _SMLO_CANDIDATES
_SLO_ERRORS_CANDIDATES  = _candidates('gdr_parameters&errors_exp_slo.dat')
_SMLO_ERRORS_CANDIDATES = _candidates('gdr_parameters&errors_exp_smlo.dat')

# ========================

SLO_FILE_HEADER = (
    "#     Experimental values of recommended GDR parameters (SLO)\n"
    "#  Z  A  El     Er1     Wr1     S1      CSp1     Er2     Wr2     S2     CSp2      S    E_{min}-E_{max} Ref.\n"
    "#-----------------------------------------------------------------------------------------------\n"
)
SMLO_FILE_HEADER = SLO_FILE_HEADER.replace('(SLO)', '(SMLO)')
MLO_FILE_HEADER  = SMLO_FILE_HEADER  # alias

# The FORMAT line is what SLO_Database.load() uses to recognise the layout.
SLO_ERRORS_FILE_HEADER = (
    "#     Experimental values and uncertainties of GDR parameters (SLO)\n"
    "#     Second line for each isotope is uncertainties of\n"
    "#     quantities in the nearest upper line.\n"
    "#     FORTRAN format of data reading for one entry is\n"
    "#     FORMAT(2I4,1X,A2,A3,1X,9F8.3,A22,/,15X,9F8.3)\n"
    "#  Z  A  El Id     Er1     Wr1     S1      CSp1     Er2     Wr2     S2     CSp2      S    E_{min}-E_{max} Ref.\n"
    "#                 dEr1    dWr1    dS1     dCSp1    dEr2    dWr2    dS2    dCSp2     dS\n"
    "#-------------------------------------------------------------------------------------------------------------\n"
)
SMLO_ERRORS_FILE_HEADER = SLO_ERRORS_FILE_HEADER.replace('(SLO)', '(SMLO)')

# Fortran-style column layout for the recommended file:
#   (2I4, 1X, A2, 2X, 9F8.3, A22)
#   Z  A  El   Er1  Wr1  S1  CSp1  Er2  Wr2  S2  CSp2  S   Emin-Emax Ref
# The errors file inserts an A3 Id between El and the 9 floats and is split
# over two lines; the second line gives uncertainties on the floats from line 1.
# Some entries leave the second (deformed) component blank, which we parse as
# zeros and convert to ``None`` so downstream code can detect "spherical only".


def _safe_float(s: str) -> float | None:
    """Parse a float from a fixed-width field, returning None for blank fields."""
    s = s.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _parse_recommended_line(line: str) -> dict:
    """Parse one line of the RIPL-4 'recommended' SLO/SMLO data files.

    Format: FORMAT(2I4, 1X, A2, 2X, 9F8.3, A22)
    """
    # Field slices match the Fortran format
    Z   = int(line[0:4])
    A   = int(line[4:8])
    El  = line[9:11].strip()  # noqa: F841 (we recompute from Z)
    Er1  = _safe_float(line[13:21])
    Wr1  = _safe_float(line[21:29])
    Sr1  = _safe_float(line[29:37])
    CSp1 = _safe_float(line[37:45])
    Er2  = _safe_float(line[45:53])
    Wr2  = _safe_float(line[53:61])
    Sr2  = _safe_float(line[61:69])
    CSp2 = _safe_float(line[69:77])
    Sr   = _safe_float(line[77:85])
    tail = line[85:].strip() if len(line) > 85 else ''

    Emin = Emax = None
    reference = None
    # Tail looks like " 17.2 - 23.7  1986Var"
    if tail:
        parts = tail.split()
        try:
            Emin = float(parts[0])
        except (ValueError, IndexError):
            pass
        if len(parts) >= 3:
            try:
                Emax = float(parts[2])
            except ValueError:
                pass
        if len(parts) >= 4:
            reference = parts[3]
        elif parts and not parts[-1].replace('.', '').replace('-', '').isdigit():
            reference = parts[-1]

    # A == 0 marks a natural-isotopic-composition row (readme). Keep it under
    # Nuclide(Z, A=0) so it cannot collide with a real isotope, and tag the Id.
    n = _c.Nuclide(Z=Z, A=A)
    Id = 'nat' if A == 0 else None

    return {
        'n': n, 'Id': Id,
        'Er1': Er1, 'Wr1': Wr1, 'Sr1': Sr1, 'CSp1': CSp1,
        'Er2': Er2, 'Wr2': Wr2, 'Sr2': Sr2, 'CSp2': CSp2,
        'Sr': Sr, 'Emin': Emin, 'Emax': Emax, 'reference': reference,
    }


def _parse_errors_first_line(line: str) -> dict:
    """Parse the first (values) line of a two-line errors record.

    Format: FORMAT(2I4, 1X, A2, A3, 1X, 9F8.3, A22)
    """
    Z   = int(line[0:4])
    A   = int(line[4:8])
    El  = line[9:11].strip()  # noqa: F841
    Id  = line[11:14].strip() or None
    Er1  = _safe_float(line[15:23])
    Wr1  = _safe_float(line[23:31])
    Sr1  = _safe_float(line[31:39])
    CSp1 = _safe_float(line[39:47])
    Er2  = _safe_float(line[47:55])
    Wr2  = _safe_float(line[55:63])
    Sr2  = _safe_float(line[63:71])
    CSp2 = _safe_float(line[71:79])
    Sr   = _safe_float(line[79:87])
    tail = line[87:].strip() if len(line) > 87 else ''

    Emin = Emax = None
    reference = None
    if tail:
        parts = tail.split()
        try:
            Emin = float(parts[0])
        except (ValueError, IndexError):
            pass
        if len(parts) >= 3:
            try:
                Emax = float(parts[2])
            except ValueError:
                pass
        if len(parts) >= 4:
            reference = parts[3]
        elif parts:
            reference = parts[-1]

    n = _c.Nuclide(Z=Z, A=A)  # A == 0: natural composition, see above
    if A == 0 and Id is None:
        Id = 'nat'

    return {
        'n': n, 'Id': Id,
        'Er1': Er1, 'Wr1': Wr1, 'Sr1': Sr1, 'CSp1': CSp1,
        'Er2': Er2, 'Wr2': Wr2, 'Sr2': Sr2, 'CSp2': CSp2,
        'Sr': Sr, 'Emin': Emin, 'Emax': Emax, 'reference': reference,
    }


def _parse_errors_second_line(line: str) -> dict:
    """Parse the second (uncertainties) line of a two-line errors record.

    Format: FORMAT(15X, 9F8.3) — column 16 onwards holds 9 fields.
    """
    # Skip leading 15 columns
    body = line[15:]
    fields = []
    for i in range(9):
        seg = body[i * 8:(i + 1) * 8]
        fields.append(_safe_float(seg))
    dEr1, dWr1, dSr1, dCSp1, dEr2, dWr2, dSr2, dCSp2, dSr = fields
    return {
        'dEr1': dEr1, 'dWr1': dWr1, 'dSr1': dSr1, 'dCSp1': dCSp1,
        'dEr2': dEr2, 'dWr2': dWr2, 'dSr2': dSr2, 'dCSp2': dCSp2,
        'dSr': dSr,
    }


def read_recommended_file(fpath: str) -> dict:
    """Read a single-line RIPL-4 recommended SLO/SMLO data file."""
    d = {}
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            try:
                pkt = _parse_recommended_line(line.rstrip('\n'))
            except (ValueError, IndexError) as exc:
                _logger.debug(f"Skipping malformed line in {fpath}: {exc}")
                continue
            d[pkt['n']] = pkt
    return d


def read_errors_records(fpath: str) -> list:
    """Read every record of a two-line RIPL-4 SLO/SMLO file, in file order.

    The errors files hold one record per nuclide and type of experimental
    data (the Id column), so a nuclide can appear several times. The readme
    calls the first record for each nuclide the recommended one.
    """
    records = []
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        # Pre-filter: drop blank and comment lines, keep order
        lines = [ln.rstrip('\n') for ln in fp if ln.strip() and not ln.lstrip().startswith('#')]
    i = 0
    while i < len(lines):
        try:
            packet = _parse_errors_first_line(lines[i])
        except (ValueError, IndexError) as exc:
            _logger.debug(f"Skipping malformed values line in {fpath}: {exc}")
            i += 1
            continue
        if i + 1 < len(lines):
            try:
                packet.update(_parse_errors_second_line(lines[i + 1]))
            except (ValueError, IndexError) as exc:
                _logger.debug(f"Skipping malformed errors line in {fpath}: {exc}")
        records.append(packet)
        i += 2
    return records


def read_errors_file(fpath: str) -> dict:
    """Read a two-line RIPL-4 SLO/SMLO file, keeping the recommended record.

    Returns one record per nuclide: the first one in the file, which the
    readme calls recommended. Use :func:`read_errors_records` for all of them.
    """
    d = {}
    for packet in read_errors_records(fpath):
        d.setdefault(packet['n'], packet)
    return d


# Backwards-compat: the original code exposed a single ``read_ascii_file``.
read_ascii_file = read_recommended_file


def _format_recommended_line(entry) -> str:
    """Render an entry in the recommended (single-line) layout."""
    def f8(x):
        return f"{x:8.3f}" if x is not None else "        "

    A = entry.n.A
    Z = entry.n.Z
    El = entry.n.element_symbol
    parts = [
        f"{Z:4d}{A:4d} {El:<2s}  ",
        f8(entry.Er1), f8(entry.Wr1), f8(entry.Sr1), f8(entry.CSp1),
        f8(entry.Er2), f8(entry.Wr2), f8(entry.Sr2), f8(entry.CSp2),
        f8(entry.Sr),
    ]
    base = "".join(parts)
    emin = f"{entry.Emin:5.1f}" if entry.Emin is not None else "     "
    emax = f"{entry.Emax:5.1f}" if entry.Emax is not None else "     "
    ref  = entry.reference if entry.reference else ""
    return f"{base} {emin} - {emax}  {ref}"


def _format_errors_lines(entry) -> str:
    """Render a fit in the errors (two-line) layout, values then 1-sigma.

    FORMAT(2I4,1X,A2,A3,1X,9F8.3,A22,/,15X,9F8.3)
    """
    def f8(x):
        return f"{x:8.3f}" if x is not None else "        "

    # The A3 Id field holds a one- or two-character Id after a blank, as in
    # the RIPL-4 files (" 4 ", " 3a", " 12").
    Id = entry.Id if entry.Id and entry.Id != 'nat' else ''
    Id3 = (' ' + Id).ljust(3) if len(Id) <= 2 else Id[:3]
    head = f"{entry.n.Z:4d}{entry.n.A:4d} {entry.n.element_symbol:<2s}{Id3} "
    values = "".join(f8(getattr(entry, k)) for k in
                     ('Er1', 'Wr1', 'Sr1', 'CSp1', 'Er2', 'Wr2', 'Sr2', 'CSp2', 'Sr'))
    errors = "".join(f8(getattr(entry, k)) for k in
                     ('dEr1', 'dWr1', 'dSr1', 'dCSp1', 'dEr2', 'dWr2', 'dSr2', 'dCSp2', 'dSr'))
    emin = f"{entry.Emin:5.1f}" if entry.Emin is not None else "     "
    emax = f"{entry.Emax:5.1f}" if entry.Emax is not None else "     "
    ref  = entry.reference if entry.reference else ""
    first = f"{head}{values} {emin} - {emax}  {ref}".rstrip()
    second = (" " * 15 + errors).rstrip()
    return f"{first}\n{second}"


def _write_errors_file(fpath: str, fits: dict, header: str) -> None:
    with open(fpath, 'w', encoding='utf-8') as fp:
        fp.write(header)
        for n, entries in fits.items():
            for entry in entries:
                fp.write(_format_errors_lines(entry) + "\n")


def write_slo_errors_ascii_file(fpath: str, fits: dict) -> None:
    """Write every fit to an errors-style (two-line) SLO ASCII file.

    Args:
        fits: Mapping nuclide -> list of entries, as ``SLO_Database.fits``.
    """
    _write_errors_file(fpath, fits, SLO_ERRORS_FILE_HEADER)


def write_smlo_errors_ascii_file(fpath: str, fits: dict) -> None:
    """Write every fit to an errors-style (two-line) SMLO ASCII file."""
    _write_errors_file(fpath, fits, SMLO_ERRORS_FILE_HEADER)


def write_slo_ascii_file(fpath: str, data: dict) -> None:
    """Write entries to a recommended-style SLO ASCII file."""
    with open(fpath, 'w', encoding='utf-8') as fp:
        fp.write(SLO_FILE_HEADER)
        for n, entry in data.items():
            fp.write(_format_recommended_line(entry) + "\n")


def write_smlo_ascii_file(fpath: str, data: dict) -> None:
    """Write entries to a recommended-style SMLO ASCII file."""
    with open(fpath, 'w', encoding='utf-8') as fp:
        fp.write(SMLO_FILE_HEADER)
        for n, entry in data.items():
            fp.write(_format_recommended_line(entry) + "\n")


# Backwards-compat alias
write_mlo_ascii_file = write_smlo_ascii_file


class Entry(_core.Experimental_GDR_Parameter, _db.NuclideDatabaseEntry):
    """An experimental GDR fit entry (SLO or SMLO)."""
    pass


class SLO_Database(_db.NuclideDatabase):
    """Experimental GDR fits, keyed by nuclide.

    ``data[n]`` (and ``get(n)``) is the recommended fit for nuclide ``n``.
    ``fits[n]`` (and ``get_fits(n)``) lists every fit for ``n`` in file order,
    recommended first. The recommended files give one fit per nuclide; the
    errors files give one per type of experimental data (``Entry.Id``).
    Natural-composition rows are keyed as ``Nuclide(Z, A=0)``.
    """

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = write_slo_ascii_file
    errors_writer: object = write_slo_errors_ascii_file

    def __init__(self, data: dict | None = None) -> None:
        super().__init__(data)
        self.fits = {n: [e] for n, e in self.data.items()}
        # True when loaded from an errors file; save() then writes that layout.
        self.errors_layout = False

    def save(self, fpath: str) -> None:
        """Save in the layout the database was loaded from.

        An errors-file database is written in the two-line layout with every
        fit and its uncertainties; otherwise the one-line recommended layout.
        """
        if self.errors_layout:
            type(self).errors_writer(fpath, self.fits)
        else:
            super().save(fpath)

    def get_fits(self, n) -> list:
        """Return every fit for nuclide ``n``, recommended first.

        Raises:
            NucleusNotFoundError: If ``n`` has no fit in this database.
        """
        if n not in self.fits:
            raise _NucleusNotFoundError(f"Nuclide {n!r} not found in database.")
        return list(self.fits[n])

    def filter(self, predicate=None, **kwargs) -> "SLO_Database":
        """As :meth:`Database.filter`, keeping every fit of the kept nuclides."""
        result = super().filter(predicate, **kwargs)
        result.fits = {n: list(self.fits[n]) for n in result.data if n in self.fits}
        result.errors_layout = self.errors_layout
        return result

    def load(self, fpath: str) -> None:
        """Load entries from either a recommended (single-line) or errors
        (two-line) data file. The format is auto-detected by header content
        and column count.
        """
        # Try the errors format first if the file has the 'dEr1' marker
        with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
            head = fp.read(2048)
        if '/' in head and 'F8.3' in head and ',/,' in head.replace(' ', ''):
            self.errors_layout = True
        elif 'errors' in head.lower() and 'second line' in head.lower():
            self.errors_layout = True
        if self.errors_layout:
            records = read_errors_records(fpath)
        else:
            records = list(read_recommended_file(fpath).values())
        for pkt in records:
            entry = self.entry(**pkt)
            self.fits.setdefault(entry.n, []).append(entry)
            # The first fit for a nuclide is the recommended one (readme).
            self.data.setdefault(entry.n, entry)


class SMLO_Database(SLO_Database):

    writer: object = write_smlo_ascii_file
    errors_writer: object = write_smlo_errors_ascii_file


# Backwards-compatible alias
MLO_Database = SMLO_Database


def _load_db(db_cls, directory: str = None, file_path: str = None,
             rel_path: str = None) -> "_db.Database":
    """Internal helper: instantiate ``db_cls`` and load the requested file."""
    db = db_cls()
    db_path = file_path
    if db_path is None:
        directory = _resolve_directory(directory)
        if rel_path is not None:
            candidates = (rel_path,) if isinstance(rel_path, str) else rel_path
            db_path = (_layout.resolve_data_file(directory, *candidates)
                       or _os.path.join(directory, candidates[0]))
    if db_path is None or not _os.path.exists(db_path):
        _logger.warning(
            f"{db_cls.__name__}: data file not found at {db_path!r}; "
            "returning empty database."
        )
        return db
    db.load(db_path)
    return db


def load_slo(directory: str = None, file_path: str = None,
             errors: bool = False) -> SLO_Database:
    """Load the RIPL-4 SLO experimental GDR parameter database.

    Args:
        directory: Path to the RIPL root directory.
        file_path: Direct path to an SLO data file (overrides ``directory``).
        errors:    If True, load the two-line file with 1-sigma uncertainties.
                   Defaults to the single-line recommended file.
    """
    rel = _SLO_ERRORS_CANDIDATES if errors else _SLO_CANDIDATES
    return _load_db(SLO_Database, directory=directory, file_path=file_path,
                    rel_path=rel)


def load_smlo(directory: str = None, file_path: str = None,
              errors: bool = False) -> SMLO_Database:
    """Load the RIPL-4 SMLO experimental GDR parameter database."""
    rel = _SMLO_ERRORS_CANDIDATES if errors else _SMLO_CANDIDATES
    return _load_db(SMLO_Database, directory=directory, file_path=file_path,
                    rel_path=rel)


def load_mlo(directory: str = None, file_path: str = None,
             errors: bool = False) -> SMLO_Database:
    """Backwards-compatible alias for :func:`load_smlo`.

    RIPL-4 no longer distributes MLO-only fits; the modern compilation uses
    the Simplified Modified Lorentzian (SMLO) approach instead.
    """
    return load_smlo(directory=directory, file_path=file_path, errors=errors)
