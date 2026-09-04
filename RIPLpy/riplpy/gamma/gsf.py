
# -*- coding: utf-8 -*-
"""Python objects which provide access to the gamma-ray strength function (GSF) database of the Reference Input Parameter Library (RIPL).

The legacy ``gamma/gamma-strength-micro`` per-nucleus directory is not shipped
in the RIPL-4 GitHub release; ``load_all``/``load_nucleus`` transparently fall
back to the RIPL-4 D1M+QRPA tables under ``gamma/d1m/``.

``fE1`` units depend on the data source actually loaded:

* Legacy ``gamma-strength-micro/z<NNN>.dat`` files (RIPL-3): ``fE1`` is the
  E1 photoabsorption cross section in **mb/MeV** (the file header reads
  ``U[MeV]  fE1[mb/MeV]`` and the .readme confirms it).
* RIPL-4 D1M+QRPA fallback: ``fE1`` is the photon strength function in
  **MeV^-3**.

The reader passes the on-disk values through verbatim; convert between the
two conventions yourself if you need apples-to-apples comparisons.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load the GSF data for a single nucleus (D1M fallback on GitHub)::
       $ import riplpy.gamma as gamma
       $ from riplpy.collections import Nuclide
       $ db = gamma.gsf.load_nucleus(Nuclide(Z=26, A=56))
       $ entry = db.get(Nuclide(Z=26, A=56))
       $ print(entry['U'][:3], entry['fE1'][:3])

"""

# ========================

# OS
import os as _os

# Logging
import logging as _logging

# Fortran format
import fortranformat as _ff

# RIPLpy
import riplpy.collections as _c
import riplpy.db as _db
from riplpy.config import resolve_directory as _resolve_directory
from . import core as _core

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = ('local_data_dir', 'FILE_HEADER', 'FORTRAN_FORMAT', 'reader', 'writer', 'read_ascii_file', 'write_ascii_file', 'Entry', 'Database', 'load')

# ========================

# Relative path to database folder
local_data_dir = _os.path.join("gamma", "gamma-strength-micro")

# ========================

FILE_HEADER = "  U[MeV]  fE1[mb/MeV]\n"
FORTRAN_FORMAT = '(f9.3,e12.3)'
reader   = _ff.FortranRecordReader(FORTRAN_FORMAT)
writer   = _ff.FortranRecordWriter(FORTRAN_FORMAT)

# ========================

def read_ascii_file(fpath: str) -> dict:
    """Read data from an ASCII file and return a dictionary. """
    d = {}
    current_nuclide = None
    package = None

    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            line_stripped = line.strip()

            # Skip empty lines and comment lines
            if not line_stripped or line_stripped.startswith('#'):
                continue

            # Check for nucleus header: "Z= XX A= YY"
            if "Z=" in line and "A=" in line:
                # Save previous isotope if exists
                if current_nuclide is not None and package is not None:
                    d[current_nuclide] = package

                # Parse new nucleus
                parts = line.split()
                Z = int(parts[1])
                A = int(parts[3])
                current_nuclide = _c.Nuclide(z=Z, a=A)
                package = {'U': [], 'fE1': []}
                continue

            # Skip header line
            if "U[MeV]" in line or "fE1" in line:
                continue

            # Parse data line
            if current_nuclide is not None:
                try:
                    data = reader.read(line)
                    package['U'].append(data[0])
                    package['fE1'].append(data[1])
                except:
                    continue

    # Save the last isotope
    if current_nuclide is not None and package is not None:
        d[current_nuclide] = package

    return d


def write_ascii_file(fpath: str, data: dict, n: _c.Nuclide = None) -> None:
    """Write ASCII file given a filename and data dictionary.

    Args:
        fpath: Path to output file.
        data: Dictionary of entries keyed by Nuclide.
        n: If provided, write only this nucleus. Otherwise write all.
    """
    with open(fpath, 'w', encoding='utf-8') as fp:
        # Determine which entries to write
        if n is not None:
            entries = [(n, data[n])]
        else:
            # Sort by Z then A
            entries = sorted(data.items(), key=lambda x: (x[0].Z, x[0].A))

        for nuclide, entry in entries:
            # Write nucleus header
            fp.write(f" Z=  {nuclide.Z:2d} A=  {nuclide.A:2d} {nuclide.element_symbol}\n")
            # Write column header
            fp.write(FILE_HEADER)
            # Write data
            if hasattr(entry, 'U'):
                # Entry object
                for u, fE1 in zip(entry.U, entry.fE1):
                    fp.write(writer.write((u, fE1)) + "\n")
            else:
                # Dict from reader
                for u, fE1 in zip(entry['U'], entry['fE1']):
                    fp.write(writer.write((u, fE1)) + "\n")


class Entry(_core.GammaStrengthFunction, _db.NuclideDatabaseEntry):
    """An object which represents an entry in the GSF database. """
    pass


class PacketEntry(_db.PacketEntry):
    """A microscopic gamma-ray strength function (GSF) packet entry."""

    _field_info = {
        'U': 'Photon-energy grid [MeV]',
        'fE1': 'E1 strength column [mb/MeV for legacy gamma-strength-micro files; MeV^-3 for the D1M+QRPA fallback]',
    }


def _write_for_db(fpath: str, data: dict) -> None:
    """Wrapper for write_ascii_file that matches the base class save() signature."""
    write_ascii_file(fpath, data)


class Database(_db.Database):

    reader: object = read_ascii_file
    entry : object = Entry
    writer: object = _write_for_db

    def load(self, fpath: str) -> None:
        """Read gamma-ray strength (GSF) information from an ASCII file for a single entry. 
           Not this can be used in successive calls, each time adding potentially new entries to the database (self.data).
        """
        _data = type(self).reader(fpath)
        for entry in _data.keys():
            self.data[entry] = PacketEntry(_data[entry])

    def load_all(self, directory: str) -> None:
        """Load the entire gamma-ray strength function (GSF) database.

        The legacy ``gamma/gamma-strength-micro`` directory is NOT shipped in
        the RIPL-4 github release. When that directory is absent we fall back
        to the new ``gamma/d1m/`` per-Z D1M+QRPA tables (whose ``U``/``fE1``
        schema we replicate) so downstream code that relied on
        ``data[n]['U']`` and ``data[n]['fE1']`` keeps working.
        """
        db_loc = _os.path.join(directory, local_data_dir)
        if _os.path.isdir(db_loc):
            for fn in _os.listdir(db_loc):
                fpath = _os.path.join(db_loc, fn)
                self.load(fpath)
            return
        # Fall back to the D1M+QRPA tables
        _logger.warning(
            f"GSF legacy directory not found ({db_loc}); falling back to "
            "D1M+QRPA tables under 'gamma/d1m/'."
        )
        from . import d1m as _d1m
        d1m_dir = _os.path.join(directory, _d1m.local_data_dir)
        if not _os.path.isdir(d1m_dir):
            _logger.warning(f"D1M fallback directory not found: {d1m_dir}")
            return
        for fn in _os.listdir(d1m_dir):
            fpath = _os.path.join(d1m_dir, fn)
            if not _os.path.isfile(fpath):
                continue
            try:
                pkts = _d1m.read_ascii_file(fpath)
            except Exception as exc:  # noqa: BLE001
                _logger.debug(f"GSF fallback: cannot parse {fpath}: {exc}")
                continue
            for n, pkt in pkts.items():
                # Keep only the schema the legacy GSF database exposed
                self.data[n] = PacketEntry({'U': pkt['U'], 'fE1': pkt['fE1']})

    def save(self, fpath: str, entry: _c.Nuclide = None) -> None:
        """Save entries to a file.

        Args:
            fpath: Path to output file.
            entry: If provided, save only this nucleus. Otherwise save all.
        """
        write_ascii_file(fpath, self.data, entry)

    def save_all(self, directory: str) -> None:
        """Save all entries to per-element files in a directory. """
        _os.makedirs(directory, exist_ok=True)
        # Group by Z
        elements = set(n.Z for n in self.data.keys())
        for Z in sorted(elements):
            element_data = {k: v for k, v in self.data.items() if k.Z == Z}
            if element_data:
                fpath = _os.path.join(directory, f"z{Z:03d}.dat")
                write_ascii_file(fpath, element_data)

    def save_entry(self, fpath: str, entry: _c.Nuclide) -> None:
        """Save a single entry to a file. """
        write_ascii_file(fpath, self.data, entry)

    def save_entries_by_Z(self, directory: str) -> None:
        """Save entries grouped by element to separate files. """
        for entry in self.data.keys():
            fn = f"gsf_Z{entry.Z:03d}A{entry.A:03d}.dat"
            fpath = _os.path.join(directory, fn)
            self.save(fpath, entry)


def load_all(directory: str = None) -> Database:
    """Load and return the GSF database. """
    directory = _resolve_directory(directory)
    d = Database()
    d.load_all(directory)
    return d


def load_nucleus(n: _c.Nuclide, directory: str = None) -> Database:
    """Load the GSF data for a single nucleus.

    Uses the legacy ``gamma/gamma-strength-micro`` per-nucleus file when
    present; otherwise falls back to the RIPL-4 D1M+QRPA per-Z table (the
    legacy directory is not shipped in the GitHub release), mirroring
    :meth:`Database.load_all`.
    """
    directory = _resolve_directory(directory)
    d = Database()
    fpath = _os.path.join(directory, local_data_dir, f'gsf_Z{n.Z:03d}A{n.A:03d}.dat')
    if _os.path.exists(fpath):
        d.load(fpath)
        return d

    # Fall back to the D1M+QRPA per-Z table and extract this nuclide.
    from . import d1m as _d1m
    d1m_file = _os.path.join(directory, _d1m.local_data_dir, f"z{n.Z:03d}_e1")
    if not _os.path.exists(d1m_file):
        _logger.warning(
            f"GSF data for {n} not available (no legacy file and no D1M "
            f"fallback at {d1m_file})."
        )
        return d
    pkts = _d1m.read_ascii_file(d1m_file)
    if n in pkts:
        d.data[n] = PacketEntry({'U': pkts[n]['U'], 'fE1': pkts[n]['fE1']})
    return d