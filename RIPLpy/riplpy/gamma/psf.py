# -*- coding: utf-8 -*-
"""Experimental Photon Strength Function (PSF) database (RIPL-4).

The PSFDatabase-v2024.1 under ``gamma/PSFDatabase-v2024.1/`` collects
experimentally measured photon strength functions from a variety of techniques:

* ``arcdrc/``       -- Average Resonance Capture / Direct Capture (DRC)
* ``nrf/``          -- Nuclear Resonance Fluorescence
* ``oslo/``         -- Oslo method
* ``pg/``           -- Primary Gamma decay
* ``photonuclear/`` -- Photonuclear (e.g. (gamma,X)) measurements
* ``pp/``           -- Polarised Photon / Photon Scattering
* ``RM/``           -- Resonance Method
* ``thc/``          -- Two-step gamma cascade / thermal capture

Each ``.dat`` file is one experimental dataset for one nucleus, in a uniform
ASCII layout where leading ``#`` lines hold the metadata (``Z``, ``A``,
authors, format, ``Col n`` labels) and the remaining lines hold a small table.
The first column is always the photon energy ``E`` [MeV], and the strength
column (typically ``f`` or ``f1`` or ``fE1`` in ``MeV^-3``) is identified by
the ``Col`` annotations. We always store at least ``E`` and ``f``; full per-row
data is exposed via ``rows``.

This is a *partial* implementation: it loads the common ``# Col`` header layout
which covers the vast majority of the database. Files that deviate from that
convention are skipped with a debug log entry.

"""

# OS
import os as _os

# Logging
import logging as _logging

# Re
import re as _re

# RIPLpy
import riplpy.collections as _c
import riplpy.db as _db
from riplpy.config import resolve_directory as _resolve_directory

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = ('local_data_dir', 'CATEGORIES', 'read_ascii_file', 'Database',
           'load', 'load_all', 'load_category')

local_data_dir = _os.path.join('gamma', 'PSFDatabase-v2024.1')

# Recognised subdirectories
CATEGORIES = ('arcdrc', 'nrf', 'oslo', 'pg', 'photonuclear', 'pp', 'RM', 'thc')

_FNAME_RE = _re.compile(r"(\d{2,3})[_ ]+(\d{2,3})")


def _z_a_from_filename(fn: str) -> tuple[int, int] | None:
    """Best-effort (Z, A) extraction from filenames like ``f1_exp_021_043_OM.dat``."""
    m = _FNAME_RE.search(fn)
    if not m:
        return None
    try:
        return int(m.group(1)), int(m.group(2))
    except ValueError:
        return None


def _parse_header(lines: list[str]) -> tuple[int | None, int | None, list[str]]:
    """Extract Z, A, and column labels from the ``#``-prefixed header block."""
    Z = A = None
    col_labels: list[str] = []
    for raw in lines:
        line = raw.lstrip()
        if not line.startswith('#'):
            break
        body = line.lstrip('#').strip()
        # Match `Z = NN, A = NN`
        if 'Z' in body and 'A' in body and '=' in body:
            cleaned = body.replace(',', ' ').replace('=', ' = ')
            tokens = cleaned.split()
            i = 0
            while i < len(tokens):
                if tokens[i] == 'Z' and i + 2 < len(tokens) and tokens[i + 1] == '=':
                    try:
                        Z = int(tokens[i + 2])
                    except ValueError:
                        pass
                    i += 3
                    continue
                if tokens[i] == 'A' and i + 2 < len(tokens) and tokens[i + 1] == '=':
                    try:
                        A = int(tokens[i + 2])
                    except ValueError:
                        pass
                    i += 3
                    continue
                i += 1
        # Collect Col n: <label> annotations
        if body.lower().startswith('col '):
            parts = body.split(':', 1)
            if len(parts) == 2:
                col_labels.append(parts[1].strip())
    return Z, A, col_labels


def read_ascii_file(fpath: str) -> dict:
    """Read one PSF dataset file. Returns ``{Nuclide: package}``."""
    # PSF files carry free-text author/reference headers that occasionally
    # contain Latin-1 bytes (e.g. a non-breaking space \xa0). Decode
    # permissively so a stray byte cannot drop an entire dataset.
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        raw_lines = fp.readlines()

    Z, A, col_labels = _parse_header(raw_lines)

    # Filename fallback
    if Z is None or A is None:
        za = _z_a_from_filename(_os.path.basename(fpath))
        if za is not None:
            Z, A = za
    if Z is None or A is None:
        raise ValueError(f"PSF: could not infer Z, A from {fpath!r}")

    rows = []
    for raw in raw_lines:
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        try:
            row = [float(x) for x in line.split()]
        except ValueError:
            continue
        if row:
            rows.append(row)

    # Identify the energy (first col) and strength column heuristically.
    E = [r[0] for r in rows] if rows else []
    f = None
    if rows:
        # Try the second column first (most files report `E dE f df` or `E f df`).
        # If the second column looks like a small bin-width (<1 typically), use
        # column index 2 (f); otherwise use index 1.
        idx_strength = 1
        if len(rows[0]) >= 3 and abs(rows[0][1]) < 1.0:
            idx_strength = 2
        f = [r[idx_strength] if idx_strength < len(r) else None for r in rows]

    n = _c.Nuclide(Z=Z, A=A)
    pkt = {
        'n': n, 'E': E, 'f': f, 'rows': rows,
        'columns': col_labels, 'source': _os.path.basename(_os.path.dirname(fpath)),
        'filename': _os.path.basename(fpath),
    }
    return {n: pkt}


class Entry(_db.PacketEntry):
    """An experimental photon strength function (PSF) dataset entry."""

    _field_info = {
        'n': 'Target nucleus',
        'E': 'Photon energy grid [MeV]',
        'f': 'Photon strength function [MeV^-3]',
        'rows': 'Raw numeric data table (one list per data line)',
        'columns': 'Column labels parsed from the file header',
        'source': 'Experimental technique category (subdirectory name)',
        'filename': 'Source data file name',
    }


class Database(_db.NuclideDatabase):
    """PSF experimental database. Keys are Nuclides; values are lists of Entry.

    Multiple datasets per nucleus are stored under ``self.data[n]`` as a list of
    :class:`Entry` packets to preserve them (datasets from different
    experiments).
    """

    reader: object = read_ascii_file
    entry : object = Entry
    local_data_dir: str = local_data_dir

    def load(self, fpath: str) -> None:
        """Merge entries from one PSF file into the database."""
        data = type(self).reader(fpath)
        entry_cls = type(self).entry
        for n, pkt in data.items():
            entry = entry_cls(pkt)
            existing = self.data.get(n)
            if existing is None:
                self.data[n] = [entry]
            elif isinstance(existing, list):
                existing.append(entry)
            else:
                self.data[n] = [existing, entry]

    def load_category(self, category: str, directory: str) -> None:
        db_loc = _os.path.join(directory, self.local_data_dir, category)
        if not _os.path.isdir(db_loc):
            _logger.debug(f"PSF category not found: {db_loc}")
            return
        for fn in sorted(_os.listdir(db_loc)):
            if not fn.endswith('.dat'):
                continue
            fpath = _os.path.join(db_loc, fn)
            if _os.path.isfile(fpath):
                try:
                    self.load(fpath)
                except Exception as exc:  # noqa: BLE001
                    _logger.debug(f"Failed to parse PSF file {fpath}: {exc}")

    def load_all(self, directory: str) -> None:
        db_loc = _os.path.join(directory, self.local_data_dir)
        if not _os.path.isdir(db_loc):
            _logger.warning(f"PSF database directory not found: {db_loc}")
            return
        for category in CATEGORIES:
            self.load_category(category, directory)


def load(directory: str = None, file_path: str = None) -> Database:
    db = Database()
    if file_path is not None:
        db.load(file_path)
    else:
        directory = _resolve_directory(directory)
        db.load_all(directory)
    return db


def load_all(directory: str = None) -> Database:
    return load(directory=directory)


def load_category(category: str, directory: str = None) -> Database:
    directory = _resolve_directory(directory)
    db = Database()
    db.load_category(category, directory)
    return db
