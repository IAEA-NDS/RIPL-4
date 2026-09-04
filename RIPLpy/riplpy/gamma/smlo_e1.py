# -*- coding: utf-8 -*-
"""SMLO E1 photon strength function tables (RIPL-4).

The github release ships per-nucleus E1 photoabsorption PSF tables under
``gamma/smlo_E1/`` (~8980 files) following the naming pattern
``fe1_the_<ZZZ>_<AAA>_photoabs_h_SMLO.dat``.

Each file uses the column layout::

    # Z = 8, A = 16
    # Parameters: Er1, Wr1, S1, beta
    # Col 1: photon energy E in MeV
    # Col 2-12: E1 PSF (fe1 in MeV^-3) at T=0.0..2.0 MeV in 0.2 MeV steps
    #   E         T=0.0       T=0.2  ...  T=2.0
    0.100  ...

For each nucleus we store ``U`` (the photon energy grid), ``fE1`` (cold-nucleus
column T=0.0), the temperature labels ``T``, the full 2-D array ``fE1_T``, and
the Lorentzian parameters ``Er1, Wr1, S1, beta`` extracted from the header.

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

__all__ = ('local_data_dir', 'read_ascii_file', 'Database', 'load',
           'load_all', 'load_nucleus')

local_data_dir = _os.path.join('gamma', 'smlo_E1')

_FNAME_RE = _re.compile(r"fe1_the_(\d+)_(\d+)_photoabs", _re.IGNORECASE)


def _z_a_from_filename(fn: str) -> tuple[int, int] | None:
    m = _FNAME_RE.search(fn)
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def _parse_zheader(line: str) -> tuple[int, int] | None:
    """Parse the ``# Z = <z>, A = <a>`` header line."""
    s = line.replace('#', '').replace(',', ' ').replace('=', ' = ')
    tokens = s.split()
    Z = A = None
    i = 0
    while i < len(tokens):
        if tokens[i] == 'Z' and i + 2 < len(tokens):
            try:
                Z = int(tokens[i + 2])
            except ValueError:
                pass
        elif tokens[i] == 'A' and i + 2 < len(tokens):
            try:
                A = int(tokens[i + 2])
            except ValueError:
                pass
        i += 1
    if Z is None or A is None:
        return None
    return Z, A


def _parse_parameters(line: str) -> dict:
    """Parse the ``# Parameters: Er1 = ..., Wr1 = ..., S1 = ..., beta = ...``."""
    out: dict = {}
    chunk = line.replace('#', '').replace(',', ' ')
    tokens = chunk.split()
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok.endswith('=') and i + 1 < len(tokens):
            key = tok[:-1].strip()
            try:
                out[key] = float(tokens[i + 1])
            except ValueError:
                pass
            i += 2
        elif tok in ('Er1', 'Wr1', 'S1', 'beta') and i + 2 < len(tokens) and tokens[i + 1] == '=':
            try:
                out[tok] = float(tokens[i + 2])
            except ValueError:
                pass
            i += 3
        else:
            i += 1
    return out


def _parse_temperature_labels(line: str) -> list[str]:
    return [tok for tok in line.split() if tok.startswith('T=')]


def read_ascii_file(fpath: str) -> dict:
    """Read a single SMLO E1 per-nucleus file."""
    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        lines = fp.readlines()

    Z = A = None
    params: dict = {}
    labels: list[str] = []
    U = []
    fE1_T = []
    for raw in lines:
        line = raw.rstrip('\n')
        if not line.strip():
            continue
        stripped = line.lstrip()
        if stripped.startswith('#'):
            zheader = _parse_zheader(line)
            if zheader is not None and Z is None:
                Z, A = zheader
            if 'Parameters' in line:
                params = _parse_parameters(line)
            if 'T=' in line:
                lbls = _parse_temperature_labels(line)
                if lbls:
                    labels = lbls
            continue
        parts = line.split()
        try:
            energy = float(parts[0])
        except (ValueError, IndexError):
            continue
        try:
            row = [float(x) for x in parts[1:]]
        except ValueError:
            continue
        U.append(energy)
        fE1_T.append(row)

    # Fallback to filename if header parsing failed
    if Z is None or A is None:
        za = _z_a_from_filename(_os.path.basename(fpath))
        if za is not None:
            Z, A = za
    if Z is None or A is None:
        raise ValueError(f"Could not infer Z, A from {fpath!r}")

    n = _c.Nuclide(Z=Z, A=A)
    fE1 = [row[0] if row else None for row in fE1_T]
    pkt = {
        'n': n, 'U': U, 'fE1': fE1, 'fE1_T': fE1_T, 'T': labels,
    }
    pkt.update({k: params[k] for k in ('Er1', 'Wr1', 'S1', 'beta') if k in params})
    return {n: pkt}


class Entry(_db.PacketEntry):
    """A SMLO E1 photoabsorption photon strength function entry."""

    _field_info = {
        'n': 'Target nucleus',
        'U': 'Photon-energy grid [MeV]',
        'fE1': 'Cold-nucleus (T=0) E1 photon strength function [MeV^-3]',
        'fE1_T': 'Per-temperature E1 strength rows (rows = U, cols = T) [MeV^-3]',
        'T': 'Temperature column labels (T=0.0..2.0 MeV)',
        'Er1': 'SMLO Lorentzian peak energy [MeV]',
        'Wr1': 'SMLO Lorentzian width [MeV]',
        'S1': 'SMLO Lorentzian strength [TRK units]',
        'beta': 'Quadrupole deformation parameter',
    }


class Database(_db.NuclideDatabase):

    reader: object = read_ascii_file
    entry : object = Entry
    local_data_dir: str = local_data_dir

    def load(self, fpath: str) -> None:
        data = type(self).reader(fpath)
        entry_cls = type(self).entry
        for n, pkt in data.items():
            self.data[n] = entry_cls(pkt)

    def load_all(self, directory: str) -> None:
        db_loc = _os.path.join(directory, self.local_data_dir)
        if not _os.path.isdir(db_loc):
            _logger.warning(f"SMLO E1 data directory not found: {db_loc}")
            return
        for fn in sorted(_os.listdir(db_loc)):
            fpath = _os.path.join(db_loc, fn)
            if _os.path.isfile(fpath):
                try:
                    self.load(fpath)
                except Exception as exc:  # noqa: BLE001
                    _logger.debug(f"Failed to parse SMLO E1 file {fpath}: {exc}")


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


def load_nucleus(n: _c.Nuclide, directory: str = None) -> Database:
    directory = _resolve_directory(directory)
    db = Database()
    fname = f"fe1_the_{n.Z:03d}_{n.A:03d}_photoabs_h_SMLO.dat"
    fpath = _os.path.join(directory, local_data_dir, fname)
    if _os.path.exists(fpath):
        db.load(fpath)
    else:
        _logger.warning(f"SMLO E1 file not found: {fpath}")
    return db
