# -*- coding: utf-8 -*-
"""Reader for ``optical/om-data/ROP2013za.dat``.

This file holds the Avrigeanu & Avrigeanu (Dec. 2013, ``iref = 9999``)
revised alpha-particle Real Optical Potential parameter set. For each
target nucleus ``(Z, A)`` it lists pre-computed Woods-Saxon parameters
and the total reaction cross section sigma_R on a fixed energy grid.

Per-nucleus structure (text):

* Header lines containing free-form formulas for ``r_R``, ``a_R``,
  ``V_R``, ``r_WV``, ``a_WV``, ``W_V``, ``r_WD``, ``a_WD``, ``W_D``.
* A boxed tabular block with twelve columns::

      E   V_R   r_R   a_R   W_V   r_WV   a_WV   W_D   r_WD   a_WD   r_C   sigma_R

Only the tabular portion is parsed; the surrounding analytic formulas
are preserved as raw text in :attr:`ROP2013Entry.header_text` for
provenance.

Example:
    >>> from riplpy.optical import rop2013
    >>> db = rop2013.load()
    >>> entry = db.get(21, 45)  # Sc-45
    >>> entry.E[:3], entry.V_R[:3]
    ([2.0, 3.0, 4.0], [164.0, 161.4, 158.8])
"""

# Standard library
import os as _os
from dataclasses import dataclass as _dataclass, field as _field
from typing import Dict as _Dict, List as _List, Optional as _Optional, Tuple as _Tuple


# =============================================================================
# Dataclasses
# =============================================================================

@_dataclass
class ROP2013Entry:
    """Tabulated revised alpha-particle ROP parameters for one (Z, A).

    Columns are stored as parallel lists. All values use units shown in
    the source file (``E`` in MeV, radii in fm, potential depths in MeV,
    ``sigma_R`` in mb).
    """

    Z: int
    A: int
    E: _List[float] = _field(default_factory=list)
    V_R: _List[float] = _field(default_factory=list)
    r_R: _List[float] = _field(default_factory=list)
    a_R: _List[float] = _field(default_factory=list)
    W_V: _List[float] = _field(default_factory=list)
    r_WV: _List[float] = _field(default_factory=list)
    a_WV: _List[float] = _field(default_factory=list)
    W_D: _List[float] = _field(default_factory=list)
    r_WD: _List[float] = _field(default_factory=list)
    a_WD: _List[float] = _field(default_factory=list)
    r_C: _List[float] = _field(default_factory=list)
    sigma_R: _List[float] = _field(default_factory=list)
    header_text: str = ''

    @property
    def n_points(self) -> int:
        """Number of tabulated energy points."""
        return len(self.E)


@_dataclass
class ROP2013Database:
    """Container for all ROP2013 entries keyed by ``(Z, A)``.

    Attributes:
        iref: The RIPL iref this dataset would be filed under (9999).
        authors: Authorship string parsed from the file header.
        reference: Citation string parsed from the file header.
        data: Mapping of ``(Z, A) -> ROP2013Entry``.
    """

    iref: int = 9999
    authors: str = ''
    reference: str = ''
    data: _Dict[_Tuple[int, int], ROP2013Entry] = _field(default_factory=dict)

    def get(self, Z: int, A: int) -> _Optional[ROP2013Entry]:
        """Return the entry for ``(Z, A)`` or ``None`` if missing."""
        return self.data.get((Z, A))

    @property
    def nuclei(self) -> _List[_Tuple[int, int]]:
        """Sorted list of ``(Z, A)`` keys in the database."""
        return sorted(self.data.keys())

    def __len__(self) -> int:
        return len(self.data)


# =============================================================================
# Parsing helpers
# =============================================================================

def _looks_like_data_row(line: str) -> bool:
    """Heuristic: a data row starts with whitespace + an integer/float E
    and has at least 12 whitespace-separated numeric tokens."""
    toks = line.split()
    if len(toks) < 12:
        return False
    try:
        # All twelve leading tokens must parse as floats (the last token
        # uses scientific notation occasionally, e.g. ``0.4795E-08``).
        for t in toks[:12]:
            float(t)
        return True
    except ValueError:
        return False


def _parse_target_header(line: str) -> _Optional[_Tuple[int, int]]:
    """Parse ``Target: Z = 21  A =  45`` -> (21, 45)."""
    stripped = line.strip()
    if not stripped.lower().startswith('target'):
        return None
    # Pull the Z value
    try:
        z_part = stripped.split('Z')[1].split('A')[0]
        z_val = int(z_part.replace('=', '').strip())
        a_part = stripped.split('A')[1]
        a_val = int(a_part.replace('=', '').strip())
        return (z_val, a_val)
    except (IndexError, ValueError):
        return None


# =============================================================================
# Public API
# =============================================================================

def read(path: str) -> ROP2013Database:
    """Read ``ROP2013za.dat`` from ``path`` and return a populated database.

    Args:
        path: Absolute path to the ``ROP2013za.dat`` file.

    Returns:
        :class:`ROP2013Database` keyed by ``(Z, A)``.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
    """
    if not _os.path.isfile(path):
        raise FileNotFoundError(
            "ROP2013za.dat not found: {p}".format(p=path))

    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()

    db = ROP2013Database()

    # Top-of-file metadata (before the first Target block).
    head_buf: _List[str] = []
    current: _Optional[ROP2013Entry] = None
    current_header: _List[str] = []
    seen_first_target = False

    # Each nucleus has a human-readable parameter table followed by a
    # "----- RIPL Library Format -----" block that mimics
    # om-parameter-u.dat. The RIPL-format block also contains a
    # ``Target:`` line we must NOT treat as a new nucleus boundary.
    in_ripl_format = False

    for raw in lines:
        if 'RIPL Library Format' in raw:
            in_ripl_format = True
            continue
        # A horizontal rule before a fresh tabular header marks the end
        # of a RIPL-format block when the next human-readable block
        # starts; we exit the mode upon seeing a new "Target:" line that
        # is not preceded by extra whitespace.

        target = _parse_target_header(raw)
        if target is not None and in_ripl_format:
            # Stay in RIPL-format mode; the matching human-readable
            # Target line for the next nucleus will reset us.
            # Heuristic: the human-readable header uses "Target: Z" with
            # exactly one space after the colon; the RIPL-format header
            # uses "Target:  Z" with two spaces.
            if raw.lstrip().startswith('Target:  '):
                continue
            # Otherwise fall through and treat as a normal header.
            in_ripl_format = False

        if target is not None:
            # Commit any in-progress entry
            if current is not None:
                current.header_text = ''.join(current_header).rstrip()
                db.data[(current.Z, current.A)] = current

            if not seen_first_target:
                # Capture authors/reference from the file preamble
                preamble = ''.join(head_buf)
                for ln in preamble.splitlines():
                    stripped = ln.strip()
                    if stripped.startswith('authors'):
                        db.authors = stripped.split('=', 1)[-1].strip()
                    elif stripped.startswith('reference'):
                        db.reference = stripped.split('=', 1)[-1].strip()
                    elif stripped.startswith('iref'):
                        try:
                            db.iref = int(stripped.split('=', 1)[-1].strip())
                        except ValueError:
                            pass
                seen_first_target = True

            current = ROP2013Entry(Z=target[0], A=target[1])
            current_header = []
            continue

        if not seen_first_target:
            head_buf.append(raw)
            continue

        # Try to parse a 12-column data row
        if _looks_like_data_row(raw):
            toks = raw.split()
            try:
                vals = [float(t) for t in toks[:12]]
            except ValueError:
                current_header.append(raw)
                continue
            if current is None:
                continue
            current.E.append(vals[0])
            current.V_R.append(vals[1])
            current.r_R.append(vals[2])
            current.a_R.append(vals[3])
            current.W_V.append(vals[4])
            current.r_WV.append(vals[5])
            current.a_WV.append(vals[6])
            current.W_D.append(vals[7])
            current.r_WD.append(vals[8])
            current.a_WD.append(vals[9])
            current.r_C.append(vals[10])
            current.sigma_R.append(vals[11])
        else:
            current_header.append(raw)

    # Commit final entry
    if current is not None:
        current.header_text = ''.join(current_header).rstrip()
        db.data[(current.Z, current.A)] = current

    return db


def load(directory: _Optional[str] = None) -> ROP2013Database:
    """Load the ROP2013 database from a RIPL directory.

    Args:
        directory: RIPL directory path; if omitted the configured path is
            used. The reader looks for ``optical/om-data/ROP2013za.dat``.

    Returns:
        Populated :class:`ROP2013Database`.

    Raises:
        FileNotFoundError: If the data file is missing.
    """
    if directory is None:
        from riplpy import config as _ripl_config
        directory = _ripl_config.get_path()
    path = _os.path.join(directory, 'optical', 'om-data', 'ROP2013za.dat')
    return read(path)
