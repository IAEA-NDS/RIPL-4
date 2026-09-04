# -*- coding: utf-8 -*-
"""Reader for the ATOMKI TALYS alpha-OMP gnu output files.

The ``optical/atomki/`` directory contains ~4359 plot-ready output files
generated with the TALYS alpha optical model potential (``alphaomp9real``)
for a wide range of target nuclei. Each file is named::

    z{Z:03d}a{A:03d}a_talys_alphaomp9real.gnu

and contains a two-line text header summarising the potential parameters
(``AP``, ``AT``, ``ZP``, ``ZT``, ``JR0``, ``r_rms``, ``lambda``, ``JR``),
followed by a single floating-point value (the energy-independent
normalization), and a multi-column tabular section listing the radial
profile of the real part of the alpha-nucleus optical potential.

Because the directory contains thousands of files, this module **does not**
preload anything at section load time. Use :func:`load_nucleus` to read a
single (Z, A) on demand.

Example:
    >>> from riplpy.optical import atomki
    >>> entry = atomki.load_nucleus(Z=26, A=42)
    >>> entry.Z, entry.A
    (26, 42)
    >>> entry.r[:3], entry.U[:3]
    ([0.025, 0.05, 0.075], [-179.21216, -179.20035, -179.18067])
"""

# Standard library
import os as _os
from dataclasses import dataclass as _dataclass, field as _field
from typing import List as _List, Optional as _Optional


# =============================================================================
# Dataclasses
# =============================================================================

@_dataclass
class AtomkiAlphaPotential:
    """Tabulated real alpha-OMP profile for a single (Z, A) target nucleus.

    Attributes:
        Z: Target atomic number.
        A: Target mass number.
        AP: Projectile mass number (typically 4 for alpha).
        ZP: Projectile atomic number (typically 2 for alpha).
        magic: True if the target is flagged as a magic nucleus.
        JR0: Volume integral JR0 in MeV*fm^3.
        r_rms: Root-mean-square matter radius in fm.
        lambda_: Normalisation constant (lambda) parsed from line 2.
        JR: Reference volume integral JR in MeV*fm^3.
        norm: The single floating-point normalisation on line 3.
        r: Radial mesh in fm.
        U_real: First potential column (real central potential, MeV).
        U_aux1: Second potential column.
        U_aux2: Third potential column (typically U_real + U_aux1).
        path: Source filename on disk.
    """

    Z: int
    A: int
    AP: int = 0
    ZP: int = 0
    magic: bool = False
    JR0: float = 0.0
    r_rms: float = 0.0
    lambda_: float = 0.0
    JR: float = 0.0
    norm: float = 0.0
    r: _List[float] = _field(default_factory=list)
    U_real: _List[float] = _field(default_factory=list)
    U_aux1: _List[float] = _field(default_factory=list)
    U_aux2: _List[float] = _field(default_factory=list)
    path: str = ''

    @property
    def n_points(self) -> int:
        """Number of radial mesh points."""
        return len(self.r)


# =============================================================================
# Filename helpers
# =============================================================================

def filename_for(Z: int, A: int) -> str:
    """Return the canonical ATOMKI filename for a given (Z, A).

    Args:
        Z: Target atomic number.
        A: Target mass number.

    Returns:
        Filename matching ``z{Z:03d}a{A:03d}a_talys_alphaomp9real.gnu``.
    """
    return 'z{Z:03d}a{A:03d}a_talys_alphaomp9real.gnu'.format(Z=Z, A=A)


def _default_directory(directory: _Optional[str] = None) -> str:
    """Resolve the ATOMKI directory under the configured RIPL path."""
    if directory is not None:
        return directory
    from riplpy import config as _ripl_config
    return _os.path.join(_ripl_config.get_path(), 'optical', 'atomki')


# =============================================================================
# Header parsing
# =============================================================================

def _parse_header(line1: str, line2: str) -> dict:
    """Parse the two-line ATOMKI header into a metadata dict.

    The header has the form::

        Potential for AP =   4 and AT =  42 (non-magic), VC for ZP =   2, ZT =  26
        JR0 =  299.102 MeV*fm^3, r_rms =  4.7774 fm, lambda =  1.2404 for JR =  371.0

    Both lines are robustly parsed by splitting on ``=`` and pulling the
    leading float token from each field.
    """
    meta: dict = {}

    def _grab(text: str, key: str) -> _Optional[float]:
        """Pull a number that follows ``key =`` in ``text``."""
        marker = key + ' ='
        idx = text.find(marker)
        if idx < 0:
            marker = key + '='
            idx = text.find(marker)
            if idx < 0:
                return None
        rest = text[idx + len(marker):].lstrip()
        token = ''
        for ch in rest:
            if ch in '0123456789+-.eE':
                token += ch
            else:
                break
        try:
            return float(token)
        except ValueError:
            return None

    ap = _grab(line1, 'AP')
    at = _grab(line1, 'AT')
    zp = _grab(line1, 'ZP')
    # ZT is the trailing token on line1 (no following text), so handle the
    # trailing newline gracefully.
    zt = _grab(line1, 'ZT')

    meta['AP'] = int(ap) if ap is not None else 0
    meta['AT'] = int(at) if at is not None else 0
    meta['ZP'] = int(zp) if zp is not None else 0
    meta['ZT'] = int(zt) if zt is not None else 0
    meta['magic'] = ('non-magic' not in line1) and ('magic' in line1)

    meta['JR0'] = _grab(line2, 'JR0') or 0.0
    meta['r_rms'] = _grab(line2, 'r_rms') or 0.0
    meta['lambda'] = _grab(line2, 'lambda') or 0.0
    meta['JR'] = _grab(line2, 'JR') or 0.0
    return meta


# =============================================================================
# Public API
# =============================================================================

def load_nucleus(Z: int, A: int, directory: _Optional[str] = None) -> AtomkiAlphaPotential:
    """Read the ATOMKI alpha-OMP profile for the target ``(Z, A)``.

    Args:
        Z: Target atomic number.
        A: Target mass number.
        directory: Optional override for the ``optical/atomki`` directory.
            When omitted, the configured RIPL path is used.

    Returns:
        :class:`AtomkiAlphaPotential` with header metadata and the radial
        mesh + potential columns.

    Raises:
        FileNotFoundError: If no file exists for the requested (Z, A).
    """
    atomki_dir = _default_directory(directory)
    path = _os.path.join(atomki_dir, filename_for(Z, A))
    if not _os.path.isfile(path):
        raise FileNotFoundError(
            "ATOMKI alpha-OMP file not found for Z={Z}, A={A}: {p}".format(
                Z=Z, A=A, p=path))

    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()

    if len(lines) < 4:
        raise ValueError("ATOMKI file too short: {p}".format(p=path))

    meta = _parse_header(lines[0], lines[1])
    try:
        norm = float(lines[2].strip())
    except ValueError:
        norm = 0.0

    r: _List[float] = []
    U_real: _List[float] = []
    U_aux1: _List[float] = []
    U_aux2: _List[float] = []
    for raw in lines[3:]:
        toks = raw.split()
        if len(toks) < 4:
            continue
        try:
            r.append(float(toks[0]))
            U_real.append(float(toks[1]))
            U_aux1.append(float(toks[2]))
            U_aux2.append(float(toks[3]))
        except ValueError:
            continue

    return AtomkiAlphaPotential(
        Z=Z, A=A,
        AP=meta['AP'], ZP=meta['ZP'], magic=meta['magic'],
        JR0=meta['JR0'], r_rms=meta['r_rms'],
        lambda_=meta['lambda'], JR=meta['JR'],
        norm=norm,
        r=r, U_real=U_real, U_aux1=U_aux1, U_aux2=U_aux2,
        path=path,
    )


def available_nuclei(directory: _Optional[str] = None) -> _List[tuple]:
    """List the (Z, A) pairs available in the ATOMKI directory.

    Args:
        directory: Optional override for the ATOMKI directory.

    Returns:
        Sorted list of ``(Z, A)`` tuples discovered from filenames.
    """
    atomki_dir = _default_directory(directory)
    if not _os.path.isdir(atomki_dir):
        return []

    out: _List[tuple] = []
    for name in _os.listdir(atomki_dir):
        if not name.endswith('_talys_alphaomp9real.gnu'):
            continue
        if not (name.startswith('z') and len(name) >= 9 and name[4] == 'a'):
            continue
        try:
            z = int(name[1:4])
            a = int(name[5:8])
        except ValueError:
            continue
        out.append((z, a))
    out.sort()
    return out
