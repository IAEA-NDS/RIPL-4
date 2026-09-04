# -*- coding: utf-8 -*-
"""Python objects which provide access to the HFB (Hartree-Fock-Bogoliubov) level densities
from the Reference Input Parameter Library (RIPL).

RIPL-3 LEGACY: The HFB spin-dependent level density tables
(``densities/total/level-densities-hfb/``) are a RIPL-3 legacy product. They
ship with the full RIPL distribution but are NOT part of the RIPL-4 GitHub
release, which provides the microscopic combinatorial tables instead
(``bsk14-comb``, ``bskg3-comb``, ``qrpabe``, ``thfb-comb``).
``riplpy.densities.load()`` skips this model with a warning when the directory
is absent.

The HFB level densities provide spin-dependent nuclear level densities for both
positive and negative parities. The data is stored in per-element files (zXXX.tab).

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Load HFB level densities for iron:
        $ import riplpy.densities.hfb as hfb
        $ db = hfb.load_element(Z=26, directory='/path/to/RIPL-4')
        $ print(db.data)

"""

# ========================

# OS
import os as _os

# Dataclasses
from dataclasses import dataclass as _dataclass
from dataclasses import field as _field

# Typing
from typing import ClassVar as _ClassVar

# Logging
import logging as _logging

# RIPLpy
import riplpy.collections as _c
import riplpy.db as _db

# Module logger
_logger = _logging.getLogger(__name__)

# ========================

__all__ = ('local_data_dir', 'MAX_SPIN', 'parse_header', 'read_ascii_file', 'ParityData', 'Entry', 'Database', 'load', 'load_all', 'load_element')

# ========================

# Relative path to database folder
local_data_dir = _os.path.join("densities", "total", "level-densities-hfb")

# Maximum spin value in the data files
MAX_SPIN = 50

# ========================


def parse_header(line: str) -> tuple:
    """Parse the header line to extract Z, A, and parity information.

    Header format: '*  Z= 26 A= 42: Positive-Parity Spin-dependent Level Density...'
    """
    # Extract Z value
    z_start = line.find('Z=') + 2
    z_end = line.find('A=') - 1
    Z = int(line[z_start:z_end].strip())

    # Extract A value
    a_start = line.find('A=') + 2
    a_end = line.find(':')
    A = int(line[a_start:a_end].strip())

    # Determine parity
    parity = 1 if 'Positive' in line else -1

    return Z, A, parity


def read_ascii_file(fpath: str) -> dict:
    """Read HFB level density data from an ASCII file and return a dictionary.

    Each file contains data for multiple isotopes of a given element, with both
    positive and negative parity blocks for each isotope.
    """
    d = {}

    current_nuclide = None
    current_parity = None
    current_data = None
    in_data_section = False

    with open(fpath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            line_stripped = line.strip()

            # Skip empty lines
            if not line_stripped:
                continue

            # Check for header block start/end
            if line_stripped.startswith('*') and 'Z=' in line_stripped and 'A=' in line_stripped:
                # Parse the header
                Z, A, parity = parse_header(line_stripped)
                current_nuclide = _c.Nuclide(Z=Z, A=A)
                current_parity = parity

                # Initialize nuclide entry if not exists
                if current_nuclide not in d:
                    d[current_nuclide] = {
                        'n': current_nuclide,
                        'positive_parity': None,
                        'negative_parity': None
                    }

                # Initialize parity data structure
                current_data = {
                    'U': [],       # Excitation energy [MeV]
                    'T': [],       # Temperature [MeV]
                    'Ncumul': [],  # Cumulative level count
                    'rho_obs': [], # Observed level density
                    'rho_tot': [], # Total level density
                    'rho_J': []    # Spin-dependent densities (list of lists)
                }
                in_data_section = False
                continue

            # Skip asterisk-only lines and column headers
            if line_stripped.startswith('*'):
                continue
            if 'U[MeV]' in line_stripped or 'RHOTOT' in line_stripped:
                in_data_section = True
                continue

            # Parse data lines
            if in_data_section and current_data is not None:
                try:
                    values = line.split()
                    if len(values) >= 5:
                        current_data['U'].append(float(values[0]))
                        current_data['T'].append(float(values[1]))
                        current_data['Ncumul'].append(float(values[2]))
                        current_data['rho_obs'].append(float(values[3]))
                        current_data['rho_tot'].append(float(values[4]))

                        # Spin-dependent densities (remaining columns)
                        spin_densities = [float(v) for v in values[5:]]
                        current_data['rho_J'].append(spin_densities)
                except (ValueError, IndexError) as e:
                    _logger.debug(f"Skipping malformed line: {line.strip()}")
                    continue

            # Check if we've reached the end of a parity block
            # This happens when we encounter a new header or end of file
            if current_data is not None and len(current_data['U']) > 0:
                if current_parity == 1:
                    d[current_nuclide]['positive_parity'] = current_data.copy()
                else:
                    d[current_nuclide]['negative_parity'] = current_data.copy()

    # Finalize the last block
    if current_data is not None and len(current_data['U']) > 0 and current_nuclide is not None:
        if current_parity == 1:
            d[current_nuclide]['positive_parity'] = current_data
        else:
            d[current_nuclide]['negative_parity'] = current_data

    return d


@_dataclass
class ParityData:
    """Level density data for a single parity.

    Contains tabulated level densities as a function of excitation energy,
    including spin-dependent values.

    Attributes:
        U: Excitation energy grid [MeV]
        T: Nuclear temperature [MeV]
        Ncumul: Cumulative level count
        rho_obs: Observed level density [MeV^-1]
        rho_tot: Total level density [MeV^-1]
        rho_J: Spin-dependent level densities [MeV^-1]
    """

    U: list = _field(default_factory=list)
    T: list = _field(default_factory=list)
    Ncumul: list = _field(default_factory=list)
    rho_obs: list = _field(default_factory=list)
    rho_tot: list = _field(default_factory=list)
    rho_J: list = _field(default_factory=list)

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'U': 'Excitation energy grid [MeV]',
        'T': 'Nuclear temperature [MeV]',
        'Ncumul': 'Cumulative level count',
        'rho_obs': 'Observed level density [MeV^-1]',
        'rho_tot': 'Total level density [MeV^-1]',
        'rho_J': 'Spin-dependent level densities [MeV^-1]',
    }

    @classmethod
    def from_dict(cls, data: dict) -> "ParityData":
        """Create a ParityData object from a dictionary."""
        if data is None:
            return None
        return cls(
            U=data.get('U', []),
            T=data.get('T', []),
            Ncumul=data.get('Ncumul', []),
            rho_obs=data.get('rho_obs', []),
            rho_tot=data.get('rho_tot', []),
            rho_J=data.get('rho_J', [])
        )

    def get_density_at_energy(self, energy: float) -> dict:
        """Get level density data at a specific excitation energy."""
        if not self.U:
            return None

        # Find the closest energy point
        closest_idx = min(range(len(self.U)), key=lambda i: abs(self.U[i] - energy))

        return {
            'U': self.U[closest_idx],
            'T': self.T[closest_idx],
            'Ncumul': self.Ncumul[closest_idx],
            'rho_obs': self.rho_obs[closest_idx],
            'rho_tot': self.rho_tot[closest_idx],
            'rho_J': self.rho_J[closest_idx] if self.rho_J else []
        }


@_dataclass
class Entry(_db.NuclideDatabaseEntry):
    """An entry in the HFB level density database.

    Contains microscopic HFB-based level densities for both parities,
    tabulated as a function of excitation energy with spin dependence.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        positive_parity: Level density data for positive parity states
        negative_parity: Level density data for negative parity states
    """

    n: _c.Nuclide = None
    positive_parity: ParityData = None
    negative_parity: ParityData = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'positive_parity': 'Level density data for positive parity states',
        'negative_parity': 'Level density data for negative parity states',
    }

    def get_total_density(self, energy: float, parity: int = None) -> float:
        """Get total level density at a given excitation energy.

        Args:
            energy: Excitation energy in MeV
            parity: +1 for positive, -1 for negative, None for sum of both

        Returns:
            Total level density in MeV^-1
        """
        result = 0.0

        if parity is None or parity == 1:
            if self.positive_parity is not None:
                data = self.positive_parity.get_density_at_energy(energy)
                if data:
                    result += data['rho_tot']

        if parity is None or parity == -1:
            if self.negative_parity is not None:
                data = self.negative_parity.get_density_at_energy(energy)
                if data:
                    result += data['rho_tot']

        return result


def write_ascii_file(fpath: str, data: dict) -> None:
    """Write HFB level density data to an ASCII file.

    Args:
        fpath: Path to the output file.
        data: Dictionary of Entry objects keyed by Nuclide.
    """
    with open(fpath, 'w', encoding='utf-8') as fp:
        # Sort entries by Z then A
        sorted_entries = sorted(data.values(), key=lambda e: (e.n.Z, e.n.A))

        for entry in sorted_entries:
            Z, A = entry.n.Z, entry.n.A
            symbol = entry.n.element_symbol

            # Write positive parity block if present
            if entry.positive_parity is not None:
                _write_parity_block(fp, Z, A, symbol, entry.positive_parity, parity=1)

            # Write negative parity block if present
            if entry.negative_parity is not None:
                _write_parity_block(fp, Z, A, symbol, entry.negative_parity, parity=-1)


def _write_parity_block(fp, Z: int, A: int, symbol: str, parity_data: "ParityData", parity: int) -> None:
    """Write a single parity block to file."""
    parity_str = "Positive" if parity == 1 else "Negative"

    # Write header
    header = f"*  Z= {Z:2d} A= {A:3d}: {parity_str}-Parity Spin-dependent Level Density [MeV-1] for {symbol:2s} {A:3d}  *"
    border = "*" * (len(header))
    fp.write(f"                    {border}\n")
    fp.write(f"                    {header}\n")
    fp.write(f"                    {border}\n")

    # Write column header
    header_cols = " U[MeV]  T[MeV]  NCUMUL   RHOOBS   RHOTOT"
    for j in range(MAX_SPIN):
        header_cols += f"     J={j:02d}"
    fp.write(header_cols + "\n")

    # Write data lines
    for i in range(len(parity_data.U)):
        line = f"{parity_data.U[i]:7.2f}{parity_data.T[i]:7.3f}"
        line += f"{parity_data.Ncumul[i]:10.2E}{parity_data.rho_obs[i]:10.2E}{parity_data.rho_tot[i]:10.2E}"

        # Spin-dependent densities
        if parity_data.rho_J and i < len(parity_data.rho_J):
            for rho in parity_data.rho_J[i]:
                line += f"{rho:10.2E}"
            # Pad with zeros if needed
            for _ in range(MAX_SPIN - len(parity_data.rho_J[i])):
                line += f"{0.0:10.2E}"
        else:
            for _ in range(MAX_SPIN):
                line += f"{0.0:10.2E}"

        fp.write(line + "\n")


class Database(_db.NuclideDatabase):
    """The HFB level density database."""

    reader: object = read_ascii_file
    entry: object = Entry
    writer: object = write_ascii_file

    def load(self, fpath: str) -> None:
        """Load HFB level density data from an ASCII file."""
        _data = type(self).reader(fpath)
        for n, data in _data.items():
            positive = ParityData.from_dict(data.get('positive_parity'))
            negative = ParityData.from_dict(data.get('negative_parity'))
            self.data[n] = Entry(n=n, positive_parity=positive, negative_parity=negative)

    def load_all(self, directory: str) -> None:
        """Load HFB level densities from all element files."""
        db_loc = _os.path.join(directory, local_data_dir)
        if not _os.path.exists(db_loc):
            _logger.warning(f"HFB level density directory not found: {db_loc}")
            return

        for fn in _os.listdir(db_loc):
            if fn.endswith('.tab'):
                fpath = _os.path.join(db_loc, fn)
                self.load(fpath)

    def load_element(self, Z: int, directory: str) -> None:
        """Load HFB level densities for a specific element."""
        fpath = _os.path.join(directory, local_data_dir, f"z{Z:03d}.tab")
        if _os.path.exists(fpath):
            self.load(fpath)
        else:
            _logger.warning(f"HFB level density file not found: {fpath}")

    def save_element(self, Z: int, directory: str) -> None:
        """Save HFB level densities for a specific element to a file."""
        fpath = _os.path.join(directory, f"z{Z:03d}.tab")
        # Filter entries for this element
        element_data = {k: v for k, v in self.data.items() if v.n.Z == Z}
        if element_data:
            type(self).writer(fpath, element_data)

    def save_all(self, directory: str) -> None:
        """Save HFB level densities to per-element files."""
        _os.makedirs(directory, exist_ok=True)
        # Group entries by Z
        elements = set(entry.n.Z for entry in self.data.values())
        for Z in sorted(elements):
            self.save_element(Z, directory)


def load(directory: str = None, file_path: str = None, Z: int = None) -> Database:
    """Load and return HFB level density data.

    Args:
        directory: Path to the RIPL directory. If provided without Z, loads all elements.
        file_path: Direct path to a specific zXXX.tab file.
        Z: Element number. If provided with directory, loads only that element.

    Returns:
        Database: An instance of the Database class with the loaded data.
    """
    db = Database()
    if file_path is not None:
        db.load(file_path)
    elif directory is not None and Z is not None:
        db.load_element(Z, directory)
    elif directory is not None:
        db.load_all(directory)
    return db


def load_all(directory: str) -> Database:
    """Load and return the entire HFB level density database."""
    return load(directory=directory)


def load_element(Z: int, directory: str) -> Database:
    """Load and return HFB level densities for a specific element."""
    return load(directory=directory, Z=Z)
