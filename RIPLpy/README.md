# RIPLpy

Python interface for the Reference Input Parameter Library (RIPL).

## Overview

RIPLpy provides access to the RIPL nuclear physics database, containing evaluated nuclear data for applications in nuclear reaction calculations, astrophysics simulations, and nuclear structure studies.

**Supported RIPL Sections:**
- **masses** - Atomic masses and mass models (AME2020, FRDM2012, HFB27, BSkG3, D1M, WS4)
- **densities** - Nuclear level density parameters (EGSM and microscopic combinatorial models)
- **fission** - Fission barrier heights and parameters
- **gamma** - Giant dipole resonance and gamma strength functions
- **levels** - Discrete nuclear level schemes
- **resonances** - Neutron resonance parameters (s-wave and p-wave)
- **optical** - Optical model potentials (~584 parameterizations)

> **RIPL-3 legacy support.** Some databases (mass models FRDM1995 / HFB14;
> level densities BSFG / CT / HFB; Mengoni-Nakajima shell corrections;
> theoretical GDR parameters) are RIPL-3 legacy products. They ship with the
> full RIPL distribution but are **not** part of the RIPL-4 GitHub release.
> RIPLpy retains their readers/writers and skips them with a warning when the
> data files are absent, so code works against both layouts.

## Installation

To install RIPLpy please perform the following actions:
1. Ensure your local environment is connected to the internet.
2. Download the code from the Git repository.
3. If necessary, unzip the compressed file and switch to the top level directory (which should include this read me file).
4. Perform the appropriate pip installation based off how you want RIPL to behave in your given environment.



```bash
# Development install (recommended)
pip install -e .
```

The above command will install the package locally to your machine. Note that the install with the "-e" option allows the package to be updated in real time, so modifications to the source code will be live as soon as you make them. This is ideal for development / research environments and is the most common way to install and use RIPLpy.

```bash
# Standard install
pip install .
```

If you want the source to stay exactly as it was on the repository, use the above command instead. Note the lack of the "-e" option. This installation method is intended for production environments.

### Optional features

`pandas` and `numpy` are only required for their respective DataFrame and
array exporters (`to_dataframe()` / `to_numpy()` and their `riplpy.*`
wrappers); everything else, including the CSV/JSON/list/ASCII exporters, works
without either. Install them via extras:

```bash
pip install -e ".[dataframe]"     # pandas, for to_dataframe()
pip install -e ".[numpy]"         # numpy, for to_numpy()
pip install -e ".[dataframe,numpy]"   # both
pip install -e ".[all]"           # all optional dependencies
```

(Drop the `-e` for a non-editable install.)

## Quick Start

```python
import riplpy

# Configure the RIPL database path
riplpy.set_path('/path/to/RIPL-4')

# Load all sections (may take time for large databases)
riplpy.load()

# Or load specific sections
import riplpy
import riplpy.masses as masses
masses.load(riplpy.get_path())
```

## Configuration

RIPLpy searches for the RIPL database in this order:
1. Path set via `riplpy.set_path()`
2. `RIPL_LOCATION` environment variable
3. `~/.riplpyrc` configuration file

```python
# Option 1: Set path programmatically
riplpy.set_path('/path/to/RIPL-4')

# Option 2: Environment variable
import os
os.environ['RIPL_LOCATION'] = '/path/to/RIPL-4'

# Check current path
print(riplpy.get_path())
```

## Core API

### Working with Nuclei

```python
from riplpy import Nuclide, Nucleus, Element

# Create nuclei using different specifications
pb208 = Nuclide(Z=82, A=208)       # By Z and A
sn132 = Nuclide(Z=50, N=82)        # By Z and N
fe56  = Nuclide(sym='Fe', A=56)    # By symbol and A

# Access properties
print(pb208.Z)              # 82
print(pb208.A)              # 208
print(pb208.N)              # 126
print(pb208.element_symbol) # 'Pb'
print(pb208.element_name)   # 'Lead'

# Nucleus is an alias for Nuclide with N-based specification
nucleus = Nucleus(Z=50, N=82)  # Sn-132
```

### Database Discovery

```python
import riplpy

riplpy.load()

# List available sections
print(riplpy.list_sections())
# ['densities', 'fission', 'gamma', 'levels', 'masses', 'optical', 'resonances']

# List databases in a section
print(riplpy.list_databases('masses'))
# ['ame20', 'frdm12', 'frdm95', 'hfb14', 'hfb27', 'bskg3', 'd1m', ...]

# Check if a nucleus is in RIPL
n = Nuclide(Z=50, A=132)
print(riplpy.in_ripl(n))
# True

# Find which sections contain a nucleus
print(riplpy.in_sections(n))
# ['densities', 'gamma', 'levels', 'masses']

# Find specific databases containing a nucleus
print(riplpy.in_dbs(n))
# [('masses', 'ame20'), ('levels', 'discrete_levels'), ...]
```

### Convenience Functions

RIPLpy provides high-level functions for common operations:

```python
import riplpy

riplpy.load()

# Mass data
mass_excess = riplpy.get_mass(82, 208)                    # Pb-208 mass excess [MeV]
mass_excess = riplpy.get_mass(82, 208, model='frdm12')    # Using FRDM2012 model
entry = riplpy.get_mass_entry(82, 208)                    # Full mass entry

# Level density parameters
ld_params = riplpy.get_level_density(50, 120, model='egsm')  # EGSM model for Sn-120

# Giant Dipole Resonance (RIPL-3 legacy theoretical GDR; full distribution only)
gdr = riplpy.get_gdr(82, 208)                             # GDR parameters for Pb-208

# Resonance spacing
res = riplpy.get_resonance(92, 238, wave='s')             # S-wave resonances for U-238

# Fission barriers
barrier = riplpy.get_fission_barrier(92, 235, model='bskg3')  # BSkG3 barrier for U-235

# Optical model potentials
omp = riplpy.get_omp(2405)                                # Get OMP by reference number
omps = riplpy.find_omp('n', 82, 208, E=14.0)              # Find OMPs for n + Pb-208 at 14 MeV
irefs = riplpy.list_omps(projectile='n')                  # List all neutron OMPs
```

### Batch Operations

Process multiple nuclei efficiently:

```python
import riplpy

riplpy.load()

# Get masses for multiple nuclei
nuclei = [(82, 208), (82, 206), (50, 132), (26, 56)]
masses = riplpy.get_masses(nuclei, model='ame20')

# Skip missing nuclei instead of raising errors
masses = riplpy.get_masses(nuclei, model='ame20', skip_missing=True)

# Batch level densities
ld_params = riplpy.get_level_densities(nuclei, model='egsm', skip_missing=True)

# Batch fission barriers
actinides = [(92, 235), (92, 238), (94, 239), (94, 241)]
barriers = riplpy.get_fission_barriers(actinides, model='bskg3', skip_missing=True)
```

### Direct Database Access

For more control, access databases directly:

```python
import riplpy
import riplpy.masses as masses
from riplpy import Nuclide

masses.load(riplpy.get_path())

# Access specific database
ame = masses.db.ame20
frdm = masses.db.frdm2012

# Get entry for a nucleus
n = Nuclide(Z=82, A=208)
entry = ame.get(n)
print(entry.Mexp)  # Experimental mass excess [MeV]
print(entry.Err)   # Uncertainty [MeV]

# Filter databases
pb_isotopes = ame.filter(Z=82)
heavy_nuclei = ame.filter(lambda e: e.n.A > 200)

# Iterate over database
for nuclide, entry in ame.data.items():
    if entry.Mexp is not None:
        print(f"{nuclide.sym}-{nuclide.A}: {entry.Mexp:.3f} MeV")
```

### RIPL-4 Array Datasets

In addition to the scalar parameter tables, the RIPL-4 release ships large
per-nucleus / per-Z **array (spectral/tabular) datasets**. These are loaded on
demand (the global `riplpy.load()` skips them by default for speed) via
per-element / per-nucleus helpers, and every entry is a uniform
`PacketEntry` supporting both `entry['U']` and `entry.U` access plus
`.as_dict()` / `to_dataframe()` export:

```python
import riplpy
from riplpy import Nuclide
riplpy.set_path('/path/to/RIPL-4')

# Gamma microscopic strength functions
import riplpy.gamma as gamma
d1m  = gamma.d1m.load_element(Z=26)                       # D1M+QRPA E1, per Z
smlo = gamma.smlo_e1.load_nucleus(Nuclide(Z=82, A=208))   # SMLO E1 photoabsorption
tlo  = gamma.tlo.load_element(Z=34)                       # TLO E1 tables
psf  = gamma.psf.load_category('oslo')                    # experimental PSF database

# Fission paths and saddle-point level densities
import riplpy.fission as fission
paths = fission.hfbpath.load_d1m()                   # HFB-D1M fission paths
rmf   = fission.rmf.load_axial()                     # RMF axial paths
nld   = fission.nld_fis.load('Max1')                 # saddle-point NLD

# Combinatorial level densities
import riplpy.densities as densities
comb  = densities.bsk14_comb.load_element(Z=26)      # BSk14 combinatorial NLD

# Optical model potentials — select by incident particle (n, p, d, t, He3, alpha)
import riplpy.optical as optical
optical.load()
omps   = optical.find_for_reaction('alpha', Z=26, A=56, E=20.0)  # alpha-OMPs for Fe-56
atomki = optical.atomki.load_nucleus(26, 56)         # ATOMKI / TALYS alphaomp9 alpha-OMP

# Every array DB exports uniformly for ML pipelines:
df = d1m.to_dataframe()
records = d1m.to_list()
```

See `CATALOG.md` for the human-readable dataset/reader map and `SCHEMA.json`
for a machine-readable field/units catalog.

### Writing & Exporting Data

Every loaded database supports the same export methods plus the native ASCII
writer, so the same recipes work for masses, level densities, fission
barriers, resonances, levels, and gamma-strength tables alike. Use these for
analysis pipelines, spreadsheet hand-off, ML ingestion, or to publish a
modified subset of RIPL back out to disk in the original format.

```python
import riplpy
riplpy.load()

db = riplpy.get_database('masses.ame20')          # any database works

# 1) Native RIPL ASCII (round-trip)
db.save('ame20_copy.dat')

# 2) Spreadsheet-friendly CSV (lists are stringified)
db.to_csv('ame20.csv')

# 3) JSON for general tooling
db.to_json('ame20.json', indent=2)

# 4) pandas DataFrame (Z, A, symbol always come first; array fields stay
#    as native Python lists). Requires the [dataframe] extra.
df = db.to_dataframe()
df.to_parquet('ame20.parquet')                    # if you have pyarrow

# 5) numpy array (Z, A, symbol first). Requires the [numpy] extra.
arr  = db.to_numpy()                              # structured (record) array
arr['Mexp']                                       # -> float64 column
X, cols = db.to_numpy(structured=False)           # dense float64 feature matrix

# 6) Plain list of dicts (Nuclide -> {Z, A, symbol})
records = db.to_list()                            # nested structures preserved
flat    = db.to_flat_list()                       # nested -> dotted keys

# Top-level convenience for ML/AI pipelines: one call, any section
df   = riplpy.to_dataframe('fission.bskg3_barriers')   # [dataframe] extra
arr  = riplpy.to_numpy('resonances.swave')             # [numpy] extra
rec  = riplpy.to_records('resonances.swave')           # no extra needed
```

For ML/AI feature engineering, `to_dataframe` and `to_numpy` are the two
primary forms — a labelled table or a raw array:

```python
import numpy as np
import riplpy
riplpy.load()

# (a) Scalar table -> DataFrame or feature matrix
db = riplpy.get_database('masses.ame20')
df = db.to_dataframe()                  # pandas: df['Mexp'], df['n.Z'], ...
X, cols = db.to_numpy(structured=False) # numpy: dense (n_nuclei, n_features)

# (b) Spectral table -> stacked (n_nuclei, n_energy) array for training
import riplpy.gamma as gamma
fe  = gamma.d1m.load_element(Z=26)      # per-Z D1M+QRPA E1 strengths
arr = fe.to_numpy()                     # arr['fE1'] is an object column of lists
spectra = np.stack(arr['fE1'])          # (n_iron_isotopes, 300) float matrix
energy  = np.asarray(arr['U'][0])       # shared 0.1-30 MeV grid
# `spectra` (and the matching DataFrame, fe.to_dataframe()) feed a model directly
```

> The DataFrame exporters need `pip install ".[dataframe]"` and the numpy
> exporters need `pip install ".[numpy]"`; the CSV/JSON/list/ASCII exporters
> have no optional dependencies.

Filter, then write a derived subset back out as ASCII (round-trips through
the section's writer):

```python
import riplpy
import riplpy.masses as masses
masses.load(riplpy.get_path())

heavy = masses.db.ame20.filter(lambda e: e.n.A > 200)
heavy.save('ame20_heavy_only.dat')                # native ASCII, same format

import riplpy.fission as fission
fission.load(riplpy.get_path())
high = fission.db.bskg3_barriers.filter(
    lambda e: e.inner['E[MeV]'] > 7.0
)
high.to_csv('bskg3_high_inner.csv')
```

The gamma-strength `gsf` database can also write per-element or per-nucleus
files (useful when staging a customised RIPL tree on disk):

```python
import riplpy.gamma as gamma
gamma.load(riplpy.get_path(), include_heavy=True)
gamma.db.gsf.save_all('out/gsf_by_Z')             # one file per element (zNNN.dat)
gamma.db.gsf.save_entries_by_Z('out/gsf_by_nuc')  # one file per nucleus
```

And the per-element loaders that back the heavy spectral datasets are also
fully exportable — they expose the same uniform `to_dataframe()` /
`to_list()` interface as the scalar tables:

```python
import riplpy.densities as densities
fe = densities.bsk14_comb.load_element(Z=26)
fe.to_dataframe().to_csv('bsk14_fe.csv', index=False)
fe.to_json('bsk14_fe.json')
```

## Section-Specific Usage

Each section has its own README with detailed examples:

- [Masses Section](riplpy/masses/README.md) - Evaluated masses and mass models
- [Densities Section](riplpy/densities/README.md) - Level density parameters
- [Fission Section](riplpy/fission/README.md) - Fission barrier data
- [Gamma Section](riplpy/gamma/README.md) - GDR and gamma strength functions
- [Levels Section](riplpy/levels/README.md) - Discrete nuclear levels
- [Resonances Section](riplpy/resonances/README.md) - Neutron resonance data
- [Optical Section](riplpy/optical/README.md) - Optical model potentials

## Examples

### Calculate Q-values

```python
import riplpy

riplpy.load()

# Neutron mass excess in MeV
NEUTRON_MASS_EXCESS = 8.071

def capture_q_value(Z_target, A_target):
    """Calculate Q-value for (n,gamma) capture reaction."""
    m_target = riplpy.get_mass(Z_target, A_target)
    m_product = riplpy.get_mass(Z_target, A_target + 1)
    # Q = M_target + M_n - M_product (gamma has zero mass)
    return m_target + NEUTRON_MASS_EXCESS - m_product

# Calculate Q-value for 208Pb(n,gamma)209Pb
q = capture_q_value(82, 208)
print(f"Q-value for Pb-208(n,g): {q:.3f} MeV")

# Calculate Q-value for a (p,n) reaction: 56Fe(p,n)56Co
def pn_q_value(Z, A):
    """Calculate Q-value for (p,n) reaction."""
    m_target = riplpy.get_mass(Z, A)
    m_proton = riplpy.get_mass(1, 1)
    m_residual = riplpy.get_mass(Z + 1, A)
    m_neutron = NEUTRON_MASS_EXCESS
    return (m_target + m_proton) - (m_residual + m_neutron)

q_pn = pn_q_value(26, 56)
print(f"Q-value for Fe-56(p,n)Co-56: {q_pn:.3f} MeV")
```

### Compare Mass Models

```python
import riplpy

riplpy.load()

# Compare mass predictions for neutron-rich Sn isotopes
models = ['ame20', 'frdm12', 'hfb27']
isotopes = [(50, A) for A in range(120, 140)]

for Z, A in isotopes:
    values = []
    for model in models:
        try:
            m = riplpy.get_mass(Z, A, model=model)
            values.append(f"{m:8.3f}")
        except:
            values.append("    N/A ")
    print(f"Sn-{A}: " + " | ".join(values))
```

### Find Optical Potentials for a Reaction

```python
import riplpy

riplpy.load()

# Find all neutron OMPs valid for U-238 at 1 MeV
potentials = riplpy.find_omp('n', Z=92, A=238, E=1.0)

print(f"Found {len(potentials)} potentials for n + U-238 at 1 MeV:")
for pot in potentials[:5]:  # Show first 5
    print(f"  iref={pot.iref}: {pot.header.author[:40]}...")
    print(f"    Model: {pot.flags.model_name}, E: {pot.validity.E_min}-{pot.validity.E_max} MeV")
```

## Requirements

- Python 3.10+
- `fortranformat`
- `pandas` (optional, for `to_dataframe()`; install via the `dataframe` extra)
- `numpy` (optional, for `to_numpy()`; install via the `numpy` extra)

Install the required dependencies with `pip install .`, and the optional
DataFrame/array exporters with `pip install ".[dataframe,numpy]"` (or
`pip install ".[all]"`). See `requirements.txt` for the required dependency
list.

## Author

Matthew Mumpower ([matthew@mumpower.net](mailto:matthew@mumpower.net))

## License

Distributed under the BSD-3 License. See `LICENSE.txt` for details.
