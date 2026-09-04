# Densities Section

Nuclear level density parameters from various models.

## Overview

The densities section provides level density model parameters used in nuclear reaction calculations. Level densities describe the number of nuclear states per unit energy interval.

**RIPL-4 Models (GitHub release + full distribution):**
- **EGSM** - Enhanced Generalized Superfluid Model
- **EGSM_NORM** - EGSM normalization parameters
- **BSk14-comb** - BSk14 microscopic combinatorial level densities
- **BSkG3-comb** - BSkG3 microscopic combinatorial level densities
- **QRPA-BE** - QRPA-based microscopic level densities
- **THFB-comb** - Temperature-dependent HFB combinatorial level densities
- **Shell Corrections (MS)** - Myers-Swiatecki shell correction energies

**RIPL-3 Legacy Models (full RIPL distribution only — NOT in the RIPL-4 GitHub release):**
- **BSFG** - Back-Shifted Fermi Gas model
- **CT** - Constant Temperature model
- **HFB** - Hartree-Fock-Bogoliubov spin-dependent level densities
- **Shell Corrections (MK)** - Mengoni-Nakajima shell corrections

Legacy models are skipped with a warning by `densities.load()` when their
data files are absent; their readers/writers remain available for users with
the full RIPL distribution.

## Quick Start

```python
import riplpy
import riplpy.densities as densities
from riplpy import Nuclide

# Load level density data
densities.load(riplpy.get_path())

# Get EGSM parameters for Fe-57 (RIPL-4; works with the GitHub release)
fe57 = Nuclide(Z=26, A=57)
entry = densities.db.egsm.get(fe57)
print(f"Level density parameter a: {entry.a:.3f} MeV^-1")
```

## Available Models

| Model | Attribute | Description | Availability |
|-------|-----------|-------------|--------------|
| EGSM | `egsm` | Enhanced Generalized Superfluid | RIPL-4 |
| EGSM Norm | `egsm_norm` | EGSM normalization | RIPL-4 |
| BSk14-comb | `bsk14_comb` | BSk14 combinatorial | RIPL-4 |
| BSkG3-comb | `bskg3_comb` | BSkG3 combinatorial | RIPL-4 |
| QRPA-BE | `qrpabe` | QRPA-based microscopic | RIPL-4 |
| THFB-comb | `thfb_comb` | Temperature-dependent HFB | RIPL-4 |
| Shell Corr. (MS) | `shellcorr_ms` | Myers-Swiatecki shell corrections | RIPL-4 |
| BSFG | `bsfg` | Back-Shifted Fermi Gas | RIPL-3 legacy |
| CT | `ct` | Constant Temperature | RIPL-3 legacy |
| HFB | `hfb` | HFB spin-dependent level densities | RIPL-3 legacy |
| Shell Corr. (MK) | `shellcorr_mk` | Mengoni-Nakajima shell corrections | RIPL-3 legacy |

## Examples

All examples below use RIPL-4 models and run against the RIPL-4 GitHub
release. RIPL-3 legacy usage is shown in a clearly-marked subsection at the
end.

### Basic Level Density Parameters (EGSM)

```python
import riplpy
import riplpy.densities as densities
from riplpy import Nuclide

densities.load(riplpy.get_path())

# Get EGSM parameters for Fe-57
fe57 = Nuclide(Z=26, A=57)
egsm = densities.db.egsm.get(fe57)

print("Fe-57 EGSM parameters:")
print(f"  Ground state spin: {egsm.Io}")
print(f"  Neutron binding energy: {egsm.Bn:.3f} MeV")
print(f"  Experimental s-wave spacing Do: {egsm.Do:.3f} keV")
print(f"  Level density parameter a: {egsm.a:.3f} MeV^-1")
print(f"  Shell correction: {egsm.Esh:.3f} MeV")
print(f"  Deformation: {egsm.Def}")
print(f"  Systematic a_sys: {egsm.a_sys} MeV^-1  (a_exp/a_sys = {egsm.a_ratio})")
```

### Apply EGSM Normalization Factors

```python
import riplpy
import riplpy.densities as densities
from riplpy import Element

densities.load(riplpy.get_path())

# EGSM normalization factors are keyed by Element, not Nuclide.
# They scale the systematic 'a' for nuclei without experimental Do.
for Z in (20, 26, 50, 82):
    ele = Element(Z=Z)
    try:
        norm = densities.db.egsm_norm.get(ele)
        print(f"{norm.symbol} (Z={norm.Z}): normalization factor = {norm.factor:.3f}")
    except KeyError:
        print(f"Z={Z}: no normalization factor")
```

### Compare EGSM Across an Isotopic Chain

```python
import riplpy
import riplpy.densities as densities
from riplpy import Nuclide

densities.load(riplpy.get_path())

# Compare the EGSM level density parameter for Sn isotopes
print("Sn isotope EGSM level density parameters:")
print(f"{'A':>4} {'a (MeV^-1)':>12} {'a_sys':>8} {'Esh (MeV)':>10}")

for A in range(112, 126, 2):
    sn = Nuclide(Z=50, A=A)
    try:
        egsm = densities.db.egsm.get(sn)
    except KeyError:
        continue
    a_sys = "n/a" if egsm.a_sys is None else f"{egsm.a_sys:8.3f}"
    print(f"{A:4d} {egsm.a:12.3f} {a_sys:>8} {egsm.Esh:10.3f}")
```

### Calculate Level Density (Fermi Gas, EGSM parameter)

```python
import riplpy
import riplpy.densities as densities
from riplpy import Nuclide
import math

densities.load(riplpy.get_path())

def fermi_gas_level_density(Z, A, U):
    """
    Estimate the level density with the Fermi gas formula using the
    EGSM level density parameter 'a'.

    Args:
        Z: Atomic number
        A: Mass number
        U: Effective excitation energy [MeV]

    Returns:
        Level density [MeV^-1]
    """
    n = Nuclide(Z=Z, A=A)
    params = densities.db.egsm.get(n)

    a = params.a  # EGSM level density parameter [MeV^-1]

    if U <= 0:
        return 0.0

    # Fermi gas formula
    sqrt_aU = math.sqrt(a * U)
    rho = math.exp(2 * sqrt_aU) / (12 * math.sqrt(2) * a**0.25 * U**1.25)

    return rho

# Calculate level density for Fe-57 at various excitation energies
print("Fe-57 Fermi gas level density vs excitation energy:")
print(f"{'U (MeV)':>10} {'rho (MeV^-1)':>15}")
for U in [1, 2, 5, 10, 15, 20]:
    rho = fermi_gas_level_density(26, 57, U)
    print(f"{U:10.1f} {rho:15.2e}")
```

### Microscopic Combinatorial Level Densities (BSk14-comb)

```python
import riplpy
import riplpy.densities as densities
from riplpy import Nuclide

densities.load(riplpy.get_path())

# The combinatorial tables (bsk14_comb, bskg3_comb, qrpabe, thfb_comb) are
# spin- and parity-dependent. densities.load() populates them automatically.
fe56 = Nuclide(Z=26, A=56)
entry = densities.db.bsk14_comb.get(fe56)

# Access the positive-parity table
pp = entry.positive_parity
print("Fe-56 BSk14-comb level density (positive parity, first few values):")
print(f"{'U (MeV)':>10} {'rho_tot':>12} {'T (MeV)':>10}")
for i in range(min(10, len(pp.U))):
    print(f"{pp.U[i]:10.2f} {pp.rho_tot[i]:12.2e} {pp.T[i]:10.3f}")

# Total level density (both parities) at a given excitation energy
rho_5 = entry.get_total_density(5.0)
print(f"\nTotal level density at U = 5 MeV: {rho_5:.3e} MeV^-1")
```

### Load a Single Element (per-Z combinatorial loader)

```python
import riplpy
import riplpy.densities as densities
from riplpy import Nuclide

# Per-element loaders avoid reading the whole library; useful for large
# combinatorial datasets. load_all(directory) loads every element.
db = densities.bsk14_comb.load_element(Z=26, directory=riplpy.get_path())
print(f"BSk14-comb iron isotopes loaded: {len(db.data)}")

fe56 = Nuclide(Z=26, A=56)
entry = db.get(fe56)
data = entry.positive_parity.get_density_at_energy(8.0)
print(f"Fe-56 near U = 8 MeV: U={data['U']:.2f} MeV, rho_tot={data['rho_tot']:.3e}")
```

### Shell Corrections (Myers-Swiatecki)

```python
import riplpy
import riplpy.densities as densities
from riplpy import Nuclide

densities.load(riplpy.get_path())

# Myers-Swiatecki shell corrections (RIPL-4). Equivalent direct loader:
#   densities.shell_corr.load_ms(riplpy.get_path())
nuclei = [
    (20, 41, "Ca-41"),
    (26, 57, "Fe-57"),
    (50, 119, "Sn-119"),
    (82, 209, "Pb-209"),
]

print("Shell correction energies (Myers-Swiatecki):")
for Z, A, name in nuclei:
    n = Nuclide(Z=Z, A=A)
    try:
        entry = densities.db.shellcorr_ms.get(n)
        print(f"  {name}: shell = {entry.shell:.3f} MeV, corr = {entry.corr:.3f} MeV, "
              f"beta2 = {entry.beta2:.3f}, beta4 = {entry.beta4:.3f}")
    except KeyError:
        print(f"  {name}: Not available")
```

### Filter Nuclei by Properties

```python
import riplpy
import riplpy.densities as densities

densities.load(riplpy.get_path())

# Find nuclei with a relatively large EGSM level density parameter
large_a = densities.db.egsm.filter(lambda e: e.a is not None and e.a > 20)
print(f"Nuclei with a > 20 MeV^-1: {len(large_a.data)}")

# Find nuclei with a negative shell correction
neg_shell = densities.db.egsm.filter(lambda e: e.Esh is not None and e.Esh < 0)
print(f"Nuclei with Esh < 0 MeV: {len(neg_shell.data)}")
```

### RIPL-3 Legacy Models (full RIPL distribution only)

> Note: The BSFG and CT models are part of the RIPL-3 legacy distribution.
> Their data files are NOT shipped in the RIPL-4 GitHub release, so
> `densities.load()` skips them with a warning and the snippets below will
> raise on the GitHub release. They only run against the full RIPL
> distribution.

```python
import riplpy
import riplpy.densities as densities
from riplpy import Nuclide

densities.load(riplpy.get_path())

# Back-Shifted Fermi Gas (BSFG) parameters for Ni-60
ni60 = Nuclide(Z=28, A=60)
bsfg = densities.db.bsfg.get(ni60)
print(f"Ni-60 BSFG: a = {bsfg.ainf:.3f} MeV^-1, "
      f"pairing = {bsfg.pairing:.3f} MeV, dW = {bsfg.dW:.3f} MeV")

# Compare BSFG and Constant Temperature (CT) for Sn isotopes
print("Sn isotope legacy parameters:")
print(f"{'A':>4} {'a_BSFG':>10} {'T_CT':>10}")
for A in range(112, 126, 2):
    sn = Nuclide(Z=50, A=A)
    try:
        bsfg = densities.db.bsfg.get(sn)
        ct = densities.db.ct.get(sn)
        print(f"{A:4d} {bsfg.ainf:10.3f} {ct.T:10.3f}")
    except KeyError:
        pass
```

The HFB spin-dependent tables (`densities.db.hfb`) are also RIPL-3 legacy;
the RIPL-4 GitHub release provides the combinatorial tables
(`bsk14_comb`, `bskg3_comb`, `qrpabe`, `thfb_comb`) instead. Both expose the
same `positive_parity` / `negative_parity` API shown above.

## Entry Fields

### EGSM Entry (RIPL-4)
| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `Io` | Ground state spin of target nucleus |
| `Bn` | Neutron binding energy / Qn [MeV] |
| `Do` | Experimental s-wave resonance spacing Dobs [keV] |
| `Derr` | Uncertainty on resonance spacing dDobs [keV] |
| `Esh` | Shell correction energy [MeV] |
| `Def` | Ground-state deformation parameter (may be `None` on older files) |
| `Dcalc` | Calculated s-wave resonance spacing [keV] (may be `None`) |
| `dap` | Upper uncertainty on level density parameter [MeV^-1] |
| `a` | Level density parameter at Bn [MeV^-1] |
| `dam` | Lower uncertainty on level density parameter [MeV^-1] |
| `a_sys` | Systematic level density parameter a_sys [MeV^-1] (may be `None`) |
| `a_ratio` | Ratio a_exp / a_sys (may be `None`) |

### EGSM Normalization Entry (RIPL-4)
| Field | Description |
|-------|-------------|
| `element` | Element object (entries are keyed by `Element`, not `Nuclide`) |
| `factor` | Normalization factor (a_exp/a_sys averaged over isotopes) |
| `Z`, `sym`, `symbol`, `name` | Convenience properties from the element |

### Combinatorial Entry (RIPL-4: BSk14-comb, BSkG3-comb, QRPA-BE, THFB-comb)
| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `positive_parity` | `ParityData` for positive-parity states (or `None`) |
| `negative_parity` | `ParityData` for negative-parity states (or `None`) |
| `get_total_density(U, parity=None)` | Total level density at energy `U` [MeV^-1] |

Each `ParityData` exposes the parallel lists `U`, `T`, `Ncumul`,
`rho_obs`, `rho_tot`, `rho_J` and the helper
`get_density_at_energy(U)`.

### Shell Correction Entry (RIPL-4: Myers-Swiatecki)
| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `shell` | Shell correction energy [MeV] |
| `corr` | Deformation correction energy [MeV] |
| `beta2` | Quadrupole deformation parameter |
| `beta4` | Hexadecapole deformation parameter |

### Legacy Entry Fields (RIPL-3, full distribution only)

**BSFG Entry**

| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `ainf` | Asymptotic level density parameter [MeV^-1] |
| `Io` | Ground state spin |
| `Bn` | Neutron binding energy [MeV] |
| `Do` | Average s-wave resonance spacing |
| `pairing` | Pairing energy shift [MeV] |
| `dW` | Shell correction energy [MeV] |

**CT Entry**

| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `T` | Nuclear temperature [MeV] |
| `E0` | Energy backshift [MeV] |
| `Ematch` | Matching energy [MeV] |

## Writing & Exporting

All level-density databases expose the standard export interface (CSV, JSON,
pandas, list-of-dicts) plus a native ASCII writer. The combinatorial tables
(`bsk14_comb`, `bskg3_comb`, `qrpabe`, `thfb_comb`) preserve their per-spin
`rho_J` arrays in DataFrames and JSON.

```python
import riplpy
import riplpy.densities as densities
densities.load(riplpy.get_path())

# Native ASCII round-trip
densities.db.egsm.save('egsm_copy.dat')

# Machine-readable exports
densities.db.egsm.to_csv('egsm.csv')
densities.db.egsm.to_json('egsm.json')
df = densities.db.egsm.to_dataframe()

# Filter + write a derived subset
big_a = densities.db.egsm.filter(lambda e: e.a is not None and e.a > 20)
big_a.to_csv('egsm_large_a.csv')

# Per-element combinatorial export (Z=26 = iron)
fe = densities.bsk14_comb.load_element(Z=26)
fe.to_dataframe().to_csv('bsk14_fe.csv', index=False)
```

## Using with Main API

The default level density model is `egsm` (RIPL-4).

```python
import riplpy

riplpy.load()

# Get EGSM level density parameters (model='egsm' is the default)
params = riplpy.get_level_density(26, 57)
print(f"a = {params.a} MeV^-1")

# Batch operations across nuclei
nuclei = [(26, 56), (28, 58), (50, 120)]
ld_params = riplpy.get_level_densities(nuclei, model='egsm', skip_missing=True)

# Combinatorial model via the main API
comb = riplpy.get_level_density(26, 56, model='bsk14_comb')
print(f"Fe-56 total rho at 5 MeV = {comb.get_total_density(5.0):.3e} MeV^-1")
```

## Physical Background

The nuclear level density ρ(U, J, π) describes the number of nuclear states per unit energy at excitation energy U with spin J and parity π. It is a crucial ingredient in statistical nuclear reaction calculations.

**BSFG Model:**
```
ρ(U) = exp(2√(aU)) / (12√2 · a^(1/4) · U^(5/4) · σ)
```
where `a` is the level density parameter and `σ` is the spin cutoff parameter.

**CT Model:**
At low energies, level density follows:
```
ρ(U) = (1/T) · exp((U - E0)/T)
```
where `T` is the nuclear temperature.
