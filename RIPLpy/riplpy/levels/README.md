# Levels Section

Discrete nuclear level schemes and level density fitting parameters.

## Overview

The levels section provides comprehensive discrete level information for nuclei across the nuclear chart, including:

- **Discrete Levels** - Complete level schemes with energies, spins, parities, and gamma-ray transitions
- **Constant Temperature** - CT model parameters fitted to discrete level data

**Available Databases:**
- **discrete_levels** - Full discrete level schemes
- **constant_temperature** - CT level density parameters

## Quick Start

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

# Load levels data (may take some time)
levels.load(riplpy.get_path())

# Get discrete levels for Fe-56
fe56 = Nuclide(Z=26, A=56)
level_scheme = levels.db.discrete_levels.get(fe56)

print(f"Fe-56 has {level_scheme.num_levels} discrete levels")
# Ground state is the first level (index 0)
gs = level_scheme.levels[0]
parity = '+' if gs.parity > 0 else '-'
print(f"Ground state: J^π = {gs.spin}{parity}")
```

## Examples

### Basic Level Information

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

levels.load(riplpy.get_path())

# Get levels for Ni-60
ni60 = Nuclide(Z=28, A=60)
scheme = levels.db.discrete_levels.get(ni60)

print(f"Ni-60 discrete levels:")
print(f"  Total levels: {scheme.num_levels}")
print(f"  Max known energy: {max(scheme.level_energies):.3f} MeV")
print(f"\nFirst 10 levels:")
print(f"{'#':>3} {'E (MeV)':>10} {'J':>6} {'π':>3}")

for i, level in enumerate(scheme.levels[:10]):
    parity = '+' if level.parity > 0 else '-'
    spin = f"{level.spin:.1f}" if level.spin is not None else "?"
    print(f"{i:3d} {level.energy:10.4f} {spin:>6} {parity:>3}")
```

### Ground State Properties

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

levels.load(riplpy.get_path())

# Get ground state spin/parity for various nuclei
nuclei = [
    (26, 56, "Fe-56"),
    (28, 60, "Ni-60"),
    (50, 120, "Sn-120"),
    (82, 208, "Pb-208"),
    (92, 238, "U-238"),
]

print("Ground state spins and parities:")
for Z, A, name in nuclei:
    n = Nuclide(Z=Z, A=A)
    try:
        scheme = levels.db.discrete_levels.get(n)
        gs = scheme.levels[0]  # Ground state is first level
        parity = '+' if gs.parity > 0 else '-'
        print(f"  {name:>10}: J^π = {gs.spin}{parity}")
    except KeyError:
        print(f"  {name:>10}: Not available")
```

### Gamma-Ray Transitions

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

levels.load(riplpy.get_path())

# Get gamma-ray transitions in Co-60 decay
co60 = Nuclide(Z=27, A=60)
scheme = levels.db.discrete_levels.get(co60)

# Gamma transitions are stored on the scheme as a dict keyed by the
# (1-based) level number. Level 1 is the ground state; gamma-emitting
# excited states start at key 2.
print(f"Co-60 gamma-ray transitions from first excited states:")
for idx, level in enumerate(scheme.levels[1:5], start=1):  # skip ground state
    transitions = scheme.gammas.get(idx + 1, [])
    if not transitions:
        continue
    print(f"\nLevel at {level.energy:.4f} MeV (J = {level.spin}):")
    for gamma in transitions:
        print(f"  → E_γ = {gamma.energy:.4f} MeV, prob = {gamma.prob_g:.4f}")
```

### Find Isomeric States

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

levels.load(riplpy.get_path())

# Find long-lived isomers in a nucleus
hf178 = Nuclide(Z=72, A=178)
scheme = levels.db.discrete_levels.get(hf178)

print("Hf-178 isomeric states (T_1/2 > 1 μs):")
for level in scheme.levels:
    if level.halflife is not None and level.halflife > 1e-6:
        parity = '+' if level.parity > 0 else '-'
        print(f"  E = {level.energy:.4f} MeV, J^π = {level.spin}{parity}, T_1/2 = {level.halflife:.2e} s")
```

### Level Density from Discrete Levels

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

levels.load(riplpy.get_path())

# Count levels vs excitation energy to estimate level density
fe56 = Nuclide(Z=26, A=56)
scheme = levels.db.discrete_levels.get(fe56)

# Bin levels by energy
bin_width = 0.5  # MeV
max_energy = 8.0  # MeV
bins = {}

for level in scheme.levels:
    if level.energy < max_energy:
        bin_idx = int(level.energy / bin_width)
        bins[bin_idx] = bins.get(bin_idx, 0) + 1

print("Fe-56 level count vs excitation energy:")
print(f"{'E (MeV)':>10} {'N_levels':>10} {'ρ (MeV^-1)':>12}")
for bin_idx in sorted(bins.keys()):
    E = (bin_idx + 0.5) * bin_width
    N = bins[bin_idx]
    rho = N / bin_width
    print(f"{E:10.2f} {N:10d} {rho:12.1f}")
```

### Constant Temperature Parameters

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

levels.load(riplpy.get_path())

# Get CT parameters fitted to discrete levels
fe56 = Nuclide(Z=26, A=56)
ct = levels.db.constant_temperature.get(fe56)

print(f"Fe-56 Constant Temperature parameters:")
print(f"  Temperature T:    {ct.T:.3f} +/- {ct.dT:.3f} MeV")
print(f"  Energy shift U0:  {ct.U0:.3f} +/- {ct.dU0:.3f} MeV")
print(f"  Matching energy Umax: {ct.Umax:.3f} MeV, levels Nmax = {ct.Nmax}")
```

### Decay Modes

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

levels.load(riplpy.get_path())

# Get decay modes for a radioactive nucleus
cs137 = Nuclide(Z=55, A=137)
scheme = levels.db.discrete_levels.get(cs137)

gs = scheme.levels[0]  # Ground state
print(f"Cs-137 ground state decay modes:")
for decay in gs.decay_modes:
    print(f"  {decay.mode}: {decay.percent:.2f}%")
```

### Filter Levels by Properties

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

levels.load(riplpy.get_path())

# Get Fe-56 levels
fe56 = Nuclide(Z=26, A=56)
scheme = levels.db.discrete_levels.get(fe56)

# Find all 2+ states
two_plus = [l for l in scheme.levels if l.spin == 2.0 and l.parity > 0]
print(f"Fe-56 has {len(two_plus)} levels with J^π = 2+")
for l in two_plus[:5]:
    print(f"  E = {l.energy:.4f} MeV")

# Find high-spin states
high_spin = [l for l in scheme.levels if l.spin is not None and l.spin >= 6]
print(f"\nFe-56 has {len(high_spin)} levels with J ≥ 6")
```

### Complete Level Information

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

levels.load(riplpy.get_path())

# Get detailed information for a specific level
pb208 = Nuclide(Z=82, A=208)
scheme = levels.db.discrete_levels.get(pb208)

# First excited state
level = scheme.levels[1]
print(f"Pb-208 first excited state:")
print(f"  Energy: {level.energy:.4f} MeV")
print(f"  Spin: {level.spin}")
print(f"  Parity: {'+' if level.parity > 0 else '-'}")
print(f"  Half-life: {level.halflife} s" if level.halflife else "  Half-life: stable")
print(f"  Spin flag: {level.jflag}")
print(f"  Number of gammas: {level.num_gammas}")
```

## Entry Fields

### Level Scheme Entry
| Field | Description |
|-------|-------------|
| `nucleus` | Nuclide object |
| `num_levels` | Total number of discrete levels |
| `levels` | List of level objects (ground state is `levels[0]`) |
| `gammas` | List of all gamma transitions |
| `extras` | Additional nuclear data |

### Discrete Level Object
| Field | Description |
|-------|-------------|
| `energy` | Excitation energy [MeV] |
| `spin` | Spin J (float or None) |
| `parity` | Parity (+1 or -1) |
| `halflife` | Half-life [s] or None |
| `jflag` | Spin estimation flag |
| `num_gammas` | Number of outgoing gamma transitions (transitions live on `scheme.gammas[level_number]`) |
| `decay_modes` | List of DecayMode objects |

### Gamma Transition Object
| Field | Description |
|-------|-------------|
| `to_level_number` | Final level index |
| `energy` | Transition energy [MeV] |
| `prob_g` | Gamma emission probability |
| `prob_em` | EM transition probability |
| `cc` | Internal conversion coefficient |

### CT Entry
| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `T`, `dT` | Nuclear temperature and 1-sigma uncertainty [MeV] |
| `U0`, `dU0` | Energy backshift and 1-sigma uncertainty [MeV] |
| `Nlev`, `Nmax`, `N0`, `Nc` | Level counts used in the CT fit |
| `Umax`, `Uc` | Matching energy and cutoff energy [MeV] |
| `Chi`, `Fit`, `Flag` | Fit quality and provenance flags |
| `sigma` | Spin cutoff [MeV] |

## Writing & Exporting

Both the discrete-level and constant-temperature databases support the
standard export interface plus native ASCII writers:

```python
import riplpy
import riplpy.levels as levels
levels.load(riplpy.get_path())

# Native ASCII round-trip (discrete-level scheme files are large)
levels.db.discrete_levels.save('discrete_levels_copy.dat')

# Machine-readable exports for the CT parameter fits
levels.db.constant_temperature.to_csv('ct.csv')
levels.db.constant_temperature.to_json('ct.json')
df = levels.db.constant_temperature.to_dataframe()

# Per-element discrete-level files (zNNN.dat) can also be loaded directly
fe = levels.discrete.load_element(Z=26, directory=riplpy.get_path())
```

## Using with Main API

The levels section is primarily accessed through direct database access since the data structures are complex.

```python
import riplpy
import riplpy.levels as levels
from riplpy import Nuclide

levels.load(riplpy.get_path())

# Access discrete levels
scheme = levels.db.discrete_levels.get(Nuclide(Z=26, A=56))

# Access CT parameters
ct = levels.db.constant_temperature.get(Nuclide(Z=26, A=56))
```

## Physical Background

Nuclear levels form a discrete spectrum at low excitation energies, transitioning to a quasi-continuum at higher energies. The level density increases approximately exponentially with excitation energy.

**Spin-Parity Assignments:**
- Measured values are considered reliable
- Estimated values are flagged with uncertainty indicators
- Common estimation methods: shell model, systematics, reaction data

**Level Density at Low Energy:**
The constant temperature model describes level density near the ground state:
```
ρ(U) = (1/T) × exp((U - E₀)/T)
```
where T is the nuclear temperature and E₀ is the energy backshift.
