# Optical Section

Optical model potentials for nuclear reaction calculations.

## Overview

The optical section provides access to the RIPL optical model potential (OMP) library, containing ~584 parameterizations for various projectile-target combinations. Optical potentials are essential for calculating elastic scattering, absorption cross sections, and transmission coefficients.

**Available Databases:**
- **index** - Metadata for all potentials (581 entries)
- **potentials** - Full OMP parameters (584 potentials)
- **deformations** - Excited-level deformation parameters (1423 entries)
- **references** - Bibliographic references (143 references)

**Supported Projectiles:**
- Neutrons (n)
- Protons (p)
- Deuterons (d)
- Tritons (t)
- Helions (³He)
- Alphas (α)

**Model Types:**
- Spherical (imodel=0): ~461 potentials
- Rigid Rotor (imodel=1): ~99 potentials
- Vibrational (imodel=2): ~10 potentials
- Soft Rotor (imodel=3): ~7 potentials
- Rigid-Soft (imodel=4): ~5 potentials
- Soft Deformed (imodel=5): ~2 potentials

## Quick Start

```python
import riplpy
import riplpy.optical as optical

# Load optical database
optical.load(riplpy.get_path())

# Get a specific potential by reference number
pot = optical.db.potentials.get(2405)
print(f"Potential iref=2405:")
print(f"  Author: {pot.header.author}")
print(f"  Projectile: {pot.projectile}")
print(f"  Valid energy: {pot.validity.E_min}-{pot.validity.E_max} MeV")
```

## Examples

### Find Potentials for a Reaction

```python
import riplpy
import riplpy.optical as optical

optical.load(riplpy.get_path())

# Find neutron OMPs for Pb-208 at 14 MeV
matches = optical.db.potentials.find_for_reaction('n', Z=82, A=208, E=14.0)

print(f"Found {len(matches)} potentials for n + Pb-208 at 14 MeV:")
for pot in matches[:5]:
    print(f"  iref={pot.iref}: {pot.header.author[:40]}")
    print(f"    Model: {pot.flags.model_name}")
    print(f"    E range: {pot.validity.E_min:.1f}-{pot.validity.E_max:.1f} MeV")
```

### Spherical OMP Parameters

```python
import riplpy
import riplpy.optical as optical

optical.load(riplpy.get_path())

# Get Koning-Delaroche global potential (iref=2405)
pot = optical.db.potentials.get(2405)

print(f"Koning-Delaroche neutron potential (iref=2405):")
print(f"  Valid for Z: {pot.validity.Z_min}-{pot.validity.Z_max}")
print(f"  Valid for A: {pot.validity.A_min}-{pot.validity.A_max}")
print(f"  Energy range: {pot.validity.E_min}-{pot.validity.E_max} MeV")
print(f"  Relativistic: {pot.flags.is_relativistic}")
print(f"  Dispersive: {pot.flags.uses_dispersion}")

# Access potential components
print(f"\nPotential has {len(pot.components)} components:")
for comp in pot.components:
    print(f"  {comp.name}: {len(comp.energy_ranges)} energy range(s)")
```

### Coupled-Channel Potentials

```python
import riplpy
import riplpy.optical as optical

optical.load(riplpy.get_path())

# Filter to rigid rotor (coupled-channel) potentials
cc_db = optical.db.potentials.filter_coupled_channel()
print(f"Found {len(cc_db.data)} coupled-channel potentials")

# Get a rigid rotor potential for actinides
rr_pots = optical.db.potentials.filter_by_model(1)  # imodel=1 = rigid rotor
print(f"\nRigid rotor potentials: {len(rr_pots.data)}")

# Look at one in detail
for iref, pot in list(rr_pots.data.items())[:1]:
    print(f"\niref={iref}: {pot.header.author}")
    print(f"  Number of isotopes: {pot.n_isotopes}")
    for iso in pot.isotopes[:2]:
        print(f"  Isotope Z={iso.Z}, A={iso.A}:")
        print(f"    Beta-2: {iso.beta2:.3f}")
        print(f"    Levels: {len(iso.levels)}")
```

### Filter by Projectile and Target

```python
import riplpy
import riplpy.optical as optical

optical.load(riplpy.get_path())

# Get all proton potentials
proton_db = optical.db.potentials.filter_by_projectile('p')
print(f"Proton potentials: {len(proton_db.data)}")

# Get potentials valid for U-238
u238_db = optical.db.potentials.filter_by_target(Z=92, A=238)
print(f"Potentials for U-238: {len(u238_db.data)}")

# Combine filters: neutron potentials for actinides
actinide_n = optical.db.potentials.filter_by_projectile('n')
actinide_n = actinide_n.filter_by_target(Z=92)
print(f"Neutron potentials for Z=92: {len(actinide_n.data)}")
```

### Access Potential Coefficients

```python
import riplpy
import riplpy.optical as optical

optical.load(riplpy.get_path())

# Get a potential and examine its Woods-Saxon parameters
pot = optical.db.potentials.get(2405)

# Real volume potential (component 0)
real_vol = pot.components[0]
print(f"{real_vol.name} parameters:")

# Each EnergyRange covers E < `epot` MeV.
for i, erange in enumerate(real_vol.energy_ranges):
    print(f"\n  Energy range {i}: E < {erange.epot} MeV")

    # Radius coefficients (R = r0 * A^(1/3) + ...); .values is the full
    # polynomial vector, .r0 is just the first (leading) term.
    rco = erange.radius
    print(f"    r0 = {rco.r0:.4f} fm")

    # Diffuseness coefficients
    aco = erange.diffuseness
    print(f"    a0 = {aco.a0:.4f} fm")

    # Potential strength (V = V0 + V1*E + ...): V0 is the leading term.
    pot_str = erange.strength
    print(f"    V0 = {pot_str.V0:.2f} MeV")
    if len(pot_str.values) > 1:
        print(f"    V1 = {pot_str.values[1]:.4f} MeV/MeV")
```

### Deformation Parameters

```python
import riplpy
import riplpy.optical as optical

optical.load(riplpy.get_path())

# Get deformation for Pu-239 first 2+ state
entries = optical.db.deformations.get_for_nucleus(94, 239)
print(f"Pu-239 deformation entries: {len(entries)}")

for entry in entries[:5]:
    par = '+' if entry.parity > 0 else '-'
    spin = entry.spin if entry.spin is not None else '?'
    print(f"  Ex={entry.Ex:.4f} MeV, J={spin}{par}, L={entry.L}, beta={entry.beta:.4f}")
```

### Bibliographic References

```python
import riplpy
import riplpy.optical as optical

optical.load(riplpy.get_path())

# Get reference for a potential
ref = optical.db.references.get(100)
print(f"Reference 100:")
print(f"  Citation: {ref.citation}")
print(f"  First author: {ref.first_author}")
print(f"  Year: {ref.year}")

# Search for references by author
koning_refs = optical.db.references.search_by_author("Koning")
print(f"\nReferences by Koning: {len(koning_refs.data)}")
```

### Modified Potentials

```python
import riplpy
import riplpy.optical as optical

# Load modified potentials (corrections to main database)
mod_db = optical.load_modified_potentials(riplpy.get_path())
print(f"Modified potentials: {len(mod_db.data)}")
print(f"Reference numbers: {sorted(mod_db.irefs)}")

# Load main database with modifications applied
db = optical.load_with_modifications(riplpy.get_path())
print(f"Total potentials (with modifications): {len(db.data)}")
```

### Database Statistics

```python
import riplpy
import riplpy.optical as optical

optical.load(riplpy.get_path())

# Get database summary
print(optical.db.potentials.info())
```

### Calculate Potential at a Point

```python
import riplpy
import riplpy.optical as optical
import math

optical.load(riplpy.get_path())

def woods_saxon(r, V0, r0, a0, A):
    """Calculate Woods-Saxon potential at radius r."""
    R = r0 * A**(1.0/3.0)
    return -V0 / (1.0 + math.exp((r - R) / a0))

# Get Koning-Delaroche potential
pot = optical.db.potentials.get(2405)

# Real volume parameters (approximate, energy-dependent)
A = 208  # Pb-208
E = 14.0  # MeV

# Get first energy range (simplified)
real_vol = pot.components[0].energy_ranges[0]
V0 = real_vol.strength.V0    # Base depth
r0 = real_vol.radius.r0      # Radius parameter
a0 = real_vol.diffuseness.a0  # Diffuseness

print(f"Woods-Saxon parameters for n + Pb-208:")
print(f"  V0 = {V0:.2f} MeV")
print(f"  r0 = {r0:.4f} fm")
print(f"  a0 = {a0:.4f} fm")
print(f"\nPotential vs radius:")
for r in [0, 2, 4, 6, 7, 8, 9, 10]:
    V = woods_saxon(r, V0, r0, a0, A)
    print(f"  r = {r:4.1f} fm: V = {V:8.2f} MeV")
```

## Entry Fields

### Spherical OMP
| Field | Description |
|-------|-------------|
| `iref` | Reference number |
| `header` | Author, title, reference info |
| `validity` | Z, A, E ranges of validity |
| `flags` | Model type, relativistic, dispersive flags |
| `components` | 6 potential components (V, W, Vso, Wso, etc.) |
| `coulomb` | Coulomb parameters |

### Coupled-Channel OMP
| Field | Description |
|-------|-------------|
| `iref` | Reference number |
| `header` | Author, title, reference info |
| `validity` | Z, A, E ranges of validity |
| `flags` | Model type flags |
| `components` | Potential components |
| `isotopes` | List of isotope-specific data |
| `n_isotopes` | Number of isotopes |

### Deformation Entry
| Field | Description |
|-------|-------------|
| `Z`, `A` | Nuclide |
| `Ex` | Excitation energy [MeV] |
| `spin`, `parity` | Level quantum numbers |
| `L` | Deformation order (2=quadrupole, 3=octupole) |
| `beta` | Deformation parameter |
| `reference` | Data source |

## Writing & Exporting

The optical-model databases (potentials, deformations, references) all
support the standard export interface:

```python
import riplpy
import riplpy.optical as optical
optical.load(riplpy.get_path())

# CSV / JSON / DataFrame exports for any optical database
optical.db.potentials.to_csv('potentials.csv')
optical.db.potentials.to_json('potentials.json')
df = optical.db.potentials.to_dataframe()
optical.db.deformations.to_csv('deformations.csv')
optical.db.references.to_json('references.json')

# Filter then export a tailored subset (e.g. all neutron potentials)
neutron = optical.db.potentials.filter_by_projectile('n')
neutron.to_csv('neutron_omps.csv')

# ATOMKI alpha-OMP per-nucleus profiles (each potential is a single
# AtomkiAlphaPotential object, not a Database, so use as_dict to export)
fe = optical.atomki.load_nucleus(26, 56)
import json
with open('atomki_fe56.json', 'w') as fp:
    json.dump({'Z': fe.Z, 'A': fe.A, 'r': fe.r, 'U_real': fe.U_real}, fp)
```

## Using with Main API

```python
import riplpy

riplpy.load()

# Convenience functions
pot = riplpy.get_omp(2405)
potentials = riplpy.find_omp('n', 82, 208, E=14.0)
irefs = riplpy.list_omps(projectile='n')
deform = riplpy.get_deformation(94, 239)
ref = riplpy.get_omp_reference(100)
```

## Physical Background

The optical model describes nuclear scattering using a complex potential:
```
U(r) = V(r) + iW(r) + V_so(r)·(l·s) + iW_so(r)·(l·s) + V_C(r)
```

**Components:**
- Real central potential V(r): Attractive, creates bound states
- Imaginary potential W(r): Accounts for absorption into compound nucleus
- Spin-orbit terms: Responsible for spin-dependent scattering
- Coulomb potential V_C: For charged projectiles

**Woods-Saxon Form:**
```
f(r) = 1 / (1 + exp((r - R)/a))
```
where R = r₀A^(1/3) and a is the diffuseness.

**Coupled-Channel:**
For deformed nuclei, the potential couples different angular momentum states, requiring a coupled-channel treatment with nuclear deformation parameters β₂, β₄, etc.
