# Fission Section

Fission barrier heights and parameters for actinide and transactinide nuclei.

## Overview

The fission section provides fission barrier data from both experimental compilations and theoretical calculations. Barriers are crucial for calculating fission cross sections and fission probabilities.

**RIPL-4 Models (GitHub release + full distribution):**
- **Empirical** - RIPL-4 empirical barrier heights (`empirical-barriers-ripl4.dat`)
- **Empirical (New)** - Alias of the RIPL-4 empirical compilation
- **EMPIRE** - EMPIRE code evaluated barriers
- **BSkG3** - Brussels-Skyrme theoretical barriers (inner/outer/isomer sections)
- **D1M** - Gogny D1M HFB low-energy fission path barriers (`barriers-d1m_lep.dat`)

**RIPL-3 Legacy Models (full RIPL distribution only — NOT in the RIPL-4 GitHub release):**
- **HFB** - Hartree-Fock-Bogoliubov theoretical barriers (`empirical-hfb-barriers.dat`)

The HFB legacy database is skipped (returns empty) with a warning by
`fission.load()` when its data file is absent; its reader/writer remain
available for users with the full RIPL distribution.

## Quick Start

```python
import riplpy
import riplpy.fission as fission
from riplpy import Nuclide

# Load fission data
fission.load(riplpy.get_path())

# Get fission barrier for U-235
u235 = Nuclide(Z=92, A=235)
barrier = fission.db.empirical_barriers.get(u235)
print(f"U-235 inner barrier: {barrier.Ea:.2f} MeV")
print(f"U-235 outer barrier: {barrier.Eb:.2f} MeV")
```

## Available Models

| Model | Attribute | Type | Availability |
|-------|-----------|------|--------------|
| Empirical | `empirical_barriers` | Experimental | RIPL-4 |
| Empirical (New) | `empirical_barriers_new` | Experimental | RIPL-4 (alias of Empirical) |
| EMPIRE | `empire_barriers` | Experimental | RIPL-4 |
| BSkG3 | `bskg3_barriers` | Theoretical | RIPL-4 |
| D1M | `d1m_barriers` | Theoretical | RIPL-4 |
| HFB | `hfb_barriers` | Theoretical | RIPL-3 legacy |

## Examples

### Basic Barrier Lookup

```python
import riplpy
import riplpy.fission as fission
from riplpy import Nuclide

fission.load(riplpy.get_path())

# Get empirical barrier for Pu-239
pu239 = Nuclide(Z=94, A=239)
entry = fission.db.empirical_barriers.get(pu239)

# Actinides have double-humped barriers
print(f"Pu-239 double-humped barrier:")
print(f"  Inner barrier (Ea): {entry.Ea:.2f} MeV")
print(f"  Outer barrier (Eb): {entry.Eb:.2f} MeV")
print(f"  Inner curvature (hwa): {entry.hwa:.2f} MeV")
print(f"  Outer curvature (hwb): {entry.hwb:.2f} MeV")
```

### Compare Barrier Models

```python
import riplpy
import riplpy.fission as fission
from riplpy import Nuclide

fission.load(riplpy.get_path())

# Compare models for U isotopes (inner barrier heights)
print("Uranium fission barriers - Inner barrier (MeV):")
print(f"{'Isotope':>10} {'Empirical':>12} {'BSkG3':>12} {'HFB(legacy)':>12}")

for A in range(232, 240):
    u = Nuclide(Z=92, A=A)
    values = []

    # Empirical barriers (Ea attribute)
    try:
        entry = fission.db.empirical_barriers.get(u)
        values.append(f"{entry.Ea:12.2f}")
    except KeyError:
        values.append("         N/A")

    # BSkG3 barriers (dict-based inner section with key 'E[MeV]')
    try:
        entry = fission.db.bskg3_barriers.get(u)
        values.append(f"{entry.inner['E[MeV]']:12.2f}")
    except KeyError:
        values.append("         N/A")

    # HFB barriers — RIPL-3 legacy, may be None
    hfb_db = fission.db.hfb_barriers
    if hfb_db is None:
        values.append("    (legacy)")
    else:
        try:
            entry = hfb_db.get(u)
            values.append(f"{entry.Bin:12.2f}")
        except KeyError:
            values.append("         N/A")

    print(f"U-{A:3d}:   " + " ".join(values))
```

### Double-Humped Barriers

```python
import riplpy
import riplpy.fission as fission
from riplpy import Nuclide

fission.load(riplpy.get_path())

# Actinides typically have double-humped barriers
actinides = [
    (92, 235, "U-235"),
    (92, 238, "U-238"),
    (94, 239, "Pu-239"),
    (94, 241, "Pu-241"),
]

print("Double-humped fission barriers (empirical):")
print(f"{'Nucleus':>10} {'Ea (MeV)':>10} {'hwa':>8} {'Eb (MeV)':>10} {'hwb':>8} {'sym_a':>6} {'sym_b':>6}")

for Z, A, name in actinides:
    n = Nuclide(Z=Z, A=A)
    try:
        entry = fission.db.empirical_barriers.get(n)
        print(f"{name:>10} {entry.Ea:10.2f} {entry.hwa:8.2f} {entry.Eb:10.2f} {entry.hwb:8.2f} {entry.syma:>6} {entry.symb:>6}")
    except KeyError:
        print(f"{name:>10} Not available")
```

### Barrier Curvatures

```python
import riplpy
import riplpy.fission as fission
from riplpy import Nuclide

fission.load(riplpy.get_path())

# Barrier curvature determines the transmission probability shape
u238 = Nuclide(Z=92, A=238)
entry = fission.db.empirical_barriers.get(u238)

print("U-238 barrier parameters:")
print(f"  Inner barrier (A):")
print(f"    Height (Ea): {entry.Ea:.2f} MeV")
print(f"    Curvature (hwa): {entry.hwa:.2f} MeV")
print(f"    Symmetry: {entry.syma}")
print(f"  Outer barrier (B):")
print(f"    Height (Eb): {entry.Eb:.2f} MeV")
print(f"    Curvature (hwb): {entry.hwb:.2f} MeV")
print(f"    Symmetry: {entry.symb}")
```

### Calculate Barrier Penetrability

```python
import riplpy
import riplpy.fission as fission
from riplpy import Nuclide
import math

fission.load(riplpy.get_path())

def hill_wheeler_penetrability(E, E_B, hw):
    """
    Calculate fission barrier penetrability using Hill-Wheeler formula.

    Args:
        E: Excitation energy [MeV]
        E_B: Barrier height [MeV]
        hw: Barrier curvature [MeV]

    Returns:
        Transmission coefficient T_f
    """
    if hw <= 0:
        return 1.0 if E > E_B else 0.0

    x = 2 * math.pi * (E - E_B) / hw
    return 1.0 / (1.0 + math.exp(-x))

# Get U-238 inner barrier
u238 = Nuclide(Z=92, A=238)
entry = fission.db.empirical_barriers.get(u238)

E_B = entry.Ea  # Inner barrier height
hw = entry.hwa   # Inner barrier curvature

print(f"U-238 fission penetrability (inner barrier Ea = {E_B:.2f} MeV):")
print(f"{'E (MeV)':>10} {'T_f':>12}")
for E in [4.0, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0]:
    T = hill_wheeler_penetrability(E, E_B, hw)
    print(f"{E:10.1f} {T:12.4f}")
```

### Filter by Barrier Height

```python
import riplpy
import riplpy.fission as fission

fission.load(riplpy.get_path())

# Find nuclei with relatively low fission barriers
# Empirical barriers use Ea (inner) / Eb (outer) — see Entry Fields below.
low_barrier = fission.db.empirical_barriers.filter(
    lambda e: e.Ea is not None and e.Ea < 5.0
)
print(f"Nuclei with empirical Ea < 5.0 MeV: {len(low_barrier.data)}")

# List them
for n, entry in list(low_barrier.data.items())[:10]:
    print(f"  {n.sym}-{n.A}: Ea = {entry.Ea:.2f} MeV")
```

## Entry Fields

### Empirical Barrier Entry
| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `Ea` | Inner barrier height [MeV] |
| `hwa` | Inner barrier curvature [MeV] |
| `syma` | Inner saddle point symmetry |
| `Eb` | Outer barrier height [MeV] |
| `hwb` | Outer barrier curvature [MeV] |
| `symb` | Outer saddle point symmetry |
| `deltaf` | Fission splitting parameter [MeV] |

### HFB Barrier Entry (RIPL-3 legacy)
| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `Bin` | Inner barrier height [MeV] |
| `hwin` | Inner barrier curvature [MeV] |
| `Bout` | Outer barrier height [MeV] |
| `hwout` | Outer barrier curvature [MeV] |
| `Bout2` | Third barrier (if present) [MeV] |

### BSkG3 Barrier Entry (RIPL-4)
| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `inner` | Inner barrier `{E[MeV], B20, B22, B30}` |
| `outer1` | First outer barrier `{E[MeV], B20, B22, B30}` |
| `outer2` | Second outer barrier `{E[MeV], B20, B22, B30}` |
| `isomer` | First shape isomer `{E[MeV], B20, B22, B30}` |

### Saddle Point Symmetries
| Code | Meaning |
|------|---------|
| `S` | Symmetric |
| `GA` | Axially asymmetric (triaxial) |
| `MA` | Mass asymmetric |

## Writing & Exporting

Each fission database (empirical, EMPIRE, BSkG3, D1M, HFB) supports the
standard export interface plus a native ASCII writer:

```python
import riplpy
import riplpy.fission as fission
fission.load(riplpy.get_path())

# Native ASCII round-trip
fission.db.bskg3_barriers.save('bskg3_copy.dat')

# Machine-readable exports (DataFrame flattens nested
# inner/outer1/outer2/isomer sections to dotted columns)
fission.db.bskg3_barriers.to_csv('bskg3.csv')
fission.db.bskg3_barriers.to_json('bskg3.json')
df = fission.db.bskg3_barriers.to_dataframe()

# Filter then save (still readable by the same loader)
high_inner = fission.db.bskg3_barriers.filter(
    lambda e: e.inner['E[MeV]'] > 7.0
)
high_inner.save('bskg3_high.dat')
```

## Using with Main API

```python
import riplpy

riplpy.load()

# Get fission barrier
barrier = riplpy.get_fission_barrier(92, 235, model='empirical')

# Batch operations
actinides = [(92, 235), (92, 238), (94, 239)]
barriers = riplpy.get_fission_barriers(actinides, model='hfb', skip_missing=True)
```

## Physical Background

Fission barriers arise from the competition between the attractive nuclear force and the repulsive Coulomb force as a nucleus deforms toward scission. Heavy nuclei (Z ≥ 90) typically have:

- **Inner barrier (A)**: First saddle point, often axially symmetric
- **Second minimum**: Shape isomer state
- **Outer barrier (B)**: Second saddle point, often mass-asymmetric

The barrier transmission probability is approximately:
```
T_f ≈ 1 / (1 + exp(-2π(E - E_B)/ℏω))
```
where E is the excitation energy, E_B is the barrier height, and ℏω is the curvature parameter.
