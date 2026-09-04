# Gamma Section

Giant Dipole Resonance (GDR) parameters and gamma-ray strength functions.

## Overview

The gamma section provides data for calculating gamma-ray emission in nuclear
reactions. RIPL-4 substantially restructured this section: the legacy
single-file products (theoretical GDR table, `gamma-strength-micro/`) were
replaced by re-derived experimental SLO/SMLO compilations and large
per-nucleus / per-Z microscopic strength-function directories.

The RIPL-4 GitHub release ships:

- **Recommended experimental GDR fits** - SLO and SMLO Lorentzian parameters
  (`experiment_slo`, `experiment_smlo`), with optional 1-sigma uncertainties
  (`experiment_slo_errors`, `experiment_smlo_errors`).
- **Experiment + systematics compilations** - ~8980 isotopes mixing
  experimental (`In=1`) and global systematics (`In=0`) values
  (`experiment_systematics_slo`, `experiment_systematics_smlo`).
- **Microscopic strength functions** - per-Z D1M+QRPA E1/M1 tables, per-nucleus
  SMLO E1 photoabsorption tables, per-Z SMLO M1 and TLO E1 tables, and the
  experimental PSFDatabase-v2024.1. These are "heavy" (thousands of files) and
  are loaded on demand (see [Heavy per-nucleus tables](#heavy-per-nucleus-tables)).

RIPL-3 legacy products (`theory_gdr`, MLO-only fits) are still exposed for
backwards compatibility but are **not** part of the GitHub release; see
[RIPL-3 legacy databases](#ripl-3-legacy-databases).

## Quick Start

```python
import riplpy
import riplpy.gamma as gamma
from riplpy import Nuclide

# Load the gamma section. include_heavy=False (default) skips the large
# per-nucleus directories so the load stays fast.
gamma.load(riplpy.get_path())

# Get recommended SLO experimental GDR parameters for Pb-208
pb208 = Nuclide(Z=82, A=208)
e = gamma.db.experiment_slo.get(pb208)
print(f"Pb-208 GDR (SLO): Er1 = {e.Er1:.2f} MeV, Wr1 = {e.Wr1:.2f} MeV")
print(f"  strength Sr1 = {e.Sr1:.3f} TRK, peak XS CSp1 = {e.CSp1:.1f} mb")
print(f"  fit range {e.Emin}-{e.Emax} MeV, ref {e.reference}")
```

## Available Databases

After `gamma.load(directory)` the following `gamma.db.*` databases are
populated. The Availability column indicates whether the underlying file is
part of the RIPL-4 GitHub release.

| Database | Attribute | Description | Availability |
|----------|-----------|-------------|--------------|
| Recommended SLO | `experiment_slo` | Recommended Single-Lorentzian experimental fits | RIPL-4 |
| Recommended SMLO | `experiment_smlo` | Recommended Simplified-Modified-Lorentzian fits | RIPL-4 |
| SLO + errors | `experiment_slo_errors` | SLO fits with 1-sigma uncertainties | RIPL-4 |
| SMLO + errors | `experiment_smlo_errors` | SMLO fits with 1-sigma uncertainties | RIPL-4 |
| Exp.+syst. SLO | `experiment_systematics_slo` | ~8980 isotopes, experiment (`In=1`) + systematics (`In=0`), SLO | RIPL-4 |
| Exp.+syst. SMLO | `experiment_systematics_smlo` | Same as above for SMLO | RIPL-4 |
| Exp. MLO (alias) | `experiment_mlo` | Backwards-compat alias of `experiment_smlo` | RIPL-3 legacy |
| Theory GDR | `theory_gdr` | Legacy theoretical (Goldhaber-Teller) systematics; empty on GitHub | RIPL-3 legacy |
| GSF (heavy) | `gsf` | Legacy tabulated GSF (D1M fallback); `None` unless `include_heavy=True` | RIPL-4 (heavy) |
| D1M+QRPA (heavy) | `gsf_d1m` | Per-Z D1M+QRPA E1/M1 PSF tables; `None` unless `include_heavy=True` | RIPL-4 (heavy) |
| SMLO E1 (heavy) | `smlo_e1` | Per-nucleus SMLO E1 photoabsorption; `None` unless `include_heavy=True` | RIPL-4 (heavy) |
| SMLO M1 (heavy) | `smlo_m1` | Per-Z SMLO M1 strength; `None` unless `include_heavy=True` | RIPL-4 (heavy) |
| TLO (heavy) | `tlo` | Per-Z TLO E1 strength; `None` unless `include_heavy=True` | RIPL-4 (heavy) |
| PSF (heavy) | `psf` | Experimental PSFDatabase-v2024.1; `None` unless `include_heavy=True` | RIPL-4 (heavy) |

`gamma.load(directory, include_heavy=False)` is the default. The five "heavy"
databases above are `None` after a normal load; access them via their direct
loaders (see [Heavy per-nucleus tables](#heavy-per-nucleus-tables)) or call
`gamma.load(directory, include_heavy=True)` to eagerly populate them (slow:
the SMLO E1 directory alone is ~8980 files).

## Examples

### Recommended experimental GDR parameters (SLO / SMLO)

```python
import riplpy
import riplpy.gamma as gamma
from riplpy import Nuclide

gamma.load(riplpy.get_path())

# Spherical nucleus: single GDR peak (Er2 is None)
sn116 = Nuclide(Z=50, A=116)
slo = gamma.db.experiment_slo.get(sn116)
smlo = gamma.db.experiment_smlo.get(sn116)

print("Sn-116 recommended GDR fits:")
print(f"  SLO : Er1 = {slo.Er1:.2f} MeV, Wr1 = {slo.Wr1:.2f} MeV, Sr1 = {slo.Sr1:.3f} TRK")
print(f"  SMLO: Er1 = {smlo.Er1:.2f} MeV, Wr1 = {smlo.Wr1:.2f} MeV, Sr1 = {smlo.Sr1:.3f} TRK")
```

### GDR for deformed nuclei (two-peak fits)

```python
import riplpy
import riplpy.gamma as gamma
from riplpy import Nuclide

gamma.load(riplpy.get_path())

# Deformed nuclei have a split GDR: the second component (Er2, Wr2) is filled
nd148 = Nuclide(Z=60, A=148)
e = gamma.db.experiment_slo.get(nd148)

print("Nd-148 GDR (SLO, deformed nucleus):")
print(f"  First peak : Er1 = {e.Er1:.2f} MeV, Wr1 = {e.Wr1:.2f} MeV")
print(f"  Second peak: Er2 = {e.Er2:.2f} MeV, Wr2 = {e.Wr2:.2f} MeV")
print(f"  Total strength Sr = {e.Sr:.3f} TRK, fit range {e.Emin}-{e.Emax} MeV")
```

### Fits with 1-sigma uncertainties

```python
import riplpy
import riplpy.gamma as gamma
from riplpy import Nuclide

gamma.load(riplpy.get_path())

pb208 = Nuclide(Z=82, A=208)
e = gamma.db.experiment_slo_errors.get(pb208)

print("Pb-208 SLO fit with uncertainties:")
print(f"  Er1 = {e.Er1:.2f} +/- {e.dEr1:.2f} MeV")
print(f"  Wr1 = {e.Wr1:.2f} +/- {e.dWr1:.2f} MeV")
print(f"  Sr1 = {e.Sr1:.3f} +/- {e.dSr1:.3f} TRK")
```

### Experiment + systematics compilation (~8980 isotopes)

```python
import riplpy
import riplpy.gamma as gamma
from riplpy import Nuclide

gamma.load(riplpy.get_path())

# This table mixes experimental (In=1) and systematics (In=0) values and uses
# the Base_GDR_Parameter fields E1/W1/E2/W2 plus Sr1/CSp1/Sr2/CSp2/Sr/In.
fe56 = Nuclide(Z=26, A=56)
s = gamma.db.experiment_systematics_smlo.get(fe56)

source = "experimental" if s.In == 1 else "systematics"
print(f"Fe-56 SMLO ({source}):")
print(f"  E1 = {s.E1:.3f} MeV, W1 = {s.W1:.2f} MeV")
print(f"  Sr1 = {s.Sr1:.3f} TRK, CSp1 = {s.CSp1:.1f} mb")

# Restrict the SLO table to the experimentally constrained subset
exp_only = gamma.db.experiment_systematics_slo.filter(lambda x: x.In == 1)
print(f"Experimentally constrained SLO entries: {len(exp_only)} of "
      f"{len(gamma.db.experiment_systematics_slo)}")
```

### Filter by GDR properties

```python
import riplpy
import riplpy.gamma as gamma

gamma.load(riplpy.get_path())

# Narrow GDRs in the recommended SLO compilation (guard None for safety)
narrow = gamma.db.experiment_slo.filter(
    lambda e: e.Wr1 is not None and e.Wr1 < 4.0
)
print(f"Recommended SLO fits with Wr1 < 4 MeV: {len(narrow)}")

# Nuclei with a resolved second (deformed) GDR component
deformed = gamma.db.experiment_slo.filter(lambda e: e.Er2 is not None)
print(f"Two-peak (deformed) SLO fits: {len(deformed)}")
for n in sorted(deformed, key=lambda x: (x.Z, x.A))[:5]:
    e = deformed.get(n)
    print(f"  {n.element_symbol}-{n.A}: Er1 = {e.Er1:.1f}, Er2 = {e.Er2:.1f} MeV")
```

### Standard Lorentzian strength function from a recommended fit

```python
import riplpy
import riplpy.gamma as gamma
from riplpy import Nuclide

gamma.load(riplpy.get_path())

def standard_lorentzian(E_gamma, CSp, E_r, Gamma_r):
    """Standard Lorentzian (SLO) photoabsorption cross section [mb].

    Args:
        E_gamma:  Gamma-ray energy [MeV]
        CSp:      Peak cross section [mb] (entry.CSp1)
        E_r:      Resonance energy [MeV] (entry.Er1)
        Gamma_r:  Resonance width [MeV] (entry.Wr1)
    """
    numer = CSp * (E_gamma * Gamma_r) ** 2
    denom = (E_gamma ** 2 - E_r ** 2) ** 2 + (E_gamma * Gamma_r) ** 2
    return numer / denom

fe54 = Nuclide(Z=26, A=54)
e = gamma.db.experiment_slo.get(fe54)

print("Fe-54 SLO photoabsorption cross section:")
print(f"{'E_gamma (MeV)':>15} {'sigma (mb)':>12}")
for E in [10, 15, e.Er1, 20, 25]:
    s = standard_lorentzian(E, e.CSp1, e.Er1, e.Wr1)
    print(f"{E:15.2f} {s:12.2f}")
```

## Heavy per-nucleus tables

The microscopic strength-function directories contain thousands of files, so
`gamma.load()` leaves `gamma.db.gsf / gsf_d1m / smlo_e1 / smlo_m1 / tlo / psf`
set to `None` by default. Access them on demand with the modules' own
`load_element(Z, directory)` / `load_nucleus(n, directory)` /
`load_category(name, directory)` helpers (preferred), or pass
`include_heavy=True` to `gamma.load()` to populate the `db.*` attributes (slow).

```python
import riplpy
import riplpy.gamma as gamma
from riplpy import Nuclide

d = riplpy.get_path()

# D1M+QRPA microscopic E1 PSF for a single element (Fe, Z=26)
fe_e1 = gamma.d1m.load_element(Z=26, directory=d, multipolarity='e1')
fe56 = Nuclide(Z=26, A=56)
pkt = fe_e1.get(fe56)
print(f"Fe-56 D1M E1: {len(pkt['U'])} energy points; "
      f"U[:3] = {pkt['U'][:3]} MeV, fE1[:3] = {pkt['fE1'][:3]} MeV^-3")
print(f"  temperature/excitation labels: {pkt['T'][:3]} ...")

# Per-nucleus SMLO E1 photoabsorption table (header carries Er1/Wr1/S1/beta)
o16 = Nuclide(Z=8, A=16)
o16_e1 = gamma.smlo_e1.load_nucleus(o16, directory=d).get(o16)
print(f"O-16 SMLO E1: Er1 = {o16_e1.get('Er1')} MeV, Wr1 = {o16_e1.get('Wr1')} "
      f"MeV, S1 = {o16_e1.get('S1')}, beta = {o16_e1.get('beta')}")

# Per-Z SMLO M1 strength table
fe_m1 = gamma.smlo_m1.load_element(Z=26, directory=d)
print(f"SMLO M1 nuclei for Z=26: {len(fe_m1)}")

# Per-Z TLO E1 strength table (deformation beta/gamma + EFF/INTER mode flag)
se = gamma.tlo.load_element(Z=34, directory=d)
se64 = Nuclide(Z=34, A=64)
tp = se.get(se64)
print(f"Se-64 TLO: beta = {tp['beta']}, gamma = {tp['gamma']}, "
      f"mode = {tp['mode']}, {len(tp['U'])} energy points")

# Experimental PSF database, one category at a time
psf = gamma.psf.load_category('photonuclear', directory=d)
print(f"PSF photonuclear datasets cover {len(psf)} nuclei")
```

Notes:

- `gamma.d1m.load_element` accepts `multipolarity='e1'` or `'m1'`.
- `gamma.smlo_e1` is per-nucleus (`load_nucleus`); `gamma.d1m`,
  `gamma.smlo_m1`, and `gamma.tlo` are per-Z (`load_element`).
- Each heavy module also exposes `load_all(directory)`, but loading every file
  is very slow; prefer the per-element / per-nucleus loaders.
- D1M / SMLO E1 packages are dicts/entries with `U` (photon-energy grid),
  `fE1` (cold-nucleus column), `T` (temperature labels), and `fE1_T` (full 2-D
  array). PSF datasets carry `E`, `f`, `rows`, `columns`, `source`,
  `filename`; multiple datasets per nucleus are stored as a list.

## RIPL-3 legacy databases

These products are part of the full RIPL distribution but are **not** shipped
in the RIPL-4 GitHub release. On the GitHub layout they load as empty databases
(with a logged warning) or resolve to a RIPL-4 alias.

- **`gamma.db.theory_gdr`** - the legacy theoretical Goldhaber-Teller GDR
  systematics (`gamma/gdr-parameters-theor.dat`). RIPL-4 replaces it with the
  per-nucleus D1M+QRPA predictions under `gamma/d1m/`. On the GitHub release
  this database has `len(gamma.db.theory_gdr) == 0`. Its entries (when the
  legacy file is present) expose `E1`, `W1`, `E2`, `W2`, and `eta`.
- **`gamma.db.experiment_mlo`** - RIPL-3 distributed MLO-only fits. RIPL-4
  dropped them in favour of SMLO, so `experiment_mlo` is now a backwards-compat
  **alias of `experiment_smlo`** (`gamma.db.experiment_mlo is
  gamma.db.experiment_smlo`).

```python
import riplpy
import riplpy.gamma as gamma

gamma.load(riplpy.get_path())

# Empty on the GitHub release (file not shipped)
print("theory_gdr entries:", len(gamma.db.theory_gdr))  # -> 0 on GitHub

# experiment_mlo is just an alias of experiment_smlo in RIPL-4
print("mlo is smlo:", gamma.db.experiment_mlo is gamma.db.experiment_smlo)
```

To use the legacy theoretical GDR table you must point RIPLpy at a full RIPL-3
distribution that still ships `gamma/gdr-parameters-theor.dat`.

## Entry Fields

### Recommended / errors experimental GDR entry (`experiment_slo`, `experiment_smlo`, `*_errors`)

| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `Er1`, `Wr1` | First component energy and width [MeV] |
| `Sr1` | First component strength [TRK units] |
| `CSp1` | First component Lorentzian peak cross section [mb] |
| `Er2`, `Wr2`, `Sr2`, `CSp2` | Second (deformed) component (None if spherical) |
| `Sr` | Total strength Sr1+Sr2 [TRK units] |
| `dEr1`, `dWr1`, `dSr1`, `dCSp1`, `dEr2`, `dWr2`, `dSr2`, `dCSp2`, `dSr` | 1-sigma uncertainties (populated only in `*_errors` databases) |
| `Emin`, `Emax` | Energy range of the fit [MeV] |
| `Id` | Identifier string (`'nat'` for natural-element rows) |
| `reference` | Reference key |

### Experiment + systematics entry (`experiment_systematics_slo`, `experiment_systematics_smlo`)

| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `E1`, `W1` | First component energy and width [MeV] |
| `E2`, `W2` | Second component energy and width [MeV] |
| `Sr1`, `CSp1` | First component strength [TRK] and peak XS [mb] |
| `Sr2`, `CSp2` | Second component strength [TRK] and peak XS [mb] |
| `Sr` | Total strength [TRK units] |
| `In` | Source flag: `1` = experimental, `0` = systematics |

### Theory GDR entry (RIPL-3 legacy `theory_gdr`)

| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `E1`, `W1` | First GDR peak energy and width [MeV] |
| `E2`, `W2` | Second GDR peak energy and width [MeV] (deformed) |
| `eta` | Deformation parameter |

### Heavy strength-function package (D1M / SMLO E1 / SMLO M1 / TLO)

| Field | Description |
|-------|-------------|
| `n` | Nuclide object |
| `U` | Photon-energy grid [MeV] |
| `fE1` | Cold-nucleus strength column [MeV^-3] |
| `T` | Temperature / excitation column labels |
| `fE1_T` | Full 2-D strength array (rows = energy, cols = `T`) |
| `Er1`, `Wr1`, `S1`, `beta` | SMLO E1 header Lorentzian parameters (SMLO E1 only) |
| `beta`, `gamma`, `mode` | TLO deformation and `EFF`/`INTER` mode flag (TLO only) |

## Writing & Exporting

All gamma databases (scalar SLO/SMLO fits and heavy strength-function
tables) support the standard export interface:

```python
import riplpy
import riplpy.gamma as gamma
gamma.load(riplpy.get_path())

# Scalar GDR fit databases
gamma.db.experiment_slo.save('exp_slo.dat')              # native ASCII
gamma.db.experiment_slo.to_csv('exp_slo.csv')
gamma.db.experiment_slo.to_json('exp_slo.json')
df = gamma.db.experiment_slo.to_dataframe()

# Experiment-only subset (drop the systematics rows)
exp_only = gamma.db.experiment_systematics_slo.filter(lambda x: x.In == 1)
exp_only.to_csv('slo_exp_only.csv')

# Heavy per-element spectral tables (load on demand)
fe = gamma.d1m.load_element(Z=26, multipolarity='e1')
fe.to_dataframe().to_csv('d1m_fe_e1.csv', index=False)   # U / fE1 / fE1_T

# The legacy GSF database can also lay itself out on disk as RIPL did:
gamma.load(riplpy.get_path(), include_heavy=True)
gamma.db.gsf.save_all('out/gsf_by_Z')                    # per-element files
gamma.db.gsf.save_entries_by_Z('out/gsf_by_nuc')         # per-nucleus files
```

## Using with Main API

```python
import riplpy

riplpy.load()

# NOTE: riplpy.get_gdr() / get_gdrs() use the RIPL-3 LEGACY theory_gdr table,
# which is EMPTY on the GitHub release (no gdr-parameters-theor.dat). On the
# GitHub layout these lookups raise / skip. Prefer the RIPL-4 databases:
import riplpy.gamma as gamma
from riplpy import Nuclide

pb208 = Nuclide(Z=82, A=208)
e = gamma.db.experiment_slo.get(pb208)
print(f"Pb-208 (SLO): Er1 = {e.Er1} MeV, Wr1 = {e.Wr1} MeV")
```

## Physical Background

The Giant Dipole Resonance (GDR) is a collective nuclear excitation where
protons oscillate against neutrons. It dominates the photoabsorption cross
section at energies around 10-25 MeV.

**Systematics:**
- Spherical nuclei: single peak near E ~ 80/A^(1/3) MeV
- Deformed nuclei: split into two peaks along the major/minor axes
- Width: Gamma ~ 4-6 MeV (spreading width + escape width)

**Gamma Strength Function:**
The E1 strength function f_E1(E_gamma) is related to the photoabsorption cross
section by:
```
f_E1(E_gamma) = sigma_abs(E_gamma) / (3 pi^2 hbar^2 c^2 E_gamma)
```
