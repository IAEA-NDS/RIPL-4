# Masses Section Release Checklist

This checklist tracks the status of the `riplpy.masses` section against the
RIPL-4 GitHub release layout (`/.../RIPL-4/github/`).

## Data files: GitHub release vs full release

The current RIPL-4 GitHub release ships a subset of the historical full
release. The table below shows which data files are present in `github/masses/`
and which are only available in the full release.

| File                        | GitHub release | Full release | Notes                                |
|-----------------------------|----------------|--------------|--------------------------------------|
| `mass-ame20.dat`            | yes            | yes          | AME 2020 atomic mass evaluation      |
| `mass-bskg3.dat`            | yes            | yes          | BSkG3 mass model                     |
| `mass-d1m.dat`              | yes            | yes          | D1M mass model                       |
| `mass-frdm12.dat`           | yes            | yes          | FRDM 2012                            |
| `mass-frdm95.dat`           | **no**         | yes          | FRDM 1995 -- skipped with a warning  |
| `mass-hfb14.dat`            | **no**         | yes          | HFB-14 -- skipped with a warning     |
| `mass-hfb27.dat`            | yes            | yes          | HFB-27 mass model                    |
| `mass-ws4.dat`              | yes            | yes          | WS4 mass model                       |
| `gs-deformations-exp.dat`   | yes            | yes          | Ground-state experimental deforms.   |
| `abundance.dat`             | yes            | yes          | Natural abundances                   |
| `Density-bskg3/z*.dat`      | yes (118 Z)    | yes          | Per-Z BSkG3 density tables           |
| `Density-d1m/z*.dat`        | yes (105 Z)    | yes          | Per-Z D1M density tables             |

## Mass model readers/writers

All eight mass-model modules are present and importable. The reader/writer
pairs are kept available even when their data files are missing so that users
with the full release continue to work.

| Module                    | DB key       | `MASS_MODELS` key(s)        | Reader/writer | Loaded by `load()` |
|---------------------------|--------------|-----------------------------|---------------|--------------------|
| `riplpy.masses.ame20`     | `ame20`      | `ame20`                     | yes           | yes                |
| `riplpy.masses.bskg3`     | `bskg3`      | `bskg3`                     | yes           | yes                |
| `riplpy.masses.d1m`       | `d1m`        | `d1m`                       | yes           | yes                |
| `riplpy.masses.frdm12`    | `frdm2012`   | `frdm12`, `frdm2012`        | yes           | yes                |
| `riplpy.masses.frdm95`    | `frdm1995`   | `frdm95`, `frdm1995`        | yes           | yes (skip if file missing) |
| `riplpy.masses.hfb14`     | `hfb14`      | `hfb14`                     | yes           | yes (skip if file missing) |
| `riplpy.masses.hfb27`     | `hfb27`      | `hfb27`                     | yes           | yes                |
| `riplpy.masses.ws4`       | `ws4`        | `ws4`                       | yes           | yes (added in this audit) |

## Auxiliary databases

| Module                          | DB key            | Source files                 | Status |
|---------------------------------|-------------------|------------------------------|--------|
| `riplpy.masses.ab`              | `natab`           | `abundance.dat`              | OK     |
| `riplpy.masses.deformations`    | `deformations`    | `gs-deformations-exp.dat`    | OK     |
| `riplpy.masses.density_bskg3`   | `density_bskg3`   | `Density-bskg3/z*.dat`       | OK     |
| `riplpy.masses.density_d1m`     | `density_d1m`     | `Density-d1m/z*.dat`         | OK     |

## Config paths

All `DATA_FILES` mass entries in `riplpy/config.py` resolve correctly against
the GitHub layout. Entries for `mass-frdm95.dat` and `mass-hfb14.dat` are kept
so the full-release path continues to work.

## Behaviour for missing files

`riplpy.masses.load_only_masses()` wraps each model load through
`_safe_load()`. When a model is marked optional (`required=False`) and its
data file is absent, `_safe_load`:

1. Catches `RiplFileNotFoundError`.
2. Logs a `WARNING` via the `riplpy.masses` logger.
3. Returns an empty `Database` instance so attribute access on
   `riplpy.masses.db.<model>` does not break downstream code
   (`in_ripl`, `list_databases`, etc.).

Currently `frdm95` and `hfb14` are marked optional; everything else is
required.

## Test verification

Run `python -m pytest tests/test_masses.py -v`:

- 20 passed
- 2 skipped (`TestFRDM95::test_load_from_directory`,
  `TestHFB14::test_load_from_directory`) -- skipped via `pytest.skip()` when
  the data file is absent.

The skipped tests will execute normally against a full RIPL-4 release.
