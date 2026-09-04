# Optical Section Release Checklist

This document tracks the release-readiness of the `riplpy.optical` section
against the RIPL-4 github data release living at
`/path/to/RIPL-4/github/optical/`.

## GitHub data layout (`optical/`)

```
optical/
├── atomki/                          # 4359 TALYS alpha-OMP gnu output files
│   └── z{Z:03d}a{A:03d}a_talys_alphaomp9real.gnu
├── code/                            # Fortran retrieval source + utilities
├── om-data/                         # Primary OMP data directory (loaded)
│   ├── mod-potentials/              # Updated/corrected OMPs (omp-XXXXXX.dat)
│   ├── om-deformations.dat          # Excited-level deformations
│   ├── om-deformations.readme
│   ├── om-index-by-Z.txt
│   ├── om-index.txt                 # 581 index entries
│   ├── om-parameter-u.dat           # Canonical 584-potential library
│   ├── om-parameter-u-032013.dat    # Dated snapshot (not loaded)
│   ├── om-parameter-u-062012.dat    # Dated snapshot (not loaded)
│   ├── om-parameter-u-062014.dat    # Dated snapshot (not loaded)
│   ├── om-parameter-u-112013.dat    # Dated snapshot (not loaded)
│   ├── om-parameter-u-112015.dat    # Dated snapshot (not loaded)
│   ├── om-parameter-u-May2014.dat   # Dated snapshot (not loaded)
│   ├── om-parameter-u.readme
│   ├── om-references.txt            # 143 bibliographic references
│   └── ROP2013za.dat                # Avrigeanu revised alpha ROP table
├── om-summary-2026/                 # Summaries, tooling, source code
├── EPJ239...pdf, JPG48...pdf, ...   # Reference PDFs (documentation)
```

## Coverage matrix

| File / directory                                | Reader                                | Status     |
| ----------------------------------------------- | ------------------------------------- | ---------- |
| `om-data/om-parameter-u.dat`                    | `optical.omp.load()`                  | ✅ 584/584 |
| `om-data/om-index.txt`                          | `optical.index.load()`                | ✅ 581/581 |
| `om-data/om-index-by-Z.txt`                     | (same content as index, not parsed)   | ➖ alias   |
| `om-data/om-references.txt`                     | `optical.references.load()`           | ✅ 143/143 |
| `om-data/om-deformations.dat`                   | `optical.deformations.load()`         | ✅ 1423/1423 |
| `om-data/mod-potentials/`                       | `optical.load_modified_potentials()`  | ✅ 7/7     |
| `om-data/ROP2013za.dat`                         | `optical.rop2013.load()` (NEW)        | ✅ 79/79   |
| `atomki/*.gnu` (4359 files)                     | `optical.atomki.load_nucleus(Z, A)`   | ✅ lazy    |

## Deferred / out-of-scope

| Item                                            | Reason                                            |
| ----------------------------------------------- | ------------------------------------------------- |
| `om-parameter-u-{date}.dat` snapshots (6 files) | Dated snapshots of the canonical file; the loader exposes the canonical `om-parameter-u.dat` only |
| `om-data/om-parameter-u.readme`                 | Documentation only; not parsed                    |
| `om-data/om-deformations.readme`                | Documentation only; not parsed                    |
| `code/` (Fortran sources, msg files, original WLH_original.py, ext_functions.f) | Original retrieval / utility source; not a Python concern |
| `om-summary-2026/` (compiled artefacts, source, `lst/`, `ncb/`, `references/`, `summary/`) | Tooling and pre-computed summaries; not parsed by riplpy |
| `EPJ239…pdf`, `JPG48…pdf`, `NDS173…pdf`, `PRC107…pdf` | Reference papers (documentation only)         |

## API surface

After `riplpy.optical.load(directory=...)` the following attributes are
populated on `optical.db`:

* `optical.db.index` — `OMPIndex` (581 entries)
* `optical.db.potentials` — `OMPDatabase` (584 potentials, spherical + CC)
* `optical.db.deformations` — `DeformationDatabase` (1423 entries)
* `optical.db.references` — `ReferenceDatabase` (143 references)
* `optical.db.rop2013` — `ROP2013Database` (79 nuclei, iref 9999)

Each of these loaders is wrapped in a `try/except FileNotFoundError` so
that a partial RIPL install only emits a `riplpy.logger.warning(...)`
rather than aborting the section. The ATOMKI directory is **not**
preloaded (it would parse ~4359 files); callers use
`optical.atomki.load_nucleus(Z, A)` on demand.

Auxiliary functions:

* `optical.load_modified_potentials(directory=None)`
* `optical.load_with_modifications(directory=None, apply_modifications=True)`
* `optical.atomki.load_nucleus(Z, A, directory=None)`
* `optical.atomki.available_nuclei(directory=None)`
* `optical.atomki.filename_for(Z, A)`
* `optical.rop2013.load(directory=None)` / `optical.rop2013.read(path)`

Top-level convenience helpers in `riplpy/__init__.py`:

* `riplpy.get_omp(iref)`
* `riplpy.list_omps(projectile=None)`
* `riplpy.find_omp(projectile, Z, A, E)`
* `riplpy.get_deformation(Z, A, Ex=None, L=None)`
* `riplpy.get_omp_reference(ref_num)`

## Test status

`pytest tests/test_optical.py` — **76 passing**:

* TestFortranFloatParser (4)
* TestOMPIndex (8)
* TestSphericalOMP (5)
* TestCoupledChannelOMP (6)
* TestOMPDatabase (14)
* TestOpticalModuleAPI (5)
* TestDeformations (6)
* TestReferences (6)
* TestDataIntegrity (6)
* TestModifiedPotentials (6)
* TestROP2013 (4)  ← new
* TestAtomki (4)   ← new

## Release notes

* `riplpy/optical/config.py` now lists `rop2013` and `atomki` paths in
  `FILE_PATHS` so callers can locate the auxiliary data files relative
  to the configured RIPL directory.
* The optical section now degrades gracefully when individual data files
  are missing rather than raising during `riplpy.load()`.
