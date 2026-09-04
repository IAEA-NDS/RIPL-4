# Gamma Section - Release Checklist

Audit of the `riplpy.gamma` section against the public RIPL-4 github release
at `/path/to/RIPL-4/github/`.

The gamma section was substantially restructured in RIPL-4: legacy single-file
products (e.g. `gdr-parameters-theor.dat`, `gamma/gamma-strength-micro/`) were
replaced by larger per-Z or per-nucleus directories, and the experimental SLO
and SMLO compilations were re-derived.

## Data files

| File / directory                                                          | Reader                                          | Status            | Notes                                                                                                  |
| ------------------------------------------------------------------------- | ----------------------------------------------- | ----------------- | ------------------------------------------------------------------------------------------------------ |
| `gamma/gdr-parameters-theor.dat`                                          | `gdr.load`                                      | missing in github | Legacy theoretical GDR file. Loader skips with a warning and returns an empty `Database`.              |
| `gamma/gamma-strength-micro/`                                             | `gsf.load_all`                                  | missing in github | Legacy GSF directory. Loader falls back to D1M+QRPA per-Z tables to preserve `U`/`fE1` schema.         |
| `gamma/gdr_parameters_exp_new/gdr_parameters_recommended_exp_slo.dat`     | `exp.load_slo` (default)                        | shipped, NEW      | Single-line records; recommended SLO GDR fits.                                                         |
| `gamma/gdr_parameters_exp_new/gdr_parameters_recommended_exp_smlo.dat`    | `exp.load_smlo` (default), `exp.load_mlo` alias | shipped, NEW      | Single-line records; recommended SMLO GDR fits. Aliased as `mlo` for back-compat.                      |
| `gamma/gdr_parameters_exp_new/gdr_parameters&errors_exp_slo.dat`          | `exp.load_slo(..., errors=True)`                | shipped, NEW      | Two-line records (values + 1-sigma uncertainties); SLO.                                                |
| `gamma/gdr_parameters_exp_new/gdr_parameters&errors_exp_smlo.dat`         | `exp.load_smlo(..., errors=True)`               | shipped, NEW      | Two-line records (values + 1-sigma uncertainties); SMLO.                                               |
| `gamma/gdr_parameters_exp&systematics/gdr-parameters_exp&systematics_slo.dat`  | `systematics.load_slo`                     | shipped, NEW      | ~8980 isotopes; mix of experimental (`In=1`) and systematics (`In=0`) values, SLO.                     |
| `gamma/gdr_parameters_exp&systematics/gdr-parameters_exp&systematics_smlo.dat` | `systematics.load_smlo`                    | shipped, NEW      | Same as above for SMLO.                                                                                |
| `gamma/d1m/z<NNN>_e1`, `z<NNN>_m1`                                        | `d1m.load_all` / `d1m.load_element`             | shipped, NEW      | D1M+QRPA microscopic E1 and M1 PSFs. Per-Z files with many nuclei concatenated.                        |
| `gamma/smlo_E1/fe1_the_<Z>_<A>_photoabs_h_SMLO.dat`                       | `smlo_e1.load_all` / `smlo_e1.load_nucleus`     | shipped, NEW      | ~8980 per-nucleus SMLO E1 photoabsorption tables.                                                      |
| `gamma/smlo_M1/z<NNN>_m1`                                                 | `smlo_m1.load_all` / `smlo_m1.load_element`     | shipped, NEW      | Per-Z SMLO M1 strength tables. Reuses the D1M multi-nucleus reader.                                    |
| `gamma/tlo/FE1_z<NN>.dat` (+ `5dch.txt`)                                  | `tlo.load_all` / `tlo.load_element`             | shipped, NEW      | Per-Z TLO E1 strength tables (CRLF line endings).                                                      |
| `gamma/PSFDatabase-v2024.1/<category>/`                                   | `psf.load_all` / `psf.load_category`            | shipped, PARTIAL  | Experimental PSF compilation across 8 categories (arcdrc, nrf, oslo, pg, photonuclear, pp, RM, thc).   |

## Readers / modules

| Module                              | Status   | Notes                                                                                              |
| ----------------------------------- | -------- | -------------------------------------------------------------------------------------------------- |
| `riplpy/gamma/core.py`              | updated  | `Experimental_GDR_Parameter` gains `CSp1, dCSp1, CSp2, dCSp2`. New `Systematics_GDR_Parameter`.    |
| `riplpy/gamma/gdr.py`               | updated  | `Database` now subclasses `NuclideDatabase`. `load()` returns empty + warning if file absent.      |
| `riplpy/gamma/gsf.py`               | updated  | Falls back to D1M tables when `gamma-strength-micro/` is absent.                                   |
| `riplpy/gamma/exp.py`               | rewritten | Parses the new RIPL-4 SLO/SMLO recommended (1-line) and errors (2-line) files. MLO is alias.       |
| `riplpy/gamma/systematics.py`       | NEW      | Parses the experiment+systematics SLO/SMLO tables (`(2I4, 9F9.3, I5)`).                            |
| `riplpy/gamma/d1m.py`               | NEW      | Multi-nucleus per-Z reader for D1M+QRPA E1 and M1 tables.                                          |
| `riplpy/gamma/smlo_e1.py`           | NEW      | Per-nucleus SMLO E1 photoabsorption reader (header carries `Er1, Wr1, S1, beta`).                  |
| `riplpy/gamma/smlo_m1.py`           | NEW      | Per-Z SMLO M1 reader; reuses `d1m.read_ascii_file`.                                                |
| `riplpy/gamma/tlo.py`               | NEW      | Per-Z TLO E1 reader; parses per-nucleus `bet`, `gam`, and `EFF`/`INTER` mode flag.                 |
| `riplpy/gamma/psf.py`               | NEW      | Experimental PSF database reader; supports all 8 PSFDatabase-v2024.1 categories.                   |
| `riplpy/gamma/__init__.py`          | updated  | `_safe_load` wraps each loader; all new databases exposed under `db.*`.                            |

## Database attributes exposed by `gamma.load(directory)`

| Attribute                          | Source                                                             |
| ---------------------------------- | ------------------------------------------------------------------ |
| `db.gsf`                           | `gamma-strength-micro/` (legacy) -> D1M fallback                   |
| `db.theory_gdr`                    | `gdr-parameters-theor.dat` (legacy, often empty in github)         |
| `db.experiment_slo`                | `gdr_parameters_recommended_exp_slo.dat`                           |
| `db.experiment_smlo`               | `gdr_parameters_recommended_exp_smlo.dat`                          |
| `db.experiment_mlo`                | alias of `experiment_smlo` (MLO == SMLO in RIPL-4)                 |
| `db.experiment_slo_errors`         | `gdr_parameters&errors_exp_slo.dat`                                |
| `db.experiment_smlo_errors`        | `gdr_parameters&errors_exp_smlo.dat`                               |
| `db.experiment_systematics_slo`    | `gdr-parameters_exp&systematics_slo.dat`                           |
| `db.experiment_systematics_smlo`   | `gdr-parameters_exp&systematics_smlo.dat`                          |
| `db.gsf_d1m`                       | `gamma/d1m/`                                                       |
| `db.smlo_e1`                       | `gamma/smlo_E1/`                                                   |
| `db.smlo_m1`                       | `gamma/smlo_M1/`                                                   |
| `db.tlo`                           | `gamma/tlo/`                                                       |
| `db.psf`                           | `gamma/PSFDatabase-v2024.1/` (multiple datasets per nucleus)       |

## Config keys (`riplpy/config.py`)

Added under the gamma section:

- `gdr_parameters_recommended_exp_slo` / `gdr_parameters_recommended_exp_smlo`
- `gdr_parameters_errors_exp_slo` / `gdr_parameters_errors_exp_smlo`
- `gdr_parameters_systematics_slo` / `gdr_parameters_systematics_smlo`
- `gsf_d1m_dir`, `gsf_smlo_e1_dir`, `gsf_smlo_m1_dir`, `gsf_tlo_dir`, `psf_database_dir`

Legacy keys (`gsf_data`, `gdr_parameters_exp`, `gdr_parameters_theor`) are
retained for backwards compatibility; loaders are tolerant to missing files.

## Tests (`tests/test_gamma.py`)

- Legacy GDR-theoretical tests now skip cleanly when the file is absent.
- The GSF tests pass against the D1M fallback (the U / fE1 schema is preserved).
- New test classes cover every new database:
  - `TestExperimentalGDR` exercises both SLO and SMLO loaders, the
    `errors=True` path, and the MLO -> SMLO alias.
  - `TestSystematicsGDR` checks the ~8980-row experiment+systematics tables and
    validates the `In` flag is in `{0, 1}`.
  - `TestD1M`, `TestSMLO_E1`, `TestSMLO_M1`, `TestTLO` exercise the new
    per-Z or per-nucleus PSF readers.
  - `TestPSFDatabase` smoke-tests the partial PSF database reader.
  - `TestGammaSectionLoad::test_section_load_does_not_raise` verifies that
    `gamma.load(directory)` populates every advertised `db.*` attribute.
- Result against the github release: **19 passed, 4 skipped** (the 4 skips are
  the legacy theoretical-GDR tests, which require a file not shipped in the
  github layout).

## Known limitations / follow-ups

- **PSFDatabase-v2024.1 partial coverage.** The reader handles the dominant
  `# Col n:` header layout used across all 8 categories, recognises a strength
  column heuristically, and stores raw rows. Specialised treatment of the
  per-experiment metadata (analysis flags, reference resolution, multi-J
  primary gamma decompositions, etc.) is deferred. The associated
  `*.readme` files are not parsed.
- **D1M load times.** Loading every per-Z file under `gamma/d1m/`, every
  per-nucleus file under `gamma/smlo_E1/`, etc. is slow at startup. The
  per-element / per-nucleus convenience loaders allow lazy access. A future
  refinement could memoise these directories.
- **Legacy theoretical GDR fits.** When the legacy
  `gdr-parameters-theor.dat` is absent, `db.theory_gdr` is empty. The
  D1M+QRPA per-nucleus predictions live in `db.gsf_d1m` / `db.smlo_e1` and
  can supply replacement GDR-like Lorentzian parameters via the SMLO E1
  header (which carries `Er1, Wr1, S1, beta` for each nucleus).
- **Natural-element entries.** Rows with `A=0` (natural isotopic composition)
  are remapped to `Nuclide(Z, A=2Z)` with `Id='nat'` so they can be stored in
  the hashable Nuclide-keyed dict. Users should filter by `entry.Id == 'nat'`
  to retrieve them.
