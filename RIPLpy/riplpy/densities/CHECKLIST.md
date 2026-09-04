# Densities Section — Release Checklist

Audit of the `riplpy.densities` section against the public RIPL-4 github
release at `/path/to/RIPL-4/github/`.

## Data files

| File / directory                                                  | Reader                                | Status              | Notes                                                                                       |
| ----------------------------------------------------------------- | ------------------------------------- | ------------------- | ------------------------------------------------------------------------------------------- |
| `densities/shellcor-ms.dat`                                       | `shell_corr.load_ms`                  | shipped             | Myers–Swiatecki shell corrections.                                                          |
| `densities/shellmengoninakajima.dat`                              | `shell_corr.load_mk`                  | missing in github   | Mengoni–Nakajima shell corrections. Loader skips with a warning when absent.                |
| `densities/total/level-densities-egsm.dat`                        | `egsm.load`                           | shipped             | Path corrected to live under `total/`.                                                       |
| `densities/total/level-densities-egsm-norm.dat`                   | `egsm_norm.load`                      | shipped             | Path corrected to live under `total/`.                                                       |
| `densities/level-densities-bfmeff.dat`                            | `bsfg.load`                           | missing in github   | Legacy file from older release. Config key retained; loader skips with warning.             |
| `densities/level-densities-ctmeff.dat`                            | `ct.load`                             | missing in github   | Legacy file from older release. Config key retained; loader skips with warning.             |
| `densities/total/level-densities-hfb/zXXX.tab`                    | `hfb.load_all` / `hfb.load_element`   | missing in github   | Per-Z HFB spin-dependent level densities. Loader skips with warning.                        |
| `densities/total/bsk14-comb/zXXX.tab` (+ `.ld`)                   | `bsk14_comb.load_*`                   | shipped, NEW reader | BSk14 + combinatorial. Reuses HFB `.tab` parser.                                            |
| `densities/total/bskg3-comb/zXXX.tab` (+ `.ld`)                   | `bskg3_comb.load_*`                   | shipped, NEW reader | BSkG3 + combinatorial. Reuses HFB `.tab` parser.                                            |
| `densities/total/qrpabe/zXXX.tab`                                 | `qrpabe.load_*`                       | shipped, NEW reader | QRPA-BE. Reuses HFB `.tab` parser.                                                          |
| `densities/total/thfb-comb/zXXX.tab` (+ `.ld`)                    | `thfb_comb.load_*`                    | shipped, NEW reader | T-HFB + combinatorial. Reuses HFB `.tab` parser.                                            |

## Readers / modules

| Module                                             | Status   | Notes                                                                                                          |
| -------------------------------------------------- | -------- | -------------------------------------------------------------------------------------------------------------- |
| `riplpy/densities/bsfg.py`                         | unchanged | Backing file missing in github release; section loader skips it with `RiplFileNotFoundError`.                  |
| `riplpy/densities/ct.py`                           | unchanged | Same as BSFG: legacy-only.                                                                                     |
| `riplpy/densities/egsm.py`                         | unchanged | Now resolves to `densities/total/...` via fixed `config.py`.                                                   |
| `riplpy/densities/egsm_norm.py`                    | unchanged | Now resolves to `densities/total/...` via fixed `config.py`.                                                   |
| `riplpy/densities/shell_corr.py`                   | unchanged | MK loader raises `FileNotFoundError` when file absent; section loader catches it.                              |
| `riplpy/densities/hfb.py`                          | unchanged | Already tolerant of missing directory (warns and returns empty database).                                      |
| `riplpy/densities/_comb_tab.py`                    | NEW       | Internal shared `.tab` directory loader used by the four combinatorial models.                                 |
| `riplpy/densities/bsk14_comb.py`                   | NEW       | Provides `load`, `load_all`, `load_element`, plus a `Database` class.                                          |
| `riplpy/densities/bskg3_comb.py`                   | NEW       | Same API as `bsk14_comb`.                                                                                      |
| `riplpy/densities/qrpabe.py`                       | NEW       | Same API as `bsk14_comb`.                                                                                      |
| `riplpy/densities/thfb_comb.py`                    | NEW       | Same API as `bsk14_comb`.                                                                                      |
| `riplpy/densities/__init__.py`                     | updated   | `_safe_load` wraps each section loader; `LEVEL_DENSITY_MODELS` extended with the four new keys.                |

## Config keys (riplpy/config.py)

- `densities_egsm` and `densities_egsm_norm` repointed under `densities/total/`.
- `densities_bsfg` / `densities_ct` retained (commented as legacy-only).
- Added directory keys:
  - `densities_bsk14_comb_dir`
  - `densities_bskg3_comb_dir`
  - `densities_qrpabe_dir`
  - `densities_thfb_comb_dir`

## Tests (`tests/test_densities.py`)

- BSFG, CT, HFB, MK shell-correction, and BSFG/CT-comparison tests are now
  guarded with `pytest.skip()` when their backing files are absent (which is the
  case for the github release).
- Added `TestCombinatorialLevelDensities`, exercising both `load_element` (per
  Z) and `load_all` for each of the four new readers.
- Added `TestDensitiesSectionLoad::test_section_load_does_not_raise` to verify
  that `riplpy.densities.load(directory)` succeeds even when legacy files are
  absent, populating the always-present databases (EGSM and the combinatorial
  models).
- Result: `16 passed, 11 skipped` against the github release.

## Known limitations / follow-ups

- The four new combinatorial readers reuse `hfb.read_ascii_file`. This parses
  the energy, temperature, cumulative count, observed/total densities, and
  the 50-spin grid that all four formats share. The associated `.ld`
  normalisation files (`bsk14-comb-ld.readme` format with `Z A Nlow Ntop alpha
  delta`) are **not yet parsed** — they are simple `(4i4,2f12.5)` records and
  could be added later as `<model>.load_norm()` if needed.
- `qrpabe` files use half-integer spin labels (`J=01/2 …`); the existing
  HFB parser stores them positionally in `rho_J`, so the indices represent
  J = 1/2, 3/2, … rather than 0, 1, 2 …. Document this if the data is plotted
  against spin.
- The legacy `level-densities-bfmeff.dat`, `level-densities-ctmeff.dat`,
  `shellmengoninakajima.dat`, and `level-densities-hfb/` directory are
  intentionally allowed to be missing. The section loader (`densities.load`)
  emits a `logger.warning` and continues.
