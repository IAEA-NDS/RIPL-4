# Levels Section Checklist

## Data Sources (GitHub Layout)

Expected files under `<RIPL_PATH>/levels/`:

| File | Purpose |
|------|---------|
| `levels-param.data` | Constant Temperature (CT) fit parameters |
| `levels-readme.html` | Documentation (not parsed) |
| `z001.dat` ... `z118.dat` | Discrete level schemes (one file per element Z) |

All files live directly in `levels/` (flat layout, no `data/discrete/` subdir).

## File-to-Reader Mapping

| File | Module | Reader | Writer | Database |
|------|--------|--------|--------|----------|
| `levels-param.data` | `riplpy/levels/ct.py` | `read_ascii_file` | `write_ascii_file` | `ct.Database` |
| `z{Z:03d}.dat` | `riplpy/levels/discrete.py` | `read_ascii_file` | `write_ascii_file` | `discrete.Database` |

## Path Configuration (`riplpy/config.py`)

| Key | Relative Path | Status |
|-----|---------------|--------|
| `levels_param` | `levels/levels-param.data` | OK |
| `discrete_levels` | `levels/z{Z:03d}.dat` | Fixed (was `levels/data/discrete/z{Z:03d}.dat`) |

`discrete.py` uses its own `local_data_dir = "levels"` constant and joins `<directory>/levels/z{Z:03d}.dat`, which is correct for the github layout.

`ct.py` uses `FILE_PATH_KEY = 'levels_param'` resolved via `config.get_data_file_path()`.

## Test Verification

Command: `python -m pytest tests/test_levels.py -xvs`

Result: **11/11 tests passing**

- `TestLevelsCT` (5 tests): load, get parameters, attribute presence, temperature sanity, write/read roundtrip
- `TestDiscreteLevels` (5 tests): load element, Fe-56, spin/parity, gammas, energy ordering
- `TestLevelsCoverage` (1 test): covers >50 elements

## Notes / Deferred

- `read_ascii_file` in `discrete.py` carries some pre-existing comments referencing legacy parsing quirks; no behavioral changes were made.
- `write_ascii_file` in `discrete.py` still has `@TODO` markers for placeholder integer/float fields in the header record; roundtrip-comparable but not field-perfect.
- HTML readme is not parsed (intentional).
