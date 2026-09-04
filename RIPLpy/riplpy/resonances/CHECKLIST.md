# Resonances Section Checklist

Status: COMPLETE
Last updated: 2026-05-14

## Source layout

The RIPL-4 github release ships the two resonance files under
`resonances/` in the data tree:

| File                     | Old (RIPL-3) name   | Description                                    |
|--------------------------|---------------------|------------------------------------------------|
| `resonances_L0.dat`      | `resonances0.dat`   | s-wave (L=0) average neutron-resonance params  |
| `resonances_L1.dat`      | `resonances1.dat`   | p-wave (L=1) average neutron-resonance params  |
| `resonances.readme`      | unchanged           | Format and provenance notes                    |

Coverage:
- s-wave: 324 nuclides
- p-wave: 248 nuclides

## Readers / Writers

| Module                          | Path key            | Reader / Writer            |
|---------------------------------|---------------------|----------------------------|
| `riplpy.resonances.swave`       | `resonances_swave`  | shared with `pwave` module |
| `riplpy.resonances.pwave`       | `resonances_pwave`  | `read_ascii_file` / `write_ascii_file` |

Both files share an identical FORTRAN format, so `swave.py` re-exports the
reader, writer and `Entry` dataclass from `pwave.py`.

`riplpy/config.py` `DATA_FILES` now points at:
- `'resonances_swave': resonances/resonances_L0.dat`
- `'resonances_pwave': resonances/resonances_L1.dat`

## Format

```
(3i4, 2x, a2, 2x, i2, f5.1, 1x, a1, f10.3, 1p, 12e12.3)
```

Each record carries: `Z, N, A, Sym, L, J, P, Sn` followed by 12 `e12.3`
fields giving `D, dD, Gg, dGg, S, dS` for each of the two evaluations
(RIPL-3 then BNL/Mughabghab 2018). Numeric fields may be blank to indicate
"no value".

## Field changes vs. previous RIPLpy

The on-disk format changed substantively for RIPL-4; the `Entry` dataclass
was extended accordingly. Backward-compatible field names that downstream
code (and tests) depend on are preserved:

| Field            | Status     | Notes                                              |
|------------------|------------|----------------------------------------------------|
| `n`              | kept       | Target `Nuclide`                                   |
| `Bn`             | kept       | Neutron binding energy [MeV]                       |
| `D`, `Derr`      | kept       | Resonance spacing (RIPL-3). UNITS: now eV (was keV)|
| `Gam`, `Gerr`    | kept       | Radiative width  (RIPL-3). UNITS: now eV (was meV) |
| `Str`, `Serr`    | kept       | Strength function (RIPL-3) [10^-4]                 |
| `Io`             | kept       | Target ground-state spin                           |
| `ref`            | removed    | No longer present in the data file                 |
| `sym`            | new        | Element symbol from the data file                  |
| `L`              | new        | Incident-neutron angular momentum (0 or 1)         |
| `parity`         | new        | Target ground-state parity ('+' / '-')             |
| `D_BNL`,    `D_BNL_err`    | new | Resonance spacing (BNL) [eV]                |
| `Gam_BNL`,  `Gam_BNL_err`  | new | Radiative width  (BNL) [eV]                 |
| `Str_BNL`,  `Str_BNL_err`  | new | Strength function (BNL) [10^-4]             |

Unit changes (eV vs keV / meV) reflect the RIPL-4 release; all `Entry`
field docstrings call out the new units explicitly.

## Test verification

`tests/test_resonances.py`: **7 / 7 pass** (one previously-skipped
roundtrip test re-enabled).

```
TestSwaveResonances::test_load_from_directory             PASSED
TestSwaveResonances::test_get_resonance_spacing           PASSED
TestSwaveResonances::test_resonance_values_reasonable     PASSED
TestSwaveResonances::test_write_and_read_roundtrip        PASSED
TestPwaveResonances::test_load_from_directory             PASSED
TestPwaveResonances::test_pwave_has_fewer_entries         PASSED
TestResonanceComparison::test_swave_and_pwave_*           PASSED
```

`tests/test_api.py::TestConvenienceFunctionsDirectLoad::test_get_resonance_direct`
also passes.
