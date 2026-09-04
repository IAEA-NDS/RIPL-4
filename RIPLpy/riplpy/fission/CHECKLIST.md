# Fission section — release checklist

Status of readers, writers, and tests against the RIPL-4 GitHub release
layout (``RIPL-4/github/fission/``).

## Data files in the RIPL-4 GitHub release

| File / directory | Reader module | Status |
|------------------|---------------|--------|
| ``empirical-barriers-ripl4.dat``           | ``riplpy.fission.empirical``     | Wired |
| ``empirical-barriers-ripl4.dat`` (alias)   | ``riplpy.fission.empirical_new`` | Wired (alias of ``empirical``) |
| ``empirical-barriers-new-EMPIRE.dat``      | ``riplpy.fission.empire``        | Wired |
| ``barriers-bskg3.dat``                     | ``riplpy.fission.bskg3``         | Wired (5-section RIPL-4 layout) |
| ``barriers-d1m_lep.dat``                   | ``riplpy.fission.d1m``           | Wired (new in RIPL-4) |
| ``hfbpath-bskg3/zXXX.dat``                 | ``riplpy.fission.hfbpath``       | Supported (opt-in via ``include_paths``) |
| ``hfbpath-d1m/zXXX.dat``                   | ``riplpy.fission.hfbpath``       | Supported (opt-in via ``include_paths``) |
| ``nld-fis-bskg3/{Max1,Max2,Max3,Min1,Min2}/zXXX.dat`` | ``riplpy.fission.nld_fis`` | Supported (opt-in via ``include_paths``) |
| ``RMF/Path_Axial/*.dat``                   | ``riplpy.fission.rmf``           | Supported (opt-in via ``include_paths``) |
| ``RMF/Path_Triaxial/*.dat``                | ``riplpy.fission.rmf``           | Supported (opt-in via ``include_paths``) |

## Legacy files removed from the RIPL-4 release

| File | Reader module | Behaviour |
|------|----------------|-----------|
| ``empirical-barriers.dat``        | ``riplpy.fission.empirical``       | Repointed to ``empirical-barriers-ripl4.dat`` |
| ``empirical-barriers-new.dat``    | ``riplpy.fission.empirical_new``   | Repointed to ``empirical-barriers-ripl4.dat`` (alias) |
| ``empirical-hfb-barriers.dat``    | ``riplpy.fission.hfb``             | Loader logs a warning and returns an empty database |
| ``HFB2007/``                      | -                                  | Deferred |
| ``leveldensities/``               | -                                  | Deferred |

## Config keys

In ``riplpy/config.py``:

| Key                                | Relative path                                |
|------------------------------------|----------------------------------------------|
| ``fission_barriers_empirical``     | ``fission/empirical-barriers-ripl4.dat``     |
| ``fission_barriers_empirical_new`` | ``fission/empirical-barriers-ripl4.dat``     |
| ``fission_barriers_empire``        | ``fission/empirical-barriers-new-EMPIRE.dat``|
| ``fission_barriers_hfb``           | ``fission/empirical-hfb-barriers.dat``       |
| ``fission_barriers_bskg3``         | ``fission/barriers-bskg3.dat``               |
| ``fission_barriers_d1m``           | ``fission/barriers-d1m_lep.dat``             |
| ``fission_path_bskg3_dir``         | ``fission/hfbpath-bskg3``                    |
| ``fission_path_d1m_dir``           | ``fission/hfbpath-d1m``                      |
| ``fission_nld_bskg3_dir``          | ``fission/nld-fis-bskg3``                    |
| ``fission_rmf_axial_dir``          | ``fission/RMF/Path_Axial``                   |
| ``fission_rmf_triaxial_dir``       | ``fission/RMF/Path_Triaxial``                |

## API

``riplpy.get_fission_barrier(Z, A, model=...)`` supports the model keys:
``empirical``, ``empirical_new``, ``empire``, ``hfb``, ``bskg3``, ``d1m``.
The ``BARRIER_MODELS`` dictionary in ``riplpy/fission/__init__.py``
maps these to the corresponding ``db.*`` attributes.

## Tests

Located in ``tests/test_fission.py``:

- ``TestEmpiricalBarriers``   — loads and validates the new RIPL-4 file.
- ``TestEmpiricalNew``        — exercises the backwards-compatibility alias.
- ``TestEmpire``              — implicit via API tests.
- ``TestHFBBarriers``         — skipped when the legacy file is absent.
- ``TestBSkG3Barriers``       — loads the 5-section RIPL-4 layout.
- ``TestD1MBarriers``         — loads the new D1M LEP file.
- ``TestFissionComparison``   — skipped if HFB data is unavailable.

## Fission-path / saddle level-density readers (opt-in)

These directories ship thousands of per-isotope tables. The readers are
implemented and wired, but **excluded from the eager ``fission.load()``
path by default** to keep ``riplpy.load()`` fast. They are opted in via
``fission.load(include_paths=True)`` (mirroring the ``include_heavy``
pattern in ``gamma`` / ``masses``). When ``include_paths`` is False,
``db.hfbpath_bskg3 = db.hfbpath_d1m = db.rmf_axial = db.rmf_triaxial =
db.nld_fis = None``. The per-element / per-model module loaders work
regardless of the flag.

### ``riplpy.fission.hfbpath`` — BSkG3 / D1M fission paths

Per-Z files ``zXXX.dat``, each concatenating several isotopes. Title line
parsed by regex (``Z=``, ``A=``, symbol, ``nbeta``, ``Egs``); ``#``/``*``
banner and column-header comment lines skipped. ``nbeta`` data lines per
Fortran ``(5f11.3,6x,i2)``; whitespace tokenisation with a **fixed-width
column fallback** for overflowed adjacent fields (e.g. d1m Np-236).

- Functions: ``load_bskg3()``, ``load_d1m()``, ``load_element(Z, model)``,
  ``load(model='bskg3')``.
- ``Entry`` schema (key = ``Nuclide``): ``n, Z, A, sym, nbeta, Egs`` plus
  parallel arrays ``beta20, beta22, beta30, dE`` (E-Egs [MeV]),
  ``MI`` (MI_ATD [hbar^2/MeV]), ``idx`` (-1 gs, 1 barrier, 2 well).
- Verified: ``load_bskg3()`` → **2491** nuclei; ``load_d1m()`` → **45**
  nuclei (7 Z files); every entry has ``len(beta20) == nbeta``.

### ``riplpy.fission.rmf`` — RMF fission paths

Per-nucleus files ``Sym<A>.dat`` under ``RMF/Path_Axial`` and
``RMF/Path_Triaxial``. ``#``-banner header parsed by regex (``Z``, ``A``,
symbol, ``nbeta``, ``Egs``, ``E(0+)``); Z trusted from the header,
filename is a cross-check. Data columns: ``beta20 beta22 beta30 E-EGS Mu``.

- Functions: ``load_axial()``, ``load_triaxial()``, ``load(layout='axial')``.
- ``Entry`` schema (key = ``Nuclide``): ``n, Z, A, sym, nbeta, Egs, E0``
  plus parallel arrays ``beta20, beta22, beta30, dE`` (E-EGS [MeV]), ``Mu``.
- Verified: ``load_axial()`` → **45** nuclei; ``load_triaxial()`` → **45**
  nuclei; every entry has ``len(beta20) == nbeta``.

### ``riplpy.fission.nld_fis`` — BSkG3 saddle/well level densities

Per-Z files ``zXXX.dat`` in subdirs ``Max1/Max2/Max3`` (inner / first-outer
/ second-outer saddles) and ``Min1/Min2`` (super-deformed wells). Two
blocks per isotope (positive then negative parity); ``*****`` banner middle
line parsed by regex (``Z``, ``A``, parity, ``b2``, ``b30``, ``b40``,
``Icr``), followed by a ``U[MeV] ... J=00..J=49`` column header and data
rows. Reuses the ``densities.hfb`` parsing approach.

- Functions: ``load(saddle)``, ``load_element(Z, saddle)`` where
  ``saddle in {Max1, Max2, Max3, Min1, Min2}``.
- ``Entry`` schema (key = ``Nuclide``): ``n, Z, A, label`` and
  ``positive_parity`` / ``negative_parity`` (``ParityData`` with
  ``b2, b30, b40, Icr, U, T, NCUMUL, RHOOBS, RHOTOT`` and ``rho_J``
  as a list-of-lists, one row per energy).
- Verified: ``load('Max1')`` → **2467** nuclei; 60-row U grid
  (0.25–200 MeV), 50 spin columns, both parities populated.

### ``fission.load()`` exposed attributes (opt-in)

``db.hfbpath_bskg3``, ``db.hfbpath_d1m``, ``db.rmf_axial``,
``db.rmf_triaxial``, ``db.nld_fis`` (Max1). Each is loaded through the
existing ``_safe_load`` wrapper so a missing directory yields an empty
database instead of crashing.

## Deferred work

``HFB2007/`` and ``leveldensities/`` legacy directories (RIPL-3) are
removed from the GitHub release; readers for these are not planned.
