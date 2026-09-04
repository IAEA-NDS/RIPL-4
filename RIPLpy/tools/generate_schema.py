# -*- coding: utf-8 -*-
"""Generate the machine-readable RIPLpy data catalog (SCHEMA.json + CATALOG.md).

This walks every loaded RIPLpy database, introspects a representative entry,
and emits:

  * ``SCHEMA.json`` -- a machine-readable catalog: for each section/database,
    the key type, entry class, entry count, and every field with its Python
    type and ``[units]`` description (sourced from each entry's
    ``_field_info``). Intended for ML/AI ingestion and tooling.
  * ``CATALOG.md``  -- a human-readable rendering of the same.

Run from the repo root with the RIPL-4 path configured, e.g.::

    RIPL_LOCATION=/path/to/RIPL-4/github python tools/generate_schema.py

Heavy per-Z / per-nucleus array datasets are sampled with one element so the
schema is captured without loading the entire (multi-thousand-file) tree.

"""

import json
import os
import sys

import riplpy
from riplpy import Nuclide

OUT_JSON = os.path.join(os.path.dirname(__file__), os.pardir, 'SCHEMA.json')
OUT_MD = os.path.join(os.path.dirname(__file__), os.pardir, 'CATALOG.md')


def _key_type(db) -> str:
    if not getattr(db, 'data', None):
        return 'Nuclide'
    k = next(iter(db.data.keys()))
    return type(k).__name__


def _sample_entry(db):
    """Return one representative entry (unwrap PSF list-valued entries)."""
    if not getattr(db, 'data', None):
        return None
    v = next(iter(db.data.values()))
    if isinstance(v, list):
        return v[0] if v else None
    return v


def _entry_schema(entry) -> dict:
    """Extract {field: {type, description}} from an entry (dataclass/Packet)."""
    if entry is None:
        return {}
    if hasattr(entry, 'field_info') and callable(entry.field_info):
        info = entry.field_info()
        # Normalize the type string.
        for name, meta in info.items():
            t = str(meta.get('type', '')).replace("<class '", "").replace("'>", "")
            meta['type'] = t.replace('riplpy.collections.', '')
        return info
    return {}


def _heavy_samples():
    """Map section.attr -> a thunk loading a single element/nucleus.

    Heavy array datasets are skipped by ``riplpy.load()``; load a small slice
    of each so its schema can still be catalogued.
    """
    import riplpy.gamma as g
    import riplpy.fission as f
    import riplpy.densities as d
    fe = Nuclide(Z=26, A=56)
    return {
        ('gamma', 'gsf_d1m'): lambda: g.d1m.load_element(Z=26),
        ('gamma', 'smlo_e1'): lambda: g.smlo_e1.load_nucleus(Nuclide(Z=8, A=16)),
        ('gamma', 'smlo_m1'): lambda: g.smlo_m1.load_element(Z=26),
        ('gamma', 'tlo'): lambda: g.tlo.load_element(Z=34),
        ('gamma', 'gsf'): lambda: g.gsf.load_nucleus(fe),
        ('gamma', 'psf'): lambda: g.psf.load_category('oslo'),
        ('fission', 'hfbpath_d1m'): lambda: f.hfbpath.load_d1m(),
        ('fission', 'rmf_axial'): lambda: f.rmf.load_axial(),
        ('fission', 'nld_fis'): lambda: f.nld_fis.load('Max1'),
    }


def main() -> int:
    riplpy.load()

    catalog = {}
    heavy = _heavy_samples()

    for section in riplpy.list_sections():
        catalog[section] = {}
        module = getattr(riplpy, section)
        for db_name in riplpy.list_databases(section):
            db = getattr(getattr(module, 'db'), db_name)
            if db is None:
                # Heavy/legacy DB not populated by load(); try a small sample.
                thunk = heavy.get((section, db_name))
                if thunk is not None:
                    try:
                        db = thunk()
                    except Exception as exc:  # noqa: BLE001
                        catalog[section][db_name] = {
                            'status': f'unavailable ({exc})'
                        }
                        continue
                else:
                    catalog[section][db_name] = {
                        'status': 'legacy/absent in this release'
                    }
                    continue

            entry = _sample_entry(db)
            catalog[section][db_name] = {
                'status': 'available',
                'key_type': _key_type(db),
                'entry_class': type(entry).__name__ if entry is not None else None,
                'n_entries': len(db.data),
                'fields': _entry_schema(entry),
            }

    with open(OUT_JSON, 'w') as fp:
        json.dump(catalog, fp, indent=2, sort_keys=True)

    lines = [
        "# RIPLpy Data Catalog",
        "",
        "Machine-generated schema for every RIPLpy database (regenerate with",
        "`python tools/generate_schema.py`). The authoritative machine-readable",
        "form is `SCHEMA.json`; this file is the human-readable rendering.",
        "",
    ]
    for section in sorted(catalog):
        lines.append(f"## {section}")
        lines.append("")
        for db_name in sorted(catalog[section]):
            meta = catalog[section][db_name]
            if meta.get('status') != 'available':
                lines.append(f"### {db_name}")
                lines.append(f"_{meta.get('status')}_")
                lines.append("")
                continue
            lines.append(
                f"### {db_name}  ({meta['n_entries']} entries, "
                f"key={meta['key_type']}, entry={meta['entry_class']})"
            )
            lines.append("")
            fields = meta.get('fields') or {}
            if fields:
                lines.append("| Field | Type | Description |")
                lines.append("|-------|------|-------------|")
                for name, fmeta in fields.items():
                    desc = fmeta.get('description', '') or ''
                    lines.append(f"| `{name}` | {fmeta.get('type','')} | {desc} |")
            else:
                lines.append("_(no field metadata)_")
            lines.append("")

    with open(OUT_MD, 'w') as fp:
        fp.write("\n".join(lines))

    print(f"Wrote {OUT_JSON} and {OUT_MD}")
    print(f"Sections: {len(catalog)}; "
          f"databases: {sum(len(v) for v in catalog.values())}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
