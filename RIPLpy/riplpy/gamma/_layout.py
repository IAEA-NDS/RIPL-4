# -*- coding: utf-8 -*-
"""Locate RIPL-4 gamma data across the layouts it is distributed in.

The IAEA RIPL-4 distribution is treated as immutable, and its gamma section
does not match the layout earlier working copies used:

* Bulk per-nucleus datasets ship as zip archives -- ``gamma/d1m.zip``,
  ``gamma/smlo_E1.zip``, ``gamma/smlo_M1.zip`` -- rather than as unpacked
  directories.
* The GDR parameter tables sit flat at ``gamma/``; unpacked working copies
  filed them under ``gamma/gdr_parameters_exp_new/`` and
  ``gamma/gdr_parameters_exp&systematics/``.

Readers call the helpers here with every layout they accept, newest first,
instead of hardcoding one. Both a pristine RIPL-4 checkout and a previously
reorganised tree therefore work, and neither is privileged.

Archives are never unpacked in place: the RIPL-4 tree may be read-only, and
writing into it would dirty a user's checkout. They are extracted once into a
cache directory -- ``$RIPLPY_CACHE``, else ``~/.cache/riplpy`` -- and reused
while the archive's size and mtime are unchanged.
"""

# OS
import os as _os

# Logging
import logging as _logging

# Zip archives
import zipfile as _zipfile

# ========================

_logger = _logging.getLogger(__name__)

__all__ = ('resolve_data_file', 'resolve_data_dir', 'cache_root')


def cache_root() -> str:
    """Return the directory used to unpack RIPL-4 archives."""
    env = _os.environ.get('RIPLPY_CACHE')
    if env:
        return _os.path.expanduser(env)
    return _os.path.join(_os.path.expanduser('~'), '.cache', 'riplpy')


def resolve_data_file(root: str, *rel_candidates: str) -> str | None:
    """Return the first candidate that exists as a file beneath ``root``.

    Args:
        root: The RIPL-4 root directory.
        rel_candidates: Relative paths to try, in preference order.

    Returns:
        The absolute path of the first match, or None when none exist.
    """
    for rel in rel_candidates:
        if not rel:
            continue
        path = _os.path.join(root, rel)
        if _os.path.isfile(path):
            return path
    return None


def _is_safe_member(name: str) -> bool:
    """Reject archive members that would escape the extraction directory."""
    if _os.path.isabs(name) or name.startswith('/'):
        return False
    parts = _os.path.normpath(name).split(_os.sep)
    return '..' not in parts


def _extract_archive(archive: str, dest_parent: str, expect: str) -> str | None:
    """Unpack ``archive`` into ``dest_parent`` and return the payload directory.

    The extraction is cached: a stamp file records the archive's size and mtime,
    and the archive is only unpacked again when those change.
    """
    try:
        stat = _os.stat(archive)
    except OSError:
        return None

    stamp_id = f"{stat.st_size}-{int(stat.st_mtime)}"
    dest = _os.path.join(dest_parent, expect)
    stamp = _os.path.join(dest_parent, f".{expect}.stamp")

    if _os.path.isdir(dest) and _os.path.isfile(stamp):
        try:
            with open(stamp, 'r', encoding='utf-8') as fh:
                if fh.read().strip() == stamp_id:
                    return dest
        except OSError:
            pass

    _logger.info(f"Unpacking {archive} into {dest_parent} (first use; cached thereafter)")
    try:
        _os.makedirs(dest_parent, exist_ok=True)
        with _zipfile.ZipFile(archive) as zf:
            members = [n for n in zf.namelist() if _is_safe_member(n)]
            skipped = len(zf.namelist()) - len(members)
            if skipped:
                _logger.warning(f"{archive}: skipped {skipped} unsafe archive member(s)")
            zf.extractall(dest_parent, members=members)
    except (OSError, _zipfile.BadZipFile) as exc:
        _logger.warning(f"Could not unpack {archive}: {exc}")
        return None

    if not _os.path.isdir(dest):
        # Archive did not contain the expected top-level folder; fall back to
        # the extraction root when it holds the payload directly.
        _logger.debug(f"{archive}: no top-level '{expect}/' directory after extraction")
        return dest_parent if _os.listdir(dest_parent) else None

    try:
        with open(stamp, 'w', encoding='utf-8') as fh:
            fh.write(stamp_id)
    except OSError:
        pass

    return dest


def resolve_data_dir(root: str, *rel_candidates: str) -> str | None:
    """Return a readable directory for the first candidate that resolves.

    A candidate resolves either as a real directory beneath ``root`` or as a
    sibling ``<candidate>.zip`` archive, which is unpacked into the cache and
    the unpacked directory returned.

    Args:
        root: The RIPL-4 root directory.
        rel_candidates: Relative directory paths to try, in preference order.

    Returns:
        A directory path whose contents the caller can list, or None.
    """
    # A real directory always wins, at any candidate position.
    for rel in rel_candidates:
        if not rel:
            continue
        path = _os.path.join(root, rel)
        if _os.path.isdir(path):
            return path

    for rel in rel_candidates:
        if not rel:
            continue
        archive = _os.path.join(root, rel + '.zip')
        if _os.path.isfile(archive):
            dest_parent = _os.path.join(cache_root(), _os.path.dirname(rel) or '.')
            resolved = _extract_archive(archive, dest_parent,
                                        _os.path.basename(rel))
            if resolved is not None:
                return resolved
    return None
