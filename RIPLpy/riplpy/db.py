# -*- coding: utf-8 -*-
"""The base objects which define access to RIPL databases.

Python 3.10+ is expected to run this code properly.

"""

# OS
import os as _os

# CSV
import csv as _csv

# JSON
import json as _json

# Dataclasses
from dataclasses import dataclass as _dataclass, fields as _fields, is_dataclass as _is_dataclass

# Typing
from typing import ClassVar as _ClassVar

# RIPLpy
from riplpy.exceptions import NucleusNotFoundError, ElementNotFoundError, FileFormatError, RiplFileNotFoundError

__all__ = ('DatabaseEntry', 'NuclideDatabaseEntry', 'PacketEntry', 'Database', 'NuclideDatabase', 'ElementDatabase', '_DbAccessor', 'load')

@_dataclass
class DatabaseEntry:
    """A generic database entry. This class is expected to be extended. """

    # Class-level metadata for field descriptions (override in subclasses)
    # Format: {'field_name': 'description with [units]'}
    _field_info: _ClassVar[dict] = {}

    @property
    def as_list(self, ) -> list:
        """Return a list representation of the object. """
        return [getattr(self, f.name) for f in _fields(self.__class__)]

    @property
    def as_tuple(self, ) -> tuple:
        """Return a tuple representation of the object. Intended for writing to file. """
        return tuple(self.as_list)

    @property
    def as_dict(self, ) -> dict:
        """Return a dictionary representation of the object. Useful packaging for other objects. """
        d = {}
        for f in _fields(self.__class__):
            d[f.name] = getattr(self, f.name)
        return d

    def __repr__(self) -> str:
        """Return a pretty string representation of the entry."""
        class_name = type(self).__name__
        field_strs = []
        for f in _fields(self.__class__):
            value = getattr(self, f.name)
            # Format based on type
            if isinstance(value, float):
                field_strs.append(f"{f.name}={value:.4g}")
            elif isinstance(value, str) and len(value) > 20:
                field_strs.append(f"{f.name}='{value[:17]}...'")
            elif hasattr(value, 'Z') and hasattr(value, 'A'):
                # Nucleus-like object
                field_strs.append(f"{f.name}={value.element_symbol}-{value.A}")
            else:
                field_strs.append(f"{f.name}={value!r}")
        return f"{class_name}({', '.join(field_strs)})"

    def summary(self) -> str:
        """Return a formatted summary of the entry with field descriptions."""
        lines = [f"Entry: {type(self).__name__}"]
        lines.append("-" * 40)

        # Get field info if available
        field_info = getattr(type(self), '_field_info', None) or {}

        for f in _fields(self.__class__):
            value = getattr(self, f.name)
            desc = field_info.get(f.name, '')

            # Format value
            if isinstance(value, float):
                val_str = f"{value:.6g}"
            elif hasattr(value, 'Z') and hasattr(value, 'A'):
                val_str = f"{value.element_symbol}-{value.A} (Z={value.Z}, N={value.N})"
            else:
                val_str = repr(value)

            # Build line
            if desc:
                lines.append(f"  {f.name}: {val_str}")
                lines.append(f"       {desc}")
            else:
                lines.append(f"  {f.name}: {val_str}")

        return "\n".join(lines)

    @classmethod
    def fields(cls) -> list:
        """Return list of field names for this entry type."""
        return [f.name for f in _fields(cls)]

    @classmethod
    def field_info(cls) -> dict:
        """Return field information including types and descriptions."""
        info = {}
        field_descriptions = getattr(cls, '_field_info', None) or {}

        for f in _fields(cls):
            info[f.name] = {
                'type': str(f.type) if f.type else 'Any',
                'description': field_descriptions.get(f.name, ''),
            }
        return info


@_dataclass
class NuclideDatabaseEntry(DatabaseEntry):

    @property
    def as_tuple(self, ) -> tuple:
        """Return a tuple representation of the object. Intended for writing to file. """
        lst = [self.n.Z, self.n.A, self.n.element_symbol]
        for f in _fields(self.__class__):
            if f.name != 'n':
                lst.append(getattr(self, f.name))
        return tuple(lst)


class PacketEntry:
    """A uniform, dict-backed database entry for array/spectral datasets.

    The scalar RIPL sections use frozen ``@dataclass`` entries. The large
    RIPL-4 array/spectral products (gamma D1M/SMLO/TLO/GSF, the combinatorial
    level densities, fission paths, …) are naturally variable-width tables and
    were historically stored as bare ``dict`` packets. ``PacketEntry`` wraps
    such a packet so that *every* database — scalar or spectral — yields an
    entry exposing the same interface (``.as_dict``/``.as_list``/``.fields()``/
    ``.field_info()``/``.summary()`` and ``__repr__``) used by the export
    machinery (``to_dataframe``/``to_json``/``to_csv``).

    Both attribute access (``entry.U``) and item access (``entry['U']``) work,
    so existing code and documentation that treated packets as plain dicts
    keeps working unchanged.

    Subclasses may set ``_field_info`` (``{name: 'description [units]'}``) for
    schema/units documentation, exactly like the dataclass entries.
    """

    _field_info: _ClassVar[dict] = {}

    def __init__(self, packet: dict | None = None, **kwargs) -> None:
        data = dict(packet) if packet else {}
        data.update(kwargs)
        # Store under a private name so attribute access cannot shadow it.
        object.__setattr__(self, '_data', data)

    # -- mapping-style access (backwards compatible with the old dict packets)
    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value

    def __contains__(self, key) -> bool:
        return key in self._data

    def get(self, key, default=None):
        return self._data.get(key, default)

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    def __iter__(self):
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    # -- attribute-style access
    def __getattr__(self, name):
        # __getattr__ is only called when normal lookup fails, so _data is safe.
        try:
            return self._data[name]
        except KeyError as exc:
            raise AttributeError(name) from exc

    # -- uniform Entry interface (parallels DatabaseEntry)
    @property
    def as_dict(self) -> dict:
        """Return the underlying packet as a plain dictionary."""
        return dict(self._data)

    @property
    def as_list(self) -> list:
        """Return the packet values in key order."""
        return list(self._data.values())

    @property
    def as_tuple(self) -> tuple:
        """Return the packet values in key order as a tuple."""
        return tuple(self._data.values())

    def fields(self) -> list:
        """Return the packet's field names."""
        return list(self._data.keys())

    def field_info(self) -> dict:
        """Return field information (type + description) for this entry."""
        info = {}
        descriptions = getattr(type(self), '_field_info', None) or {}
        for name, value in self._data.items():
            info[name] = {
                'type': type(value).__name__,
                'description': descriptions.get(name, ''),
            }
        return info

    def summary(self) -> str:
        """Return a formatted summary of the entry with field descriptions."""
        lines = [f"Entry: {type(self).__name__}", "-" * 40]
        descriptions = getattr(type(self), '_field_info', None) or {}
        for name, value in self._data.items():
            if isinstance(value, (list, tuple)):
                val_str = f"<{type(value).__name__} of {len(value)}>"
            elif hasattr(value, 'Z') and hasattr(value, 'A'):
                val_str = f"{value.element_symbol}-{value.A} (Z={value.Z}, N={value.N})"
            elif isinstance(value, float):
                val_str = f"{value:.6g}"
            else:
                val_str = repr(value)
            desc = descriptions.get(name, '')
            lines.append(f"  {name}: {val_str}")
            if desc:
                lines.append(f"       {desc}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        parts = []
        for name, value in self._data.items():
            if isinstance(value, (list, tuple)):
                parts.append(f"{name}=<{type(value).__name__} of {len(value)}>")
            elif hasattr(value, 'Z') and hasattr(value, 'A'):
                parts.append(f"{name}={value.element_symbol}-{value.A}")
            elif isinstance(value, float):
                parts.append(f"{name}={value:.4g}")
            else:
                parts.append(f"{name}={value!r}")
        return f"{type(self).__name__}({', '.join(parts)})"


class Database:
    """A generic representation of a RIPL database that accesses information from a single ASCII file. 
    
       More complicated databases may have multiple reader and writer methods.
    """

    reader: object = None # The method to read data from a specified ASCII file
    entry : object = None # The object which represents an entry in the database
    writer: object = None # The method to write data to a specified ASCII file

    def __init__(self, data: dict | None = None) -> None:
        """Initialize the object. """
        self.data = data if data is not None else {}

    def __repr__(self) -> str:
        """Return a string representation of the database."""
        class_name = type(self).__name__
        module_name = type(self).__module__

        # Try to get a short module path (e.g., 'masses.ame20' from 'riplpy.masses.ame20')
        if module_name.startswith('riplpy.'):
            module_name = module_name[7:]  # Remove 'riplpy.' prefix

        count = len(self.data)

        if count == 0:
            return f"<{class_name}: {module_name} (empty)>"

        # Try to get Z range for nuclide databases
        z_info = ""
        if self.data:
            keys = list(self.data.keys())
            if hasattr(keys[0], 'Z'):
                z_values = [k.Z for k in keys if hasattr(k, 'Z')]
                if z_values:
                    z_min, z_max = min(z_values), max(z_values)
                    z_info = f", Z={z_min}-{z_max}"

        return f"<{class_name}: {module_name} ({count} entries{z_info})>"

    def info(self) -> str:
        """Return detailed information about this database.

        Returns:
            A formatted string with database statistics and metadata.

        Examples:
            >>> print(db.info())
            Database: masses.ame20
            ----------------------------------------
            Type: Database
            Entries: 3558
            Z range: 0 - 118
            Entry type: Entry
            Fields: n, flag, Mexp, Err
        """
        lines = []
        class_name = type(self).__name__
        module_name = type(self).__module__
        if module_name.startswith('riplpy.'):
            module_name = module_name[7:]

        lines.append(f"Database: {module_name}")
        lines.append("-" * 40)
        lines.append(f"Type: {class_name}")
        lines.append(f"Entries: {len(self.data)}")

        # Z/A range for nuclide databases
        if self.data:
            keys = list(self.data.keys())
            if hasattr(keys[0], 'Z'):
                z_values = [k.Z for k in keys if hasattr(k, 'Z')]
                a_values = [k.A for k in keys if hasattr(k, 'A')]
                if z_values:
                    lines.append(f"Z range: {min(z_values)} - {max(z_values)}")
                if a_values:
                    lines.append(f"A range: {min(a_values)} - {max(a_values)}")

        # Entry type info
        if self.entry:
            entry_name = self.entry.__name__ if hasattr(self.entry, '__name__') else str(self.entry)
            lines.append(f"Entry type: {entry_name}")

            # Field names
            if hasattr(self.entry, 'fields') and callable(self.entry.fields):
                fields = self.entry.fields()
                lines.append(f"Fields: {', '.join(fields)}")

        return "\n".join(lines)

    def describe(self) -> str:
        """Return detailed field descriptions for entries in this database.

        Returns:
            A formatted string describing each field in the entry type.

        Examples:
            >>> print(db.describe())
            Entry Fields for: Entry
            ----------------------------------------
            n          : Nuclide - Target nucleus
            Mexp       : float   - Experimental mass excess [MeV]
            Err        : float   - Mass excess uncertainty [MeV]
        """
        lines = []

        if not self.entry:
            return "No entry type defined for this database."

        entry_name = self.entry.__name__ if hasattr(self.entry, '__name__') else str(self.entry)
        lines.append(f"Entry Fields for: {entry_name}")
        lines.append("-" * 40)

        # Get field info if available
        if hasattr(self.entry, 'field_info') and callable(self.entry.field_info):
            field_info = self.entry.field_info()
            for name, info in field_info.items():
                type_str = info.get('type', 'Any')
                # Simplify type strings
                type_str = type_str.replace("<class '", "").replace("'>", "")
                type_str = type_str.replace("riplpy.collections.", "")
                desc = info.get('description', '')

                if desc:
                    lines.append(f"  {name:12}: {type_str:10} - {desc}")
                else:
                    lines.append(f"  {name:12}: {type_str}")
        elif hasattr(self.entry, 'fields') and callable(self.entry.fields):
            # Fallback: just list field names
            for name in self.entry.fields():
                lines.append(f"  {name}")

        return "\n".join(lines)

    def sample(self, n: int = 1) -> object | list:
        """Return sample entry/entries from the database.

        Args:
            n: Number of samples to return (default: 1)

        Returns:
            Single entry if n=1, otherwise list of entries

        Examples:
            >>> entry = db.sample()
            >>> entries = db.sample(5)
        """
        if not self.data:
            return None if n == 1 else []

        import itertools
        keys = list(itertools.islice(self.data.keys(), n))

        if n == 1:
            return self.data[keys[0]]
        return [self.data[k] for k in keys]

    def head(self, n: int = 5) -> None:
        """Print the first n entries in a formatted table.

        Args:
            n: Number of entries to display (default: 5)

        Examples:
            >>> db.head()
            >>> db.head(10)
        """
        if not self.data:
            print("(empty database)")
            return

        keys = list(self.data.keys())[:n]
        for key in keys:
            entry = self.data[key]
            print(entry)
        if len(self.data) > n:
            print(f"... ({len(self.data) - n} more entries)")

    def load(self, fpath: str) -> None:
        """A generic function for Loading the database from an ASCII file into the data dictionary. 
           This method may need to be extended; it may not work in every case. 
           The success of loading depends sensitively on the structure of the input file.
        """
        self.data = {}
        data = type(self).reader(fpath) # This call returns the data dictionary, but not yet formatted for the db
        # Loop over the data entries
        for key in data.keys():
            # Properly format each entry for the specified db
            self.data[key] = self.entry(**data[key])

    def save(self, fpath: str) -> None:
        """Save the database to an ASCII file. """
        type(self).writer(fpath, self.data)

    def get(self, key: object) -> object:
        """Return the entry from the database given the key.

        Args:
            key (object): The key to retrieve the entry.

        Returns:
            object: The entry from the database.

        Raises:
            NucleusNotFoundError: If the key is not found in a NuclideDatabase.
            ElementNotFoundError: If the key is not found in an ElementDatabase.
            KeyError: If the key is not found in a generic Database.
        """
        try:
            return self.data[key]
        except KeyError:
            if isinstance(self, NuclideDatabase):
                raise NucleusNotFoundError(f"Nuclide with key {key} not found in database.")
            elif isinstance(self, ElementDatabase):
                raise ElementNotFoundError(f"Element with key {key} not found in database.")
            else:
                raise KeyError(f"Key {key} not found in database.")

    def set(self, key: object, value: object) -> None:
        """Set the entry in the database given the key and value. """
        self.data[key] = value

    def insert(self, key: object, value: object) -> None:
        """Insert the entry in the database given the key and value. """
        self.set(key, value)

    def remove(self, key: object) -> None:
        """Remove an existing entry in the database by providing a key (typically a nucleus-like) object. """
        self.data.pop(key)

    def clear(self, ) -> None:
        """Clear the database. """
        del self.data
        self.data = {}

    def _serialize_value(self, value: object) -> object:
        """Serialize a value for export to JSON/CSV.

        Handles Nuclide/Nucleus objects, nested dicts, lists, and dataclasses.
        """
        if value is None:
            return None
        elif isinstance(value, bool):
            return value
        elif isinstance(value, (int, float, str)):
            return value
        elif hasattr(value, 'Z') and hasattr(value, 'A') and hasattr(value, 'element_symbol'):
            # Nuclide / Nucleus / IsomericNucleus -- has the full nuclear-id
            # surface. Other Z/A-bearing objects (e.g. coupled-channel
            # ``RotationalIsotope``) fall through to the dataclass / vars()
            # path below so all their fields are serialised.
            return {'Z': value.Z, 'A': value.A, 'symbol': value.element_symbol}
        elif hasattr(value, 'item') and hasattr(value, 'dtype') and getattr(value, 'shape', None) == ():
            # 0-d numpy scalar -> native Python scalar
            return value.item()
        elif hasattr(value, 'tolist') and hasattr(value, 'dtype'):
            # numpy ndarray -> nested Python list (round-trippable)
            return self._serialize_value(value.tolist())
        elif hasattr(value, 'as_dict'):
            # DatabaseEntry / PacketEntry instance
            return {k: self._serialize_value(v) for k, v in value.as_dict.items()}
        elif _is_dataclass(value) and not isinstance(value, type):
            # Dataclass instance
            return {f.name: self._serialize_value(getattr(value, f.name)) for f in _fields(value)}
        elif isinstance(value, dict):
            return {k: self._serialize_value(v) for k, v in value.items()}
        elif isinstance(value, (list, tuple)):
            return [self._serialize_value(v) for v in value]
        else:
            return str(value)

    @staticmethod
    def _entry_as_mapping(entry: object) -> dict:
        """Return an entry's fields as a plain dict, whatever its container.

        Handles dataclass entries, ``PacketEntry``/objects with ``as_dict``,
        bare ``dict`` packets, and falls back to ``vars()`` for plain objects.
        """
        if hasattr(entry, 'as_dict'):
            return entry.as_dict
        if isinstance(entry, dict):
            return dict(entry)
        try:
            return vars(entry)
        except TypeError:
            return {'value': entry}

    def _flatten_entry(self, entry: object, prefix: str = '',
                       stringify_sequences: bool = True) -> dict:
        """Flatten an entry into a single-level dictionary for tabular export.

        Nested structures are flattened with dot notation (e.g., 'inner.E[MeV]').

        Args:
            entry: A dataclass entry, ``PacketEntry``, or dict packet.
            prefix: Key prefix used during recursion.
            stringify_sequences: If True (CSV), list/array values are stored as
                JSON strings so each cell is scalar. If False (DataFrame), they
                are kept as native Python lists so downstream code can consume
                the arrays directly.
        """
        flat = {}
        entry_dict = self._entry_as_mapping(entry)

        for key, value in entry_dict.items():
            full_key = f"{prefix}{key}" if prefix else key

            if value is None:
                flat[full_key] = None
            elif isinstance(value, bool):
                flat[full_key] = value
            elif isinstance(value, (int, float, str)):
                flat[full_key] = value
            elif hasattr(value, 'Z') and hasattr(value, 'A'):
                # Nucleus-like object - expand to Z, A, symbol
                flat[f"{full_key}.Z"] = value.Z
                flat[f"{full_key}.A"] = value.A
                flat[f"{full_key}.symbol"] = value.element_symbol
            elif hasattr(value, 'item') and hasattr(value, 'dtype') and getattr(value, 'shape', None) == ():
                # 0-d numpy scalar
                flat[full_key] = value.item()
            elif hasattr(value, 'tolist') and hasattr(value, 'dtype'):
                # numpy ndarray
                lst = value.tolist()
                flat[full_key] = _json.dumps(lst) if stringify_sequences else lst
            elif _is_dataclass(value) and not isinstance(value, type):
                # Recursively flatten nested dataclass
                nested = self._flatten_entry(value, prefix=f"{full_key}.",
                                             stringify_sequences=stringify_sequences)
                flat.update(nested)
            elif hasattr(value, 'as_dict'):
                # Nested PacketEntry / Entry
                nested = self._flatten_entry(value, prefix=f"{full_key}.",
                                             stringify_sequences=stringify_sequences)
                flat.update(nested)
            elif isinstance(value, dict):
                for k, v in value.items():
                    if isinstance(v, (list, tuple)):
                        flat[f"{full_key}.{k}"] = (
                            _json.dumps(self._serialize_value(v))
                            if stringify_sequences else self._serialize_value(v)
                        )
                    else:
                        flat[f"{full_key}.{k}"] = self._serialize_value(v)
            elif isinstance(value, (list, tuple)):
                serial = self._serialize_value(value)
                flat[full_key] = _json.dumps(serial) if stringify_sequences else serial
            else:
                flat[full_key] = str(value)

        return flat

    def _iter_export_entries(self):
        """Yield ``(key, entry)`` pairs for export.

        Databases that store *multiple* datasets per key as a ``list`` (e.g.
        the PSF experimental database, which holds several measurements per
        nucleus) are exploded into one record per dataset so every export
        format yields a uniform, flat row/record stream.
        """
        for key, value in self.data.items():
            if isinstance(value, list):
                for item in value:
                    yield key, item
            else:
                yield key, value

    def to_csv(self, fpath: str) -> None:
        """Export the database to a CSV file.

        Args:
            fpath: Path to the output CSV file.
        """
        if not self.data:
            return

        # Collect all rows and determine all column names
        rows = []
        all_columns = set()

        for key, entry in self._iter_export_entries():
            row = self._flatten_entry(entry, stringify_sequences=True)
            rows.append(row)
            all_columns.update(row.keys())

        # Sort columns for consistent output
        # Put Z, A, symbol first if present, then sort alphabetically
        priority_cols = ['n.Z', 'n.A', 'n.symbol', 'Z', 'A', 'symbol']
        sorted_cols = []
        for col in priority_cols:
            if col in all_columns:
                sorted_cols.append(col)
                all_columns.discard(col)
        sorted_cols.extend(sorted(all_columns))

        # Write CSV
        with open(fpath, 'w', newline='', encoding='utf-8') as fp:
            writer = _csv.DictWriter(fp, fieldnames=sorted_cols, extrasaction='ignore')
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    def to_json(self, fpath: str, indent: int = 2) -> None:
        """Export the database to a JSON file.

        Args:
            fpath: Path to the output JSON file.
            indent: Number of spaces for indentation (default: 2).
        """
        if not self.data:
            with open(fpath, 'w', encoding='utf-8') as fp:
                _json.dump([], fp)
            return

        # Serialize all entries
        entries = []
        for key, entry in self._iter_export_entries():
            serialized = self._serialize_value(self._entry_as_mapping(entry))
            entries.append(serialized)

        with open(fpath, 'w', encoding='utf-8') as fp:
            _json.dump(entries, fp, indent=indent)

    def to_dataframe(self) -> "pandas.DataFrame":
        """Export the database to a pandas DataFrame.

        Returns:
            pandas.DataFrame: A DataFrame with flattened entry data.

        Raises:
            ImportError: If pandas is not installed.
        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError(
                "pandas is required for to_dataframe(). Install the optional "
                "DataFrame support with: pip install 'riplpy[dataframe]' "
                "(or directly: pip install pandas)"
            )

        if not self.data:
            return pd.DataFrame()

        # Collect all rows (keep list/array fields as native Python lists so
        # downstream ML code can consume the arrays directly).
        rows = []
        for key, entry in self._iter_export_entries():
            row = self._flatten_entry(entry, stringify_sequences=False)
            rows.append(row)

        # Create DataFrame
        df = pd.DataFrame(rows)

        # Reorder columns: Z, A, symbol first
        priority_cols = ['n.Z', 'n.A', 'n.symbol', 'Z', 'A', 'symbol']
        existing_priority = [col for col in priority_cols if col in df.columns]
        other_cols = sorted([col for col in df.columns if col not in priority_cols])
        df = df[existing_priority + other_cols]

        return df

    def to_list(self) -> list:
        """Export the database as a list of serialized entry dictionaries.

        Returns:
            list: A list of dictionaries, one per entry.
        """
        entries = []
        for key, entry in self._iter_export_entries():
            serialized = self._serialize_value(self._entry_as_mapping(entry))
            entries.append(serialized)
        return entries

    def to_flat_list(self) -> list:
        """Export the database as a list of flattened entry dictionaries.

        Nested structures are flattened with dot notation (e.g., 'inner.E[MeV]').

        Returns:
            list: A list of flattened dictionaries, one per entry.
        """
        rows = []
        for key, entry in self._iter_export_entries():
            row = self._flatten_entry(entry, stringify_sequences=False)
            rows.append(row)
        return rows

    def to_numpy(self, structured: bool = True):
        """Export the database as a numpy array for ML/AI pipelines.

        Args:
            structured: If True (default), return a numpy *structured* (record)
                array mirroring :meth:`to_dataframe`: every column is preserved
                with the nuclide Z/A/symbol columns first, named fields, and
                per-column dtypes inferred (all-int -> int64, numeric with gaps
                -> float64 with NaN, list/spectral fields -> object). Access a
                column by name, e.g. ``arr['Mexp']``.

                If False, return ``(X, columns)`` where ``X`` is a dense
                float64 2-D feature matrix of the *numeric* columns only
                (missing values become NaN) and ``columns`` is the list of
                column names. Non-numeric (string) and list-valued columns are
                dropped, so ``X`` is model-ready.

        Returns:
            numpy.ndarray when ``structured`` is True; otherwise a
            ``(numpy.ndarray, list[str])`` tuple.

        Raises:
            ImportError: If numpy is not installed.

        Examples:
            >>> arr = db.to_numpy()
            >>> arr['Mexp']                     # float64 column
            >>> X, cols = db.to_numpy(structured=False)
            >>> X.shape                         # (n_rows, n_numeric_cols)
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError(
                "numpy is required for to_numpy(). Install the optional array "
                "support with: pip install 'riplpy[numpy]' "
                "(or directly: pip install numpy)"
            )

        rows = self.to_flat_list()

        # Column ordering mirrors to_dataframe: Z/A/symbol first, then sorted.
        all_cols = set()
        for row in rows:
            all_cols.update(row.keys())
        priority_cols = ['n.Z', 'n.A', 'n.symbol', 'Z', 'A', 'symbol']
        ordered = [c for c in priority_cols if c in all_cols]
        ordered += sorted(c for c in all_cols if c not in priority_cols)

        # Per-column value vectors (a missing key in a row becomes None).
        col_values = {c: [row.get(c) for row in rows] for c in ordered}

        def _classify(values):
            """Return (numpy dtype, normalized values) for one column."""
            non_none = [v for v in values if v is not None]
            is_numeric = bool(non_none) and all(
                isinstance(v, (int, float)) for v in non_none
            )
            if is_numeric:
                all_int = all(isinstance(v, int) for v in non_none)
                if all_int and None not in values:
                    return np.int64, [int(v) for v in values]
                return np.float64, [
                    float(v) if v is not None else np.nan for v in values
                ]
            return object, list(values)

        if not structured:
            numeric_cols = [
                c for c in ordered if _classify(col_values[c])[0] is not object
            ]
            X = np.empty((len(rows), len(numeric_cols)), dtype=np.float64)
            for j, c in enumerate(numeric_cols):
                X[:, j] = [
                    float(v) if v is not None else np.nan
                    for v in col_values[c]
                ]
            return X, numeric_cols

        if not ordered:
            return np.empty((len(rows),), dtype=object)

        fields = []
        normalized = {}
        for c in ordered:
            dt, vals = _classify(col_values[c])
            fields.append((c, dt))
            normalized[c] = (dt, vals)

        arr = np.empty(len(rows), dtype=fields)
        for c in ordered:
            dt, vals = normalized[c]
            if dt is object:
                # Assign element-by-element so list/array cells are stored as
                # objects rather than being broadcast into the field.
                col = np.empty(len(rows), dtype=object)
                for i, v in enumerate(vals):
                    col[i] = v
                arr[c] = col
            else:
                arr[c] = np.asarray(vals, dtype=dt)
        return arr

    def filter(self, predicate=None, **kwargs) -> "Database":
        """Filter database entries and return a new database with matching entries.

        Can filter using a predicate function or keyword arguments for common filters.

        Args:
            predicate: A callable that takes an entry and returns True/False.
                       If provided, kwargs are ignored.
            **kwargs: Keyword filters. Supported filters depend on entry type:
                - Z: Filter by atomic number (for nuclide databases)
                - A: Filter by mass number (for nuclide databases)
                - N: Filter by neutron number (for nuclide databases)
                - Any entry attribute: Filter by exact value match

        Returns:
            A new Database instance containing only matching entries

        Examples:
            >>> # Filter by Z (all lead isotopes)
            >>> pb_db = db.filter(Z=82)

            >>> # Filter by predicate (neutron-rich nuclei)
            >>> neutron_rich = db.filter(lambda e: e.n.N > e.n.Z)

            >>> # Filter by multiple criteria
            >>> heavy_pb = db.filter(Z=82, A=lambda a: a > 200)

            >>> # Chain filters
            >>> result = db.filter(Z=82).filter(lambda e: e.Mth < -10)
        """
        # Create a new database of the same type
        result = type(self)()

        for key, entry in self.data.items():
            # If predicate is provided, use it directly
            if predicate is not None:
                if callable(predicate) and predicate(entry):
                    result.data[key] = entry
                continue

            # Otherwise, apply keyword filters
            if self._matches_filters(key, entry, **kwargs):
                result.data[key] = entry

        return result

    def _matches_filters(self, key, entry, **kwargs) -> bool:
        """Check if an entry matches the given keyword filters."""
        for filter_key, filter_value in kwargs.items():
            # Handle nucleus-related filters (Z, A, N)
            if filter_key in ('Z', 'A', 'N'):
                # Get the nucleus from either the key or entry
                nucleus = None
                if hasattr(key, filter_key):
                    nucleus = key
                elif hasattr(entry, 'n') and hasattr(entry.n, filter_key):
                    nucleus = entry.n

                if nucleus is None:
                    return False

                actual_value = getattr(nucleus, filter_key)

                # Support callable filter values
                if callable(filter_value):
                    if not filter_value(actual_value):
                        return False
                elif actual_value != filter_value:
                    return False

            # Handle entry attribute filters
            elif hasattr(entry, filter_key):
                actual_value = getattr(entry, filter_key)

                # Support callable filter values
                if callable(filter_value):
                    if not filter_value(actual_value):
                        return False
                elif actual_value != filter_value:
                    return False

            else:
                # Filter key not found on entry
                return False

        return True

    def filter_by_range(self, attr: str, min_val=None, max_val=None) -> "Database":
        """Filter entries where an attribute falls within a range.

        Args:
            attr: Attribute name to filter on (can use 'Z', 'A', 'N' for nucleus attrs)
            min_val: Minimum value (inclusive). If None, no lower bound.
            max_val: Maximum value (inclusive). If None, no upper bound.

        Returns:
            A new Database instance containing only matching entries

        Examples:
            >>> # Nuclei with Z between 80 and 90
            >>> db.filter_by_range('Z', 80, 90)

            >>> # Entries with mass excess less than 0
            >>> db.filter_by_range('Mth', max_val=0)
        """
        def in_range(entry):
            # Try to get value from nucleus first, then from entry
            value = None
            if attr in ('Z', 'A', 'N') and hasattr(entry, 'n'):
                value = getattr(entry.n, attr, None)
            if value is None:
                value = getattr(entry, attr, None)

            if value is None:
                return False

            if min_val is not None and value < min_val:
                return False
            if max_val is not None and value > max_val:
                return False
            return True

        return self.filter(in_range)

    def find(self, predicate) -> object | None:
        """Find the first entry matching the predicate.

        Args:
            predicate: A callable that takes an entry and returns True/False

        Returns:
            The first matching entry, or None if no match found

        Examples:
            >>> # Find the first nucleus with A > 250
            >>> entry = db.find(lambda e: e.n.A > 250)
        """
        for entry in self.data.values():
            if predicate(entry):
                return entry
        return None

    def count(self, predicate=None, **kwargs) -> int:
        """Count entries matching the filter criteria.

        Args:
            predicate: A callable that takes an entry and returns True/False
            **kwargs: Keyword filters (same as filter())

        Returns:
            Number of matching entries

        Examples:
            >>> # Count lead isotopes
            >>> db.count(Z=82)
            41
        """
        return len(self.filter(predicate, **kwargs).data)

    def __len__(self) -> int:
        """Return the number of entries in the database."""
        return len(self.data)

    def __contains__(self, key) -> bool:
        """Check if a key exists in the database."""
        return key in self.data

    def contains(self, key) -> bool:
        """Check if a key exists in the database.

        This is an explicit method version of the `in` operator.

        Args:
            key: The key to check (typically a Nuclide or element number)

        Returns:
            True if the key exists in the database, False otherwise

        Examples:
            >>> n = Nuclide(Z=82, A=208)
            >>> db.contains(n)
            True
            >>> n in db  # Equivalent using `in` operator
            True
        """
        return key in self.data

    def __iter__(self, ) -> iter:
        """Return the list of keys for the data property as an iterable. """
        return iter(self.data.keys())


class NuclideDatabase(Database):
    """Data in this type of RIPL database are accessible via a single nucleus. """

    @property
    def nuclei(self, ) -> list:
        """Return the list of nuclei in the database. """
        return list(self.data.keys())

    def insert_by_entry(self, entry: object) -> None:
        """Insert a new entry into the database. """
        # Check that a proper entry has been provided
        if isinstance(entry, self.entry):
            # Insert the entry by its nucleus; this may overwrite an existing entry
            self.data[entry.n] = entry
        else:
            raise TypeError(f"Entry must be of type = ({type(self.entry)}); given = ({type(entry)})!")

    def remove_by_entry(self, entry: object) -> None:
        """Remove an existing entry in the database. """
        # Check that a proper entry has been provided
        if hasattr(entry, 'n'):
            try:
                self.data.pop(entry.n)
            except KeyError:
                raise NucleusNotFoundError(f"Nuclide {entry.n} not found in database.")
        else:
            raise TypeError("Entry object must have a nucleus property!")


class ElementDatabase(Database):
    """Data  in this type of RIPL database are accessible through an element number. """

    @property
    def elements(self, ) -> list:
        """Return the list of elements in the database. """
        return list(self.data.keys())

    def insert_by_entry(self, entry: object) -> None:
        """Insert a new entry into the database. """
        # Check that a proper entry has been provided
        if isinstance(entry, self.entry):
            # Insert the entry by its nucleus; this may overwrite an existing entry
            self.data[entry.Z] = entry
        else:
            raise TypeError(f"Entry must be of type = ({type(self.entry)}); given = ({type(entry)})!")

    def remove_by_entry(self, entry: object) -> None:
        """Remove an existing entry, which is keyed by element number (Z).

        The entry's ``Z`` is used when present; otherwise its element symbol
        (``sym``/``symbol``) or ``name`` is resolved to a ``Z``. Raises
        ElementNotFoundError if the element is not in the database (or is not a
        recognized symbol/name), matching :meth:`get`'s contract.
        """
        # Resolve the element key (Z) from whatever identifier the entry carries.
        if getattr(entry, 'Z', None) is not None:
            z = entry.Z
            ident = entry.Z
        else:
            from riplpy.elements import Elements
            sym = getattr(entry, 'sym', None) or getattr(entry, 'symbol', None)
            name = getattr(entry, 'name', None)
            if sym is not None:
                ident = sym
                z = Elements.SymtoZ.get(str(sym).capitalize())
            elif name is not None:
                ident = name
                z = Elements.NametoZ.get(str(name).capitalize())
            else:
                raise TypeError("Entry object must have an element number, symbol, or name!")
            if z is None:
                raise ElementNotFoundError(f"Element '{ident}' is not a recognized symbol or name.")

        try:
            self.data.pop(z)
        except KeyError:
            raise ElementNotFoundError(f"Element {ident} not found in database.")


class _DbAccessor:
    """A pass-through object which enables dynamic Database access.

    Canonical databases are stored as plain attributes. Shorthand/longhand
    aliases (e.g. ``frdm12`` -> ``frdm2012``) can be registered with
    :meth:`add_alias`; an alias resolves to its canonical attribute on access
    (``db.frdm12``) but is NOT listed by iteration, so ``list_databases`` keeps
    showing one canonical name per database.
    """

    def __init__(self, ) -> None:
        """Store dynamically assigned attributes and aliases separately. """
        super().__setattr__("_attributes", set())
        super().__setattr__("_aliases", {})

    def __setattr__(self, name, value):
        """Set attribute so long as it isn't protected. """
        if name not in ("_attributes", "_aliases") and not name.startswith("__"):
            self._attributes.add(name)
        super().__setattr__(name, value)

    def __getattr__(self, name):
        """Resolve a registered alias to its canonical attribute.

        Only invoked when normal attribute lookup fails, so canonical
        databases never reach here.
        """
        aliases = self.__dict__.get("_aliases", {})
        if name in aliases:
            return getattr(self, aliases[name])
        raise AttributeError(
            f"{type(self).__name__!r} object has no attribute {name!r}"
        )

    def add_alias(self, alias: str, canonical: str) -> None:
        """Register ``alias`` as another name for the ``canonical`` database. """
        self._aliases[alias] = canonical

    def resolve(self, name: str) -> str:
        """Return the canonical attribute name for a possibly-aliased name. """
        return self.__dict__.get("_aliases", {}).get(name, name)

    def has(self, name: str) -> bool:
        """True if ``name`` is a known database (canonical or alias). """
        return self.resolve(name) in self._attributes

    def __call__(self):
        """Iterate over assigned (canonical) attributes. """
        for attr_name in self._attributes:
            yield attr_name


def load(db_obj: Database, db_file_path: str) -> "Database":
    """Load the database into memory from a file path and return it. """
    db = db_obj()
    db.load(db_file_path)
    return db


def loader(directory: str = None, file_path: str = None, config: object = None, db_obj: object = None, file_path_key: str = None) -> "Database":
    """Load the database into memory from file path, directory, or configuration and return it.

    This is a more complex loader method. 
    The loading mechanism prioritizes the `file_path` argument first. If `file_path`
    is not provided, it attempts to construct the path using the `directory` argument.
    If neither is provided, it falls back to the default path from the configuration.

    Arguments:
        directory (str, optional): Path to the directory containing the data file.
        file_path (str, optional): Direct path to the data file.
        config (object, optional): Configuration object containing the data file path.
        database (object, optional): Database object to load the data into.
        file_path_key (str, optional): Key for the data file path in the configuration.
    Returns:
        Database: An instance of the Database class with the loaded data.

    Raises:
        FileFormatError: If the database file cannot be loaded.
        RiplFileNotFoundError: If the database file cannot be found.
    """
    db_path = None
    if file_path:
        db_path = file_path
    else:
        # Resolve the directory consistently with riplpy.load() and the
        # section loaders (set_path / RIPL_LOCATION / ~/.riplpyrc / auto).
        # ``DATA_FILES`` holds the RIPL-relative path; join it with the
        # resolved directory. (Note: config.get_data_file_path() pre-joins
        # the import-time RIPL_PATH and is therefore unsuitable here — it
        # would silently ignore an explicit ``directory`` argument.)
        try:
            resolved = config.resolve_directory(directory)
        except (ValueError, ImportError):
            resolved = None
        rel = None
        try:
            rel = config.DATA_FILES.get(file_path_key)
        except (AttributeError, KeyError):
            rel = None
        if resolved and rel:
            db_path = _os.path.join(resolved, rel)

    if db_path and _os.path.exists(db_path):
        try:
            return load(db_obj, db_path)
        except Exception as e:
            raise FileFormatError(f"Failed to load {file_path_key} data from {db_path}: {e}")

    raise RiplFileNotFoundError(f"Unable to find or load {file_path_key} data file. Please provide a valid path or set RIPL_LOCATION.")
