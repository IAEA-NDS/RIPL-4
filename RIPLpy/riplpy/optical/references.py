# -*- coding: utf-8 -*-
"""Parser for optical model references (om-references.txt).

This module provides access to the bibliographic references for the
optical model potentials in the RIPL library.
"""

import os as _os
import re
from dataclasses import dataclass as _dataclass
from typing import Dict, Optional, ClassVar as _ClassVar

from riplpy.db import Database
from . import config


@_dataclass
class Reference:
    """A bibliographic reference for an optical model potential.

    Attributes:
        ref_num: Reference number in the library
        citation: Full citation text
    """
    ref_num: int = 0
    citation: str = ''

    _field_info: _ClassVar[dict] = {
        'ref_num': 'Reference number in the library',
        'citation': 'Full bibliographic citation',
    }

    @property
    def first_author(self) -> str:
        """Extract the first author name from the citation."""
        # Try to get first author before first comma or 'and'
        text = self.citation.strip()
        if not text:
            return ''

        # Common patterns: "A.B.Name," or "A.B.Name and" or "A.B.Name et al."
        match = re.match(r'^([A-Z][a-zA-Z.\-\s]+?)(?:,|\s+and\s|\s+et\s+al)', text)
        if match:
            return match.group(1).strip()

        # Fallback: take first word-like chunk
        parts = text.split(',')[0].split(' and ')[0]
        return parts.strip()

    @property
    def year(self) -> Optional[int]:
        """Extract publication year from the citation."""
        # Look for 4-digit year in parentheses or standalone
        match = re.search(r'\((\d{4})\)', self.citation)
        if match:
            return int(match.group(1))

        # Try standalone year
        match = re.search(r'\b(19\d{2}|20\d{2})\b', self.citation)
        if match:
            return int(match.group(1))

        return None

    @property
    def journal(self) -> str:
        """Try to extract journal name from citation."""
        # Common journal patterns
        patterns = [
            r'Phys\.\s*Rev\.',
            r'Nucl\.\s*Phys\.',
            r'Nucl\.\s*Sci\.\s*Eng\.',
            r'Phys\.\s*Lett\.',
            r'J\.\s*Phys\.',
            r'At\.\s*Data',
        ]

        for pattern in patterns:
            if re.search(pattern, self.citation):
                match = re.search(pattern + r'[^,]+', self.citation)
                if match:
                    return match.group(0).strip()

        return ''

    def __repr__(self) -> str:
        short_cite = self.citation[:50] + '...' if len(self.citation) > 50 else self.citation
        return f"Reference({self.ref_num}: {short_cite})"


class ReferenceDatabase(Database):
    """Database of bibliographic references for optical model potentials."""

    entry = Reference

    def __init__(self, data: Dict[int, Reference] = None):
        """Initialize the reference database.

        Args:
            data: Optional dictionary mapping ref_num -> Reference
        """
        super().__init__(data)

    def __repr__(self) -> str:
        return f"<ReferenceDatabase: {len(self.data)} references>"

    def get(self, ref_num: int) -> Reference:
        """Get a reference by its number.

        Args:
            ref_num: Reference number

        Returns:
            Reference for the specified number

        Raises:
            EntryNotFoundError: If reference not found (a KeyError subclass,
                so ``except KeyError`` still works).
        """
        if ref_num not in self.data:
            from riplpy.exceptions import EntryNotFoundError
            raise EntryNotFoundError(f"Reference {ref_num} not found")
        return self.data[ref_num]

    def search(self, query: str) -> "ReferenceDatabase":
        """Search references by text.

        Args:
            query: Search string (case-insensitive)

        Returns:
            New ReferenceDatabase with matching references
        """
        query_lower = query.lower()
        filtered = {
            k: v for k, v in self.data.items()
            if query_lower in v.citation.lower()
        }
        return ReferenceDatabase(filtered)

    def search_by_author(self, author: str) -> "ReferenceDatabase":
        """Search references by author name.

        Args:
            author: Author name or partial name

        Returns:
            New ReferenceDatabase with matching references
        """
        author_lower = author.lower()
        filtered = {
            k: v for k, v in self.data.items()
            if author_lower in v.citation.lower().split(',')[0]
        }
        return ReferenceDatabase(filtered)

    def search_by_year(self, year: int) -> "ReferenceDatabase":
        """Search references by publication year.

        Args:
            year: Publication year

        Returns:
            New ReferenceDatabase with matching references
        """
        filtered = {
            k: v for k, v in self.data.items()
            if v.year == year
        }
        return ReferenceDatabase(filtered)

    def info(self) -> str:
        """Return detailed information about this database."""
        lines = ["Reference Database"]
        lines.append("-" * 40)
        lines.append(f"Total references: {len(self.data)}")

        # Year range
        years = [v.year for v in self.data.values() if v.year]
        if years:
            lines.append(f"Year range: {min(years)} - {max(years)}")

        return "\n".join(lines)


def read_references(filepath: str) -> ReferenceDatabase:
    """Read the optical model references file.

    The file format has reference numbers followed by citation text.
    Citations may span multiple lines (continuation lines are indented).

    Args:
        filepath: Path to om-references.txt

    Returns:
        ReferenceDatabase with all references

    Raises:
        FileNotFoundError: If the file doesn't exist
    """
    if not _os.path.exists(filepath):
        raise FileNotFoundError(f"References file not found: {filepath}")

    data = {}
    current_ref_num = None
    current_citation = []

    with open(filepath, 'r', encoding='utf-8', errors='replace') as fp:
        for line in fp:
            line = line.replace('\r', '').rstrip()

            # Skip header lines
            if 'REFERENCES' in line or not line.strip():
                continue

            # Check if this is a new reference (starts with number)
            match = re.match(r'^\s*(\d+)\.\s+(.*)$', line)
            if match:
                # Save previous reference
                if current_ref_num is not None:
                    citation = ' '.join(current_citation).strip()
                    data[current_ref_num] = Reference(
                        ref_num=current_ref_num,
                        citation=citation
                    )

                # Start new reference
                current_ref_num = int(match.group(1))
                current_citation = [match.group(2)]
            elif current_ref_num is not None and line.strip():
                # Continuation line
                current_citation.append(line.strip())

    # Save last reference
    if current_ref_num is not None:
        citation = ' '.join(current_citation).strip()
        data[current_ref_num] = Reference(
            ref_num=current_ref_num,
            citation=citation
        )

    return ReferenceDatabase(data)


def load(directory: str = None, file_path: str = None) -> ReferenceDatabase:
    """Load the optical model references database.

    Args:
        directory: RIPL directory path. If provided, constructs full path.
        file_path: Direct path to references file. Takes precedence over directory.

    Returns:
        ReferenceDatabase with loaded references

    Raises:
        FileNotFoundError: If the references file cannot be found
    """
    if file_path:
        return read_references(file_path)

    if directory:
        path = _os.path.join(directory, config.get_data_file_path('references'))
        return read_references(path)

    # Try default RIPL path
    from riplpy import config as riplpy_config
    ripl_path = riplpy_config.get_path()
    if ripl_path:
        path = _os.path.join(ripl_path, config.get_data_file_path('references'))
        return read_references(path)

    raise FileNotFoundError(
        "Cannot locate references file. Please provide directory or file_path."
    )
