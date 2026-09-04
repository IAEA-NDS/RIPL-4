# -*- coding: utf-8 -*-
"""Low-level parsing utilities for optical model parameter files.

This module provides functions for parsing the Fortran free-format data
used in the RIPL optical model files, including handling of Fortran-style
floating point notation (e.g., '1.00000-3' for 1.0e-3).
"""

import re
from typing import List, Tuple, Iterator, TextIO


def parse_fortran_float(s: str) -> float:
    """Parse a Fortran-style floating point number.

    Handles formats like:
        '1.23456'      -> 1.23456
        '1.00000-3'    -> 1.0e-3
        '1.00000+3'    -> 1.0e+3
        '-1.00000-3'   -> -1.0e-3
        '.00000'       -> 0.0
        '-.5'          -> -0.5

    Args:
        s: String representation of the number

    Returns:
        Floating point value

    Raises:
        ValueError: If the string cannot be parsed
    """
    s = s.strip()
    if not s:
        return 0.0

    # Check for Fortran-style exponent notation (e.g., '1.00000-3' or '1.00000+3')
    # This is different from standard notation because there's no 'e' or 'E'
    # Pattern: optional sign, digits, decimal, digits, then +/- and exponent digits
    # But we need to be careful: '-1.5' is not the same as '1.5-1'

    # First, try standard float parsing
    try:
        return float(s)
    except ValueError:
        pass

    # Try Fortran notation: look for +/- that's NOT at the start
    # Pattern: number followed by +/- followed by digits at end
    match = re.match(r'^([+-]?\d*\.?\d+)([+-])(\d+)$', s)
    if match:
        mantissa = float(match.group(1))
        sign = match.group(2)
        exponent = int(match.group(3))
        if sign == '-':
            exponent = -exponent
        return mantissa * (10 ** exponent)

    raise ValueError(f"Cannot parse Fortran float: '{s}'")


def tokenize_line(line: str) -> List[str]:
    """Split a line into whitespace-separated tokens.

    Args:
        line: Input line

    Returns:
        List of tokens
    """
    return line.split()


def read_floats_from_line(line: str) -> List[float]:
    """Read all floating point values from a line.

    Args:
        line: Input line with whitespace-separated values

    Returns:
        List of float values
    """
    tokens = tokenize_line(line)
    return [parse_fortran_float(t) for t in tokens]


def read_ints_from_line(line: str) -> List[int]:
    """Read all integer values from a line.

    Args:
        line: Input line with whitespace-separated values

    Returns:
        List of integer values
    """
    tokens = tokenize_line(line)
    return [int(t) for t in tokens]


def read_coefficient_array(lines_iter: Iterator[str], n_values: int, values_per_line: int = 7) -> List[float]:
    """Read a coefficient array that spans multiple lines.

    The RIPL format stores coefficient arrays with a fixed number of values
    per line (typically 7 for the first line, then continuation lines with
    leading whitespace).

    Args:
        lines_iter: Iterator over input lines
        n_values: Total number of values to read
        values_per_line: Expected values per line (default 7)

    Returns:
        List of coefficient values

    Example:
        For 13 coefficients:
        Line 1:  1.26000   .00000+0   .00000+0   .00000+0   .00000+0   .00000+0   .00000+0
        Line 2:            .00000+0   .00000+0   .00000+0   .00000+0   .00000+0   .00000+0
    """
    values = []
    while len(values) < n_values:
        line = next(lines_iter).replace('\r', '')
        line_values = read_floats_from_line(line)
        values.extend(line_values)

    # Truncate to exactly n_values (in case of over-read)
    return values[:n_values]


def read_potential_block(lines_iter: Iterator[str]) -> Tuple[float, List[float], List[float], List[float]]:
    """Read a single potential energy range block.

    Each block consists of:
        - epot (upper energy limit)
        - rco[13] (radius coefficients, 2 lines)
        - aco[13] (diffuseness coefficients, 2 lines)
        - pot[25] (potential strength, 4 lines)

    Args:
        lines_iter: Iterator over input lines

    Returns:
        Tuple of (epot, rco_values, aco_values, pot_values)
    """
    # Read epot (single value on its own line)
    epot_line = next(lines_iter).replace('\r', '')
    epot = parse_fortran_float(epot_line.strip())

    # Read rco[13] - 2 lines (7 + 6 values)
    rco = read_coefficient_array(lines_iter, 13)

    # Read aco[13] - 2 lines (7 + 6 values)
    aco = read_coefficient_array(lines_iter, 13)

    # Read pot[25] - 4 lines (7 + 7 + 7 + 4 values)
    pot = read_coefficient_array(lines_iter, 25)

    return epot, rco, aco, pot


def read_text_line(lines_iter: Iterator[str]) -> str:
    """Read a single text line, stripping whitespace.

    Args:
        lines_iter: Iterator over input lines

    Returns:
        Stripped line content
    """
    return next(lines_iter).replace('\r', '').strip()


def read_text_lines(lines_iter: Iterator[str], n_lines: int) -> List[str]:
    """Read multiple text lines.

    Args:
        lines_iter: Iterator over input lines
        n_lines: Number of lines to read

    Returns:
        List of stripped line contents
    """
    return [read_text_line(lines_iter) for _ in range(n_lines)]


def peek_line(fp: TextIO) -> str:
    """Peek at the next line without consuming it.

    Args:
        fp: File object

    Returns:
        Next line content (or empty string at EOF)
    """
    pos = fp.tell()
    line = fp.readline()
    fp.seek(pos)
    return line


def skip_to_separator(lines_iter: Iterator[str], separator: str = '+' * 80) -> bool:
    """Skip lines until a separator line is found.

    Args:
        lines_iter: Iterator over input lines
        separator: Separator string to look for

    Returns:
        True if separator was found, False if EOF reached
    """
    try:
        while True:
            line = next(lines_iter).replace('\r', '')
            if separator in line:
                return True
    except StopIteration:
        return False


def is_separator_line(line: str, separator: str = '+' * 40) -> bool:
    """Check if a line is a separator line.

    Args:
        line: Line to check
        separator: Separator pattern (default: 40+ plus signs)

    Returns:
        True if line contains separator
    """
    return separator in line.replace('\r', '')


class OMPFileReader:
    """Context manager for reading OMP parameter files.

    Provides convenient iteration over potential entries in the file,
    handling the separator-delimited format.

    Example:
        with OMPFileReader(filepath) as reader:
            for entry_lines in reader.entries():
                # Process each potential entry
                pass
    """

    def __init__(self, filepath: str):
        """Initialize the reader.

        Args:
            filepath: Path to the OMP parameter file
        """
        self.filepath = filepath
        self._fp = None

    def __enter__(self):
        self._fp = open(self.filepath, 'r', encoding='utf-8', errors='replace')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._fp:
            self._fp.close()
        return False

    def entries(self) -> Iterator[List[str]]:
        """Iterate over potential entries in the file.

        Each entry is a list of lines between separators.

        Yields:
            List of lines for each potential entry
        """
        current_entry = []
        separator = '+' * 40

        for line in self._fp:
            line = line.replace('\r', '')
            if separator in line:
                if current_entry:
                    yield current_entry
                    current_entry = []
            else:
                current_entry.append(line)

        # Don't forget the last entry if file doesn't end with separator
        if current_entry:
            yield current_entry

    def reset(self):
        """Reset file position to beginning."""
        if self._fp:
            self._fp.seek(0)
