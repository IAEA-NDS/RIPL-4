# -*- coding: utf-8 -*-
"""Custom exceptions for the RIPLpy library."""


class RIPLPyError(Exception):
    """Base class for exceptions in RIPLpy."""
    pass


class FileFormatError(RIPLPyError):
    """Raised when a file is not in the expected format."""
    pass


class RiplFileNotFoundError(RIPLPyError, FileNotFoundError):
    """Raised when a RIPL file is not found on the filesystem."""
    pass


class NucleusNotFoundError(RIPLPyError, KeyError):
    """Raised when a nucleus is not found in a database."""
    pass


class ElementNotFoundError(RIPLPyError, KeyError):
    """Raised when an element is not found in a database."""
    pass


class EntryNotFoundError(RIPLPyError, KeyError):
    """Raised when a non-nuclide-keyed entry is not found in a database.

    Used for lookups keyed by an identifier other than a nucleus/element
    (e.g. an optical-model reference number ``iref`` or a bibliographic
    reference number). Subclasses ``KeyError`` so existing ``except KeyError``
    handlers keep working.
    """
    pass


class MissingDataWarning(UserWarning):
    """Issued when databases are left empty because their data files are absent.

    A subclass of ``UserWarning`` so it is shown by default (unlike a logging
    message, which needs logging to be configured) and can be filtered or
    turned into an error with the standard ``warnings`` machinery.
    """
    pass
