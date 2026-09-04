# -*- coding: utf-8 -*-
"""Python objects that provide access to commonly used quantities that can be used in or as collections.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Create a placeholder representation of a nucleus as a Nucleus object::
       $ n = Nucleus()
    (2) Create a representation of 132-Sn as a Nucleus object::
       $ n = Nucleus(Z=50, N=82)
    (3) Create a representation of 132-Sn as a Nuclide object (which accepts more attributes upon init)::
       $ n = Nuclide(Z=50, N=82, attr='this is new')
    (4) Initialize a Nuclide via element name and mass number::
       $ n = Nuclide(sym='Sn', A=132)
    (5) Check the equivalence of a Nucleus and Nuclide object (this comparison only depends on (Z,N) - thus returns True)::
       $ n = Nucleus(Z=50, N=82)
       $ m = Nuclide(Z=50, N=82)
       $ n == m # True
    (6) Print out the symbol name of the Nucleus::
       $ n = Nucleus(Z=48, N=130)
       $ print(n.element_symbol)
    (7) Access the mass number of a Nucleus object::
       $ n = Nucleus(Z=50, N=82)
       $ n.A
    (8) Create a representation of the ground state of 132-Sn as an IsomericNucleus object which contains state information in the 'm' attribute (m=0 by default)::
       $ n = IsomericNucleus(Z=50, N=82)
    (9) Create a representation of the second excited state of 132-Sn::
       $ n = IsomericNucleus(Z=50, N=82, m=2)
    (10) Check the equivalence of a Nucleus and an IsomericNucleus object (these two are not the same)::
       $ n = Nucleus(Z=82, N=126)
       $ m = IsomericNucleus(Z=82, N=126, m=10)
       $ m == n # False
    (11) Check the equivalence of a Nucleus and an IsomericNucleus object when m == 0::
       $ n = Nucleus(Z=82, N=126)
       $ m = IsomericNucleus(Z=82, N=126)
       $ m == n # True
    (12) Create a set of Nucleus-like objects::
       $ n, m, p, q = Nucleus(Z=82, N=126), Nucleus(Z=50, N=82), IsomericNucleus(Z=82, N=126), Nucleus(Z=71, N=111)
       $ s = set([n,m,p,q]) # only contains: n, m and q (since n == p)
    (13) Create a set of Nucleus-like objects::
       $ n, m, p, q = Nucleus(Z=82, N=126), Nucleus(Z=50, N=82), IsomericNucleus(Z=82, N=126, m=2), Nucleus(Z=71, N=111)
       $ s = set([n,m,p,q]) # contains all since the isomeric nucleus (p) is distinct from the ground state (n)

"""

# ========================

# Collections
from collections import namedtuple as _namedtuple

# Dataclass
from dataclasses import dataclass as _dataclass

# RegEx
import re as _re

# Local Namespace
from .elements import Elements

# ========================

# Tuples that can be used to quickly access properties
ZNTuple              = _namedtuple('ZNTuple', 'Z, N')
ZATuple              = _namedtuple('ZATuple', 'Z, A')
NATuple              = _namedtuple('NATuple', 'N, A')
NucleusTuple         = _namedtuple('NucleusTuple', 'Z, N, A')
IsomericNucleusTuple = _namedtuple('IsomericNucleusTuple', 'Z, N, A, i')
SymNucleusTuple      = _namedtuple('SymNucleusTuple', 'Z, N, A, sym')
NamedNucleusTuple    = _namedtuple('NamedNucleusTuple', 'Z, N, A, name')
SymNamedNucleusTuple = _namedtuple('SymNamedNucleusTuple', 'Z, N, A, sym, name')


class Base(object):
    """A generic base class for collections. """

    def _subset(self, kwargs: dict) -> tuple:
        """Return a tuple of two sets: 
           1. The set of keys in kwargs that are part of the reserved kwargs.
           2. The set of keys in kwargs that are not part of the reserved kwargs.
        """
        reserved_keys = set()
        other_keys = set()
    
        for key in kwargs.keys():
            k = key.lower()  # Force to lower case to check against reserved dictionary keys
            if k in self._reserved_kwargs:
                reserved_keys.add(k)
            else:
                other_keys.add(key)
    
        return reserved_keys, other_keys


class NucleusBase(Base):
    """A base class for Nucleus-like object. 
    
    This class should not be initialized directly. Instead it contains useful methods for extended objects.
    This class assumes that several variables exist:
        cls.Z - Proton number
        cls.N - Neutron number
        cls.A - Mass number
    """

    @property
    def element_symbol(self, ) -> str:
        """Return the element symbol as a string. """
        return Elements.symbol(self.Z)
    
    @property
    def element_name(self, ) -> str:
        """Return the name of the element as a string."""
        return Elements.name(self.Z)

    @property
    def sym(self, ) -> str:
        return self.element_symbol

    @property
    def symbol(self, ) -> str:
        return self.element_symbol

    def __eq__(self, other) -> bool:
        """Return the equality (bool) of this object with another. """
        return self.Z == other.Z and self.N == other.N

    def __ne__(self, other) -> bool:
        """Return the not-equality (bool) of this object with another. """
        return not self.__eq__(other)

    def __hash__(self) -> int:
        """Return a unique (unsafe) hash of a NucleusBase object given the independent variables: (Z,N). 
        
        By default this also assumes a third component (m=0); the ground state.
        """
        return hash((self.Z, self.N, 0))

    def __str__(self, ) -> str:
        """Return a string representation of the object. """
        return f"<{self.__class__.__name__}(z={self.Z}, n={self.N})>"


class Element(Base):
    """A class representing a single element. """

    _reserved_kwargs: tuple = ('z', 'symbol', 'sym', 'name')

    def __init__(self, **kwargs) -> None:
        """Initialize the object. """

        # If no kwargs are provided, return immediately (placeholder behavior)
        if not kwargs:
            return

        # Partition the keys into reserved and other keys
        reserved_keys, other_keys = self._subset(kwargs)

        # Map of lowercase reserved keys to their original case-sensitive counterparts
        reserved_map = {k.lower(): k for k in kwargs.keys() if k.lower() in reserved_keys}

        # Retrieve reserved parameters
        z = int(kwargs.get(reserved_map.get('z'), 0))
        sym = kwargs.get(reserved_map.get('sym'))
        symbol = kwargs.get(reserved_map.get('symbol'))
        name = kwargs.get(reserved_map.get('name'))

        # Determine Z based on provided input
        if z:
            self.Z = z
        elif sym:
            self.Z = Elements.element(sym)
        elif symbol:
            self.Z = Elements.element(symbol)
        elif name:
            self.Z = Elements.element(name)
        else:
            raise Exception("Element is not valid: requires 'z', 'sym', 'symbol', or 'name'.")

        # Set any additional attributes from remaining kwargs
        for k in other_keys:
            setattr(self, k, kwargs[k])

    @property
    def symbol(self, ) -> str:
        """Return the element symbol. """
        return Elements.symbol(self.Z)
    
    @property
    def sym(self, ) -> str:
        """Return the element symbol. """
        return self.symbol

    @property
    def name(self, ) -> str:
        """Return the element name. """
        return Elements.name(self.Z)

    def __lt__(self, other: object) -> bool:
        """Return true if self is strictly smaller than the given (other) object. """
        if self.Z < other.Z:
            return True
        else:
            return False

    def __eq__(self, other) -> bool:
        """Return the equality (bool) of this object with another. """
        return self.Z == other.Z

    def __ne__(self, other) -> bool:
        """Return the equality (bool) of this object with another. """
        return not self.__eq__(other)

    def __hash__(self, ) -> int:
        """Return a unique (unsafe) hash of an Element object given the independent variables, Z. """
        return hash(self.Z)

    def __str__(self, ) -> str:
        """Return a string representation of the object.

        Falls back to a symbol-less form when ``Z`` is outside the element
        table, so formatting never crashes on an out-of-range value.
        """
        try:
            return f"<{self.__class__.__name__}(z={self.Z}, symbol={self.symbol}, name={self.name})>"
        except LookupError:
            return f"<{self.__class__.__name__}(z={self.Z})>"

    def __repr__(self, ) -> str:
        """Return a string representation of the object. """
        return f"{self.__class__.__name__}(z={self.Z})"


@_dataclass
class Nucleus(NucleusBase):
    """An object representing a nucleus defined by Z-protons and N-neutrons: (Z,N). 

    By default initialization, this object can be used as a placeholder with no (Z,N) provided.
    Otherwise, this object is intended to represent a nucleus with Z protons and N neutrons.
    Besides a simple check upon initialization, this light-weight object offers no further checking methods.
    Caution: be very careful when using hashes with this object, as it is technically mutable. 
    """

    Z: int = None # The proton number of the nucleus
    N: int = None # The neutron number of the nucleus

    @classmethod
    def from_ZA(cls, val: tuple):
        """Initialize the class from (Z,A) tuple. """
        Z, A = val
        N = A - Z
        return cls(Z,N)

    @classmethod
    def from_NA(cls, val: tuple):
        """Initialize the object from (N,A) tuple. """
        N, A = val
        Z = A - N
        return cls(Z,N)

    def __post_init__(self, ) -> None:
        """Perform post-initialization checks. """
        if self.Z is not None and self.N is not None:
            assert self.Z + self.N == self.A

    @property
    def A(self, ) -> int:
        """Return the mass number of the Nucleus. """
        return self.Z + self.N

    @A.setter
    def A(self, val: tuple) -> None:
        """Set the (dependent) mass number given both proton and neutron numbers (Z,N) as a tuple. """
        Z,N = val
        # We provide a bit of validation, mapping the input to an integer
        self.Z = int(Z)
        self.N = int(N)

    def __hash__(self) -> int:
        """Return a unique (unsafe) hash of a NucleusBase object given the independent variables: (Z,N). 
        
        In this object we have to redefine the (unsafe) hash due to it being built of a dataclass.
        By default this also assumes a third component (m=0); the ground state.
        """
        return hash((self.Z, self.N, 0))

    def __repr__(self, ) -> str:
        """Return a string representation of the object which can be used to initialize. """
        return f"Nucleus(z={self.Z}, n={self.N})"


class Nuclide(NucleusBase):
    """A Nucleus-like object which supports a more robust initialization procedure where additional attributes may be set.

       Because additional attributes may be set upon initialization, this object may perform more slowly than a Nucleus object.
       Caution: be very careful when using hashes with this object, as it is technically mutable. 
       Caution: the equivalence of this object with another depends only on the defining attributes: (Z,N).
    """

    _reserved_kwargs: tuple = ('z', 'n', 'a', 'symbol', 'sym', 'name')


    def __init__(self, **kwargs) -> None:
        """Initialize the object. """

        # If no kwargs are provided, return immediately (placeholder behavior)
        if not kwargs:
            return

        # Partition the keys into reserved and other keys
        reserved_keys, other_keys = self._subset(kwargs)

        # Map of lowercase reserved keys to their original case-sensitive counterparts
        reserved_map = {k.lower(): k for k in kwargs.keys() if k.lower() in reserved_keys}

        # Retrieve reserved parameters, default to None if not present
        z = kwargs.get(reserved_map.get('z'))
        n = kwargs.get(reserved_map.get('n'))
        a = kwargs.get(reserved_map.get('a'))
        symbol = kwargs.get(reserved_map.get('symbol'))
        sym = kwargs.get(reserved_map.get('sym'))
        name = kwargs.get(reserved_map.get('name'))

        # Convert Z, N, A to integers if they are not None
        z = int(z) if z is not None else None
        n = int(n) if n is not None else None
        a = int(a) if a is not None else None

        # Calculate Z and N based on input combinations
        if z is not None and a is not None:
            self.Z = z
            self.N = a - z
        elif z is not None and n is not None:
            self.Z = z
            self.N = n
        elif n is not None and a is not None:
            self.N = n
            self.Z = a - n
        elif sym is not None and a is not None:
            self.Z = Elements.element(sym)
            self.N = a - self.Z
        elif sym is not None and n is not None:
            self.Z = Elements.element(sym)
            self.N = n
        elif symbol is not None and a is not None:
            self.Z = Elements.element(symbol)
            self.N = a - self.Z
        elif symbol is not None and n is not None:
            self.Z = Elements.element(symbol)
            self.N = n
        elif name is not None and a is not None:
            self.Z = Elements.element(name)
            self.N = a - self.Z
        elif name is not None and n is not None:
            self.Z = Elements.element(name)
            self.N = n
        else:
            raise Exception("Nucleus is not valid: requires (Z, A), (Z, N), or (N, A).")

        # Sanity check on nucleus
        assert self.Z + self.N == self.A, "Sanity check failed: Z + N != A."

        # Set any additional attributes from remaining kwargs
        for k in other_keys:
            setattr(self, k, kwargs[k])

    @property
    def A(self, ) -> int:
        """Return the mass number of the Nuclide. """
        return self.Z + self.N

    @classmethod
    def from_string(cls, val: str):
        """Initialize the object from an isotope string.

        Accepts the common notations regardless of separator or element/mass
        order, e.g. ``'56Fe'``, ``'Fe-56'``, ``'Fe56'``, ``'Fe 56'``, and the
        full-name forms ``'Iron-56'`` / ``'56 iron'``.

        Raises:
            ValueError: If the string cannot be parsed into an element and a
                mass number, or the element is not recognized.
        """
        from .elements import isotope_str_to_tuple
        Z, N, A = isotope_str_to_tuple(val)
        return cls(Z=Z, N=N)

    def _subset(self, kwargs: dict) -> tuple:
        """Return a tuple of two sets: 
           1. The set of keys in kwargs that are part of the reserved kwargs.
           2. The set of keys in kwargs that are not part of the reserved kwargs.
        """
        reserved_keys = set()
        other_keys = set()
    
        for key in kwargs.keys():
            k = key.lower()  # Force to lower case to check against reserved dictionary keys
            if k in self._reserved_kwargs:
                reserved_keys.add(k)
            else:
                other_keys.add(key)
    
        return reserved_keys, other_keys

    def __lt__(self, other: object) -> bool:
        """Return true if self is strictly smaller than the given (other) object. """
        if (self.Z, self.A) < (other.Z, other.A):
            return True
        else:
            return False

    def __str__(self, ) -> str:
        """Return a string representation of the object.

        Falls back to a symbol-less form when ``Z`` is outside the element
        table, so formatting never crashes on an out-of-range value.
        """
        try:
            return f"<Nuclide(z={self.Z}, a={self.A}, symbol={self.element_symbol})>"
        except LookupError:
            return f"<Nuclide(z={self.Z}, a={self.A})>"

    def __repr__(self, ) -> str:
        """Return a string representation of the object which can be used to initialize. """
        return f"Nuclide(z={self.Z}, a={self.A})"


@_dataclass
class IsomericNucleus(Nucleus, NucleusBase):
    """An object representing a nucleus that also contains information about an isomeric state.
       
       The object may be initialized by providing proton number (Z), neutron number (N), and isomeric state enumeration (m).

       Caution: the equivalence of this object with another depends only on the defining attributes: (Z,N,m).
    """

    m: int = 0 # By default the state is assumed to be the ground state


    def __post_init__(self, ) -> None:
        """Perform post-initialization checks. """
        super(IsomericNucleus, self).__post_init__()
        # We insist that states are enumerated by the natural numbers
        assert self.m > -1

    def __eq__(self, other) -> bool:
        """Return the equality (bool) of this object with another. """
        test_Z = self.Z == other.Z
        test_N = self.N == other.N
        try:
            m = getattr(other, 'm')
        except AttributeError:
            # Assume ground state
            m = 0
        test_m = self.m == m

        return test_Z and test_N and test_m 

    def __hash__(self) -> int:
        """Return a unique (unsafe) hash of this object. """
        return hash((self.Z, self.N, self.m))

    def __str__(self, ) -> str:
        """Return a string representation of the object.

        Falls back to a symbol-less form when ``Z`` is outside the element
        table, so formatting never crashes on an out-of-range value.
        """
        try:
            return f"<IsomericNucleus(z={self.Z}, n={self.N}, m={self.m}, symbol={self.element_symbol})>"
        except LookupError:
            return f"<IsomericNucleus(z={self.Z}, n={self.N}, m={self.m})>"

    def __repr__(self, ) -> str:
        """Return a string representation of the object which can be used to initialize. """
        return f"IsomericNucleus(z={self.Z}, n={self.N}, m={self.m})"
