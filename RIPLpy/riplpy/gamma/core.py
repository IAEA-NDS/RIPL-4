

# Collections
from collections import namedtuple as _namedtuple

# Dataclass
from dataclasses import dataclass as _dataclass

# Typing
from typing import ClassVar as _ClassVar

# RIPLpy
import riplpy.collections as _c

# ========================

GDR_ParameterTuple = _namedtuple('GDR_ParameterTuple', 'Z, A, El, E1, W1, E2, W2')

# ========================

@_dataclass
class Base_GDR_Parameter(object):
    """A base class for GDR parameters (experiment or theory).

    Contains the fundamental GDR parameters for up to two resonance peaks
    (split GDR in deformed nuclei).

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        E1: First GDR peak energy [MeV]
        W1: First GDR peak width [MeV]
        E2: Second GDR peak energy [MeV]
        W2: Second GDR peak width [MeV]
    """

    n  : _c.Nuclide = None
    E1 : float   = None
    W1 : float   = None
    E2 : float   = None
    W2 : float   = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'E1': 'First GDR peak energy [MeV]',
        'W1': 'First GDR peak width [MeV]',
        'E2': 'Second GDR peak energy [MeV]',
        'W2': 'Second GDR peak width [MeV]',
    }

    @property
    def nucleus(self, ) -> _c.Nuclide:
        """Additional reference to the nucleus. """
        return self.n

    @property
    def symbol(self, ) -> str:
        """Return the element symbol. """
        return self.n.element_symbol

    @property
    def first(self, ) -> tuple:
        """Return the first GDR (energy, width) pair as a Python tuple. """
        return (self.E1, self.W1)

    @property
    def second(self, ) -> tuple:
        """Return the second GDR (energy, width) pair as a Python tuple. """
        return (self.E2, self.W2)


@_dataclass
class Theory_GDR_Parameter(Base_GDR_Parameter):
    """A theoretical Giant Dipole Resonance (GDR) parameter.

    Extends Base_GDR_Parameter with deformation information from
    the Goldhaber-Teller model predictions.

    Attributes:
        eta: Deformation parameter (symmetry axis / perpendicular diameter)
    """

    eta: float = None

    # Field descriptions for help/info display (extends parent)
    _field_info: _ClassVar[dict] = {
        **Base_GDR_Parameter._field_info,
        'eta': 'Deformation parameter (symmetry axis / perpendicular diameter)',
    }

    # @property
    # def as_dict(self, ) -> dict:
    #     """Return the object as a dictionary. """
    #     return {'n': self.n, 'eta': self.eta, 'E1': self.E1, 'W1': self.W1, 'E2': self.E2, 'W2': self.W2}

    # @property
    # def as_ascii_string(self, ) -> str:
    #     """Return values as a string formatted for the ASCII file. """
    #     # @TODO: Finish proper format
    #     return f"{self.n.Z:3d} {self.n.A:3d} {self.n.sym} {self.E1} {self.W1} {self.E2} {self.W2}"

    # @property
    # def as_tuple(self, ) -> tuple:
    #     """Return the object as a tuple. """
    #     return tuple(self.n, self.eta, self.E1, self.W1, self.E2, self.W2)

    # @property
    # def as_list(self, ) -> list:
    #     """Return the object as a list. """
    #     return list(self.n, self.eta, self.E1, self.W1, self.E2, self.W2)

    # def as_list(self, sort: bool = False) -> list:
    #     """Return the dataset as a list of Theory_GDR_Parameter objects. """
    #     # Attempt to sort
    #     if sort:
    #         keys = sorted(self.data.keys())
    #     else:
    #         keys = self.data.keys()
        
    #     # Placeholder return list
    #     r = []
    #     # Loop over data keys and construct list
    #     for k in keys:
    #         eta = self.data[k]['eta']
    #         E1  = self.data[k]['E1']
    #         W1  = self.data[k]['W1']
    #         E2  = self.data[k]['E2']
    #         W2  = self.data[k]['W2']
    #         r.append(Theory_GDR_Parameter(k, eta, E1, W1, E2, W2))
    #     return r


class Theory_GDR_ParameterMarshaller(object):
    """The marshaller for a Theory_GDR_Parameter object. It can perform conversions of a Theory_GDR_Parameter object to other formats. """

    @staticmethod
    def to_dict(obj: Theory_GDR_Parameter) -> dict:
        """Return the object as a dictionary. """
        return {'n': obj.n, 'eta': obj.eta, 'E1': obj.E1, 'W1': obj.W1, 'E2': obj.E2, 'W2': obj.W2}

    @staticmethod
    def to_ascii_string(obj: Theory_GDR_Parameter) -> str:
        """Return values as a string formatted for the standard ASCII file. """
        # @TODO: Finish proper format
        return f"{obj.n.Z:3d} {obj.n.A:3d} {obj.n.sym} {obj.eta} {obj.E1} {obj.W1} {obj.E2} {obj.W2}"

    @staticmethod
    def to_tuple(obj: Theory_GDR_Parameter) -> tuple:
        """Return the object as a tuple. """
        return tuple(obj.n, obj.eta, obj.E1, obj.W1, obj.E2, obj.W2)

    @staticmethod
    def to_list(obj: Theory_GDR_Parameter) -> list:
        """Return the object as a list. """
        return list(obj.n, obj.eta, obj.E1, obj.W1, obj.E2, obj.W2)

    @staticmethod
    def from_dict(d: dict) -> Theory_GDR_Parameter:
        """Return the object initialized from a dictionary. """
        return Theory_GDR_Parameter(**d)

    @staticmethod
    def from_ascii_string(s: str) -> Theory_GDR_Parameter:
        """Return the object initialized from a string"""
        lst = s.split()
        Z   = int(lst[0])
        A   = int(lst[1])
        eta = lst[3]
        E1  = float(lst[4])
        W1  = float(lst[5])
        E2  = float(lst[6]) 
        W2  = float(lst[7])
        n   = _c.Nuclide(Z=Z, A=A)
        return Theory_GDR_Parameter(n=n, E1=E1, W1=W1, E2=E2, W2=W2, eta=eta)



@_dataclass
class Experimental_GDR_Parameter(Base_GDR_Parameter):
    """An experimental GDR parameter from photoabsorption measurements.

    Contains fitted Lorentzian parameters for experimental GDR data
    with uncertainties and reference information.

    Attributes:
        Er1: First component energy [MeV]
        dEr1: Uncertainty on Er1 [MeV]
        Wr1: First component width [MeV]
        dWr1: Uncertainty on Wr1 [MeV]
        Sr1: First component strength [TRK units]
        dSr1: Uncertainty on Sr1
        Er2: Second component energy [MeV]
        dEr2: Uncertainty on Er2 [MeV]
        Wr2: Second component width [MeV]
        dWr2: Uncertainty on Wr2 [MeV]
        Sr2: Second component strength [TRK units]
        dSr2: Uncertainty on Sr2
        Sr: Total strength Sr1+Sr2 [TRK units]
        dSr: Uncertainty on Sr
        Id: Identifier string
        Emin: Lower energy limit of fit [MeV]
        Emax: Upper energy limit of fit [MeV]
        reac: Reaction flag
        reference: Reference key
    """

    Er1      : float = None
    dEr1     : float = None
    Wr1      : float = None
    dWr1     : float = None
    Sr1      : float = None
    dSr1     : float = None
    CSp1     : float = None
    dCSp1    : float = None
    Er2      : float = None
    dEr2     : float = None
    Wr2      : float = None
    dWr2     : float = None
    Sr2      : float = None
    dSr2     : float = None
    CSp2     : float = None
    dCSp2    : float = None
    Sr       : float = None
    dSr      : float = None
    Id       : str   = None
    Emin     : float = None
    Emax     : float = None
    reac     : int   = None
    reference: str   = None

    # Field descriptions for help/info display (extends parent)
    _field_info: _ClassVar[dict] = {
        **Base_GDR_Parameter._field_info,
        'Er1': 'First component energy [MeV]',
        'dEr1': 'Uncertainty on Er1 [MeV]',
        'Wr1': 'First component width [MeV]',
        'dWr1': 'Uncertainty on Wr1 [MeV]',
        'Sr1': 'First component strength [TRK units]',
        'dSr1': 'Uncertainty on Sr1',
        'CSp1': 'First component Lorentzian peak cross section [mb]',
        'dCSp1': 'Uncertainty on CSp1 [mb]',
        'Er2': 'Second component energy [MeV]',
        'dEr2': 'Uncertainty on Er2 [MeV]',
        'Wr2': 'Second component width [MeV]',
        'dWr2': 'Uncertainty on Wr2 [MeV]',
        'Sr2': 'Second component strength [TRK units]',
        'dSr2': 'Uncertainty on Sr2',
        'CSp2': 'Second component Lorentzian peak cross section [mb]',
        'dCSp2': 'Uncertainty on CSp2 [mb]',
        'Sr': 'Total strength Sr1+Sr2 [TRK units]',
        'dSr': 'Uncertainty on Sr',
        'Id': 'Identifier string',
        'Emin': 'Lower energy limit of fit [MeV]',
        'Emax': 'Upper energy limit of fit [MeV]',
        'reac': 'Reaction flag',
        'reference': 'Reference key',
    }


@_dataclass
class ExperimentalEntry(Experimental_GDR_Parameter):
    """A representation of a single experimental Giant Dipole Resonance (GDR) parameter entry in RIPL. """

    reac     : int   = None # Reaction flag (valid numbers are 1-10)
    reference: str   = None # References on the experimental data

    def check(self, ) -> None:
        """Perform sanity checks regarding the class variables that define the object. """
        assert self.reac in range(0, 11)


@_dataclass
class Systematics_GDR_Parameter(Base_GDR_Parameter):
    """Combined experimental/systematics GDR parameter (RIPL-4 SLO+SMLO systematics).

    Used for entries that may come either from experimental fits (``In=1``) or
    from the global systematics fits (``In=0``). The strength quantities follow
    the same conventions as :class:`Experimental_GDR_Parameter`.

    Attributes:
        Er1, Wr1, S1, CSp1: First-mode GDR energy, width, strength (TRK), peak XS [mb]
        Er2, Wr2, S2, CSp2: Second-mode GDR energy, width, strength (TRK), peak XS [mb]
        S: Total strength S1+S2 in units of the TRK sum rule
        In: 1 = experimental, 0 = systematics
    """

    Sr1  : float = None
    CSp1 : float = None
    Sr2  : float = None
    CSp2 : float = None
    Sr   : float = None
    In   : int   = None

    _field_info: _ClassVar[dict] = {
        **Base_GDR_Parameter._field_info,
        'Sr1':  'First component strength [TRK units]',
        'CSp1': 'First component Lorentzian peak cross section [mb]',
        'Sr2':  'Second component strength [TRK units]',
        'CSp2': 'Second component Lorentzian peak cross section [mb]',
        'Sr':   'Total strength Sr1+Sr2 [TRK units]',
        'In':   'Source flag (1 = experimental, 0 = systematics)',
    }


@_dataclass
class GammaStrengthFunction:
    """A tabulated gamma-ray strength function (GSF) for atomic nuclei.

    Contains the energy-dependent E1 photon strength function used
    in statistical model calculations.

    Attributes:
        n: Target nucleus (Nuclide object with Z, A, N properties)
        U: Excitation energy grid [MeV]
        fE1: E1 strength function values [mb/MeV]
    """

    n   : _c.Nuclide = None
    U   : list = None
    fE1 : list = None

    # Field descriptions for help/info display
    _field_info: _ClassVar[dict] = {
        'n': 'Target nucleus',
        'U': 'Excitation energy grid [MeV]',
        'fE1': 'E1 strength function values [mb/MeV]',
    }

