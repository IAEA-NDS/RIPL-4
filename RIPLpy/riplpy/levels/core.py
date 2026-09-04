

# ========================

# Dataclass
from dataclasses import dataclass as _dataclass
from dataclasses import field as _field

# ========================

__all__ = ('decay_modes', 'decay_modes_shorthand', 'decay_modes_description', 'DecayMode', 'DiscreteNuclearLevel', 'DiscreteLevelInfo', 'DiscreteGamma')

# ========================

decay_modes = {'%B-'    : 'β- decay', 
               '%EC'    : 'electron capture',
               '%EC+%B+': 'electron capture and β+ decay',
               '%N'     : 'neutron decay',
               '%A'     : 'α decay',
               '%IT'    : 'isomeric transition',
               '%P'     : 'proton decay',
               '%3HE'   : '3He decay',
               '%B+P'   : 'β+ delayed proton decay',
               '%B-N'   : 'β- delayed neutron decay',
               '%SF'    : 'spontaneous fission',
               '%ECP'   : 'electron capture delayed proton decay',
               '%ECA'   : 'electron capture delayed alpha decay',
               '%G'     : 'gamma decay',
               '%B-2N'  : 'β- delayed double neutron decay',
               '%B+2P'  : 'β+ delayed double proton decay',
              }

decay_modes_shorthand   = tuple(decay_modes.keys())
decay_modes_description = tuple(decay_modes.values())

# ========================

@_dataclass
class DecayMode(object):
    """A representation of a decay mode for the RIPL database. """

    mode      : str   = None # Decay flag as a string
    flag      : str   = None # Decay percentage modifier flag; informs a user about major uncertainties. 
                             # The modifiers are copied out of ENSDF with no modification, and can 
                             # have the following values: =, <, >, ? (unknown, but expected), 
                             # AP (approximate), GE (greater or equal), LE (less or equal), 
                             # LT (less then), SY  (value from systematics)
    percent   : str   = None # Percentage decay of different decay modes


@_dataclass
class DiscreteNuclearLevel(object):
    """A basic representation of a discrete nuclear level. """

    energy     : float = None # Energy, typically in [MeV]
    spin       : float = None # Spin as a float (with sign)
    parity     : int   = None # Either +1 or -1


@_dataclass
class DiscreteLevelInfo(DiscreteNuclearLevel):
    """An object representing an entry in the discrete nuclear level database. """

    halflife   : float = None # Half-life in [seconds]
    num_gammas : int   = None # Number of gamma-rays de-exciting the level
    jflag      : str   = None # Flag for spin estimation method
    unc        : str   = None # Flag for an uncertain level energy
    jpistr     : str   = None # String representation of spin-parity (original spins from the ENSDF file)
    decay_modes: list  = _field(default_factory=list) # Decay modes of the level
    shift      : str   = None # Value assigned to uncertain "X" in an ENSDF level
    band       : int   = None # Integer assigned to band(s) to which the level is part

    @property
    def number_decay_modes(self) -> int:
        """Return the number of decay modes. Returns None for unknown. """
        return len(self.decay_modes)


@_dataclass
class DiscreteGamma(object):
    """An object representing a gamma-ray from a discrete nuclear level. """

    to_level_number: int   = None # Final level number (unique integer in associated DiscreteLevel list)
    energy        : float = None # Energy, typically in [MeV]
    prob_g        : float = None # Probability that the level decays through photon (gamma ray) emission
    prob_em       : float = None # Probability of the electromagnetic transition (photon, conversion electron, pair creation)
                                 # The sum of the Pe gives the IT (electromagnetic transition) branching ratio of the level
                                 # This is 1 unless other decay modes are listed in the level record (see example below)
    cc            : float = None # Internal conversion coefficient of the transition
                                 # Pe, Pg, and ICC are in the following relation:
                                 # Pe = Pg * (1 + ICC)
