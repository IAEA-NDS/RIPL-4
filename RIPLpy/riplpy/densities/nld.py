
# OS
import os as _os

# Dataclasses
from dataclasses import dataclass

# ========================

# Current location
__loc__ = _os.path.dirname(_os.path.realpath(__file__))

# ========================

@dataclass
class NLD(object):
    """A simple representation of a nuclear level density (NLD). """

    energy: list
    levels: list

