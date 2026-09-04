# -*- coding: utf-8 -*-
"""Python objects that provide quick access to element information.

This module defines mappings between proton numbers (Z), element symbols, and element names. It provides a singleton `Elements` class for convenient access to this information, including methods for retrieving element properties and validating input.

Python 3.10+ is expected to run this code properly.

Examples:
    (1) Print the symbol given a proton number::
       $ print(ZtoSym[25])
    (2) The reverse of example (1)::
       $ print(SymtoZ['Mn'])
    (3) Use the Elements class, which has a bit more input protect and provides all mappings::
       $ print(Elements.symbol(25))
       $ print(Elements.name(25))
       $ print(Elements.element('mn'))
       $ print(Elements.element('manganese'))
    (4) Print out all element names using the Elements class::
       $ print(Elements.NAMES)
    (5) Print out all element symbols using the Elements class::
       $ print(Elements.SYMBOLS)

"""

# ========================

# RegEx
import re

# Local
from .exceptions import ElementNotFoundError
# Map proton number to element symbol
ZtoSym = {0  : 'Nn',
          1  : 'H',
          2  : 'He',
          3  : 'Li',
          4  : 'Be',
          5  : 'B',
          6  : 'C',
          7  : 'N',
          8  : 'O',
          9  : 'F',
          10 : 'Ne',
          11 : 'Na',
          12 : 'Mg',
          13 : 'Al',
          14 : 'Si',
          15 : 'P',
          16 : 'S',
          17 : 'Cl',
          18 : 'Ar',
          19 : 'K',
          20 : 'Ca',
          21 : 'Sc',
          22 : 'Ti',
          23 : 'V',
          24 : 'Cr',
          25 : 'Mn',
          26 : 'Fe',
          27 : 'Co',
          28 : 'Ni',
          29 : 'Cu',
          30 : 'Zn',
          31 : 'Ga',
          32 : 'Ge',
          33 : 'As',
          34 : 'Se',
          35 : 'Br',
          36 : 'Kr',
          37 : 'Rb',
          38 : 'Sr',
          39 : 'Y',
          40 : 'Zr',
          41 : 'Nb',
          42 : 'Mo',
          43 : 'Tc',
          44 : 'Ru',
          45 : 'Rh',
          46 : 'Pd',
          47 : 'Ag',
          48 : 'Cd',
          49 : 'In',
          50 : 'Sn',
          51 : 'Sb',
          52 : 'Te',
          53 : 'I',
          54 : 'Xe',
          55 : 'Cs',
          56 : 'Ba',
          57 : 'La',
          58 : 'Ce',
          59 : 'Pr',
          60 : 'Nd',
          61 : 'Pm',
          62 : 'Sm',
          63 : 'Eu',
          64 : 'Gd',
          65 : 'Tb',
          66 : 'Dy',
          67 : 'Ho',
          68 : 'Er',
          69 : 'Tm',
          70 : 'Yb',
          71 : 'Lu',
          72 : 'Hf',
          73 : 'Ta',
          74 : 'W',
          75 : 'Re',
          76 : 'Os',
          77 : 'Ir',
          78 : 'Pt',
          79 : 'Au',
          80 : 'Hg',
          81 : 'Tl',
          82 : 'Pb',
          83 : 'Bi',
          84 : 'Po',
          85 : 'At',
          86 : 'Rn',
          87 : 'Fr',
          88 : 'Ra',
          89 : 'Ac',
          90 : 'Th',
          91 : 'Pa',
          92 : 'U',
          93 : 'Np',
          94 : 'Pu',
          95 : 'Am',
          96 : 'Cm',
          97 : 'Bk',
          98 : 'Cf',
          99 : 'Es',
          100: 'Fm',
          101: 'Md',
          102: 'No',
          103: 'Lr',
          104: 'Rf',
          105: 'Db',
          106: 'Sg',
          107: 'Bh',
          108: 'Hs',
          109: 'Mt',
          110: 'Ds',
          111: 'Rg',
          112: 'Cn',
          113: 'Nh',
          114: 'Fl',
          115: 'Mc',
          116: 'Lv',
          117: 'Ts',
          118: 'Og',
          119: 'Uue',
          120: 'Ubn',
          121: 'Ubu',
          122: 'Ubb',
          123: 'Ubt',
          124: 'Ubq',
          125: 'Ubp',
          126: 'Ubh',
          127: 'Ubs',
          128: 'Ubo',
          129: 'Ube',
          130: 'Utn',
          131: 'Utu',
          132: 'Utb',
          133: 'Utt',
          134: 'Utq',
          135: 'Utp',
          136: 'Uth',
          137: 'Uts',
          138: 'Uto',
          139: 'Ute',
          140: 'Uqn'}

# Map element symbol to proton number
SymtoZ = dict(map(reversed, ZtoSym.items()))

# Map proton number to element name
ZtoName = {0  :  "Neutron",
           1  :  "Hydrogen",
           2  :  "Helium",
           3  :  "Lithium",
           4  :  "Beryllium",
           5  :  "Boron",
           6  :  "Carbon",
           7  :  "Nitrogen",
           8  :  "Oxygen",
           9  :  "Fluorine",
           10 : "Neon",
           11 : "Sodium",
           12 : "Magnesium",
           13 : "Aluminum",
           14 : "Silicon",
           15 : "Phosphorus",
           16 : "Sulfur",
           17 : "Chlorine",
           18 : "Argon",
           19 : "Potassium",
           20 : "Calcium",
           21 : "Scandium",
           22 : "Titanium",
           23 : "Vanadium",
           24 : "Chromium",
           25 : "Manganese",
           26 : "Iron",
           27 : "Cobalt",
           28 : "Nickel",
           29 : "Copper",
           30 : "Zinc",
           31 : "Gallium",
           32 : "Germanium",
           33 : "Arsenic",
           34 : "Selenium",
           35 : "Bromine",
           36 : "Krypton",
           37 : "Rubidium",
           38 : "Strontium",
           39 : "Yttrium",
           40 : "Zirconium",
           41 : "Niobium",
           42 : "Molybdenum",
           43 : "Technetium",
           44 : "Ruthenium",
           45 : "Rhodium",
           46 : "Palladium",
           47 : "Silver",
           48 : "Cadmium",
           49 : "Indium",
           50 : "Tin",
           51 : "Antimony",
           52 : "Tellurium",
           53 : "Iodine",
           54 : "Xenon",
           55 : "Cesium",
           56 : "Barium",
           57 : "Lanthanum",
           58 : "Cerium",
           59 : "Praseodymium",
           60 : "Neodymium",
           61 : "Promethium",
           62 : "Samarium",
           63 : "Europium",
           64 : "Gadolinium",
           65 : "Terbium",
           66 : "Dysprosium",
           67 : "Holmium",
           68 : "Erbium",
           69 : "Thulium",
           70 : "Ytterbium",
           71 : "Lutetium",
           72 : "Hafnium",
           73 : "Tantalum",
           74 : "Tungsten",
           75 : "Rhenium",
           76 : "Osmium",
           77 : "Iridium",
           78 : "Platinum",
           79 : "Gold",
           80 : "Mercury",
           81 : "Thallium",
           82 : "Lead",
           83 : "Bismuth",
           84 : "Polonium",
           85 : "Astatine",
           86 : "Radon",
           87 : "Francium",
           88 : "Radium",
           89 : "Actinium",
           90 : "Thorium",
           91 : "Protactinium",
           92 : "Uranium",
           93 : "Neptunium",
           94 : "Plutonium",
           95 : "Americium",
           96 : "Curium",
           97 : "Berkelium",
           98 : "Californium",
           99 : "Einsteinium",
           100: "Fermium",
           101: "Mendelevium",
           102: "Nobelium",
           103: "Lawrencium",
           104: "Rutherfordium",
           105: "Dubnium",
           106: "Seaborgium",
           107: "Bohrium",
           108: "Hassium",
           109: "Meitnerium",
           110: "Darmstadtium",
           111: "Roentgenium",
           112: "Copernicium",
           113: "Nihonium",
           114: "Flerovium",
           115: "Moscovium",
           116: "Livermorium",
           117: "Tennessine",
           118: "Oganesson",
           119: "Ununennium",
           120: "Unbinilium",
           121: "Unbiunium",
           122: "Unbibium",
           123: "Unbitrium",
           124: "Unbiquadium",
           125: "Unbipentium",
           126: "Unbihexium",
           127: "Unbiseptium",
           128: "Unbioctium",
           129: "Unbiennium",
           130: "Untrinilium",
           131: "Untriunium",
           132: "Untribium",
           133: "Untritrium",
           134: "Untripentium",
           135: "Untripentium",
           136: "Untrihexium",
           137: "Untriseptium",
           138: "Untrioctium",
           139: "Untriennium",
           140: "Unquadnilium"
           }

# Map element name to proton number
NametoZ = dict(map(reversed, ZtoName.items()))


def isotope_str_to_tuple(istr: str) -> tuple:
    """Convert an isotope string to a tuple of ``(Z, N, A)``.

    Accepts the common notations regardless of separator or of whether the
    element or the mass number comes first, e.g. ``'56Fe'``, ``'Fe-56'``,
    ``'Fe56'``, ``'Fe 56'``, and the full-name forms ``'Iron-56'`` /
    ``'56 iron'``. The element token is matched against element symbols first,
    then full names, both case-insensitively.

    Args:
        istr: The isotope string to parse.

    Returns:
        ``(Z, N, A)``.

    Raises:
        ValueError: If the string does not contain exactly one element token
            and one mass-number token, or the element is not recognized.
    """
    symbols = re.findall(r'[A-Za-z]+', istr)
    numbers = re.findall(r'\d+', istr)
    if len(symbols) != 1 or len(numbers) != 1:
        raise ValueError(
            f"Cannot parse isotope string {istr!r}; expected one element "
            f"symbol/name and one mass number (e.g. '56Fe' or 'Fe-56')."
        )

    token = symbols[0].capitalize()
    A = int(numbers[0])
    Z = SymtoZ.get(token, NametoZ.get(token))
    if Z is None:
        raise ValueError(
            f"Unrecognized element '{symbols[0]}' in isotope string {istr!r}."
        )
    return (Z, A - Z, A)


class ElementsMeta(type):
    """Generic metaclass for Elements object. Serves as a basis for the Elements singleton class.
    """

    ZtoSym : dict = ZtoSym
    SymtoZ : dict = SymtoZ
    ZtoName: dict = ZtoName
    NametoZ: dict = NametoZ

    @property
    def NAMES(cls, ) -> tuple:
        """Return all the element names.

        Returns:
            tuple: A tuple of names (str) of the elements.
        """
        return tuple(cls.NametoZ.keys())
    
    @property
    def SYMBOLS(cls, ) -> tuple:
        """Return all the element symbols.

        Returns:
            tuple: A tuple of the element symbols (str) of the elements.
        """        
        return tuple(cls.SymtoZ.keys())
    
    @property
    def ELEMENTS(cls, ) -> tuple:
        """Return all the element numbers. """
        return tuple(cls.ZtoSym.keys())


class Elements(object, metaclass=ElementsMeta):
    """A class which contains all element information accessible by straightforward methods. Built on the metaclass ElementsMeta.
    """

    @classmethod
    def name(cls, Z: int) -> str:
        """Return the element name given element proton number (Z).

        Raises:
            ElementNotFoundError: If Z is outside the element table.
        """
        try:
            return cls.ZtoName[int(Z)]
        except KeyError:
            from riplpy.exceptions import ElementNotFoundError
            raise ElementNotFoundError(f"No element with Z={Z}")

    @classmethod
    def names(cls, lst: list) -> list:
        """Return a list of element names given a list of proton numbers. """
        return [cls.name(z) for z in lst]

    @classmethod
    def symbol(cls, Z: int) -> str:
        """Return the element symbol given element proton number (Z).

        Raises:
            ElementNotFoundError: If Z is outside the element table.
        """
        try:
            return cls.ZtoSym[int(Z)]
        except KeyError:
            from riplpy.exceptions import ElementNotFoundError
            raise ElementNotFoundError(f"No element with Z={Z}")


    @classmethod
    def symbols(cls, lst: list) -> list:
        """Return a list of element symbols given a collection of proton numbers. """
        return [cls.symbol(z) for z in lst]

    @classmethod
    def element_from_symbol(cls, sym: str) -> int:
        """Return the element proton number (Z) given the symbol. """
        return cls.SymtoZ[str(sym).capitalize()]

    @classmethod
    def elements_from_symbols(cls, syms: tuple) -> list:
        """Return the element proton number (Z) given the symbols (syms). """
        return [cls.element_from_symbol(_sym) for _sym in syms]

    @classmethod
    def element_from_name(cls, name: str) -> int:
        """Return the element proton number (Z) given the element name. """
        return cls.NametoZ[str(name).capitalize()]

    @classmethod
    def elements_from_names(cls, names: tuple) -> list:
        """Return the element proton number (Z) given the element names (names). """
        return [cls.element_from_name(_name) for _name in names]

    @classmethod
    def element(cls, txt: str) -> int:
        """Attempt to infer the element proton number (Z) from symbol or name. """
        try:
            # First attempt to decipher the symbol
            return cls.element_from_symbol(txt)
        except KeyError:
            # Next try to decipher the name
            try:
                return cls.element_from_name(txt)
            except KeyError:
                raise ElementNotFoundError(f"Element {txt} proton number not found!")

    @classmethod
    def elements(cls, lst: list) -> list:
        """Return a list of element proton numbers given a list of names or symbols (or combination of both). """
        return [cls.element(e) for e in lst]

    @classmethod
    def valid_symbol(cls, symbol: str) -> bool:
        """Check if the provided symbol is in the list of known element symbols. """
        symbol = str(symbol).lower().capitalize()
        if symbol in cls.SYMBOLS:
            return True
        else:
            return False

    @classmethod
    def valid_name(cls, name: str) -> bool:
        """Check if the provided name is in the list of known element names. """
        name = str(name).lower().capitalize()
        if name in cls.NAMES:
            return True
        else:
            return False

    @classmethod
    def valid_element(cls, Z: int) -> bool:
        """Check if the provided element name is in the list of known elements. """
        if Z in cls.ELEMENTS:
            return True
        else:
            return False

    @classmethod
    def validate(cls, input: object, type: str) -> bool:
        """Validate the provided input as a ('name', 'symbol' or 'element'). """
        type = str(type).lower()
        match type:
            case 'name':
                return cls.valid_name(input)
            case 'symbol':
                return cls.valid_symbol(input)
            case 'element':
                return cls.valid_element(input)
            case _:
                raise ValueError("Invalid input!")
