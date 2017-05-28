from indigoX.exception import IndigoSearchError

class _Element(object):
    def __init__(self, name, symbol, group, period, number, mass, atomic_radii,
                 covalent_radii, vdw_radii, chi, hyper=None):
        self.name = name
        self.symbol = symbol
        self.group = group
        self.period = period
        self.number = number
        self.mass = mass
        self.atomic_radii = atomic_radii
        self.covalent_radii = covalent_radii
        self.vdw_radii = vdw_radii
        self.chi = chi
        if group < 13:
            self.valence = group
        elif group == 18 and symbol == 'He':
            self.valence = 2
        else:
            self.valence = group - 10
            
        if group == 1 or symbol == 'He':
            self.octet = 2
        elif group == 2:
            self.octet = 4
        elif group == 13:
            self.octet = 6
        else:
            self.octet = 8
        
        if hyper is None:
            self.hyper = self.octet
        else:
            self.hyper = hyper
    
    def __str__(self):
        return self.symbol
    
    def __eq__(self, c):
        if self.symbol == c or self.number == c or self.name == c:
            return True
        return object.__eq__(self, c)
    
class _PeriodicTable(object):
    
    def __init__(self, elements):
        self.elements_number = dict()
        self.elements_name = dict()
        self.elements_symbol = dict()
        for e in elements:
            self.elements_number[e.number] = e
            self.elements_name[e.name] = e
            self.elements_symbol[e.symbol] = e
            
    def __getattr__(self, name):
        return self[name]
            
    def __getitem__(self, name):
        try:
            int(name)
        except ValueError:
            intable = False
        else:
            intable = True
            
        if not intable and name.title() in self.elements_name:
            return self.elements_name[name.title()]
        elif not intable and name.title() in self.elements_symbol:
            return self.elements_symbol[name.title()]
        elif intable and int(name) in self.elements_number:
            return self.elements_number[int(name)]
        else:
            raise IndigoSearchError('Unknown element type: {}'
                                    .format(name.title()))            
    
_elements = [
        # Name Symbol Group Period AtomicNumber Mass Radius Cvradius VdWradius chi
        _Element("Nullium",       "X",  18, 9,   0,   0.0000,    0.00,  0.0,  0.0,  0.0 ),
        _Element("Actinium",      "Ac",  3, 7,  89, 227.0278,    1.88,  0.0,  0.0,  1.3 ),
        _Element("Aluminum",      "Al", 13, 3,  13,  26.981539,  1.43,  1.25, 2.05, 1.61),
        _Element("Americium",     "Am",  0, 7,  95, 243.0614,    1.73,  0.0,  0.0,  1.3 ),
        _Element("Antimony",      "Sb", 15, 5,  51, 121.76,      1.82,  1.41, 2.2,  2.05),
        _Element("Argon",         "Ar", 18, 3,  18,  39.948,     1.74,  0.0,  1.91, 0.0 ),
        _Element("Arsenic",       "As", 15, 4,  33,  74.92159,   1.25,  1.21, 2.0,  2.18),
        _Element("Astatine",      "At", 17, 6,  85, 209.9871,    0.0,   0.0,  0.0,  1.96),
        _Element("Barium",        "Ba",  2, 6,  56, 137.327,     2.17,  1.98, 0.0,  0.89),
        _Element("Berkelium",     "Bk",  0, 7,  97, 247.0703,    1.70,  0.0,  0.0,  1.3 ),
        _Element("Beryllium",     "Be",  2, 2,   4,   9.012182,  1.13,  0.89, 0.0,  1.57),
        _Element("Bismuth",       "Bi", 15, 6,  83, 208.98037,   1.55,  1.52, 2.4,  2.0 ),
        _Element("Bohrium",       "Bh",  7, 7, 107, 262.12,      0.0,   0.0,  0.0,  0.0 ),
        _Element("Boron",         "B" , 13, 2,   5,  10.811,     0.83,  0.88, 2.08, 2.04),
        _Element("Bromine",       "Br", 17, 4,  35,  79.904,     0.0,   1.14, 1.95, 2.96, 12),
        _Element("Cadmium",       "Cd", 12, 5,  48, 112.411,     1.49,  1.41, 0.0,  1.69),
        _Element("Caesium",       "Cs",  1, 6,  55, 132.90543,   2.654, 2.35, 2.62, 0.79),
        _Element("Calcium",       "Ca",  2, 4,  20,  40.078,     1.97,  1.74, 0.0,  1.0 ),
        _Element("Californium",   "Cf",  0, 7,  98, 251.0796,    1.69,  0.0,  0.0,  1.3 ),
        _Element("Carbon",         "C", 14, 2,   6,  12.011,     0.77,  0.77, 1.85, 2.55),
        _Element("Cerium",        "Ce",  0, 6,  58, 140.115,     1.825, 1.65, 0.0,  1.12),
        _Element("Chlorine",      "Cl", 17, 3,  17,  35.4527,    0.0,   0.99, 1.81, 3.16),
        _Element("Chromium",      "Cr",  6, 4,  24,  51.9961,    1.25,  0.0,  0.0,  1.66),
        _Element("Cobalt",        "Co",  9, 4,  27,  58.9332,    1.25,  1.16, 0.0,  1.88),
        _Element("Copper",        "Cu", 11, 4,  29,  63.546,     1.28,  1.17, 0.0,  1.9 ),
        _Element("Curium",        "Cm",  0, 7,  96, 247.0703,    1.74,  0.0,  0.0,  1.3 ),
        _Element("Dubnium",       "Db",  4, 7, 104, 261.11,      0.0,   0.0,  0.0,  0.0 ),
        _Element("Dysprosium",    "Dy",  0, 6,  66, 162.5,       1.77,  1.59, 0.0,  1.23),
        _Element("Einsteinium",   "Es",  0, 7,  99, 252.083,     2.03,  0.0,  0.0,  1.3 ),
        _Element("Erbium",        "Er",  0, 6,  68, 167.26,      1.76,  1.57, 0.0,  1.25),
        _Element("Europium",      "Eu",  0, 6,  63, 151.965,     2.04,  1.85, 0.0,  1.2 ),
        _Element("Fermium",       "Fm",  0, 7, 100, 257.0951,    0.0,   0.0,  0.0,  1.3 ),
        _Element("Fluorine",      "F" , 17, 2,   9,  18.9984032, 0.709, 0.58, 1.35, 3.98),
        _Element("Francium",      "Fr",  1, 7,  87, 223.0197,    2.7,   0.0,  0.0,  0.7 ),
        _Element("Gadolinium",    "Gd",  0, 6,  64, 157.25,      1.8,   1.61, 0.0,  0.94),
        _Element("Gallium",       "Ga", 13, 4,  31,  69.723,     1.22,  1.25, 0.0,  1.81),
        _Element("Germanium",     "Ge", 14, 4,  32,  72.61,      1.23,  1.22, 0.0,  2.01),
        _Element("Gold",          "Au", 11, 6,  79, 196.96654,   1.44,  1.34, 0.0,  2.0 ),
        _Element("Hafnium",       "Hf",  4, 6,  72, 178.49,      1.56,  1.44, 0.0,  1.5 ),
        _Element("Hahnium",       "Hn",  8, 7, 108,   0.0,       0.0,   0.0,  0.0,  0.0 ),
        _Element("Helium",        "He", 18, 1,   2,   4.002602,  1.28,  0.0,  1.22, 0.0 ),
        _Element("Holmium",       "Ho",  0, 6,  67, 164.93032,   1.77,  1.58, 0.0,  1.24),
        _Element("Hydrogen",      "H" ,  1, 1,   1,   1.00797,   0.78,  0.3,  1.2,  2.2 ),
        _Element("Indium",        "In", 13, 5,  49, 114.818,     1.63,  1.5,  0.0,  1.78),
        _Element("Iodine",        "I" , 17, 5,  53, 126.90447,   0.0,   1.33, 2.15, 2.66),
        _Element("Iridium",       "Ir",  9, 6,  77, 192.217,     1.36,  1.26, 0.0,  2.28),
        _Element("Iron",          "Fe",  8, 4,  26,  55.845,     1.24,  1.16, 0.0,  1.83),
        _Element("Joliotium",     "Jl",  5, 7, 105, 262.114,     0.0,   0.0,  0.0,  0.0 ),
        _Element("Krypton",       "Kr", 18, 4,  36,  83.80,      0.0,   1.89, 1.98, 0.0 ),
        _Element("Lanthanum",     "La",  3, 6,  57, 138.9055,    1.88,  1.69, 0.0,  1.1 ),
        _Element("Lawrencium",    "Lr",  3, 7, 103, 262.11,      0.0,   0.0,  0.0,  0.0 ),
        _Element("Lead",          "Pb", 14, 6,  82, 207.2,       1.75,  1.54, 0.0,  2.02),
        _Element("Lithium",       "Li",  1, 2,   3,   6.941,     1.52,  1.23, 0.0,  0.98),
        _Element("Lutetium",      "Lu",  3, 6,  71, 174.967,     1.72,  1.56, 0.0,  1.3 ),
        _Element("Magnesium",     "Mg",  2, 3,  12,  24.30506,   1.6,   1.36, 0.0,  1.31),
        _Element("Manganese",     "Mn",  7, 4,  25,  54.93805,   1.24,  1.77, 0.0,  1.55),
        _Element("Meitnerium",    "Mt",  9, 7, 109,   0.0,       0.0,   0.0,  0.0,  0.0 ),
        _Element("Mendelevium",   "Md",  0, 7, 101, 258.1,       0.0,   0.0,  0.0,  1.3 ),
        _Element("Mercury",       "Hg", 12, 6,  80, 200.59,      1.60,  1.44, 0.0,  1.8 ),
        _Element("Molybdenum",    "Mo",  6, 5,  42,  95.94,      1.36,  1.29, 0.0,  2.16),
        _Element("Neodymium",     "Nd",  0, 6,  60, 144.24,      1.82,  1.64, 0.0,  1.14),
        _Element("Neon",          "Ne", 18, 2,  10,  20.1797,    0.0,   0.0,  1.6,  0.0 ),
        _Element("Neptunium",     "Np",  0, 7,  93, 237.0482,    1.5,   0.0,  0.0,  1.28),
        _Element("Nickel",        "Ni", 10, 4,  28,  58.6934,    1.25,  1.15, 0.0,  1.91),
        _Element("Niobium",       "Nb",  5, 5,  41,  92.90638,   1.43,  1.34, 0.0,  1.6 ),
        _Element("Nitrogen",      "N" , 15, 2,   7,  14.00674,   0.71,  0.7,  1.54, 3.04),
        _Element("Nobelium",      "No",  0, 7, 102, 259.1009,    0.0,   0.0,  0.0,  0.0 ),
        _Element("Osmium",        "Os",  8, 6,  76, 190.23,      1.35,  1.26, 0.0,  2.2 ),
        _Element("Oxygen",        "O" , 16, 2,   8,  15.9994,    0.6,   0.66, 1.4,  3.44),
        _Element("Palladium",     "Pd", 10, 5,  46, 106.42,      1.38,  1.28, 0.0,  2.2 ),
        _Element("Phosphorus",    "P" , 15, 3,  15,  30.973762,  1.15,  1.10, 1.9,  2.19, 10),
        _Element("Platinum",      "Pt", 10, 6,  78, 195.08,      1.38,  1.29, 0.0,  2.54),
        _Element("Plutonium",     "Pu",  7, 0,  94, 244.0642,    0.0,   0.0,  0.0,  1.3 ),
        _Element("Polonium",      "Po", 16, 6,  84, 208.9824,    1.67,  1.53, 0.0,  2.2 ),
        _Element("Potassium",     "K" ,  1, 4,  19,  39.0983,    2.27,  2.03, 2.31, 0.82),
        _Element("Praseodymium",  "Pr",  0, 6,  59, 140.90765,   1.83,  1.65, 0.0,  1.13),
        _Element("Promethium",    "Pm",  0, 6,  61, 144.9127,    1.81,  0.0,  0.0,  0.94),
        _Element("Protactinium",  "Pa",  0, 7,  91, 231.03588,   1.61,  0.0,  0.0,  1.38),
        _Element("Radium",        "Ra",  2, 7,  88, 226.0254,    2.23,  0.0,  0.0,  0.89),
        _Element("Radon",         "Rn", 18, 6,  86, 222.0176,    0.0,   0.0,  0.0,  0.7 ),
        _Element("Rhenium",       "Re",  7, 6,  75, 186.207,     1.37,  1.28, 0.0,  2.2 ),
        _Element("Rhodium",       "Rh",  9, 5,  45, 102.9055,    1.34,  1.25, 0.0,  2.28),
        _Element("Rubidium",      "Rb",  1, 5,  37,  85.4678,    1.475, 0.0,  2.44, 0.82),
        _Element("Ruthenium",     "Ru",  8, 5,  44, 101.07,      1.34,  1.24, 0.0,  2.2 ),
        _Element("Rutherfordium", "Rf",  6, 7, 106, 263.118,     0.0,   0.0,  0.0,  0.0 ),
        _Element("Samarium",      "Sm",  0, 6,  62, 150.36,      1.8,   1.66, 0.0,  1.17),
        _Element("Scandium",      "Sc",  3, 4,  21,  44.95591,   1.61,  1.44, 0.0,  1.36),
        _Element("Selenium",      "Se", 16, 4,  34,  78.96,      2.15,  1.17, 2.0,  2.55),
        _Element("Silicon",       "Si", 14, 3,  14,  28.0855,    1.17,  1.17, 2.0,  1.9 ),
        _Element("Silver",        "Ag", 11, 5,  47, 107.8682,    1.44,  1.34, 0.0,  1.93),
        _Element("Sodium",        "Na",  1, 3,  11,  22.989768,  1.54,  0.0,  2.31, 0.93),
        _Element("Strontium",     "Sr",  2, 5,  38,  87.62,      2.15,  1.92, 0.0,  0.95),
        _Element("Sulfur",        "S" , 16, 3,  16,  32.066,     1.04,  1.04, 1.85, 2.58, 12),
        _Element("Tantalum",      "Ta",  5, 6,  73, 180.9479,    1.43,  1.34, 0.0,  2.36),
        _Element("Technetium",    "Tc",  7, 5,  43,  98.9072,    1.36,  0.0,  0.0,  1.9 ),
        _Element("Tellurium",     "Te", 16, 5,  52, 127.6,       1.43,  1.37, 2.2,  2.1 ),
        _Element("Terbium",       "Tb",  0, 6,  65, 158.92534,   1.78,  1.59, 0.0,  1.22),
        _Element("Thallium",      "Tl", 13, 6,  81, 204.3833,    1.7,   1.55, 0.0,  2.33),
        _Element("Thorium",       "Th",  0, 7,  90, 232.0381,    1.80,  0.0,  0.0,  0.0 ),
        _Element("Thulium",       "Tm",  0, 6,  69, 168.93421,   1.75,  1.56, 0.0,  0.96),
        _Element("Tin",           "Sn", 14, 5,  50, 118.71,      1.41,  1.4,  2.0,  1.96),
        _Element("Titanium",      "Ti",  4, 4,  22,  47.867,     1.45,  1.32, 0.0,  1.54),
        _Element("Tungsten",      "W" ,  6, 6,  74, 183.84,      1.37,  1.3,  0.0,  1.9 ),
        _Element("Uranium",       "U" ,  0, 7,  92, 238.0289,    1.54,  0.0,  0.0,  1.26),
        _Element("Vanadium",      "V" ,  5, 4,  23,  50.9415,    1.32,  0.0,  0.0,  1.63),
        _Element("Xenon",         "Xe", 18, 5,  54, 131.29,      2.18,  2.09, 2.16, 2.6 ),
        _Element("Ytterbium",     "Yb",  0, 6,  70, 173.04,      1.94,  1.7,  0.0,  1.27),
        _Element("Yttrium",       "Y" ,  3, 5,  39,  88.90585,   1.81,  1.62, 0.0,  1.22),
        _Element("Zinc",          "Zn", 12, 4,  30,  65.39,      1.33,  1.25, 0.0,  1.65),
        _Element("Zirconium",     "Zr",  4, 5,  40,  91.224,     1.6,   1.45, 0.0,  1.3 )
        ]

PeriodicTable = _PeriodicTable(_elements)
