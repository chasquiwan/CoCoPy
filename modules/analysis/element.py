#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   ChemPy - A chemistry toolkit for Python
#
#   Copyright (c) 2010 by Joshua W. Allen (jwallen@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

################################################################################
#
# CoCoPy - A python toolkit for rotational spectroscopy
#
# Copyright (c) 2013 by David Schmitz (david.schmitz@chasquiwan.de).
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of 
# this software and associated documentation files (the “Software”), to deal in the 
# Software without restriction, including without limitation the rights to use, 
# copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
# Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN 
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH 
# THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# 
# MIT Licence (http://mit-license.org/)
#
################################################################################


"""
This module contains information about the chemical elements. Information for
each element is stored as attributes of an object of the :class:`Element` 
class. 

Element objects for each chemical element (1-112) have also been declared as 
module-level variables, using each element's symbol as its variable name. These
should be used in most cases to conserve memory.
"""

################################################################################

class Element:
    """
    A chemical element. The attributes are:

    =========== =============== ================================================
    Attribute   Type            Description
    =========== =============== ================================================
    `number`    ``int``         The atomic number of the element
    `symbol`    ``str``         The symbol used for the element
    `name`      ``str``         The IUPAC name of the element
    `mass`      ``float``       The mass of the element in kg/mol
    =========== =============== ================================================
    
    This class is specifically for properties that all atoms of the same element
    share. Ideally there is only one instance of this class for each element.
    """
    
    def __init__(self, number, symbol, name, mass):
        self.number = number
        self.symbol = intern(symbol)
        self.name = name
        self.mass = mass
    
    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return self.symbol
    
    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "Element(%s, '%s', '%s', %s)" % (self.number, self.symbol, self.name, self.mass)
    
################################################################################

def getElement(number=0, symbol=''):
    """
    Return the :class:`Element` object with attributes defined by the given
    parameters. Only the parameters explicitly given will be used, so you can
    search by atomic `number` or by `symbol` independently.
    """
    for element in elementList:
        if (number == 0 or element.number == number) and (symbol == '' or element.symbol == symbol):
            return element
    # If we reach this point that means we did not find an appropriate element,
    # so we raise an exception
    # TODO

################################################################################

# Declare an instance of each element (1 to 112)
# The variable names correspond to each element's symbol
# The elements are sorted by increasing atomic number and grouped by period
# Recommended IUPAC nomenclature is used throughout (including 'aluminium' and 
# 'caesium')

# Period 1
H  = Element(1,   'H' , 'hydrogen'      , 1.0078250)
He = Element(2,   'He', 'helium'        , 4.002602)

# Period 2
Li = Element(3,   'Li', 'lithium'       , 6.941)
Be = Element(4,   'Be', 'beryllium'     , 9.012182)
B  = Element(5,   'B',  'boron'         , 10.811)
C  = Element(6,   'C' , 'carbon'        , 12.0)
N  = Element(7,   'N' , 'nitrogen'      , 14.0030740)
O  = Element(8,   'O' , 'oxygen'        , 15.994914619)
F  = Element(9,   'F' , 'fluorine'      , 18.998403)
Ne = Element(10,  'Ne', 'neon'          , 20.1797)

# Period 3
Na = Element(11,  'Na', 'sodium'        , 22.989770)
Mg = Element(12,  'Mg', 'magnesium'     , 24.3050)
Al = Element(13,  'Al', 'aluminium'     , 26.981538)
Si = Element(14,  'Si', 'silicon'       , 28.0855)
P  = Element(15,  'P' , 'phosphorus'    , 30.973761)
S  = Element(16,  'S' , 'sulfur'        , 32.065)
Cl = Element(17,  'Cl', 'chlorine'      , 35.453)
Ar = Element(18,  'Ar', 'argon'         , 39.348)

# Period 4
K  = Element(19,  'K' , 'potassium'     , 39.0983)
Ca = Element(20,  'Ca', 'calcium'       , 40.078)
Sc = Element(21,  'Sc', 'scandium'      , 44.955910)
Ti = Element(22,  'Ti', 'titanium'      , 47.867)
V  = Element(23,  'V' , 'vanadium'      , 50.9415)
Cr = Element(24,  'Cr', 'chromium'      , 51.9961)
Mn = Element(25,  'Mn', 'manganese'     , 54.938049)
Fe = Element(26,  'Fe', 'iron'          , 55.845)
Co = Element(27,  'Co', 'cobalt'        , 58.933200)
Ni = Element(28,  'Ni', 'nickel'        , 58.6934)
Cu = Element(29,  'Cu', 'copper'        , 63.546)
Zn = Element(30,  'Zn', 'zinc'          , 65.409)
Ga = Element(31,  'Ga', 'gallium'       , 69.723)
Ge = Element(32,  'Ge', 'germanium'     , 73.9211788)
As = Element(33,  'As', 'arsenic'       , 74.92160)
Se = Element(34,  'Se', 'selenium'      , 78.96)  
Br = Element(35,  'Br', 'bromine'       , 78.9183361)
Kr = Element(36,  'Kr', 'krypton'       , 83.798)

# Period 5
Rb = Element(37,  'Rb', 'rubidium'      , 85.4678)
Sr = Element(38,  'Sr', 'strontium'     , 87.62)
Y  = Element(39,  'Y' , 'yttrium'       , 88.90585)
Zr = Element(40,  'Zr', 'zirconium'     , 91.224)
Nb = Element(41,  'Nb', 'niobium'       , 92.90638)
Mo = Element(42,  'Mo', 'molybdenum'    , 95.94)
Tc = Element(43,  'Tc', 'technetium'    , 98.00)
Ru = Element(44,  'Ru', 'ruthenium'     , 101.07)
Rh = Element(45,  'Rh', 'rhodium'       , 102.90550)
Pd = Element(46,  'Pd', 'palladium'     , 106.42)
Ag = Element(47,  'Ag', 'silver'        , 107.8682)
Cd = Element(48,  'Cd', 'cadmium'       , 112.411)
In = Element(49,  'In', 'indium'        , 114.818)
Sn = Element(50,  'Sn', 'tin'           , 118.710)
Sb = Element(51,  'Sb', 'antimony'      , 121.760)
Te = Element(52,  'Te', 'tellurium'     , 127.60)
I  = Element(53,  'I' , 'iodine'        , 126.90447)
Xe = Element(54,  'Xe', 'xenon'         , 131.293)

# Period 6
Cs = Element(55,  'Cs', 'caesium'       , 132.90545)
Ba = Element(56,  'Ba', 'barium'        , 137.327)
La = Element(57,  'La', 'lanthanum'     , 138.9055)
Ce = Element(58,  'Ce', 'cerium'        , 140.116)
Pr = Element(59,  'Pr', 'praesodymium'  , 140.90765)
Nd = Element(60,  'Nd', 'neodymium'     , 144.24)
Pm = Element(61,  'Pm', 'promethium'    , 145.00)
Sm = Element(62,  'Sm', 'samarium'      , 150.36)
Eu = Element(63,  'Eu', 'europium'      , 151.964)
Gd = Element(64,  'Gd', 'gadolinium'    , 157.25)
Tb = Element(65,  'Tb', 'terbium'       , 158.92534)
Dy = Element(66,  'Dy', 'dysprosium'    , 162.500)
Ho = Element(67,  'Ho', 'holmium'       , 164.93032)
Er = Element(68,  'Er', 'erbium'        , 167.259)
Tm = Element(69,  'Tm', 'thulium'       , 168.93421)
Yb = Element(70,  'Yb', 'ytterbium'     , 173.04)
Lu = Element(71,  'Lu', 'lutetium'      , 174.967)
Hf = Element(72,  'Hf', 'hafnium'       , 178.49)
Ta = Element(73,  'Ta', 'tantalum'      , 180.9479)
W  = Element(74,  'W' , 'tungsten'      , 183.84)
Re = Element(75,  'Re', 'rhenium'       , 186.207)
Os = Element(76,  'Os', 'osmium'        , 190.23)
Ir = Element(77,  'Ir', 'iridium'       , 192.217)
Pt = Element(78,  'Pt', 'platinum'      , 195.078)
Au = Element(79,  'Au', 'gold'          , 196.96655)
Hg = Element(80,  'Hg', 'mercury'       , 200.59)
Tl = Element(81,  'Tl', 'thallium'      , 204.3833)
Pb = Element(82,  'Pb', 'lead'          , 207.2)
Bi = Element(83,  'Bi', 'bismuth'       , 208.98038)
Po = Element(84,  'Po', 'polonium'      , 209.0)
At = Element(85,  'At', 'astatine'      , 210.0)
Rn = Element(86,  'Rn', 'radon'         , 222.0)

# Period 7
Fr = Element(87,  'Fr', 'francium'      , 223.00)
Ra = Element(88,  'Ra', 'radium'        , 226.0)
Ac = Element(89,  'Ac', 'actinum'       , 227.0)
Th = Element(90,  'Th', 'thorium'       , 232.0381)
Pa = Element(91,  'Pa', 'protactinum'   , 231.03588)
U  = Element(92,  'U' , 'uranium'       , 238.02891)
Np = Element(93,  'Np', 'neptunium'     , 237.00)
Pu = Element(94,  'Pu', 'plutonium'     , 244.0642)
Am = Element(95,  'Am', 'americium'     , 243.00)
Cm = Element(96,  'Cm', 'curium'        , 247.00)
Bk = Element(97,  'Bk', 'berkelium'     , 247.00)
Cf = Element(98,  'Cf', 'californium'   , 251.00)
Es = Element(99,  'Es', 'einsteinium'   , 252.00)
Fm = Element(100, 'Fm', 'fermium'       , 257.00)
Md = Element(101, 'Md', 'mendelevium'   , 258.00)
No = Element(102, 'No', 'nobelium'      , 259.00)
Lr = Element(103, 'Lr', 'lawrencium'    , 262.00)
Rf = Element(104, 'Rf', 'rutherfordium' , 261.00)
Db = Element(105, 'Db', 'dubnium'       , 262.00)
Sg = Element(106, 'Sg', 'seaborgium'    , 266.00)
Bh = Element(107, 'Bh', 'bohrium'       , 264.00)
Hs = Element(108, 'Hs', 'hassium'       , 277.00)
Mt = Element(109, 'Mt', 'meitnerium'    , 268.00)
Ds = Element(110, 'Ds', 'darmstadtium'  , 281.00)
Rg = Element(111, 'Rg', 'roentgenium'   , 272.00)
Cn = Element(112, 'Cn', 'copernicum'    , 285.00)

# A list of the elements, sorted by increasing atomic number
elementList = [
    H, He,
    Li, Be, B, C, N, O, F, Ne,
    Na, Mg, Al, Si, P, S, Cl, Ar,
    K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
    Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe,
    Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
    Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn
]
