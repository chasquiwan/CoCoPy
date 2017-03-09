#!/usr/bin/python
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
# Copyright (c) 2016 by David Schmitz (david.schmitz@chasquiwan.de).
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

'''
Todo:
    - 1 Write comments
    - 2
    - 3
'''

"""
This module provides classes and methods for working with molecules and
molecular configurations. A molecule is represented internally using a graph
data type, where atoms correspond to vertices and bonds correspond to edges.
Both :class:`Atom` and :class:`Bond` objects store semantic information that
describe the corresponding atom or bond.
"""

import numpy as np
import element as elements
from graph import Vertex, Edge, Graph
import scipy.constants as constants

#from pattern import AtomPattern, BondPattern, MoleculePattern, AtomType
#from pattern import getAtomType, fromAdjacencyList, toAdjacencyList

################################################################################

class Atom(Vertex):
    """
    An atom. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `element`           :class:`Element`    The chemical element the atom represents
    `radicalElectrons`  ``short``           The number of radical electrons
    `spinMultiplicity`  ``short``           The spin multiplicity of the atom
    `implicitHydrogens` ``short``           The number of implicit hydrogen atoms bonded to this atom
    `charge`            ``short``           The formal charge of the atom
    `label`             ``str``             A string label that can be used to tag individual atoms
    =================== =================== ====================================

    Additionally, the ``mass``, ``number``, and ``symbol`` attributes of the
    atom's element can be read (but not written) directly from the atom object,
    e.g. ``atom.symbol`` instead of ``atom.element.symbol``.
    """

    def __init__(self, element=None, coord=np.array([0.0, 0.0, 0.0]), radicalElectrons=0, spinMultiplicity=1, implicitHydrogens=0, charge=0, label=''):
        Vertex.__init__(self)
        if isinstance(element, str):
            self.element = elements.__dict__[element]
        else:
            self.element = element
        self.coord = coord
        self.radicalElectrons = radicalElectrons
        self.spinMultiplicity = spinMultiplicity
        self.implicitHydrogens = implicitHydrogens
        self.charge = charge
        self.label = label
        self.atomType = None
        self.realMass = self.element.mass
        self.mullikenCharge = charge

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<Atom '%s'>" % (
            str(self.element) +
            ''.join(['.' for i in range(self.radicalElectrons)]) +
            ''.join(['+' for i in range(self.charge)]) +
            ''.join(['-' for i in range(-self.charge)])
        )

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "Atom(element='%s', radicalElectrons=%s, spinMultiplicity=%s, implicitHydrogens=%s, charge=%s, label='%s')" % (self.element, self.radicalElectrons, self.spinMultiplicity, self.implicitHydrogens, self.charge, self.label)

    @property
    def mass(self): return self.element.mass

    @property
    def number(self): return self.element.number

    @property
    def symbol(self): return self.element.symbol

    def equivalent(self, other):
        """
        Return ``True`` if `other` is indistinguishable from this atom, or
        ``False`` otherwise. If `other` is an :class:`Atom` object, then all
        attributes except `label` must match exactly. If `other` is an
        :class:`AtomPattern` object, then the atom must match any of the
        combinations in the atom pattern.
        """
        if isinstance(other, Atom):
            atom = other
            return (self.element is atom.element and
                self.radicalElectrons == atom.radicalElectrons and
                self.spinMultiplicity == atom.spinMultiplicity and
                self.implicitHydrogens == atom.implicitHydrogens and
                self.charge == atom.charge)

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise. If `other` is an :class:`Atom` object, then this is the same
        as the :meth:`equivalent()` method. If `other` is an
        :class:`AtomPattern` object, then the atom must match or be more
        specific than any of the combinations in the atom pattern.
        """
        if isinstance(other, Atom):
            return self.equivalent(other)

    def copy(self):
        """
        Generate a deep copy of the current atom. Modifying the
        attributes of the copy will not affect the original.
        """
        a = Atom(self.element, self.radicalElectrons, self.spinMultiplicity, self.implicitHydrogens, self.charge, self.label)
        a.atomType = self.atomType
        return a

    def isHydrogen(self):
        """
        Return ``True`` if the atom represents a hydrogen atom or ``False`` if
        not.
        """
        return self.element.number == 1

    def isNonHydrogen(self):
        """
        Return ``True`` if the atom does not represent a hydrogen atom or
        ``False`` if not.
        """
        return self.element.number > 1

    def isCarbon(self):
        """
        Return ``True`` if the atom represents a carbon atom or ``False`` if
        not.
        """
        return self.element.number == 6

    def isOxygen(self):
        """
        Return ``True`` if the atom represents an oxygen atom or ``False`` if
        not.
        """
        return self.element.number == 8

    def incrementRadical(self):
        """
        Update the atom pattern as a result of applying a GAIN_RADICAL action,
        where `radical` specifies the number of radical electrons to add.
        """
        # Set the new radical electron counts and spin multiplicities
        self.radicalElectrons += 1
        self.spinMultiplicity += 1

    def decrementRadical(self):
        """
        Update the atom pattern as a result of applying a LOSE_RADICAL action,
        where `radical` specifies the number of radical electrons to remove.
        """
        # Set the new radical electron counts and spin multiplicities
        if self.radicalElectrons - 1 < 0:
            print 'Error'
        self.radicalElectrons -= 1
        if self.spinMultiplicity - 1 < 0:
            self.spinMultiplicity -= 1 - 2
        else:
            self.spinMultiplicity -= 1

    def applyAction(self, action):
        """
        Update the atom pattern as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        # Invalidate current atom type
        self.atomType = None
        # Modify attributes if necessary
        if action[0].upper() in ['CHANGE_BOND', 'FORM_BOND', 'BREAK_BOND']:
            # Nothing else to do here
            pass
        elif action[0].upper() == 'GAIN_RADICAL':
            for i in range(action[2]): self.incrementRadical()
        elif action[0].upper() == 'LOSE_RADICAL':
            for i in range(abs(action[2])): self.decrementRadical()
        else:
            print 'error'

################################################################################

class Bond(Edge):
    """
    A chemical bond. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `order`             ``str``             The bond order (``S`` = single, `D`` = double, ``T`` = triple, ``B`` = benzene)
    =================== =================== ====================================

    """

    def __init__(self, order=1):
        Edge.__init__(self)
        self.order = order

    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<Bond '%s'>" % (self.order)

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "Bond(order='%s')" % (self.order)

    def equivalent(self, other):
        """
        Return ``True`` if `other` is indistinguishable from this bond, or
        ``False`` otherwise. `other` can be either a :class:`Bond` or a
        :class:`BondPattern` object.
        """
        if isinstance(other, Bond):
            bond = other
            return (self.order == bond.order)

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is a specific case of `other`, or ``False``
        otherwise. `other` can be either a :class:`Bond` or a
        :class:`BondPattern` object.
        """
        # There are no generic bond types, so isSpecificCaseOf is the same as equivalent
        return self.equivalent(other)

    def copy(self):
        """
        Generate a deep copy of the current bond. Modifying the
        attributes of the copy will not affect the original.
        """
        return Bond(self.order)

    def isSingle(self):
        """
        Return ``True`` if the bond represents a single bond or ``False`` if
        not.
        """
        return self.order == 'S'

    def isDouble(self):
        """
        Return ``True`` if the bond represents a double bond or ``False`` if
        not.
        """
        return self.order == 'D'

    def isTriple(self):
        """
        Return ``True`` if the bond represents a triple bond or ``False`` if
        not.
        """
        return self.order == 'T'

    def isBenzene(self):
        """
        Return ``True`` if the bond represents a benzene bond or ``False`` if
        not.
        """
        return self.order == 'B'

    def incrementOrder(self):
        """
        Update the bond as a result of applying a CHANGE_BOND action to
        increase the order by one.
        """
        if self.order == 'S': self.order = 'D'
        elif self.order == 'D': self.order = 'T'
        else:
            print 'Error'

    def decrementOrder(self):
        """
        Update the bond as a result of applying a CHANGE_BOND action to
        decrease the order by one.
        """
        if self.order == 'D': self.order = 'S'
        elif self.order == 'T': self.order = 'D'
        else:
            print 'Error'

    def __changeBond(self, order):
        """
        Update the bond as a result of applying a CHANGE_BOND action,
        where `order` specifies whether the bond is incremented or decremented
        in bond order, and should be 1 or -1.
        """
        if order == 1:
            if self.order == 'S': self.order = 'D'
            elif self.order == 'D': self.order = 'T'
            else:
                print 'Error'
        elif order == -1:
            if self.order == 'D': self.order = 'S'
            elif self.order == 'T': self.order = 'D'
            else:
                print 'Error'
        else:
            print 'Error'

    def applyAction(self, action):
        """
        Update the bond as a result of applying `action`, a tuple
        containing the name of the reaction recipe action along with any
        required parameters. The available actions can be found
        :ref:`here <reaction-recipe-actions>`.
        """
        if action[0].upper() == 'CHANGE_BOND':
            if action[2] == 1:
                self.incrementOrder()
            elif action[2] == -1:
                self.decrementOrder()
            else:
                print 'Error'
        else:
            print 'Error'

################################################################################

class Molecule(Graph):
    """
    A representation of a molecular structure using a graph data type, extending
    the :class:`Graph` class. The `atoms` and `bonds` attributes are aliases
    for the `vertices` and `edges` attributes. Corresponding alias methods have
    also been provided.
    """

    def __init__(self, atoms=None, bonds=None, SMILES='', InChI='', implicitH=False):
        Graph.__init__(self, atoms, bonds)
        self.implicitHydrogens = False
        self.properties = dict()
        if SMILES != '': self.fromSMILES(SMILES, implicitH)
        elif InChI != '': self.fromInChI(InChI, implicitH)


    def __str__(self):
        """
        Return a human-readable string representation of the object.
        """
        return "<Molecule '%s'>" % (self.toSMILES())

    def __repr__(self):
        """
        Return a representation that can be used to reconstruct the object.
        """
        return "Molecule(SMILES='%s')" % (self.toSMILES())

    def __getAtoms(self): return self.vertices
    def __setAtoms(self, atoms): self.vertices = atoms
    atoms = property(__getAtoms, __setAtoms)

    def __getBonds(self): return self.edges
    def __setBonds(self, bonds): self.edges = bonds
    bonds = property(__getBonds, __setBonds)

    def addAtom(self, atom):
        """
        Add an `atom` to the graph. The atom is initialized with no bonds.
        """
        return self.addVertex(atom)

    def addBond(self, atom1, atom2, bond):
        """
        Add a `bond` to the graph as an edge connecting the two atoms `atom1`
        and `atom2`.
        """
        return self.addEdge(atom1, atom2, bond)

    def getBonds(self, atom):
        """
        Return a list of the bonds involving the specified `atom`.
        """
        return self.getEdges(atom)

    def getBond(self, atom1, atom2):
        """
        Returns the bond connecting atoms `atom1` and `atom2`.
        """
        return self.getEdge(atom1, atom2)

    def hasAtom(self, atom):
        """
        Returns ``True`` if `atom` is an atom in the graph, or ``False`` if
        not.
        """
        return self.hasVertex(atom)

    def hasBond(self, atom1, atom2):
        """
        Returns ``True`` if atoms `atom1` and `atom2` are connected
        by an bond, or ``False`` if not.
        """
        return self.hasEdge(atom1, atom2)

    def removeAtom(self, atom):
        """
        Remove `atom` and all bonds associated with it from the graph. Does
        not remove atoms that no longer have any bonds as a result of this
        removal.
        """
        return self.removeVertex(atom)

    def removeBond(self, atom1, atom2):
        """
        Remove the bond between atoms `atom1` and `atom2` from the graph.
        Does not remove atoms that no longer have any bonds as a result of
        this removal.
        """
        return self.removeEdge(atom1, atom2)

    def sortAtoms(self):
        """
        Sort the atoms in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """
        return self.sortVertices()

    def getFormula(self):
        """
        Return the molecular formula for the molecule.
        """
        import pybel
        mol = pybel.Molecule(self.toOBMol())
        return mol.formula

    def getMolecularWeight(self):
        """
        Return the molecular weight of the molecule in kg/mol.
        """
        return sum([atom.realMass for atom in self.vertices])

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        g = Graph.copy(self, deep)
        other = Molecule(g.vertices, g.edges)
        return other

    def merge(self, other):
        """
        Merge two molecules so as to store them in a single :class:`Molecule`
        object. The merged :class:`Molecule` object is returned.
        """
        g = Graph.merge(self, other)
        molecule = Molecule(atoms=g.vertices, bonds=g.edges)
        return molecule

    def split(self):
        """
        Convert a single :class:`Molecule` object containing two or more
        unconnected molecules into separate class:`Molecule` objects.
        """
        graphs = Graph.split(self)
        molecules = []
        for g in graphs:
            molecule = Molecule(atoms=g.vertices, bonds=g.edges)
            molecules.append(molecule)
        return molecules

    def makeHydrogensImplicit(self):
        """
        Convert all explicitly stored hydrogen atoms to be stored implicitly.
        An implicit hydrogen atom is stored on the heavy atom it is connected
        to as a single integer counter. This is done to save memory.
        """
        # Check that the structure contains at least one heavy atom
        for atom in self.vertices:
            if not atom.isHydrogen():
                break
        else:
            # No heavy atoms, so leave explicit
            return

        # Count the hydrogen atoms on each non-hydrogen atom and set the
        # `implicitHydrogens` attribute accordingly
        hydrogens = []
        for atom in self.vertices:
            if atom.isHydrogen():
                neighbor = self.edges[atom].keys()[0]
                neighbor.implicitHydrogens += 1
                hydrogens.append(atom)

        # Remove the hydrogen atoms from the structure
        for atom in hydrogens:
            self.removeAtom(atom)

        # Set implicitHydrogens flag to True
        self.implicitHydrogens = True

    def makeHydrogensExplicit(self):
        """
        Convert all implicitly stored hydrogen atoms to be stored explicitly.
        An explicit hydrogen atom is stored as its own atom in the graph, with
        a single bond to the heavy atom it is attached to. This consumes more
        memory, but may be required for certain tasks (e.g. subgraph matching).
        """


        # Create new hydrogen atoms for each implicit hydrogen
        hydrogens = []
        for atom in self.vertices:
            while atom.implicitHydrogens > 0:
                H = Atom(element='H')
                bond = Bond(order='S')
                hydrogens.append((H, atom, bond))
                atom.implicitHydrogens -= 1

        # Add the hydrogens to the graph
        numAtoms = len(self.vertices)
        for H, atom, bond in hydrogens:
            self.addAtom(H)
            self.addBond(H, atom, bond)
            # If known, set the connectivity information
            H.connectivity1 = 1
            H.connectivity2 = atom.connectivity1
            H.connectivity3 = atom.connectivity2
            H.sortingLabel = numAtoms
            numAtoms += 1

        # Set implicitHydrogens flag to False
        self.implicitHydrogens = False

    def clearLabeledAtoms(self):
        """
        Remove the labels from all atoms in the molecule.
        """
        for atom in self.vertices:
            atom.label = ''

    def containsLabeledAtom(self, label):
        """
        Return :data:`True` if the molecule contains an atom with the label
        `label` and :data:`False` otherwise.
        """
        for atom in self.vertices:
            if atom.label == label: return True
        return False

    def getLabeledAtom(self, label):
        """
        Return the atoms in the molecule that are labeled.
        """
        for atom in self.vertices:
            if atom.label == label: return atom
        return None

    def getLabeledAtoms(self):
        """
        Return the labeled atoms as a ``dict`` with the keys being the labels
        and the values the atoms themselves. If two or more atoms have the
        same label, the value is converted to a list of these atoms.
        """
        labeled = {}
        for atom in self.vertices:
            if atom.label != '':
                if atom.label in labeled:
                    labeled[atom.label] = [labeled[atom.label]]
                    labeled[atom.label].append(atom)
                else:
                    labeled[atom.label] = atom
        return labeled

    def fromCML(self, cmlstr, implicitH=False):
        """
        Convert a string of CML `cmlstr` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import pybel
        cmlstr = cmlstr.replace('\t', '')
        mol = pybel.readstring('cml', cmlstr)
        self.fromOBMol(mol.OBMol, implicitH)
        return self

    def fromInChI(self, inchistr, implicitH=False):
        """
        Convert an InChI string `inchistr` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import pybel
        mol = pybel.readstring('inchi', inchistr)
        self.fromOBMol(mol.OBMol, implicitH)
        return self

    def fromSMILES(self, smilesstr, implicitH=False):
        """
        Convert a SMILES string `smilesstr` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import pybel
        mol = pybel.readstring('smiles', smilesstr)
        self.fromOBMol(mol.OBMol, implicitH)
        return self

    def fromFile(self, fname='', ftype='g03', calctype='opt', implicitH=False):
        """
        Convert a SMILES string `smilesstr` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        if fname != '':
            import pybel
            mol = pybel.readfile(ftype, fname).next()
            self.fromOBMol(mol.OBMol, implicitH)
            if ftype == 'g03' or ftype == 'g09':
                import gaussian
                self.properties.update(gaussian.readout(fname, calctype))
            return self

    def fromString(self, strg='', ftype='g03', implicitH=False):
        """
        Convert a SMILES string `smilesstr` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        if strg != '':
            import pybel
            mol = pybel.readstring(ftype, strg)
            self.fromOBMol(mol.OBMol, implicitH)

            return self

    def fromOBMol(self, obmol, implicitH=False):
        """
        Convert an OpenBabel OBMol object `obmol` to a molecular structure. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """

        self.vertices = []
        self.edges = {}

        # Add hydrogen atoms to complete molecule if needed
        obmol.AddHydrogens()

        # Iterate through atoms in obmol
        for i in range(0, obmol.NumAtoms()):
            obatom = obmol.GetAtom(i + 1)

            # Use atomic number as key for element
            number = obatom.GetAtomicNum()
            element = elements.getElement(number=number)


            # Process spin multiplicity
            radicalElectrons = 0
            spinMultiplicity = obatom.GetSpinMultiplicity()
            if spinMultiplicity == 0:
                radicalElectrons = 0; spinMultiplicity = 1
            elif spinMultiplicity == 1:
                radicalElectrons = 2; spinMultiplicity = 1
            elif spinMultiplicity == 2:
                radicalElectrons = 1; spinMultiplicity = 2
            elif spinMultiplicity == 3:
                radicalElectrons = 2; spinMultiplicity = 3

            # Process charge
            charge = obatom.GetFormalCharge()
            coord = np.array([obatom.x(), obatom.y(), obatom.z()])
            label = obatom.GetIdx()

            atom = Atom(element, coord, radicalElectrons, spinMultiplicity, 0, charge, str(label))
            self.vertices.append(atom)
            self.edges[atom] = {}


            # Add bonds by iterating again through atoms
            for j in range(0, i):
                obatom2 = obmol.GetAtom(j + 1)
                obbond = obatom.GetBond(obatom2)
                if obbond is not None:
                    order = 0

                    # Process bond type
                    if obbond.IsSingle(): order = 'S'
                    elif obbond.IsDouble(): order = 'D'
                    elif obbond.IsTriple(): order = 'T'
                    elif obbond.IsAromatic(): order = 'B'

                    bond = Bond(order)
                    atom1 = self.vertices[i]
                    atom2 = self.vertices[j]
                    self.edges[atom1][atom2] = bond
                    self.edges[atom2][atom1] = bond

        # Make hydrogens implicit to conserve memory
        if implicitH: self.makeHydrogensImplicit()

        return self

    def toCML(self):
        """
        Convert the molecular structure to CML. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import pybel
        mol = pybel.Molecule(self.toOBMol())
        cml = mol.write('cml').strip()
        return '\n'.join([l for l in cml.split('\n') if l.strip()])

    def toInChI(self):
        """
        Convert a molecular structure to an InChI string. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import openbabel
        # This version does not write a warning to stderr if stereochemistry is undefined
        obmol = self.toOBMol()
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat('inchi')
        obConversion.SetOptions('w', openbabel.OBConversion.OUTOPTIONS)
        return obConversion.WriteString(obmol).strip()

    def toSMILES(self):
        """
        Convert a molecular structure to an SMILES string. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import pybel
        mol = pybel.Molecule(self.toOBMol())
        return mol.write('smiles').strip()

    def toFile(self, ftype):
        """
        Convert a molecular structure to an SMILES string. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """
        import pybel
        mol = pybel.Molecule(self.toOBMol())
        return mol.write(ftype)

    def toOBMol(self):
        """
        Convert a molecular structure to an OpenBabel OBMol object. Uses
        `OpenBabel <http://openbabel.org/>`_ to perform the conversion.
        """

        import openbabel

        # Make hydrogens explicit while we perform the conversion
        implicitH = self.implicitHydrogens
        self.makeHydrogensExplicit()

        # Sort the atoms before converting to ensure output is consistent
        # between different runs
        self.sortAtoms()

        atoms = self.vertices
        bonds = self.edges

        obmol = openbabel.OBMol()
        for atom in atoms:
            a = obmol.NewAtom()
            a.SetAtomicNum(atom.number)
            a.SetFormalCharge(atom.charge)
            a.SetVector(atom.coord[0], atom.coord[1], atom.coord[2])
        orders = {'S': 1, 'D': 2, 'T': 3, 'B': 5}
        for atom1, bonds in bonds.iteritems():
            for atom2, bond in bonds.iteritems():
                index1 = atoms.index(atom1)
                index2 = atoms.index(atom2)
                if index1 < index2:
                    order = orders[bond.order]
                    obmol.AddBond(index1+1, index2+1, order)

        obmol.AssignSpinMultiplicity(True)

        # Restore implicit hydrogens if necessary
        if implicitH: self.makeHydrogensImplicit()

        return obmol

    def toJson(self, fname='molecule.json'):
        """ Converts an OpenBabel molecule to json for use in Blender """
        import json
        from json import encoder
        encoder.FLOAT_REPR = lambda o: format(o, '.3f')

        def dumps(object):
            """Outputs json with small formatting edits."""
            # Pretty print json string with truncated floats
            json_string = json.dumps(object, indent=4, sort_keys=True)
            # Make all lists of floats one line and return
            return make_one_line_lists(json_string)


        def make_one_line_lists(json_string):
            """Display float lists as one line in json. Useful for vectors."""
            json_string = json_string.split("\n")
            for i, row in enumerate(json_string):

                # Iterate through all rows that start a list
                if row[-1] != "[" or not has_next_float(json_string, i):
                    continue

                # Move down rows until the list ends, deleting and appending.
                while has_next_float(json_string, i):
                    row += " " + json_string[i + 1].strip()
                    del json_string[i + 1]

                # Finish off with the closing bracket
                json_string[i] = row + " " + json_string[i + 1].strip()
                del json_string[i + 1]

            # Recombine the list into a string and return
            return "\n".join(json_string)

        def has_next_float(json_string, i):
            """Tests if the next row in a split json string is a float."""
            try:
                float(json_string[i + 1].strip().replace(",", ""))
                return True
            except:
                return False

        # Get centroid to center molecule at (0, 0, 0)
        centroid = [0, 0, 0]
        for atom in self.atoms:
            centroid = [c + a for c, a in zip(centroid, atom.coord)]
        centroid = [c / float(len(self.atoms)) for c in centroid]

        # Openbabel atom types have valence ints. Remove those.
        # There are other flags on common atoms (aromatic, .co, am, etc.)
        parse_type = lambda t: t[0] if len(t) > 2 else re.sub("(\d|\W)", "", t)

        # Save atom element type and 3D location.
        atoms = [{"element": atom.element.symbol,
                  "location": [a - c for a, c in zip(atom.coord, centroid)]}
                 for atom in self.atoms]

        # Save number of bonds and indices of endpoint atoms
        # Switch from 1-index to 0-index counting
        bonds = list()
        for x in self.bonds.iteritems():
            for y in x[1].iteritems():
                bonds.append({'source': int(x[0].label)-1, 'target': int(y[0].label)-1, 'order': y[1].order})
        f = open(fname, 'w+')
        f.write(dumps({"atoms": atoms, "bonds": bonds}))
        f.close()

    def getDihedralAngle(self, atoms=[]):
        if len(atoms) == 4:
            atoms = np.array(atoms)
            atoms = atoms - 1
            #http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
            b1 = self.atoms[atoms[1]].coord - self.atoms[atoms[0]].coord
            b2 = self.atoms[atoms[1]].coord - self.atoms[atoms[2]].coord
            b3 = self.atoms[atoms[3]].coord - self.atoms[atoms[2]].coord

            n1 = np.cross(b1, b2) / np.linalg.norm(np.cross(b1, b2))
            n2 = np.cross(b2, b3) / np.linalg.norm(np.cross(b2, b3))

            n_b2 = b2 / np.linalg.norm(b2)
            m1 = np.cross(n1, n_b2)

            x = np.dot(n1, n2)
            y = np.dot(m1, n2)

            rad = np.arctan2(y,x)
            grad = rad * 180.0 / np.pi

            return grad, rad

    def getAngle(self, atoms=[]):
        if len(atoms) == 3:
            atoms = np.array(atoms)
            atoms = atoms - 1

            b1 = self.atoms[atoms[0]].coord - self.atoms[atoms[1]].coord
            b2 = self.atoms[atoms[2]].coord - self.atoms[atoms[1]].coord

            rad = np.arccos(np.dot(b1, b2)/(np.linalg.norm(b1)*np.linalg.norm(b2)))

            grad = rad * 180.0 / np.pi
            return grad, rad

    def getBondLength(self, atoms=[]):
        if len(atoms) == 2:
            atoms = np.array(atoms)
            atoms = atoms - 1

            b1 = self.atoms[atoms[0]].coord - self.atoms[atoms[1]].coord

            return np.linalg.norm(b1)



    def getCenterOfMass(self, atoms=[]):
        """
        Calculate and return the [three-dimensional] position of the center of
        mass of the current geometry. If a list `atoms` of atoms is specified,
        only those atoms will be used to calculate the center of mass.
        Otherwise, all atoms will be used.
        """

        atoms = np.array(atoms)
        center = np.zeros(3, np.float64); mass = 0.0
        if len(atoms) == 0:
            atoms = np.arange(0, len(self.atoms), 1)
        else:
            atoms = atoms - 1

        for i in atoms:
            center += self.atoms[i].realMass * self.atoms[i].coord
            mass += self.atoms[i].realMass
        center /= mass

        return center

    def moveToCenterOfMass(self):
        center = self.getCenterOfMass()
        for atom in self.atoms:
            atom.coord -= center

    def getMomentOfInertiaTensor(self, atoms=[]):
        """
        Calculate and return the moment of inertia tensor for the current
        geometry in kg*m^2. If the coordinates are not at the center of mass,
        they are temporarily shifted there for the purposes of this calculation.
        """
        atoms = np.array(atoms)
        I = np.zeros((3,3), np.float64)
        centerOfMass = self.getCenterOfMass(atoms)
        if len(atoms) == 0:
            atoms = np.arange(0, len(self.atoms), 1)
        else:
            atoms = atoms - 1

        for i in atoms:
            mass = self.atoms[i].realMass * constants.u
            coord = self.atoms[i].coord - centerOfMass
            coord *= 1.0E-10
            I[0,0] += mass * (coord[1] ** 2 + coord[2] ** 2)
            I[1,1] += mass * (coord[0] ** 2 + coord[2] ** 2)
            I[2,2] += mass * (coord[0] ** 2 + coord[1] ** 2)
            I[0,1] -= mass * coord[0] * coord[1]
            I[0,2] -= mass * coord[0] * coord[2]
            I[1,2] -= mass * coord[1] * coord[2]
        I[1,0] = I[0,1]
        I[2,0] = I[0,2]
        I[2,1] = I[1,2]

        return I

    def getPrincipalMomentsOfInertia(self, atoms=[]):
        """
        Calculate and return the principal moments of inertia and corresponding
        principal axes for the current geometry. The moments of inertia are in
        kg*m^2, while the principal axes have unit length.
        """
        I0 = self.getMomentOfInertiaTensor(atoms)
        # Since I0 is real and symmetric, diagonalization is always possible
        I, V = np.linalg.eigh(I0)
        return I, V

    def getRotationalConstants(self, atoms=[]):
        """
        Calculate and return the principal moments of inertia and corresponding
        principal axes for the current geometry. The moments of inertia are in
        kg*m^2, while the principal axes have unit length.
        """
        I, V = self.getPrincipalMomentsOfInertia(atoms)
        rot = constants.h/(8.0*np.pi**2.0*I)
        return rot

    def rotateToPrincipleAxisSystem(self):
        self.moveToCenterOfMass()
        I, V = self.getPrincipalMomentsOfInertia()
        for atom in self.atoms:
            atom.coord = np.dot(V.T, atom.coord)

    def rotateMolecule(self, alpha=0.0, beta=0.0, gamma=0.0):
#        self.moveToCenterOfMass()
        alpha = float(alpha); beta = float(beta); gamma = float(gamma)
        R_x = np.zeros((3,3)); R_y = np.zeros((3,3)); R_z = np.zeros((3,3))

        R_x[0,0] = 1.; R_x[1,1] = np.cos(alpha/180.*np.pi); R_x[2,2] = np.cos(alpha/180.*np.pi)
        R_x[1,2] = -np.sin(alpha/180*np.pi); R_x[2,1] = np.sin(alpha/180.*np.pi)

        R_y[0,0] = np.cos(beta/180.*np.pi); R_y[1,1] = 1.; R_y[2,2] = np.cos(beta/180.*np.pi)
        R_y[0,2] = np.sin(beta/180.*np.pi); R_y[2,0] = -np.sin(beta/180.*np.pi)

        R_z[0,0] = np.cos(gamma/180.*np.pi); R_z[1,1] = np.cos(gamma/180.*np.pi); R_z[2,2] = 1.
        R_z[0,1] = -np.sin(gamma/180.*np.pi); R_z[1,0] = np.sin(gamma/180.*np.pi)

        R = np.dot(R_x,R_y)
        R = np.dot(R, R_z)

        for atom in self.atoms:
            atom.coord = np.dot(R, atom.coord)


    def mirrorMolecule(self, plane='yz'):
        for atom in self.atoms:
            if plane == 'yz':
                atom.coord[0] = -atom.coord[0]
            elif plane == 'xz':
                atom.coord[1] = -atom.coord[1]
            elif plane == 'xy':
                atom.coord[2] = -atom.coord[2]
            else:
                atom.coord[0] = -atom.coord[0]
                atom.coord[1] = -atom.coord[1]
                atom.coord[2] = -atom.coord[2]


    def rotateMoleculeAroundBond(self, bond, alpha=0.0):
        self.moveToCenterOfMass()
        # TODO

    def rotateTopAroundBond(self, bond, alpha=0.0):
        self.moveToCenterOfMass()
        # TODO

    def getBondOrientation(self, atom1, atom2):
        lam = np.zeros(3); vec = np.zeros(3)
        self.rotateToPrincipleAxisSystem()
        if np.linalg.norm(self.atoms[atom1].coord) > np.linalg.norm(self.atoms[atom2].coord):
            vec = self.atoms[atom1].coord - self.atoms[atom2].coord
        else:
            vec = self.atoms[atom2].coord - self.atoms[atom1].coord

        lam = vec/np.linalg.norm(vec)

        return lam, vec

    def getTopOrientation(self, conAtom=1, topAtoms=[]):
        conAtom -= 1
        lam = np.zeros(3); vec = np.zeros(3)
        self.rotateToPrincipleAxisSystem()
        centerOfMass = self.getCenterOfMass(topAtoms)
        norm = np.linalg.norm(self.atoms[conAtom].coord - centerOfMass)

        for i in np.arange(0,3,1):
            vec[i] = centerOfMass[i] - self.atoms[conAtom].coord[i]
        lam = vec / norm

        return lam, vec

    def getIntRotParams(self, conAtom=1, onAxisAtom=1, offAxisAtoms=[], rep='prolate', cof='yes'):

        orientation = dict()
        rho = np.zeros(3)
        self.rotateToPrincipleAxisSystem()

        I_alpha, V = self.getPrincipalMomentsOfInertia(np.append(offAxisAtoms, onAxisAtom))
        rot_alpha = self.getRotationalConstants(np.append(offAxisAtoms, onAxisAtom))
        I_alpha = I_alpha[2]; rot_alpha = rot_alpha[2]

        I, V = self.getPrincipalMomentsOfInertia()
        if cof == 'yes':
            lam, vec = self.getTopOrientation(conAtom, np.append(offAxisAtoms, onAxisAtom))
        else:
            lam, vec = self.getBondOrientation(conAtom-1, onAxisAtom-1)
        rot = self.getRotationalConstants()

        for i in np.arange(0,3,1):
            rho[i] = lam[i]*I_alpha/I[i]
        if rep == 'prolate':
            x = 1
            y = 2
            z = 0
        else:
            x = 1
            y = 0
            z = 2

        orientation['gamma'] = [np.arccos(rho[x] / np.sqrt(rho[x]**2 + rho[y]**2))]
        orientation['gamma'].append(orientation['gamma'][0]*180/np.pi)
        orientation['beta'] = [np.arccos(rho[z] / np.linalg.norm(rho))]
        orientation['beta'].append(orientation['beta'][0]*180/np.pi)
        orientation['delta'] = [np.arccos(lam[z])]
        orientation['delta'].append(orientation['delta'][0]*180/np.pi)
        orientation['epsilon'] = [np.arccos((np.sign(lam[y])*lam[x])/np.sqrt(lam[x]**2 + lam[y]**2))]
        orientation['epsilon'].append(orientation['epsilon'][0]*180/np.pi)
        orientation['lam'] = lam
        orientation['lam_rad'] = np.arccos(lam)
        orientation['lam_deg'] = orientation['lam_rad']*180/np.pi
        orientation['I_alpha'] = I_alpha
        orientation['I'] = I
        orientation['vec'] = vec
        orientation['rho'] = rho
        orientation['rot'] = rot
        orientation['F0'] = rot_alpha
        orientation['norm_rho'] = np.linalg.norm(rho)
        orientation['r'] = 1. - lam[y]**2*I_alpha/I[y] - lam[z]**2*I_alpha/I[z]
        orientation['F'] = constants.h/(8*np.pi**2*orientation['r']*I_alpha)

        conAtom -= 1

        return orientation

    def getIsotopologues(self, rot_meas=np.array([])):

        rot_diff = dict()
        rot_shifted = dict()
        rot_iso = dict()
        iso = dict(C=[13.], O=[18.], N=[15.], Ge=[70.,72.,73.,76.])
        rot_normal = self.getRotationalConstants()

        for i in iso.iterkeys():
            for j in self.atoms:
                if j.symbol == i:
                    for k in iso[i]:
                        j.realMass = k
                        rot_diff[str(j.symbol) + str(j.label) + '_{0:.2f}'.format(k)] = self.getRotationalConstants() - rot_normal
                        rot_iso[str(j.symbol) + str(j.label) + '_{0:.2f}'.format(k)] = self.getRotationalConstants()
                        if len(rot_meas > 0):
                            rot_shifted[str(j.symbol) + str(j.label) + '_{0:.2f}'.format(k)] = self.getRotationalConstants() - rot_normal + rot_meas
                        j.realMass = j.element.mass

        return rot_iso, rot_diff, rot_shifted
