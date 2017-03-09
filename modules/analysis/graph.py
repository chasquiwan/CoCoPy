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

"""
This module contains an implementation of a graph data structure (the 
:class:`Graph` class) and functions for manipulating that graph, including 
efficient isomorphism functions.
"""


################################################################################

class Vertex(object):
    """
    A base class for vertices in a graph. Contains several connectivity values
    useful for accelerating isomorphism searches, as proposed by
    `Morgan (1965) <http://dx.doi.org/10.1021/c160017a018>`_.

    ==================  ========================================================
    Attribute           Description
    ==================  ========================================================
    `connectivity1`     The number of nearest neighbors
    `connectivity2`     The sum of the neighbors' `connectivity1` values
    `connectivity3`     The sum of the neighbors' `connectivity2` values
    `sortingLabel`      An integer used to sort the vertices
    ==================  ========================================================

    """

    def __init__(self):
        self.resetConnectivityValues()

    def equivalent(self, other):
        """
        Return :data:`True` if two vertices `self` and `other` are semantically
        equivalent, or :data:`False` if not. You should reimplement this
        function in a derived class if your vertices have semantic information.
        """
        return True

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is semantically more specific than `other`,
        or ``False`` if not. You should reimplement this function in a derived
        class if your edges have semantic information.
        """
        return True

    def resetConnectivityValues(self):
        """
        Reset the cached structure information for this vertex.
        """
        self.connectivity1 = -1
        self.connectivity2 = -1
        self.connectivity3 = -1
        self.sortingLabel = -1

def getVertexConnectivityValue(vertex):
    """
    Return a value used to sort vertices prior to poposing candidate pairs in
    :meth:`__VF2_pairs`. The value returned is based on the vertex's
    connectivity values (and assumes that they are set properly).
    """
    return ( -256*vertex.connectivity1 - 16*vertex.connectivity2 - vertex.connectivity3 )

def getVertexSortingLabel(vertex):
    """
    Return a value used to sort vertices prior to poposing candidate pairs in
    :meth:`__VF2_pairs`. The value returned is based on the vertex's
    connectivity values (and assumes that they are set properly).
    """
    return vertex.sortingLabel

################################################################################

class Edge(object):
    """
    A base class for edges in a graph. This class does *not* store the vertex
    pair that comprises the edge; that functionality would need to be included
    in the derived class.
    """

    def __init__(self):
        pass

    def equivalent(self, other):
        """
        Return ``True`` if two edges `self` and `other` are semantically
        equivalent, or ``False`` if not. You should reimplement this
        function in a derived class if your edges have semantic information.
        """
        return True

    def isSpecificCaseOf(self, other):
        """
        Return ``True`` if `self` is semantically more specific than `other`,
        or ``False`` if not. You should reimplement this function in a derived
        class if your edges have semantic information.
        """
        return True

################################################################################

class Graph:
    """
    A graph data type. The vertices of the graph are stored in a list
    `vertices`; this provides a consistent traversal order. The edges of the
    graph are stored in a dictionary of dictionaries `edges`. A single edge can
    be accessed using ``graph.edges[vertex1][vertex2]`` or the :meth:`getEdge`
    method; in either case, an exception will be raised if the edge does not
    exist. All edges of a vertex can be accessed using ``graph.edges[vertex]``
    or the :meth:`getEdges` method.
    """

    def __init__(self, vertices=None, edges=None):
        self.vertices = vertices or []
        self.edges = edges or {}
        
    def addVertex(self, vertex):
        """
        Add a `vertex` to the graph. The vertex is initialized with no edges.
        """
        self.vertices.append(vertex)
        self.edges[vertex] = dict()
        return vertex

    def addEdge(self, vertex1, vertex2, edge):
        """
        Add an `edge` to the graph as an edge connecting the two vertices
        `vertex1` and `vertex2`.
        """
        self.edges[vertex1][vertex2] = edge
        self.edges[vertex2][vertex1] = edge
        return edge

    def getEdges(self, vertex):
        """
        Return a list of the edges involving the specified `vertex`.
        """
        return self.edges[vertex]

    def getEdge(self, vertex1, vertex2):
        """
        Returns the edge connecting vertices `vertex1` and `vertex2`.
        """
        return self.edges[vertex1][vertex2]

    def hasVertex(self, vertex):
        """
        Returns ``True`` if `vertex` is a vertex in the graph, or ``False`` if
        not.
        """
        return vertex in self.vertices

    def hasEdge(self, vertex1, vertex2):
        """
        Returns ``True`` if vertices `vertex1` and `vertex2` are connected
        by an edge, or ``False`` if not.
        """
        return vertex2 in self.edges[vertex1] if vertex1 in self.edges else False

    def removeVertex(self, vertex):
        """
        Remove `vertex` and all edges associated with it from the graph. Does
        not remove vertices that no longer have any edges as a result of this
        removal.
        """
        for vertex2 in self.vertices:
            if vertex2 is not vertex:
                if vertex in self.edges[vertex2]:
                    del self.edges[vertex2][vertex]
        del self.edges[vertex]
        self.vertices.remove(vertex)

    def removeEdge(self, vertex1, vertex2):
        """
        Remove the edge having vertices `vertex1` and `vertex2` from the graph.
        Does not remove vertices that no longer have any edges as a result of
        this removal.
        """
        del self.edges[vertex1][vertex2]
        del self.edges[vertex2][vertex1]

    def copy(self, deep=False):
        """
        Create a copy of the current graph. If `deep` is ``True``, a deep copy
        is made: copies of the vertices and edges are used in the new graph.
        If `deep` is ``False`` or not specified, a shallow copy is made: the
        original vertices and edges are used in the new graph.
        """
        other = Graph()
        for vertex in self.vertices:
            other.addVertex(vertex.copy() if deep else vertex)
        for vertex1 in self.vertices:
            for vertex2 in self.edges[vertex1]:
                if deep:
                    index1 = self.vertices.index(vertex1)
                    index2 = self.vertices.index(vertex2)
                    other.addEdge(other.vertices[index1], other.vertices[index2],
                    self.edges[vertex1][vertex2].copy())
                else:
                    other.addEdge(vertex1, vertex2, self.edges[vertex1][vertex2])
        return other

    def merge(self, other):
        """
        Merge two graphs so as to store them in a single Graph object.
        """

        # Create output graph
        new = Graph()

        # Add vertices to output graph
        for vertex in self.vertices:
            new.addVertex(vertex)
        for vertex in other.vertices:
            new.addVertex(vertex)

        # Add edges to output graph
        for v1 in self.vertices:
            for v2 in self.edges[v1]:
                new.edges[v1][v2] = self.edges[v1][v2]
        for v1 in other.vertices:
            for v2 in other.edges[v1]:
                new.edges[v1][v2] = other.edges[v1][v2]

        return new

    def split(self):
        """
        Convert a single Graph object containing two or more unconnected graphs
        into separate graphs.
        """

        # Create potential output graphs
        new1 = self.copy()
        new2 = Graph()

        if len(self.vertices) == 0:
            return [new1]

        # Arbitrarily choose last atom as starting point
        verticesToMove = [ self.vertices[-1] ]

        # Iterate until there are no more atoms to move
        index = 0
        while index < len(verticesToMove):
            for v2 in self.edges[verticesToMove[index]]:
                if v2 not in verticesToMove:
                    verticesToMove.append(v2)
            index += 1

        # If all atoms are to be moved, simply return new1
        if len(new1.vertices) == len(verticesToMove):
            return [new1]

        # Copy to new graph
        for vertex in verticesToMove:
            new2.addVertex(vertex)
        for v1 in verticesToMove:
            for v2, edge in new1.edges[v1].iteritems():
                new2.edges[v1][v2] = edge

        # Remove from old graph
        for v1 in new2.vertices:
            for v2 in new2.edges[v1]:
                if v1 in verticesToMove and v2 in verticesToMove:
                    del new1.edges[v1][v2]
        for vertex in verticesToMove:
            new1.removeVertex(vertex)

        new = [new2]
        new.extend(new1.split())
        return new

    def resetConnectivityValues(self):
        """
        Reset any cached connectivity information. Call this method when you
        have modified the graph.
        """
        for vertex in self.vertices: vertex.resetConnectivityValues()
        
    def updateConnectivityValues(self):
        """
        Update the connectivity values for each vertex in the graph. These are
        used to accelerate the isomorphism checking.
        """

        assert str(self.__class__) != 'chempy.molecule.Molecule' or not self.implicitHydrogens, "%s has implicit hydrogens" % self

        for vertex1 in self.vertices:
            count = len(self.edges[vertex1])
            vertex1.connectivity1 = count
        for vertex1 in self.vertices:
            count = 0
            edges = self.edges[vertex1]
            for vertex2 in edges: count += vertex2.connectivity1
            vertex1.connectivity2 = count
        for vertex1 in self.vertices:
            count = 0
            edges = self.edges[vertex1]
            for vertex2 in edges: count += vertex2.connectivity2
            vertex1.connectivity3 = count
        
    def sortVertices(self):
        """
        Sort the vertices in the graph. This can make certain operations, e.g.
        the isomorphism functions, much more efficient.
        """
        # Only need to conduct sort if there is an invalid sorting label on any vertex
        for vertex in self.vertices:
            if vertex.sortingLabel < 0: break
        else:
            return
        self.vertices.sort(key=getVertexConnectivityValue)
        for index, vertex in enumerate(self.vertices):
            vertex.sortingLabel = index
