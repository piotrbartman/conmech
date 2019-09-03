"""
Created at 22.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from simulation.grid.edge import Edge


class Grid:
    DIRICHLET = "Dirichlet"
    NEUMANN = "Neumann"
    CONTACT = "Contact"

    def __init__(self):
        self.Points = np.zeros((0, 3))
        self.points = np.zeros((0, 2))
        self.Edges = np.zeros((0, 3), dtype=np.int)  # TODO optimize (remove and use only edges)
        self.edges = np.zeros((0, 4, 2), dtype=np.int)
        # TODO: bad practice
        #  i, j, type: (always i<j on plane)
        #  0 - no edge
        #  1 - from normal go right to normal, 2 - from normal go up to normal,
        #  3 - from normal go right and up to cross, 4 - from cross go right and up to normal,
        #  5 - from normal go right and down to cross, 6 - from cross go right and down to normal
        #
        self.borders = {"Dirichlet": 0, "Neumann": 0, "Contact": 0}
        self.Height = 0
        self.Length = 0
        self.SizeH = 0
        self.SizeL = 0
        self.longTriangleSide = 0
        self.halfLongTriangleSide = 0
        self.shortTriangleSide = 0
        self.halfShortTriangleSide = 0
        self.TriangleArea = 0

    @property
    def ind_num(self):
        return len(self.Points) - self.borders["Dirichlet"] - 1

    def getPoint(self, x, y):
        i = 0
        while i < len(self.Points):
            if self.Points[i][0] == x and self.Points[i][1] == y:
                return i
            else:
                i += 1
        return -1

    # TODO: order of args still matters
    def get_edge(self, i, j):
        result = None
        for k in range(len(self.edges[i])):
            if self.edges[i, k, 0] == j:
                result = (i, j, self.edges[i, k, 1])
                break

        return result

    def get_independent_points(self):
        return self.get_points(dirichlet=False, neumann=True, contact=True, inside=True)

    def get_points(self, *, dirichlet: bool, neumann: bool, contact: bool, inside: bool):
        if inside:
            start = 0
        elif contact:
            start = len(self.Points) \
                    - self.borders["Dirichlet"] \
                    - self.borders["Neumann"] \
                    - self.borders["Contact"]
        elif neumann:
            start = len(self.Points) \
                    - self.borders["Dirichlet"] \
                    - self.borders["Neumann"]
        elif dirichlet:
            start = len(self.Points) \
                    - self.borders["Dirichlet"]
        else:
            start = len(self.Points)

        if not dirichlet:
            stop = len(self.Points) \
                   - self.borders["Dirichlet"] \
                   - 1
        else:
            stop = len(self.Points)

        return self.Points[start: stop]
