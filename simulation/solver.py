"""
Created at 21.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from simulation.matrices import Matrices
from simulation.f import F


class Solver:

    def __init__(self, grid, F0, FN, alpha, regular_dphi, b, rho):

        self.grid = grid

        self.B = Matrices.construct_B(grid)
        self.F = F(grid, F0, FN)

        self.alpha = alpha
        self.regular_dphi = regular_dphi
        self.b = b
        self.rho = rho

        self.u = np.empty(self.grid.ind_num)

    @staticmethod
    def distance(p1, p2):
        return np.sqrt(np.sum((p1 - p2) ** 2))

    def n_down(self, e):
        # [0,-1]
        e1 = int(self.grid.Edges[e][0])
        e2 = int(self.grid.Edges[e][1])
        dx = self.grid.Points[e2][0] - self.grid.Points[e1][0]
        dy = self.grid.Points[e2][1] - self.grid.Points[e1][1]
        norm = np.sqrt(dx * dx + dy * dy)
        n = np.array([float(dy) / norm, float(-dx) / norm])
        if n[1] > 0:
            n = -n
        return n

    @staticmethod
    def u_at_middle(e1, e2, ind_num, u):
        umL = 0  # u at mL
        if e1 < ind_num:
            umL += u[e1] * 0.5
        if e2 < ind_num:
            umL += u[e2] * 0.5
        return umL

    def JZu(self):
        JZu = np.zeros(self.grid.ind_num)

        for i in range(self.grid.ind_num):
            for e in range(-self.grid.BorderEdgesD - self.grid.BorderEdgesN - self.grid.BorderEdgesC,
                           -self.grid.BorderEdgesD - self.grid.BorderEdgesN):
                e1 = self.grid.Edges[e][0]
                e2 = self.grid.Edges[e][1]
                if i == e1 or i == e2:
                    umL = Solver.u_at_middle(e1, e2, self.grid.ind_num, self.u)

                    p1 = self.grid.Points[e1][:2]
                    p2 = self.grid.Points[e2][:2]

                    L = Solver.distance(p1, p2)

                    JZu[i] += L * 0.5 * self.regular_dphi(umL, self.b, self.rho)

        JZu *= self.alpha

        return JZu

    def Bu1(self):
        result = np.dot((self.B[(1, 1)] + self.B[(2, 2)]), self.u)
        return result

    def f(self, u_vector):
        self.u = u_vector

        X = self.Bu1() \
            + self.JZu() \
            - self.F.Zero

        return 100000000 * X  # 10000000000

    ########################################################

    @staticmethod
    def jtZ(uT, vT, rho=0.0000001):  # uT, vT - vectors; REGULARYZACJA Coulomba
        M = 1 / np.sqrt(float(uT[0] * uT[0] + uT[1] * uT[1]) + float(rho ** 2))
        result = M * float(uT[0]) * float(vT[0]) + M * float(uT[1]) * float(vT[1])
        return result
