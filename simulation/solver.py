"""
Created at 21.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from simulation.matrices import Matrices
from simulation.f import F
import numba


class Solver:

    def __init__(self, grid, F0, FN, alpha, regular_dphi, b, rho):

        self.grid = grid

        B = Matrices.construct_B(grid)
        self.B = B[(1, 1)] + B[(2, 2)]
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

    # TODO cleanup
    @staticmethod
    @numba.jit()
    def numba_JZu(ind_num, Edges, Points, u, alpha, BorderEdgesD, BorderEdgesN, BorderEdgesC):
        JZu = np.zeros(ind_num)

        for i in range(ind_num):
            for e in range(-BorderEdgesD - BorderEdgesN - BorderEdgesC,
                           -BorderEdgesD - BorderEdgesN):
                e1 = Edges[e][0]
                e2 = Edges[e][1]
                if i == e1 or i == e2:
                    # umL = Solver.u_at_middle(e1, e2, ind_num, u)

                    p1 = Points[e1][:2]
                    p2 = Points[e2][:2]

                    L = np.sqrt(np.sum((p1 - p2) ** 2)) #Solver.distance(p1, p2)
                    x1 = (p1[0] + p2[0]) * 0.5
                    # TODO
                    JZu[i] += L * 0.5 * (-x1) #self.regular_dphi(umL, b, rho)

        JZu *= alpha
        return JZu

    def JZu(self):
        # JZu = np.zeros(self.grid.ind_num)
        #
        # for i in range(self.grid.ind_num):
        #     for e in range(-self.grid.BorderEdgesD - self.grid.BorderEdgesN - self.grid.BorderEdgesC,
        #                    -self.grid.BorderEdgesD - self.grid.BorderEdgesN):
        #         e1 = self.grid.Edges[e][0]
        #         e2 = self.grid.Edges[e][1]
        #         if i == e1 or i == e2:
        #             umL = Solver.u_at_middle(e1, e2, self.grid.ind_num, self.u)
        #
        #             p1 = self.grid.Points[e1][:2]
        #             p2 = self.grid.Points[e2][:2]
        #
        #             L = Solver.distance(p1, p2)
        #             x1 = (p1[0] + p2[0]) * 0.5
        #             JZu[i] += L * 0.5 * (-x1)#(-np.pi * np.sin(np.pi * (p1[0] + p2[0]) * 0.5)) # self.regular_dphi(umL, self.b, self.rho)
        #
        # JZu *= self.alpha
        JZu = Solver.numba_JZu(self.grid.ind_num, self.grid.Edges, self.grid.Points, self.u, self.alpha,
                               self.grid.BorderEdgesD, self.grid.BorderEdgesN, self.grid.BorderEdgesC)
        return JZu

    # TODO cleanup
    @staticmethod
    @numba.jit()
    def numba_Bu1(B, u):
        result = np.dot(B, u)
        return result

    def Bu1(self):
        # result = np.dot((self.B[(1, 1)] + self.B[(2, 2)]), self.u)
        result = Solver.numba_Bu1(self.B, self.u)
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
