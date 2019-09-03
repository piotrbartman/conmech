"""
Created at 29.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from simulation.grid.point import Point
from simulation.grid.edge import Edge
from simulation.f import F
from simulation.solver.solver import Solver


class SolverFactory:
    @staticmethod
    def construct(grid, setup):
        solver = Solver(grid, setup.alpha, setup.regular_dphi, setup.b, setup.rho)

        solver.u = np.empty(grid.ind_num)

        B = SolverFactory.construct_B(grid)
        solver.B = B[(1, 1)] + B[(2, 2)]

        solver.F = F(grid, setup.F0, setup.FN)

        solver.ub = np.empty(grid.ind_num)
        for i in range(grid.ind_num):
            p = grid.Points[i]
            solver.ub[i] = setup.ub(p[0], p[1], setup.b)

        return solver

    @staticmethod
    def construct_B(grid):
        AX = np.zeros((grid.ind_num, 8))  # area with dx
        AY = np.zeros((grid.ind_num, 8))  # area with dy

        grid.get_points(dirichlet=False, neumann=True, contact=True, inside=True)
        for i in range(grid.ind_num):
            p = grid.Points[i]
            AX[i], AY[i] = Point.gradients(p)

        W11 = SolverFactory.multiply(grid, AX, AX)
        W22 = SolverFactory.multiply(grid, AY, AY)

        B = {(1, 1): W11,
             (2, 2): W22
             }

        return B

    @staticmethod
    def multiply(grid, AK, AL):
        W = np.zeros([grid.ind_num, grid.ind_num])

        # grid.get_points()
        for i in range(grid.ind_num):
            W[i][i] = np.sum(AK[i] * AL[i])

            # c - contacting triangles numbers
            for j in range(grid.ind_num):
                edge = grid.get_edge(i, j)

                if edge is not None:  # edge was found
                    c1i, c1j, c2i, c2j = Edge.c(edge)
                    W[i][j] = AK[i][c1i] * AL[j][c1j] + AK[i][c2i] * AL[j][c2j]
                    W[j][i] = AL[i][c1i] * AK[j][c1j] + AL[i][c2i] * AK[j][c2j]

        return W
