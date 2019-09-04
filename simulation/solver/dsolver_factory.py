"""
Created at 04.09.2019

@author: Piotr Bartman
"""

import numpy as np
from simulation.mesh.point import Point
from simulation.mesh.edge import Edge
from simulation.mesh.dmesh import DMesh
from simulation.f import F
from simulation.solver.solver import Solver


class SolverFactory:
    @staticmethod
    def construct(mesh, setup):
        solver = Solver(mesh, setup.alpha, setup.regular_dphi, setup.b, setup.rho)

        solver.u = np.empty(mesh.ind_num)

        B = SolverFactory.construct_B(mesh)
        solver.B = B[(1, 1)] + B[(2, 2)]

        solver.F = F(mesh, setup.F0, setup.FN)

        solver.ub = np.empty(mesh.ind_num)
        for i, p in enumerate(mesh.get_independent_points()):
            solver.ub[i] = setup.ub(p[0], p[1], setup.b)

        return solver

    @staticmethod
    def construct_B(mesh: DMesh):
        AX = np.zeros((mesh.ind_num, 16))  # area with dx
        AY = np.zeros((mesh.ind_num, 16))  # area with dy

        for i, p in enumerate(mesh.get_independent_points()):
            AX[i], AY[i] = Point.gradients(p)

        W11 = SolverFactory.multiply(mesh, AX, AX)
        W22 = SolverFactory.multiply(mesh, AY, AY)

        B = {(1, 1): W11,
             (2, 2): W22
             }

        return B

    @staticmethod
    def multiply(mesh, AK, AL):
        W = np.zeros([mesh.ind_num, mesh.ind_num])

        # mesh.get_points()
        for i in range(mesh.ind_num):
            W[i][i] = np.sum(AK[i] * AL[i])

            # c - contacting triangles numbers
            for j in range(mesh.ind_num):
                edge = mesh.get_edge(i, j)

                if edge is not None:  # edge was found
                    c1i, c1j, c2i, c2j = Edge.c(edge)
                    W[i][j] = AK[i][c1i] * AL[j][c1j] + AK[i][c2i] * AL[j][c2j]
                    W[j][i] = AL[i][c1i] * AK[j][c1j] + AL[i][c2i] * AK[j][c2j]

        return W
