"""
Created at 29.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from simulation.mesh.point import Point
from simulation.mesh.edge import Edge
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
        for i in range(mesh.ind_num):
            p = mesh.Points[i]
            solver.ub[i] = setup.ub(p[0], p[1], setup.b)

        return solver

    @staticmethod
    def construct_B(mesh):
        AX = np.zeros((mesh.ind_num, 8))  # area with dx
        AY = np.zeros((mesh.ind_num, 8))  # area with dy

        mesh.get_points(dirichlet=False, neumann=True, contact=True, inside=True)
        for i in range(mesh.ind_num):
            p = mesh.Points[i]
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

    @staticmethod
    def todo(mesh):
        Q = np.zeros((mesh.ind_num, mesh.ind_num))

        # TODO fix mesh.get_independent_points()
        # for i, p in enumerate(mesh.get_independent_points()):
        aside = 1  # mesh.shortTriangleSide
        for i, p in enumerate(mesh.Points):
            Q[i][i] = (1 / 12) * (aside ** 2) * Point.triangles(p).sum()

            # c - contacting triangles numbers
            for j in range(mesh.ind_num):
                edge = mesh.get_edge(i, j)

                if edge is not None:  # edge was found
                    Q[i][j] = (1 / 12) * (aside ** 2)
                    if not Point.is_inside(mesh.Points[i]) and not Point.is_inside(mesh.Points[j]):
                        Q[i][j] /= 2
                    Q[j][i] = Q[i][j]

        return Q
