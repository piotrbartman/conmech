"""
Created at 21.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from simulation.mesh.point import Point


class F:
    def __init__(self, mesh, F0, FN, t=None):
        self.F = F.construct_f_vector(mesh, F0, FN, t)

    @staticmethod
    def construct_f_vector(mesh, F0, FN, t):
        neighborhood = np. array([[-1, 0], [- .5, .5], [0, 1], [.5, .5], [1, 0], [.5, -.5], [0, -1], [-.5, -.5]])
        neighborhood *= mesh.halfLongTriangleSide

        F = np.zeros(mesh.ind_num)

        for i, p in enumerate(mesh.get_independent_points()):
            x = p[0]
            y = p[1]
            p_type = p[2]

            if p_type != Point.CROSS:  # normal point
                # TODO: ...
                idxs = ()
                if p_type == Point.LEFT_TOP_CORNER:
                    idxs = (4, 5)
                if p_type == Point.LEFT_SIDE:
                    idxs = (2, 3, 4, 5)
                if p_type == Point.LEFT_BOTTOM_CORNER:
                    idxs = (2, 3)
                if p_type == Point.TOP:
                    idxs = (4, 5, 6, 7)
                if p_type == Point.RIGHT_TOP_CORNER:
                    idxs = (6, 7)
                if p_type == Point.RIGHT_SIDE:
                    idxs = (0, 1, 6, 7)
                if p_type == Point.RIGHT_BOTTOM_CORNER:
                    idxs = (0, 1)
                if p_type == Point.BOTTOM:
                    idxs = (0, 1, 2, 3)
                if p_type == Point.NORMAL_MIDDLE:
                    idxs = (0, 1, 2, 3, 4, 5, 6, 7)

                for idx in idxs:
                    for vit_j in range(2):
                        n_idx = (idx + vit_j) % 8
                        if t is None:
                            F[i] += F0(x + neighborhood[n_idx][0], y + neighborhood[n_idx][1])
                        else:
                            F[i] += F0(x + neighborhood[n_idx][0], y + neighborhood[n_idx][1], t)

                F[i] = (float(mesh.TriangleArea) / 6) * F[i]

            else:  # cross point

                idxs = (1, 3, 5, 7)

                for idx in idxs:
                    for vit_j in range(2):
                        n_idx = (idx + 2 * vit_j) % 8
                        if t is None:
                            F[i] += F0(x + neighborhood[n_idx][0], y + neighborhood[n_idx][1])
                        else:
                            F[i] += F0(x + neighborhood[n_idx][0], y + neighborhood[n_idx][1], t)

                F[i] = (float(mesh.TriangleArea) / 6) * F[i]

        # mesh.get_points()
        for i in range(mesh.ind_num):
            for e in range(-mesh.borders["Dirichlet"] - mesh.borders["Neumann"], -mesh.borders["Dirichlet"]):
                e1 = mesh.Edges[e][0]
                e2 = mesh.Edges[e][1]
                p1 = mesh.Points[e1][:2]
                p2 = mesh.Points[e2][:2]
                x = (p1 + p2) * 0.5
                if t is not None:
                    x = (*x, t)
                if i == e1 or i == e2:
                    F[i] += (mesh.longTriangleSide * 0.5) * FN(*x)

        return F
