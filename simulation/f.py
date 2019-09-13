"""
Created at 21.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from simulation.mesh.point import Point


class F:
    def __init__(self, mesh, F0, FN):
        self.F = F.construct_f_vector(mesh, F0, FN)

    @staticmethod
    def construct_f_vector(mesh, F0, FN):
        half_long_triangle_side = mesh.halfLongTriangleSide
        half_short_triangle_side = mesh.halfShortTriangleSide

        F = np.zeros(mesh.ind_num)

        for i, p in enumerate(mesh.get_independent_points()):
            x = p[0]
            y = p[1]
            t = p[2]

            if t != Point.CROSS:  # normal point
                # TODO: clean up
                values_in_triangle = np.zeros(8)

                values_in_triangle[0] += F0(x - half_long_triangle_side, y)
                values_in_triangle[0] += F0(x - half_short_triangle_side, y + half_short_triangle_side)

                values_in_triangle[1] += F0(x - half_short_triangle_side, y + half_short_triangle_side)
                values_in_triangle[1] += F0(x, y + half_long_triangle_side)

                values_in_triangle[2] += F0(x, y + half_long_triangle_side)
                values_in_triangle[2] += F0(x + half_short_triangle_side, y + half_short_triangle_side)

                values_in_triangle[3] += F0(x + half_short_triangle_side, y + half_short_triangle_side)
                values_in_triangle[3] += F0(x + half_long_triangle_side, y)

                values_in_triangle[4] += F0(x + half_long_triangle_side, y)
                values_in_triangle[4] += F0(x + half_short_triangle_side, y - half_short_triangle_side)

                values_in_triangle[5] += F0(x + half_short_triangle_side, y - half_short_triangle_side)
                values_in_triangle[5] += F0(x, y - half_long_triangle_side)

                values_in_triangle[6] += F0(x, y - half_long_triangle_side)
                values_in_triangle[6] += F0(x - half_short_triangle_side, y - half_short_triangle_side)

                values_in_triangle[7] += F0(x - half_short_triangle_side, y - half_short_triangle_side)
                values_in_triangle[7] += F0(x - half_long_triangle_side, y)

                # TODO: ...
                idxs = ()
                if t == Point.LEFT_TOP_CORNER:  # TODO Point.LEFT_TOP_CORNER:
                    idxs = (4, 5)
                if t == Point.LEFT_SIDE:
                    idxs = (2, 3, 4, 5)
                if t == Point.LEFT_BOTTOM_CORNER:
                    idxs = (2, 3)
                if t == Point.TOP:
                    idxs = (4, 5, 6, 7)
                if t == Point.RIGHT_TOP_CORNER:  # TODO Point.LEFT_TOP_CORNER:
                    idxs = (6, 7)
                if t == Point.RIGHT_SIDE:
                    idxs = (0, 1, 6, 7)
                if t == Point.RIGHT_BOTTOM_CORNER:
                    idxs = (0, 1)
                if t == Point.BOTTOM:
                    idxs = (0, 1, 2, 3)
                if t == Point.NORMAL_MIDDLE:
                    idxs = (0, 1, 2, 3, 4, 5, 6, 7)

                for idx in idxs:
                    F[i] += values_in_triangle[idx]

                F[i] = (float(mesh.TriangleArea) / 6) * F[i]

            else:  # cross point

                values_in_triangle = np.zeros(4)

                values_in_triangle[0] += F0(x - half_short_triangle_side, y - half_short_triangle_side)
                values_in_triangle[0] += F0(x - half_short_triangle_side, y + half_short_triangle_side)

                values_in_triangle[1] += F0(x - half_short_triangle_side, y + half_short_triangle_side)
                values_in_triangle[1] += F0(x + half_short_triangle_side, y + half_short_triangle_side)

                values_in_triangle[2] += F0(x + half_short_triangle_side, y + half_short_triangle_side)
                values_in_triangle[2] += F0(x + half_short_triangle_side, y - half_short_triangle_side)

                values_in_triangle[3] += F0(x + half_short_triangle_side, y - half_short_triangle_side)
                values_in_triangle[3] += F0(x - half_short_triangle_side, y - half_short_triangle_side)

                F[i] += values_in_triangle[0]
                F[i] += values_in_triangle[1]
                F[i] += values_in_triangle[2]
                F[i] += values_in_triangle[3]

                F[i] = (float(mesh.TriangleArea) / 6) * F[i]

        # mesh.get_points()
        for i in range(mesh.ind_num):
            for e in range(-mesh.borders["Dirichlet"] - mesh.borders["Neumann"], -mesh.borders["Dirichlet"]):
                e1 = mesh.Edges[e][0]
                e2 = mesh.Edges[e][1]
                p1 = mesh.Points[e1][:2]
                p2 = mesh.Points[e2][:2]
                x = (p1 + p2) * 0.5
                if i == e1 or i == e2:
                    F[i] += (mesh.longTriangleSide * 0.5) * FN(*x)

        return F
