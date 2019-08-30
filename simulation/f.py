"""
Created at 21.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from simulation.grid.point import Point


class F:
    def __init__(self, grid, F0, FN):
        self.F = F.construct_f_vector(grid, F0, FN)

    @staticmethod
    def construct_f_vector(grid, F0, FN):
        half_long_triangle_side = grid.halfLongTriangleSide
        half_short_triangle_side = grid.halfShortTriangleSide

        F = np.zeros(grid.ind_num)

        # grid.get_points()
        for i in range(grid.ind_num):
            x = grid.Points[i][0]
            y = grid.Points[i][1]
            t = grid.Points[i][2]

            if t != grid.CROSS:  # normal point
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
                if t == Point.TOP:
                    idxs = (4, 5, 6, 7)
                if t == Point.LEFT_TOP_CORNER:
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

                F[i] = (float(grid.TriangleArea) / 6) * F[i]

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

                F[i] = (float(grid.TriangleArea) / 6) * F[i]

        # grid.get_points()
        for i in range(grid.ind_num):
            for e in range(-grid.BorderEdgesD - grid.BorderEdgesN, -grid.BorderEdgesD):
                e1 = int(grid.Edges[e][0])
                e2 = int(grid.Edges[e][1])
                p1 = grid.Points[int(e1)][0:2]
                p2 = grid.Points[int(e2)][0:2]
                x = (p1 + p2) * 0.5
                if i == e1 or i == e2:
                    F[i] += (grid.longTriangleSide * 0.5) * FN(*x)

        return F
