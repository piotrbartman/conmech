"""
Created at 21.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np


class F:
    def __init__(self, grid, F0, FN):
        self.F0 = F0
        self.FN = FN
        self.grid = grid
        self.F = None
        self.Zero = np.zeros(self.grid.ind_num)
        self.One = np.zeros(self.grid.ind_num)

    # TODO: inject?
    ########################################################

    def f0(self, x):
        return self.F0(*x)

    def fN(self, x):
        return self.FN(*x)

    ########################################################

    def setF(self):
        half_long_triangle_side = self.grid.halfLongTriangleSide
        half_short_triangle_side = self.grid.halfShortTriangleSide

        self.F = np.zeros(self.grid.ind_num)

        for i in range(self.grid.ind_num):
            x = self.grid.Points[i][0]
            y = self.grid.Points[i][1]
            t = self.grid.Points[i][2]

            # TODO: static class Point?
            if t != self.grid.CROSS:  # normal point
                # TODO: clean up
                values_in_triangle = np.zeros(8)

                values_in_triangle[0] += self.f0([x - half_long_triangle_side, y])
                values_in_triangle[0] += self.f0([x - half_short_triangle_side, y + half_short_triangle_side])

                values_in_triangle[1] += self.f0([x - half_short_triangle_side, y + half_short_triangle_side])
                values_in_triangle[1] += self.f0([x, y + half_long_triangle_side])

                values_in_triangle[2] += self.f0([x, y + half_long_triangle_side])
                values_in_triangle[2] += self.f0([x + half_short_triangle_side, y + half_short_triangle_side])

                values_in_triangle[3] += self.f0([x + half_short_triangle_side, y + half_short_triangle_side])
                values_in_triangle[3] += self.f0([x + half_long_triangle_side, y])

                values_in_triangle[4] += self.f0([x + half_long_triangle_side, y])
                values_in_triangle[4] += self.f0([x + half_short_triangle_side, y - half_short_triangle_side])

                values_in_triangle[5] += self.f0([x + half_short_triangle_side, y - half_short_triangle_side])
                values_in_triangle[5] += self.f0([x, y - half_long_triangle_side])

                values_in_triangle[6] += self.f0([x, y - half_long_triangle_side])
                values_in_triangle[6] += self.f0([x - half_short_triangle_side, y - half_short_triangle_side])

                values_in_triangle[7] += self.f0([x - half_short_triangle_side, y - half_short_triangle_side])
                values_in_triangle[7] += self.f0([x - half_long_triangle_side, y])

                # TODO: ...
                idxs = ()
                if t == self.grid.TOP:
                    idxs = (4, 5, 6, 7)
                if t == self.grid.LEFT_TOP_CORNER:
                    idxs = (6, 7)
                if t == self.grid.RIGHT_SIDE:
                    idxs = (0, 1, 6, 7)
                if t == self.grid.RIGHT_BOTTOM_CORNER:
                    idxs = (0, 1)
                if t == self.grid.BOTTOM:
                    idxs = (0, 1, 2, 3)
                if t == self.grid.NORMAL_MIDDLE:
                    idxs = (0, 1, 2, 3, 4, 5, 6, 7)

                for idx in idxs:
                    self.F[i] += values_in_triangle[idx]

                self.F[i] = (float(self.grid.TriangleArea) / 6) * self.F[i]

            else:  # cross point

                values_in_triangle = np.zeros(4)

                values_in_triangle[0] += self.f0([x - half_short_triangle_side, y - half_short_triangle_side])
                values_in_triangle[0] += self.f0([x - half_short_triangle_side, y + half_short_triangle_side])

                values_in_triangle[1] += self.f0([x - half_short_triangle_side, y + half_short_triangle_side])
                values_in_triangle[1] += self.f0([x + half_short_triangle_side, y + half_short_triangle_side])

                values_in_triangle[2] += self.f0([x + half_short_triangle_side, y + half_short_triangle_side])
                values_in_triangle[2] += self.f0([x + half_short_triangle_side, y - half_short_triangle_side])

                values_in_triangle[3] += self.f0([x + half_short_triangle_side, y - half_short_triangle_side])
                values_in_triangle[3] += self.f0([x - half_short_triangle_side, y - half_short_triangle_side])

                self.F[i] += values_in_triangle[0]
                self.F[i] += values_in_triangle[1]
                self.F[i] += values_in_triangle[2]
                self.F[i] += values_in_triangle[3]

                self.F[i] = (float(self.grid.TriangleArea) / 6) * self.F[i]

        for i in range(self.grid.ind_num):
            for e in range(-self.grid.BorderEdgesD - self.grid.BorderEdgesN, -self.grid.BorderEdgesD):
                e1 = int(self.grid.Edges[e][0])
                e2 = int(self.grid.Edges[e][1])
                p1 = self.grid.Points[int(e1)][0:2]
                p2 = self.grid.Points[int(e2)][0:2]
                x = (p1 + p2) * 0.5
                if i == e1 or i == e2:
                    self.F[i] += ((self.grid.longTriangleSide * 0.5) * self.fN(x))

        self.Zero = self.F[:]
