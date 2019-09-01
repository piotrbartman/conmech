"""
Created at 22.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from simulation.grid.grid import Grid
from simulation.grid.point import Point
from simulation.grid.edge import Edge


class GridFactory:

    @staticmethod
    def construct(sizeH, sizeL, height, left, top, right, bottom):
        grid = Grid()
        grid.SizeH = sizeH
        grid.SizeL = sizeL
        grid.Height = float(height)
        grid.longTriangleSide = grid.Height / sizeH
        grid.Length = grid.longTriangleSide * sizeL

        grid.halfLongTriangleSide = grid.longTriangleSide * 0.5
        grid.shortTriangleSide = grid.longTriangleSide * np.sqrt(2.) * 0.5
        grid.halfShortTriangleSide = grid.shortTriangleSide * 0.5
        grid.TriangleArea = (grid.longTriangleSide * grid.longTriangleSide) / 4.

        GridFactory.startBorder(grid, 0, 0)
        for i in range(sizeH):
            x, y = 0, (i + 1) * grid.longTriangleSide
            GridFactory.add_border(grid, x, y, point=Point.LEFT_SIDE, edge_inversion=False, edge=Edge.VERTICAL)
        grid.Points[0, 2] = Point.LEFT_TOP_CORNER
        grid.borders[left] += sizeH

        for i in range(sizeL):
            x, y = (i + 1) * grid.longTriangleSide, grid.Height
            GridFactory.add_border(grid, x, y, point=Point.TOP, edge_inversion=False, edge=Edge.HORIZONTAL)
        grid.Points[0, 2] = Point.RIGHT_TOP_CORNER
        grid.borders[top] += sizeL

        for i in range(sizeH, 0, -1):
            x, y = grid.Length, (i - 1) * grid.longTriangleSide
            GridFactory.add_border(grid, x, y, point=Point.RIGHT_SIDE, edge_inversion=True, edge=Edge.VERTICAL)
        grid.Points[0, 2] = Point.RIGHT_BOTTOM_CORNER
        grid.borders[right] += sizeH

        for i in range(sizeL, 1, -1):
            x, y = (i - 1) * grid.longTriangleSide, 0
            GridFactory.add_border(grid, x, y, point=Point.BOTTOM, edge_inversion=True, edge=Edge.HORIZONTAL)
        GridFactory.stopBorder(grid)
        grid.borders[bottom] += sizeL

        for i in range(sizeL):
            for j in range(1, sizeH):
                x1 = i * grid.longTriangleSide
                x2 = (i + 1) * grid.longTriangleSide
                y = j * grid.longTriangleSide
                GridFactory.addPoint(grid, x1, y, Point.NORMAL_MIDDLE)
                GridFactory.addPoint(grid, x2, y, Point.NORMAL_MIDDLE)
                a = grid.getPoint(x1, y)
                b = grid.getPoint(x2, y)
                GridFactory.addEdge(grid, a, b, Edge.HORIZONTAL)

        for i in range(1, sizeL):
            for j in range(sizeH):
                x = i * grid.longTriangleSide
                y1 = j * grid.longTriangleSide
                y2 = (j + 1) * grid.longTriangleSide
                a = grid.getPoint(x, y1)
                b = grid.getPoint(x, y2)
                GridFactory.addEdge(grid, a, b, Edge.VERTICAL)

        for i in range(sizeL):
            for j in range(sizeH):
                x = (i + 0.5) * grid.longTriangleSide
                y = (j + 0.5) * grid.longTriangleSide
                GridFactory.addPoint(grid, x, y, Point.CROSS)
                a = grid.getPoint(x, y)
                b = grid.getPoint(i * grid.longTriangleSide, (j + 1) * grid.longTriangleSide)
                GridFactory.addEdge(grid, a, b, Edge.BOTTOM)
                b = grid.getPoint((i + 1) * grid.longTriangleSide, (j + 1) * grid.longTriangleSide)
                GridFactory.addEdge(grid, a, b, Edge.X_TOP)
                b = grid.getPoint((i + 1) * grid.longTriangleSide, j * grid.longTriangleSide)
                GridFactory.addEdge(grid, a, b, Edge.X_BOTTOM)
                b = grid.getPoint(i * grid.longTriangleSide, j * grid.longTriangleSide)
                GridFactory.addEdge(grid, a, b, Edge.TOP)

        max_edges = 4
        grid.edges = np.zeros((np.max(grid.Edges) + 1, max_edges, 2), dtype=np.int)
        grid.edges -= 1
        for i in range(len(grid.Edges)):
            edges = grid.edges[grid.Edges[i][0]]
            for j in range(max_edges):
                if edges[j][0] == -1:
                    edges[j][0] = grid.Edges[i][1]
                    edges[j][1] = grid.Edges[i][2]
                    break
        return grid


    @staticmethod
    def startBorder(grid, x, y):
        GridFactory.addPoint(grid, x, y, Point.LEFT_BOTTOM_CORNER)

    @staticmethod
    def add_border(grid, x, y, point, edge_inversion, edge):
        i, j = 1, 0
        if edge_inversion:
            i, j = j, i
        GridFactory.addPoint(grid, x, y, point)
        GridFactory.addEdge(grid, i, j, edge)

    @staticmethod
    def stopBorder(grid):
        GridFactory.addEdge(grid, len(grid.Points) - 1, 0, 1)

    @staticmethod
    def addPoint(grid, x, y, t):
        i = 0
        while i < len(grid.Points):
            if grid.Points[i][0] == x and grid.Points[i][1] == y:
                return
            else:
                i += 1
        grid.Points = np.append([[x, y, t]], grid.Points, axis=0)
        for i in range(len(grid.Edges)):
            grid.Edges[i][0] += 1
            grid.Edges[i][1] += 1

    @staticmethod
    def addEdge(grid, i, j, t):  # (i,j): xi <= xj and yi < yj
        a = i
        b = j
        if (grid.Points[j][0] < grid.Points[i][0] or
                (grid.Points[j][0] == grid.Points[i][0] and grid.Points[j][1] < grid.Points[i][1])):
            a = j
            b = i
        grid.Edges = np.append([[a, b, t]], grid.Edges, axis=0)
