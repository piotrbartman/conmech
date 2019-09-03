"""
Created at 22.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from simulation.mesh.mesh import Mesh
from simulation.mesh.point import Point
from simulation.mesh.edge import Edge


class MeshFactory:

    @staticmethod
    def construct(sizeH, sizeL, height, left, top, right, bottom):
        mesh = Mesh()
        mesh.SizeH = sizeH
        mesh.SizeL = sizeL
        mesh.Height = float(height)
        mesh.longTriangleSide = mesh.Height / sizeH
        mesh.Length = mesh.longTriangleSide * sizeL

        mesh.halfLongTriangleSide = mesh.longTriangleSide * 0.5
        mesh.shortTriangleSide = mesh.longTriangleSide * np.sqrt(2.) * 0.5
        mesh.halfShortTriangleSide = mesh.shortTriangleSide * 0.5
        mesh.TriangleArea = (mesh.longTriangleSide * mesh.longTriangleSide) / 4.

        MeshFactory.startBorder(mesh, 0, 0)
        for i in range(sizeH):
            x, y = 0, (i + 1) * mesh.longTriangleSide
            MeshFactory.add_border(mesh, x, y, point=Point.LEFT_SIDE, edge_inversion=False, edge=Edge.VERTICAL)
        mesh.Points[0, 2] = Point.LEFT_TOP_CORNER
        mesh.borders[left] += sizeH

        for i in range(sizeL):
            x, y = (i + 1) * mesh.longTriangleSide, mesh.Height
            MeshFactory.add_border(mesh, x, y, point=Point.TOP, edge_inversion=False, edge=Edge.HORIZONTAL)
        mesh.Points[0, 2] = Point.RIGHT_TOP_CORNER
        mesh.borders[top] += sizeL

        for i in range(sizeH, 0, -1):
            x, y = mesh.Length, (i - 1) * mesh.longTriangleSide
            MeshFactory.add_border(mesh, x, y, point=Point.RIGHT_SIDE, edge_inversion=True, edge=Edge.VERTICAL)
        mesh.Points[0, 2] = Point.RIGHT_BOTTOM_CORNER
        mesh.borders[right] += sizeH

        for i in range(sizeL, 1, -1):
            x, y = (i - 1) * mesh.longTriangleSide, 0
            MeshFactory.add_border(mesh, x, y, point=Point.BOTTOM, edge_inversion=True, edge=Edge.HORIZONTAL)
        MeshFactory.stopBorder(mesh)
        mesh.borders[bottom] += sizeL

        for i in range(sizeL):
            for j in range(1, sizeH):
                x1 = i * mesh.longTriangleSide
                x2 = (i + 1) * mesh.longTriangleSide
                y = j * mesh.longTriangleSide
                MeshFactory.addPoint(mesh, x1, y, Point.NORMAL_MIDDLE)
                MeshFactory.addPoint(mesh, x2, y, Point.NORMAL_MIDDLE)
                a = mesh.getPoint(x1, y)
                b = mesh.getPoint(x2, y)
                MeshFactory.addEdge(mesh, a, b, Edge.HORIZONTAL)

        for i in range(1, sizeL):
            for j in range(sizeH):
                x = i * mesh.longTriangleSide
                y1 = j * mesh.longTriangleSide
                y2 = (j + 1) * mesh.longTriangleSide
                a = mesh.getPoint(x, y1)
                b = mesh.getPoint(x, y2)
                MeshFactory.addEdge(mesh, a, b, Edge.VERTICAL)

        for i in range(sizeL):
            for j in range(sizeH):
                x = (i + 0.5) * mesh.longTriangleSide
                y = (j + 0.5) * mesh.longTriangleSide
                MeshFactory.addPoint(mesh, x, y, Point.CROSS)
                a = mesh.getPoint(x, y)
                b = mesh.getPoint(i * mesh.longTriangleSide, (j + 1) * mesh.longTriangleSide)
                MeshFactory.addEdge(mesh, a, b, Edge.BOTTOM)
                b = mesh.getPoint((i + 1) * mesh.longTriangleSide, (j + 1) * mesh.longTriangleSide)
                MeshFactory.addEdge(mesh, a, b, Edge.X_TOP)
                b = mesh.getPoint((i + 1) * mesh.longTriangleSide, j * mesh.longTriangleSide)
                MeshFactory.addEdge(mesh, a, b, Edge.X_BOTTOM)
                b = mesh.getPoint(i * mesh.longTriangleSide, j * mesh.longTriangleSide)
                MeshFactory.addEdge(mesh, a, b, Edge.TOP)

        max_edges = 4
        mesh.edges = np.zeros((np.max(mesh.Edges) + 1, max_edges, 2), dtype=np.int)
        mesh.edges -= 1
        for i in range(len(mesh.Edges)):
            edges = mesh.edges[mesh.Edges[i][0]]
            for j in range(max_edges):
                if edges[j][0] == -1:
                    edges[j][0] = mesh.Edges[i][1]
                    edges[j][1] = mesh.Edges[i][2]
                    break
        return mesh


    @staticmethod
    def startBorder(mesh, x, y):
        MeshFactory.addPoint(mesh, x, y, Point.LEFT_BOTTOM_CORNER)

    @staticmethod
    def add_border(mesh, x, y, point, edge_inversion, edge):
        i, j = 1, 0
        if edge_inversion:
            i, j = j, i
        MeshFactory.addPoint(mesh, x, y, point)
        MeshFactory.addEdge(mesh, i, j, edge)

    @staticmethod
    def stopBorder(mesh):
        MeshFactory.addEdge(mesh, len(mesh.Points) - 1, 0, 1)

    @staticmethod
    def addPoint(mesh, x, y, t):
        i = 0
        while i < len(mesh.Points):
            if mesh.Points[i][0] == x and mesh.Points[i][1] == y:
                return
            else:
                i += 1
        mesh.Points = np.append([[x, y, t]], mesh.Points, axis=0)
        for i in range(len(mesh.Edges)):
            mesh.Edges[i][0] += 1
            mesh.Edges[i][1] += 1

    @staticmethod
    def addEdge(mesh, i, j, t):  # (i,j): xi <= xj and yi < yj
        a = i
        b = j
        if (mesh.Points[j][0] < mesh.Points[i][0] or
                (mesh.Points[j][0] == mesh.Points[i][0] and mesh.Points[j][1] < mesh.Points[i][1])):
            a = j
            b = i
        mesh.Edges = np.append([[a, b, t]], mesh.Edges, axis=0)
