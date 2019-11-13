"""
Created at 29.10.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

from new.mesh import Mesh
from new.mesh_factory import Regular2D
import numpy as np
from new.constants import *
import pylab
from matplotlib import pyplot as plt


class TestMesh:

    def test(self):
        mesh = Regular2D.construct(6, 2, left=DIRICHLET, top=DIRICHLET, right=DIRICHLET, bottom=DIRICHLET)

        assert len(mesh.point) == 25
        assert len(mesh.edge) == 60
        assert len(mesh.element) == 36
        assert not np.isnan(mesh.point.coordinates).any()
        assert not np.isnan(mesh.edge.points).any()
        assert not np.isnan(mesh.element.edges).any()

    def test_draw(self):
        mesh = Regular2D.construct(6, 2, left=DIRICHLET, top=DIRICHLET, right=DIRICHLET, bottom=DIRICHLET)

        pylab.axes().set_aspect('equal', 'box')
        plt.scatter(mesh.point[:, 0], mesh.point[:, 1])

        for edge in mesh.edge.points:
            x1, y1 = mesh.point[int(edge[0])][0], mesh.point[int(edge[0])][1]
            x2, y2 = mesh.point[int(edge[1])][0], mesh.point[int(edge[1])][1]
            plt.plot([x1, x2], [y1, y2], 'k-', lw=0.5)
        plt.show()

    def test_draw_border(self):
        mesh = Regular2D.construct(6, 2, left=DIRICHLET, top=DIRICHLET, right=DIRICHLET, bottom=DIRICHLET)

        pylab.axes().set_aspect('equal', 'box')
        plt.scatter(mesh.point[:, 0], mesh.point[:, 1])

        for border in mesh.subarea['border'][DIRICHLET]:
            for edge in border:
                x1, y1 = mesh.point[int(edge[0])][0], mesh.point[int(edge[0])][1]
                x2, y2 = mesh.point[int(edge[1])][0], mesh.point[int(edge[1])][1]
                plt.plot([x1, x2], [y1, y2], 'k-', lw=0.5)
        plt.show()
