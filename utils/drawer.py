"""
Created at 21.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import matplotlib.pyplot as plt
import pylab
import numpy as np


class Drawer:
    @staticmethod
    def draw(solver, setup):
        mesh = solver.mesh
        txt = 'CROSS EQUATION GR' + str(mesh.SizeH) + ' ' + str(mesh.SizeL) \
              + ') F0[' + str(setup.F0) \
              + '] FN[' + str(setup.FN) + ']'

        plt.close()
        pylab.axes().set_aspect('equal', 'box')

        shadow = 0.1
        thickness1 = thickness2 = 2

        i = len(mesh.Edges) - 1
        j = len(mesh.Edges) - mesh.borders["Dirichlet"] - 1
        while j < i:
            x1 = mesh.Points[int(mesh.Edges[i, 0])][0]
            y1 = mesh.Points[int(mesh.Edges[i, 0])][1]
            x2 = mesh.Points[int(mesh.Edges[i, 1])][0]
            y2 = mesh.Points[int(mesh.Edges[i, 1])][1]
            plt.plot([x1, x2], [y1, y2], 'k-', alpha=shadow, lw=thickness1)
            i -= 1
        j -= mesh.borders["Neumann"]
        while j < i:
            x1 = mesh.Points[int(mesh.Edges[i, 0])][0]
            y1 = mesh.Points[int(mesh.Edges[i, 0])][1]
            x2 = mesh.Points[int(mesh.Edges[i, 1])][0]
            y2 = mesh.Points[int(mesh.Edges[i, 1])][1]
            plt.plot([x1, x2], [y1, y2], 'k-', alpha=shadow, lw=thickness1)
            i -= 1
        j -= mesh.borders["Contact"]
        while j < i:
            x1 = mesh.Points[int(mesh.Edges[i, 0])][0]
            y1 = mesh.Points[int(mesh.Edges[i, 0])][1]
            x2 = mesh.Points[int(mesh.Edges[i, 1])][0]
            y2 = mesh.Points[int(mesh.Edges[i, 1])][1]
            plt.plot([x1, x2], [y1, y2], 'k-', alpha=shadow, lw=thickness1)
            i -= 1
        while -1 < i:
            x1 = mesh.Points[int(mesh.Edges[i, 0])][0]
            y1 = mesh.Points[int(mesh.Edges[i, 0])][1]
            x2 = mesh.Points[int(mesh.Edges[i, 1])][0]
            y2 = mesh.Points[int(mesh.Edges[i, 1])][1]
            plt.plot([x1, x2], [y1, y2], 'k-', alpha=shadow, lw=thickness1)
            i -= 1

            # ------------
        u = np.concatenate((solver.u, np.zeros(solver.mesh.borders["Dirichlet"] + 1)))
        plt.scatter(mesh.Points[:, 0], mesh.Points[:, 1], marker='o', c=u, cmap="Reds")

        plt.colorbar()

        # i = len(mesh.Edges) - 1
        # j = len(mesh.Edges) - mesh.borders["Dirichlet"] - 1
        # while j < i:
        #     x1 = solver.DisplacedPoints[int(mesh.Edges[i, 0])][0]
        #     y1 = solver.DisplacedPoints[int(mesh.Edges[i, 0])][1]
        #     x2 = solver.DisplacedPoints[int(mesh.Edges[i, 1])][0]
        #     y2 = solver.DisplacedPoints[int(mesh.Edges[i, 1])][1]
        #     plt.plot([x1, x2], [y1, y2], 'r-', lw=thickness2)
        #     i -= 1
        # j -= mesh.borders["Neumann"]
        # while j < i:
        #     x1 = solver.DisplacedPoints[int(mesh.Edges[i, 0])][0]
        #     y1 = solver.DisplacedPoints[int(mesh.Edges[i, 0])][1]
        #     x2 = solver.DisplacedPoints[int(mesh.Edges[i, 1])][0]
        #     y2 = solver.DisplacedPoints[int(mesh.Edges[i, 1])][1]
        #     plt.plot([x1, x2], [y1, y2], 'b-', lw=thickness2)
        #     i -= 1
        # j -= mesh.borders["Contact"]
        # while j < i:
        #     x1 = solver.DisplacedPoints[int(mesh.Edges[i, 0])][0]
        #     y1 = solver.DisplacedPoints[int(mesh.Edges[i, 0])][1]
        #     x2 = solver.DisplacedPoints[int(mesh.Edges[i, 1])][0]
        #     y2 = solver.DisplacedPoints[int(mesh.Edges[i, 1])][1]
        #     plt.plot([x1, x2], [y1, y2], 'y-', lw=thickness2)
        #     i -= 1
        # while -1 < i:
        #     x1 = solver.DisplacedPoints[int(mesh.Edges[i, 0])][0]
        #     y1 = solver.DisplacedPoints[int(mesh.Edges[i, 0])][1]
        #     x2 = solver.DisplacedPoints[int(mesh.Edges[i, 1])][0]
        #     y2 = solver.DisplacedPoints[int(mesh.Edges[i, 1])][1]
        #     plt.plot([x1, x2], [y1, y2], 'k-', lw=thickness2)
        #     i -= 1

            # ------------

        # plt.savefig(txt + '.png', transparent=True, bbox_inches='tight', pad_inches=0, dpi=300)  # DPI 500
        # print(txt + '.png')
        plt.show()
        plt.close()
