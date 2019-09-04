"""
Created at 03.09.2019

@author: Piotr Bartman
"""

import numpy as np
from simulation.mesh.dmesh import DMesh
from simulation.mesh.dcell import DCell
from simulation.mesh.dborders import DBorders


class DMeshFactory:
    @staticmethod
    def mesh(mesh_shape, physical_shape, left, top, right, bottom):
        mesh = DMesh()

        mesh.shape = physical_shape

        mesh.points = np.empty((mesh_shape[0] * 2, mesh_shape[1] * 2))

        cell_x_height = physical_shape[0] / (mesh_shape[0] * 2)
        cell_y_height = physical_shape[1] / (mesh_shape[1] * 2)
        mesh.cell = DCell(cell_x_height, cell_y_height)

        mesh.borders = DBorders(left, top, right, bottom)

        mesh.ind_num = DMeshFactory.calculate_independent_number(mesh.x_shape, mesh.y_shape, mesh.borders)

    @staticmethod
    def calculate_independent_number(x, y, borders):
        dirichlet = DBorders.DIRICHLET

        result = (x - 2) * (y - 2)  # inside

        result += y - 2 if borders.left != dirichlet else 0
        result += y - 2 if borders.right != dirichlet else 0
        result += x - 2 if borders.top != dirichlet else 0
        result += x - 2 if borders.bottom != dirichlet else 0

        result += 1 if borders.left != dirichlet and borders.top != dirichlet else 0
        result += 1 if borders.left != dirichlet and borders.bottom != dirichlet else 0
        result += 1 if borders.right != dirichlet and borders.top != dirichlet else 0
        result += 1 if borders.right != dirichlet and borders.bottom != dirichlet else 0
        
        return result


