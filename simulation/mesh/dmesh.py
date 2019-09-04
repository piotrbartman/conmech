"""
Created at 03.09.2019

@author: Piotr Bartman
"""

import numpy as np
from simulation.mesh.dcell import DCell
from simulation.mesh.dborders import DBorders


class DMesh:

    def __init__(self):
        self.shape: tuple = None
        self.points: np.ndarray = None
        self.cell: DCell = None
        self.borders: DBorders = None
        self.ind_num = 0

    def get_independent_points(self):
        x_start = 0 if self.borders.left != DBorders.DIRICHLET else 1
        x_stop = self.x_shape if self.borders.right != DBorders.DIRICHLET else self.x_shape - 1
        y_start = 0 if self.borders.top != DBorders.DIRICHLET else 1
        y_stop = self.y_shape if self.borders.bottom != DBorders.DIRICHLET else self.y_shape - 1

        result = self.points[x_start:x_stop, y_start:y_stop].reshape((-1,))

        return result

    @property
    def x_shape(self):
        return self.points.shape[0]

    @property
    def y_shape(self):
        return self.points.shape[1]
