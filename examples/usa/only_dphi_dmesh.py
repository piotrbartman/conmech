"""
Created at 21.08.2019

@author: Piotr Bartman
"""

import numpy as np
from simulation.mesh.mesh import Mesh
from simulation.mesh.dmesh_factory import DMeshFactory
from simulation.mesh.mesh_factory import MeshFactory
from simulation.solver.solver_factory import SolverFactory
from simulation.solver.solver import Solver
from utils.drawer import Drawer
import numba


class Setup:
    gridHeight = 1
    cells_number = (16, 16)  # number of triangles per aside
    grid_left_border = Mesh.DIRICHLET
    grid_top_border = Mesh.DIRICHLET
    grid_right_border = Mesh.DIRICHLET
    grid_bottom_border = Mesh.CONTACT

    alpha = 18

    b = 3
    rho = 1e-3

    @staticmethod
    def F0(x1, x2):
        result = 0
        return result

    @staticmethod
    def FN(x1, x2):
        result = 0
        return result

    @staticmethod
    @numba.jit()
    def regular_dphi(r, b, rho):
        x = r - b
        result = x / np.sqrt(x ** 2 + rho ** 2)
        return result

    @staticmethod
    @numba.jit()
    def ub(x1, x2, b):
        if x2 == 0:
            result = b
        else:
            result = 0
        return result


def u_alpha():
    setup = Setup()
    mesh = DMeshFactory.construct(setup.cells_number[0],
                                  setup.cells_number[1],
                                  setup.gridHeight,
                                  left=setup.grid_left_border,
                                  top=setup.grid_top_border,
                                  right=setup.grid_right_border,
                                  bottom=setup.grid_bottom_border)
    solver = SolverFactory.construct(mesh=mesh, setup=setup)

    solver.solve(verbose=True)

    Drawer.draw(solver, setup)


if __name__ == '__main__':
    u_alpha()
