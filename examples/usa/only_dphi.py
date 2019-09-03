"""
Created at 21.08.2019

@author: Piotr Bartman
"""

import numpy as np
from simulation.mesh.mesh import Mesh
from simulation.mesh.mesh_factory import MeshFactory
from simulation.solver.solver_factory import SolverFactory
from simulation.solver.solver import Solver
from utils.drawer import Drawer
import numba


class Setup:
    gridHeight = 1
    cells_number = (20, 20)  # number of triangles per aside
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
        result = x / np.sqrt(x**2 + rho**2)
        return result

    @staticmethod
    @numba.jit()
    def ub(x1, x2, b):
        if x2 == 0:
            result = b
        else:
            result = 0
        return result


def u_infinity():
    setup = Setup()
    mesh = MeshFactory.construct(setup.cells_number[0],
                                 setup.cells_number[1],
                                 setup.gridHeight,
                                 left=setup.grid_left_border,
                                 top=setup.grid_top_border,
                                 right=setup.grid_right_border,
                                 bottom=setup.grid_bottom_border)
    B = SolverFactory.construct_B(mesh)
    B = B[(1, 1)] + B[(2, 2)]
    ub = np.empty(mesh.ind_num)
    for i in range(mesh.ind_num):
        p = mesh.Points[i]
        ub[i] = setup.ub(p[0], p[1], setup.b)
    Bu = Solver.numba_Bu1(B, ub)

    mesh = MeshFactory.construct(setup.cells_number[0],
                                 setup.cells_number[1],
                                 setup.gridHeight,
                                 left=setup.grid_left_border,
                                 top=setup.grid_top_border,
                                 right=setup.grid_right_border,
                                 bottom=Mesh.DIRICHLET)  # !!!
    solver = SolverFactory.construct(mesh=mesh, setup=setup)
    solver.condition = lambda: Bu[:mesh.ind_num]

    solver.solve(verbose=True)
    solver.u += solver.ub

    Drawer.draw(solver, setup)


def u_alpha():
    setup = Setup()
    mesh = MeshFactory.construct(setup.cells_number[0],
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
    u_infinity()
