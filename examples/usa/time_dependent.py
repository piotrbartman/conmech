"""
Created at 21.08.2019

@author: Piotr Bartman
"""

import numpy as np
from pathlib import Path
from simulation.mesh.mesh import Mesh
from simulation.mesh.mesh_factory import MeshFactory
from simulation.solver.solver_factory import SolverFactory
from simulation.solver.solver import Solver
from utils.drawer import Drawer
import numba

SQUARE = 24


class Setup:
    def __init__(self):
        self.gridHeight = 1
        self.cells_number = (SQUARE, SQUARE)  # number of triangles per aside
        self.grid_left_border = Mesh.NEUMANN
        self.grid_top_border = Mesh.NEUMANN
        self.grid_right_border = Mesh.NEUMANN
        self.grid_bottom_border = Mesh.NEUMANN

        self.alpha_0 = 0

        self.time_start = 0
        self.time_stop = 1 / 0.5
        self.time_step = 1 / 32

        # TODO
        self.alpha = None
        self.regular_dphi = None
        self.b = None
        self.rho = None
        self.ub = lambda a, b, c: np.NaN

    @staticmethod
    def F0(x1, x2):
        if 0.4 < x1 < 0.6 and 0.4 < x2 < 0.6:
            result = 2
        else:
            result = 0
        # result = 0.0
        return result

    @staticmethod
    def FN(x1, x2):
        result = 0
        return result


def u_alpha(setup, path=None):
    mesh = MeshFactory.construct(setup.cells_number[0],
                                 setup.cells_number[1],
                                 setup.gridHeight,
                                 left=setup.grid_left_border,
                                 top=setup.grid_top_border,
                                 right=setup.grid_right_border,
                                 bottom=setup.grid_bottom_border)

    start_vector = np.full(mesh.ind_num, setup.alpha_0)
    solver = SolverFactory.construct(mesh=mesh, setup=setup)
    solver.u = start_vector
    time_steps = int((setup.time_stop - setup.time_start) // setup.time_step)

    A = SolverFactory.todo(mesh)
    for _ in range(time_steps):
        b = np.dot(A, solver.u)
        # print(b)
        b = b - setup.time_step * np.dot(2 * solver.B, solver.u)
        b = b + setup.time_step * solver.F.F
        solver.u = np.linalg.solve(A, b)
        print(np.allclose(np.dot(A, solver.u), b))

        Drawer.draw(solver, setup, path=path)


if __name__ == '__main__':
    sim_setup = Setup()
    sim_path = str(Path(__file__).parent.absolute()) + "/results"
    u_alpha(sim_setup)
