"""
Created at 13.09.2019

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

SQUARE = 12


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
    def F0(x1, x2, t):
        result = (- np.e ** -t) * ((x1**3)/3 - (x1**2)/2) - 1 * (2*x1 - 1)
        return result

    @staticmethod
    def FN(x1, x2, t):
        result = 0
        return result

    
def exact_solution(x1, x2, t):
    return (np.e ** -t) * ((x1 ** 3) / 3 - (x1 ** 2) / 2)


def u_alpha(setup, path=None):
    mesh = MeshFactory.construct(setup.cells_number[0],
                                 setup.cells_number[1],
                                 setup.gridHeight,
                                 left=setup.grid_left_border,
                                 top=setup.grid_top_border,
                                 right=setup.grid_right_border,
                                 bottom=setup.grid_bottom_border)

    start_vector = np.full(mesh.ind_num, setup.alpha_0)
    solver = SolverFactory.construct(mesh=mesh, setup=setup, t=0)
    solver.u = start_vector
    time_steps = int((setup.time_stop - setup.time_start) // setup.time_step)

    A = SolverFactory.todo(mesh)
    # A = np.eye(mesh.ind_num)
    for curr_time_step in range(time_steps):
        b = np.dot(A, solver.u)
        # print(b)
        b = b - 1 * setup.time_step * np.dot(solver.B, solver.u)
        if curr_time_step < 6:
            b = b + setup.time_step * solver.F.F
        solver.u = np.linalg.solve(A, b)
        print(np.allclose(np.dot(A, solver.u), b))

        Drawer.draw(solver, setup, path=path)


if __name__ == '__main__':
    sim_setup = Setup()
    sim_path = str(Path(__file__).parent.absolute()) + "/results"
    u_alpha(sim_setup)
