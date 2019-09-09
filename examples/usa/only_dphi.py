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


class Setup:
    gridHeight = 1
    cells_number = (16, 16)  # number of triangles per aside
    grid_left_border = Mesh.DIRICHLET
    grid_top_border = Mesh.DIRICHLET
    grid_right_border = Mesh.DIRICHLET
    grid_bottom_border = Mesh.CONTACT

    alpha = 0

    b = 3
    rho = 1e-3

    @staticmethod
    def F0(x1, x2):
        if 0.4 < x1 < 0.6 and 0.4 < x2 < 0.6:
            result = 128
        else:
            result = 0
        # result = 0
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


def u_infinity(setup):
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

    mesh_special = MeshFactory.construct(setup.cells_number[0],
                                         setup.cells_number[1],
                                         setup.gridHeight,
                                         left=setup.grid_left_border,
                                         top=setup.grid_top_border,
                                         right=setup.grid_right_border,
                                         bottom=Mesh.DIRICHLET)  # !!!
    solver = SolverFactory.construct(mesh=mesh_special, setup=setup)
    solver.condition = lambda: Bu[:mesh_special.ind_num]

    solver.solve(verbose=True)
    solver.u += solver.ub

    solver.mesh = mesh
    path = str(Path(__file__).parent.absolute()) + "/results"
    Drawer.draw(solver, setup, path=path, fixed_contact=True)


def u_alpha(setup, alphas=None, quality=None):
    mesh = MeshFactory.construct(setup.cells_number[0],
                                 setup.cells_number[1],
                                 setup.gridHeight,
                                 left=setup.grid_left_border,
                                 top=setup.grid_top_border,
                                 right=setup.grid_right_border,
                                 bottom=setup.grid_bottom_border)

    if alphas is None:
        alphas = [setup.alpha]
    if quality is None:
        quality = [1000]
    else:
        assert len(quality) == len(alphas)
    start_vector = None
    for i in range(len(alphas)):
        setup.alpha = alphas[i]
        solver = SolverFactory.construct(mesh=mesh, setup=setup)

        solver.solve(start_vector=start_vector, quality=quality[i], verbose=True)

        path = str(Path(__file__).parent.absolute()) + "/results"
        Drawer.draw(solver, setup, path=path)
        start_vector = solver.u


if __name__ == '__main__':
    setup = Setup()
    sim_length = 10
    sim_alphas = [2 ** (i + 0) for i in range(sim_length)]
    sim_quality = [max(8 ** (sim_length/2 - i/2), 1) for i in range(sim_length)]
    u_alpha(setup, alphas=sim_alphas, quality=sim_quality)
    u_infinity(setup)
