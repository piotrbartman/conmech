"""
Created at 21.08.2019

@author: Piotr Bartman
"""

import numpy as np
from simulation.simulation_runner import SimulationRunner
from simulation.grid.grid import Grid
from simulation.grid.grid_factory import GridFactory
from simulation.solver.solver_factory import SolverFactory
from simulation.solver.solver import Solver
from utils.drawer import Drawer
import numba


class Setup:
    gridHeight = 1
    cells_number = (20, 20)  # number of triangles per aside
    grid_left_border = Grid.DIRICHLET
    grid_top_border = Grid.DIRICHLET
    grid_right_border = Grid.DIRICHLET
    grid_bottom_border = Grid.DIRICHLET

    alpha = 15

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


if __name__ == '__main__':
    setup = Setup()
    grid = GridFactory.construct(setup.cells_number[0],
                                 setup.cells_number[1],
                                 setup.gridHeight,
                                 left=setup.grid_left_border,
                                 top=setup.grid_top_border,
                                 right=setup.grid_right_border,
                                 bottom=Grid.CONTACT)
    B = SolverFactory.construct_B(grid)
    B = B[(1, 1)] + B[(2, 2)]
    ub = np.empty(grid.ind_num)
    for i in range(grid.ind_num):
        p = grid.Points[i]
        ub[i] = setup.ub(p[0], p[1], setup.b)
    Bu = Solver.numba_Bu1(B, ub)

    grid = GridFactory.construct(setup.cells_number[0],
                                 setup.cells_number[1],
                                 setup.gridHeight,
                                 left=setup.grid_left_border,
                                 top=setup.grid_top_border,
                                 right=setup.grid_right_border,
                                 bottom=setup.grid_bottom_border)
    solver = SolverFactory.construct(grid=grid, setup=setup)
    print(Bu[:grid.ind_num])
    solver.Bu = Bu

    solver.solve(verbose=True)

    # TODO draw "Contact" border as b
    Drawer.draw(solver, setup)
