"""
Created at 21.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
import scipy.optimize

from simulation.grid_factory import GridFactory
from utils.drawer import Drawer
from simulation.solver import Solver


class SimulationRunner:
    @staticmethod
    def run(setup):
        grid = GridFactory.construct(setup.cells_number[0],
                                     setup.cells_number[1],
                                     setup.gridHeight)
        solver = Solver(grid, setup.f, setup.g, setup.alpha, setup.regular_dphi, setup.b, setup.rho)

        u_vector = np.zeros(grid.ind_num)

        solver.F.setF()

        while True:
            u_vector = scipy.optimize.fsolve(solver.f, u_vector)
            quality_inv = np.linalg.norm(solver.f(u_vector))
            if quality_inv < 1:
                break
            else:
                print(f"Quality = {quality_inv**-1} is too low, trying again...")

        # solver.set_u_and_displaced_points(u_vector)

        # TODO: separate drawing
        Drawer.draw(solver, setup)
