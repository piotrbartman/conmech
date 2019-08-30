"""
Created at 21.08.2019

@author: Piotr Bartman
"""

import numpy as np
from simulation.simulation_runner import SimulationRunner


class Setup:
    gridHeight = 1
    cells_number = (30, 30)  # number of triangles per aside

    alpha = 1000

    b = 0
    rho = 1e-8

    @staticmethod
    def f(x1, x2):
        result = 2 * np.pi ** 2 + np.sin(np.pi * x1) * np.sin(np.pi * x2)
        return result

    @staticmethod
    def g(x1, x2):
        if x1 == 1:
            result = - np.pi * np.sin(np.pi * x2)
        if x2 == 1:
            result = - np.pi * np.sin(np.pi * x1)
        return result

    @staticmethod
    def regular_dphi(r, b, rho):
        # x = r - b
        # result = x / np.sqrt(x**2 + rho**2)
        result = 0
        return result


if __name__ == '__main__':
    setup = Setup()
    SimulationRunner.run(setup)
