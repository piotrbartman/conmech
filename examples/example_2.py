"""
Created at 21.08.2019

@author: Piotr Bartman
"""

import numpy as np
from simulation.simulation_runner import SimulationRunner


class Setup:
    gridHeight = 1

    cells_number = (2, 5)  # number of triangles per aside
    F0 = np.array([-0.2, -0.2])  # inner forces
    FN = np.array([0, 0])  # outer forces


if __name__ == '__main__':
    setup = Setup()
    SimulationRunner.run(setup)
