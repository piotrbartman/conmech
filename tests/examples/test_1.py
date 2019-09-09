"""
Created at 21.08.2019

@author: Piotr Bartman
"""

from simulation.simulation_runner import SimulationRunner
from examples.usa.only_dphi import Setup, u_infinity, u_alpha


def test_simple_run():
    setup = Setup()
    SimulationRunner.run(setup)


def test_alpha_loop_and_infinity():
    setup = Setup()
    setup.cells_number = (8, 8)
    sim_length = 10
    sim_alphas = [2 ** (i + 0) for i in range(sim_length)]
    sim_quality = [max(8 ** (sim_length / 2 - i / 2), 1) for i in range(sim_length)]
    u_alpha(setup, alphas=sim_alphas, quality=sim_quality)
    u_infinity(setup)
