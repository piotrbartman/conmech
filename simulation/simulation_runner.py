"""
Created at 21.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

from simulation.grid.grid_factory import GridFactory
from utils.drawer import Drawer
from simulation.solver.solver_factory import SolverFactory
import time


class SimulationRunner:
    @staticmethod
    def run(setup):
        print("Setups...")
        t0 = time.perf_counter()
        t = time.perf_counter()

        # Create grid and solver
        grid = GridFactory.construct(setup.cells_number[0],
                                     setup.cells_number[1],
                                     setup.gridHeight)
        solver = SolverFactory.construct(grid, setup.f, setup.g, setup.alpha, setup.regular_dphi, setup.b, setup.rho)

        print(f"Done after {(time.perf_counter() - t)} s")
        print("Solving...")
        t = time.perf_counter()

        # Solve problem
        solver.solve(verbose=True)

        print(f"Done after {(time.perf_counter() - t)} s")
        # TODO: separate drawing
        print("Drawing...")
        t = time.perf_counter()

        # Draw results
        # solver.set_u_and_displaced_points(u_vector)
        Drawer.draw(solver, setup)

        print(f"Done after {(time.perf_counter() - t)} s")
        print(f"All done after {(time.perf_counter() - t0)} s")
