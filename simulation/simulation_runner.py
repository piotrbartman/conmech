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
                                     setup.gridHeight,
                                     left=setup.grid_left_border,
                                     top=setup.grid_top_border,
                                     right=setup.grid_right_border,
                                     bottom=setup.grid_bottom_border)
        solver = SolverFactory.construct(grid=grid, setup=setup)

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
