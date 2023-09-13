import pickle
from dataclasses import dataclass
from typing import Optional, Type

import numpy as np
from conmech.helpers.config import Config
from conmech.mesh.boundaries_description import BoundariesDescription
from conmech.plotting.drawer import Drawer
from conmech.scenarios.problems import PoissonProblem, ContactLaw, DummyContactLaw
from conmech.simulations.problem_solver import PoissonSolver
from examples.Bartman_Ochal_Tarzia_2023 import make_slope_contact_law, simulate, draw
from examples.error_estimates_tarzia import compare


@dataclass()
class StaticPoissonSetup(PoissonProblem):
    grid_height: ... = 1
    elements_number: ... = (4, 8)

    contact_law: ... = make_slope_contact_law(alpha=1000, b=8, example="13")

    @staticmethod
    def internal_temperature(x: np.ndarray, t: Optional[float] = None) -> np.ndarray:
        return np.array([0.0])

    @staticmethod
    def outer_temperature(x: np.ndarray, t: Optional[float] = None) -> np.ndarray:
        if x[0] == 2.0:
            return np.array([-1])
        return np.array([0])

    boundaries: ... = BoundariesDescription(
        contact=lambda x: x[0] == 0.0,
    )


def main(config: Config):
    """
    Entrypoint to example.

    To see result of simulation you need to call from python `main(Config().init())`.
    """
    setup = StaticPoissonSetup(mesh_type="cross")
    alphas = [1e1, 1e3, np.inf]

    for alpha in alphas:
        try:
            if config.force:
                simulate(config, alpha, None, _set_alpha, setup)
            draw(config, alpha, None)
        except FileNotFoundError:
            simulate(config, alpha, None, _set_alpha, setup)
            draw(config, alpha, None)


def _set_alpha(setup, alpha):
    setup.contact_law = make_slope_contact_law(alpha=alpha, b=8, example="13")
    if alpha == np.inf:
        setup.boundaries = BoundariesDescription(
            dirichlet=(
                lambda x: x[0] == 0.0,
                lambda x: np.full(x.shape[0], 8),
            ),
        )


if __name__ == "__main__":
    main(Config(outputs_path="./output/Tarzia", force=False).init())

