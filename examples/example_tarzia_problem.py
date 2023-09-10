import pickle
from dataclasses import dataclass
from typing import Optional, Type

import numpy as np
from conmech.helpers.config import Config
from conmech.mesh.boundaries_description import BoundariesDescription
from conmech.plotting.drawer import Drawer
from conmech.scenarios.problems import PoissonProblem, ContactLaw
from conmech.simulations.problem_solver import PoissonSolver
from examples.error_estimates_tarzia import compare


def make_slope_contact_law(slope: float) -> Type[ContactLaw]:
    class TarziaContactLaw(ContactLaw):
        @staticmethod
        def potential_normal_direction(u_nu: float) -> float:
            b = 5
            r = u_nu
            # EXAMPLE 11
            # if r < b:
            #     result = (r - b) ** 2
            # else:
            #     result = 1 - np.exp(-(r-b))
            # EXAMPLE 13
            result = 0.5 * (r - b) ** 2
            result *= slope
            return result

        @staticmethod
        def subderivative_normal_direction(u_nu: float, v_nu: float) -> float:
            raise NotImplementedError()

        @staticmethod
        def regularized_subderivative_tangential_direction(
            u_tau: np.ndarray, v_tau: np.ndarray, rho=1e-7
        ) -> float:
            """
            Coulomb regularization
            """
            raise NotImplementedError()

    return TarziaContactLaw


@dataclass()
class StaticPoissonSetup(PoissonProblem):
    grid_height: ... = 1
    elements_number: ... = (6, 12)

    contact_law: ... = make_slope_contact_law(slope=1000)

    @staticmethod
    def internal_temperature(x: np.ndarray, t: Optional[float] = None) -> np.ndarray:
        if 0.4 <= x[0] <= 0.6 and 0.4 <= x[1] <= 0.6:
            return np.array([-10.0])
        return np.array([2 * np.pi**2 * np.sin(np.pi * x[0]) * np.sin(np.pi * x[1])])

    @staticmethod
    def outer_temperature(x: np.ndarray, t: Optional[float] = None) -> np.ndarray:
        if x[1] > 0.5:
            return np.array([1.0])
        return np.array([-1.0])

    boundaries: ... = BoundariesDescription(
        dirichlet=(
            lambda x: x[0] == 0.0,
            lambda x: np.full(x.shape[0], 5),
        ),
        contact=lambda x: x[0] == 2.0,
    )


def main(config: Config):
    """
    Entrypoint to example.

    To see result of simulation you need to call from python `main(Config().init())`.
    """
    alphas = [0.01, 0.1, 1, 10, 100, 1000, 10000]
    ihs = [4, 8, 16, 32, 64, 128, 256]
    alphas = alphas[:]
    ihs = ihs[:]  # TODO

    for ih in ihs:
        for alpha in alphas:
            try:
                if config.force:
                    simulate(config, alpha, ih)
                draw(config, alpha, ih)
            except FileNotFoundError:
                simulate(config, alpha, ih)
                draw(config, alpha, ih)
    convergence(config, alphas, ihs)


def simulate(config, alpha, ih):
    print(f"Simulate {alpha=}, {ih=}")
    setup = StaticPoissonSetup(mesh_type="cross")
    setup.contact_law = make_slope_contact_law(slope=alpha)
    setup.elements_number = (1 * ih, 2 * ih)

    runner = PoissonSolver(setup, "schur")

    state = runner.solve(verbose=True)

    if config.outputs_path:
        with open(
            f"{config.outputs_path}/alpha_{alpha}_ih_{ih}",
            "wb+",
        ) as output:
            # Workaround
            state.body.dynamics = None
            state.body.properties.relaxation = None
            state.setup = None
            state.constitutive_law = None
            pickle.dump(state, output)


def draw(config, alpha, ih):
    with open(f"{config.outputs_path}/alpha_{alpha}_ih_{ih}", "rb") as output:
        state = pickle.load(output)
    max_ = max(max(state.temperature), 1)
    min_ = min(min(state.temperature), 0)
    drawer = Drawer(state=state, config=config)
    drawer.cmap = "plasma"
    drawer.field_name = "temperature"
    # drawer.draw(
    #     show=config.show, save=config.save, foundation=False, field_max=max_, field_min=min_
    # )


def convergence(config, alphas, ihs):
    cvgs = {
        "hn_ac": [(ihs[0], a) for a in alphas],
        "hx_ac": [(ihs[-1], a) for a in alphas],
        "hc_an": [(ih, alphas[0]) for ih in ihs],
        "hc_ax": [(ih, alphas[-1]) for ih in ihs],
        "hx_ax": [(ihs[i], alphas[i]) for i in range(len(ihs))]
    }

    for cvg, params in cvgs.items():
        ih = params[-1][0]
        alpha = params[-1][1]
        with open(f"{config.outputs_path}/alpha_{alpha}_ih_{ih}", "rb") as output:
            ref = pickle.load(output)
        for ih, alpha in params[:-1]:
            with open(f"{config.outputs_path}/alpha_{alpha}_ih_{ih}", "rb") as output:
                sol = pickle.load(output)
            print(cvg, compare(ref, sol))


if __name__ == "__main__":
    main(Config(outputs_path="./output/BOT2023", force=False).init())

    "cd ~/devel/conmech && git pull; PYTHONPATH=/home/prb/devel/conmech venv/bin/python3.11 examples/example_tarzia_problem.py &"
    # hn_ac(0, 1675.6377901503872)
    # hn_ac(0, 1424.4089747346386)
    # hn_ac(0, 570.3035428534573)
    # hn_ac(0, 82.4145654649069)
    # hn_ac(0, 9.024423274856195)
    # hn_ac(0, 0.8570601649719317)

    # hx_ac(0, 396671.37946907803)
    # hx_ac(0, 362252.1828834345)
    # hx_ac(0, 196957.26862225882)
    # hx_ac(0, 41166.897184751695)
    # hx_ac(0, 6405.289959436339)
    # hx_ac(0, 969.8479417923099)

    # hc_an(0, 34810.85013895854)
    # hc_an(0, 22577.794240169485)
    # hc_an(0, 5505.7223168591445)
    # hc_an(0, 6102.412144549319)

    # hc_ax(0, 70704.1114522041)
    # hc_ax(0, 70893.28055019652)
    # hc_ax(0, 63141.28381364561)
    # hc_ax(0, 6887.845735312516)

    # hx_ax(0, 369667.3839689251)
    # hx_ax(0, 305495.01151754276)
    # hx_ax(0, 98587.20060468125)
    # hx_ax(0, 42172.21447780231)
