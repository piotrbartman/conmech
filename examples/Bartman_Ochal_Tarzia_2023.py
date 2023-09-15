# CONMECH @ Jagiellonian University in Kraków
#
# Copyright (C) 2023  Piotr Bartman <piotr.bartman@uj.edu.pl>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
# USA.
import pickle
from dataclasses import dataclass
from typing import Optional, Type

import matplotlib.pyplot as plt
import numba
import numpy as np
from conmech.helpers.config import Config
from conmech.mesh.boundaries_description import BoundariesDescription
from conmech.plotting.drawer import Drawer
from conmech.scenarios.problems import PoissonProblem, ContactLaw, DummyContactLaw
from conmech.simulations.problem_solver import PoissonSolver
from examples.error_estimates_tarzia import compare


def make_slope_contact_law(alpha: float, b: float, example: str = "11") -> Type[ContactLaw]:
    """
    Potentials based on https://doi.org/10.48550/arXiv.2106.04702
    """
    if alpha == np.inf:
        return DummyContactLaw

    @numba.njit()
    def example_11(r):
        """
        EXAMPLE 11.
        """
        if r < b:
            result = (r - b) ** 2
        else:
            result = 1 - np.exp(-(r-b))
        return result

    @numba.njit()
    def example_13(r):
        """
        EXAMPLE 13.
        """
        result = 0.5 * (r - b) ** 2
        return result

    examples = {"11": example_11, "13": example_13}
    example_func = examples[example]

    class TarziaContactLaw(ContactLaw):
        @staticmethod
        def potential_normal_direction(u_nu: float) -> float:
            result = alpha * example_func(u_nu)
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

    contact_law: ... = None

    @staticmethod
    def internal_temperature(x: np.ndarray, t: Optional[float] = None) -> np.ndarray:
        return np.array([0.0])

    @staticmethod
    def outer_temperature(x: np.ndarray, t: Optional[float] = None) -> np.ndarray:
        if x[0] == 2.0:
            return np.array([ - 1 ])
        return np.array([ + 1 ])

    boundaries: ... = BoundariesDescription(
        dirichlet=(
            lambda x: x[1] == 1.0,
            lambda x: np.full(x.shape[0], 0),
        ),
        contact=lambda x: x[1] == 0.0
    )


def main(config: Config):
    """
    Entrypoint to example.

    To see result of simulation you need to call from python `main(Config().init())`.
    """
    alphas = [1e-2, 1, 1e1, 1e3, 1e4, np.inf]
    ihs = [4, 8, 16, 32, 64, 128, 256]
    alphas = alphas[:-1]
    ihs = ihs[:5]  # TODO

    for ih in ihs:
        for alpha in alphas:
            try:
                if config.force:
                    simulate(config, alpha, ih)
                draw(config, alpha, ih)
            except FileNotFoundError:
                simulate(config, alpha, ih)
                draw(config, alpha, ih)
    try:
        if config.force:
            convergence(config, alphas, ihs)
        draw_convergence(config, alphas, ihs)
    except (FileNotFoundError, KeyError):
        convergence(config, alphas, ihs)
        draw_convergence(config, alphas, ihs)


def _set_alpha(setup, alpha):
    setup.contact_law = make_slope_contact_law(alpha=alpha, b=0, example="11")
    if alpha == np.inf:
        setup.boundaries = BoundariesDescription(
            dirichlet=(
                lambda x: x[1] == 1.0 or x[1] == 0.0,
                lambda x: np.full(x.shape[0], 0),
            ),
        )


def simulate(config, alpha, ih=None, alpha_setter=_set_alpha, setup=None):
    print(f"Simulate {alpha=}, {ih=}")
    setup = setup or StaticPoissonSetup(mesh_type="cross")
    alpha_setter(setup, alpha)
    if ih is not None:
        setup.elements_number = (1 * ih, 2 * ih)

    runner = PoissonSolver(setup, solving_method="auto")

    state = runner.solve(verbose=True, disp=True, method="BFGS")

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
    if not config.show:
        return
    max_ = max(max(state.temperature), 0.5)
    min_ = min(min(state.temperature), -0.5)
    drawer = Drawer(state=state, config=config)
    drawer.cmap = "plasma"
    drawer.field_name = "temperature"
    drawer.original_mesh_color = None
    drawer.deformed_mesh_color = None
    drawer.draw(
        show=config.show, save=config.save, foundation=False, field_max=max_, field_min=min_
    )


def convergence(config, alphas, ihs):
    cvgs = {
        "hn_ac": {(ihs[0], a): None for a in alphas},
        "hx_ac": {(ihs[-1], a): None for a in alphas},
        "hc_an": {(ih, alphas[0]): None for ih in ihs},
        "hc_ax": {(ih, alphas[-1]): None for ih in ihs},
        "hc_ac": {(ihs[i], alphas[0]): None for i in range(len(ihs))}  # TODO alphas[0]
    }

    for cvg, map_ in cvgs.items():
        params = list(map_.keys())
        ih = params[-1][0]
        alpha = params[-1][1]
        with open(f"{config.outputs_path}/alpha_{alpha}_ih_{ih}", "rb") as output:
            ref = pickle.load(output)
        for ih, alpha in params[:-1]:
            with open(f"{config.outputs_path}/alpha_{alpha}_ih_{ih}", "rb") as output:
                sol = pickle.load(output)
            cvgs[cvg][(ih, alpha)] = compare(ref, sol)
    if config.outputs_path:
        with open(
            f"{config.outputs_path}/convergences",
            "wb+",
        ) as output:
            pickle.dump(cvgs, output)

def plot_helper(X, Y, title, xlabel):
    plt.plot(X, Y)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(r"error")
    plt.grid()
    plt.semilogx()
    plt.show()

def draw_convergence(config, alphas, ihs):
    with open(f"{config.outputs_path}/convergences", "rb") as file:
        cvgs = pickle.load(file)
    if not config.show:
        return
    convs = {
        "L4": [cvgs["hn_ac"][(ihs[0], alpha)] for alpha in alphas],
        "L2": [cvgs["hx_ac"][(ihs[-1], alpha)] for alpha in alphas],
        "L3": [cvgs["hc_an"][(ih, alphas[0])] for ih in ihs],
        "L1": [cvgs["hc_ax"][(ih, alphas[-1])] for ih in ihs],
        "L5": [cvgs["hc_ac"][(ihs[i], alphas[0])] for i in range(len(ihs))]  # TODO alphas[0]
    }

    plot_helper(ihs[:-1], convs["L1"][:-1], "L1", r"$h$")
    plot_helper(alphas[:-1], convs["L2"][:-1], "L2", r"$\alpha$")
    plot_helper(ihs[:-1], convs["L3"][:-1], "L3", r"$h$")
    plot_helper(alphas[:-1], convs["L4"][:-1], "L4", r"$\alpha$")
    plot_helper(ihs[:-1], convs["L5"][:-1], "L5", r"$(h, \alpha)$")


def example_11_plot():
    b = 0
    X = np.linspace(b - 1, b + 1, num=1000)

    def dj(r):
        if r < b:
            return 2 * (r - b)
        else:
            return np.exp(-(r - b))

    Y = np.empty_like(X)
    for i, x in enumerate(X):
        Y[i] = dj(x)
    plt.plot(X, Y)
    plt.plot(X, X * 0, color="black")
    plt.title("EXAMPLE 11: $b=0$")
    plt.xlabel("$r$")
    plt.ylabel(r"$\partial j(r)$")
    plt.grid()
    # plt.semilogx()
    plt.show()


if __name__ == "__main__":
    config = Config(outputs_path="output/BOT2023", force=True, show=False).init()
    main(config)

    if config.show:
        example_11_plot()

    "PYTHONPATH=/home/prb/devel/conmech venv/bin/python3.11 examples/Bartman_Ochal_Tarzia_2023.py &"
