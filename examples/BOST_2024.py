# CONMECH @ Jagiellonian University in Kraków
#
# Copyright (C) 2023  Piotr Bartman-Szwarc <piotr.bartman@uj.edu.pl>
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
from typing import Iterable

import numpy as np
from conmech.helpers.config import Config
from conmech.mesh.boundaries_description import BoundariesDescription
from conmech.plotting.drawer import Drawer
from conmech.scenarios.problems import ContactLaw, StaticDisplacementProblem
from conmech.simulations.problem_solver import StaticSolver as StaticProblemSolver
from conmech.properties.mesh_description import (
    BOST2023MeshDescription,
)


GAP = 0.05


class BOST23(ContactLaw):
    @staticmethod
    def potential_normal_direction(u_nu: float) -> float:
        u_nu -= GAP
        # EXAMPLE 10
        a = 0.1
        b = 0.1
        if u_nu <= 0:
            return 0.0
        if u_nu < b:
            return (a + np.exp(-b)) / (2 * b) * u_nu ** 2
        return a * u_nu - np.exp(- u_nu) + ((b + 2) * np.exp(-b) - a * b) / 2

    @staticmethod
    def potential_tangential_direction(u_tau: np.ndarray) -> float:
        return np.sum(u_tau * u_tau) ** 0.5

    @staticmethod
    def subderivative_normal_direction(u_nu: float, v_nu: float) -> float:
        return 0

    @staticmethod
    def regularized_subderivative_tangential_direction(
        u_tau: np.ndarray, v_tau: np.ndarray, rho=1e-7
    ) -> float:
        """
        Coulomb regularization
        """
        return 0


@dataclass
class StaticSetup(StaticDisplacementProblem):
    mu_coef: ... = 40
    la_coef: ... = 4
    contact_law: ... = BOST23

    @staticmethod
    def inner_forces(x, t=None):
        return np.array([-1.2, -0.8])

    @staticmethod
    def outer_forces(x, t=None):
        return np.array([0, 0])

    @staticmethod
    def friction_bound(u_nu: float) -> float:
        if u_nu < 0:
            return 0
        return u_nu

    boundaries: ... = BoundariesDescription(
        contact=lambda x: x[1] <= 0.1, dirichlet=lambda x: x[0] == 0
    )


def main(config: Config, igs: Iterable[int]):
    """
    Entrypoint to example.

    To see result of simulation you need to call from python `main(Config().init())`.
    """
    simulate = config.force
    try:
        for ig in igs:
            with open(f"{config.outputs_path}/BOST_ig_{ig}", "rb") as output:
                _ = pickle.load(output)
    except IOError:
        simulate = True

    if simulate:
        print("Simulating...")
        for ig in igs:
            mesh_descr = BOST2023MeshDescription(
                initial_position=None,
                max_element_perimeter=1 / 32,
            )
            setup = StaticSetup(mesh_descr)

            def potential_normal_direction(u_nu: float) -> float:
                u_nu -= GAP
                # EXAMPLE 10
                a = 0.1
                b = 0.1
                if u_nu <= 0:
                    result = 0.0
                elif u_nu < b:
                    result = (a + np.exp(-b)) / (2 * b) * u_nu ** 2
                else:
                    result = a * u_nu - np.exp(- u_nu) + ((b + 2) * np.exp(-b) - a * b) / 2
                return ig * result

            setup.contact_law.potential_normal_direction = potential_normal_direction
            if config.test:
                setup.elements_number = (2, 4)
            runner = StaticProblemSolver(setup, "schur")

            state = runner.solve(
                verbose=True,
                fixed_point_abs_tol=0.001,
                initial_displacement=setup.initial_displacement,
            )
            with open(f"{config.outputs_path}/BOST_ig_{ig}", "wb+") as output:
                state.body.dynamics.force.outer.source = None
                state.body.dynamics.force.inner.source = None
                state.body.properties.relaxation = None
                state.setup = None
                state.constitutive_law = None
                pickle.dump(state, output)

    print(f"Plotting {igs=}")
    for ig in igs:
        with open(f"{config.outputs_path}/BOST_ig_{ig}", "rb") as output:
            state = pickle.load(output)

        print(f"{ig=}")

        state.displaced_nodes[:, 1] += GAP
        state.body.mesh.nodes[:, 1] += GAP
        drawer = Drawer(state=state, config=config)
        drawer.colorful = False
        drawer.outer_forces_scale = 1
        drawer.draw(show=config.show, save=config.save)
    print("Done")



if __name__ == "__main__":
    show = True
    main(
        Config(save=not show, show=show, force=False).init(),
        [0] + [2**i for i in range(32)] + [2**64]
    )
