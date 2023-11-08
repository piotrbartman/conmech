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

from dataclasses import dataclass

import numpy as np
from conmech.helpers.config import Config
from conmech.mesh.boundaries_description import BoundariesDescription
from conmech.plotting.drawer import Drawer
from conmech.scenarios.problems import ContactLaw, StaticDisplacementProblem
from conmech.simulations.problem_solver import StaticSolver as StaticProblemSolver
from conmech.properties.mesh_description import (
    CrossMeshDescription,
    BOST2023MeshDescription,
)


GAP = 0.05


class BOST23(ContactLaw):
    @staticmethod
    def potential_normal_direction(u_nu: float) -> float:
        u_nu -= GAP
        if u_nu <= 0:
            return 0.0
        if u_nu < 0.1:
            return 10 * u_nu * u_nu
        return 0.1

    @staticmethod
    def potential_tangential_direction(u_tau: np.ndarray) -> float:
        return np.log(np.sum(u_tau * u_tau) ** 0.5 + 1)

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
        if u_nu < 0.1:
            return 8 * u_nu
        return 0.8

    boundaries: ... = BoundariesDescription(
        contact=lambda x: x[1] == 0, dirichlet=lambda x: x[0] == 0
    )


def main(config: Config):
    """
    Entrypoint to example.

    To see result of simulation you need to call from python `main(Config().init())`.
    """
    mesh_descr = BOST2023MeshDescription(
        initial_position=None,
        max_element_perimeter=1 / 8,
    )
    setup = StaticSetup(mesh_descr)
    if config.test:
        setup.elements_number = (2, 4)
    runner = StaticProblemSolver(setup, "schur")

    state = runner.solve(
        verbose=True,
        fixed_point_abs_tol=0.001,
        initial_displacement=setup.initial_displacement,
    )
    state.displaced_nodes[:, 1] += GAP
    state.body.mesh.nodes[:, 1] += GAP
    drawer = Drawer(state=state, config=config)
    drawer.colorful = False
    drawer.outer_forces_scale = 1
    drawer.draw(show=config.show, save=config.save)


if __name__ == "__main__":
    main(Config().init())
