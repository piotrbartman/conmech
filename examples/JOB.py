"""
Created at 21.08.2019
"""
from dataclasses import dataclass

import numpy as np
from conmech.helpers.config import Config
from conmech.mesh.boundaries_description import BoundariesDescription
from conmech.plotting.drawer import Drawer
from conmech.scenarios.problems import ContactLaw, Quasistatic
from conmech.simulations.problem_solver import TimeDependent as QuasistaticProblemSolver


class JureczkaOchal2018(ContactLaw):
    @staticmethod
    def potential_normal_direction(u_nu: float) -> float:
        if u_nu <= 0:
            return 0.0
        if u_nu < 0.1:
            return 30 * u_nu * u_nu
        return 3

    @staticmethod
    def potential_tangential_direction(u_tau: np.ndarray) -> float:
        return -0.3 * np.exp(-np.linalg.norm(u_tau)) + 0.7 * np.linalg.norm(u_tau)

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


@dataclass()
class QuasistaticSetup(Quasistatic):
    grid_height: ... = 1.0
    elements_number: ... = (16, 16)
    mu_coef: ... = 2
    la_coef: ... = 2
    th_coef: ... = 4
    ze_coef: ... = 4
    time_step = 1/32
    contact_law: ... = JureczkaOchal2018

    @staticmethod
    def inner_forces(x):
        return np.array([-2.5, -0.5])

    @staticmethod
    def outer_forces(x):
        return np.array([0, 0])

    @staticmethod
    def friction_bound(u_nu: float) -> float:
        if u_nu < 0:
            return 0
        if u_nu < 0.1:
            return 30 * u_nu
        return 3

    boundaries: ... = BoundariesDescription(
        contact=lambda x: x[1] == 0, dirichlet=lambda x: x[0] == 0
    )


def main(show: bool = True, save: bool = False):
    setup = QuasistaticSetup(mesh_type="cross")
    runner = QuasistaticProblemSolver(setup, "schur")

    states = runner.solve(
        n_steps=32,
        output_step=(0, 32),
        verbose=True,
        initial_displacement=setup.initial_displacement,
        initial_velocity=setup.initial_velocity,
    )
    config = Config()
    for state in states:
        Drawer(state=state, config=config).draw(show=show, save=save)


if __name__ == "__main__":
    main(show=True)
