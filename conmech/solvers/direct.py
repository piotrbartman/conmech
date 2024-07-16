"""
Created at 18.02.2021
"""

from typing import Optional, Callable

import numpy as np
import scipy.optimize

from conmech.dynamics.statement import Statement
from conmech.scenarios.problems import ContactLaw
from conmech.scene.body_forces import BodyForces
from conmech.solvers._solvers import SolversRegistry
from conmech.solvers.solver import Solver
from conmech.solvers.solver_methods import make_equation


@SolversRegistry.register("dynamic", "direct")
@SolversRegistry.register("static", "direct")
class Direct(Solver):
    def __init__(
        self,
        statement: Statement,
        body: BodyForces,
        time_step: float,
        contact_law: Optional[ContactLaw] = None,
        friction_bound: Optional[Callable[[float], float]] = None,
        driving_vector: bool = False,
    ):
        super().__init__(
            statement,
            body,
            time_step,
            contact_law,
            friction_bound,
            driving_vector,
        )
        self.equation: Optional[Callable] = None

        if contact_law is not None:
            self.equation = make_equation(
                jn=contact_law.subderivative_normal_direction,
                jt=contact_law.regularized_subderivative_tangential_direction,
                contact=(
                    contact_law.general_contact_condition
                    if hasattr(contact_law, "general_contact_condition")
                    else None
                ),
                h_functional=friction_bound,
            )

    def __str__(self) -> str:
        return "direct"

    @property
    def node_relations(self) -> np.ndarray:
        return self.statement.left_hand_side.data

    @property
    def node_forces(self) -> np.ndarray:
        return self.statement.right_hand_side

    def _solve_impl(
            self,
            initial_guess: np.ndarray,
            *,
            variable_old: np.ndarray,
            displacement: np.ndarray,
            **kwargs
    ) -> np.ndarray:
        displacement = np.squeeze(displacement.copy().reshape(1, -1))
        if self.equation is not None:
            result = scipy.optimize.fsolve(
                self.equation,
                initial_guess,
                args=(
                    variable_old,
                    self.body.mesh.nodes,
                    self.body.mesh.contact_boundary,
                    self.body.mesh.boundaries.contact_normals,
                    self.node_relations,
                    self.node_forces,
                    displacement,
                    self.body.dynamics.acceleration_operator.SM1.data,
                    self.time_step
                ),
            )
        else:
            result = np.linalg.solve(self.node_relations, self.node_forces)
            result_len = len(result)
            var_len = len(initial_guess.ravel())
            if result_len < var_len:
                result_ = np.zeros(var_len)
                result_[:result_len] = result[:]
                result = result_
        return np.asarray(result)
