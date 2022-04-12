"""
Created at 18.02.2021
"""
import numpy as np

from conmech.dynamics.statement import Statement, Variables


class Solver:
    def __init__(
        self,
        statement,
        mesh,
        body_prop,
        variables,
        contact_law,
        friction_bound,
    ):
        self.body_prop = body_prop
        self.contact_law = contact_law
        self.friction_bound = friction_bound

        self.mesh = mesh
        self.statement: Statement = statement

        self.current_time = 0
        self.u_vector = np.zeros(self.mesh.independent_nodes_count * 2)
        self.v_vector = np.zeros(self.mesh.independent_nodes_count * 2)
        self.t_vector = np.zeros(self.mesh.independent_nodes_count)
        self.var = variables

        self.elasticity = mesh.elasticity

        self.statement.update(
            Variables(
                displacement=self.u_vector,
                velocity=self.v_vector,
                temperature=self.t_vector,
                time_step=self.var.time_step,
            )
        )

    def __str__(self):
        raise NotImplementedError()

    def iterate(self, velocity):
        self.v_vector = velocity.reshape(-1)
        self.u_vector = self.u_vector + self.var.time_step * self.v_vector

    def solve(self, initial_guess, *, velocity: np.ndarray, **kwargs):
        raise NotImplementedError()
