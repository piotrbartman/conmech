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
        self.var = variables

        self.elasticity = mesh.elasticity

        self.statement.update(self.var)

    def __str__(self):
        raise NotImplementedError()

    def iterate(self, velocity):
        self.var.velocity = velocity.reshape(-1)
        self.var.displacement = self.var.displacement + self.var.time_step * self.var.velocity

    def solve(self, initial_guess, *, velocity: np.ndarray, **kwargs):
        raise NotImplementedError()
