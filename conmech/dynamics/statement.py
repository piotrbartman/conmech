from dataclasses import dataclass
from typing import Optional

import numpy as np


class Variables:
    def __init__(self, *, displacement=None, velocity=None, temperature=None, electric_potential=None, time_step=0):
        self.displacement: Optional[np.ndarray] = displacement
        self.velocity: Optional[np.ndarray] = velocity
        self.temperature: Optional[np.ndarray] = temperature
        self.electric_potential: Optional[np.ndarray] = electric_potential
        self.time_step: float = time_step


class Statement:
    def __init__(self, dynamics, dimension):
        self.dynamics = dynamics
        self.dimension = dimension
        self.left_hand_side = None
        self.right_hand_side = None

    def update_left_hand_side(self, var: Variables):
        raise NotImplementedError()

    def update_right_hand_side(self, var: Variables):
        raise NotImplementedError()

    def update(self, var: Variables):
        self.update_left_hand_side(var)
        self.update_right_hand_side(var)


class StaticDisplacementStatement(Statement):
    def __init__(self, dynamics):
        super().__init__(dynamics, 2)

    def update_left_hand_side(self, var: Variables):
        self.left_hand_side = self.dynamics.elasticity

    def update_right_hand_side(self, var: Variables):
        self.right_hand_side = self.dynamics.forces.forces_vector


class QuasistaticVelocityStatement(Statement):
    def __init__(self, dynamics):
        super().__init__(dynamics, 2)

    def update_left_hand_side(self, var: Variables):
        self.left_hand_side = self.dynamics.viscosity

    def update_right_hand_side(self, var: Variables):
        assert var.displacement is not None

        self.right_hand_side = (
            self.dynamics.forces.forces_vector - self.dynamics.elasticity @ var.displacement.T
        )


class DynamicVelocityStatement(Statement):
    def __init__(self, dynamics):
        super().__init__(dynamics, 2)

    def update_left_hand_side(self, var):
        assert var.time_step

        self.left_hand_side = (
            self.dynamics.viscosity + (1 / var.time_step) * self.dynamics.acceleration_operator
        )

    def update_right_hand_side(self, var):
        assert var.displacement is not None
        assert var.velocity is not None
        assert var.time_step

        A = -1 * self.dynamics.elasticity @ var.displacement

        A += (1 / var.time_step) * self.dynamics.acceleration_operator @ var.velocity

        self.right_hand_side = self.dynamics.forces.forces_vector + A


class DynamicVelocityWithTemperatureStatement(DynamicVelocityStatement):
    def update_right_hand_side(self, var):
        super().update_right_hand_side(var)

        assert var.temperature is not None

        A = self.dynamics.thermal_expansion.T @ var.temperature

        self.right_hand_side += A


class TemperatureStatement(Statement):
    def __init__(self, dynamics):
        super().__init__(dynamics, 1)

    def update_left_hand_side(self, var):
        assert var.time_step

        ind = self.dynamics.independent_nodes_count

        self.left_hand_side = (1 / var.time_step) * self.dynamics.acceleration_operator[
            :ind, :ind
        ] + self.dynamics.thermal_conductivity[:ind, :ind]

    def update_right_hand_side(self, var):
        assert var.velocity is not None
        assert var.temperature is not None
        assert var.time_step

        rhs = (-1) * self.dynamics.thermal_expansion @ var.velocity

        ind = self.dynamics.independent_nodes_count

        rhs += (
            (1 / var.time_step) * self.dynamics.acceleration_operator[:ind, :ind] @ var.temperature
        )
        self.right_hand_side = rhs
        # self.right_hand_side = self.inner_temperature.F[:, 0] + Q1 - C2Xv - C2Yv  # TODO #50


class PiezoelectricStatement(Statement):
    def __init__(self, dynamics):
        super().__init__(dynamics, 1)

    def update_left_hand_side(self, var):
        ind = self.dynamics.independent_nodes_count

        self.left_hand_side = self.dynamics.permittivity[:ind, :ind]

    def update_right_hand_side(self, var):
        assert var.velocity is not None
        assert var.time_step
        assert var.temperature is not None

        rhs = (-1) * self.dynamics.thermal_expansion @ var.velocity

        ind = self.dynamics.independent_nodes_count

        rhs += (
            (1 / var.time_step) * self.dynamics.acceleration_operator[:ind, :ind] @ var.temperature
        )
        self.right_hand_side = rhs
        # self.right_hand_side = self.inner_temperature.F[:, 0] + Q1 - C2Xv - C2Yv  # TODO #50
