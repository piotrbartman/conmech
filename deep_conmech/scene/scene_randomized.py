import numpy as np

from conmech.helpers import nph
from conmech.properties.body_properties import DynamicBodyProperties
from conmech.properties.mesh_properties import MeshProperties
from conmech.properties.obstacle_properties import ObstacleProperties
from conmech.properties.schedule import Schedule
from conmech.scene.scene import Scene


class SceneRandomized(Scene):
    def __init__(
        self,
        mesh_prop: MeshProperties,
        body_prop: DynamicBodyProperties,
        obstacle_prop: ObstacleProperties,
        schedule: Schedule,
        normalize_by_rotation: bool,
        create_in_subprocess: bool,
        with_schur: bool = True,
    ):
        super().__init__(
            mesh_prop=mesh_prop,
            body_prop=body_prop,
            obstacle_prop=obstacle_prop,
            schedule=schedule,
            normalize_by_rotation=normalize_by_rotation,
            create_in_subprocess=create_in_subprocess,
            with_schur=with_schur,
        )
        self.velocity_in_random_factor = 0
        self.displacement_in_random_factor = 0
        self.displacement_to_velocity_noise = 0
        self.velocity_randomization = np.zeros_like(self.initial_nodes)
        self.displacement_randomization = np.zeros_like(self.initial_nodes)
        # printer.print_setting_internal(self, f"output/setting_{helpers.get_timestamp()}.png", None, "png", 0)

    # def remesh(self):
    #    super().remesh()
    #    self.set_randomization(self.randomized_inputs)

    @property
    def randomized_inputs(self):
        return (
            self.velocity_in_random_factor != 0
            or self.displacement_in_random_factor != 0
            or self.displacement_to_velocity_noise != 0
        )

    def unset_randomization(self):
        self.velocity_in_random_factor = 0
        self.displacement_in_random_factor = 0
        self.displacement_to_velocity_noise = 0
        self.regenerate_randomization()

    def set_randomization(self, config):
        self.velocity_in_random_factor = config.td.velocity_in_random_factor
        self.displacement_in_random_factor = config.td.displacement_in_random_factor
        self.displacement_to_velocity_noise = config.td.displacement_to_velocity_noise
        self.regenerate_randomization()

    def regenerate_randomization(self):
        self.velocity_randomization = nph.generate_normal(
            rows=self.nodes_count,
            columns=self.dimension,
            scale=self.velocity_in_random_factor,
        )
        self.displacement_randomization = nph.generate_normal(
            rows=self.nodes_count,
            columns=self.dimension,
            scale=self.displacement_in_random_factor,
        )
        # Do not randomize boundaries
        self.displacement_randomization[self.boundary_indices] = 0.0
        self.velocity_randomization[self.boundary_indices] = 0.0

    @property
    def normalized_velocity_randomization(self):
        return self.normalize_rotate(self.velocity_randomization)

    @property
    def normalized_displacement_randomization(self):
        return self.normalize_rotate(self.displacement_randomization)

    @property
    def randomized_velocity_old(self):
        return self.velocity_old + self.velocity_randomization

    @property
    def randomized_displacement_old(self):
        return self.displacement_old + self.displacement_randomization

    @property
    def a_correction(self):
        u_correction = self.displacement_to_velocity_noise * (
            self.displacement_randomization / (self.time_step**2)
        )
        v_correction = (
            (1.0 - self.displacement_to_velocity_noise)
            * self.velocity_randomization
            / self.time_step
        )
        return -1.0 * (u_correction + v_correction)

    @property
    def normalized_a_correction(self):
        return self.normalize_rotate(self.a_correction)

    def make_dirty(self):
        self.velocity_old = self.randomized_velocity_old
        self.displacement_old = self.randomized_displacement_old

        self.unset_randomization()

    def iterate_self(self, acceleration, temperature=None):
        _ = temperature
        if self.randomized_inputs:
            self.regenerate_randomization()
        super().iterate_self(acceleration)

    @property
    def input_velocity_old(self):  # normalized_randomized_velocity_old
        return self.normalized_velocity_old + self.normalized_velocity_randomization

    @property
    def input_displacement_old(self):  # normalized_randomized_displacement_old
        return self.normalized_displacement_old + self.normalized_displacement_randomization

    @property
    def input_forces(self):
        return self.normalized_inner_forces

    @property
    def boundary_forces(self):
        return self.normalized_inner_forces[self.boundary_indices]