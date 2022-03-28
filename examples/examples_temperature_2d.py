import deep_conmech.scenarios as scenarios
from deep_conmech.scenarios import *
from conmech.helpers import cmh
from deep_conmech.common.plotter import plotter_mapper
from deep_conmech.scenarios import *
from deep_conmech.simulator.setting.setting_temperature import SettingTemperature
from deep_conmech.simulator.solver import Solver
import os


def main(mesh_density=4, final_time=4, plot_animation=True):
    all_scenarios = []

    polygon_scenario = lambda i, forces_function, obstacle : TemperatureScenario(
                id=f"polygon_{i}",
                mesh_data=MeshData(
                    dimension=2,
                    mesh_type=m_polygon,
                    scale=[1],
                    mesh_density=[mesh_density],
                    is_adaptive=False,
                ),
                body_prop=default_temp_body_prop,
                obstacle_prop=default_obstacle_prop,
                schedule=Schedule(final_time=final_time),
                forces_function=forces_function,
                obstacles=o_front,
                heat_function=np.array([0]),
            )

    all_scenarios.extend(
        [
            polygon_scenario(i=0, forces_function=f_rotate_fast, obstacle=o_side),
            polygon_scenario(i=1, forces_function=f_rotate_fast, obstacle=o_front),
            polygon_scenario(i=2, forces_function=np.array([1,0]), obstacle=o_front)
        ]
    )

    C_temp_body_prop = [
        default_temp_body_prop,
        get_temp_body_prop(
            C_coeff=np.array([[1.5, 0], [0, 0.5]]), K_coeff=default_K_coeff,
        ),
        get_temp_body_prop(
            C_coeff=np.array([[1.0, -0.5], [-0.5, 1.0]]), K_coeff=default_K_coeff,
        ),
        # not allowed in physical law
        get_temp_body_prop(
            C_coeff=np.array([[1.0, 0.5], [-0.5, 1.0]]), K_coeff=default_K_coeff,
        ),
    ]
    all_scenarios.extend(
        [
            TemperatureScenario(
                id=f"C_{i}",
                mesh_data=MeshData(
                    dimension=2,
                    mesh_type=m_rectangle,
                    scale=[1],
                    mesh_density=[mesh_density],
                    is_adaptive=False,
                ),
                body_prop=temp_body_prop,
                obstacle_prop=default_obstacle_prop,
                schedule=Schedule(final_time=final_time),
                forces_function=np.array([0, 0]),
                obstacles=o_side,
                heat_function=np.array([2]),
            )
            for i, temp_body_prop in enumerate(C_temp_body_prop)
        ]
    )

    K_temp_body_prop = [
        default_temp_body_prop,
        get_temp_body_prop(
            C_coeff=default_C_coeff, K_coeff=np.array([[0.5, 0], [0, 0.5]]),
        ),
        get_temp_body_prop(
            C_coeff=default_C_coeff, K_coeff=np.array([[0.1, 0.1], [0.1, 0.1]]),
        ),
        get_temp_body_prop(
            C_coeff=default_C_coeff, K_coeff=np.array([[0.1, -0.1], [-0.1, 0.1]]),
        ),
        # not allowed in physical law
        get_temp_body_prop(
            C_coeff=default_C_coeff, K_coeff=np.array([[0.1, -0.1], [0.1, 0.1]]),
        ),
    ]
    all_scenarios.extend(
        [
            TemperatureScenario(
                id=f"K_{i}",
                mesh_data=MeshData(
                    dimension=2,
                    mesh_type=m_rectangle,
                    scale=[1],
                    mesh_density=[mesh_density],
                    is_adaptive=False,
                ),
                body_prop=temp_body_prop,
                obstacle_prop=default_obstacle_prop,
                schedule=Schedule(final_time=final_time),
                forces_function=np.array([0, 0]),
                obstacles=o_side,
                heat_function=h_corner,
            )
            for i, temp_body_prop in enumerate(K_temp_body_prop)
        ]
    )


    for scenario in all_scenarios:
        print("-----")
        plotter_mapper.print_one_dynamic(
            Solver.solve_with_temperature,
            scenario,
            SettingTemperature.get_setting,
            catalog="EXAMPLES TEMPERATURE 2D",  # os.path.basename(__file__),
            simulate_dirty_data=False,
            plot_base=False,
            plot_detailed=True,
            plot_animation=plot_animation,
        )


if __name__ == "__main__":
    main()