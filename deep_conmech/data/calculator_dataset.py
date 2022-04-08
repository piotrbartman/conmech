from typing import List

from conmech.scenarios.scenarios import Scenario
from conmech.solvers.calculator import Calculator
from deep_conmech.data.scenario_dataset import ScenariosDataset
from deep_conmech.training_config import TrainingConfig


class CalculatorDataset(ScenariosDataset):
    def __init__(
            self,
            description: str,
            all_scenarios: List[Scenario],
            load_to_ram: bool,
            config: TrainingConfig,
    ):
        super().__init__(
            description=f"{description}_calculator", all_scenarios=all_scenarios,
            solve_function=Calculator.solve_all,
            load_to_ram=load_to_ram, config=config
        )