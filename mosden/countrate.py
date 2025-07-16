from mosden.utils.input_handler import InputHandler
from mosden.utils.csv_handler import CSVHandler
from mosden.base import BaseClass
import os


class CountRate(BaseClass):
    """
    Class to handle the delayed neutron count rate calculations
    """
    def __init__(self, input_path: str) -> None:
        super().__init__(input_path)
        self.data_dir: str = self.input_data['file_options']['processed_data_dir']
        self.output_dir: str = self.input_data['file_options']['output_dir']

        self.parent_feed: bool = self.input_data['modeling_options']['parent_feeding']
        self.num_times: int = self.input_data['modeling_options']['num_decay_times']
        self.decay_time: float = self.input_data['modeling_options']['decay_time']
        return None
    
    def calculate_count_rate(self) -> None:
        """
        Calculate the delayed neutron count rate from
        concentrations using various methods
        """
        # Get data (conc.csv and emission probability processed data)
        # Calculate count rate (method (depletion linear interp parent feeding or not), decay time, decay time steps, depletion time steps)
        # Generate CSV with count rate over time
        emission_prob_path: str = self.processed_data_paths['emission_probability']
        emission_prob_data = CSVHandler(emission_prob_path, create=False)

        return