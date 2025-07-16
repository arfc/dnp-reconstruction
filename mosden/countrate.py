from mosden.utils.input_handler import InputHandler
from mosden.utils.csv_handler import CSVHandler
import os


class CountRate:
    """
    Class to handle the delayed neutron count rate calculations
    """
    def __init__(self, input_path: str) -> None:
        self.input_path: str = input_path
        self.input_handler: InputHandler = InputHandler(input_path)
        self.input_data: dict = self.input_handler.read_input()
        self.data_dir: str = self.input_data['file_options']['processed_data_dir']
        self.output_dir: str = self.input_data['file_options']['output_dir']

        self.parent_feed: bool = self.input_data['modeling_options']['parent_feeding']
        self.num_times: int = self.input_data['modeling_options']['num_decay_times']
        self.decay_time: float = self.input_data['modeling_options']['decay_time']
        
        self.half_life_dir: str = self.input_data['data_options']['decay_constant']
        self.cross_section_dir: str = self.input_data['data_options']['cross_section']
        self.emission_probability_dir: str = self.input_data['data_options']['emission_probability']
        self.fission_yield_dir: str = self.input_data['data_options']['fission_yield']['data'] 
        return None
    
    def calculate_count_rate(self) -> None:
        """
        Calculate the delayed neutron count rate from
        concentrations using various methods
        """
        # Get data (conc.csv and emission probability processed data)
        # Calculate count rate (method (depletion linear interp parent feeding or not), decay time, decay time steps, depletion time steps)
        # Generate CSV with count rate over time
        emission_prob_path: str = os.path.join(self.data_dir, self.emission_probability_dir, 'eval.csv')
        emission_prob_data = CSVHandler(emission_prob_path, create=False)

        return