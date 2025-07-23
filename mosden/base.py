from mosden.utils.input_handler import InputHandler
from pathlib import Path
from mosden.utils.csv_handler import CSVHandler
import os

class BaseClass:
    def __init__(self, input_path: str) -> None:
        self.input_path: str = input_path
        self.input_handler: InputHandler = InputHandler(input_path)
        self.input_data: dict = self.input_handler.read_input()
        self.name: str = self.input_data['name']
        self.energy_MeV: float = self.input_data['data_options']['energy_MeV']
        self.fissiles: dict[str, float] = self.input_data['data_options']['fissile_fractions']
        
        self.data_types: list[str] = ['fission_yield', 'half_life', 'cross_section', 'emission_probability']

        self.processed_data_dir: str = self.input_data['file_options']['processed_data_dir']
        self.concentration_path: str = os.path.join(self.input_data['file_options']['output_dir'], 'concentrations.csv')
        self.countrate_path: str = os.path.join(self.input_data['file_options']['output_dir'], 'count_rate.csv')
        self.group_path: str = os.path.join(self.input_data['file_options']['output_dir'], 'group_parameters.csv')
        return None
    

    def _read_processed_data(self, data_type: str) -> dict[str: dict[str: float]]:
        """
        Read the processed data for a given fissile nuclide.

        Parameters
        ----------
        data_type : str
            The type of data to read (e.g., "fission_yield", "half_life", "cross_section", "emission_probability").

        Returns
        -------
        data : dict[str: dict[str: float]]
            The processed data for the fissile nuclide.
        """
        data_path = os.path.join(self.processed_data_dir, f'{data_type}.csv')
        csv_handler = CSVHandler(data_path, create=False)
        if not csv_handler._file_exists():
            raise FileNotFoundError(f"Processed data file {data_path} does not exist.") 
        data = csv_handler.read_csv()
        return data