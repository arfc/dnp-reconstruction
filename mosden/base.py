from mosden.utils.input_handler import InputHandler
from pathlib import Path
from mosden.utils.csv_handler import CSVHandler
import os
import logging
import json
from typing import Any
from time import time

class BaseClass:
    def __init__(self, input_path: str) -> None:
        self.omc_data_words: list[str] = ['omcchain']
        self.endf_data_words: list[str] = ['nfy']
        self.iaea_data_words: list[str] = ['iaea']

        self.input_path: str = input_path
        self.input_handler: InputHandler = InputHandler(input_path)
        self.input_data: dict = self.input_handler.read_input()

        self.log_file: str = self.input_data.get('file_options', {}).get('log_file', 'log.log')
        self.log_level: int = self.input_data.get('file_options', {}).get('log_level', 1)
        logger_overwrite: bool = self.input_data.get('file_options', {}).get('overwrite', {}).get('logger', False)

        self.logger: logging.Logger = logging.getLogger(__name__)
        if logger_overwrite:
            log_mode = 'w'
        else:
            log_mode = 'a'
        logging.basicConfig(filename=self.log_file,
                            level=self.log_level,
                            filemode=log_mode)
        
        self.name: str = self.input_data['name']
        self.energy_MeV: float = self.input_data['data_options']['energy_MeV']
        self.fissiles: dict[str, float] = self.input_data['data_options']['fissile_fractions']
        
        self.data_types: list[str] = ['fission_yield', 'half_life', 'cross_section', 'emission_probability']

        self.processed_data_dir: str = self.input_data['file_options']['processed_data_dir']
        self.concentration_path: str = os.path.join(self.input_data['file_options']['output_dir'], 'concentrations.csv')
        self.countrate_path: str = os.path.join(self.input_data['file_options']['output_dir'], 'count_rate.csv')
        self.group_path: str = os.path.join(self.input_data['file_options']['output_dir'], 'group_parameters.csv')
        self.postproc_path: str = os.path.join(self.input_data['file_options']['output_dir'], 'postproc.json')

        self.post_data: dict[str: float|str|list] = dict()
        if Path(self.postproc_path).exists():
            with open(self.postproc_path, 'r') as f:
                self.post_data = json.load(f)
        else:
            self.post_data = dict()
        
        self.names: dict[str: str] = {
            'countsMC': 'countsMC',
            'groupfitMC': 'groupfitMC'
        }
        return None
    
    def time_track(self, starttime: float, modulename: str ='') -> None:
        self.logger.info(f'{modulename} took {round(time()-starttime, 3)}s')
        return None

    
    def save_postproc(self) -> None:
        if Path(self.postproc_path).exists():
            with open(self.postproc_path, 'r') as f:
                existing_data = json.load(f)
            existing_data.update(self.post_data)
            data_to_write = existing_data
        else:
            data_to_write = self.post_data
        with open(self.postproc_path, 'w') as f:
            json.dump(data_to_write, f, indent=4)
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