from utils.input_handler import InputHandler
from utils.csv_handler import CSVHandler
import os

class Preprocess:
    def __init__(self, input_path: str) -> None:
        """
        This class handles preprocessing of files into a common format.

        Parameters
        ----------
        input_path : str
            Path to the input file.
        """
        self.input_path = input_path
        self.input_handler = InputHandler(input_path)
        self.input_data = self.input_handler.read_input()

        self.data_dir = self.input_data['data_directory']
        return None
    
    def openmc_preprocess(self) -> None:
        """
        Processes OpenMC all chain_* files and cross section data using OpenMC
        """
        self._openmc_chain_preprocess()
        return None
    
    def _openmc_chain_preprocess(self) -> None:
        """
        Processes OpenMC all chain_* files
        """
        for file in os.listdir(self.data_dir + '/chain/'):
            file_data = self._process_chain_file(file)
            csv_path = self.data_dir + '/processed/'
            CSVHandler(csv_path).write_csv(file_data)
        return None
    
    def _process_chain_file(self, file: str) -> dict[str, dict[str, float]]:
        """
        Processes a single OpenMC chain file and returns the data as a dictionary.

        Parameters
        ----------
        file : str
            Name of the chain file to process.

        Returns
        -------
        data : dict[str, dict[str, float]]
            Dictionary containing the processed data.
        """
        # Placeholder for actual processing logic
        data = {}
        # Implement the logic to read and process the chain file here
        return data