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
        self.overwrite = self.input_data['overwrite']
        self.out_dir = self.input_data['output_directory']
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
            full_path = os.path.join(self.data_dir + '/chain/', file)
            file_data = self._process_chain_file(full_path)
            csv_path = self.data_dir + self.out_dir + file.split('.')[0] + '.csv'
            CSVHandler(csv_path, self.overwrite).write_csv(file_data)
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
        import openmc.deplete
        chain = openmc.deplete.Chain.from_xml(file)
        nuclides: list[openmc.deplete.Nuclide] = chain.nuclides
        for nuc in nuclides:
            data[nuc.name] = {}
            data[nuc.name]['half_life'] = nuc.half_life
            input(nuc.yield_data)
            input(nuc.yield_energies)
            help(nuc)
            input('...')
        fission_yields: list[dict[str: {str: float}]] = chain.fission_yields
        # Placeholder for actual processing logic
        data = {}
        # Implement the logic to read and process the chain file here
        return data


preproc = Preprocess('../examples/keepin_1957/pre_input.json')
preproc.openmc_preprocess()