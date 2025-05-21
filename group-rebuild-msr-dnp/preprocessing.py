from utils.input_handler import InputHandler

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
        return None
    
    def openmc_preprocess(self) -> None:
        """
        Processes OpenMC all chain_* files and cross section data using OpenMC
        """
        return None
    
    def _openmc_chain_preprocess(self) -> None:
        """
        Processes OpenMC all chain_* files
        """
        return None