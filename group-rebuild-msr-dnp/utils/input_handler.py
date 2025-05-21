import json

class InputHandler:
    def __init__(self, input_path: str) -> None:
        """
        This class handles the input files for all processing.

        Parameters
        ----------
        input_path : str
            Path to the input file.
        """
        self.input_path = input_path
        return None
    
    def read_input(self) -> dict:
        """
        Read the input file and return the data as a dictionary.

        Returns
        -------
        output : dict
            Dictionary containing settings and data selections.
        """
        with open(self.input_path, 'r') as file:
            output = json.load(file)
        return output
