from utils.nuclide_handler import get_all_nuclides
import numpy as np
import mypy
import warnings
import logging

class Concentrations:
    """
    Class to handle the formation of concentrations from data
    """
    def __init__(self, input_data: dict) -> None:
        self.input_data = input_data
        self.all_nuclides = get_all_nuclides(input_data)
        return None
    
    def data_concentrations(self) -> dict[str, list[float]]:
        """
        Generate the concentrations of each nuclide from nuclear data based on
        irradiation of the sample for the irradiation times.

        Returns
        -------
        concs : dict[str, list[float]]
            Dictionary of concentrations for each nuclide

        """
        raise NotImplementedError("This method should be implemented")
        return {}
    
if __name__ == "__main__":
    # Example usage
    input_data = {
        "nuclide1": [1.0, 2.0, 3.0],
        "nuclide2": [4.0, 5.0, 6.0]
    }
    conc = Concentrations(input_data)
    print(conc.data_concentrations())