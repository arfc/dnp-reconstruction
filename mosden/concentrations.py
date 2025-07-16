import numpy as np
from mosden.utils.input_handler import InputHandler
from mosden.utils.csv_handler import CSVHandler
from pathlib import Path
from uncertainties import ufloat
from mosden.base import BaseClass

class Concentrations(BaseClass):
    """
    Class to handle the formation of concentrations from data
    """
    def __init__(self, input_path: str) -> None:
        super().__init__(input_path)
        self.output_dir: str = self.input_data['file_options']['output_dir']

        self.model_method: str = self.input_data['modeling_options']['concentration_handling']
        self.fissiles: dict[str, float] = self.input_data['data_options']['fissile_fractions']
        self.data_path: str = self.input_data['file_options']['processed_data_dir']
        self.overwrite: bool = self.input_data['file_options']['overwrite']

        self.reprocessing: dict[str: float] = self.input_data['modeling_options']['reprocessing']
        self.reprocess: bool = (sum(self.reprocessing.values()) > 0)
        self.reprocess_locations: list[str] = self.input_data['modeling_options']['reprocessing_locations']
        self.t_in: float = self.input_data['modeling_options']['incore_s']
        self.t_ex: float = self.input_data['modeling_options']['excore_s']
        self.t_net: float = self.input_data['modeling_options']['net_irrad_s']
        self.irrad_type: str = self.input_data['modeling_options']['irrad_type']

        return None
    
    def generate_concentrations(self) -> None:
        """
        Generate the concentrations of each nuclide from based on
        irradiation of the sample for the irradiation times.
        """
        data: dict[str: dict[str: float]] = dict()
        if self.model_method == 'CFY':
            if self.t_ex > 0.0:
                raise NotImplementedError('Excore residence not available for CFY')
            if self.reprocess:
                raise NotImplementedError('Reprocessing not available for CFY')
            data = self.CFY_concentrations()
        else:
            raise NotImplementedError(
                f"Concentration handling method '{self.model_method}' is not implemented"
            )

        CSVHandler(self.concentration_path, self.overwrite).write_csv(data)
        return

    def CFY_concentrations(self) -> None:
        """
        Generate the concentrations of each nuclide using the CFY method.
        """
        concentrations: dict[str: dict[str: ufloat]] = dict()
        all_nucs: set[str] = set()
        for fissile, fraction in self.fissiles.items():
            CFY_data = self._read_processed_data(fissile, 'fission_yield')
            for nuclide in CFY_data.keys():
                concs = ufloat(CFY_data[nuclide]['CFY'], CFY_data[nuclide]['sigma CFY'])
                concentrations[nuclide] = concentrations.get(nuclide, 0) + concs * fraction
                all_nucs.add(nuclide)


        data: dict[str: dict[str: float]] = dict()
        for nuc in all_nucs:
            data[nuc] = {}
            data[nuc]['Concentration'] = concentrations[nuc].n
            data[nuc]['sigma Concentration'] = concentrations[nuc].s
            
        return data




if __name__ == "__main__":
    input_path = "../examples/keepin_1957/input.json"
    concentrations = Concentrations(input_path)
    concentrations.generate_concentrations()