import numpy as np
from mosden.utils.input_handler import InputHandler
from mosden.utils.csv_handler import CSVHandler
from pathlib import Path
from uncertainties import ufloat

class Concentrations:
    """
    Class to handle the formation of concentrations from data
    """
    def __init__(self, input_path: str) -> None:
        self.input_path: str = input_path
        self.input_handler: InputHandler = InputHandler(input_path)
        self.input_data: dict = self.input_handler.read_input()
        self.output_dir: str = self.input_data['file_options']['output_dir']

        self.model_method: str = self.input_data['modeling_options']['concentration_handling']
        self.fissiles: dict[str, float] = self.input_data['data_options']['fissile_fractions']
        self.data_path: str = self.input_data['file_options']['processed_data_dir']
        self.energy: float = self.input_data['data_options']['energy_MeV']
        self.overwrite: bool = self.input_data['file_options']['overwrite']
        self.FY_file_name: str = self.input_data['data_options']['fission_yield']['type']

        self.reprocessing: dict[str: float] = self.input_data['modeling_options']['reprocessing']
        self.reprocess: bool = (sum(self.reprocessing.values()) > 0)
        self.reprocess_locations: list[str] = self.input_data['modeling_options']['reprocessing_locations']
        self.t_in: float = self.input_data['modeling_options']['incore_s']
        self.t_ex: float = self.input_data['modeling_options']['excore_s']
        self.t_net: float = self.input_data['modeling_options']['net_irrad_s']
        self.irrad_type: str = self.input_data['modeling_options']['irrad_type']
        
        self.half_life_dir: str = self.input_data['data_options']['decay_constant']
        self.cross_section_dir: str = self.input_data['data_options']['cross_section']
        self.emission_probability_dir: str = self.input_data['data_options']['emission_probability']
        self.fission_yield_dir: str = self.input_data['data_options']['fission_yield']['data']
        return None
    
    def generate_concentrations(self) -> None:
        """
        Generate the concentrations of each nuclide from based on
        irradiation of the sample for the irradiation times.
        """
        if self.model_method == 'CFY':
            if self.t_ex > 0.0:
                raise NotImplementedError('Excore residence not available for CFY')
            if self.reprocess:
                raise NotImplementedError('Reprocessing not available for CFY')
            self.CFY_concentrations()
        else:
            raise NotImplementedError(
                f"Concentration handling method '{self.model_method}' is not implemented"
            )

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
            
        output_path = Path(self.output_dir) / f"concentrations.csv"
        CSVHandler(output_path, self.overwrite).write_csv(data)
        return

    def _read_processed_data(self, fissile: str, data_type: str) -> dict[str, dict[str, float]]:
        """
        Read the processed data for a given fissile nuclide.

        Parameters
        ----------
        fissile : str
            The name of the fissile nuclide.
        data_type : str
            The type of data to read (e.g., "fission_yield", "half_life", "cross_section", "emission_probability").

        Returns
        -------
        dict
            The processed data for the fissile nuclide.
        """
        if data_type == 'fission_yield':
            use_path = self.fission_yield_dir
        elif data_type == 'half_life':
            use_path = self.half_life_dir
        elif data_type == 'cross_section':
            use_path = self.cross_section_dir
        elif data_type == 'emission_probability':
            use_path = self.emission_probability_dir
        else:
            raise ValueError(f"Unknown data type: {data_type}")

        data_path = Path(self.data_path) / use_path / fissile / f"{self.energy}MeV" / self.FY_file_name
        csv_handler = CSVHandler(data_path, create=False)
        if not csv_handler._file_exists():
            raise FileNotFoundError(f"Processed data file {data_path} does not exist.") 
        data = csv_handler.read_csv()
        return data



if __name__ == "__main__":
    input_path = "../examples/keepin_1957/input.json"
    concentrations = Concentrations(input_path)
    concentrations.generate_concentrations()