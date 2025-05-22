from utils.input_handler import InputHandler
from utils.csv_handler import CSVHandler
import os
import numpy as np

class Preprocess:
    def __init__(self, input_path: str) -> None:
        """
        This class handles preprocessing of files into a common format.

        Parameters
        ----------
        input_path : str
            Path to the input file.
        """
        self.input_path: str = input_path
        self.input_handler: InputHandler = InputHandler(input_path)
        self.input_data: dict = self.input_handler.read_input()

        self.data_dir: str = self.input_data['data_directory']
        self.overwrite: bool = self.input_data['overwrite']
        self.out_dir: str = self.input_data['output_directory']
        self.fissile_targets: list = self.input_data['fissile_targets']
        self.energy_MeV: float = self.input_data['energy_MeV']
        return None
    
    def openmc_preprocess(self) -> None:
        """
        Processes OpenMC all chain_* files and cross section data using OpenMC
        """
        for fissile in self.fissile_targets:
            self._openmc_chain_preprocess(fissile)
        return None
    
    def _openmc_chain_preprocess(self, fissile: str) -> None:
        """
        Processes OpenMC all chain_* files

        Parameters
        ----------
        fissile : str
            Name of the fissile target to process.
        """
        for file in os.listdir(self.data_dir + '/chain/'):
            full_path: str = os.path.join(self.data_dir + '/chain/', file)
            file_data: dict[str: dict[str, float]] = self._process_chain_file(full_path, fissile)
            csv_path: str = self.out_dir + '/' + file.split('.')[0] + '.csv'
            CSVHandler(csv_path, self.overwrite).write_csv(file_data)
        return None
    
    def _process_chain_file(self, file: str, fissile: str) -> dict[str, dict[str, float]]:
        """
        Processes a single OpenMC chain file and returns the data as a dictionary.

        Parameters
        ----------
        file : str
            Name of the chain file to process.
        fissile : str
            Name of the fissile target to process.

        Returns
        -------
        data : dict[str, dict[str, float]]
            Dictionary containing the processed data.
        """
        import openmc.deplete
        chain: openmc.deplete.chain = openmc.deplete.Chain.from_xml(file)
        nuclides: list[openmc.deplete.Nuclide] = chain.nuclides
        nuc_dict: dict[str, int] = chain.nuclide_dict
        target_index: int = nuc_dict[fissile]
        target_nuc: openmc.deplete.Nuclide = nuclides[target_index]
        FY_data: openmc.deplete.FissionYieldDistribution = target_nuc.yield_data
        energies = FY_data.energies
        products = FY_data.products
        FY_chain: dict[str: dict[float, float]] = {}
        for product in products:
            FY_chain[product] = {}
            for energy in energies:
                FY_chain[product][energy] = FY_data[energy][product]

        fit_FY_chain: dict[str: float] = self._fit_fy_chain(FY_chain, order=1)
        chain_nucs: list[str] = list(fit_FY_chain.keys())
        nuclide_halflives: dict[str: float] = {nuc.name: nuc.half_life for nuc in nuclides}

        data: dict[str: dict[str: float]] = dict()
        for nuc in chain_nucs:
            data[nuc] = {}
            data[nuc]['half_life'] = nuclide_halflives[nuc]
            data[nuc]['FY'] = fit_FY_chain[nuc]
        
        return data
    
    def _fit_fy_chain(self, FY_chain: dict[str: dict[float, float]], order:int=1) -> dict[str: float]:
        """
        Fit the fission yield chain data.

        Parameters
        ----------
        FY_chain : dict[str: dict[float, float]]
            Dictionary containing the fission yield chain data.
        order : int
            Order of the polynomial fit.

        Returns
        -------
        fit_FY_chain : dict[str: float]
            Dictionary containing the fit fission yield chain data.
        """
        fit_FY_chain: dict[str: float] = {}
        for product in FY_chain:
            xs = list(FY_chain[product].keys())
            ys = list(FY_chain[product].values())
            coeffs = np.polyfit(xs, ys, order)
            fit_FY_chain[product] = np.polyval(coeffs, self.energy_MeV)
        return fit_FY_chain


preproc = Preprocess('../examples/keepin_1957/pre_input.json')
preproc.openmc_preprocess()