from mosden.utils.input_handler import InputHandler
from mosden.utils.csv_handler import CSVHandler
import os
import numpy as np
import re
from uncertainties import ufloat
import pandas as pd

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

        self.data_dir: str = self.input_data['file_options']['unprocessed_data_dir']
        self.overwrite: bool = self.input_data['file_options']['overwrite']['preprocessing']
        self.out_dir: str = self.input_data['file_options']['processed_data_dir']
        self.fissile_targets: list = list(self.input_data['data_options']['fissile_fractions'].keys())
        self.energy_MeV: float = self.input_data['data_options']['energy_MeV']

        self.half_life_dir: str = self.input_data['data_options']['decay_constant']
        self.cross_section_dir: str = self.input_data['data_options']['cross_section']
        self.emission_probability_dir: str = self.input_data['data_options']['emission_probability']
        self.fission_yield_dir: str = self.input_data['data_options']['fission_yield']['data']
        return None
    
    def run(self) -> None:
        """
        Run the preprocessing steps.
        """
        self.openmc_preprocess()
        self.endf_preprocess()
        self.iaea_preprocess()
        return None
    
    def openmc_preprocess(self) -> None:
        """
        Processes OpenMC all chain_* files and cross section data using OpenMC
        """
        for fissile in self.fissile_targets:
            self._openmc_chain_preprocess(fissile)
        return None
    
    def endf_preprocess(self) -> None:
        """
        Processes ENDF data
        """
        for fissile in self.fissile_targets:
            self._endf_nfy_preprocess(fissile)
        return None
    
    def iaea_preprocess(self) -> None:
        """
        Processes IAEA data for emission probabilities and half-lives.
        """
        self._iaea_dn_preprocess()
        return None

    def _iaea_dn_preprocess(self) -> None:
        """
        Processes IAEA data for emission probabilities and half-lives.

        Parameters
        ----------
        fissile : str
            Name of the fissile target to process, not needed here but kept for pathing.
        """
        full_path: str = os.path.join(self.data_dir, 'iaea', 'eval.csv')
        df = pd.read_csv(full_path, header=1)
        print(df.columns)
        # nuc, half life (uncertainty), emission probability (uncertainty)
        data: dict[str, dict[str, float]] = {}
        for _, row in df.iterrows():
            iaea_nuc = row['nucid']
            nuc = self._iaea_to_mosden_nuc(iaea_nuc)
            half_life = row[' T1/2 [s] ']
            half_life_uncertainty = row[' D T1/2 [s]']
            if row['D pn1'] < 0:
                row['D pn1'] = 0
            if row['D pn2 '] < 0:
                row['D pn2 '] = 0
            if row['D pn3'] < 0:
                row['D pn3'] = 0
            P1 = ufloat(row[' pn1 % ']/100, row['D pn1']/100)
            P2 = ufloat(row[' pn2 % ']/100, row['D pn2 ']/100)
            P3 = ufloat(row[' pn3 % ']/100, row['D pn3']/100)
            emission_prob = 1*P1 + 2*P2 + 3*P3
            data[nuc] = {}
            data[nuc]['half_life'] = half_life
            data[nuc]['sigma half_life'] = half_life_uncertainty
            data[nuc]['emission probability'] = emission_prob.n
            data[nuc]['sigma emission probability'] = emission_prob.s

        csv_path: str = os.path.join(self.out_dir, 'iaea', 'eval.csv')
        CSVHandler(csv_path, self.overwrite).write_csv(data)
        return None
    
    def _iaea_to_mosden_nuc(self, iaea_nuc: str) -> str:
        """
        Converts IAEA nuclide format to MoSDeN format.

        Parameters
        ----------
        iaea_nuc : str
            IAEA nuclide identifier.

        Returns
        -------
        str
            MoSDeN formatted nuclide identifier.
        """
        i = 0
        while i < len(iaea_nuc) and iaea_nuc[i].isdigit():
            i += 1
        mass = iaea_nuc[:i]
        element = iaea_nuc[i:].capitalize()
        return f"{element}{mass}"

    def _openmc_chain_preprocess(self, fissile: str) -> None:
        """
        Processes OpenMC all chain_* files

        Parameters
        ----------
        fissile : str
            Name of the fissile target to process.
        """
        data_dir: str = os.path.join(self.data_dir, self.half_life_dir, 'omcchain')
        out_dir: str = os.path.join(self.out_dir, self.half_life_dir, fissile, f'{self.energy_MeV}MeV')
        for file in os.listdir(data_dir):
            full_path: str = os.path.join(data_dir, file)
            file_data: dict[str: dict[str: float]] = self._process_chain_file(full_path, fissile)
            csv_path: str = os.path.join(out_dir, file.split('.')[0] + '.csv')
            CSVHandler(csv_path, self.overwrite).write_csv(file_data)
        return None
    
    def _endf_nfy_preprocess(self, fissile: str) -> None:
        """
        Processes ENDF data for the specified fissile target.

        Parameters
        ----------
        fissile : str
            Name of the fissile target to process.
        """
        data_dir: str = os.path.join(self.data_dir, self.fission_yield_dir, 'nfy')
        out_dir: str = os.path.join(self.out_dir, self.fission_yield_dir, fissile, f'{self.energy_MeV}MeV')
        for file in os.listdir(data_dir):
            adjusted_fissile: str = self._endf_fissile_name(fissile)
            if not file.startswith(f'nfy-{adjusted_fissile}'):
                continue
            full_path: str = os.path.join(data_dir, file)
            file_data : dict[str: dict[str: float]] = self._process_endf_nfy_file(full_path)
            csv_path: str = os.path.join(out_dir, 'nfy.csv')
            CSVHandler(csv_path, self.overwrite).write_csv(file_data) 
        return None
    
    def _endf_fissile_name(self, fissile: str) -> str:
        """
        Adjusts the fissile target name for ENDF processing.

        Parameters
        ----------
        fissile : str
            Name of the fissile target to adjust.

        Returns
        -------
        adjusted_fissile : str
            Adjusted fissile target name.
        """
        match = re.match(r'([A-Za-z]+)(\d+)', fissile)
        if not match:
            raise ValueError(f"Invalid nuclide format: {fissile}")

        symbol, A = match.groups()
        A = int(A)
        
        periodic_table = {
            'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
            'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
            'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26,
            'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
            'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
            'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
            'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58,
            'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66,
            'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74,
            'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
            'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90,
            'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
            'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105,
            'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112,
            'Fl': 114, 'Lv': 116
        }

        symbol = symbol.capitalize() if len(symbol) == 1 else symbol[0].upper() + symbol[1:].lower()

        Z = periodic_table.get(symbol)
        if Z is None:
            raise ValueError(f"Unknown element symbol: {symbol}")
        
        return f"{Z:03d}_{symbol}_{A}"
    
    def _process_endf_nfy_file(self, file: str) -> dict[str, dict[str: float]]:
        """
        Processes a single ENDF NFY file and returns the data as a dictionary.

        Parameters
        ----------
        file : str
            Name of the NFY file to process.
        fissile : str
            Name of the fissile target to process.

        Returns
        -------
        data : dict[str, dict[str: float]]
            Dictionary containing the processed data.
        """
        import openmc.data
        fpys = openmc.data.FissionProductYields(openmc.data.endf.Evaluation(file))
        energies = fpys.energies
        fys = fpys.cumulative
        endf_nucs: list = list(fys[0].keys())
        fit_FY_nfy = self._fit_fy_endf(energies, fys)

        data: dict[str: dict[str: float]] = dict()
        for nuc in endf_nucs:
            data[nuc] = {}
            data[nuc]['CFY'] = fit_FY_nfy[nuc].n
            data[nuc]['sigma CFY'] = fit_FY_nfy[nuc].s
        return data
        
    
    def _process_chain_file(self, file: str, fissile: str) -> dict[str, dict[str: float]]:
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
        FY_chain: dict[str: dict[float: float]] = {}
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
            if nuclide_halflives[nuc] != None:
                data[nuc]['half_life'] = nuclide_halflives[nuc]
            else:
                data[nuc]['half_life'] = np.inf
            
            data[nuc]['IFY'] = fit_FY_chain[nuc]
        
        return data

    def _fit_fy_endf(self, energies: list[float], fys: list[dict[str: ufloat]]) -> dict[str: ufloat]:
        """
        Fit the fission yield data from ENDF files.
        Uses the closest fit, not an interpolation.

        Parameters
        ----------
        energies : list[float]
            List of energy values.
        ifys : list[dict[str: ufloat]]
            Dictionary of fission yields or uncertainties at energy indices.

        Returns
        -------
        dict[str: ufloat]
            Dictionary containing the fitted fission yield data.
        """
        fit_FY_endf: dict[str: float] = {}
        endf_nucs: list[str] = list(fys[0].keys())
        energy_index: int = np.argmin(np.abs(np.array(energies) - self.energy_MeV * 1e6))
        for i, nuc in enumerate(endf_nucs):
            fission_yield = fys[energy_index][nuc]
            fit_FY_endf[nuc] = ufloat(fission_yield.n, fission_yield.s)
        return fit_FY_endf

    def _fit_fy_chain(self, FY_chain: dict[str: dict[float: float]], order:int=1) -> dict[str: float]:
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
            fit_FY_chain[product] = self._energy_fit(energies=list(FY_chain[product].keys()),
                                                     values=list(FY_chain[product].values()), 
                                                     order=order)
        return fit_FY_chain
    
    def _energy_fit(self, energies: list[float], values: list[float], order: int = 1) -> float:
        """
        Fit the energy values to a polynomial of the specified order.
        Evaluates at `self.energy_MeV`.

        Parameters
        ----------
        energies : list[float]
            List of energy values.
        values : list[float]
            List of values corresponding to the energies.
        order : int
            Order of the polynomial fit.

        Returns
        -------
        fit_value : float
            The fitted value at `self.energy_MeV`.
        """
        xs = [energy*1e-6 for energy in energies]
        if order == 1:
            fit_value = np.interp(self.energy_MeV, xs, values)
        else:
            coeffs = np.polyfit(xs, values, order)
            fit_value = np.polyval(coeffs, self.energy_MeV)
        return fit_value


if __name__ == "__main__":
    preproc = Preprocess('../examples/keepin_1957/input.json')
    preproc.run()