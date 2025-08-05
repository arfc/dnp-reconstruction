import numpy as np
from mosden.utils.csv_handler import CSVHandler
from uncertainties import ufloat
from mosden.base import BaseClass
from time import time


class Concentrations(BaseClass):
    def __init__(self, input_path: str) -> None:
        """
        This class handles the formation of concentrations from data.

        Parameters
        ----------
        input_path : str
            Path to the input file
        """
        super().__init__(input_path)
        self.output_dir: str = self.input_data['file_options']['output_dir']
        modeling_options: dict = self.input_data.get('modeling_options', {})
        file_options: dict = self.input_data.get('file_options', {})
        overwrite: dict = file_options.get('overwrite', {})

        self.model_method: str = modeling_options.get(
            'concentration_handling', 'CFY')
        self.overwrite: bool = overwrite.get('concentrations', False)

        self.reprocessing: dict[str: float] = modeling_options.get(
            'reprocessing', {})
        self.reprocess: bool = (sum(self.reprocessing.values()) > 0)
        self.reprocess_locations: list[str] = modeling_options.get(
            'reprocessing_locations', [])
        self.t_in: float = modeling_options.get('incore_s', 0.0)
        self.t_ex: float = modeling_options.get('excore_s', 0.0)
        self.t_net: float = modeling_options.get('net_irrad_s', 0.0)
        self.irrad_type: str = modeling_options.get('irrad_type', 'saturation')

        return None

    def generate_concentrations(self) -> None:
        """
        Generate the concentrations of each nuclide from based on
        irradiation of the sample for the irradiation times.
        """
        start = time()
        data: dict[str: dict[str: float]] = dict()
        if self.model_method == 'CFY':
            if self.t_ex > 0.0:
                raise NotImplementedError(
                    'Excore residence not available for CFY')
            if self.reprocess:
                raise NotImplementedError('Reprocessing not available for CFY')
            if self.irrad_type != 'saturation':
                self.logger.error(
                    'CFY is intended for a saturation irradiation')
            data = self.CFY_concentrations()
        elif self.model_method == 'IFY':
            if self.t_ex > 0.0:
                raise NotImplementedError(
                    'Excore residence not available for IFY')
            if self.reprocess:
                raise NotImplementedError('Reprocessing not available for IFY')
            if self.irrad_type != 'pulse':
                self.logger.error(
                    'IFY method does not use cumulative fission yields')
            self.logger.error(
                'IFY method has not been verified. Use with caution')
            data = self.IFY_concentrations()
        else:
            raise NotImplementedError(
                f"Concentration handling method '{
                    self.model_method}' is not implemented")

        CSVHandler(self.concentration_path, self.overwrite).write_csv(data)
        self.save_postproc()
        self.time_track(start, 'Concentrations')
        return

    def CFY_concentrations(self) -> None:
        """
        Generate the concentrations of each nuclide using the CFY method.
        """
        concentrations: dict[str: dict[str: ufloat]] = dict()
        all_nucs: set[str] = set()
        CFY_data = self._read_processed_data('fission_yield')
        half_life_data = self._read_processed_data('half_life')
        for nuclide in CFY_data.keys():
            concs = ufloat(
                CFY_data[nuclide]['CFY'],
                CFY_data[nuclide]['sigma CFY'])

            try:
                hl = ufloat(
                    half_life_data[nuclide]['half_life'],
                    half_life_data[nuclide]['sigma half_life'])
            except KeyError:
                continue
            lam = np.log(2) / hl
            concentrations[nuclide] = concs / lam
            all_nucs.add(nuclide)

        data: dict[str: dict[str: float]] = dict()
        for nuc in all_nucs:
            data[nuc] = {}
            data[nuc]['Concentration'] = concentrations[nuc].n
            data[nuc]['sigma Concentration'] = concentrations[nuc].s

        return data

    def IFY_concentrations(self) -> None:
        """
        Generate the concentrations of each nuclide using the IFY method.
        """
        concentrations: dict[str: dict[str: ufloat]] = dict()
        all_nucs: set[str] = set()
        IFY_data = self._read_processed_data('fission_yield')
        for nuclide in IFY_data.keys():
            concs = IFY_data[nuclide]['IFY']
            concentrations[nuclide] = concs
            all_nucs.add(nuclide)

        data: dict[str: dict[str: float]] = dict()
        for nuc in all_nucs:
            data[nuc] = {}
            data[nuc]['Concentration'] = concentrations[nuc]
            data[nuc]['sigma Concentration'] = 0.0

        return data


if __name__ == "__main__":
    input_path = "../examples/keepin_1957/input.json"
    concentrations = Concentrations(input_path)
    concentrations.generate_concentrations()
