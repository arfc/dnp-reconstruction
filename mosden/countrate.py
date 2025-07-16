from mosden.utils.input_handler import InputHandler
from mosden.utils.csv_handler import CSVHandler
from mosden.base import BaseClass
import os
import numpy as np
from uncertainties import ufloat, unumpy

class CountRate(BaseClass):
    """
    Class to handle the delayed neutron count rate calculations
    """
    def __init__(self, input_path: str) -> None:
        super().__init__(input_path)
        self.data_dir: str = self.input_data['file_options']['processed_data_dir']
        self.output_dir: str = self.input_data['file_options']['output_dir']
        self.overwrite: bool = self.input_data['file_options']['overwrite']['count_rate']

        self.parent_feed: bool = self.input_data['modeling_options']['parent_feeding']
        self.num_times: int = self.input_data['modeling_options']['num_decay_times']
        self.decay_time: float = self.input_data['modeling_options']['decay_time']
        self.decay_times: np.ndarray = np.linspace(0, self.decay_time, self.num_times)
        return None
    
    def calculate_count_rate(self) -> None:
        """
        Calculate the delayed neutron count rate from
        concentrations using various methods
        """
        data: dict[str: list[float]] = dict()

        if self.input_data['modeling_options']['parent_feeding']:
            raise NotImplementedError('Requires OpenMC depletion')
        else:
            data = self._count_rate_from_data()

        CSVHandler(self.countrate_path, self.overwrite).write_count_rate_csv(data)
        return
    
    def _count_rate_from_data(self) -> dict[str: list[float]]:
        """
        Calculate the delayed neutron count rate from existing data
        """
        data: dict[str: list[float]] = dict()
        count_rate: np.ndarray = np.zeros(len(self.decay_times))
        sigma_count_rate: np.ndarray = np.zeros(len(self.decay_times))

        any_fissile: str = list(self.fissiles.keys())[0]
        processed_data_paths= self._get_data_paths(processed=True, directory=False, fissile=any_fissile)
        emission_prob_data = CSVHandler(processed_data_paths['emission_probability'], create=False).read_csv()
        half_life_data = CSVHandler(processed_data_paths['half_life'], create=False).read_csv()
        concentration_data = CSVHandler(self.concentration_path, create=False).read_csv()

        emission_nucs = list(emission_prob_data.keys())
        half_life_nucs = list(half_life_data.keys())
        conc_nucs = list(concentration_data.keys())
        net_unique_nucs = list(set(emission_nucs+half_life_nucs+conc_nucs))
        net_similar_nucs = list(set(emission_nucs) & set(half_life_nucs) & set(conc_nucs))

        for nuc in net_similar_nucs:
            Pn_data = emission_prob_data[nuc]
            Pn = ufloat(Pn_data['emission probability'], Pn_data['sigma emission probability'])

            hl_data = half_life_data[nuc]
            try:
                halflife = ufloat(hl_data['half_life'], hl_data['sigma half_life'])
            except KeyError:
                halflife = hl_data['half_life']
            decay_const = np.log(2) / halflife

            conc_data = concentration_data[nuc]
            conc = ufloat(conc_data['Concentration'], conc_data['sigma Concentration'])

            counts = Pn * decay_const * conc * unumpy.exp(-decay_const * self.decay_times)
            count_rate += unumpy.nominal_values(counts)
            sigma_count_rate += unumpy.std_devs(counts)

        data = {
            'times': self.decay_times,
            'counts': count_rate,
            'sigma counts': sigma_count_rate
        }
        return data
    

if __name__ == '__main__':
    delayed_neutrons = CountRate('../examples/keepin_1957/input.json')
    delayed_neutrons.calculate_count_rate()