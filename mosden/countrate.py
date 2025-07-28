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
        self.processed_data_dir: str = self.input_data['file_options']['processed_data_dir']
        self.output_dir: str = self.input_data['file_options']['output_dir']
        self.overwrite: bool = self.input_data['file_options']['overwrite']['count_rate']

        self.parent_feed: bool = self.input_data['modeling_options']['parent_feeding']
        self.num_times: int = self.input_data['modeling_options']['num_decay_times']
        self.decay_time: float = self.input_data['modeling_options']['decay_time']
        self.decay_times: np.ndarray = np.linspace(0, self.decay_time, self.num_times)
        self.method = self.input_data['modeling_options']['count_rate_handling']
        return None
    
    def calculate_count_rate(self) -> None:
        """
        Calculate the delayed neutron count rate from
        concentrations using various methods
        """
        data: dict[str: list[float]] = dict()
        if self.method == 'data':
            data = self._count_rate_from_data()
        else:
            raise NotImplementedError(f'{self.method} not available in countrate')

        CSVHandler(self.countrate_path, self.overwrite).write_count_rate_csv(data)
        self.save_postproc()
        return
    
    def _count_rate_from_data(self) -> dict[str: list[float]]:
        """
        Calculate the delayed neutron count rate from existing data
        """
        data: dict[str: list[float]] = dict()
        count_rate: np.ndarray = np.zeros(len(self.decay_times))
        sigma_count_rate: np.ndarray = np.zeros(len(self.decay_times))

        emission_prob_data = CSVHandler(os.path.join(self.processed_data_dir, 'emission_probability.csv'), create=False).read_csv()
        half_life_data = CSVHandler(os.path.join(self.processed_data_dir, 'half_life.csv'), create=False).read_csv()
        concentration_data = CSVHandler(self.concentration_path, create=False).read_csv()

        emission_nucs = list(emission_prob_data.keys())
        half_life_nucs = list(half_life_data.keys())
        conc_nucs = list(concentration_data.keys())
        net_unique_nucs = list(set(emission_nucs+half_life_nucs+conc_nucs))
        net_similar_nucs = list(set(emission_nucs) & set(half_life_nucs) & set(conc_nucs))

        nuc_data: dict[dict[str: list[str]]] = {
            "emission_nucs": emission_nucs,
            "half_life_nucs": half_life_nucs,
            "conc_nucs": conc_nucs,
            "net_unique_nucs": net_unique_nucs,
            "net_similar_nucs": net_similar_nucs
        }

        self.post_data.append(nuc_data)

        if len(net_similar_nucs) == 0:
            raise Exception('Error: no data exists for given emission, half life, and concentration data')

        for nuc in net_similar_nucs:
            Pn_data = emission_prob_data[nuc]
            Pn = ufloat(Pn_data['emission probability'], Pn_data['sigma emission probability'])

            hl_data = half_life_data[nuc]
            try:
                halflife = ufloat(hl_data['half_life'], hl_data['sigma half_life'])
            except KeyError:
                self.logger.warning('Half-life does not have uncertainties')
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