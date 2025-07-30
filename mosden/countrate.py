from mosden.utils.input_handler import InputHandler
from mosden.utils.csv_handler import CSVHandler
from mosden.base import BaseClass
import os
import numpy as np
from uncertainties import ufloat, unumpy
from time import time
from typing import Callable

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

        self.irrad_type: str = self.input_data['modeling_options']['irrad_type']

        return None
    
    def calculate_count_rate(self, MC_run: bool=False, sampler_func: str=None) -> dict[str: list[float]]:
        """
        Calculate the delayed neutron count rate from
        concentrations using various methods
        """
        start = time()
        data: dict[str: list[float]] = dict()
        if self.method == 'data':
            self.emission_prob_data = CSVHandler(os.path.join(self.processed_data_dir, 'emission_probability.csv'), create=False).read_csv()
            self.half_life_data = CSVHandler(os.path.join(self.processed_data_dir, 'half_life.csv'), create=False).read_csv()
            self.concentration_data = CSVHandler(self.concentration_path, create=False).read_csv()
            data = self._count_rate_from_data(MC_run, sampler_func)
        elif self.method == 'groupfit':
            self.group_params = CSVHandler(self.group_path, create=False).read_vector_csv()
            data = self._count_rate_from_groups()
        else:
            raise NotImplementedError(f'{self.method} not available')

        if not MC_run:
            CSVHandler(self.countrate_path, self.overwrite).write_count_rate_csv(data)
            self.save_postproc()
            self.time_track(start, 'Countrate')
        return data
    
    def _count_rate_from_groups(self) -> dict[str: list[float]]:
        """
        Calculate the delayed neutron count rate from group parameters
        """
        from mosden.groupfit import Grouper
        data: dict[str: list[float]] = dict()
        count_rate: np.ndarray = np.zeros(len(self.decay_times))
        sigma_count_rate: np.ndarray = np.zeros(len(self.decay_times))
        grouper = Grouper(self.input_path)
        if self.irrad_type == 'pulse':
            fit_function = grouper._pulse_fit_function
        elif self.irrad_type == 'saturation':
            fit_function = grouper._saturation_fit_function
        else:
            raise NotImplementedError(f'{self.irrad_type} not supported in nonlinear least squares solver')

        parameters = np.zeros(grouper.num_groups*2, dtype=object)
        for i in range(grouper.num_groups):
            yield_val = ufloat(self.group_params['yield'][i], self.group_params['sigma yield'][i])
            half_life = ufloat(self.group_params['half_life'][i], self.group_params['sigma half_life'][i])
            parameters[i] = yield_val
            parameters[grouper.num_groups + i] = half_life

        counts = fit_function(self.decay_times, parameters)
        count_rate = np.asarray(unumpy.nominal_values(counts), dtype=float)
        sigma_count_rate = np.asarray(unumpy.std_devs(counts), dtype=float)

        data = {
            'times': self.decay_times,
            'counts': count_rate,
            'sigma counts': sigma_count_rate
        }
        return data
    
    def _count_rate_from_data(self, MC_run: bool=False, sampler_func: str=None) -> dict[str: list[float]]:
        """
        Calculate the delayed neutron count rate from existing data
        """
        data: dict[str: list[float]] = dict()
        count_rate: np.ndarray = np.zeros(len(self.decay_times))
        sigma_count_rate: np.ndarray = np.zeros(len(self.decay_times))

        emission_nucs = list(self.emission_prob_data.keys())
        half_life_nucs = list(self.half_life_data.keys())
        conc_nucs = list(self.concentration_data.keys())
        net_unique_nucs = list(set(emission_nucs+half_life_nucs+conc_nucs))
        net_similar_nucs = list(set(emission_nucs) & set(half_life_nucs) & set(conc_nucs))

        if len(net_similar_nucs) == 0:
            raise Exception('Error: no data exists for given emission, half life, and concentration data')

        for nuc in net_similar_nucs:
            Pn_data = self.emission_prob_data[nuc]
            Pn = ufloat(Pn_data['emission probability'], Pn_data['sigma emission probability'])

            hl_data = self.half_life_data[nuc]
            try:
                halflife = ufloat(hl_data['half_life'], hl_data['sigma half_life'])
            except KeyError:
                self.logger.warning('Half-life does not have uncertainties')
                halflife = hl_data['half_life']
            decay_const = np.log(2) / halflife

            conc_data = self.concentration_data[nuc]
            conc = ufloat(conc_data['Concentration'], conc_data['sigma Concentration'])

            if MC_run and sampler_func:
                if sampler_func == 'normal':
                    Pn = np.random.normal(Pn.n, Pn.s)
                    decay_const  = np.random.normal(decay_const.n, decay_const.s)
                    conc = np.random.normal(conc.n, conc.s)
                elif sampler_func == 'uniform':
                    Pn = np.random.uniform(Pn.n-Pn.s, Pn.n+Pn.s)
                    decay_const = np.random.uniform(decay_const.n-decay_const.s,
                                                    decay_const.n+decay_const.s)
                    conc = np.random.uniform(conc.n-conc.s,
                                             conc.n+conc.s)
                else:
                    raise NotImplementedError(f'{sampler_func} not available')
                if conc < 0.0:
                    conc = 0.0
                if decay_const < 0.0:
                    decay_const = 0.0
                if Pn < 0.0:
                    Pn = 0.0

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