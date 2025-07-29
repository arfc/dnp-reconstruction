import numpy as np
from mosden.utils.input_handler import InputHandler
from mosden.utils.csv_handler import CSVHandler
from pathlib import Path
from uncertainties import ufloat, unumpy
from mosden.base import BaseClass
from scipy.optimize import least_squares, curve_fit
from typing import Callable
from math import ceil
from time import time
import warnings
from logging import INFO

class Grouper(BaseClass):
    """
    Class to handle the formation of groups from count rates
    """
    def __init__(self, input_path: str) -> None:
        super().__init__(input_path)
        self.output_dir: str = self.input_data['file_options']['output_dir']

        self.model_method: str = self.input_data['group_options']['method']
        self.num_groups: int = self.input_data['group_options']['num_groups']
        self.overwrite: bool = self.input_data['file_options']['overwrite']['group_fitting']
        self.MC_samples: int = self.input_data['group_options']['samples']

        self.t_in: float = self.input_data['modeling_options']['incore_s']
        self.t_ex: float = self.input_data['modeling_options']['excore_s']
        self.t_net: float = self.input_data['modeling_options']['net_irrad_s']
        self.irrad_type: str = self.input_data['modeling_options']['irrad_type']
        self.sample_func: str = self.input_data['group_options']['sample_func']
        return None
    
    def generate_groups(self) -> None:
        """
        Generate some number of groups based on the selected method
        """
        start = time()
        data: dict[str: dict[str: float]] = dict()
        if self.model_method == 'nlls':
            data = self._nonlinear_least_squares()
        else:
            raise NotImplementedError(f'{self.model_method} is not implemented')
        CSVHandler(self.group_path, self.overwrite).write_groups_csv(data, sortby='half_life')
        self.save_postproc()
        self.time_track(start, 'Groupfit')
        return None
    
    def _residual_function(self, parameters: np.ndarray[float], times: np.ndarray[float], counts: np.ndarray[float],
                           fit_func : Callable) -> float:
        """
        Calculate the residual of the current set of parameters

        Parameters
        ----------
        parameters : np.ndarray[float]
            Half life and yield parameters [yield1, yield2, ..., h1, h2, ...]
        times : np.ndarray[float]
            List of times
        counts : np.ndarray[float]
            List of nominal times
        fit_func : Callable
            Function that takes times and parameters to return list of counts

        Returns
        -------
        residual : float
            Value of the residual
        """
        residual = (counts - fit_func(times, parameters)) / (counts + 1e-12)
        return residual
    
    def _pulse_fit_function(self, times: np.ndarray[float], parameters: np.ndarray[float]) -> np.ndarray[float|object]:
        yields = parameters[:self.num_groups]
        half_lives = parameters[self.num_groups:]
        counts: np.ndarray[float] = np.zeros(len(times))
        for group in range(self.num_groups):
            lam = np.log(2) / half_lives[group]
            a = yields[group]
            try:
                counts += (a * lam * np.exp(-lam * times))
            except TypeError:
                if group == 0:
                    counts: np.ndarray[object] = np.zeros(len(times), dtype=object)
                counts += (a * lam * unumpy.exp(-lam * times))
        return counts

    def _saturation_fit_function(self, times: np.ndarray[float], parameters: np.ndarray[float]) -> np.ndarray[float|object]:
        yields = parameters[:self.num_groups]
        half_lives = parameters[self.num_groups:]
        counts: np.ndarray[float] = np.zeros(len(times))
        tot_cycles: int = ceil(self.t_net / (self.t_in + self.t_ex))
        for group in range(self.num_groups):
            lam = np.log(2) / half_lives[group]
            a = yields[group]
            cycle_sum = 0
            for j in range(1, tot_cycles+1):
                try:
                    cycle_sum += np.exp(-lam * (self.t_net - j*self.t_in - (j-1)*self.t_ex))
                except TypeError:
                    cycle_sum += unumpy.exp(-lam * (self.t_net - j*self.t_in - (j-1)*self.t_ex))
            try:
                counts += a * np.exp(-lam * times) * (1 - np.exp(-lam * self.t_net + (1 - np.exp(lam * self.t_ex) * cycle_sum)))
            except TypeError:
                if group == 0:
                    counts: np.ndarray[object] = np.zeros(len(times), dtype=object)
                counts += a * unumpy.exp(-lam * times) * (1 - unumpy.exp(-lam * self.t_net + (1 - unumpy.exp(lam * self.t_ex) * cycle_sum)))
        return counts
    
    def _nonlinear_least_squares(self, count_data: dict[str: np.ndarray[float]] = None) -> dict[str: dict[str: float]]:
        """
        Run nonlinear least squares fit on the delayed neutron count rate curve
        to generate group half-lives and yields
        """
        from mosden.countrate import CountRate
        initial_parameter_guess = np.ones(self.num_groups*2)
        if count_data == None:
            count_data = CSVHandler(self.countrate_path).read_vector_csv()
        times = np.asarray(count_data['times'])
        counts = np.asarray(count_data['counts'])
        count_err = np.asarray(count_data['sigma counts'])
        if self.irrad_type == 'pulse':
            fit_function = self._pulse_fit_function
        elif self.irrad_type == 'saturation':
            fit_function = self._saturation_fit_function
        else:
            raise NotImplementedError(f'{self.irrad_type} not supported in nonlinear least squares solver')
        
        min_half_life = 1e-3
        max_half_life = 1e3
        max_yield = 1.0
        lower_bounds = np.concatenate((np.zeros(self.num_groups), np.ones(self.num_groups) * min_half_life))
        upper_bounds = np.concatenate((np.ones(self.num_groups) * max_yield, np.ones(self.num_groups) * max_half_life))

        bounds = (lower_bounds, upper_bounds) 
        result = least_squares(self._residual_function,
                               initial_parameter_guess,
                               bounds=bounds,
                               method='trf',
                               ftol=1e-12,
                               gtol=1e-12,
                               xtol=1e-12,
                               verbose=0,
                               max_nfev=1e5,
                               args=(times, counts, fit_function))

        sampled_params: list[float] = list()
        tracked_counts: list[float] = list()
        sampled_params.append(result.x)
        countrate = CountRate(self.input_path)
        self.logger.info(f'Currently using {self.sample_func} sampling')
        for _ in range(self.MC_samples):
            if self.MC_samples == 1:
                break
            data = countrate.calculate_count_rate(MC_run=True, sampler_func=self.sample_func)
            count_sample = data['counts']
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                result = least_squares(self._residual_function,
                                    result.x,
                                    bounds=bounds,
                                    method='trf',
                                    ftol=1e-12,
                                    gtol=1e-12,
                                    xtol=1e-12,
                                    verbose=0,
                                    max_nfev=1e3,
                                    args=(times, count_sample, fit_function))
            tracked_counts.append([i for i in count_sample])
            sampled_params.append(result.x)
        sampled_params: np.ndarray[float] = np.asarray(sampled_params)
        param_means: np.ndarray[float] = np.mean(sampled_params, axis=0)
        param_stds: np.ndarray[float] = np.std(sampled_params, axis=0)

        param_data = dict()
        param_data['groupfitMC'] = list()
        for iterval in range(self.MC_samples):
            param_data['groupfitMC'].append([i for i in sampled_params[iterval]])
        self.post_data['groupfitMC'] = param_data
        self.post_data['countsMC'] = tracked_counts

        data: dict[str: dict[str: float]] = dict()
        for group in range(self.num_groups):
            data[group] = dict()
            data[group]['yield'] = param_means[group]
            data[group]['sigma yield'] = param_stds[group]
            data[group]['half_life'] = param_means[self.num_groups+group]
            data[group]['sigma half_life'] = param_stds[self.num_groups+group]
        return data

if __name__ == "__main__":
    input_path = "../examples/keepin_1957/input.json"
    groupcalc = Grouper(input_path)
    groupcalc.generate_groups()
    data = CSVHandler(groupcalc.group_path).read_vector_csv()
    print(data)