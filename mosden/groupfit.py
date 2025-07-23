import numpy as np
from mosden.utils.input_handler import InputHandler
from mosden.utils.csv_handler import CSVHandler
from pathlib import Path
from uncertainties import ufloat
from mosden.base import BaseClass
from scipy.optimize import least_squares, curve_fit
from typing import Callable

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
        return None
    
    def generate_groups(self) -> None:
        """
        Generate some number of groups based on the selected method
        """
        data: dict[str: dict[str: float]] = dict()
        if self.model_method == 'nlls':
            data = self._nonlinear_least_squares()
        else:
            raise NotImplementedError(f'{self.model_method} is not implemented')
        CSVHandler(self.group_path, self.overwrite).write_groups_csv(data, sortby='half_life')
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
        residual = (counts - fit_func(times, parameters)) / counts
        return residual
    
    def _pulse_fit_function(self, times: np.ndarray[float], parameters: np.ndarray[float]) -> np.ndarray[float]:
        yields = parameters[:self.num_groups]
        half_lives = parameters[self.num_groups:]
        counts: np.ndarray[float] = np.zeros(len(times))
        for group in range(self.num_groups):
            lam = np.log(2) / half_lives[group]
            a = yields[group]
            counts += (a * lam * np.exp(-lam * times))
        return counts

    def _saturation_fit_function(self, times: np.ndarray[float], parameters: np.ndarray[float]) -> np.ndarray[float]:
        yields = parameters[:self.num_groups]
        half_lives = parameters[self.num_groups:]
        counts: np.ndarray[float] = np.zeros(len(times))
        for group in range(self.num_groups):
            lam = np.log(2) / half_lives[group]
            a = yields[group]
            counts += (a * np.exp(-lam * times))
        return counts
 
    def _nonlinear_least_squares(self) -> dict[str: dict[str: float]]:
        """
        Run nonlinear least squares fit on the delayed neutron count rate curve
        to generate group half-lives and yields
        """
        initial_parameter_guess = np.ones(self.num_groups*2)
        count_data = CSVHandler(self.countrate_path).read_countrate_csv()
        times = np.asarray(count_data['times'])
        counts = np.asarray(count_data['counts'])
        count_err = np.asarray(count_data['sigma counts'])
        if self.irrad_type == 'pulse':
            fit_function = self._pulse_fit_function
        elif self.irrad_type == 'saturation':
            fit_function = self._saturation_fit_function
        else:
            raise NotImplementedError(f'{self.irrad_type} not supported in nonlinear least squares solver')
        
        result = least_squares(self._residual_function,
                               initial_parameter_guess,
                               bounds=(0, 1000),
                               method='trf',
                               gtol=1e-8,
                               verbose=1,
                               args=(times, counts, fit_function))
        # Add uncertainty calculation
        data: dict[str: dict[str: float]] = dict()
        group_params = result.x
        print(group_params)
        for group in range(self.num_groups):
            data[group] = dict()
            data[group]['yield'] = result.x[group]
            data[group]['half_life'] = result.x[self.num_groups+group]
        return data

if __name__ == "__main__":
    input_path = "../examples/keepin_1957/input.json"
    groupcalc = Grouper(input_path)
    groupcalc.generate_groups()