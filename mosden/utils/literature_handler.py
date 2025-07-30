from mosden.base import BaseClass
from uncertainties import ufloat
import numpy as np

class Literature(BaseClass):
    def __init__(self, input_path) -> None:
        """
        This class holds data from the literature for passing into MoSDeN post processing.
        """
        super().__init__(input_path)
        return None
    
    def get_group_data(self) -> dict[dict[str: list[float]]]:
        """
        Get countrate data from various sources and compile into a dictionary
        """
        data_holder = dict()
        self.logger.warning('Current literature energy binning could be improved')
        if self.energy_MeV < 1e-3:
            energy = 'thermal'
        else:
            energy = 'fast'
        names = ['keepin']

        for name in names:
            data_holder[name] = dict()
            for fiss, frac in self.fissiles.items():
                data_holder[name][fiss] = self._get_group_data_helper(fiss, frac, name, energy)
        data_holder = self._merge_fiss(data_holder)
        return data_holder
    
    def _merge_fiss(self, data: dict[str: dict[str: list[float]]]) -> dict[dict[str: list[float]]]:
        """
        Merge same name fissile data into a single dictionary
        """
        merged_data = dict()
        for name, name_data in data.items():
            merged_data[name] = dict()
            yields = dict()
            halflives = dict()
            for fiss, params in name_data.items():
                for i in range(len(params['yield'])):
                    yields[i] = yields.get(i, 0.0) + ufloat(params['yield'][i], params['sigma yield'][i])
                    halflives[i] = halflives.get(i, 0.0) + ufloat(params['half_life'][i], params['sigma half_life'][i])
            merged_data[name]['yield'] = [yields[i].n for i in yields]
            merged_data[name]['sigma yield'] = [yields[i].s for i in yields]
            merged_data[name]['half_life'] = [halflives[i].n for i in halflives]
            merged_data[name]['sigma half_life'] = [halflives[i].s for i in halflives]

        return merged_data


    def _get_group_data_helper(self, fiss: str, frac: float, name: str, energy: str) -> dict[str: list[float]]:
        if name == 'keepin':
            if fiss == 'U235':
                if energy == 'thermal':
                    # Parish et al. 1999
                    net_yield = ufloat(0.0158, 0.0011)
                    yields = [a*net_yield*frac for a in [ufloat(0.033, 0.003),
                                                    ufloat(0.219, 0.005),
                                                    ufloat(0.196, 0.022),
                                                    ufloat(0.395, 0.011),
                                                    ufloat(0.115, 0.009),
                                                    ufloat(0.042, 0.008)]]
                    decay_constants = [ufloat(0.0124, 0.0003),
                                    ufloat(0.0305, 0.0009),
                                    ufloat(0.111, 0.004),
                                    ufloat(0.301, 0.012),
                                    ufloat(1.14, 0.15),
                                    ufloat(3.01, 0.29)]
                    half_lives = [np.log(2)/lam*frac for lam in decay_constants]
                elif energy == 'fast':
                    # Keepin et al. 1957
                    net_yield = ufloat(0.0165, 0.0005)
                    yields = [a*net_yield*frac for a in [ufloat(0.038, 0.003),
                                                    ufloat(0.213, 0.005),
                                                    ufloat(0.188, 0.016),
                                                    ufloat(0.407, 0.007),
                                                    ufloat(0.128, 0.008),
                                                    ufloat(0.026, 0.003)]]
                    half_lives = [frac*hl for hl in [ufloat(54.51, 0.94),
                                    ufloat(21.84, 0.54),
                                    ufloat(6.00, 0.17),
                                    ufloat(2.23, 0.06),
                                    ufloat(0.496, 0.029),
                                    ufloat(0.179, 0.017)]]
            elif fiss == 'U238':
                if energy == 'thermal':
                    # Keepin et al. 1957
                    net_yield = ufloat(0.0412, 0.0017)
                    yields = [a*net_yield*frac for a in [ufloat(0.013, 0.001),
                                                    ufloat(0.137, 0.002),
                                                    ufloat(0.162, 0.020),
                                                    ufloat(0.388, 0.012),
                                                    ufloat(0.225, 0.013),
                                                    ufloat(0.075, 0.005)]]
                    half_lives = [frac*hl for hl in [ufloat(52.38, 1.29),
                                    ufloat(21.58, 0.39),
                                    ufloat(4.53, 0.19),
                                    ufloat(1.93, 0.07),
                                    ufloat(0.490, 0.023),
                                    ufloat(0.172, 0.009)]]

        try:
            group_params = {
                'yield': [a.n for a in yields],
                'sigma yield': [a.s for a in yields],
                'half_life': [hl.n for hl in half_lives],
                'sigma half_life': [hl.s for hl in half_lives]
            }
        except UnboundLocalError:
            raise KeyError(f"Data for {fiss} in {name} at {energy} energy not found")
        return group_params

            
