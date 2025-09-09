from mosden.base import BaseClass
from uncertainties import ufloat
import numpy as np

class Literature(BaseClass):
    def __init__(self, input_path) -> None:
        """
        This class holds data from the literature for passing into MoSDeN post processing.

        Parameters
        ----------
        input_path : str
            Path to the input file containing fissile data.
        """
        super().__init__(input_path)
        self.available_names: list[str] = ['keepin', 'charlton', 'endfb6', 'mills', 'saleh', 'synetos', 'tuttle', 'waldo', 'brady', 'Modified 0D Scaled']
        return None
    
    def get_group_data(self, names:list[str]=None) -> dict[dict[str: list[float]]]:
        """
        Get countrate data from various sources and compile into a dictionary

        Parameters
        ----------
        names : list[str]
            List of names for which to retrieve data. Default is ['keepin'].
        """
        data_holder = dict()
        self.logger.warning('Current literature energy binning could be improved')
        if self.energy_MeV < 1e-3:
            energy = 'thermal'
        else:
            energy = 'fast'

        if not names:
            names = self.available_names
        for name in names:
            data_holder[name] = dict()
            for fiss, frac in self.fissiles.items():
                data = self._get_group_data_helper(fiss, frac, name, energy)
                if data is None:
                    del data_holder[name]
                    continue
                data_holder[name][fiss] = data
        data_holder = self._merge_fiss(data_holder)
        return data_holder
    
    def _merge_fiss(self, data: dict[str: dict[str: list[float]]]) -> dict[dict[str: list[float]]]:
        """
        Merge same name fissile data into a single dictionary.

        Parameters
        ----------
        data : dict[str: dict[str: list[float]]]
            Dictionary containing group data for different names, fissile nuclides, and energies.
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
        """
        Helper function to retrieve group data for a specific fissile nuclide and energy.
        Parameters
        ----------
        fiss : str
            Fissile nuclide for which to retrieve data (e.g., 'U235').
        frac : float
            Fraction of the fissions from the given fissile nuclide.
        name : str
            Name of the literature source (e.g., 'keepin').
        energy : str
            Energy level ('thermal' or 'fast').

        Returns
        -------
        group_params : dict[str: list[float]]
            Group data for the specified fissile nuclide and energy.
        
        Raises
        ------
        UnboundLocalError
            If the fissile nuclide or energy level is not found in the specified literature source
        """
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
        
        elif name == 'charlton':
            if fiss == 'U235':
                if energy == 'thermal':
                    # Parish et al. 1999
                    net_yield = ufloat(0.0158, 0.0)
                    yields = [a * net_yield * frac for a in [
                        ufloat(0.034, 0.002),
                        ufloat(0.218, 0.003),
                        ufloat(0.202, 0.015),
                        ufloat(0.384, 0.008),
                        ufloat(0.118, 0.009),
                        ufloat(0.042, 0.009)]]
                    decay_constants = [ufloat(0.0125, 0.0001), ufloat(0.0309, 0.0007), ufloat(0.115, 0.002),
                                    ufloat(0.317, 0.009), ufloat(1.41, 0.09), ufloat(3.01, 0.29)]
                    half_lives = [np.log(2)/lam * frac for lam in decay_constants]

        elif name == 'endfb6':
            if fiss == 'U235':
                if energy == 'thermal':
                    # Parish et al. 1999
                    net_yield = ufloat(0.0167, 0.0)
                    yields = [a * net_yield * frac for a in [
                        ufloat(0.035, 0.0),
                        ufloat(0.1807, 0.0),
                        ufloat(0.1725, 0.0),
                        ufloat(0.3868, 0.0),
                        ufloat(0.1586, 0.0),
                        ufloat(0.0664, 0.0)]]
                    decay_constants = [ufloat(0.0134, 0.0), ufloat(0.0327, 0.0), ufloat(0.1208, 0.0),
                                    ufloat(0.3028, 0.0), ufloat(0.8495, 0.0), ufloat(2.853, 0.0)]
                    half_lives = [np.log(2)/lam * frac for lam in decay_constants]

        elif name == 'mills':
            if fiss == 'U235':
                if energy == 'thermal':
                    # Parish et al. 1999
                    #  (this actually has a typo in group 5 0.115 instead of 0.1745)
                    #  (several typos actually, group 1 uncert)
                    # Pulling from Mills et al. 1992
                    # Microscopic
                    net_yield = ufloat(0.0164484, 0.00114)
                    yields = [a * net_yield * frac for a in [
                        ufloat(0.0343, 0.0002909),
                        ufloat(0.1974, 0.002048),
                        ufloat(0.1193, 0.009091),
                        ufloat(0.4002, 0.01407),
                        ufloat(0.1745, 0.01748),
                        ufloat(0.0742, 0.005143)]]
                    decay_constants = [ufloat(0.01254, 0.0000225), ufloat(0.0304, 0.0001563), ufloat(0.09027, 0.004451),
                                    ufloat(0.2501, 0.01016), ufloat(0.6455, 0.05122), ufloat(2.46, 0.08903)]
                    half_lives = [np.log(2)/lam * frac for lam in decay_constants]

        elif name == 'saleh':
            if fiss == 'U235':
                if energy == 'thermal':
                    # Parish et al. 1999
                    net_yield = ufloat(0.0154, 0.0004)
                    yields = [a * net_yield * frac for a in [
                        ufloat(0.0354, 0.003),
                        ufloat(0.236, 0.02),
                        ufloat(0.193, 0.017),
                        ufloat(0.385, 0.067),
                        ufloat(0.109, 0.009),
                        ufloat(0.0412, 0.008)]]
                    decay_constants = [ufloat(0.0125, 0.0009), ufloat(0.0306, 0.0002), ufloat(0.111, 0.007),
                                    ufloat(0.3, 0.005), ufloat(1.10, 0.03), ufloat(3.01, 0.29)]
                    half_lives = [np.log(2)/lam * frac for lam in decay_constants]

        elif name == 'synetos':
            if fiss == 'U235':
                if energy == 'thermal':
                    # Parish et al. 1999
                    net_yield = ufloat(0.0151, 0.0007)
                    yields = [a * net_yield * frac for a in [
                        ufloat(0.034, 0.002),
                        ufloat(0.217, 0.009),
                        ufloat(0.142, 0.023),
                        ufloat(0.386, 0.017),
                        ufloat(0.179, 0.017),
                        ufloat(0.042, 0.008)]]
                    decay_constants = [ufloat(0.0123, 0.0001), ufloat(0.0310, 0.0004), ufloat(0.099, 0.011),
                                    ufloat(0.254, 0.02), ufloat(0.9, 0.11), ufloat(3.01, 0.29)]
                    half_lives = [np.log(2)/lam * frac for lam in decay_constants]

        elif name == 'tuttle':
            if fiss == 'U235':
                if energy == 'thermal':
                    # Parish et al. 1999
                    net_yield = ufloat(0.0170, 0.0002)
                    yields = [a * net_yield * frac for a in [
                        ufloat(0.038, 0.004),
                        ufloat(0.213, 0.007),
                        ufloat(0.188, 0.024),
                        ufloat(0.407, 0.010),
                        ufloat(0.128, 0.012),
                        ufloat(0.026, 0.004)]]
                    decay_constants = [ufloat(0.0127, 0.0003), ufloat(0.0317, 0.0012), ufloat(0.115, 0.004),
                                    ufloat(0.311, 0.012), ufloat(1.40, 0.12), ufloat(3.87, 0.55)]
                    half_lives = [np.log(2)/lam * frac for lam in decay_constants]

        elif name == 'waldo':
            if fiss == 'U235':
                if energy == 'thermal':
                    # Parish et al. 1999
                    net_yield = ufloat(0.0167, 0.0007)
                    yields = [a * net_yield * frac for a in [
                        ufloat(0.033, 0.003),
                        ufloat(0.219, 0.023),
                        ufloat(0.196, 0.023),
                        ufloat(0.395, 0.016),
                        ufloat(0.115, 0.009),
                        ufloat(0.042, 0.005)]]
                    decay_constants = [ufloat(0.0127, 0.0003), ufloat(0.0317, 0.0012), ufloat(0.115, 0.004),
                                    ufloat(0.311, 0.012), ufloat(1.40, 0.12), ufloat(3.87, 0.55)]
                    half_lives = [np.log(2)/lam * frac for lam in decay_constants]

        elif name == 'brady':
            if fiss == 'U235':
                if energy == 'thermal':
                    # Brady and England 1989
                    # Microscopic
                    net_yield = ufloat(0.0178, 0.001)
                    yields = [a * net_yield * frac for a in [
                        ufloat(0.038, 0.0),
                        ufloat(0.1918, 0.0),
                        ufloat(0.1638, 0.0),
                        ufloat(0.3431, 0.0),
                        ufloat(0.1744, 0.0),
                        ufloat(0.0890, 0.0)]]
                    decay_constants = [ufloat(0.0133, 0.0), ufloat(0.0325, 0.0), ufloat(0.1219, 0.0),
                                    ufloat(0.3169, 0.0), ufloat(0.9886, 0.0), ufloat(2.9544, 0.0)]
                    half_lives = [np.log(2)/lam * frac for lam in decay_constants]

        elif name == 'Modified 0D Scaled':
            if fiss == 'U235':
                if energy == 'thermal':
                    # Microscopic
                    net_yield = ufloat(0.0158, 0.0011)
                    yields = [a*net_yield*frac for a in [
                                ufloat(0.033, 0.003),
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

                    yields[0] = ufloat(0.00052, 0.00002)
                    yields[1] = ufloat(0.00297, 0.00009)
                    half_lives[0] = ufloat(55.52977, 0.15676)
                    half_lives[1] = ufloat(23.54278, 0.12560)



        try:
            group_params = {
                'yield': [a.n for a in yields],
                'sigma yield': [a.s for a in yields],
                'half_life': [hl.n for hl in half_lives],
                'sigma half_life': [hl.s for hl in half_lives]
            }
        except UnboundLocalError:
            self.logger.warning(f"Data for {fiss} in {name} at {energy} energy not found")
            return None
        return group_params

            
if __name__ == "__main__":
    input_path = "../../examples/huynh_2014/input.json"
    lit = Literature(input_path)
    data = lit.get_group_data(lit.available_names)
    target_key = 'half_life'
    target_name = None
    for name, val in data.items():
        if name != target_name and target_name is not None:
            continue
        print(name)
        for key, item_val in val.items():
            if key != target_key and target_key is not None:
                continue
            print(key)
            yield_val = sum([ufloat(val['yield'][i], val['sigma yield'][i])*100 for i in range(len(item_val))])
            halflife_val = 1/yield_val * sum([ufloat(val['yield'][i], val['sigma yield'][i])*100*ufloat(val['half_life'][i], val['sigma half_life'][i]) for i in range(len(item_val))])
            if target_key == 'yield':
                print(f'{yield_val = }')
            elif target_key == 'half_life':
                print(f'{halflife_val = }')
            for item in item_val:
                print(f'{round(item, 5):.5f} & ', end='')
            print()
        print('\n')
