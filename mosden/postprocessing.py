from mosden.base import BaseClass
import matplotlib.pyplot as plt
from mosden.utils.csv_handler import CSVHandler
from mosden.countrate import CountRate
import os
import numpy as np
from uncertainties import ufloat
from logging import INFO

class PostProcess(BaseClass):
    """
    Class to handle the postprocessing from the postprocessing json file
    """
    def __init__(self, input_path: str) -> None:
        super().__init__(input_path)
        self.processed_data_dir: str = self.input_data['file_options']['processed_data_dir']
        self.output_dir: str = self.input_data['file_options']['output_dir']
        self.overwrite: bool = self.input_data['file_options']['overwrite']['postprocessing']
        return None
    
    def run(self) -> None:
        self.compare_yields()
        self.MC_NLLS_analysis()
        return None
    
    def MC_NLLS_analysis(self) -> None:
        """
        Analyze Monte Carlo Non-linear Least Squares results
        """
        self._plot_counts()
        self._plot_MC_group_params()
        return None

    
    def compare_yields(self) -> None:
        """
        Compare the total DN yields from summing individuals and from group parameters
        """
        summed_yield: float = self._get_summed_yield()
        group_yield: float = self._get_group_yield()
        self.logger.info(f'Summed yield: {summed_yield}')
        self.logger.info(f'Group yield {group_yield}')
        return None
    
    def _plot_counts(self) -> None:
        counts = self.post_data[self.names['countsMC']]
        countrate = CountRate(self.input_path)
        times = countrate.decay_times
        for MC_iterm, count_val in enumerate(counts):
            label = 'Sampled' if MC_iterm == 0 else None
            plt.plot(times, count_val, alpha=0.1, color='r', label=label)
        count_data = CSVHandler(self.countrate_path).read_vector_csv()['counts']
        plt.plot(times, count_data, color='black', linestyle='', marker='x', label='Mean', markersize=5, markevery=5)
        countrate.method = 'groupfit'
        group_counts = countrate.calculate_count_rate()
        plt.plot(times, group_counts['counts'], color='blue', alpha=0.75, label='Group Fit', linestyle='--')
        plt.xlabel('Time [s]')
        plt.ylabel(r'Count Rate $[n \cdot s^{-1}]$')
        plt.yscale('log')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}MC_counts.png')
        plt.close() 
        return None
    
    def _get_group_yield(self) -> float:
        group_data = CSVHandler(self.group_path, create=False).read_vector_csv()
        yields = [ufloat(y, std) for y, std in zip(group_data['yield'], group_data['sigma yield'])]
        net_yield = sum(yields)
        net_yield = ufloat(round(net_yield.n, 5), round(net_yield.s, 5))
        return net_yield

    
    def _get_summed_yield(self, num_top: int=10) -> float:
        nuc_yield: dict[str: float] = dict()
        emission_prob_data = CSVHandler(os.path.join(self.processed_data_dir, 'emission_probability.csv'), create=False).read_csv()
        concentration_data = CSVHandler(self.concentration_path, create=False).read_csv()
        emission_nucs = list(emission_prob_data.keys())
        conc_nucs = list(concentration_data.keys())
        net_nucs = list(set(emission_nucs) & set(conc_nucs))

        for nuc in net_nucs:
            emission_data = emission_prob_data[nuc]
            Pn = ufloat(emission_data['emission probability'],
                        emission_data['sigma emission probability'])
            conc_data = concentration_data[nuc]
            N = ufloat(conc_data['Concentration'],
                       conc_data['sigma Concentration'])
            nuc_yield[nuc] = Pn * N
        
        sorted_yields = dict(sorted(nuc_yield.items(), key=lambda item: item[1].n, reverse=True))
        net_yield = sum([i for i in sorted_yields.values()])
        extracted_vals = dict()
        running_sum = 0
        sizes = list()
        labels = list()
        if self.log_level <= INFO:
            counter = 0
            self.logger.info(f'Writing individual nuclide emission times concentration (net yield)')
            for nuc, yield_val in sorted_yields.items():
                self.logger.info(f'{nuc} - {round(yield_val.n, 3)} +/- {round(yield_val.s, 3)}')
                sizes.append(yield_val.n)
                labels.append(nuc)
                running_sum += yield_val
                counter += 1
                extracted_vals[nuc] = yield_val
                if counter > num_top:
                    break
            self.logger.info(f'Finished writing individual nuclide emission times concentration (net yield)')
            remainder = net_yield.n - running_sum.n
            sizes.append(remainder)
            labels.append('Other')
            colormap = plt.cm.rainbow
            colors = [colormap(i) for i in np.linspace(0.15, 1, num_top+2)]
            fig, ax = plt.subplots()
            ax.pie(sizes, labels=labels, autopct='%1.1f%%',
                pctdistance=0.7, labeldistance=1.1,
                colors=colors)
            ax.axis('equal')
            plt.tight_layout()
            fig.savefig(f'{self.output_dir}dnp_yield.png')
            plt.close()

            labels = [i.capitalize() for i in self.fissiles.keys()]
            sizes = list(self.fissiles.values())
            remainder = 1 - sum(sizes)
            if remainder > 0.0:
                labels.append('Other')
                sizes.append(remainder)
            colors = [colormap(i) for i in np.linspace(0.15, 1, len(labels))]
            fig, ax = plt.subplots()
            _, _, _ = ax.pie(sizes, labels=labels, autopct='%1.1f%%',
                pctdistance=0.7, labeldistance=1.1,
                colors=colors, textprops={'fontsize': 12})

            ax.axis('equal')
            
            plt.tight_layout()
            fig.savefig(f'{self.output_dir}fission_fraction.png')
            plt.close()
            
        return net_yield




