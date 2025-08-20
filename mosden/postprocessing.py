from logging import INFO
from uncertainties import ufloat, unumpy
import numpy as np
import os
from mosden.utils.literature_handler import Literature
from mosden.countrate import CountRate
from mosden.utils.csv_handler import CSVHandler
from mosden.base import BaseClass
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
plt.style.use('mosden.plotting')


class PostProcess(BaseClass):
    def __init__(self, input_path: str) -> None:
        """
        This class handles the postprocessing from the postprocessing json file

        Parameters
        ----------
        input_path : str
            Path to the input file
        """
        super().__init__(input_path)
        file_options: dict = self.input_data.get('file_options', {})
        modeling_options: dict = self.input_data.get('modeling_options', {})
        data_options: dict = self.input_data['data_options']
        overwrite: dict = file_options.get('overwrite', {})
        self.processed_data_dir: str = file_options.get('processed_data_dir',
                                                        '')
        self.output_dir: str = os.path.join(file_options.get('output_dir', ''),
                                            'images/')
        self.overwrite: bool = overwrite.get('postprocessing', False)
        self.num_groups: int = self.input_data['group_options']['num_groups']
        self.MC_samples: int = self.input_data['group_options']['samples']
        self.irrad_type: str = self.input_data['modeling_options']['irrad_type']
        self.use_data: list[str] = ['keepin', 'charlton', 'endfb6', 'mills']#, 'saleh', 'synetos', 'tuttle', 'waldo']
        self.nuclides: list[str] = ['Br87', 'I137', 'Br88', 'Br89', 'I138', 'Rb94', 'Rb93', 'Te136']
        self.markers: list[str] = ['v', 'o', 'x', '^', 's', 'D']
        self.linestyles: list[str] = ['--', '..', '-.']
        self.load_post_data()
        self.decay_times: np.ndarray[float] = CountRate(input_path).decay_times
        self.num_decay_times = modeling_options['num_decay_times']
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.t_in: float = modeling_options.get('incore_s', 0.0)
        self.t_ex: float = modeling_options.get('excore_s', 0.0)
        self.decay_time_spacing: str = data_options['decay_time_spacing']
        self.total_decay_time: float = modeling_options['decay_time']
        self.group_data = None

        self.MC_yields, self.MC_half_lives = self._get_MC_group_params()

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
        if self.MC_samples > 1:
            self._plot_MC_group_params()
            self._plot_sensitivities(off_nominal=True)
        return None

    def _plot_sensitivities(self, off_nominal: bool = True) -> None:
        """
        Plot the sensitivities of emission probabilities, concentrations,
          and half-lives

        Parameters
        ----------
        off_nominal : bool, optional
            Whether to plot off-nominal sensitivities, by default True
        """
        def scatter_helper(
                data: dict,
                group_params: np.ndarray[float],
                xlab: str,
                ylab: str,
                savename: str,
                savedir: str,
                off_nominal: bool = True) -> None:
            markers = self.markers
            nuclides = self.nuclides or list(data[0].keys())
            xlabel_replace = {
                "Half-life": fr"$T_i [s]$",
                "Decay Constant": fr"$\lambda_i [s^{{-1}}]$",
                "Concentration": fr"$N_i [-]$",
                "Emission Probability": fr'$P_{{n, i}} [-]$'
            }
            ylabel_replace = {
                "Half-life": fr"$T_k [s]$",
                "Decay Constant": fr"$\lambda_k [s^{{-1}}]$",
                "Yield": fr"$\bar{{\nu}}_{{d, k}} [-]$",
            }
            offnom_ylabel_replace = {
                "Half-life": fr"$\Delta T_k [s]$",
                "Decay Constant": fr"$\Delta \lambda_k [s^{{-1}}]$",
                "Yield": fr"$\Delta \bar{{\nu}}_{{d, k}} [-]$",
            }
            xlab = xlabel_replace[xlab]
            if off_nominal:
                ylab = offnom_ylabel_replace[ylab]
            else:
                ylab = ylabel_replace[ylab]
            for nuc in nuclides:
                for group in range(self.num_groups):
                    data_vals = [data[nuc] for data in data]
                    group_vals = group_params[group, 1:]
                    plot_val = group_vals
                    if off_nominal:
                        plot_val = group_vals - np.mean(group_vals)
                    plt.scatter(
                        data_vals,
                        plot_val,
                        label=f'Group {group + 1}',
                        alpha=0.5,
                        s=4,
                        marker=markers[group])
                plt.legend(markerscale=2)
                plt.xlabel(xlab)
                plt.ylabel(ylab)
                plt.savefig(f'{savedir}{savename}_{nuc}.png')
                plt.close()
            return None

        pn_save_dir = os.path.join(self.output_dir, 'sens_pn/')
        if not os.path.exists(pn_save_dir):
            os.makedirs(pn_save_dir)
        lam_save_dir = os.path.join(self.output_dir, 'sens_hl/')
        if not os.path.exists(lam_save_dir):
            os.makedirs(lam_save_dir)
        conc_save_dir = os.path.join(self.output_dir, 'sens_conc/')
        if not os.path.exists(conc_save_dir):
            os.makedirs(conc_save_dir)
        Pn_data = self.post_data['PnMC']
        hl_data = self.post_data['hlMC']
        conc_data = self.post_data['concMC']
        scatter_helper(
            Pn_data,
            self.MC_yields,
            'Emission Probability',
            'Yield',
            'sens_pn_yield',
            pn_save_dir,
            off_nominal=off_nominal)
        scatter_helper(
            hl_data,
            self.MC_yields,
            'Half-life',
            'Yield',
            'sens_lam_yield',
            lam_save_dir,
            off_nominal=off_nominal)
        scatter_helper(
            conc_data,
            self.MC_yields,
            'Concentration',
            'Yield',
            'sens_conc_yield',
            conc_save_dir,
            off_nominal=off_nominal)
        scatter_helper(
            Pn_data,
            self.MC_half_lives,
            'Emission Probability',
            'Half-life',
            'sens_pn_halflife',
            pn_save_dir,
            off_nominal=off_nominal)
        scatter_helper(
            hl_data,
            self.MC_half_lives,
            'Half-life',
            'Half-life',
            'sens_lam_halflife',
            lam_save_dir,
            off_nominal=off_nominal)
        scatter_helper(
            conc_data,
            self.MC_half_lives,
            'Concentration',
            'Half-life',
            'sens_conc_halflife',
            conc_save_dir,
            off_nominal=off_nominal)
        return None

    def compare_yields(self) -> None:
        """
        Compare the total DN yields from summing individuals and from 
          group parameters
        """
        num_top = 7
        num_stack = 2
        summed_yield, summed_avg_halflife = self._get_summed_params(num_top)
        group_yield, group_avg_halflife = self._get_group_params()

        self.summed_yield = summed_yield
        self.summed_avg_halflife = summed_avg_halflife
        self.group_yield = group_yield
        self.group_avg_halflife = group_avg_halflife

        self._plot_nuclide_count_rates(num_stack)
        self.logger.info(f'Summed yield: {summed_yield}')
        self.logger.info(f'Summed average half-life: {summed_avg_halflife} s')
        self.logger.info(f'Group yield {group_yield}')
        self.logger.info(f'Group average half-life: {group_avg_halflife} s')
        return None
    
    def _plot_nuclide_count_rates(self, num_stack: int=1):
        """
        Plot the most important nuclide (by delayed neutron counts emitted) at
            each time step.

        Parameters
        ----------
        num_stack : int, optional
            number of nuclides to plot stacked at each time, by default 1
        """
        data_dict = self._get_data()
        net_nucs = data_dict['net_nucs']
        count_rates = dict()

        for nuc in net_nucs:
            Pn = data_dict['nucs'][nuc]['emission_probability']
            N = data_dict['nucs'][nuc]['concentration']
            hl = data_dict['nucs'][nuc]['half_life']
            lam = np.log(2) / hl
            count_rates[nuc] = list()
            counts = Pn * lam * N * unumpy.exp(-lam * self.decay_times)
            count_rates[nuc] = counts

        biggest_nucs_list = list()
        for ti in range(len(self.decay_times)):
            cur_t_counts = dict()
            for nuc in net_nucs:
                cur_t_counts[nuc] = count_rates[nuc][ti]
            for nuc in range(num_stack):
                try:
                    max_nuc = max(cur_t_counts, key=cur_t_counts.get)
                except ValueError:
                    self.logger.warning("Max nuc evaluation failed")
                    break
                biggest_nucs_list.append(max_nuc)
                del cur_t_counts[max_nuc]
        biggest_nucs = list(dict.fromkeys(biggest_nucs_list))

        for nuci, nuc in enumerate(biggest_nucs):
            rate_n = unumpy.nominal_values(count_rates[nuc])
            rate_s = unumpy.std_devs(count_rates[nuc])
            upper = rate_n + rate_s
            lower = rate_n - rate_s
            plt.fill_between(self.decay_times, lower, upper, color=f'C{nuci}',
                             alpha=0.5)
            plt.plot(self.decay_times, rate_n, color=f'C{nuci}', label=f'{nuc}',
                     linestyle='--', marker=self.markers[nuci%len(self.markers)], markevery=5,
                     markersize=3)
        
        plt.xlabel('Time [s]')
        plt.ylabel(r'Delayed Neutron Rate $[s^{-1}]$')
        plt.xscale('log')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}individual_nuclide_rates.png')
        plt.close()

        for nuci, nuc in enumerate(biggest_nucs):
            rate_n = unumpy.nominal_values(count_rates[nuc])
            rate_s = unumpy.std_devs(count_rates[nuc])
            upper = cumulative_trapezoid(rate_n + rate_s, self.decay_times, initial=0)
            lower = cumulative_trapezoid(rate_n - rate_s, self.decay_times, initial=0)
            rate_n = cumulative_trapezoid(rate_n, self.decay_times, initial=0)
            rate_s = cumulative_trapezoid(rate_s, self.decay_times, initial=0)

            plt.fill_between(self.decay_times, lower, upper, color=f'C{nuci}', alpha=0.5)
            plt.plot(self.decay_times, rate_n, color=f'C{nuci}', label=f'{nuc}',
                     linestyle='--', marker=self.markers[nuci%len(self.markers)], markevery=5,
                     markersize=3)
        plt.xlabel('Time [s]')
        plt.ylabel('Total Delayed Neutron Counts')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}individual_nuclide_counts.png')
        plt.close()
 
        return None

    def _plot_MC_group_params(self) -> None:
        """
        Plot the group parameters from the Monte Carlo NLLS analysis
        """
        def helper_func(name: str,
                        items: np.ndarray[float],
                        group_data: dict[str: object]):
            if name == 'yield':
                label_name = 'Yield'
                xlabel = 'Yield'
                scaler = 1
                scale_label = ''
            elif name == 'half_life':
                label_name = 'Half-life'
                xlabel = 'Half-life [s]'
                scaler = 1
                scale_label = ''
            else:
                raise NotImplementedError(f'{name} not defined')

            group_item = [scaler * ufloat(y,
                                          std) for y,
                          std in zip(group_data[name],
                                     group_data[f'sigma {name}'])]
            for group, item in enumerate(items):
                item = item * scaler
                bins = np.linspace(
                    min(item), max(item), int(
                        np.sqrt(
                            len(item))))
                counts, edges = np.histogram(item, bins=bins)
                normalized_counts = counts
                bin_centers = 0.5 * (edges[:-1] + edges[1:])
                plt.bar(
                    bin_centers,
                    normalized_counts,
                    width=np.diff(edges),
                    label=f'Sampled {label_name}',
                    alpha=0.5,
                    color='red',
                    edgecolor='black')

                plt.axvline(
                    group_item[group].n,
                    color='blue',
                    linestyle='--',
                    label=fr'Group {label_name} ± $1\sigma$')
                plt.axvspan(
                    group_item[group].n -
                    group_item[group].s,
                    group_item[group].n +
                    group_item[group].s,
                    color='blue',
                    alpha=0.25)

                plt.axvline(items[group, 0], color='black',
                            linestyle='-', label=f'Nominal {label_name}')

                plt.xlabel(xlabel + scale_label)
                plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(
                    useMathText=True))
                plt.ticklabel_format(style='sci', axis='x', scilimits=(-2, 2))
                plt.ylabel('Frequency')
                plt.legend()
                plt.tight_layout()
                plt.savefig(f'{self.output_dir}MC_group{group + 1}_{name}.png')
                plt.close()
            return None

        self.group_data = CSVHandler(
            self.group_path,
            create=False).read_vector_csv()

        helper_func('yield', self.MC_yields, self.group_data)
        helper_func('half_life', self.MC_half_lives, self.group_data)

        return None

    def _get_MC_group_params(
            self) -> tuple[np.ndarray[float], np.ndarray[float]]:
        """
        Get the Monte Carlo group parameters from the postprocessing data
        Returns yields and half-lives as numpy arrays

        Returns
        -------
        yields, half_lives : tuple[np.ndarray[float], np.ndarray[float]]
            Tuple containing the yields and half-lives as numpy arrays
        """
        parameters = self.post_data[self.names['groupfitMC']]
        yields = np.zeros((self.num_groups, self.MC_samples))
        half_lives = np.zeros((self.num_groups, self.MC_samples))
        for MC_i, params in enumerate(parameters):
            yield_val = params[:self.num_groups]
            half_life_val = params[self.num_groups:]
            sort_idx = np.argsort(half_life_val)[::-1]
            yields[:, MC_i] = np.asarray(yield_val)[sort_idx]
            half_lives[:, MC_i] = np.asarray(half_life_val)[sort_idx]

        return yields, half_lives

    def _plot_counts(self) -> None:
        """
        Plot the counts from all sources
        """
        sample_color = 'red'
        mean_color = 'black'
        group_color = 'blue'

        counts = self.post_data[self.names['countsMC']]
        countrate = CountRate(self.input_path)
        times = countrate.decay_times
        alpha_MC: float = 1 / np.sqrt(self.MC_samples)
        for MC_iterm, count_val in enumerate(counts):
            label = 'Sampled' if MC_iterm == 0 else None
            plt.plot(
                times,
                count_val,
                alpha=alpha_MC,
                color=sample_color,
                label=label)
        count_data = CSVHandler(self.countrate_path).read_vector_csv()
        plt.errorbar(
            times,
            count_data['counts'],
            count_data['sigma counts'],
            color=mean_color,
            linestyle='',
            marker='x',
            label='Mean',
            markersize=5,
            markevery=5)
        countrate.method = 'groupfit'
        group_counts = countrate.calculate_count_rate(write_data=False)
        plt.plot(
            times,
            group_counts['counts'],
            color=group_color,
            alpha=0.75,
            label='Group Fit',
            linestyle='--',
            zorder=3)
        plt.fill_between(
            times,
            group_counts['counts'] -
            group_counts['sigma counts'],
            group_counts['counts'] +
            group_counts['sigma counts'],
            color=group_color,
            alpha=0.3,
            zorder=2,
            edgecolor='black')
        literature_data = Literature(
            self.input_path).get_group_data(
            self.use_data)
        first: bool = True
        for name, lit_data in literature_data.items():
            if name == 'endfb6':
                name = 'ENDF/B-VI'
            else:
                name = name.capitalize()
            countrate.group_params = lit_data
            data = countrate._count_rate_from_groups()
            plt.plot(times, data['counts'], label=f'{name} 6-Group Fit')
            plt.fill_between(
                times,
                data['counts'] - data['sigma counts'],
                data['counts'] + data['sigma counts'],
                alpha=0.3,
                zorder=2,
                edgecolor='black')
            if first:
                base_name = name
                base_counts = data['counts']
                first = False

        plt.xlabel('Time [s]')
        plt.ylabel(r'Count Rate $[n \cdot s^{-1}]$')
        plt.yscale('log')
        leg = plt.legend()
        for line in leg.legend_handles:
            if line.get_label() == 'Sampled':
                line.set_alpha(0.5)
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}MC_counts.png')
        plt.close()

        for MC_iterm, count_val in enumerate(counts):
            label = 'Sampled' if MC_iterm == 0 else None
            plt.plot(
                times,
                count_val /
                base_counts,
                alpha=alpha_MC,
                color=sample_color,
                label=label,
                linestyle='',
                marker='o',
                markersize=3)
        plt.errorbar(
            times,
            count_data['counts'] /
            base_counts,
            count_data['sigma counts'],
            color=mean_color,
            linestyle='',
            marker='x',
            label='Mean',
            markersize=5,
            markevery=5)

        plt.plot(
            times,
            group_counts['counts'] /
            base_counts,
            color=group_color,
            alpha=0.75,
            label='Group Fit',
            linestyle='--',
            zorder=3)
        plt.fill_between(
            times,
            (group_counts['counts'] - group_counts['sigma counts']) /
            base_counts,
            (group_counts['counts'] + group_counts['sigma counts']) /
            base_counts,
            color=group_color,
            alpha=0.3,
            zorder=2,
            edgecolor='black')
        for name, lit_data in literature_data.items():
            if name == 'endfb6':
                name = 'ENDF/B-VI'
            else:
                name = name.capitalize()
            countrate.group_params = lit_data
            data = countrate._count_rate_from_groups()
            plt.plot(
                times,
                data['counts'] /
                base_counts,
                label=f'{name} 6-Group Fit')
            plt.fill_between(
                times,
                (data['counts'] - data['sigma counts']) / base_counts,
                (data['counts'] + data['sigma counts']) / base_counts,
                alpha=0.3,
                zorder=2,
                edgecolor='black')
        plt.xlabel('Time [s]')
        plt.ylabel(fr'{base_name} Normalized Count Rate')
        leg = plt.legend()
        for line in leg.legend_handles:
            if line.get_label() == 'Sampled':
                line.set_alpha(0.5)
        plt.xscale('log')
        plt.tight_layout()
        plt.savefig(f'{self.output_dir}{base_name}_counts.png')
        plt.close()

        return None

    def _get_group_params(self) -> tuple[float, float]:
        """
        Get the group parameters from the postprocessing data

        returns
        -------
        net_yield, avg_half_life : tuple[float, float]
            net yield and average half-life of the group.

        """
        group_data = CSVHandler(
            self.group_path,
            create=False).read_vector_csv()
        yields = [
            ufloat(
                y, std) for y, std in zip(
                group_data['yield'], group_data['sigma yield'])]
        halflives = [
            ufloat(
                hl, std) for hl, std in zip(
                group_data['half_life'], group_data['sigma half_life'])]
        net_yield = sum(yields)
        if net_yield.n <= 0.0:
            net_yield = ufloat(1e-12, 1e-12)
        lam_vals = np.log(2) / halflives
        # Equation based on Parish 1999 Status of Six-group DN data
        avg_halflife = sum(yields / (net_yield * lam_vals))
        net_yield = ufloat(round(net_yield.n, 5),
                           round(net_yield.s, 5))
        avg_halflife = ufloat(round(avg_halflife.n, 5),
                              round(avg_halflife.s, 5))
        return net_yield, avg_halflife
    
    def _get_data(self) -> dict[str: dict]:
        """
        Collect the data from the processed data files and add them to a
            dictionary

        Returns
        -------
        data_dict : dict[str: dict]
            Dictionary of processed data name to its data
        """
        data_dict = dict()
        emission_prob_data = CSVHandler(
            os.path.join(
                self.processed_data_dir,
                'emission_probability.csv'),
            create=False).read_csv()
        data_dict['emission_probability'] = emission_prob_data
        halflife_data = CSVHandler(
            os.path.join(
                self.processed_data_dir,
                'half_life.csv'),
            create=False).read_csv()
        data_dict['half_life'] = halflife_data
        concentration_data = CSVHandler(
            self.concentration_path,
            create=False).read_csv()
        data_dict['concentration'] = concentration_data

        emission_nucs = list(emission_prob_data.keys())
        conc_nucs = list(concentration_data.keys())
        net_nucs = list(set(emission_nucs) & set(conc_nucs))
        data_dict['net_nucs'] = net_nucs
        data_dict['nucs'] = {}

        for nuc in net_nucs:
            data_dict['nucs'][nuc] = {}
            emission_data = emission_prob_data[nuc]
            Pn = ufloat(emission_data['emission probability'],
                        emission_data['sigma emission probability'])
            conc_data = concentration_data[nuc]
            N = ufloat(conc_data['Concentration'],
                       conc_data['sigma Concentration'])
            hl_data = halflife_data[nuc]
            uncert = hl_data.get('sigma half_life', 1e-12)
            hl = ufloat(hl_data['half_life'], uncert)
            data_dict['nucs'][nuc]['emission_probability'] = Pn
            data_dict['nucs'][nuc]['concentration'] = N
            data_dict['nucs'][nuc]['half_life'] = hl
        return data_dict


    def _get_summed_params(self, num_top: int = 10) -> tuple[float, float]:
        """
        Get the summed parameters from the postprocessing data

        Parameters
        ----------
        num_top : int, optional
            Number of top contributors to consider, by default 10

        returns
        -------
        net_yield, avg_half_life : tuple[float, float]
            net yield and average half-life of the group.
        """
        nuc_yield: dict[str: float] = dict()
        data_dict = self._get_data()
        net_nucs = data_dict['net_nucs']

        halflife_times_yield: dict[str: float] = dict()

        for nuc in net_nucs:
            Pn = data_dict['nucs'][nuc]['emission_probability']
            N = data_dict['nucs'][nuc]['concentration']
            hl = data_dict['nucs'][nuc]['half_life']
            lam_val = np.log(2) / hl
            nuc_yield[nuc] = Pn * N
            halflife_times_yield[nuc] = nuc_yield[nuc] / lam_val

        sorted_yields = dict(
            sorted(
                nuc_yield.items(),
                key=lambda item: item[1].n,
                reverse=True))
        net_yield = sum([i for i in sorted_yields.values()])
        if net_yield.n <= 0.0:
            net_yield = ufloat(1e-12, 1e-12)
        # Parish 1999 uses relative alpha_i values, not yields
        avg_halflife = sum(
            [i / net_yield for i in halflife_times_yield.values()])
        extracted_vals = dict()
        running_sum = 0
        sizes = list()
        labels = list()
        if self.log_level <= INFO:
            counter = 0
            self.logger.info(
                f'Writing nuclide emission times concentration (net yield)')
            for nuc, yield_val in sorted_yields.items():
                self.logger.info(
                    f'{nuc} - {round(yield_val.n, 3)} +/- {round(yield_val.s, 3)}')
                sizes.append(yield_val.n)
                labels.append(nuc)
                running_sum += yield_val
                counter += 1
                extracted_vals[nuc] = yield_val
                if counter > num_top:
                    break
            self.logger.info(
                f'Finished nuclide emission times concentration (net yield)')
            remainder = net_yield.n - running_sum.n
            sizes.append(remainder)
            labels.append('Other')
            colormap = plt.cm.rainbow
            colors = [colormap(i) for i in np.linspace(0.15, 1, num_top + 2)]
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
            wedges, _, _ = ax.pie(sizes, autopct='%1.1f%%',
                             pctdistance=0.7, labeldistance=1.1,
                             colors=colors, textprops={'fontsize': 12})

            ax.legend(
                wedges,
                labels,
                title="Relative Fission Rates",
                loc="center left",
                bbox_to_anchor=(1, 0, 0.5, 1)
            )

            ax.axis('equal')

            plt.tight_layout()
            fig.savefig(f'{self.output_dir}fission_fraction.png')
            plt.close()

        return net_yield, avg_halflife
