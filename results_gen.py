import scale_handler
import ensdf_handler
import keepin_handler
import prk_handler
import misc_funcs
import linear_least_squares as lls
import spectra_handler
import matplotlib.pyplot as plt
import time
import numpy as np
import sys




class RESULTS:
    """
    This class handles generation of results by saving specific settings
        and options directly into itself.
    """

    def __init__(self,
                 times,
                 fissions,
                 efficiency,
                 normalize_value,
                 volume,
                 mass_normalize,
                 ensdf_fname='./ensdf_data/eval_net.xlsx',
                 ensdf_sheet='Sheet1',
                 ORIGEN_out = './scale_outputs/godiva_irrad_post_pulse.out',
                 TRITON_out = './scale_outputs/godiva_3d_depl.out'):
        """
        Initialize. Runs necessary datasets depending on what is
            requested.

        Parameters
        ----------
        times : array-like
            Time values to evaluate at
        fissions : float
            Number of fission events within the sample
        efficiency : float
            Efficiency of the detector for delayed neutrons
        normalize_value : float
            Value ORIGEN multiplies its concentrations by
        volume : float
            Volume of the sample
        mass_normalize : float
            Mass normalization ORIGEN uses for activity   
        ensdf_fname : str
            Path to IAEA .xlsx data
        ensdf_sheet : str
            Name of sheet with data in .xlsx file
        ORIGEN_out : str
            Path to ORIGEN output file
        TRITON_out : str
            Path to TRITON output file

        Returns
        -------
        None
        """
        self.times = times
        self.fissions = fissions
        self.efficiency = efficiency
        self.normalize_value = normalize_value
        self.volume = volume
        self.mass_normalize = mass_normalize
        self.ensdf_fname = ensdf_fname
        self.ensdf_sheet = ensdf_sheet
        self.triout = TRITON_out
        self.oriout = ORIGEN_out
        return


    def delnu_gen(self,
                  scale_output,
                  activity='ensdf',
                  timestep=0,
                  show_iso='all',
                  custom_time=None,
                  errs=True):
        """
        Generate delayed neutron data based on input specifications.

        Parameters
        ----------
        scale_output : str
            Path to SCALE TRITON or ORIGEN output file
        activity : str
            Which method to use for decay constants ('ensdf', 'origen', 'debug')
        timestep : int
            Step at which to evaluate ORIGEN concentrations and activities
        show_iso : str
            Which isotope to show ('all' shows all isotopes) (ex/ xe135)
        custom_time : vector
            A different set of times to feed in for a separate calculation
        

        Returns
        -------
        delnu_counts : array-like
            Detected delay neutron counts from sample at each time
        """
        ensdf_gen = ensdf_handler.ENSDF(self.ensdf_fname,
                                        self.ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        SCALE_gen = scale_handler.SCALE(scale_output, fissions,
                      efficiency, normalize_value, volume,
                       mass_normalize)
        if type(custom_time) != type(None):
            time_use = custom_time
        else:
            time_use = self.times
        delnu_counts, delnu_errs = SCALE_gen.simulate_ensdf_SCALE(time_use, ensdf_dict,
                                                     timestep, show_iso,
                                                     activity, errs=errs)
        n_per_f = misc_funcs.delnu_per_fiss(time_use,
                                            delnu_counts,
                                            self.fissions,
                                            self.efficiency)
        return np.array(delnu_counts), np.array(delnu_errs)

    def compare_delnu_perc_diff(self,
                                counts_1,
                                counts_2,
                                savename='delnu_compare',
                                split=20):
        """
        Compares the delayed neutron results from one decay data
            and another decay data.
        Saves plots.

        Parameters
        ----------
        split : float
            Approximate time to split plot (None means no split)
        counts_1 : array-like
            Counts at each time value
        counts_2 : array-like
            Counts at each time value
        savename : str
            Name of file to save as
            
        Returns
        -------
        None
        """
        split_index = 0
        for ind, t in enumerate(self.times):
            if t < split:
                split_index += 1
            else:
                break
        if split:
            perc_diff = list()
            for ind in range(split_index):
                perc_diff.append((counts_1[ind]-counts_2[ind])/counts_1[ind] * 100)
            plt.plot(self.times[:int(split_index)], perc_diff)
            plt.ylabel('Delayed Neutron Count % Diff')
            plt.xlabel('Time [s]')
            plt.tight_layout()
            plt.savefig(f'{savename}-0-{split}.png')
            plt.close()

            perc_diff = list()
            for ind in range(int(split_index), len(counts_1)):
                perc_diff.append((counts_1[ind]-counts_2[ind])/counts_1[ind] * 100)
            plt.plot(self.times[int(split_index):], perc_diff)
            plt.ylabel('Delayed Neutron Count % Diff')
            plt.xlabel('Time [s]')
            plt.tight_layout()
            plt.savefig(f'{savename}-{split}-end.png')
            plt.close()
            
        perc_diff = list()
        for ind in range(len(counts_1)):
            perc_diff.append((counts_1[ind]-counts_2[ind])/counts_1[ind] * 100)
        plt.plot(self.times, perc_diff)
        plt.ylabel('Delayed Neutron Count % Diff')
        plt.xlabel('Time [s]')
        plt.tight_layout()
        plt.savefig(f'{savename}-full.png')
        plt.close()
        return

    def subtime_delnu_perc_diff(self,
                                time_1,
                                counts_1,
                                time_2,
                                counts_2,
                                savename='delnu_compare',
                                split=20):
        """
        Compares the delayed neutron results from one decay data
            and another decay data with different times.
        Saves plots.

        Parameters
        ----------
        time_1 : vector
            Time vector
        counts_1 : array-like
            Counts at each time value
        time_2 : vector
            Time vector
        counts_2 : array-like
            Counts at each time value
        savename : str
            Name of file to save as
        split : float
            Approximate time to split plot (None means no split)
            
        Returns
        -------
        None
        """
        split_index = 0
        use_time = list()
        use_counts_1 = list()
        use_counts_2 = list()
        for ind1, t1 in enumerate(time_1):
            if round(t1, 5) in np.round(time_2, 5):
                use_counts_1.append(counts_1[ind1])
                use_time.append(t1)
                if t1 < split:
                    split_index += 1
        for ind2, t2 in enumerate(time_2):
            if round(t2, 5) in np.round(time_1, 5):
                use_counts_2.append(counts_2[ind2])
        if split:
            perc_diff = list()
            for ind in range(int(split_index)):
                perc_diff.append((use_counts_1[ind]-use_counts_2[ind])/use_counts_1[ind] * 100)
            plt.plot(use_time[:int(split_index)], perc_diff, marker='.')
            plt.ylabel('Delayed Neutron Count Difference [%]')
            plt.xlabel('Time [s]')
            plt.tight_layout()
            plt.savefig(f'{savename}-0-{split}.png')
            plt.close()

            perc_diff = list()
            for ind in range(int(split_index), len(use_counts_1)):
                perc_diff.append((use_counts_1[ind]-use_counts_2[ind])/use_counts_1[ind] * 100)
            plt.plot(use_time[int(split_index):], perc_diff, marker='.')
            plt.ylabel('Delayed Neutron Count Difference [%]')
            plt.xlabel('Time [s]')
            plt.tight_layout()
            plt.savefig(f'{savename}-{split}-end.png')
            plt.close()

        perc_diff = list()
        for ind in range(len(use_counts_1)):
            perc_diff.append((use_counts_1[ind]-use_counts_2[ind])/use_counts_1[ind] * 100)
        plt.plot(use_time, perc_diff, marker='.')
        plt.ylabel('Delayed Neutron Count Difference [%]')
        plt.xlabel('Time [s]')
        plt.tight_layout()
        plt.savefig(f'{savename}-full-pcnt.png')
        plt.close()

        perc_diff = list()
        for ind in range(len(use_counts_1)):
            perc_diff.append(abs(use_counts_1[ind]-use_counts_2[ind]))#/use_counts_1[ind] * 100)
        plt.plot(use_time, perc_diff, marker='.')
        plt.ylabel('Absolute Delayed Neutron Count Difference [#/s]')
        plt.xlabel('Time [s]')
        plt.yscale('log')
        plt.tight_layout()
        plt.savefig(f'{savename}-full-abs.png')
        plt.close()

        return

    def ILLS_fit(self,
                 counts,
                 numgroups,
                 decay_dims,
                 nodes,
                 count_errs,
                 halflife_base=None,
                 percent_variance=0.3,
                 times=None,
                 ref_counts=None,
                 ref_times=None):
        """
        Use the iterative linear least squares method to create
            groups to fit the given data.

        Parameters
        ----------
        counts : vector
            Delayed neutrons detected at a given time
        numgroups : int
            Number of groups to use
        decay_dims : int
            Number of dimensions to use for decays (1 or 2)
        nodes : int
            Number of nodes to use (should be odd)
        halflife_base : vector
            None if not known, otherwise vector containing the expected
                group half lives
        percent_variance : float
            If using 2d lambda_values and orig values are given,
                this term represents how far away from the given value
                to search
        times : vector
            Used if there is a custom time set that needs to be used
        ref_counts : vector
            Used for refining the abundance results
        ref_times : vector
            Used for refining the abundance results

        Returns
        -------
        abund : vector
            Abundance values for each group
        lami : vector
            Decay constant values for each group
        soln_vec : vector
            a_i*lambda_i, lambda_i for each group
        
        """
        if nodes % 2 == 0:
            print(f'Nodes value {nodes} is not odd') 
            raise Exception
        if decay_dims == 2:
            run_type = '2d'
            max_iters = nodes ** numgroups
        elif decay_dims == 1:
            n = nodes + numgroups - 1
            k = numgroups
            max_iters = math.factorial(n) / (math.factorial(k) * math.factorial(n-k))
            run_type = '1d'
        elif decay_dims == 0:
            max_half = 60
            print(f'Using randomly generated decay constants: 0 - {max_half}s')
            max_iters = 5E4
            print(f'Running {max_iters} iterations')
            run_type = '0d'
        else:
            print('Bad shape {np.shape(np.shape(lambda_values))[0]}')
            raise Exception

        if run_type == '2d' and type(halflife_base) != type(None):
            min_halflives = halflife_base * (1-percent_variance)
            max_halflives = halflife_base * (1+percent_variance)
            halflife_values = np.zeros((numgroups, nodes))
            half_errs = list()
            for group in range(numgroups):
                halflife_values[group, :] = np.linspace(min_halflives[group], max_halflives[group], nodes)
                half_errs.append((max_halflives[group]-min_halflives[group])/(nodes-1))
        elif run_type == '2d':
            print(f'WARNING: Using default min and max values')
            halflife_values = np.zeros((numgroups, nodes))
            min_halflives = [0.01, 0.2, 1.0, 5.0, 10, 40 ]
            max_halflives = [0.20, 1.0, 5.0, 10 , 40, 100]
            half_errs = list()
            print(f'WARNING: Currently set for {len(min_halflives)} groups')
            for group in range(numgroups):
                halflife_values[group, :] = np.linspace(min_halflives[group], max_halflives[group], nodes)
                half_errs.append((max_halflives[group]-min_halflives[group])/(nodes-1))
        elif run_type == '1d':
            print(f'Using default min and max values')
            minhalf = 0.1
            maxhalf = 60
            dhalf = (maxhalf-minhalf)/(num_nodes-1)
            halflife_values = np.arange(minhalf, maxhalf+dhalf*0.1, dhalf)
        elif run_type == '0d':
            halflife_values = max_half
            half_errs = np.sqrt(1 / max_iters) * np.ones(numgroups)
        else:
            print(f'Error in selecting run')
            raise Exception
        lam_errs = list()
        try:
            lambda_values = np.log(2) / halflife_values
            lambda_values.sort()
        except np.AxisError:
            pass
        if type(times) != type(None):
            use_time = times
        else:
            use_time = self.times
        lin_solver = lls.LLS(use_time,
                         counts,
                         self.fissions,
                         self.efficiency,
                         numgroups)
        abund, lami, covariance_mat, residual = lin_solver.lami_solver(lambda_values,
                                                                       max_iters,
                                                                       counts,
                                                                       count_errs,
                                                                       fissions_error)
        if type(abund) == type(None):
            print('ILLS did not converge to a solution')
            raise Exception
##        a_errs = list()
##        for ind in range(numgroups):
##            a_errs.append(residual / (len(self.times) - numgroups) * covariance_mat[ind][ind])
        half_life_final_vals = np.log(2)/lami
        for ind, each in enumerate(half_life_final_vals):
            errval = np.log(2)/each**2 * half_errs[ind]
            lam_errs.append(errval)
        soln_vec = list()
        for each in range(len(abund)):
            soln_vec.append(abund[each] * lami[each])
            soln_vec.append(lami[each])
        #print('Using TEST values '*5)
        #lam_errs = 0 * np.array(lam_errs)
        if group_abundance_err:
            a_errs = lin_solver.MC_abun_err(lami, lam_errs,
                                            counts, count_errs,
                                            tot_iters=abund_iters,
                                            fiserr=fissions_error)
        else:
            print(f'Abundance errors approximated')
        print(f'\nGroup n/F: {np.sum(abund)}')
        print(f'Decay Constants: {lami}')
        print(f'Decay Uncertainties: {lam_errs}')
        print(f'Half Lives: {half_life_final_vals}')
        print(f'Half Life Uncertainties: {half_errs}')
        print(f'Group Yields: {abund.T}\n')
        print(f'Group Yield Uncertainties: {a_errs}')

        print('Formatted Half Lives')
        for ind, half in enumerate(np.log(2) / lami):
            if ind != len(lami) - 1:
                print(f'{np.round(half, 3)} $\pm$ {np.round(half_errs[ind], 3)} & ', end='')
            else:
                print(f'{np.round(half, 3)} $\pm$ {np.round(half_errs[ind], 3)}')
        print('\n')
        print('Formatted Abundancies')
        for ind, abun in enumerate(abund.T*100):
            if ind != len(abund) - 1:
                print(f'{np.round(abun, 3)} $\pm$ {np.round(a_errs[ind]*100, 3)} & ', end='')
            else:
                print(f'{np.round(abun, 3)} $\pm$ {np.round(a_errs[ind]*100, 3)}')
        print('\n')
        print(f'Next iteration half lives:')
        print('np.array([', end='')
        for ind, half in enumerate(np.log(2) / lami):
            if ind != len(lami) - 1:
                print(f'{np.round(half, 5)}', end=', ')
            else:
                print(f'{np.round(half, 5)}', end='])\n')
        if type(ref_counts) != type(None) and type(ref_times) != type(None):
            print(f'Refining abundances')
            abund_new, lami, covariance_mat, residual = lin_solver.abund_refine(lami,
                                                      ref_times,
                                                      ref_counts)
            print(f'Refined yields: {abund_new.T[0]}\n')
            print(f'Refined n/f: {sum(abund_new)}')


        keepin_response = keepin_handler.KEEPIN()
        #delnu_counts, errors = keepin_response.simulate_lin_solve(use_time,
        #                                                  soln_vec,
        #                                                  self.fissions,
        #                                                  self.efficiency,
        #                                                          a_errs,
        #                                                          lam_errs)
        #response = list()
        #for cnt in delnu_counts:
        #    response.append(cnt[0])
        #n_per_f = misc_funcs.delnu_per_fiss(use_time,
        #                                    response,
        #                                    self.fissions,
        #                                    self.efficiency)
        #for index in range(len(delnu_counts)):
        #    avg_pcnt_diff = np.mean(abs(delnu_counts[index] - counts[index]) / delnu_counts[index]) * 100
        #print(f'Average % difference: {avg_pcnt_diff}')

        return abund, lami, soln_vec, a_errs, lam_errs


def group_results(pathname,
                  fit_groups,
                  decay_nodes,
                  results_generator,
                  counts,
                  count_errs,
                  decay_dims,
                  imdir,
                  halflife_base,
                  percent_variation,
                  times,
                  fissions,
                  efficiency,
                  use_errorbars):
    """
    Generate group results for given settings

    Parameters
    ----------
    pathname : str
        Path to where to put outputs
    fit_groups : int
        Number of groups to use in fit
    decay_nodes : int
        Number of nodes to use for generating decay constants
    results_generator : RESULTS class
        Results class
    counts : vector
        Count data to use to generate results
    count_errs : vector
        Associated count errors
    decay_dims : int
        Number of dimensions to use to generate decay matrix or vector
    imdir : str
        Name of default directory to use
    halflife_base : vector
        Default values of half lives used for group fit
    percent_variation : List
        List of percent variation to apply to base (0.1 -> 10%)
    times : vector
        Times to generate group counts at
    fissions : float
        Number of fission events in sample
    efficiency : float
        Efficiency of the neutron detector
    use_errorbars : bool
        Whether to plot errors or not
        

    Returns
    -------
    None

    """
    name = pathname.replace('/', '')
    current_path = imdir+pathname
    misc_funcs.dir_handle(current_path)
    for pcnt in percent_variation:
        print('-'*40)
        numgroups = fit_groups
        nodes = decay_nodes
        print(f'{name} {numgroups} group fit, {nodes} nodes with {pcnt*100}% variance')
        abund, lami, soln_vec, a_errs, lam_errs = results_generator.ILLS_fit(counts,
                                                 numgroups,
                                                 decay_dims,
                                                 nodes,
                                                                             count_errs,
                                                 halflife_base=halflife_base,
                                                 percent_variance=pcnt,
                                                                             times=times)
        keepin_response = keepin_handler.KEEPIN()
        group_counts, group_errs = keepin_response.simulate_lin_solve(times, soln_vec, fissions, efficiency,
                                                                                    a_errs, lam_errs)
        if use_errorbars:
            misc_funcs.multplt(times, group_counts, label=f'{numgroups} group fit', alpha=alpha,
                           errors=group_errs)
            misc_funcs.multplt(times, counts, label=f'Data', alpha=alpha,
                           errors=count_errs)
        else:
            misc_funcs.multplt(times, group_counts, label=f'{numgroups} group fit', alpha=alpha)
            misc_funcs.multplt(times, counts, label=f'{name}', alpha=alpha)

        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'{name}_{numgroups}_group_fit.png')
        plt.close()

        if len(times) == len(results_generator.times) and np.all(np.isclose(times, results_generator.times)):
            results_generator.compare_delnu_perc_diff(counts,
                                                      group_counts,
                                                      savename=current_path+'perc_diff',
                                                      split=5)
        else:
            results_generator.subtime_delnu_perc_diff(times,
                                                        counts,
                                                        times,
                                                      group_counts,
                                                      savename=current_path+'perc_diff',
                                                      split=5)

        prcnt_diff = list()
        for index in range(len(counts)):
            prcnt_diff.append((abs(counts[index] - group_counts[index]) /
                               counts[index]) * 100)
        avg_pcnt_diff = np.mean(prcnt_diff)
        print(f'Average % difference: {avg_pcnt_diff}')
        if not np.all(np.isclose(np.log(2)/lami, halflife_base)):
            print(f'New fit for {pcnt*100}% variation')
            print('-'*50)
            group_results(pathname,
                              fit_groups,
                              decay_nodes,
                              results_generator,
                              counts,
                              count_errs,
                              decay_dims,
                              imdir,
                              np.log(2)/lami,
                              percent_variation,
                              times,
                              fissions,
                              efficiency,
                              use_errorbars)
        else:
            print('Converged for given criteria')

    return

def prefab_group_delnu(times,
                       lam_vec,
                       lam_err,
                       abun_vec,
                       abun_err,
                       fissions,
                       efficiency):
    """
    Useful for when the group data is already known and can be directly modeled

    Parameters
    ----------
    times : vector
        Time values to generate counts at
    lam_vec : vector
        Decay constants for each group
    lam_err : vector
        Errors for each decay constant for each group
    abun_vec : vector
        Group yield (abundance) for each group
    abun_err : vector
        Group yield error for each group
    fissions : float
        Number of fission events in target sample
    efficiency : float
        Efficiency of the neutron detector

    Returns
    -------
    delnu : vector
        Number of counts at each time step
    errors : vector
        Uncertainty in counts at each time step
    """
    delnu = list()
    errors = list()
    err_solve = False
    if type(lam_err) != type(None) and type(abun_err) != type(None):
        err_solve = True
    
    for t in times:
        detect = 0
        err = 0
        for ind in range(len(lam_vec)):
            lami = lam_vec[ind]
            ai = abun_vec[ind]
            if irradiation == 'infinite':
                a_val = ai
            elif irradiation == 'pulse':
                a_val = lami * ai
            detect += (a_val * np.exp(-lami * t))
            if err_solve:
                ai_err = abun_err[ind]
                lami_err = lam_err[ind]
                if irradiation == 'pulse':
                    err += ((lami * np.exp(-lami * t) * ai_err)**2 +
                        (ai*(1-lami*t)*np.exp(-lami*t)*lami_err)**2)
                elif irradiation == 'infinite':
                    err += ((np.exp(-lami * t) * ai_err)**2 +
                        (ai*t*np.exp(-lami*t)*lami_err)**2)
        detect = fissions * efficiency * detect
        delnu.append(detect)
        if err_solve:
            err = np.sqrt(err) * fissions * efficiency
            errors.append(err)
        else:
            errors.append(0)
    return delnu, errors

if __name__ == '__main__':
    from settings import *
    misc_funcs.dir_handle(imdir)

    # Set time to be universal
    time_construct = scale_handler.SCALE(ORIGEN_out,
                 fissions,
                 efficiency,
                 normalize_value)
    times, _ = time_construct.origen_delnu_parser('all')

    # Set energy mesh to be universal
    ORIGEN_result = scale_handler.SCALE(ORIGEN_out,
                                        fissions,
                                        efficiency,
                                        normalize_value,
                                        volume,
                                        mass_normalize)
    _, energy_data, _, _ = ORIGEN_result.origen_spectra_parser()
    energy_mesh = energy_data.copy()
    energy_bin_kev = np.round(np.mean(np.diff(energy_mesh)), 3) * 1E3

    # Run definitions
    results_generator = RESULTS(times,
                                fissions,
                                efficiency,
                                normalize_value,
                                volume,
                                mass_normalize,
                                ensdf_fname=ensdf_fname,
                                ensdf_sheet=ensdf_sheet,
                                ORIGEN_out=ORIGEN_out,
                                TRITON_out=TRITON_out)

    # Determine minimum delayed neutron generation required
        # Concentrations-DecayConstants
        # ORIGEN-ORIGEN, ORIGEN-ENSDF, TRITON-ORIGEN, TRITON-ENSDF
        # Pure ORIGEN
        # Keepin, Brady/England
    pathname = 'delnu-plain/'
    current_path = imdir+pathname
    misc_funcs.dir_handle(current_path)
    # Pure ORIGEN
    if ori_pure_ensdf_comp or \
       ori_pure_group_fit or \
       view_pn_ori_pure or \
       iaea_ori_pure or \
       collect_data or \
       ori_ensdf_counts_cmp:
        print('-'*40)
        print('ORIGEN Pure')
        scale_builder = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        t_ori, pure_delnu_ori = scale_builder.origen_delnu_parser(target)
        n_per_f = misc_funcs.delnu_per_fiss(t_ori, pure_delnu_ori, fissions, efficiency)
        plt.plot(t_ori, pure_delnu_ori)
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.savefig(current_path+f'pure_delnu_ori_{target}.png')
        plt.close()
        
    # ORIGEN-ENSDF
    if run_compare_decay or \
       run_keep_brad_scale or \
       ori_ensdf_group_fit or \
       ori_pure_ensdf_comp or \
       run_ori_tri_compare or \
       ori_ensdf_keep_err or \
       iaea_ori_pure or \
       collect_data or \
       ori_ensdf_counts_cmp:
        origen_ensdf_dn, origen_ensdf_errs = results_generator.delnu_gen(ORIGEN_out,
                                                      activity='ensdf',
                                                      timestep=0,
                                                      show_iso=target)
        plt.plot(times, origen_ensdf_dn)
        origen_ensdf_errs = np.asarray(origen_ensdf_errs)
        if use_errorbars:
            plt.fill_between(times, origen_ensdf_dn+origen_ensdf_errs, origen_ensdf_dn-origen_ensdf_errs, alpha=alpha/2)
        else:
            pass
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.savefig(current_path+f'origen_ensdf_dn_{target}.png')
        plt.close()
    # ORIGEN-ORIGEN
    if run_compare_decay or \
       ori_pure_ensdf_comp or \
       view_pn_ori_pure or \
       collect_data:
        origen_origen_dn, origen_origen_errs = results_generator.delnu_gen(ORIGEN_out,
                                                       activity='origen',
                                                       timestep=0,
                                                       show_iso=target,
                                                       errs=False)
        plt.plot(times, origen_origen_dn)
        if use_errorbars:
            plt.fill_between(times, origen_origen_dn+origen_origen_errs, origen_origen_dn-origen_origen_errs, alpha=alpha/2)
        else:
            pass
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.savefig(current_path+f'origen_origen_dn_{target}.png')
        plt.close()
    # TRITON-ENSDF
    if False:#run_ori_tri_compare or \
       #tri_ensdf_group_fit or \
       #collect_data:
        triton_ensdf_dn, triton_ensdf_errs = results_generator.delnu_gen(TRITON_out,
                                              activity='ensdf',
                                              timestep=0,
                                              show_iso=target)
        plt.plot(times, triton_ensdf_dn)
        if use_errorbars:
            plt.fill_between(times, triton_ensdf_dn+triton_ensdf_errs, triton_ensdf_dn-triton_ensdf_errs, alpha=alpha/2)
        else:
            pass
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.savefig(current_path+f'triton_ensdf_dn_{target}.png')
        plt.close()
    # TRITON-ORIGEN
    if False:
        triton_origen_dn, triton_origen_errs = results_generator.delnu_gen(TRITON_out,
                                              activity='origen',
                                              timestep=0,
                                              show_iso=target)
        plt.plot(times, triton_origen_dn)
        if use_errorbars:
            #plt.errorbar(times, triton_origen_dn, yerr=triton_origen_errs, elinewidth=1)
            plt.fill_between(times, triton_origen_dn+triton_origen_errs, triton_origen_dn-triton_origen_errs, alpha=alpha/2)
        else:
            pass
            #plt.plot(times, triton_origen_dn)
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.savefig(current_path+f'triton_origen_dn_{target}.png')
        plt.close()
    # Keepin
    if run_keep_brad_scale or \
       test_group_fit or \
       ori_ensdf_keep_err or \
       keepin_pure or \
       collect_data or \
       keepin_pure_ori_fit:
        #name = '6keepin235fast'
        #keepin_response = keepin_handler.KEEPIN(name)
        #keepin_delnu, keepin_errs = keepin_response.simulate_instant(times, fissions, efficiency)
        #keepin_delnu = np.array(keepin_delnu)
        #keepin_errs = np.array(keepin_errs)
        keepin_delnu, keepin_errs = prefab_group_delnu(times,
                       keepin_lamvec,
                       keepin_lamerr,
                       keepin_abuvec,
                       keepin_abuerr,
                       fissions,
                       efficiency)
        plt.plot(times, keepin_delnu)
        keepin_delnu = np.asarray(keepin_delnu)
        keepin_errs = np.asarray(keepin_errs)
        if use_errorbars:
            #plt.errorbar(times, keepin_delnu, yerr=keepin_errs, elinewidth=1)
            plt.fill_between(times, keepin_delnu+keepin_errs, keepin_delnu-keepin_errs, alpha=alpha/2)
        else:
            pass
            #plt.plot(times, keepin_delnu)
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.savefig(current_path+'keepin_delnu.png')
        plt.close()

    # Brady/England
    if run_keep_brad_scale or \
       collect_data or \
       keepin_pure_ori_fit:
        #name = '6brengland235fast'
        #brady_england_response = keepin_handler.KEEPIN(name)
        #brady_england_delnu, be_errs = brady_england_response.simulate_instant(times, fissions, efficiency)
        brady_england_delnu, be_errs = prefab_group_delnu(times,
                       be_lamvec,
                       be_lamerr,
                       be_abuvec,
                       be_abuerr,
                       fissions,
                       efficiency)
        plt.plot(times, brady_england_delnu)
        brady_england_delnu = np.array(brady_england_delnu)
        be_errs = np.array(be_errs)
        if use_errorbars:
            #plt.errorbar(times, brady_england_delnu, yerr=be_errs, elinewidth=1)
            plt.fill_between(times, brady_england_delnu+be_errs, brady_england_delnu-be_errs, alpha=alpha/2)
        else:
            pass
            #plt.plot(times, brady_england_delnu)
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.savefig(current_path+'brady_england_delnu.png')
        plt.close()

    # Spectral matrix generation
    if alt_spc_lstsq_oriaea or \
            iaea_ori_2d_spectra or \
            spectra_lstsq_oriaea or \
            spec_compare_oriaea or \
            spectra_puori_fit:
        ORIGEN_result = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        ORIGEN_dict = ORIGEN_result.ensdf_matcher(ensdf_dict,
                                                  0,
                                                  target)
        spectrum_dealer = spectra_handler.SPECTRA(energy_mesh, times)
        #spectral_matrix, valid_list = spectrum_dealer.spectral_matrix_constructor(ORIGEN_dict, ensdf_handler)
        spectral_matrix, valid_list, count_matrix = spectrum_dealer.spectral_matrix_constructor(ORIGEN_dict, ensdf_handler, alt_norm=True)







    # Analysis

    if run_compare_decay:
        # Decay compare IAEA ORIGEN with delayed neutrons
        pathname = 'decay-compare/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        results_generator.compare_delnu_perc_diff(origen_ensdf_dn,
                                                  origen_origen_dn,
                                                  savename=current_path+'perc_diff',
                                                  split=2)
        misc_funcs.multplt(times, origen_origen_dn, label='origen-origen', alpha=alpha)
        misc_funcs.multplt(times, origen_ensdf_dn, label='origen-iaea', alpha=alpha)
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.savefig(current_path+'delnu_tot_compare.png')
        plt.close()

    if run_ori_tri_compare:
        # Compare ORIGEN and TRITON concentration results
        pathname = 'ori-tri-conc/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        results_generator.compare_delnu_perc_diff(origen_ensdf_dn,
                                                  triton_ensdf_dn,
                                                  savename=current_path+'perc_diff',
                                                  split=2)
        misc_funcs.multplt(times, triton_ensdf_dn, label='triton-iaea', alpha=alpha)
        misc_funcs.multplt(times, origen_ensdf_dn, label='origen-iaea', alpha=alpha)
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.savefig(current_path+'delnu_tot_compare.png')
        plt.close()

    if run_keep_brad_scale:
        # Plot SCALE delnu (ORIGEN-IAEA) with groups and data
        pathname = 'ori-iaea-data/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        plt.plot(times, keepin_delnu, label='Keepin')
        
        plt.plot(keepin_response.true_data_time, keepin_response.true_data_resp,
                 label='Keepin True', linestyle='', marker='.')


        misc_funcs.multplt(times, keepin_delnu, label='Brady-England', alpha=alpha)
        misc_funcs.multplt(times, origen_ensdf_dn, label='origen-iaea', alpha=alpha)

        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.savefig(current_path+'delnu_tot_compare.png')
        plt.close()












    if test_group_fit:
        # 6-group fit to Keepin 6-group data
        pathname = f'test-{fit_groups}-fit/'
        group_results(pathname,
                  fit_groups,
                  decay_nodes,
                  results_generator,
                  keepin_delnu,
                  keepin_errs,
                  decay_dims,
                  imdir,
                  halflife_base,
                  percent_variation,
                  times,
                  fissions,
                  efficiency,
                  use_errorbars)


    if ori_ensdf_group_fit:
        # 6-group fit to SCALE delnu (ORIGEN-IAEA)
        pathname = f'origen-iaea-{fit_groups}-fit/'
        group_results(pathname,
                  fit_groups,
                  decay_nodes,
                  results_generator,
                  origen_ensdf_dn,
                  origen_ensdf_errs,
                  decay_dims,
                  imdir,
                  halflife_base,
                  percent_variation,
                  times,
                  fissions,
                  efficiency,
                  use_errorbars)

    if tri_ensdf_group_fit:
        # 6-group fit to SCALE delnu (TRITON-IAEA)
        pathname = f'triton-iaea-{fit_groups}-fit/'
        group_results(pathname,
                  fit_groups,
                  decay_nodes,
                  results_generator,
                  triton_ensdf_dn,
                  triton_ensdf_errs,
                  decay_dims,
                  imdir,
                  halflife_base,
                  percent_variation,
                  times,
                  fissions,
                  efficiency,
                  use_errorbars)

    if ori_pure_group_fit:
        # 6-group fit to SCALE delnu (ORIGEN-Pure)
        pathname = f'origen-pure-{fit_groups}-fit/'
        pure_errs_ori = np.array(pure_delnu_ori) * 0
        group_results(pathname,
                  fit_groups,
                  decay_nodes,
                  results_generator,
                  pure_delnu_ori,
                  pure_errs_ori,
                  decay_dims,
                  imdir,
                  halflife_base,
                  percent_variation,
                  t_ori,
                  fissions,
                  efficiency,
                  use_errorbars)


    if ori_pure_ensdf_comp:
        # Compare ORIGEN delayed neutron output with
        #    ORIGEN using its own decay constants but
        #    ENSDF Pn values
        pathname = 'origen-pure-origen/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        misc_funcs.multplt(t_ori, pure_delnu_ori, label='origen_pure_dn', alpha=alpha)
        misc_funcs.multplt(times, origen_origen_dn, label=f'origen_origen_dn', alpha=alpha)
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.tight_layout()
        plt.savefig(current_path+'Pn-compare.png')
        plt.close()

        results_generator.subtime_delnu_perc_diff(times,
                                                  origen_ensdf_dn,
                                                  t_ori,
                                                  pure_delnu_ori,
                                                  savename=current_path+'perc_diff',
                                                  split=5)
    if ori_ensdf_counts_cmp:
        # Compare ORIGEN delayed neutron output with
        #    ORIGEN using its own decay constants but
        #    ENSDF Pn and lambda values
        pathname = 'compare-origen-pure-origen/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        misc_funcs.multplt(t_ori, pure_delnu_ori, label='Pure ORIGEN', alpha=alpha)
        print(f'Pure ORIGEN delnu 0s: {pure_delnu_ori[0]}')
        misc_funcs.multplt(times, origen_ensdf_dn, label=f'IAEA ORIGEN', alpha=alpha)
        print(f'ORIGEN-IAEA delnu 0s: {origen_ensdf_dn[0]}')
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.tight_layout()
        plt.savefig(current_path+'total-compare.png')
        plt.close()

        results_generator.subtime_delnu_perc_diff(times,
                                                  origen_ensdf_dn,
                                                  t_ori,
                                                  pure_delnu_ori,
                                                  savename=current_path+'perc_diff',
                                                  split=5)

    if ori_ensdf_keep_err:
        # Compare the errors of IAEA in ORIGEN and Keepin 6 group fit
        pathname = 'origen-iaea-keepin-err/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        misc_funcs.multplt(times, keepin_errs, label='keepin_errs', alpha=alpha)
        misc_funcs.multplt(times, origen_ensdf_errs, label=f'origen_ensdf_errs', alpha=alpha)
        plt.yscale('log')
        plt.ylabel('Uncertainty in Counts [#/s]')
        plt.xlabel('Time [s]')
        plt.tight_layout()
        plt.savefig(current_path+'direct-err-compare.png')
        plt.close()



    if test_custom_fit:
        # Custom fit defined here
        pathname = 'custom-fit-compare/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        lam_vec = [0.0132458,  0.03051683, 0.11325934, 0.32377951, 1.34372515, 3.79640257]
        lam_err = [0.00013797703653062053, 0.00029343110454323917, 0.001110385717928915,
                   0.0033727032561812215, 0.012920434119403924, 0.03721963300320383]
        abun_vec = [0.0003788,  0.00274314, 0.00290058, 0.00770585, 0.0026248,  0.00154029]
        abun_err = [1.3046296078626844e-09, 1.4356865542789333e-09, 5.113537184258931e-10,
                    1.7170475634739636e-10, 3.6320708014947495e-11, 2.8076518098563682e-12]
        group_delnu, group_errs = prefab_group_delnu(times,
                       lam_vec,
                       lam_err,
                       abun_vec,
                       abun_err,
                       fissions,
                       efficiency)
        delnu, errs = results_generator.delnu_gen(ORIGEN_out,
                                                      activity='ensdf',
                                                      timestep=0,
                                                      show_iso=target)

        results_generator.compare_delnu_perc_diff(origen_ensdf_dn,
                                                  origen_origen_dn,
                                                  savename=current_path+'perc_diff',
                                                  split=2)
        misc_funcs.multplt(times, group_delnu, label='origen-origen', alpha=alpha, errors=group_errs)
        misc_funcs.multplt(times, delnu, label='origen-ensdf', alpha=alpha, errors=errs)
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.savefig(current_path+'delnu_tot_compare.png')
        plt.close()

    if triton_no_ori_over:
        # Generate results for TRITON data after removing all
        #   isotopes which are also present in ORIGEN
        pathname = 'tri-ori-iso-compare/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        TRITON_gen = scale_handler.SCALE(TRITON_out, fissions, efficiency, normalize_value)
        ORIGEN_gen = scale_handler.SCALE(ORIGEN_out, fissions, efficiency, normalize_value)
        TRITON_IAEA_dict = TRITON_gen.ensdf_matcher(ensdf_dict, 0)
        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        ORIGEN_IAEA_dict = ORIGEN_gen.ensdf_matcher(ensdf_dict, 0)
        isotopes_unique = list()
        for isotope in TRITON_IAEA_dict:
            if isotope in ORIGEN_IAEA_dict:
                pass
            else:
                isotopes_unique.append(isotope)
        contribution = list()
        counts_list_list = list()
        print(f'{len(isotopes_unique)} isotopes to check')
        for isotope in isotopes_unique:
            counts, errs = TRITON_gen.simulate_ensdf_SCALE(times, ensdf_dict, 0, detect_isotope=isotope)
            #misc_funcs.multplt(times, counts, errors=errs, label=f'{isotope}', alpha=alpha)
            counts_list_list.append(counts)
            contribution.append(sum(counts))
        print('Sum value of each isotope')
        for ind, isotope in enumerate(isotopes_unique):
            print(f'{isotope} : {contribution[ind]}')
        print(f'Overall initial count difference: {sum(contribution)}')
        # Re-order lists to be highest contribution first
        orig_contribution = contribution.copy()
        contribution.sort(reverse=True)
        counts_list_list_new = list()
        isotopes_new = list()
        num_contributors = 5
        for ind, cur_cont in enumerate(contribution):
            if ind > num_contributors:
                break
            else:
                move_from = orig_contribution.index(cur_cont)
                isotopes_new.append(isotopes_unique[move_from])
                counts_list_list_new.append(counts_list_list[move_from])
        # Plot top 5
        for ind, scaleval in enumerate(contribution):
            if ind == num_contributors:
                break
            else:
                print(f'Contributor #{ind+1}: {isotopes_new[ind]} : Initial {counts_list_list_new[ind][0]}')
                #misc_funcs.multplt(times, counts_list_list_new[ind], errors=errs, label=f'{isotopes_new[ind]}', alpha=alpha)
        plt.stackplot(times, counts_list_list_new[0:num_contributors], labels=isotopes_new[0:ind])
            
        #plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'{num_contributors}_max_contributors.png')
        plt.close()

    if keepin_pure:
        # Compare Keepins group results with actual data given by Keepin
        #keepin_delnu, keepin_errs, keepin_response
        pathname = 'keepin-keepin/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        results_generator.subtime_delnu_perc_diff(times,
                                                  keepin_delnu,
                                                  keepin_response.true_data_time,
                                                  keepin_response.true_data_resp,
                                                  savename=current_path+'perc_diff',
                                                  split=5)
        misc_funcs.multplt(times, keepin_delnu, label='Keepin Fit', alpha=alpha, errors=keepin_errs)
        plt.plot(keepin_response.true_data_time, keepin_response.true_data_resp,
                           label='Keepin Data', alpha=alpha, linestyle='', marker='.')
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.savefig(current_path+'delnu_tot_compare.png')
        plt.close()
        keepin_delnu, keepin_errs

    if targets_iso_iaea_ori:
        # Designate specific isotopes to investigate
        pathname = 'lam_target_investigate/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        counts_list_list_new = list()
        for cur_target in target_list_iaea_ori:
            diff_counts = list()
            origen_ensdf_dn, origen_ensdf_errs = results_generator.delnu_gen(ORIGEN_out,
                                              activity='ensdf',
                                              timestep=0,
                                              show_iso=cur_target,
                                                                             errs=True)
            origen_origen_dn, origen_ensdf_errs = results_generator.delnu_gen(ORIGEN_out,
                                              activity='origen',
                                              timestep=0,
                                              show_iso=cur_target,
                                                                             errs=False)
            for each in range(len(times)):
                diff_val = (origen_origen_dn[each] - origen_ensdf_dn[each])
                diff_counts.append(diff_val)

            plt.plot(times, diff_counts)
            plt.xlabel('Time [s]')
            plt.ylabel('Pure - IAEA Count Rate [#/s]')
            plt.tight_layout()
            plt.savefig(current_path+f'{cur_target}.png')
            plt.close()

            counts_list_list_new.append(diff_counts)
            print(f'{cur_target} peak: {max(diff_counts)}')
##            plt.plot(times, origen_ensdf_dn, label=f'{target}')
##            if use_errorbars:
##                plt.fill_between(times, origen_ensdf_dn+origen_ensdf_errs, origen_ensdf_dn-origen_ensdf_errs, alpha=alpha/2)
##            else:
##                pass
        label_use = [label_dict[i] for i in target_list_iaea_ori]
        #plt.stackplot(times, np.abs(counts_list_list_new), labels=label_use)
        for i, each in enumerate(counts_list_list_new):
            plt.plot(times, np.abs(each), label=label_use[i])
        plt.yscale('log')
        plt.xscale('log')
        plt.ylim((1, 1e8))
        #plt.ylabel('Pure ORIGEN - IAEA ORIGEN Count Rate [#/s]')
        plt.ylabel('Count Rate Difference [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'origen_ensdf_dn_targets_lam.png')
        plt.close()

    if targets_pure_ori_ori:
        # Designate specific isotopes to investigate
        pathname = 'pn_target_investigate/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        counts_list_list_new = list()
        scale_builder = scale_handler.SCALE(ORIGEN_out,
                                    fissions,
                                    efficiency,
                                    normalize_value,
                                    volume,
                                    mass_normalize)
        for cur_target in target_list_ori_ori:
            diff_counts = list()
            t_ori, pure_delnu_ori = scale_builder.origen_delnu_parser(target)
            origen_origen_dn, origen_ensdf_errs = results_generator.delnu_gen(ORIGEN_out,
                                              activity='origen',
                                              timestep=0,
                                              show_iso=cur_target,
                                                                             errs=False)
            use_time = list()
            use_counts_1 = list()
            use_counts_2 = list()
            for ind1, t1 in enumerate(times):
                if round(t1, 5) in np.round(t_ori, 5):
                    use_counts_1.append(origen_origen_dn[ind1])
                    use_time.append(t1)
            for ind2, t2 in enumerate(t_ori):
                if round(t2, 5) in np.round(times, 5):
                    use_counts_2.append(pure_delnu_ori[ind2])

            #print(f'{target}')
            for each in range(len(use_time)):
                #input(f'Pure / ORI counts: {use_counts_2[each] / use_counts_1[each]}')
                diff_counts.append(use_counts_1[each] - use_counts_2[each])
            counts_list_list_new.append(diff_counts)
            print(f'{cur_target} peak: {max(diff_counts)}')
            #plt.plot(use_time, diff_counts, label=f'{target}')
        label_use = [label_dict[i] for i in target_list_ori_ori]
        #plt.stackplot(use_time, np.abs(counts_list_list_new), labels=label_use)
        for i, each in enumerate(counts_list_list_new):
            plt.plot(times, np.abs(each), label=label_use[i])
        plt.yscale('log')
        plt.xscale('log')
        #plt.yscale('log')
        plt.ylim((1, 1e8))
        plt.ylabel('Count Rate Difference [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'origen_origen_dn_targets_pn.png')
        plt.close()

    if targets_net_pure_ori:
        # Designate specific isotopes to investigate
        pathname = 'net_target_investigate/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        counts_list_list_new_pre = list()
        counts_list_list_new_post = list()
        counts_list_list_new_full = list()

        split = 150
        scale_builder = scale_handler.SCALE(ORIGEN_out,
                                    fissions,
                                    efficiency,
                                    normalize_value,
                                    volume,
                                    mass_normalize)

#        print('Temporary debug measure; line ~1313')
#        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
#        ensdf_dict = ensdf_gen.parse_file()
#
#        ORIGEN_dict = scale_builder.ensdf_matcher(ensdf_dict,
#                                                  0,
#                                                  target)
#        target_list_pure_tot = ORIGEN_dict.keys()
#        # Remove from print to here

        for cur_target in target_list_pure_tot:
            diff_counts = list()
            t_ori, pure_delnu_ori = scale_builder.origen_delnu_parser(cur_target)
            origen_iaea_dn, origen_ensdf_errs = results_generator.delnu_gen(ORIGEN_out,
                                              activity='ensdf',
                                              timestep=0,
                                              show_iso=cur_target,
                                                                             errs=False,
                                                                             custom_time=t_ori)
#            use_time = list()
#            use_counts_1 = list()
#            use_counts_2 = list()
#            for ind1, t1 in enumerate(times):
#                if round(t1, 5) in np.round(t_ori, 5):
#                    use_counts_1.append(origen_iaea_dn[ind1])
#                    use_time.append(t1)
#            for ind2, t2 in enumerate(t_ori):
#                if round(t2, 5) in np.round(times, 5):
#                    use_counts_2.append(pure_delnu_ori[ind2])
            use_time = t_ori
            use_counts_1 = origen_iaea_dn
            use_counts_2 = pure_delnu_ori

#            if len(use_counts_2) < 1:
#                continue

            use_time = np.asarray(use_time)
            split_index = (np.abs(use_time - split)).argmin()
            
            #print(f'{target}')
            for each in range(len(use_time)):
                #input(f'Pure / ORI counts: {use_counts_2[each] / use_counts_1[each]}')
                diff_counts.append(use_counts_2[each] - use_counts_1[each])
#            plt.plot(use_time[:split_index], diff_counts[:split_index])
#            plt.xlabel('Time [s]')
#            plt.ylabel('Pure - IAEA Count Rate [#/s]')
#            plt.tight_layout()
#            plt.savefig(current_path+f'{cur_target}-pre{split}.png')
#            plt.close()
#
#            plt.plot(use_time[split_index:], diff_counts[split_index:])
#            plt.xlabel('Time [s]')
#            plt.ylabel('Pure - IAEA Count Rate [#/s]')
#            plt.tight_layout()
#            plt.savefig(current_path+f'{cur_target}-post{split}.png')
#            plt.close()

            plt.plot(use_time, diff_counts)
            plt.xlabel('Time [s]')
            plt.ylabel('Pure - IAEA Count Rate [#/s]')
            plt.tight_layout()
            plt.savefig(current_path+f'{cur_target}.png')
            plt.close()

            #counts_list_list_new_pre.append(diff_counts[:split_index])
            #counts_list_list_new_post.append(diff_counts[split_index:])
            counts_list_list_new_full.append(diff_counts)
            biggest = max(diff_counts)
            negativest = -min(diff_counts)
            print(f'{cur_target} 0s diff: {diff_counts[0]}')
            #print(f'{cur_target} peak: {max(biggest, negativest)}')
            #plt.plot(use_time, diff_counts, label=f'{target}')

        label_use = [label_dict[i] for i in target_list_pure_tot]
#        plt.stackplot(use_time[:split_index], np.abs(counts_list_list_new_pre), labels=label_use)
#        plt.yscale('log')
#        plt.xscale('log')
#        #plt.yscale('log')
#        #plt.ylim((1e1, 1e7))
#        plt.ylabel('Count Rate Difference [#/s]')
#        plt.xlabel('Time [s]')
#        plt.legend()
#        plt.tight_layout()
#        plt.savefig(current_path+f'origen_iaea_dn_targets-pre{split}.png')
#        plt.close()
#
#        plt.stackplot(use_time[split_index:], np.abs(counts_list_list_new_post), labels=label_use)
#        plt.yscale('log')
#        plt.xscale('log')
#        #plt.yscale('log')
#        #plt.ylim((1e1, 1e7))
#        plt.ylabel('Count Rate Difference [#/s]')
#        plt.xlabel('Time [s]')
#        plt.legend()
#        plt.tight_layout()
#        plt.savefig(current_path+f'origen_iaea_dn_targets-post{split}.png')
#        plt.close()

        #plt.stackplot(use_time, np.abs(counts_list_list_new_full), labels=label_use)
        for i, each in enumerate(counts_list_list_new_full):
            plt.plot(times, np.abs(each), label=label_use[i])
        plt.yscale('log')
        plt.xscale('log')
        plt.ylim((1, 1e8))
        #plt.yscale('log')
        #plt.ylim((1e1, 1e7))
        plt.ylabel('Count Rate Difference [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'origen_iaea_dn_targets_net.png')
        plt.close()


        t_ori, pure_delnu_ori = scale_builder.origen_delnu_parser(target)
        origen_iaea_dn, origen_ensdf_errs = results_generator.delnu_gen(ORIGEN_out,
                                          activity='ensdf',
                                          timestep=0,
                                          show_iso=target,
                                                                         errs=False)
        plt.plot(t_ori, pure_delnu_ori, label='Pure ORIGEN')
        plt.plot(times, origen_iaea_dn, label='ORIGEN-IAEA')
        plt.ylabel('Delayed Neutrons [#/s]')
        plt.xlabel('Time [s]')
        plt.yscale('log')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'net_counts.png')
        plt.close()

    if find_worst_lam_isos:
        origen_origen_dn, origen_origen_errs = results_generator.delnu_gen(ORIGEN_out,
                                                       activity='lamdebug',
                                                       timestep=0,
                                                       show_iso=target,
                                                       errs=False)
        
    if find_worst_pn_isos:
        origen_origen_dn, origen_origen_errs = results_generator.delnu_gen(ORIGEN_out,
                                                       activity='pndebug',
                                                       timestep=0,
                                                       show_iso=target,
                                                       errs=False)
    if find_worst_tot_isos:
        origen_origen_dn, origen_origen_errs = results_generator.delnu_gen(ORIGEN_out,
                                                       activity='puredebug',
                                                       timestep=0,
                                                       show_iso=target,
                                                       errs=False)

    if view_pn_ori_pure:
        # Designate specific isotopes to investigate
        pathname = 'pn-generic-plots/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        plt.plot(times, origen_origen_dn)
        if use_errorbars:
            plt.fill_between(times, origen_origen_dn+origen_origen_errs, origen_origen_dn-origen_origen_errs, alpha=alpha/2)
        else:
            pass
        misc_funcs.multplt(times, origen_origen_dn, errors=origen_origen_errs, alpha=alpha,
                           label='origen_origen_dn')
        misc_funcs.multplt(t_ori, pure_delnu_ori, alpha=alpha,
                           label='pure_origen_dn')
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.savefig(current_path+f'pure_delnu_ori_{target}.png')
        plt.close()

    if iaea_ori_pure:
        # Compare Pure ORIGEN and ORIGEN-IAEA
        # t_ori, pure_delnu_ori
        # origen_ensdf_dn, origen_ensdf_errs
        pathname = 'pure-iaea-compare/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        misc_funcs.multplt(times, origen_ensdf_dn, label='IAEA-ORIGEN',
                           alpha=alpha)
        misc_funcs.multplt(t_ori, pure_delnu_ori, label='Pure ORIGEN',
                           alpha=alpha)
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.savefig(current_path+f'pure_iaea_ori_{target}.png')
        plt.close()
        results_generator.subtime_delnu_perc_diff(times,
                                                  origen_ensdf_dn,
                                                  t_ori,
                                                  pure_delnu_ori,
                                                  savename=current_path+'perc_diff',
                                                  split=5)

    if ori_iaea_keepin_prk:
        # Compare the point reactor kinetics response for
        #   ORIGEN-IAEA and Keepin group fits
        pathname = 'prk-ori-iaea-keepin/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        
        for reactivity_magnitude in reactivity_magnitudes:
            group_solve = prk_handler.GPRK(ori_iaea_lamvec, ori_iaea_abuvec,
                                           prk_times, reactivity_magnitude,
                                           gen_time=gen_time,
                                           lambda_errs=ori_iaea_lamerr,
                                           abundance_errs=ori_iaea_abuerr)
            soln_matrix, err_matrix = group_solve.group_prk()
            misc_funcs.multplt(prk_times, soln_matrix[0, :], errors=err_matrix[0, :],
                               label=f'ORIGEN-IAEA Step Insertion of {reactivity_magnitude}$',
                               alpha=alpha)

            group_solve = prk_handler.GPRK(keepin_lamvec, keepin_abuvec,
                                           prk_times, reactivity_magnitude,
                                           gen_time=gen_time,
                                           lambda_errs=keepin_lamerr,
                                           abundance_errs=keepin_abuerr)
            soln_matrix, err_matrix = group_solve.group_prk()
            misc_funcs.multplt(prk_times, soln_matrix[0, :], errors=err_matrix[0, :],
                               label=f'Keepin [6] step insertion of {reactivity_magnitude}$',
                               alpha=alpha)
        
        plt.yscale('log')
        #plt.yticks([0.1, 1, 10])
        plt.xlabel('Time [s]')
        plt.ylabel('Realtive Neutron Density')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'response.png')
        plt.close()

    if puori_iaea_ori_prk:
        # Compare the point reactor kinetics response for
        #   ORIGEN-IAEA and Pure ORIGEN group fits
        pathname = 'prk-ori-iaea-pure/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        for reactivity_magnitude in reactivity_magnitudes:
            group_solve = prk_handler.GPRK(ori_iaea_lamvec, ori_iaea_abuvec,
                                           prk_times, reactivity_magnitude,
                                           gen_time=gen_time,
                                           lambda_errs=ori_iaea_lamerr,
                                           abundance_errs=ori_iaea_abuerr)
            soln_matrix, err_matrix = group_solve.group_prk()
            misc_funcs.multplt(prk_times, soln_matrix[0, :], errors=err_matrix[0, :],
                               label=f'ORIGEN-IAEA step insertion of {reactivity_magnitude}$',
                               alpha=alpha)

            group_solve = prk_handler.GPRK(pure_ori_lamvec, pure_ori_abuvec,
                                           prk_times, reactivity_magnitude,
                                           gen_time=gen_time,
                                           lambda_errs=pure_ori_lamerr,
                                           abundance_errs=pure_ori_abuerr)
            soln_matrix, err_matrix = group_solve.group_prk()
            misc_funcs.multplt(prk_times, soln_matrix[0, :], errors=err_matrix[0, :],
                               label=f'Pure ORIGEN step insertion of {reactivity_magnitude}$',
                               alpha=alpha)
        
        
        plt.yscale('log')
        #plt.yticks([0.1, 1, 10])
        plt.xlabel('Time [s]')
        plt.ylabel('Relative Neutron Density')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'response.png')
        plt.close()

    if keepin_pure_iaea_prk:
        # Compare the point reactor kinetics response for
        #   ORIGEN-IAEA, Keepin, and Pure ORIGEN group fits
        pathname = 'prk-keep-ori-iaea-pure/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        for reactivity_magnitude in reactivity_magnitudes:
            group_solve = prk_handler.GPRK(ori_iaea_lamvec, ori_iaea_abuvec,
                                           prk_times, reactivity_magnitude,
                                           gen_time=gen_time,
                                           lambda_errs=ori_iaea_lamerr,
                                           abundance_errs=ori_iaea_abuerr)
            soln_matrix, err_matrix = group_solve.group_prk()
            misc_funcs.multplt(prk_times, soln_matrix[0, :], errors=err_matrix[0, :],
                               label=f'ORIGEN-IAEA step insertion of {reactivity_magnitude}$',
                               alpha=alpha)

            group_solve = prk_handler.GPRK(pure_ori_lamvec, pure_ori_abuvec,
                                           prk_times, reactivity_magnitude,
                                           gen_time=gen_time,
                                           lambda_errs=pure_ori_lamerr,
                                           abundance_errs=pure_ori_abuerr)
            soln_matrix, err_matrix = group_solve.group_prk()
            misc_funcs.multplt(prk_times, soln_matrix[0, :], errors=err_matrix[0, :],
                               label=f'Pure ORIGEN step insertion of {reactivity_magnitude}$',
                               alpha=alpha)

            group_solve = prk_handler.GPRK(keepin_lamvec, keepin_abuvec,
                                           prk_times, reactivity_magnitude,
                                           gen_time=gen_time,
                                           lambda_errs=keepin_lamerr,
                                           abundance_errs=keepin_abuerr)
            soln_matrix, err_matrix = group_solve.group_prk()
            misc_funcs.multplt(prk_times, soln_matrix[0, :], errors=err_matrix[0, :],
                               label=f'Keepin [6] step insertion of {reactivity_magnitude}$',
                               alpha=alpha)
        
        plt.yscale('log')
        #plt.yticks([0.1, 1, 10])
        plt.xlabel('Time [s]')
        plt.ylabel('Relative Neutron Density')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'response.png')
        plt.close()


    if important_pure_ori:
        # Print the important isotopes at each time step from Pure ORIGEN
        origen_origen_dn, origen_origen_errs = results_generator.delnu_gen(ORIGEN_out,
                                                       activity='purecheck',
                                                       timestep=0,
                                                       show_iso=target,
                                                       errs=False)
    if keepin_pure_ori_fit:
        # Manually entered best fit
        pathname = 'fits-compare/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        group_counts, group_errs = prefab_group_delnu(times,
                                                      pure_ori_lamvec,
                                                      pure_ori_lamerr,
                                                      pure_ori_abuvec,
                                                      pure_ori_abuerr,
                                                      fissions,
                                                      efficiency)
        misc_funcs.multplt(times, group_counts, errors=group_errs, alpha=alpha,
                           label='Pure ORIGEN 6 Group Fit')

        #group_counts, group_errs = prefab_group_delnu(times,
        #                                              ori_iaea_lamvec,
        #                                              ori_iaea_lamerr,
        #                                              ori_iaea_abuvec,
        #                                              ori_iaea_abuerr,
        #                                              fissions,
        #                                              efficiency)
        #misc_funcs.multplt(times, group_counts, errors=group_errs, alpha=alpha,
        #                   label='ORIGEN-IAEA 6 Group Fit')
        misc_funcs.multplt(times, np.array(keepin_delnu), errors=np.array(keepin_errs),
                           alpha=alpha,
                           label='Keepin 6 Group Fit')
        misc_funcs.multplt(times, brady_england_delnu, errors=be_errs, alpha=alpha,
                           label='Brady/England 6 Group Fit')
        plt.xlabel('Time [s]')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.tight_layout()
        plt.savefig(current_path+f'fits.png')
        plt.close()


        # Difference between Pure fit and (Keepin and Brady/England)

        results_generator.subtime_delnu_perc_diff(times,
                                                  group_counts,
                                                  times,
                                                  np.asarray(keepin_delnu),
                                                  savename=current_path+'perc-diff-Pure-Keepin',
                                                  split=5)
        
        results_generator.subtime_delnu_perc_diff(times,
                                                  group_counts,
                                                  times,
                                                  np.asarray(brady_england_delnu),
                                                  savename=current_path+'perc-diff-Pure-BE',
                                                  split=5)

        # End Difference

    if collect_data:
        # Continue generating useful data
        print('Pure ORIGEN 6-group fit')
        group_counts, group_errs = prefab_group_delnu(times,
                                                      pure_ori_lamvec,
                                                      pure_ori_lamerr,
                                                      pure_ori_abuvec,
                                                      pure_ori_abuerr,
                                                      fissions,
                                                      efficiency)
        misc_funcs.delnu_per_fiss(times, group_counts, fissions, efficiency)
        print('ORIGEN-IAEA 6-group fit')
        group_counts, group_errs = prefab_group_delnu(times,
                                                      ori_iaea_lamvec,
                                                      ori_iaea_lamerr,
                                                      ori_iaea_abuvec,
                                                      ori_iaea_abuerr,
                                                      fissions,
                                                      efficiency)
        misc_funcs.delnu_per_fiss(times, group_counts, fissions, efficiency)

    if pure_ori_t_spectra:
        # Generate spectra for pure ORIGEN output
        pathname = 'spectra-t-pure-ori/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        scale_builder = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        

        time_data, energy_data, spectra_matrix, bin_data = scale_builder.origen_spectra_parser()
        print(f'Pure ORIGEN Spectrum Generation for {len(time_data)} times')
        # Spectra for different times
        n = 1
        cur_begin_time = time.time()
        for tind, t in enumerate(time_data):
            if tind/len(time_data) >= 0.1 * n:
                n += 1
                net_time = time.time() - cur_begin_time
                full_complete_time = net_time / (0.1 * (n-1))
                print(f'    Progress: {round(tind/len(time_data) * 100)}% in {round(net_time, 0)}s')
                print(f'        Estimated completion in {round(full_complete_time - net_time, 0)}s')
            #norm_factor = spectra_normalize / np.sum(spectra_matrix[:, tind])
            if spectra_normalized:
                normalize = np.sum(spectra_matrix[:, tind])
                name = 'Probability'
            else:
                normalize = 1/efficiency
                name = 'Neutron Intensity'
            y = spectra_matrix[:, tind] / normalize
            plt.step(energy_data, y, where='mid')
            plt.title(f'{t}s spectra')
            plt.xlabel('Energy [MeV]')
            plt.ylabel(f'{name} / {energy_bin_kev} keV')
            plt.savefig(current_path+f'{tind}-{name}.png')
            plt.close()

    if pure_ori_2d_spectra:
        # Generate 2D spectra for pure ORIGEN output
        pathname = 'spectra-2d-pure-ori/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        scale_builder = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)

        time_data, energy_data, spectra_matrix, bin_data = scale_builder.origen_spectra_parser()
        
        # Heatmap - 10^4 counts
        fig, ax = plt.subplots()
        x, y = np.meshgrid(time_data, energy_data)
        # Column is energy spectra for given energy (i.e. time index)
        z = np.zeros(np.shape(spectra_matrix))
        for tind, t in enumerate(time_data):
            if spectra_normalized:
                normalize = np.sum(spectra_matrix[:, tind])
                name = 'Probability'
            else:
                normalize = 1
                name = 'Neutron Intensity'
            #norm_factor = spectra_normalize / np.sum(spectra_matrix[:, tind])
            z[:, tind] = spectra_matrix[:, tind] / normalize #norm_factor * spectra_matrix[:, tind]
        c = ax.pcolormesh(x, y, z, cmap='magma')
        cbar = fig.colorbar(c, ax=ax)
        cbar.set_label(name)
        plt.xlabel('Time [s]')
        plt.ylabel('Energy [MeV]')
        plt.tight_layout()
        plt.savefig(current_path+f'{name}.png')
        plt.close()

    if iaea_ori_t_spectra:
        # Generate time dependent spectra using IAEA spectra data and ORIGEN concentrations
        pathname = 'spectra-t-iaea-ori/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        print(f'IAEA-ORIGEN Spectrum Generation for {len(times)} times')
        # Concentrations 
        ORIGEN_result = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        ORIGEN_dict = ORIGEN_result.ensdf_matcher(ensdf_dict,
                                                  0,
                                                  target)
        # Spectra
        valid_list = list()
        for isotope in ORIGEN_dict.keys():
            if isotope in iaea_spectra:
                bins, values = ensdf_handler.spectra_analysis(iaea_spectra[isotope],
                                                              path='spectra/spectra/',
                                                              display=False)
                ORIGEN_dict[isotope]['spectrum_bins'] = bins
                ORIGEN_dict[isotope]['spectrum_values'] = values
                valid_list.append(isotope)


        # Now iterate through time
        max_energy = 2
        n = 1
        impactful_iso = ''
        prev_impactful_iso = None
        cur_begin_time = time.time()
        for tind, t in enumerate(times):
            if tind/len(times) >= 0.1 * n:
                n += 1
                net_time = time.time() - cur_begin_time
                full_complete_time = net_time / (0.1 * (n-1))
                print(f'    Progress: {round(tind/len(times) * 100)}% in {round(net_time, 0)}s')
                print(f'        Estimated completion in {round(full_complete_time - net_time, 0)}s')

            activity_err = 0
            impactful_activity = 0
            net_spectra = np.zeros(len(energy_mesh))
            net_spec_err = np.zeros(len(energy_mesh))
            # sum over each isotope
            for isotope in valid_list:
                atoms = ORIGEN_dict[isotope]['conc'][0]
                Pn = ORIGEN_dict[isotope]['emission'][0]
                lam = np.log(2) / ORIGEN_dict[isotope]['halflife'][0]
                activity = (efficiency * ORIGEN_dict[isotope]['emission'][0] *
                            atoms[tind] * lam)
                if activity > impactful_activity:
                    impactful_iso = isotope
                    impactful_activity = activity
                Pn_err = ORIGEN_dict[isotope]['emission'][1]
                lam_err = (np.log(2) / ORIGEN_dict[isotope]['halflife'][0]**2 *
                           ORIGEN_dict[isotope]['halflife'][1])
                atom_err = ORIGEN_dict[isotope]['conc'][1]

                activity_err += np.sqrt(((lam * atoms[tind] * np.exp(-lam * t) * Pn_err)**2 +
                            (Pn * atoms[tind] * (1-lam*t) * np.exp(-lam*t) * lam_err)**2 +
                            (Pn * lam * np.exp(-lam*t) * atom_err)**2))
                for eind, e in enumerate(energy_mesh):
                    # spectra uses closest datapoint
                    closest_bindex = (np.abs(ORIGEN_dict[isotope]['spectrum_bins'] -
                                          e)).argmin()
                    net_spectra[eind] += (activity *
                                          ORIGEN_dict[isotope]['spectrum_values'][closest_bindex])
                    cur_spectra_err = 0
                    err_val = (efficiency * activity_err / max_energy) ** 2 + cur_spectra_err**2
                    net_spec_err[eind] += err_val
            net_spec_err = np.sqrt(net_spec_err)
            if spectra_normalized:
                name = 'Probability'
                normalize = np.sum(net_spectra)
                # new error calc
                use_err = np.zeros(len(energy_mesh))
                sum_err = np.sqrt(np.sum(net_spec_err**2))
                for eind, enerr in enumerate(net_spec_err):
                    use_err[eind] = np.sqrt( (enerr / normalize)**2 + (-net_spectra[eind]/normalize**2 * sum_err)**2 )
            else:
                normalize = 1
                name = 'Neutron Intensity'
                use_err = net_spec_err
            if impactful_iso != prev_impactful_iso:
                print(f'{round(t, 2)}s : {impactful_iso} most impactful with {impactful_activity} +/- {net_spec_err[0]} counts')
                prev_impactful_iso = impactful_iso
            y = net_spectra/normalize
            plt.step(energy_mesh, y, where='mid')
            plt.fill_between(energy_mesh, y+use_err,
                             y-use_err, alpha=alpha/2,
                             step='mid')
            plt.title(f'{t}s')
            plt.xlabel('Energy [MeV]') 
            plt.ylabel(f'{name} / {energy_bin_kev} keV')
            plt.savefig(current_path+f'{t}s-{name}.png')
            plt.close()
        
    if iaea_ori_2d_spectra:
        # Generate 2D spectra for ORIGEN-IAEA data
        pathname = 'spectra-2d-ori-iaea/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        # Concentrations 
        ORIGEN_result = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        ORIGEN_dict = ORIGEN_result.ensdf_matcher(ensdf_dict,
                                                  0,
                                                  target)
        # Spectra
        ORIGEN_dict, valid_list = spectra_handler.SPECTRA(energy_mesh, times).update_dict(ORIGEN_dict, ensdf_handler)


        # Now iterate through time

        fig, ax = plt.subplots()
        x, y = np.meshgrid(times, energy_mesh)
        # Column is energy spectra for given energy
#        z = np.zeros((len(energy_mesh), len(times)))
#        for tind, t in enumerate(times):
#            net_spectra = np.zeros(len(energy_mesh))
#            # sum over each isotope
#            for isotope in valid_list:
#                atoms = ORIGEN_dict[isotope]['conc'][0]
#                Pn = ORIGEN_dict[isotope]['emission'][0]
#                lam = np.log(2) / ORIGEN_dict[isotope]['halflife'][0]
#                activity = (efficiency * ORIGEN_dict[isotope]['emission'][0] *
#                            atoms * lam * np.exp(-lam * t))
#                for eind, e in enumerate(energy_mesh):
#                    # spectra uses closest datapoint
#                    closest_bindex = (np.abs(ORIGEN_dict[isotope]['spectrum_bins'] -
#                                          e)).argmin()
#                    net_spectra[eind] += (activity *
#                                          ORIGEN_dict[isotope]['spectrum_values'][closest_bindex])
#
#
#            if spectra_normalized:
#                normalize = np.sum(net_spectra)
#                name = 'Probability'
#            else:
#                normalize = 1
#                name = 'Neutron Intensity'
#            z[:, tind] = net_spectra / normalize
        z = spectral_matrix
        c = ax.pcolormesh(x, y, z, cmap='magma')
        cbar = fig.colorbar(c, ax=ax)
        cbar.set_label(name)
        plt.xlabel('Time [s]')
        plt.ylabel('Energy [MeV]')
        plt.tight_layout()
        plt.savefig(current_path+f'{name}.png')
        plt.close()

    if spectra_puoriaea_com:
        # Compare the spectra results of Pure ORIGEN output and
        #   ORIGEN-IAEA output.
        # 1D and 2D
        pathname = 'spectra-puoriaea-compare/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        # Need data for both but only for matching time steps
        # Concentrations 
        ORIGEN_result = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        ORIGEN_dict = ORIGEN_result.ensdf_matcher(ensdf_dict,
                                                  0,
                                                  target)
        # Spectra
        ORIGEN_dict, valid_list = spectra_handler.SPECTRA(energy_mesh, times).update_dict(ORIGEN_dict, ensdf_handler)

#        print('REMOVE')
#        net_activity_per_t = 0
#        for isotope in ORIGEN_dict.keys():
#            t = 0
#            atoms = ORIGEN_dict[isotope]['conc'][0]
#            Pn = ORIGEN_dict[isotope]['emission'][0]
#            lam = np.log(2) / ORIGEN_dict[isotope]['halflife'][0]
#            activity = (efficiency * Pn * atoms * lam * np.exp(-lam * t))
#            net_activity_per_t += activity
#        print(f'Total ORIGEN-IAEA activity at {t}s: {net_activity_per_t}')
#        print('END REMOVE')
        

        time_data, energy_data, spectra_matrix, bin_data = ORIGEN_result.origen_spectra_parser()

        ORIGEN_dict, valid_list = spectra_handler.SPECTRA(energy_data, times).update_dict(ORIGEN_dict, ensdf_handler)

        basex, basey = np.meshgrid(time_data, energy_data)
        spec_puorigen = np.zeros(np.shape(spectra_matrix))
        spec_oriaea   = np.zeros(np.shape(spectra_matrix))
        print(f'Pure ORIGEN and ORIGEN-IAEA Spectra Generation for {len(time_data)} times')
        print('Activity error calculation still uses exponential decay, use Sampler uncertainty over time')
        # Spectra for different times
        n = 1
        max_energy = max(energy_data)
        cur_begin_time = time.time()
        l2_norms = list()
        diff_avgs = list()
        for tind, t in enumerate(time_data):
            # ORIGEN_IAEA Section
            activity_err = 0
            net_spectra = np.zeros(len(energy_data))
            net_spec_err = np.zeros(len(energy_data))
            net_activity_per_t = 0
            # sum over each isotope
            #print(f'{t}s')
            for isotope in valid_list:
                atoms = ORIGEN_dict[isotope]['conc'][0]
                Pn = ORIGEN_dict[isotope]['emission'][0]
                lam = np.log(2) / ORIGEN_dict[isotope]['halflife'][0]
                activity = (efficiency * Pn * atoms[tind] * lam)
                net_activity_per_t += activity
                Pn_err = ORIGEN_dict[isotope]['emission'][1]
                lam_err = (np.log(2) / ORIGEN_dict[isotope]['halflife'][0]**2 *
                           ORIGEN_dict[isotope]['halflife'][1])
                atom_err = ORIGEN_dict[isotope]['conc'][1]
                activity_err += np.sqrt(((lam * atoms[0] * np.exp(-lam * t) * Pn_err)**2 +
                            (Pn * atoms[0] * (1-lam*t) * np.exp(-lam*t) * lam_err)**2 +
                            (Pn * lam * np.exp(-lam*t) * atom_err)**2))
                #print(isotope)
                #print(activity_err)
                #print(activity)
                #print(activity_err/activity)
                prev_bindex = 0
                bin_indeces = np.digitize(ORIGEN_dict[isotope]['spectrum_bins'], energy_data)
                for eind, e in enumerate(energy_data):
                    if activity < activity_err:
                        err_frac = 0
                    else:
                        err_frac = activity_err / activity
                    net_spectra[eind] += activity * ORIGEN_dict[isotope]['spectrum_values'][eind]
                    cur_spectra_err = 0
                    # Constant error to each bin
                    #net_spec_err[eind] += (efficiency * activity_err / max_energy) ** 2 + cur_spectra_err**2
                    net_spec_err[eind] += (efficiency * activity_err / len(energy_data)) ** 2 + cur_spectra_err**2
                    # Constant ratio of error to each bin
                    #net_spec_err[eind] += (err_frac * ORIGEN_dict[isotope]['spectrum_values'][eind]) ** 2 + cur_spectra_err**2 
            net_spec_err = np.sqrt(net_spec_err)
            if spectra_normalized:
                name = 'Probability'
                normalize = np.sum(net_spectra)
                # new error calc
                use_err = np.zeros(len(energy_data))
                sum_err = np.sqrt(np.sum(net_spec_err**2))
                for eind, enerr in enumerate(net_spec_err):
                    use_err[eind] = np.sqrt( (enerr / normalize)**2 + (-net_spectra[eind]/normalize**2 * sum_err)**2 )
            else:
                normalize = 1
                name = 'Neutron Intensity'
                use_err = net_spec_err
                #print(f'ORIGEN-IAEA net counts at {t}s: {np.sum(net_spectra)}')
            y = net_spectra/normalize
            spec_oriaea[:, tind] = net_spectra/normalize
            plt.step(energy_data, y, where='mid', label='IAEA-ORIGEN')
            if spectra_uncertainty:
                plt.fill_between(energy_data, y+use_err,
                                 y-use_err, alpha=alpha/2,
                                 step='mid')


            # Pure ORIGEN Section
            if tind/len(time_data) >= 0.1 * n:
                n += 1
                net_time = time.time() - cur_begin_time
                full_complete_time = net_time / (0.1 * (n-1))
                #print(f'    Progress: {round(tind/len(time_data) * 100)}% in {round(net_time, 0)}s')
                #print(f'        Estimated completion in {round(full_complete_time - net_time, 0)}s')
            #norm_factor = spectra_normalize / np.sum(spectra_matrix[:, tind])
            if spectra_normalized:
                normalize = np.sum(spectra_matrix[:, tind])
                name = 'Probability / 2.0 keV'
            else:
                normalize = mass_normalize/(efficiency)
                name = 'Neutron Intensity'
                #input(f'Pure ORIGEN net counts at {t}s: {np.sum(spectra_matrix[:, tind] / normalize)}\n')
            y = spectra_matrix[:, tind] / normalize
            spec_puorigen[:, tind] = spectra_matrix[:, tind] / normalize
            plt.step(energy_data, y, where='mid', label='Pure ORIGEN')
            plt.xlabel('Energy [MeV]')
            plt.ylabel(name)
            plt.legend()
            plt.tight_layout()
            plt.savefig(current_path+f'{tind}.png')
            plt.close()

            cur_pcnt_diff = abs(spec_oriaea[:, tind] - spec_puorigen[:, tind]) #* 100 / (spectra_matrix[:, tind]/normalize)
            #plt.step(energy_data, cur_pcnt_diff, where='mid')
            #plt.xlabel('Energy [MeV]')
            #plt.ylabel('Absolute Difference')
            #plt.title(f'{t}s')
            #plt.savefig(current_path+f'{tind}-{name}-diff.png')
            #plt.close()

            # Do difference of average energy over time
            avg_pure = [spec_puorigen[i, tind] * e for i, e in enumerate(energy_data)]
            avg_iaea = [spec_oriaea[i, tind] * e for i, e in enumerate(energy_data)]
            #actual_cur_pcnt_diff = abs(spec_oriaea[:, tind] - spec_puorigen[:, tind]) / (0.5*spec_oriaea[:, tind] + 0.5*spec_puorigen[:, tind]) * 100
            diff_avgs.append(abs(sum(avg_pure) - sum(avg_iaea)))
            #l2_norms.append(np.linalg.norm(cur_pcnt_diff))
        plt.close()
        plt.plot(time_data, diff_avgs)
        plt.xlabel('Time [s]')
        plt.ylabel('Difference of Average Energies [MeV]')
        plt.savefig(current_path+f'l2normdiff.png')
        plt.close()
        print('generating movie')
        movie_start = time.time()
        misc_funcs.movie_gen(current_path, len(time_data))
        print(f'movie took {time.time() - movie_start} s')


        # 2D plot gen
        fig, ax = plt.subplots()
        z = abs(spec_oriaea - spec_puorigen)
        c = ax.pcolormesh(basex, basey, z, cmap='magma')
        cbar = fig.colorbar(c, ax=ax)
        cbar.set_label(f'Absoulte Difference in {name}')
        plt.xlabel('Time [s]')
        plt.ylabel('Energy [MeV]')
        plt.tight_layout()
        plt.savefig(current_path+f'{name}-diff.png')
        plt.close()
 
        fig, ax = plt.subplots()
        z = spec_puorigen
        c = ax.pcolormesh(basex, basey, z, cmap='magma')
        cbar = fig.colorbar(c, ax=ax)
        cbar.set_label(f'{name}')
        plt.xlabel('Time [s]')
        plt.ylabel('Energy [MeV]')
        plt.tight_layout()
        plt.savefig(current_path+f'{name}-Pure-ORIGEN.png')
        plt.close()

        fig, ax = plt.subplots()
        z = spec_oriaea
        c = ax.pcolormesh(basex, basey, z, cmap='magma')
        cbar = fig.colorbar(c, ax=ax)
        cbar.set_label(f'{name}')
        plt.xlabel('Time [s]')
        plt.ylabel('Energy [MeV]')
        plt.tight_layout()
        plt.savefig(current_path+f'{name}-ORIGEN-IAEA.png')
        plt.close()


    if spectra_expstrp_oria:
        # Perform the exponential stripping method
        #   to generate group spectral fits.
        #   Only available with isotope data.

        pathname = 'spectra-expstrp-ORIAEA/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        # Concentrations 
        ORIGEN_result = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        ORIGEN_dict = ORIGEN_result.ensdf_matcher(ensdf_dict,
                                                  0,
                                                  target)
        # Spectra
        valid_list = list()
        max_energy = 2
        ORIGEN_dict, valid_list = spectra_handler.SPECTRA(energy_mesh, times).update_dict(ORIGEN_dict, ensdf_handler)

        lst_sqrs = lls.LLS(times, None, fissions, efficiency, fit_groups)
        norm_six_spectra = lst_sqrs.spectra_abun_exp_strip(ORIGEN_dict, len(energy_mesh), ori_iaea_abuvec, valid_list)
        for group_ind, group_result in enumerate(norm_six_spectra):
            if spectra_normalized:
                normalize = np.sum(np.abs(group_result))
                name = 'Probability'
            else:
                normalize = 1
                name = 'Neutron Intensity'
            plt.step(energy_data, group_result/normalize, where='mid')
            plt.xlabel('Energy [MeV]')
            plt.ylabel(f'{name} / {energy_bin_kev} keV')
            plt.tight_layout()
            plt.savefig(current_path+f'group{group_ind+1}-spectra.png')
            plt.close()
            

    if spectra_lstsq_oriaea:
        # Use least squares fractional method to determine group spectra

        pathname = 'spectra-oriaea-lstsq-fit/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        lamvec = ori_iaea_lamvec
        abuvec = ori_iaea_abuvec

        ORIGEN_result = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        ORIGEN_dict = ORIGEN_result.ensdf_matcher(ensdf_dict,
                                                  0,
                                                  target)
        spectrum_dealer = spectra_handler.SPECTRA(energy_mesh, times)
        ORIGEN_dict, valid_list = spectrum_dealer.update_dict(ORIGEN_dict, ensdf_handler)
        group_spectra = spectrum_dealer.least_squares_spectrum(ori_iaea_lamvec, ORIGEN_dict, valid_list)
        for group_ind, group_result in enumerate(group_spectra):
            if spectra_normalized:
                normalize = np.sum(np.abs(group_result))
                name = 'Probability'
            else:
                normalize = 1
                name = 'Neutron Intensity'
            plt.step(energy_mesh, group_result/normalize, where='mid')
            plt.xlabel('Energy [MeV]')
            plt.ylabel(f'{name} / {energy_bin_kev} keV')
            plt.tight_layout()
            plt.savefig(current_path+f'group{group_ind+1}-spectra.png')
            plt.close()


        # Compare group spectra results with actual spectra
        for tind, t in enumerate(times):
            # Actual spectra
            plt.step(energy_mesh, spectral_matrix[:, tind], where='mid',
                    label='Data')
            plt.xlabel('Energy [MeV]')
            # Fit spectra
            fit_spectra = np.zeros(len(energy_mesh))
            for group_ind, group_result in enumerate(group_spectra):
                lam = lamvec[group_ind]
                abu = abuvec[group_ind]
                if irradiation == 'pulse':
                    a_val = lam * abu
                elif irradiation == 'infinite':
                    a_val = abu
                fit_spectra += a_val * np.exp(-lam * t) * group_result
            if spectra_normalized:
                name = 'Probability'
                normalize = np.sum(np.abs(fit_spectra))
            else:
                name = 'Counts'
                normalize = 1/(fissions * efficiency)
            plt.ylabel(f'{name} / {energy_bin_kev} keV')
            plt.step(energy_mesh, fit_spectra/normalize, where='mid', 
                    label=f'{fit_groups} Fit')
            #plt.title(f'{t}s')
            plt.legend()
            plt.savefig(current_path+f'{tind}.png')
            plt.close()
        misc_funcs.movie_gen(current_path, len(times))

    if alt_spc_lstsq_oriaea:
        # Use least squares to fit the spectral matrix data
        pathname = 'spectra-alt-oriaea-lstsq-fit/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        lamvec = ori_iaea_lamvec
        abuvec = ori_iaea_abuvec
        lamerr = ori_iaea_lamerr
        abuerr = ori_iaea_abuerr

        print('Performing least squares with count matrix')
        group_spectra, group_resid, group_errs = spectrum_dealer.pcnt_data_recon_lstsq(count_matrix, lamvec, lamerr, abuvec, abuerr)
        print(group_errs)
        for group_ind, group_result in enumerate(group_spectra):
            if spectra_normalized:
                normalize = np.sum(np.abs(group_result))
                name = 'Probability'
                sum_err = np.sqrt(np.sum(group_errs[group_ind]**2))
                use_err = np.sqrt( (group_errs[group_ind] / normalize)**2 + (-group_result/normalize**2 * sum_err)**2 )
            else:
                normalize = 1
                name = 'Neutron Intensity'
                use_err = group_errs[group_ind]
            plt.step(energy_mesh, group_result/normalize, where='mid')
            plt.fill_between(energy_mesh, group_result/normalize + use_err, group_result/normalize - use_err, alpha=alpha/2, step='mid')
            plt.xlabel('Energy [MeV]') 
            plt.ylabel(f'{name} / {energy_bin_kev} keV')
            plt.tight_layout()
            plt.savefig(current_path+f'group{group_ind+1}-spectra.png')
            plt.close()

        # Compare group spectra results with actual spectra
        n = 1
        cur_begin_time = time.time()
        for tind, t in enumerate(times):
            # Only print later time values
            if tind < len(times) - 10:
                continue
            if tind/len(times) >= 0.1 * n:
                n += 1
                net_time = time.time() - cur_begin_time
                full_complete_time = net_time / (0.1 * (n-1))
                print(f'    Progress: {round(tind/len(times) * 100)}% in {round(net_time, 0)}s')
                print(f'        Estimated completion in {round(full_complete_time - net_time, 0)}s')
            # Actual spectra
            plt.step(energy_mesh, spectral_matrix[:, tind], where='mid',
                    label='Data')
            plt.xlabel('Energy [MeV]')
            # Fit spectra
            fit_spectra = np.zeros(len(energy_mesh))
            fit_err = np.zeros(len(energy_mesh))
            for group_ind, group_result in enumerate(group_spectra):
                lam = lamvec[group_ind]
                abu = abuvec[group_ind]
                if irradiation == 'pulse':
                    a_val = lam * abu
                elif irradiation == 'infinite':
                    a_val = abu
                fit_spectra += a_val * np.exp(-lam * t) * group_result
                if irradiation == 'pulse':
                    fit_err += ( (abuerr[group_ind] * lam * np.exp(-lam * t) * group_result)**2 +
                        (lamerr[group_ind] * np.exp(-lam * t) * (abu * group_result - abu * group_result * t * lam))**2 + 
                        (group_errs[group_ind] * abu * lam * np.exp(-lam * t))**2 )**2
                elif irradiation == 'infinite':
                    fit_err += ( (abuerr[group_ind] * np.exp(-lam * t) * group_result)**2 +
                        (lamerr[group_ind] * np.exp(-lam * t) * (abu * group_result * t))**2 + 
                        (group_errs[group_ind] * abu * np.exp(-lam * t))**2 )**2
            use_err = np.zeros(len(energy_mesh))
            if spectra_normalized:
                name = 'Probability'
                normalize = np.sum(np.abs(fit_spectra))
                sum_err = np.sqrt(np.sum(fit_err**2))
                use_err = np.sqrt( (fit_err / normalize)**2 + (-fit_spectra/normalize**2 * sum_err)**2 )
            else:
                name = 'Counts'
                normalize = 1/(fissions * efficiency) 
                use_err = np.sqrt(fit_err)
            plt.ylabel(f'{name} / {energy_bin_kev} keV')
            plt.step(energy_mesh, fit_spectra/normalize, where='mid', 
                    label=f'{fit_groups} Fit') 
            plt.fill_between(energy_mesh, fit_spectra/normalize - use_err, fit_spectra/normalize + use_err, alpha=alpha, step='mid')
            #plt.title(f'{t}s')
            plt.legend()
            plt.savefig(current_path+f'{tind}.png')
            plt.close()
        misc_funcs.movie_gen(current_path, len(times))

    if spec_compare_oriaea:
        # Compare the group fits and plot against each other

        pathname = 'spectra-compare-oriaea-lstsq-fits/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        lamvec = ori_iaea_lamvec
        lamerr = ori_iaea_lamerr
        abuvec = ori_iaea_abuvec
        abuerr = ori_iaea_abuerr

        ORIGEN_result = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        ORIGEN_dict = ORIGEN_result.ensdf_matcher(ensdf_dict,
                                                  0,
                                                  target)
        spectrum_dealer = spectra_handler.SPECTRA(energy_mesh, times)

        ORIGEN_dict, valid_list = spectra_handler.SPECTRA(energy_mesh, times).update_dict(ORIGEN_dict, ensdf_handler)

        print('Performing least squares with spectral matrix')
        alt_group_spectra, group_resid, group_errs = spectrum_dealer.pcnt_data_recon_lstsq(count_matrix, lamvec, lamerr, abuvec, abuerr)
        q_group_spectra = spectrum_dealer.least_squares_spectrum(ori_iaea_lamvec, ORIGEN_dict, valid_list)
        for group_ind, group_result in enumerate(q_group_spectra):
            if spectra_normalized:
                normalize_q = np.sum(np.abs(q_group_spectra[group_ind]))
                normalize_alt = np.sum(np.abs(alt_group_spectra[group_ind]))
                name = 'Probability'
            else:
                normalize_q = 1
                normalize_alt = 1
                name = 'Neutron Intensity'
            plt.step(energy_mesh, q_group_spectra[group_ind]/normalize_q, where='mid',
                    label = 'Fractional Fit')
            plt.step(energy_mesh, alt_group_spectra[group_ind]/normalize_alt, where='mid',
                    label = 'Data Fit')
            plt.xlabel('Energy [MeV]') 
            plt.ylabel(f'{name} / {energy_bin_kev} keV')
            plt.legend()
            plt.tight_layout()
            plt.savefig(current_path+f'group{group_ind+1}-spectra.png')
            plt.close()

    if spectra_puori_fit:
        # Create spectral group fits for the pure ORIGEN data

        pathname = 'spectra-alt-puorigen-lstsq-fit/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)
        lamvec = pure_ori_lamvec
        lamerr = pure_ori_lamerr
        abuvec = pure_ori_abuvec
        abuerr = pure_ori_abuerr

        scale_builder = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)

        time_data, energy_data, spectra_matrix, bin_data = scale_builder.origen_spectra_parser()
        print('Performing least squares with spectra matrix')
        spectrum_dealer = spectra_handler.SPECTRA(energy_data, times)
        group_spectra, group_resid, group_err = spectrum_dealer.pcnt_data_recon_lstsq(count_matrix, lamvec, lamerr, abuvec, abuerr)
        if spectra_normalized:
            for tind in range(len(time_data)):
                spectra_matrix[:, tind] = spectra_matrix[:, tind] / np.sum(spectra_matrix[:, tind])
        for group_ind, group_result in enumerate(group_spectra):
            if spectra_normalized:
                normalize = np.sum(np.abs(group_result))
                name = 'Probability'
            else:
                normalize = 1
                name = 'Neutron Intensity'
            plt.step(energy_data, group_result/normalize, where='mid')
            plt.xlabel('Energy [MeV]') 
            plt.ylabel(f'{name} / {energy_bin_kev} keV')
            plt.tight_layout()
            plt.savefig(current_path+f'group{group_ind+1}-spectra.png')
            plt.close()

        # Compare group spectra results with actual spectra
        n = 1
        cur_begin_time = time.time()
        print('2518 and 2504 removed plotting')
#        for tind, t in enumerate(time_data):
#            if tind/len(time_data) >= 0.1 * n:
#                n += 1
#                net_time = time.time() - cur_begin_time
#                full_complete_time = net_time / (0.1 * (n-1))
#                print(f'    Progress: {round(tind/len(time_data) * 100)}% in {round(net_time, 0)}s')
#                print(f'        Estimated completion in {round(full_complete_time - net_time, 0)}s')
#            # Actual spectra
#            plt.step(energy_data, spectra_matrix[:, tind], where='mid',
#                    label='Data')
#            plt.xlabel('Energy [MeV]')
#            # Fit spectra
#            fit_spectra = np.zeros(len(energy_data))
#            for group_ind, group_result in enumerate(group_spectra):
#                lam = lamvec[group_ind]
#                abu = abuvec[group_ind]
#                fit_spectra += lam * abu * np.exp(-lam * t) * group_result
#            if spectra_normalized:
#                name = 'Probability'
#                normalize = np.sum(np.abs(fit_spectra))
#            else:
#                name = 'Counts'
#                normalize = 1/(fissions * efficiency)
#            plt.ylabel(f'{name} / {energy_bin_kev} keV')
#            plt.step(energy_data, fit_spectra/normalize, where='mid', 
#                    label=f'{fit_groups} Fit')
#            plt.title(f'{t}s')
#            plt.legend()
#            plt.savefig(current_path+f'{tind}.png')
#            plt.close()
#        misc_funcs.movie_gen(current_path, len(time_data))

    if keepin_data_origen:
        # Generate plot comparing Keepin data and ORIGEN data
        pathname = 'request-713/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        name = '6keepin235fast'
        keepin_response = keepin_handler.KEEPIN(name)
        plt.plot(keepin_response.true_data_time, keepin_response.true_data_resp,
                 label='Approx Keepin Data', linestyle='', marker='.')
        
        scale_builder = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        t_ori, pure_delnu_ori = scale_builder.origen_delnu_parser(target)
        n_per_f = misc_funcs.delnu_per_fiss(t_ori, pure_delnu_ori, fissions, efficiency)
        plt.plot(t_ori, pure_delnu_ori, label='ORIGEN Output')
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'keepin-origen-data.png')
        plt.close()

#    if movie_gen:
#        # Create a gif
#        pathname = 'movie/'
#        current_path = imdir+pathname
#        misc_funcs.dir_handle(current_path)
#
#        filenames = list()
#        for i in range(903):
#            name = f'{i}.png'
#            filenames.append(movie_target_path + name)
#
#        import imageio
#        images = []
#        for filename in filenames:
#            images.append(imageio.v2.imread(filename))
#        imageio.mimsave(f'{current_path}movie.gif', images) 
    
    if sanity_check:
        # Check that Keepin group fit aligns with Keepin data
        pathname = 'Keepin-validation/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        name = '6keepin235fast'
        keepin_response = keepin_handler.KEEPIN(name)
        plt.plot(keepin_response.true_data_time, keepin_response.true_data_resp,
                 label='Approx Keepin Data', linestyle='', marker='.')

        group_times_use = np.arange(0, 330, 1)
        group_delnu, group_errs = prefab_group_delnu(group_times_use,
                       keepin_lamvec,
                       keepin_lamerr,
                       keepin_abuvec,
                       keepin_abuerr,
                       fissions,
                       efficiency)
        misc_funcs.delnu_per_fiss(group_times_use, group_delnu, fissions, efficiency)
        plt.plot(group_times_use, group_delnu, label='Keepin Group Fit')
        print(f'Efficiency needed @330s: {keepin_response.true_data_resp[-1] / group_delnu[-1]}')
        print(f'Efficiency needed @0s: {keepin_response.true_data_resp[0] / group_delnu[0]}')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Count Rate [#/s]')
        plt.yscale('log')
        plt.tight_layout()
        plt.savefig(current_path+f'count-compare.png')
        plt.close()

    if keepin_base_compare:
        # Keepin normalized group fits comparison with uncertainties
        pathname = 'request-713/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        keepin_delnu, keepin_errs = prefab_group_delnu(times,
                       keepin_lamvec,
                       keepin_lamerr,
                       keepin_abuvec,
                       keepin_abuerr,
                       fissions,
                       efficiency)

        keepin_delnu = np.asarray(keepin_delnu)
        keepin_errs = np.asarray(keepin_errs)
        brady_england_delnu, be_errs = prefab_group_delnu(times,
                       be_lamvec,
                       be_lamerr,
                       be_abuvec,
                       be_abuerr,
                       fissions,
                       efficiency)
        brady_england_delnu = np.asarray(brady_england_delnu)
        be_errs = np.asarray(be_errs)

        ori_group_delnu, ori_group_errs = prefab_group_delnu(times,
                       pure_ori_lamvec,
                       pure_ori_lamerr,
                       pure_ori_abuvec,
                       pure_ori_abuerr,
                       fissions,
                       efficiency)

        ori_group_delnu = np.asarray(ori_group_delnu)
        ori_group_errs  = np.asarray(ori_group_errs)


        iaea_group_delnu, iaea_group_errs = prefab_group_delnu(times,
                       ori_iaea_lamvec,
                       ori_iaea_lamerr,
                       ori_iaea_abuvec,
                       ori_iaea_abuerr,
                       fissions,
                       efficiency)

        iaea_group_delnu = np.asarray(iaea_group_delnu)
        iaea_group_errs  = np.asarray(iaea_group_errs)


        #scale_group_delnu, _ = prefab_group_delnu(times,
        #               scale_lamvec,
        #               pure_ori_lamerr,
        #               scale_abuvec,
        #               pure_ori_abuerr,
        #               fissions,
        #               efficiency)


        #plt.plot(times, scale_group_delnu / keepin_delnu, label='SCALE Fit')
        plt.plot(times, keepin_delnu / keepin_delnu, label='Keepin [6]')
        plt.plot(times, brady_england_delnu / keepin_delnu, label='Brady/England [11]')
        plt.plot(times, ori_group_delnu / keepin_delnu, label='Pure ORIGEN')
        plt.plot(times, iaea_group_delnu / keepin_delnu, label='IAEA-ORIGEN')
        plt.xlabel('Time [s]')
        plt.ylabel('Keepin Normalized Counts')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'keepin-normalized-counts-clean.png')
        plt.close()

        
        #keepin_err = np.zeros(len(keepin_delnu))
        keepin_err = np.sqrt( (keepin_delnu/keepin_delnu**2 * keepin_errs)**2 + (keepin_errs/keepin_delnu)**2 )
        ori_err = np.sqrt( (ori_group_delnu/keepin_delnu**2 * keepin_errs)**2 + (ori_group_errs/keepin_delnu)**2 )
        bradeng_err = np.sqrt( (brady_england_delnu/keepin_delnu**2 * keepin_errs)**2 + (be_errs/keepin_delnu)**2 )

        misc_funcs.multplt(times, keepin_delnu / keepin_delnu, label='Keepin [6]', errors=keepin_err, alpha=alpha, cnst_line=True)
        misc_funcs.multplt(times, brady_england_delnu / keepin_delnu, label='Brady/England [11]', errors=bradeng_err, alpha=alpha, cnst_line=True)
        misc_funcs.multplt(times, ori_group_delnu / keepin_delnu, label='Pure ORIGEN', errors=ori_err, alpha=alpha, cnst_line=True)
        misc_funcs.multplt(times, iaea_group_delnu / keepin_delnu, label='IAEA-ORIGEN', errors=ori_err, alpha=alpha, cnst_line=True)
        plt.xscale('log')

        #misc_funcs.multplt(times, scale_group_delnu / keepin_delnu, label='SCALE Fit', alpha=alpha)

        plt.xlabel('Time [s]')
        plt.ylabel('Keepin Normalized Counts')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'keepin-normalized-counts.png')
        plt.close()

    if integral_keep_ori_be:
        # Integral as a function of time for the Keepin, ORIGEN, and BradyEngland group fits
        pathname = 'request-713/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        sf = 3 # Start From

        #name = '6keepin235fast'
        #keepin_response = keepin_handler.KEEPIN(name)
        #keepin_delnu, keepin_errs = keepin_response.simulate_instant(times, fissions, efficiency)
        keepin_delnu, keepin_errs = prefab_group_delnu(times,
                       keepin_lamvec,
                       keepin_lamerr,
                       keepin_abuvec,
                       keepin_abuerr,
                       fissions,
                       efficiency)

        keepin_delnu = np.asarray(keepin_delnu)
        keepin_errs = np.asarray(keepin_errs)
        brady_england_delnu, be_errs = prefab_group_delnu(times,
                       be_lamvec,
                       be_lamerr,
                       be_abuvec,
                       be_abuerr,
                       fissions,
                       efficiency)
        #name = '6brengland235fast'
        #brady_england_response = keepin_handler.KEEPIN(name)
        #brady_england_delnu, be_errs = brady_england_response.simulate_instant(times, fissions, efficiency)

        brady_england_delnu = np.asarray(brady_england_delnu)
        be_errs = np.asarray(be_errs)



        ori_group_delnu, ori_group_errs = prefab_group_delnu(times,
                       pure_ori_lamvec,
                       pure_ori_lamerr,
                       pure_ori_abuvec,
                       pure_ori_abuerr,
                       fissions,
                       efficiency)

        ori_group_delnu = np.asarray(ori_group_delnu)
        ori_group_errs  = np.asarray(ori_group_errs)

        keep_int = np.zeros(len(times))
        ori_int = np.zeros(len(times))
        be_int = np.zeros(len(times))
        keep_int_err = np.zeros(len(times))
        ori_int_err = np.zeros(len(times))
        be_int_err = np.zeros(len(times))
        dt = np.mean(np.diff(times))
        for tind, t in enumerate(times):
            # Generate integral for each group at each time
            keep_int[tind] = np.trapz(keepin_delnu[:tind], x=times[:tind])
            ori_int[tind] = np.trapz(ori_group_delnu[:tind], x=times[:tind])
            be_int[tind] = np.trapz(brady_england_delnu[:tind], x=times[:tind])
            # Uncertainties
            if tind > 2:
                keep_add = np.sqrt(3) * keepin_errs[1:tind-1]
                ori_add = np.sqrt(3) * ori_group_errs[1:tind-1]
                be_add = np.sqrt(3) * be_errs[1:tind-1]
            else:
                keep_add = 0
                ori_add = 0
                be_add = 0
            keep_int_err[tind] = dt/2 * np.sqrt( np.sum(keepin_errs[:tind]**2) + np.sum(keep_add**2))
            ori_int_err[tind] = dt/2 * np.sqrt( np.sum(ori_group_errs[:tind]**2) + np.sum(ori_add**2))
            be_int_err[tind] = dt/2 * np.sqrt( np.sum(be_errs[:tind]**2) + np.sum(be_add**2))


        misc_funcs.multplt(times[sf:], keep_int[sf:], label='Keepin Integral', errors=keep_int_err[sf:])
        misc_funcs.multplt(times[sf:], ori_int[sf:], label='ORIGEN Integral', errors=ori_int_err[sf:])
        misc_funcs.multplt(times[sf:], be_int[sf:], label='Brady/England Integral', errors=be_int_err[sf:])
        plt.xlabel('Time [s]')
        plt.ylabel('Net Counts')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'Integral-fits.png')
        plt.close()


        misc_funcs.multplt(times[sf:], keep_int[sf:] / fissions / efficiency, label='Keepin Integral', errors=keep_int_err[sf:]/fissions/efficiency)
        misc_funcs.multplt(times[sf:], ori_int[sf:] / fissions / efficiency, label='ORIGEN Integral', errors=ori_int_err[sf:]/fissions/efficiency)
        misc_funcs.multplt(times[sf:], be_int[sf:] / fissions / efficiency, label='Brady/England Integral', errors=be_int_err[sf:]/fissions/efficiency)
        plt.xlabel('Time [s]')
        plt.ylabel('Yield [delnu per fiss]')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'Yield-fits.png')
        plt.close()

        keepin_err = np.zeros(len(keepin_delnu))
        ori_err = np.sqrt( (ori_int/keep_int**2 * keep_int_err)**2 + (ori_int_err/keep_int)**2 )
        bradeng_err = np.sqrt( (be_int/keep_int**2 * keep_int_err)**2 + (be_int_err/keep_int)**2 )

        misc_funcs.multplt(times[sf:], keep_int[sf:] / keep_int[sf:], label='Keepin Integral', errors=keepin_err[sf:])
        misc_funcs.multplt(times[sf:], ori_int[sf:] / keep_int[sf:], label='ORIGEN Integral', errors=ori_err[sf:])
        misc_funcs.multplt(times[sf:], be_int[sf:] / keep_int[sf:], label='Brady/England Integral', errors=bradeng_err[sf:])
        plt.xlabel('Time [s]')
        plt.ylabel('Keepin Normalized Counts')
        plt.legend()
        plt.tight_layout()
        plt.savefig(current_path+f'Integral-fits-normalized.png')
        plt.close()


    if display_endf_spectra:
        # Generate the group spectra from ENDF data
        pathname = 'spectra-endf-display/'
        current_path = imdir+pathname
        misc_funcs.dir_handle(current_path)

        endf_spectra = spectra_handler.SPECTRA(energy_mesh, times).read_endf_excel(endf_spectra_filename, endf_spectra_sheetname)
        for gind, group in enumerate(endf_spectra):
            plt.plot(energy_mesh, endf_spectra[gind, :])
            plt.xlabel('Energy [MeV]')
            if spectra_normalized:
                normalize = np.sum(endf_spectra[gind, :])
                name = 'Probability'
            else:
                normalize = 1
                # The ENDF data is already normalized to probability
                name = 'Probability'
            plt.ylabel(f'{name} / {energy_bin_kev} keV')
            plt.tight_layout()
            plt.savefig(current_path+f'group-{gind+1}.png')
            plt.close()

    if yield_contributions:
        # Display the percentage contribution of yield from each isotope
        #   as well as which isotope have spectral data.
        
        
        ORIGEN_result = scale_handler.SCALE(ORIGEN_out,
                                            fissions,
                                            efficiency,
                                            normalize_value,
                                            volume,
                                            mass_normalize)
        ensdf_gen = ensdf_handler.ENSDF(ensdf_fname, ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        ORIGEN_dict = ORIGEN_result.ensdf_matcher(ensdf_dict,
                                                  0,
                                                  target)
        # Spectra
        ORIGEN_dict, valid_list = spectra_handler.SPECTRA(energy_mesh, times).update_dict(ORIGEN_dict, ensdf_handler)
        yield_list = list()
        iso_list = list()
        spectra_bool_list = list()
        # Sort by yield
        for isotope in ORIGEN_dict:
            Pn = ORIGEN_dict[isotope]['emission'][0]
            atom0 = ORIGEN_dict[isotope]['conc'][0][0]
            iso_yield = Pn * atom0 / fissions
            yield_list.append(iso_yield)
            iso_list.append(isotope)
            if 'spectrum_values' in ORIGEN_dict[isotope]:
                spectra_bool_list.append(True)
            else:
                spectra_bool_list.append(False)
        # Print
        net_yield = np.sum(yield_list)
        yield_list, iso_list, spectra_bool_list = zip(*sorted(zip(yield_list, iso_list, spectra_bool_list)))
        yield_list = list(yield_list)
        iso_list = list(iso_list)
        spectra_bool_list = list(spectra_bool_list)
        net_present = 0
        for index, yld in enumerate(yield_list):
            yld_pcnt = round(yld / net_yield * 100, 3)
            print(f'{iso_list[index]} : {spectra_bool_list[index]} : {yld_pcnt}%')
            if spectra_bool_list[index]:
                net_present += yld_pcnt
        print(f'Total spectra accounted for: {net_present}%')
        print(f'Predicted total yield: {net_yield}')

        print(f'{len(ORIGEN_dict) / len(ensdf_dict) * 100}% of isotopes in IAEA database in ORIGEN')

        





    # Done
    end = time.time()
    print(f'Took {round(end - begin, 0)}s')
