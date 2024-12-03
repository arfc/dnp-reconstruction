import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid
import misc_funcs
import settings


class KEEPIN:
    """
    Handles 1957 Keepin data usage
    """
    def __init__(self,
                 name=None):
        """
        Initialize

        Parameters
        ----------
        name : str
            Name of group fit

        Returns
        -------
        None
        """
        self.name = name
        self.true_data_time = [0.1,   0.3,  1.2, 1.7, 21.5, 190,  330]
        self.true_data_resp = [2.6E6, 2E6, 1E6, 8E5, 4E4, 4E2, 5E1]
        return

    def data_store(self):
        """
        Stores the Keepin data. Could be reorganized to be a separate file.

        Parameters
        ----------
        None

        Returns
        -------
        group_data : dict
            key : string
                Name of group
            value : list
                Half life, error, rel_abundance, err, yield, err
        total_yield : list
            yield, err
        """
        group_data = dict()
        if self.name == '6keepin235fast':
            group_data['g1'] = [54.51, 0.94, 0.038, 0.003, 0.00063, 0.00005]
            group_data['g2'] = [21.84, 0.54, 0.213, 0.005, 0.00351, 0.00011]
            group_data['g3'] = [6.00, 0.17,  0.188, 0.016, 0.00310, 0.00028]
            group_data['g4'] = [2.23, 0.06, 0.407, 0.007, 0.00672, 0.00023]
            group_data['g5'] = [0.496, 0.029, 0.128, 0.008, 0.00211, 0.00015]
            group_data['g6'] = [0.179, 0.017, 0.026, 0.003, 0.00043, 0.00005]
            self.groups = 6
            total_yield = [0.0165, 0.0005]
            pm = u'\u00b1'
            print(f'{self.name} n/f: {total_yield[0]} {pm} {total_yield[1]}\n')
        elif self.name == '6brengland235fast':
            group_data['g1'] = [52.116, 0, 0.0350, 0, 0.000721, 1.442E-6]
            group_data['g2'] = [21.197, 0, 0.1807, 0, 0.00372242, 7.4E-6]
            group_data['g3'] = [5.7380, 0, 0.1725, 0, 0.0035535, 7.1E-6]
            group_data['g4'] = [2.2891, 0, 0.3868, 0, 0.00796808, 1.6E-5]
            group_data['g5'] = [0.8159, 0, 0.1586, 0, 0.00326716, 6.5E-6]
            group_data['g6'] = [0.2430, 0, 0.0664, 0, 0.00136784, 2.7E-6]
            self.groups = 6
            total_yield = [0.0206, 0.002]
            print(f'{self.name} n/f: {total_yield[0]}\n')
        elif self.name == 'test_1':
            group_data['g1'] = [10, 0, 1, 0, 0.1, 0]
            self.groups = 1
            total_yield = [0.1, 0]
            print(f'{self.name} n/f: {total_yield[0]}\n')
        elif self.name == 'test_2':
            group_data['g1'] = [10, 0, 1/3, 0, 0.3, 0]
            group_data['g2'] = [1, 0, 2/3, 0, 0.6, 0]
            self.groups = 2
            total_yield = [0.9, 0]
            print(f'{self.name} n/f: {total_yield[0]}\n')
        elif self.name == 'test_3':
            group_data['g1'] = [10, 0, 1/4, 0, 10/4, 0]
            group_data['g2'] = [1, 0, 1/4, 0, 10/4, 0]
            group_data['g3'] = [0.1, 0, 1/2, 0, 10/2, 0]
            self.groups = 3
            total_yield = [10, 0]
            print(f'{self.name} n/f: {total_yield[0]}\n')
        elif self.name == 'test_6':
            group_data['g1'] = [54.51, 0.94, 0.038, 0.003, 0.063, 0.005]
            group_data['g2'] = [21.84, 0.54, 0.213, 0.005, 0.351, 0.011]
            group_data['g3'] = [6.00, 0.17,  0.188, 0.016, 0.310, 0.028]
            group_data['g4'] = [2.23, 0.06, 0.407, 0.007, 0.672, 0.023]
            group_data['g5'] = [0.496, 0.029, 0.128, 0.008, 0.211, 0.015]
            group_data['g6'] = [0.179, 0.017, 0.026, 0.003, 0.043, 0.005]
            self.groups = 6
            total_yield = [0.0165, 0.0005]
            print(f'{self.name} n/f: {total_yield[0]}\n')
        else:
            print(f'{self.name} not available')
            raise Exception('Unavailable name')

        return group_data, total_yield
        
    def simulate_instant(self, times, fissions, efficiency):
        """
        Simulate irradiation followed by measurement based
            on group values. Creates the net delayed neutron output over time.

        Parameters
        ----------
        times : list
            List of times to generate data points
        fissions : float
            Number of fissions in the sample

        Returns
        -------
        delnu : list
            Delayed neutrons at given times
        cur_err : list
            Error of the delayed neutron count at given times
        
        """
        gdata, tot_yield = self.data_store()
        delnu = list()
        groups = list()
        errs = list()
        for each in range(self.groups):
            groups.append('g' + str(each + 1))
        lead_term = fissions * efficiency
        for t in times:
            detect = 0
            cur_err = 0
            for g in groups:
                a = gdata[g][4]
                lam = np.log(2)/gdata[g][0]
                err_a = gdata[g][5]
                err_lam = np.log(2) / gdata[g][0]**2 * gdata[g][1]
                if settings.irradiation == 'pulse':
                    term_a = lam * np.exp(-lam * t)
                    term_lam = a * (1 - lam*t) * np.exp(-lam * t)
                    detect += (a * lam * np.exp(-lam * t))
                elif settings.irradiation == 'infinite':
                    term_a = np.exp(-lam * t)
                    term_lam = (a * t * np.exp(-lam * t))
                    detect == (a * np.exp(-lam * t))
                cur_err += (term_a * err_a)**2 + (term_lam * err_lam)**2
            detect = detect * lead_term
            cur_err = np.sqrt(cur_err) * lead_term
            delnu.append(detect)
            errs.append(cur_err)
        n_per_f = misc_funcs.delnu_per_fiss(times,
                                    delnu,
                                    fissions,
                                    efficiency)
        if type(delnu) == type(list):
            delnu = np.array(delnu)
        if type(errs) == type(list):
            delnu = np.array(errs)
        return delnu, errs

    def simulate_lin_solve(self, times, soln_vec, fissions, efficiency,
                           a_errs=None, lam_errs=None):
        """
        Simulate instantaneous irradiation using the groups generated by the
            linear least squares method

        Parameters
        ----------
        times : list
            List of times to generate data points
        soln_vec : numpy array
            12x1 a_i*lambda_i, lambda_i for 6 groups
            
        Returns
        -------
        delnu : list
            Delayed neutrons at given times
        errors : list
            Errors at given times
        """
        delnu = list()
        groups = list()
        errors = list()
        normalize = 0
        err_solve = False
        if type(a_errs) != type(None) and \
           type(lam_errs) != type(None):
            err_solve = True
        for index, g in enumerate(soln_vec):
            if index % 2 == 0:
                normalize += g
                groups.append(index)
        for t in times:
            detect = 0
            err = 0
            for index, g in enumerate(groups):
                lami = soln_vec[g+1]
                ai = soln_vec[g] / soln_vec[g+1]
                #detect += (soln_vec[g] * np.exp(-soln_vec[g+1] * t))
                if settings.irradiation == 'pulse':
                    a_val = lami * ai
                elif settings.irradiation == 'infinite':
                    a_val = ai
                detect += (a_val * np.exp(-lami * t))
                if err_solve:
                    ai_err = a_errs[index]
                    lami_err = lam_errs[index]
                    if settings.irradiation == 'pulse':
                        err += ((lami * np.exp(-lami * t) * ai_err)**2 +
                            (ai*(1-lami*t)*np.exp(-lami*t)*lami_err)**2)
                    elif settings.irradiation == 'infinite':
                        err += ((np.exp(-lami * t) * ai_err)**2 +
                            (ai*t*np.exp(-lami*t)*lami_err)**2)
            detect = fissions * efficiency * float(detect)
            delnu.append(detect)
            if err_solve:
                err = float(np.sqrt(err)) * fissions * efficiency
                errors.append(err)
            else:
                errors.append(0)
        return delnu, errors

def debug_run(deb_group):
    data_name = 'test_' + str(deb_group)
    keepin_response = KEEPIN(data_name)
    keepin_group_data, keepin_net_data = keepin_response.data_store()
    keepin_delnu = keepin_response.simulate_instant(times, fissions, efficiency)
    plt.plot(times, keepin_delnu, label=f'{deb_group} keepin')
    int_keepin_cnt = cumulative_trapezoid(keepin_delnu, x=times)
    tot_keepin_cnt = int_keepin_cnt[-1] - int_keepin_cnt[0]
    print(f'Max Keepin counts: {max(keepin_delnu)}')
    print(f'Total Kepein Counts: {tot_keepin_cnt}')
    # test simulate_lin_solve
    soln_vec = list()
    for g in keepin_group_data.keys():
        lami = np.log(2)/keepin_group_data[g][0]
        ai = keepin_group_data[g][4]
        soln_vec.append(lami * ai)
        soln_vec.append(lami)
        print(f'Lambda: {lami}')
        print(f'ai: {ai}')
    lls_delnu = keepin_response.simulate_lin_solve(times, soln_vec, fissions, efficiency)
    print(f'Soln: {soln_vec}')
    #print('Hardcoding group 2 test')
    #plt.plot(times,
    #         0.3*np.log(2)/10*np.exp(-np.log(2)/10 * times) + 0.6*np.log(2)*np.exp(-np.log(2) * times),
    #         label='Hard coded')
    plt.plot(times, lls_delnu, label=f'{deb_group} LLS')
    plt.yscale('log')
    plt.ylabel('Delayed Neutron Count Rate [#/s]')
    plt.xlabel('Time [s]')
    plt.legend()
    plt.show()
    return
        


if __name__ == '__main__':
    import ensdf_handler
    import misc_funcs
    dt = 0.1
    end_time = 330
    fissions = 1E16
    efficiency = 1.650637878787879e-07
    times = np.arange(0, end_time+dt, dt)
    default = True
    debug = False
    debug_group = 6

    if debug:
         debug_run(debug_group)

    if default:
        name = '6keepin235fast'
        keepin_response = KEEPIN(name)
        keepin_group_data, keepin_net_data = keepin_response.data_store()
        keepin_delnu, keepin_errs = keepin_response.simulate_instant(times, fissions, efficiency)
        ensdf_keepin_sim = ensdf_handler.ENSDF('./ensdf_data/eval_net.xlsx',
                                                 'Sheet1')
        #ensdf_normalized_eff = 0.006272693743682443 * efficiency #keepin_max / ensdf_max
        ensdf_keepin_delnu, ensdf_keepin_data = ensdf_keepin_sim.simulate_keepin_group_abun(keepin_group_data,
                                                                         times,
                                                                                            fissions,
                                                                                            efficiency)
        plt.plot(times, keepin_delnu, label='keepin')
        int_keepin_cnt = cumulative_trapezoid(keepin_delnu, x=times)
        tot_keepin_cnt = int_keepin_cnt[-1] - int_keepin_cnt[0]
        true_keepin_cnt = cumulative_trapezoid(keepin_response.true_data_resp,
                                   x=keepin_response.true_data_time)
        true_tot_keepin_cnt = true_keepin_cnt[-1] - true_keepin_cnt[0]
        print(f'Keepin n/f: {misc_funcs.delnu_per_fiss(times, keepin_delnu, fissions, efficiency)}')
        print(f'True Total Keepin Counts: {true_tot_keepin_cnt}')
        print(f'Calculated eff: {true_tot_keepin_cnt / (fissions * keepin_net_data[0])}')
        print(f'Max Keepin counts: {max(keepin_delnu)}')
        print(f'Total Kepein Counts: {tot_keepin_cnt}')
        predict_keepin_cnt = fissions * efficiency * keepin_net_data[0]
        print(f'Predicted Total Keepin Counts: {predict_keepin_cnt}')
        plt.plot(times, ensdf_keepin_delnu, label='ensdf-keepin')
        int_ensdf_cnt = cumulative_trapezoid(ensdf_keepin_delnu, x=times)
        tot_ensdf_cnt = int_ensdf_cnt[-1] - int_ensdf_cnt[0]
        print(f'Max ensdf counts: {max(ensdf_keepin_delnu)}')
        print(f'Total ensdf Counts: {tot_ensdf_cnt}')
        plt.plot(keepin_response.true_data_time, keepin_response.true_data_resp,
                label='True', linestyle='', marker='.')
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate [#/s]')
        plt.xlabel('Time [s]')
        plt.legend()
        plt.show()
