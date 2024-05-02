import numpy as np
import keepin_handler as keepin
import matplotlib.pyplot as plt
import scipy
import itertools
import math
import time
import settings
import scipy.optimize as scp
import statistics

class LLS:
    """
    Linear least squares implementation
    """

    def __init__(self, times, counts, fissions, efficiency, ngroups):
        """
        Initialize

        Parameters
        ----------
        times : list
            List of time values
        counts : list
            List of counts at each time

        Returns
        -------
        None
        """
        self.times = times
        self.counts = counts
        self.fissions = fissions
        self.efficiency = efficiency
        self.numgroups = ngroups
        return

    def spectra_abun_exp_strip(self, isotope_dict, num_energy_bins,
            group_yields, valid_isos):
        """
        Applies the exponential stripping method to generated 
            group spectra.
            Starting from the longest lived group, the isotopes
            are iterated through until the current sum of isotopic
            yields are equal or greater than the current group yield.
            If greater, than a fraction of the isotope is applied to that group
            and the next group.
            The next group is then iterated through, continuing from
            the current isotope.

        Parameters
        ----------
        isotope_dict : dict
            keys : str
                Isotope names
            values : dict
                keys : str
                    Categories of data for extraction (should include spectra)
                values : list
                    Associated values and uncertainties
        num_energy_bins : int
            Number of energy bins with datapoints
        group_yields : numpy array
            1D array of the yield for each group
        valid_isos : list
            List of strings containing names of isotopes with spectral data

        Returns
        -------
        norm_group_spectra : 2D numpy array
            Normalized matrix where rows are the groups and columns
                are the energy indeces. Normalized such that each group provides
                the probability of an emitted neutron having a given energy.
        """
        print(f'Total group yield: {np.sum(group_yields)}')
        norm_six_spectra = np.zeros((self.numgroups, num_energy_bins))
        running_isotope_yield = 0
        for cur_group in range(self.numgroups):
            print(f'Desired group {cur_group+1} yield: {group_yields[cur_group]}')
            yield_sum = 0
            while yield_sum < group_yields[cur_group]:
                # Determine next longest lived isotope
                longest_life = 0
                longest_lived = ''
                for isotope in valid_isos:
                    if isotope_dict[isotope]['halflife'][0] > longest_life:
                        longest_life = isotope_dict[isotope]['halflife'][0]
                        longest_lived = isotope
                if longest_life == 0:
                    print('Isotope yield sums did not allow for complete convergence')
                    running_isotope_yield += yield_sum
                    print(f'Total isotope yield: {running_isotope_yield}')
                    print(f'Total group yield needed: {np.sum(group_yields)}')
                    print(f'% Completeness: {running_isotope_yield / np.sum(group_yields) * 100}%')
                    print('Insufficient spectral data to account for total yield')
                    raise Exception
                # Scale its contribution according to yield after normalizing
                temp_spectra = isotope_dict[longest_lived]['spectrum_values'] 
                norm_spectra = temp_spectra / np.sum(temp_spectra)
                iso_yield = isotope_dict[longest_lived]['emission'][0] * isotope_dict[longest_lived]['conc'][0][0] / self.fissions
                yield_sum += iso_yield
                # Fractional contribution
                fraction = 1
                if yield_sum > group_yields[cur_group]:
                    yield_sum -= iso_yield
                    fraction = (group_yields[cur_group] - yield_sum) / iso_yield
                    norm_six_spectra[cur_group+1, :] += (1-fraction) * iso_yield / group_yields[cur_group] * norm_spectra
                    yield_sum += fraction * iso_yield
                norm_six_spectra[cur_group, :] += fraction * iso_yield / group_yields[cur_group] * norm_spectra
                isotope_dict.pop(longest_lived)
                valid_isos.pop(np.where(np.array(valid_isos) == longest_lived)[0][0])
                print(f'Current iso {longest_lived} yield: {fraction*iso_yield}')
                
            running_isotope_yield += yield_sum
            print(f'Group {cur_group+1} Yield: {yield_sum}')
            if fraction != 1:
                print(f'Contributed iso {longest_lived} yield: {(1-fraction)*iso_yield}')
        return norm_six_spectra


    def MC_abun_err(self,
                 lamvec,
                 lamerr,
                 cntvec,
                 cnterr=None,
                 fiserr=6.464567E13,
                    tot_iters=1000):
        """
        Uses a Monte Carlo approach to generate the abundance uncertainties

        Parameters
        ----------
        lamvec : array
            Decay constants for each group
        lamerr : array
            Associated errors for each decay constant
        cntvec : array
            Counts at each time step
        cnterr : array
            Uncertainty in counts at each time step
        fiserr : float
            Uncertainty in total fissions in sample
        tot_iters : int
            Total number of iterations to perform before stopping

        Returns
        -------
        abundance_errors : list
            Uncertainty in group yields
        """
        totyield = 0
        totyield_err = 0
        # Update b vector of errors
        stored_abundances = dict()
        naming_convention = 'a'
        abundance_errors = list()
        if type(cnterr) == type(None):
            cnterr = np.zeros(len(cntvec))
        berr = np.zeros((1, len(cnterr)))
        for errind, err in enumerate(cnterr):
            nd_err = (err * 1 / (self.efficiency * self.fissions))**2
            fs_err = (fiserr * (cntvec[errind] / (self.efficiency * self.fissions**2)))**2
            berr[0][errind] = np.sqrt(nd_err + fs_err)
        print(f'Performing abundance error search {tot_iters} times')
        n = 1
        lami_begin_time = time.time()
        for cur_iter in range(tot_iters):
            if cur_iter == 0:
                # Base solution
                print(f'Base Decays: {lamvec}')
                print(f'Errors: {lamerr}')
                A = np.zeros((len(self.times), self.numgroups))
                for index, t in enumerate(self.times):
                    for coldex, lam in enumerate(lamvec):
                        A[index, coldex] = lam * np.exp(-lam * t)
                b = np.array(cntvec) / (self.efficiency * self.fissions)
                #x, res = scp.nnls(A, b, maxiter=None)
                x, res, rank, s = np.linalg.lstsq(A, b.T, rcond=None)
                for ai_ind, ai in enumerate(x):
                    name = naming_convention + str(ai_ind + 1)
                    stored_abundances[name] = list()
                    stored_abundances[name].append(ai)

            if cur_iter/tot_iters >= 0.1 * n:
                n += 1
                net_time = time.time() - lami_begin_time
                full_complete_time = net_time / (0.1 * (n-1))
                print(f'    Progress: {round(cur_iter/tot_iters * 100)}% in {round(net_time, 0)}s')
                print(f'        Estimated completion in {round(full_complete_time - net_time, 0)}s')
            # Generate random values (each lam and all counts)
            num_rand_vals = self.numgroups + 1
            rand_vals = np.random.rand(num_rand_vals)
            
            # Assign decay constants and counts based on uncertainties and
            #   random values
            use_lam = np.zeros((1, self.numgroups))
            for each in range(self.numgroups):
                use_lam[0][each] = lamvec[each] + 2 * (rand_vals[each] - 0.5) * lamerr[each]

            b = np.zeros((1, len(cntvec)))
            for each in range(len(cntvec)):
                b[0][each] = (cntvec[each] / (self.efficiency * self.fissions) +
                           2 * (rand_vals[-1] - 0.5) * berr[0][each])

            # Build A mat
            A = np.zeros((len(self.times), self.numgroups))
            for index, t in enumerate(self.times):
                for coldex, lam in enumerate(use_lam[0]):
                    A[index, coldex] = lam * np.exp(-lam * t)
            # Solve LLS problem
            x, res, rank, s = np.linalg.lstsq(A, b.T, rcond=None)
            #x, res = scp.nnls(A, b, maxiter=None)

            # Check if valid LLS solution
            #if not np.all(x >= 0):
            #    continue
            
            # Record each abundance value in a dictionary
            for ai_ind, ai in enumerate(x):
                name = naming_convention + str(ai_ind + 1)
                stored_abundances[name].append(ai[0])
        # Determine mean value and uncertainty (depends on spread of data)
        for each in stored_abundances.keys():
            histval = stored_abundances[each]
            n, bins, patches = plt.hist(histval, bins=int(tot_iters/100), density=False)

            # Calculate 1 stnd dev
            n_sum = 0
            done = False
            for ind, eachval in enumerate(n):
                n_sum += eachval / len(stored_abundances[each])
                if n_sum >= 0.5 + 0.682689492/2 and not done:
                    # Set up this way to collect left half and right part of std
                    ai_err = abs(stored_abundances[each][0] - bins[ind + 1])
                    done = True
            if done == False:
                ai_err = abs(stored_abundances[each][0] - bins[-1])
            # https://stackoverflow.com/questions/50786699/
            #   how-to-calculate-the-standard-deviation-from-a-histogram-python-matplotlib
            mids = 0.5*(bins[1:] + bins[:-1])
            #dev = np.sqrt(np.average((mids - stored_abundances[each][0])**2, weights=n))
            dev = statistics.stdev(stored_abundances[each])
            #plt.vlines(stored_abundances[each][0] + ai_err, 0, max(n), linestyle='dotted', label='1\u03C3', color='red')
##            plt.axvspan(stored_abundances[each][0] - ai_err,
##                        stored_abundances[each][0] + ai_err,
##                        0, max(n), label='1\u03C3', color='red',
##                        alpha=0.25)
            # Print what percentage is in range
##            counter = 0
##            for ind, subeach in enumerate(histval):
##                if subeach > stored_abundances[each][0] - ai_err and \
##                   subeach < stored_abundances[each][0] + ai_err:
##                    counter += 1
##            pcnt_within = counter / len(stored_abundances[each]) * 100
##            print(f'Sigma within: {pcnt_within}%')
            plt.vlines(stored_abundances[each][0], 0, max(n), linestyle='dotted', label='Solution',
                       color='green')
            #plt.vlines(np.mean(histval), 0, max(n), linestyle='dotted', label='Mean', color='black')
            #plt.vlines(dev, 0, max(n), linestyle='dotted', label='Mean 1\u03C3', color='yellow')
            plt.axvspan(stored_abundances[each][0] - dev,
                        stored_abundances[each][0] + dev,
                        0, max(n), label='\u03C3', color='yellow',
                        alpha=0.25)
            # Print what percentage is in range
##            counter = 0
##            for ind, subeach in enumerate(histval):
##                if subeach > stored_abundances[each][0] - dev and \
##                   subeach < stored_abundances[each][0] + dev:
##                    counter += 1
##            pcnt_within = counter / len(stored_abundances[each]) * 100
            print(f'{stored_abundances[each][0]} +/- {dev}')
            
            

            plt.xlabel('Group Yield')
            plt.ylabel('Frequency')
            plt.legend()
            plt.title(f'{each}')
            #plt.show()
            plt.savefig(settings.imdir + f'g{each}-yield-MC.png')
            plt.close()
            abundance_errors.append(dev) #ai_err
            totyield += stored_abundances[each][0]
            totyield_err += dev ** 2
        totyield_err = np.sqrt(totyield_err)
        print(f'Yield: {totyield} +/- {totyield_err}')
            
        
        return abundance_errors

    def calc_cycle(self, vector, minlam, maxlam, dlam):
        """
        Iterates through all combinations for vector

        Returns
        -------
        vector : arraylike
            next valid combination
        complete : bool
            if all values are maxlam
        """
        complete = False
        for ind, val in enumerate(vector):
            if np.all(np.isclose(vector, maxlam)):
                complete = True
                break
            elif val < maxlam:
                vector[ind] += dlam
                break
            elif val >= maxlam:
                # Cycle subsequent max, increment next nonmax then stop
                for subind, subval in enumerate(vector):
                    if subval >= maxlam:
                        vector[subind] = minlam
                    else:
                        vector[subind] += dlam
                        break
                break
            else:
                print('Unknown Outcome')
                raise Exception
        return vector, complete

    def lami_solver(self,
                    lambda_values,
                    max_iters,
                    cntvec,
                    cnterr,
                    fiserr):
        """
        Solves LLS while assuming known decay constants.
            Performs a simplistic search for optimal decay constants.


        Returns
        -------
        abund_vec : 1D numpy array
            Abundance values
        lami_vec : 1D numpy array
            Best fit decay constants
        """
        print(np.log(2) / lambda_values)
        A = np.zeros((len(self.times), self.numgroups))
        b = np.zeros(len(self.times))
        x = np.zeros(self.numgroups)
        lambda_vector = None
        minimum_res = np.inf
        min_pcnt_diff = np.inf
        min_chi_res = np.inf
        abund_vec = None
        lami_vec = None
        cur_iter = 0
        n = 1
        print(f'Performing linear least squares {max_iters} times')
        if np.shape(np.shape(lambda_values))[0] == 2:
            run_type = '2d'
        elif np.shape(np.shape(lambda_values))[0] == 1:
            run_type = '1d'
            #comb_replace_index(lambda_vector, minlam, maxlam, dlam, self.numgroups)
        elif np.shape(np.shape(lambda_values))[0] == 0:
            run_type = '0d'
        else:
            print('Bad shape {np.shape(np.shape(lambda_values))[0]}')
            raise Exception
        lami_begin_time = time.time()
        for index, t in enumerate(self.times):
            b[index] = self.counts[index] / (self.fissions * self.efficiency)
        counts = b * (self.fissions * self.efficiency)
        berr = np.zeros(len(cnterr))
        for errind, err in enumerate(cnterr):
            nd_err = (err * 1 / (self.efficiency * self.fissions))**2
            fs_err = (fiserr * (cntvec[errind] / (self.efficiency * self.fissions**2)))**2
            berr[errind] = np.sqrt(nd_err + fs_err)
        # Starting from min val, go through every combination of decay constants
        complete = False
        while not complete:
            cur_iter += 1
            if cur_iter/max_iters >= 0.1 * n:
                n += 1
                net_time = time.time() - lami_begin_time
                full_complete_time = net_time / (0.1 * (n-1))
                print(f'    Progress: {round(cur_iter/max_iters * 100)}% in {round(net_time, 0)}s')
                print(f'        Estimated completion in {round(full_complete_time - net_time, 0)}s')
            #lambda_vector, complete = self.calc_cycle(lambda_vector, minlam, maxlam, dlam)
            if run_type == '2d':
                lambda_vector, complete = comb_replace_index_2d(lambda_vector, lambda_values, self.numgroups)
            elif run_type == '1d':
                lambda_vector, complete = comb_replace_index(lambda_vector, lambda_values, self.numgroups)
                #comb_replace_index(lambda_vector, minlam, maxlam, dlam, self.numgroups)
            elif run_type == '0d':
                lambda_vector = np.log(2) / lambda_values * np.random.rand(self.numgroups, 1)
                lambda_vector = np.log(2) / lambda_vector
                if cur_iter == max_iters:
                    complete = True
            else:
                print('Bad shape {np.shape(np.shape(lambda_values))[0]}')
                raise Exception
            #print(lambda_vector)
            #print(lambda_vector, complete)
            # Build A mat
            for index, t in enumerate(self.times):
                for coldex, lam in enumerate(lambda_vector):
                    if settings.irradiation == 'pulse':
                        A[index, coldex] = lam * np.exp(-lam * t)
                    elif settings.irradiation == 'infinite':
                        A[index, coldex] = np.exp(-lam * t)
            # Solve LLS problem
            #x, res, rank, s = np.linalg.lstsq(A, b, rcond=None)
            x, res = scp.nnls(A, b, maxiter=None) #Non-negative only
            if np.any(x < 0):
                continue
            #if res.size == 0:
            #    res = np.linalg.norm(A @ x - b)
            # CALCULATING % DIFF AND USING THAT ####################
            #if res < minimum_res and np.all(x >= 0):
                #minimum_res = min(minimum_res, res)

            
            keepin_response = keepin.KEEPIN()
            soln_vec = list()
            for each in range(len(x)):
                soln_vec.append(x.copy()[each] * lambda_vector.copy()[each])
                soln_vec.append(lambda_vector.copy()[each])
            group_counts, group_errs = keepin_response.simulate_lin_solve(self.times,
                                                                          soln_vec,
                                                                          self.fissions,
                                                                          self.efficiency)
##            # Instead of using pcnt diff, try Minkowski
##            avg_pcnt_diff = scipy.spatial.distance.minkowski(counts, group_counts, p=10)

            # Instead of using pcnt diff, try chi-squared
            group_counts = np.array(group_counts)
            chi_square = np.zeros(len(group_counts))
            for ind, grpcnt in enumerate(group_counts):
                chi_square[ind] = ((grpcnt - counts[ind]) / berr[ind])**2
            chi_square_res = sum(chi_square)
            
            prcnt_diff = list()
            for index in range(len(counts)):
                prcnt_diff.append((abs(counts[index] - group_counts[index]) /
                                   counts[index]))
            avg_pcnt_diff = np.mean(prcnt_diff)

            


##            if avg_pcnt_diff < min_pcnt_diff and np.all(x >= 0):
            if chi_square_res < min_chi_res and np.all(x >= 0):



                
            #######################################
            
                min_pcnt_diff = avg_pcnt_diff
                min_chi_res = chi_square_res
                lami_vec = lambda_vector.copy()
                abund_vec = x.copy()
##                print(f'Chi Square: {chi_square}')
##                print(f'% Diff: {min_pcnt_diff}')
                #print(minimum_res)
                #print(abund_vec)
                #input(lami_vec)
                #print(np.log(2) / lami_vec)
                #print(res)
        covariance_mat = np.linalg.inv(A.T @ A)
        print(f'Residual: {res}')
        print(f'% Diff: {min_pcnt_diff*100}')
        print(f'Chi Square: {min_chi_res}')
        print(f'Abund: {abund_vec}')
        print(f'Lami: {lami_vec}')
        return abund_vec, lami_vec, covariance_mat, res


    def abund_refine(self,
                     lambda_vector,
                     times,
                     counts):
        """
        Refine abundances by applying the selected decay vector
            with a larger data pool
        """
        A = np.zeros((len(times), self.numgroups))
        b = np.zeros((len(times), 1))
        x = np.zeros((self.numgroups, 1))
        abund_vec = None
        lami_vec = None
        for index, cnt in enumerate(counts):
            b[index] = cnt / (self.fissions * self.efficiency)
        # Build A mat
        for index, t in enumerate(times):
            for coldex, lam in enumerate(lambda_vector):
                A[index, coldex] = lam * np.exp(-lam * t)
        # Solve LLS problem
        x, res, rank, s = np.linalg.lstsq(A, b, rcond=None)
        if res.size == 0:
            res = np.linalg.norm(A @ x - b)
        lami_vec = lambda_vector.copy()
        abund_vec = x.copy()
        print(f'Residual: {res[0]}')
        return abund_vec, lami_vec

        

    def solver(self):
        """
        Solves the linear least squares problem

        Parameters
        ----------
        None

        Returns
        -------
        soln_vec : 1D numpy array
            12x1 list of a_i, lam_i values
        """
        # Build matrices
        n = len(self.times)
        ngroups = self.numgroups
        nvals = 2 * ngroups
        A = np.zeros((n, nvals))
        b = np.zeros((n, 1))
        x = np.zeros((nvals, 1))
        soln_vec = np.zeros((nvals, 1))
        for index, t in enumerate(self.times):
            A[index, :] = [1, -t] * ngroups #, 1, -t, 1, -t, 1, -t, 1, -t, 1, -t]
            b[index] = np.log(self.counts[index])
        x, res, rank, s = np.linalg.lstsq(A, b, rcond=None)
        for index, each in enumerate(x):
            if index % 2 == 0:
                soln_vec[index] = np.exp(x[index])
            else:
                soln_vec[index] = x[index]
        self.data_out(soln_vec)
        #print(np.linalg.lstsq(A, b))
        #print(f'Residual: {res}')
        return soln_vec

    def data_out(self, soln_vec):
        ai = list()
        lami = list()
        for ind, each in enumerate(soln_vec):
            if ind % 2 == 0:
                ai.append(each / soln_vec[ind+1] /
                          self.fissions / self.efficiency)
            else:
                lami.append(np.log(2) / each)
        print(f'n/F: {sum(ai)}')
        print(f'Half Lives: {lami}')
        print(f'Abundances: {ai}')
        print(f'Rel abundances: {ai / sum(ai)}')
        return

def debug_run_main(deb_group):
    data_name = 'test_' + str(deb_group)
    keepin_response = keepin.KEEPIN(data_name, 'fast')
    counts = keepin_response.simulate_instant(times, fissions, efficiency) 
    lin_solver = LLS(times, counts, fissions, efficiency, deb_group)
    soln_vec = lin_solver.solver()
    new_counts = keepin_response.simulate_lin_solve(times, soln_vec, fissions, efficiency)
    plt.plot(times, counts, label = f'{deb_group} Keepin')
    plt.plot(times, new_counts, label = f'{deb_group} LLS')
    plt.yscale('log')
    plt.ylabel('Delayed Neutron Count Rate [#/s]')
    plt.xlabel('Time [s]')
    plt.legend()
    plt.show()
    plt.close()
    return

def debug_run_reduced(deb_group, lambda_combos):
    data_name = 'test_' + str(deb_group)
    keepin_response = keepin.KEEPIN(data_name)
    counts = keepin_response.simulate_instant(times, fissions, efficiency) 
    lin_solver = LLS(times, counts, fissions, efficiency, deb_group)
    abund, lami = lin_solver.lami_solver()
    soln_vec = list()
    for each in range(len(abund)):
        soln_vec.append(abund[each] * lami[each])
        soln_vec.append(lami[each])
    new_counts = keepin_response.simulate_lin_solve(times, soln_vec, fissions, efficiency)
    plt.plot(times, counts, label = f'{deb_group} Keepin')
    plt.plot(times, new_counts, label = f'{deb_group} LLS')
##    soln_vec = list()
##    abund = fissions * efficiency * np.array([0.6, 0.3])
##    for each in range(len(abund)):
##        soln_vec.append(abund[each] * lami[each])
##        soln_vec.append(lami[each])
##    new_counts = keepin_response.simulate_lin_solve(times, soln_vec, fissions, efficiency)
##    plt.plot(times, new_counts, label = f'{deb_group} True? LLS')
    plt.yscale('log')
    plt.ylabel('Delayed Neutron Count Rate [#/s]')
    plt.xlabel('Time [s]')
    plt.legend()
    plt.show()
    plt.close()
    print(f'\nn/F: {sum(abund)}')
    print(f'Lams: {lami}')
    print(f'Half Lifes: {np.log(2) / lami}')
    print(f'Abundances: {abund}')
    print()
    return

def comb_replace_index(lam_vec, lambda_values, group_count):
    """
    Generate the combinations_with_replacement value on the fly
        for a given previous lam_vec using
        
    Parameters
    ----------
    lam_vec : vector
        Previously implemented lambda values
    lambda_values : vector
        Possible values of lambda for each group
    group_count : int
        Number of DNP groups to build

    Returns
    -------
    lam_vec : vector
        Next iteration of vector
    """
    complete = False
    if lam_vec is None:
        lam_vec = minlam * np.ones((1, group_count))[0]
    else:
        # If all are same number, increase first and reset others or finish if max
        maxlam = max(lambda_values)
        if np.all(np.isclose(lam_vec, lam_vec[0])) and np.isclose(lam_vec[0], maxlam):
            complete = True
        elif np.all(np.isclose(lam_vec, lam_vec[0])):
            lam_vec[0] = lambda_values[np.where(np.isclose(lambda_values, lam_vec[0]))[0][0] + 1]#dlam
            #lam_vec[0] += dlam
            lam_vec[1:] = minlam
        # If first is largest value in array, increase smallest value (closest first if tie)
        elif np.isclose(lam_vec[0], np.amax(lam_vec)):
            #lam_vec[np.argmin(lam_vec)] += dlam
            lam_vec[1:], trash = comb_replace_index(lam_vec[1:], lambda_values, group_count)
            #lam_vec[1:], trash = comb_replace_index(lam_vec[1:], minlam, maxlam, dlam, groups)
    if np.all(np.isclose(lam_vec, lam_vec[0])) and np.isclose(lam_vec[0], maxlam):
        complete = True
    #print(lam_vec)
    return lam_vec, complete

def comb_replace_index_2d(lam_vec, lambda_values, group_count):
    """
    Generate the combinations_with_replacement value on the fly
        for a given previous lam_vec using

    Parameters
    ----------
    lam_vec : vector
        Previously implemented lambda values
    lambda_values : matrix
        Possible values of lambda for each group
    group_count : int
        Number of DNP groups to build

    Returns
    -------
    lam_vec : vector
        Next iteration of vector
        
    """
    #print(f'inlet vec: {lam_vec}')
    complete = False
    #input(lam_vec)
    if lam_vec is None:
        lam_vec = lambda_values[:, 0].copy()
    else:
        val = lam_vec[0]
        if np.all(np.isclose(lam_vec, lambda_values[:, -1])):
            complete = True
        elif val < lambda_values[0, -1]:
            cur_col = np.where(np.isclose(val, lambda_values[0]))[0][0]
            lam_vec[0] = lambda_values[0, cur_col + 1]
        elif val >= lambda_values[0, -1]:
            #print(f'{val} >= {lambda_values[0, -1]}')
            # Cycle subsequent max, increment next nonmax then stop
            lam_vec[0] = lambda_values[0, 0].copy()
            comb_replace_index_2d(lam_vec[1:], lambda_values[1:], group_count)
        else:
            print('Unknown Outcome')
            raise Exception
    return lam_vec, complete


def generic_MC_lstsq_err(A, A_err, b, b_err, tot_iters = 1000):
    """
    Generates the solution distribution given a least squares problem.
        Uses a stochastic approach.

    Parameters
    ----------
    A : numpy array
        2D array of coefficients
    A_err : numpy array
        2D array same shape as A with uncertainties
    b : numpy array
        1D numpy vector of observed values
    b_err : numpy array
        1D numpy vector same shape as b with uncertainties
    tot_iters : int
        Number of iterations to perform

    Returns
    -------
    x : numpy array
        Parameter 1D vector of solutions to least squares problem
    x_errs : numpy_array
        Uncertainties of x
    
    """

    #print(f'Performing abundance error search {tot_iters} times')
    import spectra_handler
    spectra_class = spectra_handler.SPECTRA(None, None)
    stored_abundances = dict()
    totyield = 0
    totyield_err = 0
    n = 1
    lami_begin_time = time.time()
    for cur_iter in range(tot_iters):
        use_A = A.copy()
        use_b = b.copy()
        if cur_iter == 0:
            # Base solution
            #x_true, res_true, rank, s = spectra_class.lsqnonneg(A, b)
            x_true, res_true, res_many = spectra_class.lsqnonneg(A, b)
            x_errs = list()
            for xi_ind, xi in enumerate(x_true):
                name = 'x' + str(xi_ind + 1)
                stored_abundances[name] = list()
                stored_abundances[name].append(xi)

        #if cur_iter/tot_iters >= 0.1 * n:
        #    n += 1
        #    net_time = time.time() - lami_begin_time
        #    full_complete_time = net_time / (0.1 * (n-1))
        #    print(f'    Progress: {round(cur_iter/tot_iters * 100)}% in {round(net_time, 0)}s')
        #    print(f'        Estimated completion in {round(full_complete_time - net_time, 0)}s')
        # Generate random values (each A col and b)
        num_rand_vals = len(A[0, :]) + 1
        rand_vals = np.random.rand(num_rand_vals)
        
        # Assign decay constants and counts based on uncertainties and
        #   random values
        use_lam = np.zeros((1, len(A[0, :])))
        for row in range(len(A[:, 0])):
            for col in range(len(A[0, :])):
                use_A[row, col] = A[row, col] + 2 * (rand_vals[col] - 0.5) * A_err[row, col]

        for each in range(len(b)):
            use_b[each] = (b[each] + 2 * (rand_vals[-1] - 0.5) * b_err[each])

        # Solve LLS problem
        x, r_sq, res = spectra_class.lsqnonneg(use_A, use_b)
        
        # Record each abundance value in a dictionary
        for xi_ind, xi in enumerate(x):
            name = 'x' + str(xi_ind + 1)
            stored_abundances[name].append(xi)
    # Determine mean value and uncertainty (depends on spread of data)
    for each in stored_abundances.keys():
        histval = stored_abundances[each]
        n, bins, patches = plt.hist(histval, bins=int(tot_iters/100), density=False)

        # Calculate 1 stnd dev
        n_sum = 0
        done = False
        for ind, eachval in enumerate(n):
            n_sum += eachval / len(stored_abundances[each])
            if n_sum >= 0.5 + 0.682689492/2 and not done:
                # Set up this way to collect left half and right part of std
                ai_err = abs(stored_abundances[each][0] - bins[ind + 1])
                done = True
        if done == False:
            ai_err = abs(stored_abundances[each][0] - bins[-1])
        # https://stackoverflow.com/questions/50786699/
        #   how-to-calculate-the-standard-deviation-from-a-histogram-python-matplotlib
        mids = 0.5*(bins[1:] + bins[:-1])
       #dev = np.sqrt(np.average((mids - stored_abundances[each][0])**2, weights=n))
        dev = statistics.stdev(stored_abundances[each])
        plt.vlines(stored_abundances[each][0], 0, max(n), linestyle='dotted', label='Solution',
                   color='green')
        plt.axvspan(stored_abundances[each][0] - dev,
                    stored_abundances[each][0] + dev,
                    0, max(n), label='\u03C3', color='yellow',
                    alpha=0.25)
        #print(f'{stored_abundances[each][0]} +/- {dev}')
        
        

        #plt.xlabel('Parameter Solution')
        #plt.ylabel('Frequency')
        #plt.legend()
        #plt.title(f'{each}')
        #plt.show()
        #plt.savefig(settings.imdir + f'g{each}-yield-MC.png')
        plt.close()
        x_errs.append(dev) # MC err

        totyield += stored_abundances[each][0]
        totyield_err += dev ** 2
    totyield_err = np.sqrt(totyield_err)
    x_errs = np.asarray(x_errs)
    #print(f'Sum: {totyield} +/- {totyield_err}')

    # Fit Error
    #x_errs += max(np.abs(A@x_true - b))

    #temp_time = range(len(b))
    #errors_fit = A_err@x + A@x_errs
    #plt.plot(temp_time, A@x_true, label='Fit')
    #plt.fill_between(temp_time, A@x_true+errors_fit, A@x_true-errors_fit, alpha=0.25)
    #plt.plot(temp_time, b, label='Data')
    #plt.fill_between(temp_time, b+b_err, b-b_err, alpha=0.25)
    #plt.legend()
    #plt.show()

    return x_true, x_errs


if __name__ == '__main__':
    begin = time.time()
    dt = 0.1
    tmax = 33 #330
    times = np.arange(0, tmax+dt, dt)
    #fissions = 1#1E16
    #efficiency = 1#5.874643361662476e-08
    fissions = 1.013343827616795e+16
    volume = 0.1583105694
    efficiency = 1.650637878787879e-07
    default = True
    debug = False
    numgroups = 6
    debug_group = 2
    run_type = '2d' #'1d'
    percent_variance = 0.3 # 30%

    #A = np.array([[1, 0], [1, 0], [0, 1]])
    A = np.array([[2], [1], [1]])
    #b = np.array([2, 1, 1])
    b = np.array([1, 1, 1])
    #A_err = np.array([[0, 0], [0, 0], [0, 0]])
    A_err = np.zeros(np.shape(A))
    #b_err = np.array([0, 0, 0])
    b_err = np.zeros(np.shape(b))

    generic_MC_lstsq_err(A, A_err, b, b_err, tot_iters = 10000)
    quit()


    name = '6keepin235fast'
    keepin_response = keepin.KEEPIN(name)
    counts = keepin_response.simulate_instant(times, fissions, efficiency) 
    lin_solver = LLS(times, counts, fissions, efficiency, numgroups)
    lamvec = np.log(2) / ([54.51, 21.84, 6.00, 2.23, 0.496, 0.179])
    halerr = [0.94, 0.54, 0.17, 0.06, 0.029, 0.017]
    lamerr = list(np.log(2) / np.array([54.51, 21.84, 6.00, 2.23, 0.496, 0.179])**2 * np.array(halerr))
    #lamerr = [0.0025431927373323992, 0.006347501653479353, 0.024159887785289134,
    #          0.07606762112101241, 0.3549140709472327, 1.2049494664231994]
    fiss_err = 0#6.464567E13
    cnterr = None
    lin_solver.MC_abun_err(lamvec,
                 lamerr,
                 counts[0],
                 cnterr,
                 fiserr=fiss_err,
                    tot_iters=500)


##    print(f'Times 0 - {tmax}; steps {dt}')
##    if run_type == '2d':
##        avg_halflives = np.array([54.51, 21.84, 6.00, 2.23, 0.496, 0.179])
##        #min_halflives = [0.01, 0.2, 1.0, 5.0, 10, 40 ]
##        #max_halflives = [0.20, 1.0, 5.0, 10 , 40, 100]
##        num_nodes = 3
##        min_halflives = avg_halflives * (1-percent_variance)
##        max_halflives = avg_halflives * (1+percent_variance)
##        halflife_values = np.zeros((numgroups, num_nodes))
##        for group in range(numgroups):
##            halflife_values[group, :] = np.linspace(min_halflives[group], max_halflives[group], num_nodes)
##        print(f'Half lives: {halflife_values}')
##        lambda_values = np.log(2) / halflife_values
##        lambda_values.sort()
##        print(lambda_values)
##    elif run_type == '1d':
##        minhalf = 0.1
##        maxhalf = 60
##        half_divis = 5
##        dhalf = (maxhalf-minhalf)/(half_divis-1)
##        halflife_values = np.arange(minhalf, maxhalf+dhalf*0.1, dhalf)
##
##        minlam = np.log(2)/maxhalf
##        maxlam = np.log(2)/minhalf
##    
##        dlam = (maxlam-minlam)/(lam_divis-1)
##        lambda_values = np.sort(np.log(2) / halflife_values)
##        #lambda_values = np.sort(np.log(2) / np.array([54.51, 21.84, 6.00, 2.23, 0.496, 0.179]))
##        #halflife_values = np.log(2) / lambda_values

    


    
    # Start coarse, find best region, set bounds such that previous best is
        #avg iterate
    ######################
##    lin_solver = LLS(0, 0, fissions, efficiency, numgroups)
##    vec = None
##    potential_vals = np.array([[0, 0.1, 0.3],
##                               [1, 2, 3],
##                               [10, 20, 80]])
##    complete = False
##    while not complete:
##        #vec, complete = lin_solver.calc_cycle(vec, 0, 2, 1)
##        vec, complete = comb_replace_index_2d(vec, potential_vals, 3)
##        input(vec)
##        
##    input()
    ###########################
##    
##    if debug:
##        #lambda_combos = [x for x in itertools.combinations_with_replacement(lambda_values, r=debug_group)]
##        if run_type == '1d':
##            n = len(lambda_values) + debug_group - 1
##            k = debug_group
##            max_iters = math.factorial(n) / (math.factorial(k) * math.factorial(n-k))
##        elif run_type == '2d':
##            max_iters = num_nodes ** debug_group
##        print(f'Total iterations: {int(max_iters)}')
##        #print(f' {((maxlam-minlam)/dlam + 1) ** debug_group} iterations\n')
##        #debug_run_reduced(debug_group, lambda_values)
##        #debug_run_main(debug_group)
##
##
##
##    if default:
##        #lambda_combos = [x for x in itertools.combinations_with_replacement(lambda_values, r=numgroups)]
##        if run_type == '1d':
##            n = len(lambda_values) + numgroups - 1
##            k = numgroups
##            max_iters = math.factorial(n) / (math.factorial(k) * math.factorial(n-k))
##        elif run_type == '2d':
##            max_iters = num_nodes ** numgroups
##        print(f'Total iterations: {int(max_iters)}')
##        name = '6keepin235fast'
##        keepin_response = keepin.KEEPIN(name)
##        counts = keepin_response.simulate_instant(times, fissions, efficiency)
##        lin_solver = LLS(times, counts, fissions, efficiency, numgroups)
##        abund, lami = lin_solver.lami_solver(lambda_values,
##                                             max_iters)
##        soln_vec = list()
##        for each in range(len(abund)):
##            soln_vec.append(abund[each] * lami[each])
##            soln_vec.append(lami[each])
##        new_counts = keepin_response.simulate_lin_solve(times, soln_vec, fissions, efficiency)
##        #soln_vec = lin_solver.solver()
##        #new_counts = keepin_response.simulate_lin_solve(times, soln_vec, fissions, efficiency)
##        plt.plot(times, counts, label = 'Keepin')
##        #plt.plot(keepin_response.true_data_time, keepin_response.true_data_resp,
##        #         label='True', linestyle='', marker='.')
##        plt.plot(times, new_counts, label = 'LLS')
##        plt.yscale('log')
##        plt.ylabel('Delayed Neutron Count Rate [#/s]')
##        plt.xlabel('Time [s]')
##        plt.legend()
##        plt.savefig('images/6group_fit_run.png')
##        print(f'\nn/F: {sum(abund)}')
##        print(f'Lams: {lami}')
##        print(f'Half Lifes: {np.log(2) / lami}')
##        print(f'Abundances: {abund}')
##        print()
##    end = time.time()
##    print(f'Finished in {end - begin}s')
##    print(f'Completed {max_iters/(end-begin)} iterations per second')
