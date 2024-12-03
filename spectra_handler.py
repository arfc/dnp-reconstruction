import settings
import numpy as np
import time
import scipy.optimize as scp
import matplotlib.pyplot as plt
from pandas import read_excel

class SPECTRA:
    """
    Handles spectral data from different sources,
        normalization, and adding spectral data
        into a given dictionary.
    """

    def __init__(self,
            energy_mesh,
            times):
        """
        Initializes the class

        Parameters
        ----------
        energy_mesh : numpy array
            1D midpoints of the energy bins

        Returns
        -------
        None
        """
        self.energy_mesh = energy_mesh
        self.times = times

        return


    def read_endf_excel(self, filename, sheetname):
        """
        Read the ENDF U-235 spectral 6 group results from Excel

        Parameters
        ----------
        filename : str
            Path to the Excel file as well as the file's name
        sheetname : str
            Name of the sheet with the data in Excel

        Returns
        -------
        endf_spectral_matrix : numpy array
            2D numpy array where each row is a unique group and each column
                is a different energy bin based on the class' energy bin

        """
        df = read_excel(filename,
                        sheet_name=sheetname)
        ev_group_ind = -1
        prob_group_ind = -1
        temp_store = dict()
        for key in df:
            if 'eV' in key:
                ev_group_ind = int(key.replace('g', '').replace('eV', '')) - 1
                energies = np.asarray(df[key])
                energy_list = list()
                for ind, each in enumerate(energies):
                    if str(each) == 'nan':
                        pass
                    else:
                        if '+' in each:
                            splitval = '+'
                        elif '-' in each:
                            splitval = '-'

                        data = str(each).split(splitval)
                        value = data[0] + 'E' + splitval + data[1]
                        energy_list.append(float(value) / 1E6)
                energies = np.asarray(energy_list)
            elif 'prob' in key:
                prob_group_ind = int(key.replace('g', '').replace('prob', '')) - 1
                probs = np.asarray(df[key])
                probs_list = list()
                for ind, each in enumerate(probs):
                    if str(each) == 'nan':
                        pass
                    else:
                        if '+' in each:
                            splitval = '+'
                        elif '-' in each:
                            splitval = '-'
                        data = str(each).split(splitval)
                        value = data[0] + 'E' + splitval + data[1]
                        probs_list.append(float(value) / 1E6)
                probs = np.asarray(probs_list)
            else:
                raise Exception('Unknown Column')



            if ev_group_ind == prob_group_ind:
                temp_store[ev_group_ind] = {'energies': energies, 'probs': probs}

        endf_spectral_matrix = np.zeros((ev_group_ind+1, len(self.energy_mesh)))

        # Convert ENDF spectral data to current energy mesh
        # Find where energies should be binned, then put probs in the correct bins
        for key in temp_store:
            bin_indeces = np.digitize(temp_store[key]['energies'], self.energy_mesh[:-1])
            for dictdex, bin_index in enumerate(bin_indeces):
                endf_spectral_matrix[key, bin_index] += temp_store[key]['probs'][dictdex]
 
        return endf_spectral_matrix


    def update_dict(self, feed_dict, ensdf_handler):
        """
        Given a dictionary with isotopes as keys (i.e. xe135),
            add normalized spectral data (probability normalization,
            or 1 neutron spectrum)

        Parameters
        ----------
        feed_dict : dict
            keys : str
                Name of isotope
            values : unknown
                unknown
        ensdf_handler : module
            Module by the same name

        Returns
        -------
        feed_dict : dict
            keys : str
                Name of isotope
            values : numpy array
                Same as previous, but adds the normalized spectral results
                    and bins
        valid_list : list
            List of isotopes which contain spectral data
        """
        valid_list = list()
        ensdf_gen = ensdf_handler.ENSDF(settings.ensdf_fname, settings.ensdf_sheet)
        ensdf_dict = ensdf_gen.parse_file()
        for isotope in feed_dict.keys():
            # Check each spectral library and handle appropriately
            if isotope in settings.iaea_spectra:
                bins, values = ensdf_handler.spectra_analysis(settings.iaea_spectra[isotope],
                                                              path='spectra/spectra/',
                                                              display=False)

                feed_dict[isotope]['spectrum_bins'] = bins
                feed_dict[isotope]['spectrum_values'] = values
                valid_list.append(isotope)

            elif isotope in settings.endf_spectra:
                # https://www-nds.iaea.org/exfor/endf.htm
                df = read_excel(settings.endf_spectra_filename,
                        sheet_name=isotope)
                for key in df:
                    if 'eV' in key:
                        energies = np.asarray(df[key])
                        energy_list = list()
                        for ind, each in enumerate(energies):
                            energy_list.append(float(each) / 1E6)
                        energies = np.asarray(energy_list)
                    elif 'prob' in key:
                        probs = np.asarray(df[key])
                        probs_list = list()
                        for ind, each in enumerate(probs):
                            probs_list.append(float(each) / 1E6)
                        probs = np.asarray(probs_list)
               
                feed_dict[isotope]['spectrum_bins'] = energies
                feed_dict[isotope]['spectrum_values'] = probs
                valid_list.append(isotope)


        for isotope in valid_list:
            meshed_spectra = np.zeros(len(self.energy_mesh))
            for eind, e in enumerate(self.energy_mesh):
                # spectra uses closest datapoint
                closest_bindex = (np.abs(feed_dict[isotope]['spectrum_bins'] -
                                      e)).argmin()
                meshed_spectra[eind] = (feed_dict[isotope]['spectrum_values'][closest_bindex])

            feed_dict[isotope]['spectrum_bins'] = self.energy_mesh
            # These spectra need to be normalize no matter what
            normalize = np.sum(np.abs(meshed_spectra))
            feed_dict[isotope]['spectrum_values'] = meshed_spectra/normalize

        return feed_dict, valid_list

    def least_squares_spectrum(self, group_lambdas, isotope_dict, valid_list):
        """
        Generates the group spectra by using least squares to determine
            the contribution of each isotope to each groups' yield.

        Parameters
        ----------
        group_lambdas : numpy array
            1D numpy array for each group the decay constant in seconds^-1
        isotope_dict : dict
            key : str
                Name of isotope (xe135)
            values : dict
                key : str
                    Name of category (needs to contain spectra data)
                values : varies
                    Lists of values and uncertainties for most
        valid_list : list
            List of strings of isotopes with spectra

        Returns
        -------
        group_spectra : numpy array
            Rows are each group spectra, columns are energy midpoint values
        """
        group_spectra = np.zeros((settings.fit_groups, len(self.energy_mesh)))
        normalizer = np.zeros(settings.fit_groups)
        for isotope in valid_list:
            single_iso = False
            A = np.zeros((len(self.times)+1, settings.fit_groups))
            b = np.zeros(len(self.times)+1)
            # Group decays should sum to isotope decay
            A[-1, :] = np.ones(settings.fit_groups)
            b[-1] = 1
            lami = np.log(2) / isotope_dict[isotope]['halflife'][0]
            conc = isotope_dict[isotope]['conc'][0]
            Pn   = isotope_dict[isotope]['emission'][0]
            spec = isotope_dict[isotope]['spectrum_values']
            
            for tind, t in enumerate(self.times):
                for gind in range(settings.fit_groups):
                    A[tind, gind] = group_lambdas[gind] * np.exp(-group_lambdas[gind] * t)
                b[tind] = lami * np.exp(-lami * t)
            

            # Depending on the isotope, can make certain columns 0
            #   in order to match lam < lam < lam setup with only m groups
            bigger_g = 5
            smaller_g = 0
            for gind, g in enumerate(group_lambdas):
                if g < lami:
                    smaller_g = max(smaller_g, gind)
                if g > lami:
                    bigger_g = min(bigger_g, gind)
            if bigger_g == 0:
                single_iso = True
            if smaller_g == 5:
                single_iso = True
            A_use = A[:, smaller_g:bigger_g+1]
            if not single_iso:
                #x, res, rank, s = np.linalg.lstsq(A_use, b.T, rcond=None)
                x, res = scp.nnls(A_use, b.T, maxiter=None)
            else:
                x = [1]
            frac_vec = np.zeros(settings.fit_groups)
            frac_vec[smaller_g] = x[0]
            if not single_iso:
                frac_vec[bigger_g] = x[1]

            # x vector is fraction of isotope to each group
            # Sum contribution to each group
            for gind in range(settings.fit_groups):
                spec_i_to_g = frac_vec[gind] * conc[0] * Pn * spec
                group_spectra[gind, :] += spec_i_to_g
                normalizer[gind] += frac_vec[gind] * conc[0] * Pn

        for gind in range(settings.fit_groups):
            group_spectra[gind, :] = group_spectra[gind, :] / normalizer[gind]
            


        return group_spectra

    def lsqnonneg(self, C, d, x0=None, tol=None, itmax_factor=3):
        '''Linear least squares with nonnegativity constraints.
        (x, resnorm, residual) = lsqnonneg(C,d) returns the vector x that minimizes norm(d-C*x)
        subject to x >= 0, C and d must be real

        A Python implementation of NNLS algorithm
        References:
        [1]  Lawson, C.L. and R.J. Hanson, Solving Least-Squares Problems, Prentice-Hall, Chapter 23, p. 161, 1974.
        Contributed by Klaus Schuch (schuch@igi.tugraz.at)
        based on MATLAB's lsqnonneg function
        https://gist.github.com/jdelafon/b7fdc7a0bae42af56366fc7786cc5d54
        '''

        eps = 2.22e-16    # from matlab
        def msize(x, dim):
            s = x.shape
            if dim >= len(s):
                return 1
            else:
                return s[dim]

        if tol is None:
            norm_c = abs(C).sum().max()
            tol = 10*eps*norm_c*(max(C.shape)+1)

        C = np.asarray(C)

        (m,n) = C.shape
        P = np.zeros(n)
        Z = np.arange(1, n+1)

        if x0 is None: x=P
        else:
            if any(x0 < 0): x=P
            else: x=x0

        ZZ = Z

        resid = (d - np.dot(C, x))
        #resid = (d - np.dot(C, x)) / d
        w = np.dot(C.T, resid)

        outeriter=0; it=0
        itmax=itmax_factor*n
        exitflag=1
        dont_break = True

        # outer loop to put variables into set to hold positive coefficients
        while np.any(Z) and np.any(w[ZZ-1] > tol) and dont_break:
            outeriter += 1

            t = w[ZZ-1].argmax()
            t = ZZ[t]

            P[t-1]=t
            Z[t-1]=0

            PP = np.where(P != 0)[0]+1
            ZZ = np.where(Z != 0)[0]+1

            CP = np.zeros(C.shape)

            CP[:, PP-1] = C[:, PP-1]
            CP[:, ZZ-1] = np.zeros((m, msize(ZZ, 1)))

            #z=np.dot(np.linalg.pinv(CP), d)
            D = np.diag(1/d)
            D_square = np.square(D)
            z = np.linalg.pinv(CP.T @ D_square @ CP) @ CP.T @ D_square @ d

            z[ZZ-1] = np.zeros((msize(ZZ,1), msize(ZZ,0)))

            # inner loop to remove elements from the positve set which no longer belong
            while np.any(z[PP-1] <= tol):
                it += 1

                if it > itmax:
                    max_error = z[PP-1].max()
                    #print(f'Iterations {it} exceeded, error up to {max_error}')
                    dont_break = False
                    #return (z, sum(resid*resid), resid)
                    #raise Exception(f'Exiting: Iteration count {it} exceeded\n  Try raising the tolerance tol. {max_error}')

                QQ = np.where((z <= tol) & (P != 0))[0]
                alpha = min(x[QQ]/(x[QQ] - z[QQ]))
                x = x + alpha*(z-x)

                ij = np.where((abs(x) < tol) & (P != 0))[0]+1
                Z[ij-1] = ij
                P[ij-1] = np.zeros(max(ij.shape))
                PP = np.where(P != 0)[0]+1
                ZZ = np.where(Z != 0)[0]+1

                CP[:, PP-1] = C[:, PP-1]
                CP[:, ZZ-1] = np.zeros((m, msize(ZZ, 1)))

                z=np.dot(np.linalg.pinv(CP), d)
                D = np.diag(1/d)
                D_square = np.square(D)
                #z = np.linalg.pinv(CP.T @ D_square @ CP) @ CP.T @ D_square @ d
                z[ZZ-1] = np.zeros((msize(ZZ,1), msize(ZZ,0)))

            x = z
            #resid = (d - np.dot(C, x))
            resid = (d - np.dot(C, x)) / d
            w = np.dot(C.T, resid)

        return (x, sum(resid * resid), resid)

    def pcnt_data_recon_lstsq(self, spectral_matrix, lam_vec, lam_err, abu_vec, abu_err, uncertainty=settings.spectra_uncertainty):
        """
        This method also uses least squares, but takes in data of the 
            spectra evolution over time in order to generate the results.
            Uses the least square percentage regression method.

        Parameters
        ----------
        spectral_matrix : numpy array
            2D matrix where columns are times, rows are energies, and values are counts
        lam_vec : numpy array
            1D array of group decay constants
        abu_vec : numpy array
            1D array of group yields
        uncertainty : Bool
            True to generate group spectra uncertainties

        Returns
        -------
        group_spectra : numpy array
            Rows are each group spectra, columns are energy midpoint values
        res : numpy array
            Residual at each term
        group_uncertainties : numpy array
            Uncertainties, same shape as group spectra
        """
        A = np.zeros((len(self.times), settings.fit_groups))
        n = 1
        cur_begin_time = time.time()
        group_spectra = np.zeros((settings.fit_groups, len(self.energy_mesh)))
        group_uncertainties = group_spectra.copy()
        A_err = A.copy()
        print('The spectral matrix does not have uncertainty considered yet due to lacking spectral uncertainty')
        if uncertainty:
            from linear_least_squares import generic_MC_lstsq_err
        for gind in range(settings.fit_groups):
            lam = lam_vec[gind]
            yld = abu_vec[gind]
            if settings.irradiation == 'pulse':
                a_val = lam * yld
            elif settings.irradiation == 'infinite':
                a_val = yld
            for tind, t in enumerate(self.times):
                A[tind, gind] = a_val * np.exp(-lam * t) * (settings.efficiency * settings.fissions)
                A_err[tind, gind] = (settings.efficiency * settings.fissions * (abu_err[gind] * lam * np.exp(-lam * t))**2 +
                        (lam_err[gind] * np.exp(-lam * t) * (yld - yld*t * lam))**2 )
        cur_begin_time = time.time()
        for eind, energy in enumerate(self.energy_mesh):
            if eind/len(self.energy_mesh) >= 0.1 * n:
                n += 1
                net_time = time.time() - cur_begin_time
                full_complete_time = net_time / (0.1 * (n-1))
                print(f'    Progress: {round(eind/len(self.energy_mesh) * 100)}% in {round(net_time, 0)}s')
                print(f'        Estimated completion in {round(full_complete_time - net_time, 0)}s')
            b = np.zeros(len(self.times))
            b_err = b.copy()
            for tind, t in enumerate(self.times):
                b[tind] = spectral_matrix[eind, tind] #/ (settings.efficiency * settings.fissions)
            x, res_sq, res = self.lsqnonneg(A, b)
            if uncertainty:
                _, spec_errs = generic_MC_lstsq_err(A, A_err, b, b_err, tot_iters = settings.spectra_iters)
                group_uncertainties[:, eind] = spec_errs

            group_spectra[:, eind] = x
            if np.any(x < 0):
                print('Negative probability in group spectra generation')
                print(f'Energy: {energy}')
                #print(f'x: {x}')
                #print(f'Counts: {b * (settings.efficiency * settings.fissions)}')
                #print(f'B mat: {b}')
                #print(f'Counts at 0s: {np.sum(spectral_matrix[:, 0])}')
                #input(f'A mat: {A}')

            #plt.plot(self.times, A@x, label='LLS')
            #plt.plot(self.times, b, label='Counts')
            #plt.xlabel('Time [s]')
            #plt.ylabel('Spectra')
            #plt.yscale('log')
            #plt.legend()
            #plt.title(f'{energy} [MeV]')
            #plt.tight_layout()
            #plt.savefig(f'./{energy}.png')
            #plt.show()
            #plt.close()

        return group_spectra, res, group_uncertainties


    def data_recon_lstsq_deprecated(self, spectral_matrix, lam_vec, abu_vec):
        """
        This method also uses least squares, but takes in data of the 
            spectra evolution over time in order to generate the results

        Parameters
        ----------
        spectral_matrix : numpy array
            2D matrix where columns are times, rows are energies, and values are counts
        lam_vec : numpy array
            1D array of group decay constants
        abu_vec : numpy array
            1D array of group yields

        Returns
        -------
        group_spectra : numpy array
            Rows are each group spectra, columns are energy midpoint values
        """
        #new_times = np.linspace(settings.start_time, settings.end_time, settings.fit_groups)
        #old_times = self.times.copy()
        #print(f'Using custom times in data_recon_lstsq {new_times}')
        #self.times = new_times
        A = np.zeros((len(self.times), settings.fit_groups))
        n = 1
        cur_begin_time = time.time()
        group_spectra = np.zeros((settings.fit_groups, len(self.energy_mesh)))
        for gind in range(settings.fit_groups):
            lam = lam_vec[gind]
            yld = abu_vec[gind]
            if settings.irradiation == 'pulse':
                a_val = lam * yld
            elif settings.irradiation == 'infinite':
                a_val = yld
            for tind, t in enumerate(self.times):
                A[tind, gind] = a_val * np.exp(-lam * t) * (settings.efficiency * settings.fissions)
        A_inv = np.linalg.pinv(A)
        for eind, energy in enumerate(self.energy_mesh):
            b = np.zeros(len(self.times))
            for tind, t in enumerate(self.times):
                b[tind] = spectral_matrix[eind, tind] #/ (settings.efficiency * settings.fissions)

            #x, res, rank, s = np.linalg.lstsq(A, b.T, rcond=None)
            x = A_inv @ b
            #x, res = scp.nnls(A, b)

            group_spectra[:, eind] = x
            if np.any(x < 0):
                print('Negative probability in group spectra generation')
                print(f'Energy: {energy}')
                #print(f'x: {x}')
                #print(f'Counts: {b * (settings.efficiency * settings.fissions)}')
                #print(f'B mat: {b}')
                #print(f'Counts at 0s: {np.sum(spectral_matrix[:, 0])}')
                #input(f'A mat: {A}')

            plt.plot(self.times, A@x, label='LLS')
            plt.plot(self.times, b, label='Counts')
            plt.xlabel('Time [s]')
            plt.ylabel('Spectra')
            plt.yscale('log')
            plt.legend()
            plt.title(f'{energy} [MeV]')
            plt.tight_layout()
            #plt.savefig(f'./{energy}.png')
            plt.show()
            plt.close()

        return group_spectra
    
    def spectral_matrix_constructor(self, ORIGEN_dict, ensdf_handler, alt_norm=False):
        """
        Build a spectral matrix where each row represents a new energy
            and each column is a new time

        Parameters
        ----------
        ORIGEN_dict : dict
            keys : str
                Name of isotopes
            values : dict
                keys : str
                    Name of categories
                values : numpy array
                    spectral data, atoms, etc.
        ensdf_handler : module
            Module with the same name

        Returns
        -------
        spectral_matrix : numpy array
            2D numpy array of probabilitiy of energy emission normalized
                at each time step
        valid_list : list
            List of strings; isotopes with spectra
        """
        ORIGEN_dict, valid_list = self.update_dict(ORIGEN_dict, ensdf_handler)
        print('Constructing spectral matrix')
        n = 1
        cur_begin_time = time.time()
        # Column is energy spectra for given energy
        spectral_matrix = np.zeros((len(self.energy_mesh), len(self.times)))
        if alt_norm:
            count_matrix = np.zeros((len(self.energy_mesh), len(self.times)))
        for tind, t in enumerate(self.times):
            if tind/len(self.times) >= 0.1 * n:
                n += 1
                net_time = time.time() - cur_begin_time
                full_complete_time = net_time / (0.1 * (n-1))
                print(f'    Progress: {round(tind/len(self.times) * 100)}% in {round(net_time, 0)}s')
                print(f'        Estimated completion in {round(full_complete_time - net_time, 0)}s')
            net_spectra = np.zeros(len(self.energy_mesh))
            # sum over each isotope
            for isotope in valid_list:
                atoms = ORIGEN_dict[isotope]['conc'][0]
                Pn = ORIGEN_dict[isotope]['emission'][0]
                lam = np.log(2) / ORIGEN_dict[isotope]['halflife'][0]
                activity = (settings.efficiency * ORIGEN_dict[isotope]['emission'][0] *
                            atoms[tind] * lam)
                for eind, e in enumerate(self.energy_mesh):
                    # spectra uses closest datapoint
                    closest_bindex = (np.abs(ORIGEN_dict[isotope]['spectrum_bins'] -
                                          e)).argmin()
                    net_spectra[eind] += (activity *
                                          ORIGEN_dict[isotope]['spectrum_values'][closest_bindex])


            if settings.spectra_normalized:
                normalize = np.sum(np.abs(net_spectra))
                name = 'Probability'
            else:
                normalize = 1
                name = 'Neutron Intensity'
            spectral_matrix[:, tind] = net_spectra / normalize
            if alt_norm:
                count_matrix[:, tind] = net_spectra
        if alt_norm:
            return spectral_matrix, valid_list, count_matrix

        return spectral_matrix, valid_list

