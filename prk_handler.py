import numpy as np
import matplotlib.pyplot as plt



class GPRK:
    """
    This class handles the point reactor kinetic equations.
        There are two primary cases considered:
            A set group input and individual precursors input
        This class handles the group input variant
        The PRKEs are solved one time step at a time.
        To determine the value at the next time step, linear
            intrpolation is used. For sufficiently small
            time steps, this method is fine.
    """

    def __init__(self,
                 lambda_vector,
                 abundance_vector,
                 times,
                 react_magnitude,
                 lambda_errs=None,
                 abundance_errs=None,
                 gen_time=5.56122E-9,
                 nubar=2.60340,
                 n0=1,
                 reactivity_type='step',
                 nubar_err=2.16512E-5):
        """
        Initialize the class

        Parameters
        ----------
        lambda_vector : array-like
            1D vector array of decay constants
        abundance_vector : array-like
            1D vector array of group yields
        times : array-like
            1D vector array of time values
        react_magnitude : float
            How many dollars of reactivity to insert (1 == 1$)
        gen_time : float
            Neutron generation time
        nubar : float
            Average number of neutrons produced per fission event
        n0 : float
            Initial neutron density solution to PRKEs
        reactivity_type : str
            Type of reactivity insertion/withdrawal to simulate

        Returns
        -------
        None
        """
        self.lambda_vector = lambda_vector
        self.beta_i = abundance_vector / nubar
        self.beta = sum(self.beta_i)
        self.n0 = n0
        self.times = times
        self.gen_time = gen_time
        self.r_type = reactivity_type
        self.r_mag = react_magnitude
        self.errors = False
        if type(lambda_errs) != type(None) and \
           type(abundance_errs) != type(None):
            self.lam_err = lambda_errs
            self.abu_err = abundance_errs
            self.errors = True
            self.beta_i_err = np.zeros((len(abundance_vector), 1))
            for errind, ai in enumerate(abundance_vector):
                delai = 1/nubar * self.abu_err[errind]
                delnubar = abundance_vector[errind] / nubar**2 * nubar_err
                self.beta_i_err[errind] = np.sqrt(delai**2 + delnubar**2)
            self.beta_err = np.linalg.norm(self.beta_i_err)
        return

    def reactivity(self, t):
        """
        Reactivity as a function of time

        Parameters
        ----------
        t : float
            Time value currently evaluating

        Returns
        -------
        r : float
            Reactivity at current time
        """
        if self.r_type == 'step':
            if t < 0:
                r = 0
            else:
                r = self.r_mag * self.beta
        

        return r


    
    def group_prk(self):
        """
        Builds and runs the PRK matrix problem
            for a set number of DNP groups

        Parameters
        ----------
        None

        Returns
        -------
        soln_matrix : matrix
            2D matrix; each row is for each variable, column for time step
        err_matrix : matrix
            2D matrix; each row is for each variable, column for time step
        """
        num_groups = len(self.lambda_vector)
        n = num_groups+1
        A = np.zeros((n, n))
        soln_matrix = np.zeros((n, len(self.times)))
        err_matrix = np.zeros((n, len(self.times)))
        prev_t = 0
        for row in range(n):
            if row == 0:
                for col in range(1, n):
                    A[row, col] = self.lambda_vector[col - 1]
            else:
                A[row, 0] = self.beta_i[row - 1] / self.gen_time
                for col in range(1, n):
                    if row == col:
                        A[row, col] = -self.lambda_vector[col - 1]
        for cur_col, t in enumerate(self.times):
            r = self.reactivity(t)
            A[0, 0] = (r - self.beta) / self.gen_time
            if t == 0:
                n0 = self.n0
                ci0 = list()
                for ind in range(n-1):
                    # Initial conc wrong?
                    prec_val = self.beta_i[ind] * n0 / (self.lambda_vector[ind] *
                                                   self.gen_time)
                    ci0.append(prec_val)
                soln_matrix[0, cur_col] = n0
                for cur_row, each in enumerate(ci0):
                    soln_matrix[cur_row+1, cur_col] = each

                if self.errors:
                    err_matrix[:, cur_col] = 0
            else:
                n0 += b[0] * (t - prev_t)
                #print(n0)
                ci0_new = list()
                for ind in range(n-1):
                    ci0_new.append(ci0[ind] + b[ind+1] * (t - prev_t))
                ci0 = ci0_new.copy()
                soln_matrix[0, cur_col] = n0
                for cur_row, each in enumerate(ci0):
                    soln_matrix[cur_row+1, cur_col] = each

                
                if self.errors:
                    dt = t - prev_t
                    # n err
                    dn = 0
                    #   n0
                    dn += (1 + (r - self.beta)/self.gen_time * dt +
                           dt**2/self.gen_time *
                           sum(self.beta_i * self.lambda_vector) *
                           err_matrix[0, cur_col-1])**2
                    #   beta
                    dn += (n0 * dt / self.gen_time * self.beta_err)**2
                    for ind in range(len(self.lambda_vector)):
                        #   lami
                        dn += (self.lam_err[ind] * (ci0[ind] * dt +
                                               self.beta_i[ind] / self.gen_time * n0 * dt**2 -
                                               2 * self.lambda_vector[ind] * ci0[ind] * dt**2)
                               )**2
                        #   betai
                        dn += (self.beta_i_err[ind] * (n0 * self.lambda_vector[ind] *
                                                  dt**2 / self.gen_time))**2
                        #   Ci0
                        dn += (err_matrix[ind+1, cur_col-1] * (self.lambda_vector[ind] *
                                                               dt -
                                                               self.lambda_vector[ind]**2 *
                                                               dt**2))**2
                    n_err = np.sqrt(dn)
                    err_matrix[0, cur_col] = n_err
                    for ind in range(len(self.lambda_vector)):
                        # Ci err
                        dci = 0
                        #   Ci0
                        dci += (err_matrix[ind+1, cur_col-1] * (1 -
                                                                self.lambda_vector[ind] *
                                                                dt))**2
                        #   betai
                        dci += (self.beta_i_err[ind] * (n0 * dt / self.gen_time))**2
                        # n0
                        dci += (err_matrix[0, cur_col-1] * (self.beta_i[ind] * dt /
                                                            self.gen_time))**2
                        # lami
                        dci += (self.lam_err[ind] * (dt * ci0[ind]))**2
                        ci_err = np.sqrt(dci)
                        err_matrix[ind+1, cur_col] = ci_err
                    
                    
            x_vec = np.zeros((n, 1))
            for ind in range(n):
                if ind == 0:
                    x_vec[ind] = n0
                else:
                    x_vec[ind] = ci0[ind-1]
            b = A @ x_vec
            prev_t = t

        return soln_matrix, err_matrix




if __name__ == '__main__':
    lambda_vector = np.array([0.01271596, 0.03081311, 0.11215974,
                              0.32044158, 1.35677102, 3.75954429])
    abunda_vector = np.array([0.00037125, 0.0027802,  0.00277704,
                              0.00778691, 0.00266239, 0.0015291 ])
    lam_err = np.array([0.00019073945529993002, 0.0004487346818842039,
                        0.001633394242058498, 0.004955282156812542,
                        0.0197588013152238, 0.05475064498520127])
    abu_err = np.array([1.5116938474450134e-11, 1.655471987831139e-11,
                        5.674897275809026e-12, 1.687695117506137e-12,
                        2.2820378889822704e-13, 3.3426026441962236e-14])
    times = np.arange(0, 20, 1E-3)
    #gen_time = 5.56122E-9

    gen_time = 1E-5
    reactivity_magnitude = 0.5
    group_solve = GPRK(lambda_vector, abunda_vector,
                       times, reactivity_magnitude,
                       gen_time=gen_time,
                       lambda_errs=lam_err,
                       abundance_errs=abu_err)
    soln_matrix, err_matrix = group_solve.group_prk()

    plt.plot(times, soln_matrix[0, :]/soln_matrix[0, :][0],
             label=f'Step Insertion of {reactivity_magnitude}$')
    plt.fill_between(times, soln_matrix[0, :] + err_matrix[0, :],
                     soln_matrix[0, :] - err_matrix[0, :],
                     alpha=0.5)
    plt.yscale('log')
    plt.yticks([0.1, 1, 10, 100])
    plt.xlabel('Time [s]')
    plt.ylabel('Relative Neutron Density')
    plt.legend()
    plt.show()
    plt.close()

    # Precursor plots
    for each in range(len(sol_matrix[:, 0])):
        plt.plot(times, soln_matrix[each+1, :],
                 label=f'Step Insertion of {reactivity_magnitude}$')
        plt.fill_between(times, soln_matrix[each+1, :] + err_matrix[each+1, :],
                         soln_matrix[each+1, :] - err_matrix[each+1, :],
                         alpha=0.5)
        plt.yscale('log')
        plt.yticks([0.1, 1, 10, 100])
        plt.xlabel('Time [s]')
        plt.ylabel(f'Group {each+1} Precursor Conc')
        plt.legend()
        plt.show()
        plt.close()
        


