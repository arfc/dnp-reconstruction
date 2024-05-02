from pandas import read_excel
import numpy as np
import matplotlib.pyplot as plt
import settings


class ENSDF:
    """
    This script handles extraction of IAEA/ENSDF Evaluation data from .xlsx files
    """

    def __init__(self,
                 filename,
                 sheetname):
        """
        Initialize

        Parameters
        ----------
        filename : str
            Name of file containing ENSDF data
        sheetname : str
            Name of sheet in Excel

        Returns
        -------
        None
        """
        self.file = filename
        self.sheet = sheetname
        return

    def parse_file(self):
        """
        Goes through each file line by line,
            extracts isotope name, half life,
            and average delayed neutron emission quantity

        Parameters
        ----------
        None

        Returns
        -------
        ensdf_data : dict
            key : str
                Name of isotope (i.e. xe135)
            value : list
                List of half life (s) followed by average emission per decay
        """
        df = read_excel(self.file, sheet_name = self.sheet)
        skip_next_nuc = False
        ensdf_data = dict()
        element = None
        atomic_mass = None
        half_life = None
        half_present = False
        num_emit = None
        high_state = False
        low_state = False
        for line_index in range(3, len(df)):
            if skip_next_nuc:
                skip_next_nuc = False
                continue
            # Only record data when halflife is present
            half_emiss = list()
            # atomic mass eval
            nuc_val = df['Nuclide'][line_index]
            isomer = df['Isomer'][line_index]
            isomer_state = False
            low_state = False
            high_state = False
            if isomer == 'ISO':
                isomer_state = True
            elif isomer == 'HS':
                high_state = True
            elif isomer == 'LS':
                low_state = True
            
            if type(nuc_val) is int and not skip_next_nuc:
                atomic_mass = nuc_val
                if isomer_state or low_state:
                    skip_next_nuc = False
                else:
                    skip_next_nuc = True
            elif type(nuc_val) is int and skip_next_nuc:
                skip_next_nuc = False

            # element eval
            ele_guess = df['Unnamed: 1'][line_index]
            if type(ele_guess) is str:
                element = ele_guess
            elif isomer_state or low_state:
                element = df['Unnamed: 1'][line_index - 2]
            if isomer_state or low_state:
                isotope = str(element.lower().strip('*')) + str(atomic_mass) + 'm'
            else:
                isotope = str(element.lower().strip('*')) + str(atomic_mass)

            # half life eval
            half_guess = df['T1/2'][line_index]
            if type(half_guess) is str:
                half_present = True
                halflife = half_guess.partition('(')[0]
                half_err = half_guess.partition('(')[2].partition(')')[0]
                divisor = len(halflife.partition('.')[-1])
                try:
                    halflife = float(halflife)
                    half_err = float(half_err) / 10**(divisor + len(half_err) - 1)
                except ValueError:
                    half_present = False
                    halflife = None
                    units = 's'
                units = half_guess.partition(')')[2].replace(' ', '')
                if units == 'ms' and half_present:
                    halflife = halflife / 1000
                    half_err = half_err / 1000
                    
                elif units == 's' and half_present:
                    pass
                elif units == 'm ms' and half_present:
                    halflife = halflife / 1000
                    half_err = half_err / 1000
                elif not half_present:
                    pass
                else:
                    print(f'Unknown units [{units}] for {isotope}')
                    print(f'half_guess is {half_guess}')
                    raise Exception

            else:
                half_present = False


            # Average emissions per decay
            num_emit = 0
            emit_err = 0
            one_val = 0
            two_val = 0
            three_val = 0

            one_emit = df['%P(1n)'][line_index]
            # Check "?" values, 34Na, "~" values not considered
            if type(one_emit) is str:
                one_val = one_emit.partition('(')[0].replace(' ', '').replace('<', '').replace('>', '').replace('~', '').replace('≤', '')
                one_err = one_emit.partition('(')[2].partition(')')[0]
                divisor = len(one_emit.partition('.')[-1].partition('(')[0].replace(' ', ''))
                if one_err == '':
                    one_err = 0
                elif '-' in one_err:
                    one_err = one_err.replace('+', '').split('-')
                    one_err_choice = max(one_err)
                    one_err = float(one_err_choice) / 10**(divisor + len(one_err_choice) - 1)
                else:
                    one_err = float(one_err) / 10**(divisor + len(one_err) - 1)

                try:
                    one_val = float(one_val)
                except ValueError:
                    one_val = 0

            two_emit = df['%P(2n)'][line_index]
            if type(two_emit) is str:
                two_val = two_emit.partition('(')[0].replace(' ', '').replace('<', '').replace('>', '').replace('~', '').replace('≤', '')
                two_err = two_emit.partition('(')[2].partition(')')[0]
                divisor = len(two_emit.partition('.')[-1].partition('(')[0].replace(' ', ''))
                if two_err == '':
                    two_err = 0
                elif '-' in two_err:
                    two_err = two_err.replace('+', '').split('-')
                    two_err_choice = max(two_err)
                    two_err = float(two_err_choice) / 10**(divisor + len(two_err_choice) - 1)
                else:
                    two_err = float(two_err) / 10**(divisor + len(two_err) - 1)
                
                try:
                    two_val = float(two_val)
                except ValueError:
                    two_val = 0

            three_emit = df['%P(3n)'][line_index]
            if type(three_emit) is str:
                three_val = three_emit.partition('(')[0].replace(' ', '').replace('<', '').replace('>', '').replace('~', '').replace('≤', '')
                three_err = three_emit.partition('(')[2].partition(')')[0]
                divisor = len(three_emit.partition('.')[-1].partition('(')[0].replace(' ', ''))
                if three_err == '':
                    three_err = 0
                elif '-' in three_err:
                    three_err = three_err.replace('+', '').split('-')
                    three_err_choice = max(three_err)
                    three_err = float(three_err_choice) / 10**(divisor + len(three_err_choice) - 1)
                else:
                    three_err = float(three_err) / 10**(divisor + len(three_err) - 1)

                try:
                    three_val = float(three_val)
                except ValueError:
                    three_val = 0


            num_emit = (one_val + 2*two_val + 3*three_val) / 100
            emit_err = np.sqrt(one_err**2 + 4*two_err**2 + 9*three_err**2) / 100
            
            #decay_const = np.log(2) / halflife
            pm = u'\u00b1'
            #print(f'{isotope} : {halflife} {pm} {half_err} : {num_emit} {pm} {emit_err}')
            if half_present and num_emit > 1E-13:
                ensdf_data[isotope] = {'halflife' : [halflife, half_err],
                                       'emission' : [num_emit, emit_err]}
                #[halflife, num_emit, half_err, emit_err]
        return ensdf_data


    def simulate_keepin_group_abun(self, groupdata, times, fissions, efficiency):
        """
        Generate the ENSDF dataset, evaluate which group each precursor is in
            and then assign it an abundance based on that group

        Using group abundance values means that the # delayed neutrons per decay
            is not considered. A better approach would be to determine the
            actual concentration of each isotope.

        Parameters
        ----------
        groupdata : dict
            key : string
                Name of group (g1, g2, ...)
            value : list
                Half life, error, rel_abundance, err, yield, err
        times : list
            List of times to simulate
        fissions : float
            Number of fissions in sample
        efficiency : float
            Efficiency of detector (or normalization factor)

        Returns
        -------
        detector_data : list
            Delayed neutrons at each time

        fulldata : dict
            key : str
                Name of isotope (i.e. xe135)
            value : list
                List of half life (s), average emission per decay, and
                abundance
        """
        ensdf = self.parse_file()
        fulldata = dict()
        for isotope in ensdf:
            data_list = ensdf[isotope]
            #input(data_list)
            data_halflife = data_list['halflife'][0]
            prev_min_dist = 1E100
            if data_halflife is None:
                continue
            for group in groupdata:
                # Assign to group based on smallest halflife distance
                cur_dist = abs(groupdata[group][0]**2 - data_halflife**2)
                min_dist = min(cur_dist, prev_min_dist)
                if min_dist == cur_dist:
                    min_group = group
                prev_min_dist = min_dist
                
            abundance = groupdata[min_group][4]/100 # Group i delnu/fiss
            fulldata[isotope] = [data_halflife, data_list['emission'][0], abundance]
        detector_data = list()
        lead_term = fissions * efficiency
        for t in times:
            counts = 0
            for iso in fulldata:
                lamb = np.log(2) / fulldata[iso][0]
                ai = fulldata[iso][2]
                if settings.irradiation == 'pulse':
                    a_val = ai * lamb
                elif settings.irradiation == 'infinite':
                    a_val = ai
                counts += lead_term * a_val * np.exp(-lamb * t)
            detector_data.append(counts)
        return detector_data, fulldata


def spectra_analysis(textfilenames,
                     path='spectra/spectra/',
                     display=False):
    """
    Reads in the IAEA spectra data from given file sets,
        normalizes the data such that the sum is one. This allows
        the y-axis to represent the probability a given neutron will
        reside in the energy bin.

    Parameters
    ----------
    textfilenames : list
        List of text files which go together
    path : str
        Path to file
    display : bool (optional)
        Generate a plot of the spectrum

    Returns
    -------
    bins : numpy 1D array
        center bin values
    values : numpy 1D array
        normalized spectrum values
    """
    bins = list()
    values = list()
    for datafile in textfilenames:
        temp_bins = list()
        temp_values = list()
        with open(path+datafile) as f:
            lines = f.readlines()
            for line in lines:
                if line[0] != '#':
                    data = line.split()
                    bins.append(float(data[0]))
                    values.append(float(data[1]))
                    temp_bins.append(float(data[0]))
                    temp_values.append(float(data[1]))
        temp_bins = np.array(temp_bins)
        temp_values = np.array(temp_values)  
        temp_inds = temp_bins.argsort()
        sort_temp_bins = temp_bins[temp_inds]
        sort_temp_values = temp_values[temp_inds]
        
        if display:
            plt.step(sort_temp_bins, sort_temp_values,
                     label=f'{datafile}', where='mid')
            plt.xlabel('Energy [MeV]')
            plt.ylabel('Provided Value')
            plt.legend()
            plt.tight_layout()
            plt.show()
            plt.close()
    bins = np.array(bins)
    values = np.array(values)             
    inds = bins.argsort()
    bins = bins[inds]
    values = values[inds]
    # Now normalize into probability form
    norm = 1 / sum(values)
    values = values * norm
    

    if display:
        plt.step(bins, values, label='Total', where='mid')
        plt.xlabel('Energy [MeV]')
        plt.ylabel('Probability')
        plt.legend()
        plt.tight_layout()
        plt.show()
        plt.close()
        
    return bins, values




if __name__ == '__main__':
    import keepin_handler
##    fissions = 1#1E16
##    efficiency = 1
##    name = '6keepin235fast'
##    keepin_response = keepin_handler.KEEPIN(name)
##    times = np.arange(0, 10, 0.01)
##    keepin_group_data, keepin_net_data = keepin_response.data_store()
##    print('Uncertainties not incorporated')
##    print('Isomers not incorporated')
##    ensdf_keepin_sim = ENSDF('./ensdf_data/eval_net.xlsx',
##                                             'Sheet1')
##    ensdf_data = ensdf_keepin_sim.parse_file()
##    #for each in ensdf_data:
##        #print(f'{each} : {ensdf_data[each]}')
##    print(f'Total isotopes: {len(ensdf_data)}')
##    ensdf_inserted_isos = ['na34', 'na35', 'si35', 'v61', 'v63', 'co71', 'co72', 'ag125', 'cd133']
##    # na35 Pn = ?, same for ag125
##    # si35 < 5.3 Pn? same for co71, co72
##    # v61 ~= 6.0 Pn? same for v63, cd133, na34
##    ensdf_inserted_halflives = [0.0055, 0.0015, 0.780, 0.047, 0.017, 0.080, 0.062, 0.166, 0.064]
##    ensdf_inserted_emiss_dec = [0.1500, 0.0000, 0.053, 0.060, 0.350, 0.036, 0.080, 0.000, 1.000]
##    
##    ensdf_keepin_delnu, ensdf_keepin_data = ensdf_keepin_sim.simulate_keepin_group_abun(keepin_group_data,
##                                                                                        times,
##                                                                                        fissions,
##                                                                                        efficiency)
##    plt.plot(times, ensdf_keepin_delnu, label='ensdf-keepin')
##    keepin_delnu, keepin_err = keepin_response.simulate_instant(times, fissions, efficiency)
##    plt.plot(times, keepin_delnu, label='keepin')
##    plt.yscale('log')
##    plt.ylabel('Delayed Neutron Count Rate [#/s]')
##    plt.xlabel('Time [s]')
##    plt.legend()
##    plt.show()
    text_files = ['fig27_full.dat']
    bins, vals = spectra_analysis(text_files, display=True)
