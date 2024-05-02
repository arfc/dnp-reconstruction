import numpy as np
import matplotlib.pyplot as plt
import time
from pandas import read_excel
import settings

class SCALE:
    """
    This class handles the SCALE output, including extracting data,
        sorting into a form similar to ENSDF, exporting concentrations,
        plotting with given decay constants for a given isotope,
        determining decay constants given two data points,
        and extracting delayed neutron count data from ORIGEN
    """

    def __init__(self,
                 filename,
                 fissions,
                 efficiency,
                 normalize_value,
                 volume=0.1583105694,
                 mass_normalize=21.90177):
        """
        Initialize

        Parameters
        ----------
        filename : str
            Name of file containing SCALE data
            (This should be a .txt file which has been pre-organized
            to only contain
        ftype : str ['triton', 'origen']
            Type of file
        fissions : float
            Number of fissions
        efficiency : float
            Efficiency of delayed neutron detector
        normalize_value : float
            TRITON mass multiplier
        volume : float
            Volume of the sample

        Returns
        -------
        None
        """
        self.fname = filename
        self.vol = volume
        self.fiss = fissions
        self.eff = efficiency
        self.norm = normalize_value
        self.mnorm = mass_normalize
        self.ENSDF_data = dict()
        self.ENSDF_data['na34'] = {'halflife': [0.0055, 0.0001], 'emission': [1.15, 0.2]}
        self.ENSDF_data['na35'] = {'halflife': [0.0015, 0.00005], 'emission': [0.57, 0.57]}
        self.ENSDF_data['si35'] = {'halflife': [0.78, 0.012], 'emission': [0.053, 0]}
        self.ENSDF_data['v61'] = {'halflife': [0.047, 0.0012], 'emission': [0.06, 0]}
        self.ENSDF_data['v63'] = {'halflife': [0.017, 0.0003], 'emission': [0.35, 0]}
        self.ENSDF_data['co71'] = {'halflife': [0.008, 0.0003], 'emission': [0.036, 0.009]}
        self.ENSDF_data['co72'] = {'halflife': [0.0062, 0.0003], 'emission': [0.08, 0.02]}
        #self.ENSDF_data['ag125'] = {'halflife': [0.0015, 0.00005], 'emission': [0, 0]}
        self.ENSDF_data['cd133'] = {'halflife': [0.0064, 0.0008], 'emission': [1, 0.5]} # Fabricated uncertainty

        return

    def origen_activity_parser_deprecated(self, timestep):
        """
        Parses the ORIGEN .out style of file. At a particular timestep,
            gathers activity data. Combined with concentration data,
            can generate decay constant values.

        Parameters
        ----------
        timestep : int
            Time index at which to evaluate

        Returns
        -------
        act_data : dict
            key : str
                Name of isotope (i.e. xe135) (hyphens are removed)
            value : float
                Value of activity immediately after pulse 
        """
        check_val = 'Nuclide concentrations in becquerels'
        end_val = 'totals'
        #decay_check_val = 'Nuclide concentrations in becquerels'
        act_data = dict()
        marker_indeces = list()
        end_indeces = list()
        recent_start = False
        net_composition = 0
        with open(self.fname, 'r') as f:
            lines = f.readlines()
        # Need to find each instance instead of only one
        for ind, line in enumerate(lines):
            if check_val in line:
                marker_indeces.append(ind)
                recent_start = True
            elif recent_start and end_val in line:
                recent_start = False
                end_indeces.append(ind)
            else:
                pass
        marker_indeces = marker_indeces[:-1]
        end_indeces = end_indeces[:-1]
        start_offset = 6
        end_offset = -1
        for type_index in range(len(marker_indeces)):
            for line in lines[marker_indeces[type_index]+start_offset:end_indeces[type_index] + end_offset]:
                data = line.split()
                element = data[0].split('-')[0]
                weight = data[0].split('-')[1]
                isotope = str(element) + str(weight)
                activ = float(data[1+timestep]) / self.mnorm
                net_composition += activ
                if isotope in act_data.keys():
                    act_data[isotope] += activ
                else:
                    act_data[isotope] = activ
        return act_data


    def origen_activity_parser(self, timestep):
        """
        Parses the ORIGEN .out style of file. At a particular timestep,
            gathers activity data. Combined with concentration data,
            can generate decay constant values.

        Parameters
        ----------
        timestep : int
            Time index at which to evaluate

        Returns
        -------
        act_data : dict
            key : str
                Name of isotope (i.e. xe135) (hyphens are removed)
            value : float
                Value of activity immediately after pulse 
        """
        # So far have branch_frac * lambda * atoms
        #   print isotopes that dont have branch_frac data
        check_val = 'nuclide     atoms         (1/s)     fraction     n/s         MeV'
        end_val = 'total'
        act_data = dict()
        marker_indeces = list()
        end_indeces = list()
        recent_start = False
        with open(self.fname, 'r') as f:
            lines = f.readlines()
        # Need to find each instance instead of only one
        for ind, line in enumerate(lines):
            if check_val in line:
                marker_indeces.append(ind)
                recent_start = True
            elif recent_start and end_val in line:
                recent_start = False
                end_indeces.append(ind)
            else:
                pass
        start_offset = 1
        end_offset = -1
        for type_index in range(len(marker_indeces)):
            for line in lines[marker_indeces[type_index]+start_offset:end_indeces[type_index] + end_offset]:
                data = line.split()
                element = data[0].split('-')[0]
                weight = data[0].split('-')[1]
                isotope = str(element) + str(weight)
                if isotope in act_data.keys():
                    if isotope[-1] != 'm':
                        isotope += 'm'
                        act_data[isotope] = float(data[2])
                    else:
                        print(f'{isotope} already exists')
                        raise Exception
                else:
                    act_data[isotope] = float(data[2])
        return act_data

    def origen_delnu_parser(self,
                            collect):
        """
        Parses the ORIGEN .out style of file. Collects the delayed neutrons
            from each isotope as well as the time values for each

        Parameters
        ----------
        collect : str
            Either a specific isotope ('xe135'), or 'all' for the total

        Returns
        -------
        times : vector
            Times (seconds) used in ORIGEN
        counts : vector
            Total delayed neutrons emitted in ORIGEN
        """
        check_val = 'Delayed neutron intensity by nuclide (neutrons/sec) for case'
        end_val = 'total'
        marker_indeces = list()
        time_indeces = list()
        end_indeces = list()
        recent_start = False
        net_composition = 0
        self.pure_dict = dict()
        with open(self.fname, 'r') as f:
            lines = f.readlines()
        # Need to find each instance instead of only one
        for ind, line in enumerate(lines):
            if check_val in line:
                marker_indeces.append(ind)
                recent_start = True
            elif recent_start and end_val in line:
                recent_start = False
                end_indeces.append(ind)
            elif 'time' in line:
                time_indeces.append(ind)
            else:
                pass
        marker_index = marker_indeces[0]
        end_index = end_indeces[0]
        times_index = time_indeces[-3]
        start_offset = 3
        end_offset = 1
        time_line = True
        times = list()
        counts = list()
        t_begin = 2
        for line in lines[times_index+t_begin:]:
            data = line.split()
            try:
                times.append(float(data[1]))
            except ValueError:
                break
        
        for line in lines[marker_index+start_offset:end_index + end_offset]:
            temp_counts = list()
            data = line.split()
            if len(data) == 1:
                continue
            elif data[0] == 'total' and collect == 'all':
                for cnt in data[1:]:
                    counts.append(float(cnt) * self.eff / self.mnorm)
            elif data[0] == 'total':
                pass
            else:
                try:
                    element = data[0].split('-')[0]
                    weight = data[0].split('-')[1]
                except IndexError:
                    continue
                isotope = str(element) + str(weight)
                for cnt in data[1:]:
                    try:
                        temp_counts.append(float(cnt) * self.eff / self.mnorm)
                    except ValueError:
                        temp_counts.append(0)
                self.pure_dict[isotope] = temp_counts
                if isotope == collect:
                    counts = temp_counts.copy()
        return times, counts

    def origen_spectra_parser(self):
        """
        Parses the ORIGEN .out style of file for spectra.

        Parameters
        ----------
        None

        Returns
        -------
        time_data : 1D numpy array
            Time values used
        energy_data : 1D numpy array
            Energy bin midpoints
        spectra_matrix : 2D numpy array
            Rows are the energies, columns are times, values are counts
        bin_data : 1D numpy array
            Bins used
        """
        check_val = 'Delayed neutron source intensity (1/s) as a function of time'
        end_val = 'total'
        recent_start = False
        energy_data = list()
        spectra_matrix = list()
        bin_data = list()
        with open(self.fname, 'r') as f:
            lines = f.readlines()
        # Need to find relative end
        for ind, line in enumerate(lines):
            if check_val in line:
                marker_indeces = ind
                recent_start = True
            elif recent_start and end_val in line:
                end_indeces = ind
                recent_start = False
            else:
                pass
        start_offset = 3
        end_offset = -1
        for line in lines[marker_indeces+start_offset : end_indeces+end_offset]:
            cur_row = list()
            data = line.split()
            if data[0] == 'boundaries':
                time_data = [float(x.replace('sec', '')) for x in data[2:]]
            if data[1] == '-':
                curbin_data = [float(data[0]), float(data[2])]
                cur_mid_energy = (float(data[0]) + float(data[2])) / 2
                energy_data.append(cur_mid_energy)
                try:
                    cur_row = [float(x.replace('sec', '')) for x in data[3:]]
                except ValueError:
                    # ORIGEN sometimes removes the E in sci notation
                    list_row = list()
                    for val_inspect in data[3:]:
                        if 'E' in val_inspect:
                            list_row.append(float(val_inspect))
                        else:
                            if '+' in val_inspect:
                                usesplitter = '+'
                            elif '-' in val_inspect:
                                usesplitter = '-'
                            val_split = val_inspect.split(usesplitter)
                            list_row.append(float(val_split[0] + 'E' +
                                                  usesplitter + val_split[1]))
                    cur_row.append(list_row)
                    cur_row = cur_row[0]
                bin_data.append(curbin_data)
                spectra_matrix.append(cur_row)
        energy_data = np.array(energy_data)
        time_data = np.array(time_data)
        spectra_matrix = np.array(spectra_matrix)
        bin_data = np.array(bin_data)
        # Energy bins large to small; reverse
        spectra_matrix = np.flip(spectra_matrix, 0)
        energy_data = np.flip(energy_data, 0)
        bin_data = np.flip(bin_data, 0)
        return time_data, energy_data, spectra_matrix, bin_data

    def origen_parser(self):
        """
        Parses the ORIGEN .out style of file. Volume is normalized to 1 and
            atoms need to be divided by normalization factor

        Parameters
        ----------
        None

        Returns
        -------
        comp_data : dict
            key : str
                Name of isotope (i.e. xe135) (hyphens are removed)
            value : numpy array
                Values of atom/barn-cm immediately after pulse 
        """
        check_val = 'Nuclide concentrations in atoms/barn-cm'
        end_val = 'totals'
        comp_data = dict()
        marker_indeces = list()
        end_indeces = list()
        recent_start = False
        net_composition = 0
        with open(self.fname, 'r') as f:
            lines = f.readlines()
        # Need to find each instance instead of only one
        for ind, line in enumerate(lines):
            if check_val in line:
                marker_indeces.append(ind)
                recent_start = True
            elif recent_start and end_val in line:
                recent_start = False
                end_indeces.append(ind)
            else:
                pass
        marker_indeces = marker_indeces[:-1]
        end_indeces = end_indeces[:-1]
        start_offset = 6
        end_offset = -1
        for type_index in range(len(marker_indeces)):
            for line in lines[marker_indeces[type_index]+start_offset:end_indeces[type_index] + end_offset]:
                data = line.split()
                element = data[0].split('-')[0]
                weight = data[0].split('-')[1]
                isotope = str(element) + str(weight)
                iso_data = list()
                for data_point in data[1:]:
                    try:
                        float(data_point)
                    except ValueError:
                        if '-' in data_point:
                            choice = '-'
                        elif '+' in data_point:
                            choice = '+'
                        new_data = data_point.split(choice)
                        new_point = new_data[0] + 'E' + choice + new_data[1]
                        data_point = float(new_point)
                    iso_data.append(float(data_point) / self.norm)
                #conc = float(data[1+timestep])
                #net_composition += conc / self.norm
                iso_data = np.asarray(iso_data)
                if isotope in comp_data.keys():
                    comp_data[isotope] += iso_data#conc / self.norm
                else:
                    comp_data[isotope] = iso_data#conc / self.norm
        #print(f'Net atoms: {net_composition * 1E24}')
        return comp_data
        

    def triton_parser(self):
        """
        Parses the TRITON .out style of file

        Parameters
        ----------
        None

        Returns
        -------
        comp_data : dict
            key : str
                Name of isotope (i.e. xe135) (hyphens are removed)
            value : float
                Value of atom/barn-cm immediately after pulse       
        """
        check_val = 'end-of-step 1 isotopics'
        comp_data = dict()
        with open(self.fname, 'r') as f:
            lines = f.readlines()
        for ind, line in enumerate(lines):
            if check_val in line:
                marker_index = ind
                break
            else:
                pass
        data_start_index = marker_index + 13
        num_isotopes = 2237
        net_composition = 0
        # Values are lost to duplicate isotopes
        for line in lines[data_start_index:data_start_index+num_isotopes]:
            data = line.split()
            base = data[0].split(':')[0]
            element = base.split('-')[1]
            weight = base.split('-')[2]
            isotope = str(element) + str(weight)
            conc = float(data[-1])
            net_composition += conc
            #####
            #ensdf_inserted_isos = ['na34', 'na35', 'si35', 'v61', 'v63', 'co71', 'co72', 'ag125', 'cd133']
            #if isotope in ensdf_inserted_isos:
            #    input(f'{isotope} : {conc}')
            #####
            if isotope in comp_data.keys():
                comp_data[isotope] += conc
            else:
                comp_data[isotope] = conc
        #print(f'Net atoms: {net_composition * self.vol * 1E24}')
        return comp_data

    def gen_comp_data(self, timestep):
        """
        Parse file for each isotope and concentration for each time step.
            Because only TRITON files mention TRITON, can search for that phrase
            to determine file type.

        Paramters
        ---------
        None

        Returns
        -------
        comp_data : dict
            key : str
                Name of isotope (i.e. xe135) (hyphens are removed)
            value : numpy array
                Value of atom/barn-cm at each time step (initial for TRITON, deprecated)
        """
        check_val = 'TRITON'
        with open(self.fname, 'r') as f:
            lines = f.readlines()
        ftype = 'origen'
        for ind, line in enumerate(lines):
            if check_val in line:
                ftype = 'triton'
                break
            else:
                pass
        print('-'*40)
        if ftype == 'triton':
            print('\nTRITON Concentrations')
            comp_data = self.triton_parser()
            raise Exception('TRITON concentrations not time dependent')
        elif ftype == 'origen':
            print('\nORIGEN Concentrations')
            comp_data = self.origen_parser()
        else:
            raise Exception
        self.ftype = ftype

        return comp_data

    def ensdf_matcher(self, ensdf_dict, timestep, target='all'):
        """
        Create a formated set of atom counts from directly after
            the pulse occurs. These concentrations are then paired with
            their associated Pn and lambda values from ENSDF.

        Parameters
        ----------
        ensdf_data : dict
            key : str
                Name of isotope (i.e. xe135)
            value : dict
                key : str
                    Identifier (emissions, halflife, conc)
                value : list
                    value, uncertainty

        Returns
        -------
        net_data : dict
            key : str
                Name of isotope (i.e. xe135)
            value : list
                key : str
                    Identifier (emissions, halflife, conc)
                value : list
                    value, uncertainty
        """
        timestep = 0
        scale_data = self.gen_comp_data(timestep)
        ensdf_data_copy = ensdf_dict.copy()

        # Add ENSDF data to IAEA data
        ensdf_data_copy['na34'] = self.ENSDF_data['na34']
        ensdf_data_copy['na35'] = self.ENSDF_data['na35']
        ensdf_data_copy['si35'] = self.ENSDF_data['si35']
        ensdf_data_copy['v61'] = self.ENSDF_data['v61']
        ensdf_data_copy['v63'] = self.ENSDF_data['v63']
        ensdf_data_copy['co71'] = self.ENSDF_data['co71']
        ensdf_data_copy['co72'] = self.ENSDF_data['co72']
        #ensdf_data_copy['ag125'] = self.ENSDF_data['ag125']
        ensdf_data_copy['cd133'] = self.ENSDF_data['cd133']
        #####
        
        net_data = dict()
        in_ensdf_not_origen_count = 0
        print(f'    Number {self.ftype.upper()} isos: {len(scale_data)}')
        print(f'    Number IAEA isos: {len(ensdf_dict)}')
        print('Using ORIGEN concentration uncertainties, metastable same as base')
        for isotope in ensdf_data_copy:
            if target == 'all':
                pass
            elif target == isotope:
                pass
            else:
                continue
            try:
                atoms_barn_cm = scale_data[isotope]
            except KeyError:
                #print(f'Isotope {isotope} not found in origen data')
                in_ensdf_not_origen_count += 1
                continue
            atoms = atoms_barn_cm * 1E24 * self.vol
            atom_err = self.conc_uncert(isotope)
            ensdf_data_copy[isotope]['conc'] = [atoms, atom_err]
            net_data[isotope] = ensdf_data_copy[isotope]
        print(f'    {in_ensdf_not_origen_count} isotopes in IAEA but not {self.ftype}')
        #print(f'3g == {scale_data["u235"] * 1E24 * self.vol / 6.022E23 * 235}g')
        return net_data

    def conc_uncert(self,
                    isotope,
                    filename='./scale_outputs/response_table.1.stddev.xlsx',
                    sheetname='response_table.1.stddev'):
        """
        Extract uncertainty in ORIGEN concentrations from response table
            csv file.

        Parameters
        ----------
        isotope : str
            Name of isotope to get uncertainty data for
        
        Returns
        -------
        uncertainty : float
            Uncertainty in concentration
        """
        df = read_excel(filename,
                        sheet_name=sheetname)
        if isotope[-1] == 'm':
            isotope = isotope[:-1]
        search_name = f'irrad:oriout.{isotope}'
        if isotope == 'sb134':
            uncertainty = 0
        else:
            try:
                uncertainty = df[search_name][0]
            except KeyError:
                uncertainty = 0
                #print(f'        {isotope} not found')
        return uncertainty

    def simulate_ensdf_SCALE(self, times, ensdf_dict, timestep, detect_isotope='all',
                             activity='ENSDF', errs=True):
        """
        Simulate the delayed neutron response based on ENSDF
            using SCALE composition data.

        Parameters
        ----------
        times : list
            Times at which to evaluate delayed neutron emissions
        ensdf_data : dict
            key : str
                Name of isotope (i.e. xe135)
            value : list
                List of half life (s) followed by average emission per decay
        detect_isotope : str
            Name of isotope to detect delayed neutrons from. Can be set to 'all'
        activity : str
            Where to pull decay constant data from (either 'ENSDF', 'ORIGEN', 'LAMDEBUG', 'PNDEBUG')
        errs : bool
            Whether or not to calculate the uncertainty for each time

        Returns
        ------
        counts : list
            List of counts evaluated at each time provided
        """
        net_data = self.ensdf_matcher(ensdf_dict, timestep, target=detect_isotope)
        counts = list()
        iso_list = list()
        lam_list = list()
        Pn_list = list()
        atoms_list = list()
        lam_err_list = list()
        errors = list()

        #input('Temporary debug measure (Press Continue)')
        for each in settings.DEBUG_IGNORE_ISOTOPES:
            print(f'DEBUG IGNORE {each}')
            net_data.pop(each)
        #times, _ = self.origen_delnu_parser('all')

        prev_max_iso = ''
        prev_debug_worst_iso = ''
        print(f'    Using {len(net_data)} isotopes')
        print(f'{activity.upper()} decay constants')
        if activity.upper() == 'PNDEBUG' or \
           activity.upper() == 'PUREDEBUG' or \
           activity.upper() == 'PURECHECK':
            # Generate self.pure_dict (counts)
            pure_time_data, discard = self.origen_delnu_parser(detect_isotope)
        if activity.upper() == 'ORIGEN' or \
           activity.upper() == 'LAMDEBUG' or \
           activity.upper() == 'PNDEBUG' or \
           activity.upper() == 'PUREDEBUG':
            activity_data = self.origen_activity_parser(timestep)
            # Trim net data to only contain isos also in activity data
            trimmed = dict()
            print(f'Removed isotopes due to lacking emission/decay data')
            for conc_iso in net_data.keys():
                if conc_iso in activity_data.keys():
                    trimmed[conc_iso] = net_data[conc_iso]
            net_data = trimmed.copy()
        elif activity.upper() == 'ENSDF' or \
             activity.upper() == 'PURECHECK':
            pass
        else:
            print(f'Activity {activity} not recognized')
            raise Exception
        for isotope in net_data:
            if net_data[isotope]['halflife']:
                # DEFAULT DATA IS IAEA WITH ORIGEN CONCENTRATIONS
                Pn = net_data[isotope]['emission'][0]
                atoms = net_data[isotope]['conc'][0]
                lam = np.log(2) / net_data[isotope]['halflife'][0]
    
                if activity.upper() == 'ENSDF':
                    if errs:
                        lam_err_list.append(np.log(2) /
                                            net_data[isotope]['halflife'][0]**2 *
                                            net_data[isotope]['halflife'][1])
                elif activity.upper() == 'ORIGEN':
                    lam = activity_data[isotope]
#                elif activity.upper() == 'ORIGEN' or \
#                     activity.upper() == 'PNDEBUG':
#                    puori_times, puori_counts = self.origen_delnu_parser(isotope)
#                    lam = np.log(puori_counts[-1] / puori_counts[0]) / (puori_times[0] - puori_times[-1]) #activity_data[isotope]
#                    if errs:
#                        raise Exception
#
#                    if activity.upper() == 'PNDEBUG':
#                        ENSDF_Pn = Pn
#                        ORIGEN_Pn = self.pure_dict[isotope][0] / (lam * atoms * self.eff) # at t=0, so no exp
#                        Pn = ORIGEN_Pn
#                        if round(abs(ENSDF_Pn - ORIGEN_Pn), 2) != 0.0:
#                            #print(f'ENSDF {isotope}: {lam}')
#                            #print(f'ORIGEN {isotope}: {lam2}')
#                            print(f'{isotope} % diff: {abs(ORIGEN_Pn - ENSDF_Pn) / ENSDF_Pn * 100}')
#                elif activity.upper() == 'LAMDEBUG':
#                    # Decay constants
#                    lam2 = np.log(2) / net_data[isotope]['halflife'][0] # IAEA lam
#                    puori_times, puori_counts = self.origen_delnu_parser(isotope)
#                    lam = np.log(puori_counts[-1] / puori_counts[0]) / (puori_times[0] - puori_times[-1]) #activity_data[isotope] # Pure lam
#                    if errs:
#                        raise Exception
#                    if round(abs(lam - lam2), 2) != 0.0:
#                        #print(f'ENSDF {isotope}: {lam}')
#                        #print(f'ORIGEN {isotope}: {lam2}')
#                        print(f'{isotope} % diff: {abs(lam - lam2) / lam2 * 100}')
#                        #input()
#                        pass
#                    # Pn values
#                    #ENSDF_Pn = net_data[isotope]['emission'][0]
#                    #ORIGEN_Pn = self.pure_dict[iso] / (lam * atoms) # at t=0, so no exp
#                else:
#                    print(f'Activity {activity}')
#                    raise Exception('Unknown activity')
                iso_list.append(isotope)
                lam_list.append(lam)
                Pn_list.append(Pn)
                atoms_list.append(atoms)
        print(f'    Times : Most impactful isotope during that time : halflife')
        saved_max_iso = ''
        for cur_t_ind, t in enumerate(times):
            detect = 0
            max_count = -1
            max_ind = 0
            max_iso = ''
            max_half = list()
            cur_err = 0
            worst_atoms = 0
            worst_counts = 0
            worst_count_val = 0
            debug_max_diff = 0
            debug_worst_iso = ''
            emiss_max_diff = 0
            emiss_worst_iso = ''
            if activity.upper() == 'PNDEBUG' or activity.upper() == 'PUREDEBUG':
                pn_checked = False
            else:
                pn_checked = True
            for ind, isotope in enumerate(iso_list):
                if isotope == detect_isotope or detect_isotope == 'all':
                    lam = lam_list[ind]
                    Pn = Pn_list[ind]
                    atoms = atoms_list[ind]
                    count_val = self.eff * Pn * lam * atoms[cur_t_ind] # * np.exp(-lam * t) #Using prev defined vals
                    if activity.upper() == 'LAMDEBUG':
                        # Calculate PURE ORIGEN lambda
                        #puori_times, puori_counts = self.origen_delnu_parser(isotope)
                        #if np.isclose(puori_counts[0], puori_counts[1]):
                        #        target_index = -1
                        #else:
                        #        target_index = 1
                        #lam2_top = np.log(puori_counts[target_index] / puori_counts[0])
                        lam2 = activity_data[isotope]#lam2_top / (puori_times[0] - puori_times[target_index]) #activity_data[isotope]
                        other_count_val = self.eff * Pn * lam2 * atoms[cur_t_ind]# * np.exp(-lam2 * t) #PURE lambda counts
                        rel_diff = abs(count_val - other_count_val)
                        if rel_diff > debug_max_diff:
                            worst_atoms = atoms[cur_t_ind]
                            debug_worst_iso = isotope
                            debug_max_diff = rel_diff
                            worst_iaea_lam = lam
                            worst_ori_lam = lam2
                            worst_counts = count_val # IAEA
                            worst_count_val = other_count_val # PURE
                    elif activity.upper() == 'PNDEBUG' or activity.upper() == 'PUREDEBUG':
                        if round(t, 5) in np.round(pure_time_data, 5):
                            pn_checked = True
                            #puori_times, puori_counts = self.origen_delnu_parser(isotope)
                            #if np.isclose(puori_counts[0], puori_counts[1]):
                            #    target_index = -1
                            #else:
                            #    target_index = 1
                            lam2 = activity_data[isotope] #np.log(puori_counts[target_index] / puori_counts[0]) / (puori_times[0] - puori_times[target_index]) #activity_data[isotope]
                            Pn2 = self.pure_dict[isotope][0] / (lam2 * atoms[0] * self.eff)
                            if activity.upper() == 'PUREDEBUG':
                                lam2 = lam2 # Change Pn and lambda
                            else:
                                lam2 = lam # Keeping IAEA lam
                            time_index = np.where(np.isclose(t, pure_time_data))[0][0]
                            other_count_val = self.eff * Pn2 * lam2 * atoms[cur_t_ind] #* np.exp(-lam2 * t)#self.pure_dict[isotope][time_index] # IAEA
                            rel_diff = abs(count_val - other_count_val)
                            if rel_diff > debug_max_diff:
                                worst_counts = count_val # IAEA
                                worst_count_val = other_count_val # PURE
                                worst_atoms = atoms[cur_t_ind]
                                debug_worst_iso = isotope
                                debug_max_diff = rel_diff
                                worst_iaea_lam = np.log(2) / net_data[isotope]['halflife'][0]
                                worst_ori_lam = lam2 #activity_data[isotope]
                                worst_iaea_pn = Pn
                                ori_lam = lam2 #activity_data[isotope]
                                worst_ori_pn = Pn2
                    elif activity.upper() == 'PURECHECK':
                        if round(t, 5) in np.round(pure_time_data, 5) and isotope in self.pure_dict:
                            time_index = np.where(np.isclose(t, pure_time_data))[0][0]
                            count_val = self.pure_dict[isotope][time_index]
                    detect += count_val
                    if errs:
                        Pn_err = net_data[isotope]['emission'][1]
                        lam_err = lam_err_list[ind]
                        atom_err = net_data[isotope]['conc'][1]
                        cur_err += ((lam * atoms[0] * np.exp(-lam * t) * Pn_err)**2 +
                                    (Pn * atoms[0] * (1-lam*t) * np.exp(-lam*t) * lam_err)**2 +
                                    (Pn * lam * np.exp(-lam*t) * atom_err)**2)
                    if count_val > max_count:
                        max_count = count_val
                        max_iso = isotope
                        max_ind = ind
                        max_half = np.log(2) / lam
            if debug_worst_iso != prev_debug_worst_iso:
                if not pn_checked:
                    pass
                else:
                    print(f'    {np.round(t, 4)}s: Worst iso {debug_worst_iso}')
                    print(f'        Atoms : {worst_atoms}')
                    if activity.upper() == 'LAMDEBUG':
                        print(f'        Lambda : IAEA {worst_iaea_lam} : ORIGEN {worst_ori_lam}')
                    elif activity.upper() == 'PNDEBUG':
                        print(f'        Pn : IAEA {worst_iaea_pn} : ORIGEN {worst_ori_pn}')
                    elif activity.upper() == 'PUREDEBUG':
                        print(f'        Pn : IAEA {worst_iaea_pn} : ORIGEN {worst_ori_pn}')
                        print(f'        Lambda : IAEA {worst_iaea_lam} : ORIGEN {worst_ori_lam}')
                    print(f'        Counts : IAEA {worst_counts}    : ORIGEN   {worst_count_val}')
                    prev_debug_worst_iso = debug_worst_iso
            #print(f'{max_iso} : {prev_max_iso}')
            if max_iso != prev_max_iso:
                if activity.upper() != 'PURECHECK':
                    print(f'    {np.round(t, 4)}s: {max_iso} : {max_half} s')
                    print(f'          Max count: {max_count}')
                    prev_max_iso = max_iso
                elif activity.upper() == 'PURECHECK':
                    #if round(t, 5) in np.round(pure_time_data, 5) and max_iso != saved_max_iso:
                    print(f'    {np.round(t, 4)}s: {max_iso} : {max_half} s')
                    print(f'          Max count: {max_count}')
                    saved_max_iso = max_iso
                    prev_max_iso = max_iso
                #print(f'Lambda : Pn : Atoms')
                #print(f'{lam_list[max_ind]} : {Pn_list[max_ind]} : {atoms_list[max_ind]}')
                
            counts.append(detect)
            errors.append(np.sqrt(cur_err) * self.eff)
            #print(atoms)
        return counts, errors

        

        
if __name__ == '__main__':
    # INITIALIZE
    begin = time.time()
    import ensdf_handler
    import keepin_handler
    import misc_funcs
    from settings import *
    

    # GENERATE DATA

    ensdf_gen = ensdf_handler.ENSDF('./ensdf_data/eval_net.xlsx',
                                             'Sheet1')
    ensdf_dict = ensdf_gen.parse_file()


    activity = 'puredebug'
    errs = False
    filename = './scale_outputs/godiva_irrad_post_pulse.out'
    runname = 'ENSDF-ORIGEN'
    timestep = 0
    ORIGEN_gen = SCALE(filename, fissions,
                      efficiency, normalize_value, volume,
                       mass_normalize)
    #generic_data = ORIGEN_gen.origen_parser(show_iso)
    ORI_counts, ORI_err = ORIGEN_gen.simulate_ensdf_SCALE(times, ensdf_dict,
                                                 timestep, show_iso, activity,
                                                          errs=errs)
    time_data, energy_data, spectra_matrix, bin_data = ORIGEN_gen.origen_spectra_parser()

    # Spectra for different times
##    for tind, t in enumerate(time_data):
##        norm_factor = spectra_normalize / np.sum(spectra_matrix[:, tind])
##        plt.step(energy_data, spectra_matrix[:, tind] * norm_factor)
##        plt.title(f'Spectra at {t} s')
##        plt.xlabel('Energy [MeV]')
##        plt.ylabel(f'Relative Intensity')
##        plt.show()
##        plt.close()

    # Counts over time for given energy
##    for eind, e in enumerate(energy_data):
##        norm_factor = spectra_normalize / np.sum(spectra_matrix[eind, :])
##        plt.step(time_data, spectra_matrix[eind, :] * norm_factor)
##        plt.title(f'Spectra at {e} MeV')
##        plt.xlabel('Time [s]')
##        plt.ylabel(f'Relative Intensity')
##        plt.show()
##        plt.close()


    # Heatmap - 10^4 counts
    fig, ax = plt.subplots()
    x, y = np.meshgrid(time_data, energy_data)
    # Column is energy spectra for given energy
    z = np.zeros(np.shape(spectra_matrix))
    for tind, t in enumerate(time_data):
        norm_factor = spectra_normalize / np.sum(spectra_matrix[:, tind])
        z[:, tind] = norm_factor * spectra_matrix[:, tind]
    c = ax.pcolormesh(x, y, z, cmap='magma')
    cbar = fig.colorbar(c, ax=ax)
    cbar.set_label('Relative Counts')
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [MeV]')
    plt.tight_layout()
    #ax.set_zlabel('Relative Intensity')
    plt.show()
    plt.close()


    # 3D counts, times, energy
##    from matplotlib import cm
##    x, y = np.meshgrid(time_data, energy_data)
##    # Column is energy spectra for given energy
##    z = np.zeros(np.shape(spectra_matrix))
##    for tind, t in enumerate(time_data):
##        norm_factor = spectra_normalize / np.sum(spectra_matrix[:, tind])
##        z[:, tind] = norm_factor * spectra_matrix[:, tind]
##    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
##    surf = ax.plot_surface(x, y, z, cmap=cm.magma,
##                       linewidth=1, antialiased=True)
##    fig.colorbar(surf, shrink=0.5, aspect=5)
##    plt.xlabel('Time [s]')
##    plt.ylabel('Energy [MeV]')
##    ax.set_zlabel('Relative Intensity')
##    plt.show()



    
    #print(times)
    #print(ORI_counts)
    #plt.errorbar(times, ORI_counts, yerr=ORI_err)
    #plt.ylabel('log')
    #plt.show()

##    ORIGEN_gen.origen_delnu_parser('all')
##
##
##    ensdf_dict = ensdf_gen.parse_file()
##    filename = './scale_outputs/godiva_3d_depl.out'
##    runname = 'ENSDF-TRITON'
##    TRITON_gen = SCALE(filename, fissions,
##                      efficiency, normalize_value, volume)
##    #generic_data = scale_gen.triton_parser()
##    TRI_counts = TRITON_gen.simulate_ensdf_SCALE(times, ensdf_dict, timestep, show_iso, activity)
##    plt.plot(times, TRI_counts, label=f'{show_iso} {runname}')
##    print(f'{runname} n/f: {misc_funcs.delnu_per_fiss(times, TRI_counts, fissions, efficiency)}\n')
##
##
##
##    activity = 'origen' #ensdf, origen, debug
##    filename = './scale_outputs/godiva_irrad_post_pulse.out'
##    runname = 'ENSDF-ORIGEN-ACT'
##    timestep = 0
##    ORIGEN_gen = SCALE(filename, fissions,
##                      efficiency, normalize_value, volume)
##    #generic_data = ORIGEN_gen.origen_parser(0)
##    ORI_ACT_counts = ORIGEN_gen.simulate_ensdf_SCALE(times, ensdf_dict, timestep, show_iso, activity)
##    plt.plot(times, ORI_ACT_counts, label=f'{show_iso} {runname}')
##    print(f'{runname} n/f: {misc_funcs.delnu_per_fiss(times, ORI_ACT_counts, fissions, efficiency)}\n')
##
##    
##    name = '6keepin235fast'
##    keepin_response = keepin_handler.KEEPIN(name)
##    plt.plot(keepin_response.true_data_time, keepin_response.true_data_resp,
##             label='Keepin True', linestyle='', marker='.')
##    #print(f'Keepin True n/f: {misc_funcs.delnu_per_fiss(keepin_response.true_data_time, keepin_response.true_data_resp, fissions, efficiency)}\n')
##
##    
##    keepin_delnu = keepin_response.simulate_instant(times, fissions, efficiency)
##    print(f'Keepin Fit n/f: {misc_funcs.delnu_per_fiss(times, keepin_delnu, fissions, efficiency)}\n')
##    plt.plot(times, keepin_delnu, label='Keepin')
##
##
##    name = '6brengland235fast'
##    keepin_response = keepin_handler.KEEPIN(name)
##    brady_england_delnu = keepin_response.simulate_instant(times, fissions, efficiency)
##    print(f'Brady-England Fit n/f: {misc_funcs.delnu_per_fiss(times, brady_england_delnu, fissions, efficiency)}\n')
##    plt.plot(times, keepin_delnu, label='Brady-England')
##    
##    plt.yscale('log')
##    plt.ylabel('Delayed Neutron Count Rate [#/s]')
##    plt.xlabel('Time [s]')
##    plt.legend()
##    plt.savefig('delnu_tot_compare.png')
##    #plt.show()
##    plt.close()


    
    end = time.time()
    print(f'Took {end - begin}s')
