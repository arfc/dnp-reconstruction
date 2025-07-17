import json

class InputHandler:
    def __init__(self, input_path: str) -> None:
        """
        This class handles the input files for all processing.

        Parameters
        ----------
        input_path : str
            Path to the input file.
        """
        self.input_path = input_path
        return None
    
    def read_input(self, check=True, adjust_data=True) -> dict:
        """
        Read the input file and return the data as a dictionary.

        Parameters
        ----------
        check : bool, optional
            If True, checks the behaviour of the input data. Default is True.
        adjust_data : bool, optional
            If True, adjusts the data to fit desired formatting. Default is True.

        Returns
        -------
        output : dict
            Dictionary containing settings and data selections.
        """
        with open(self.input_path, 'r') as file:
            output = json.load(file)
        if check:
            self._check_behaviour(output)
        if adjust_data:
            self._adjust_data(output)
        return output
    
    def _check_behaviour(self, data: dict) -> None:
        """
        Check the behaviour of the input data.

        Parameters
        ----------
        data : dict
            Dictionary containing settings and data selections.
        
        Raises
        ------
        ValueError
            If the behaviour is not supported.
        """
        def _check_if_in_options(item, options):
            if item not in options:
                raise ValueError(f'Option {item} not supported in {options}')
            return None

        if data['file_options']['unprocessed_data_dir'] == data['file_options']['processed_data_dir']:
            raise ValueError("Unprocessed and processed data directories cannot be the same.")
        if data['file_options']['unprocessed_data_dir'] == data['file_options']['output_dir']:
            raise ValueError("Unprocessed data directory cannot be the same as the output directory.")
        if data['file_options']['processed_data_dir'] == data['file_options']['output_dir']:
            raise ValueError("Processed data directory cannot be the same as the output directory.")
        if data['modeling_options']['parent_feeding'] and not data['modeling_options']['concentration_handling'] == 'depletion':
            raise ValueError("Parent feeding option requires depletion method for concentration handling")
        
        possible_half_life_options = ['test-data', 'endfb71', 'iaea']
        _check_if_in_options(data['data_options']['half_life']['data'], possible_half_life_options)
        possible_cross_sections_options = ['test-data', 'endfb71']
        _check_if_in_options(data['data_options']['cross_section']['data'], possible_cross_sections_options)
        possible_fission_yields_options = ['test-data', 'endfb71']
        _check_if_in_options(data['data_options']['fission_yield']['data'], possible_fission_yields_options)
        possible_emission_probabilities = ['test-data', 'iaea']
        _check_if_in_options(data['data_options']['emission_probability']['data'], possible_emission_probabilities)
        possible_fy_names = ['nfy', 'chain']
        _check_if_in_options(data['data_options']['fission_yield']['name'], possible_fy_names)        
        possible_emission_names = ['eval']
        _check_if_in_options(data['data_options']['emission_probability']['name'], possible_emission_names)
        possible_xs_names = ['notyetavailable']
        _check_if_in_options(data['data_options']['cross_section']['name'], possible_xs_names) 
        possible_decay_names = ['chain', 'eval']
        _check_if_in_options(data['data_options']['half_life']['name'], possible_decay_names)
        possible_decay_spacings = ['linear']
        _check_if_in_options(data['data_options']['decay_time_spacing'], possible_decay_spacings)
        possible_concentration_options = ['CFY']
        _check_if_in_options(data['modeling_options']['concentration_handling'], possible_concentration_options)
        possible_irradiation_options = ['saturation', 'pulse']
        _check_if_in_options(data['modeling_options']['irrad_type'], possible_irradiation_options)

        if sum(data['data_options']['fissile_fractions'].values()) != 1.0:
            raise ValueError("Fissile fractions must sum to 1.0. Current sum: "
                             f"{sum(data['data_options']['fissile_fractions'].values())}")
        return
    
    def _adjust_data(self, data: dict) -> None:
        """
        Adjust the input data to fit desired formatting.

        Parameters
        ----------
        data : dict
            Dictionary containing settings and data selections.
        """
        # Adjust fission yield type and chain file suffix based on energy
        energy = data['data_options']['energy_MeV']
        # Two energy group division for selection of chain file
        if energy <= 1e-3:
            chain_suffix = 'pwr'
        elif energy > 1e-3:
            chain_suffix = 'sfr'
        chain_middle = data['data_options']['fission_yield']['data']
        
        data_selection = data['data_options']
        data_type = 'fission_yield'
        cur_data = data_selection[data_type]
        if cur_data['name'] == 'nfy':
            cur_data['name'] = 'nfy.csv'
        elif cur_data['name'] == 'chain':
            cur_data['name'] = 'chain_' + chain_middle + '_' + chain_suffix + '.csv'

        data_type = 'emission_probability'
        cur_data = data_selection[data_type]
        if cur_data['name'] == 'eval':
            cur_data['name'] = 'eval.csv'
        
        data_type = 'half_life'
        cur_data = data_selection[data_type]
        if cur_data['name'] == 'nfy':
            cur_data['name'] = 'nfy.csv'
        elif cur_data['name'] == 'chain':
            cur_data['name'] = 'chain_' + chain_middle + '_' + chain_suffix + '.csv'
        elif cur_data['name'] == 'eval':
            cur_data['name'] = 'eval.csv'


    

if __name__ == "__main__":
    handler = InputHandler("../../examples/keepin_1957/input.json")
    input_data = handler.read_input()
    handler.check_behaviour(input_data)