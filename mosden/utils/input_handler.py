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
        if data['file_options']['unprocessed_data_dir'] == data['file_options']['processed_data_dir']:
            raise ValueError("Unprocessed and processed data directories cannot be the same.")
        if data['file_options']['unprocessed_data_dir'] == data['file_options']['output_dir']:
            raise ValueError("Unprocessed data directory cannot be the same as the output directory.")
        if data['file_options']['processed_data_dir'] == data['file_options']['output_dir']:
            raise ValueError("Processed data directory cannot be the same as the output directory.")
        if data['modeling_options']['parent_feeding'] and not data['modeling_options']['concentration_handling'] == 'depletion':
            raise ValueError("Parent feeding option requires depletion method for concentration handling")
        
        possible_data_options = ['test-data', 'endfb71']
        if data['data_options']['half_life']['data'] not in possible_data_options:
            raise ValueError(f"Decay data option '{data['data_options']['half_life']['data']}' is not supported. "
                             f"Supported options are: {possible_data_options}")
        if data['data_options']['cross_section']['data'] not in possible_data_options:
            raise ValueError(f"Cross section option '{data['data_options']['cross_section']['data']}' is not supported. "
                             f"Supported options are: {possible_data_options}")
        if data['data_options']['fission_yield']['data'] not in possible_data_options:
            raise ValueError(f"Fission yield option '{data['data_options']['fission_yield']['data']}' is not supported. "
                             f"Supported options are: {possible_data_options}")

        possible_emission_probabilities = ['test-data', 'iaea']
        if data['data_options']['emission_probability']['data'] not in possible_emission_probabilities:
            raise ValueError(f"Emission probability option '{data['data_options']['emission_probability']['data']}' is not supported. "
                             f"Supported options are: {possible_emission_probabilities}")

        possible_fy_names = ['nfy', 'chain']
        if data['data_options']['fission_yield']['name'] not in possible_fy_names:
            raise ValueError(f"Fission yield name '{data['data_options']['fission_yield']['name']}' is not supported. "
                             f"Supported options are: {possible_fy_names}")
        
        possible_emission_names = ['eval']
        if data['data_options']['emission_probability']['name'] not in possible_emission_names:
            raise ValueError(f"Emission name '{data['data_options']['emission_probability']['name']}' is not supported. "
                             f"Supported options are: {possible_emission_names}")

        possible_xs_names = ['notyetavailable']
        if data['data_options']['cross_section']['name'] not in possible_xs_names:
            raise ValueError(f"Cross section name '{data['data_options']['cross_section']['name']}' is not supported. "
                             f"Supported options are: {possible_xs_names}")
        
        possible_decay_names = ['chain', 'eval']
        if data['data_options']['half_life']['name'] not in possible_decay_names:
            raise ValueError(f"Half life name '{data['data_options']['half_life']['name']}' is not supported. "
                             f"Supported options are: {possible_decay_names}")


        
        possible_decay_spacings = ['linear']
        if data['data_options']['decay_time_spacing'] not in possible_decay_spacings:
            raise ValueError(f"Decay spacing option '{data['data_options']['decay_time_spacing']}' is not supported. "
                             f"Supported options are: {possible_decay_spacings}")
        
        if sum(data['data_options']['fissile_fractions'].values()) != 1.0:
            raise ValueError("Fissile fractions must sum to 1.0. Current sum: "
                             f"{sum(data['data_options']['fissile_fractions'].values())}")
        

        possible_concentration_options = ['CFY']
        if data['modeling_options']['concentration_handling'] not in possible_concentration_options:
            raise ValueError(f"Concentration handling option '{data['modeling_options']['concentration_handling']}' is not supported. "
                             f"Supported options are: {possible_concentration_options}")
        
        possible_irradiation_options = ['saturation', 'pulse']
        if data['modeling_options']['irrad_type'] not in possible_irradiation_options:
            raise ValueError(f"Irradiation type '{data['modeling_options']['irrad_type']}' is not supported. "
                             f"Supported options are: {possible_irradiation_options}")
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
            
        if data['data_options']['fission_yield']['name'] == 'nfy':
            data['data_options']['fission_yield']['name'] = 'nfy.csv'
        elif data['data_options']['fission_yield']['name'] == 'chain':
            data['data_options']['fission_yield']['name'] = 'chain_' + chain_middle + '_' + chain_suffix + '.csv'
        
        if data['data_options']['emission_probability']['name'] == 'eval':
            data['data_options']['emission_probability']['name'] = 'eval.csv'

    

if __name__ == "__main__":
    handler = InputHandler("../../examples/keepin_1957/input.json")
    input_data = handler.read_input()
    handler.check_behaviour(input_data)