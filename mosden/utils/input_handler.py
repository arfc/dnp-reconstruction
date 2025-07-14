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
    
    def read_input(self, check=True) -> dict:
        """
        Read the input file and return the data as a dictionary.

        Parameters
        ----------
        check : bool, optional
            If True, checks the behaviour of the input data. Default is True.

        Returns
        -------
        output : dict
            Dictionary containing settings and data selections.
        """
        with open(self.input_path, 'r') as file:
            output = json.load(file)
        if check:
            self._check_behaviour(output)
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
        
        possible_data_options = ['endfb71']
        if data['data_options']['decay_constant'] not in possible_data_options:
            raise ValueError(f"Decay data option '{data['data_options']['decay_constant']}' is not supported. "
                             f"Supported options are: {possible_data_options}")

        if data['data_options']['cross_section'] not in possible_data_options:
            raise ValueError(f"Cross section option '{data['data_options']['cross_section']}' is not supported. "
                             f"Supported options are: {possible_data_options}")

        if data['data_options']['fission_yield']['data'] not in possible_data_options:
            raise ValueError(f"Fission yield option '{data['data_options']['fission_yield']['data']}' is not supported. "
                             f"Supported options are: {possible_data_options}")
        possible_fy_types = ['nfy', 'chain']
        if data['data_options']['fission_yield']['type'] not in possible_fy_types:
            raise ValueError(f"Fission yield type '{data['data_options']['fission_yield']['type']}' is not supported. "
                             f"Supported options are: {possible_fy_types}")

        
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

        return

    
if __name__ == "__main__":
    handler = InputHandler("../../examples/keepin_1957/input.json")
    input_data = handler.read_input()
    handler.check_behaviour(input_data)