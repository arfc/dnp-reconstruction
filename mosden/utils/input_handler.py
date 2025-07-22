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
        self.preproc_choices: dict = dict()
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
        if data['modeling_options']['parent_feeding'] and not data['modeling_options']['concentration_handling'] == 'depletion':
            raise ValueError("Parent feeding option requires depletion method for concentration handling")
        
        possible_decay_spacings = ['linear']
        _check_if_in_options(data['data_options']['decay_time_spacing'], possible_decay_spacings)
        possible_concentration_options = ['CFY', 'IFY']
        _check_if_in_options(data['modeling_options']['concentration_handling'], possible_concentration_options)
        possible_irradiation_options = ['saturation', 'pulse']
        _check_if_in_options(data['modeling_options']['irrad_type'], possible_irradiation_options)
        possible_count_rate_options = ['data']
        _check_if_in_options(data['modeling_options']['count_rate_handling'], possible_count_rate_options)

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
        return None




    

if __name__ == "__main__":
    handler = InputHandler("../../examples/keepin_1957/input.json")
    input_data = handler.read_input()
    handler.check_behaviour(input_data)