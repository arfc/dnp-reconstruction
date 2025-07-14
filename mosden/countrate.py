from mosden.utils.input_handler import InputHandler
from mosden.utils.csv_handler import CSVHandler


class CountRate:
    """
    Class to handle the delayed neutron count rate calculations
    """
    def __init__(self, input_path: str) -> None:
        self.input_path: str = input_path
        self.input_handler: InputHandler = InputHandler(input_path)
        self.input_data: dict = self.input_handler.read_input()
        self.output_dir: str = self.input_data['file_options']['output_dir']
        
        return None
    
    