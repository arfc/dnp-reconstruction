from mosden.utils.input_handler import InputHandler
from pathlib import Path
from mosden.utils.csv_handler import CSVHandler

class BaseClass:
    def __init__(self, input_path: str) -> None:
        self.input_path: str = input_path
        self.input_handler: InputHandler = InputHandler(input_path)
        self.input_data: dict = self.input_handler.read_input()
        
        self.data_types: list[str] = ['fission_yield', 'half_life', 'cross_section', 'emission_probability']
        return None
    
    def _get_data_pathing(self, data_type: str) -> str:
        """
        Get the path for the specific type of data of interest

        Parameters
        ----------
        data_type : str
            The type of data to read (e.g., "fission_yield", "half_life", "cross_section", "emission_probability").

        Returns
        -------
        use_path : str
            Portion of the path that leads to the data
        use_name : str
            Name of the csv file
        """
        if data_type in self.data_types:
            use_data = self.input_data['data_options'][data_type]
        else:
            raise ValueError(f"Unknown data type: {data_type}")
        
        use_path = use_data['data']
        use_name = use_data['name']
        return use_path, use_name

    
    def _read_processed_data(self, fissile: str, data_type: str) -> dict[str, dict[str, float]]:
        """
        Read the processed data for a given fissile nuclide.

        Parameters
        ----------
        fissile : str
            The name of the fissile nuclide.
        data_type : str
            The type of data to read (e.g., "fission_yield", "half_life", "cross_section", "emission_probability").

        Returns
        -------
        dict
            The processed data for the fissile nuclide.
        """
        use_path, use_name = self._get_data_pathing(data_type)

        data_path = Path(self.data_path) / use_path / fissile / f"{self.energy}MeV" / use_name
        csv_handler = CSVHandler(data_path, create=False)
        if not csv_handler._file_exists():
            raise FileNotFoundError(f"Processed data file {data_path} does not exist.") 
        data = csv_handler.read_csv()
        return data