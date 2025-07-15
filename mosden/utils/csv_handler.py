import pandas as pd
import os

class CSVHandler:
    def __init__(self, file_path: str, overwrite: bool=False, create=True) -> None:
        self.file_path = file_path
        if create:
            self._create_directory() 
        self.overwrite = overwrite
        return None
    
    def _create_directory(self) -> None:
        """
        Create the directory for the file path if it does not exist.
        """
        directory = os.path.dirname(self.file_path)
        if not os.path.exists(directory):
            os.makedirs(directory)
        return None
    
    def _file_exists(self) -> bool:
        """
        Check if the file exists at the specified path.

        Returns
        -------
        bool
            True if the file exists, False otherwise.
        """
        return os.path.isfile(self.file_path)

    def read_csv(self) -> dict[str, dict[str, float]]:
        df = pd.read_csv(self.file_path, index_col=0)
        data = df.to_dict(orient='index')
        return data

    def write_csv(self, data: dict[str, dict[str, float]]) -> None:
        if not self.overwrite and self._file_exists():
            raise FileExistsError(f"File {self.file_path} already exists. Set overwrite=True to overwrite.")
        df = pd.DataFrame.from_dict(data, orient='index')
        df.index.name = 'Nuclide'
        df.to_csv(self.file_path, index=True)
        return None