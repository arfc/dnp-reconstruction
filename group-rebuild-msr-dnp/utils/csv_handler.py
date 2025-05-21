import pandas as pd

class CSVHandler:
    def __init__(self, file_path: str) -> None:
        self.file_path = file_path
        return None

    def read_csv(self) -> dict[str, dict[str, float]]:
        df = pd.read_csv(self.file_path)
        data = df.to_dict(orient='index')
        return data

    def write_csv(self, data: dict[str, dict[str, float]]) -> None:
        df = pd.DataFrame.from_dict(data, orient='index')
        df.to_csv(self.file_path, index=False)
        return None