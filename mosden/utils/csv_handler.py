import pandas as pd
import os
from uncertainties import ufloat

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

    def read_csv(self, raw_iaea=False) -> dict[str, dict[str, float]]:
        if raw_iaea:
            return self._read_iaea_csv()
        df = pd.read_csv(self.file_path, index_col=0)
        data = df.to_dict(orient='index')
        return data
    
    def _read_iaea_csv(self) -> dict[str, dict[str, float]]:
        """
        Read the IAEA CSV file and return the data as a dictionary.

        Returns
        -------
        dict
            The data read from the IAEA CSV file.
        """
        data: dict[str, dict[str, float]] = {}
        df = pd.read_csv(self.file_path, header=1)
        for _, row in df.iterrows():
            iaea_nuc = row['nucid']
            nuc = self._iaea_to_mosden_nuc(iaea_nuc)
            half_life = row[' T1/2 [s] ']
            half_life_uncertainty = row[' D T1/2 [s]']
            if row['D pn1'] < 0:
                row['D pn1'] = 0
            if row['D pn2 '] < 0:
                row['D pn2 '] = 0
            if row['D pn3'] < 0:
                row['D pn3'] = 0
            P1 = ufloat(row[' pn1 % ']/100, row['D pn1']/100)
            P2 = ufloat(row[' pn2 % ']/100, row['D pn2 ']/100)
            P3 = ufloat(row[' pn3 % ']/100, row['D pn3']/100)
            prob_beta = ufloat(row['  beta- %']/100, row[' D beta-']/100)
            emission_prob = (1*P1 + 2*P2 + 3*P3) * prob_beta
            data[nuc] = {}
            data[nuc]['half_life'] = half_life
            data[nuc]['sigma half_life'] = half_life_uncertainty
            data[nuc]['emission probability'] = emission_prob.n
            data[nuc]['sigma emission probability'] = emission_prob.s
        return data
    
    def _iaea_to_mosden_nuc(self, iaea_nuc: str) -> str:
        """
        Converts IAEA nuclide format to MoSDeN format.

        Parameters
        ----------
        iaea_nuc : str
            IAEA nuclide identifier.

        Returns
        -------
        str
            MoSDeN formatted nuclide identifier.
        """
        i = 0
        while i < len(iaea_nuc) and iaea_nuc[i].isdigit():
            i += 1
        mass = iaea_nuc[:i]
        element = iaea_nuc[i:].capitalize()
        return f"{element}{mass}"

    def write_csv(self, data: dict[str, dict[str, float]]) -> None:
        if not self.overwrite and self._file_exists():
            raise FileExistsError(f"File {self.file_path} already exists. Set overwrite=True to overwrite.")
        df = pd.DataFrame.from_dict(data, orient='index')
        df.index.name = 'Nuclide'
        df.to_csv(self.file_path, index=True)
        return None