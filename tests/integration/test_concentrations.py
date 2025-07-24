from mosden.concentrations import Concentrations
from pathlib import Path
from mosden.utils.csv_handler import CSVHandler
import numpy as np
import pytest

@pytest.mark.parametrize("input_path, reference_output_path", [
    ("tests/integration/test-data/input1.json", "tests/integration/test-data/reference/test1"),
    ("tests/integration/test-data/input2.json", "tests/integration/test-data/reference/test2"),
    ("tests/integration/test-data/input3.json", "tests/integration/test-data/reference/test3"),
    ("tests/integration/test-data/input4.json", "tests/integration/test-data/reference/test4")
] )
def test_generate_concentrations(input_path, reference_output_path):
    """
    Test the concentration generation method.
    """
    concentrations = Concentrations(input_path)
    concentrations.processed_data_dir = reference_output_path
    
    concentrations.generate_concentrations()
    
    output_path = Path(concentrations.output_dir) / "concentrations.csv"
    assert output_path.exists(), f"Output file {output_path} does not exist."

    data = CSVHandler(output_path).read_csv()
    assert data, "Output file is empty."
    assert "Concentration" in data['Xe135'], "Concentration data is missing."
    assert "sigma Concentration" in data['Xe135'], "Sigma Concentration data is missing."
    
    for nuclide, values in data.items():
        assert isinstance(values['Concentration'], float), f"Concentration for {nuclide} is not a float."
        assert isinstance(values['sigma Concentration'], float), f"Sigma Concentration for {nuclide} is not a float."

    reference_path = Path(reference_output_path) / "concentrations.csv"
    reference_data = CSVHandler(reference_path).read_csv()

    assert data == reference_data, "Output concentrations do not match the expected reference concentrations."

    return