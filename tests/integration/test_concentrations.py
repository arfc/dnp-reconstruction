from mosden.concentrations import Concentrations
from pathlib import Path
from mosden.utils.csv_handler import CSVHandler
import numpy as np
import pytest

@pytest.mark.parametrize("input_path, reference_output_path", [
    ("tests/integration/test-data/input1.json", "tests/integration/test-data/output/test1/")
] )
def test_generate_CFY_concentrations(input_path, reference_output_path):
    """
    Test the CFY concentration generation method.
    """
    concentrations = Concentrations(input_path)
    
    # Generate concentrations
    concentrations.generate_concentrations()
    
    # Check if output file is created
    output_path = Path(concentrations.output_dir) / "concentrations.csv"
    assert output_path.exists(), f"Output file {output_path} does not exist."

    # Check if the output file has the expected content
    data = CSVHandler(output_path).read_csv()
    assert data, "Output file is empty."
    assert "Concentration" in data['Xe135'], "Concentration data is missing."
    assert "sigma Concentration" in data['Xe135'], "Sigma Concentration data is missing."
    
    # Check if the concentrations are calculated correctly
    for nuclide, values in data.items():
        assert isinstance(values['Concentration'], float), f"Concentration for {nuclide} is not a float."
        assert isinstance(values['sigma Concentration'], float), f"Sigma Concentration for {nuclide} is not a float."
    
    # Check if the concentrations are correctly calculated for Xe135
    assert data['Xe135']['Concentration'] > 0, "Concentration for Xe135 should be greater than 0."
    assert data['Xe135']['sigma Concentration'] >= 0, "Sigma Concentration for Xe135 should be non-negative."

    assert np.isclose(data['Xe135']['Concentration'], 0.06624316), "Concentration for Xe135 should be 0.06624316"
    assert np.isclose(data['Xe135']['sigma Concentration'], 0.000391777), "Sigma Concentration for Xe135 should be 0.000391777"

    # Check all concentrations
    reference_path = Path(reference_output_path) / "concentrations.csv"
    reference_data = CSVHandler(reference_path).read_csv()

    assert data == reference_data, "Output concentrations do not match the expected reference concentrations."

    return