from mosden.countrate import CountRate
from pathlib import Path
from mosden.utils.csv_handler import CSVHandler
import numpy as np
import pytest
import os

@pytest.mark.parametrize("input_path, reference_output_path", [
    ("tests/integration/test-data/input1.json", "tests/integration/test-data/reference/test1"),
    ("tests/integration/test-data/input2.json", "tests/integration/test-data/reference/test2"),
    ("tests/integration/test-data/input3.json", "tests/integration/test-data/reference/test3")
] )
def test_calculate_count_rate(input_path, reference_output_path):
    """
    Test the count rate calculation.
    """
    countrate = CountRate(input_path)
    countrate.processed_data_dir = reference_output_path
    countrate.concentration_path = os.path.join(countrate.input_data['file_options']['output_dir'], 'concentrations.csv')

    # Generate count rate
    countrate.calculate_count_rate()
    
    # Check if output file is created
    output_path = Path(countrate.output_dir) / "count_rate.csv"
    assert output_path.exists(), f"Output file {output_path} does not exist."

    # Check if the output file has the expected content
    data = CSVHandler(output_path).read_countrate_csv()
    assert data, "Output file is empty."

    # Check all concentrations
    reference_path = Path(reference_output_path) / "count_rate.csv"
    reference_data = CSVHandler(reference_path).read_countrate_csv()

    assert data == reference_data, "Output count rate does not match the expected reference count rate."

    return