from pathlib import Path
from mosden.preprocessing import Preprocess
import pytest
from mosden.utils.csv_handler import CSVHandler

@pytest.mark.parametrize("input_path, reference_output_path", [
    ("tests/integration/test-data/input1.json", "tests/integration/test-data/reference/test1/")
] )
def test_openmc_preprocess(input_path, reference_output_path):
    """
    Test the OpenMC preprocessing functionality.
    """
    preproc = Preprocess(str(input_path))
    preproc.openmc_preprocess()

    # Check if the output directory exists
    output_dir = Path(preproc.out_dir)
    assert output_dir.exists(), f"Output directory {output_dir} does not exist."
    
    # Check if the expected files are created
    for fissile in preproc.fissile_targets:
        path_nuc_energy = f"{fissile}/{preproc.energy_MeV}MeV"
        expected_files = list(output_dir.glob(f"{path_nuc_energy}/chain_*.csv"))
        assert expected_files, f"No CSV files found in {output_dir}/{path_nuc_energy}/chain_*.csv."

        # Check if the files contain expected data
        for file in expected_files:
            with open(file, 'r') as f:
                content = f.read()
                if fissile == 'U235':
                    assert 'Xe135' in content, f"'Xe135' not found in {file}."
                    assert '1.0' in content, f"'1.0' not found in {file}."
                    assert '100.0' in content, f"'100.0' not found in {file}."
                elif fissile == 'U238':
                    assert 'Xe135' in content, f"'Xe135' not found in {file}."
                    assert '0.123' in content, f"'0.123' not found in {file}."
                    assert '100.0' in content, f"'100.0' not found in {file}."
        
                # Check if files match expected output
                reference_file = Path(reference_output_path) / fissile / f"{preproc.energy_MeV}MeV" / file.name
                assert reference_file.exists(), f"Reference file {reference_file} does not exist."
                reference_data = CSVHandler(reference_file).read_csv()
                data = CSVHandler(file).read_csv()
                assert data == reference_data, f"Output data for {file} does not match reference data."

    return

@pytest.mark.parametrize("input_path, reference_output_path", [
    ("tests/integration/test-data/input1.json", "tests/integration/test-data/reference/test1/")
] )
def test_endf_preprocess(input_path, reference_output_path):
    """
    Test the ENDF preprocessing functionality.
    """
    preproc = Preprocess(str(input_path))
    preproc.endf_preprocess()

    # Check if the output directory exists
    output_dir = Path(preproc.out_dir)
    assert output_dir.exists(), f"Output directory {output_dir} does not exist."

    # Check if the expected files are created
    for fissile in preproc.fissile_targets:
        path_nuc_energy = f"{fissile}/{preproc.energy_MeV}MeV"
        expected_files = list(output_dir.glob(f"{path_nuc_energy}/nfy*.csv"))
        assert expected_files, f"No CSV files found in {output_dir}/{path_nuc_energy}/nfy*.csv."

        # Check if the files contain expected data
        for file in expected_files:
            with open(file, 'r') as f:
                lines = f.readlines()
                if fissile == 'U235':
                    line_num = 867
                    assert 'Xe135' in lines[line_num], f"'Xe135' not found on line {line_num + 1} in {file}."
                    assert '0.065385' in lines[line_num], f"'0.065385' not found on line {line_num + 1} in {file}."
                    assert '0.000457695' in lines[line_num], f"'0.000457695' not found on line {line_num + 1} in {file}."
                elif fissile == 'U238':
                    line_num = 911
                    assert 'Xe135' in lines[line_num], f"'Xe135' not found on line {line_num + 1} in {file}."
                    assert '0.0696758' in lines[line_num], f"'0.0696758' not found on line {line_num + 1} in {file}."
                    assert '0.000696758' in lines[line_num], f"'0.000696758' not found on line {line_num + 1} in {file}."

                # Check if files match expected output
                reference_file = Path(reference_output_path) / fissile / f"{preproc.energy_MeV}MeV" / file.name
                assert reference_file.exists(), f"Reference file {reference_file} does not exist."
                reference_data = CSVHandler(reference_file).read_csv()
                data = CSVHandler(file).read_csv()
                assert data == reference_data, f"Output data for {file} does not match reference data." 
    return

@pytest.mark.parametrize("input_path, reference_output_path", [
    ("tests/integration/test-data/input1.json", "tests/integration/test-data/reference/test1/")
] )
def test_iaea_preprocess(input_path, reference_output_path):
    """
    Test the IAEA preprocessing functionality.
    """
    preproc = Preprocess(str(input_path))
    preproc.iaea_preprocess()

    # Check if the output directory exists
    output_dir = Path(preproc.out_dir)
    assert output_dir.exists(), f"Output directory {output_dir} does not exist."

    # Check if the expected files are created
    iaea_dir = output_dir / 'iaea'
    assert iaea_dir.exists(), f"IAEA directory {iaea_dir} does not exist."
    
    csv_path = iaea_dir / 'eval.csv'
    assert csv_path.exists(), f"CSV file {csv_path} does not exist."

    # Check if the CSV file contains expected data
    data = CSVHandler(csv_path).read_csv()
    assert '87BR' in data, "'87BR' not found in IAEA data."
    assert '11XX' in data, "'11XX' not found in IAEA data."
    assert '12XX' in data, "'12XX' not found in IAEA data."
    assert '13XX' in data, "'13XX' not found in IAEA data."
    assert '14XX' in data, "'14XX' not found in IAEA data."
    
    # Check if the emission probabilities are calculated correctly
    assert 'emission probability' in data['11XX'], "'emission probability' not found for 11XX."
    assert 'sigma emission probability' in data['11XX'], "'sigma emission probability' not found for 11XX."
    assert 'emission probability' in data['12XX'], "'emission probability' not found for 12XX."
    assert 'sigma emission probability' in data['12XX'], "'sigma emission probability' not found for 12XX."
    assert 'emission probability' in data['13XX'], "'emission probability' not found for 13XX."
    assert 'sigma emission probability' in data['13XX'], "'sigma emission probability' not found for 13XX."
    assert 'emission probability' in data['14XX'], "'emission probability' not found for 14XX."
    assert 'sigma emission probability' in data['14XX'], "'sigma emission probability' not found for 14XX."

    # Check if the values match the expected output
    assert data['11XX']['emission probability'] == 1.0, "Emission probability for 11XX does not match expected value."
    assert data['12XX']['emission probability'] == 100.0, "Emission probability for 12XX does not match expected value."
    assert data['13XX']['emission probability'] == 55.64, "Emission probability for 13XX does not match expected value."
    assert data['14XX']['emission probability'] == 55.64, "Emission probability for 14XX does not match expected value."

    return

