from pathlib import Path
from mosden.preprocessing import Preprocess
import shutil

def test_openmc_preprocess():
    """
    Test the OpenMC preprocessing functionality.
    """
    input_path = Path(__file__).parent / 'inputs' / 'pre_input.json'
    preproc = Preprocess(str(input_path))
    preproc.openmc_preprocess()
    # Check if the output directory exists
    output_dir = Path(preproc.out_dir)
    assert output_dir.exists(), f"Output directory {output_dir} does not exist."
    
    # Check if the expected files are created
    for fissile in preproc.fissile_targets:
        path_nuc_energy = f"{fissile}/{preproc.energy_MeV}MeV"
        expected_files = list(output_dir.glob(f"{path_nuc_energy}/*.csv"))
        assert expected_files, f"No CSV files found in {path_nuc_energy}."

    # Check if the files contain expected data
    for file in expected_files:
        with open(file, 'r') as f:
            content = f.read()
            assert 'Xe135' in content, f"'Xe135' not found in {file}."
            assert '1.0' in content, f"'1.0' not found in {file}."
            assert '100.0' in content, f"'100.0' not found in {file}."

    # Clean up the output directory after the test
    for file in expected_files:
        file.unlink()
    shutil.rmtree(output_dir, ignore_errors=True)
    print("OpenMC preprocessing test passed successfully.")
    return