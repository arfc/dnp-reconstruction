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
        expected_files = list(output_dir.glob(f"{path_nuc_energy}/chain_*.csv"))
        assert expected_files, f"No CSV files found in {path_nuc_energy}."

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

    shutil.rmtree(output_dir, ignore_errors=True)
    return

def test_endf_preprocess():
    """
    Test the ENDF preprocessing functionality.
    """
    input_path = Path(__file__).parent / 'inputs' / 'pre_input.json'
    preproc = Preprocess(str(input_path))
    preproc.endf_preprocess()

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
                if fissile == 'U235':
                    assert 'Xe135' in content, f"'Xe135' not found in {file}."
                    assert '0.065385' in content, f"'0.065385' not found in {file}."
                    assert '0.000457695' in content, f"'0.000457695' not found in {file}."
                elif fissile == 'U238':
                    assert 'Xe135' in content, f"'Xe135' not found in {file}."
                    assert '0.123' in content, f"'0.123' not found in {file}."
                    assert '100.0' in content, f"'100.0' not found in {file}."

    shutil.rmtree(output_dir, ignore_errors=True)
    return
