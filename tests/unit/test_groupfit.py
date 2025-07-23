from mosden.groupfit import Grouper

def test_grouper_init():
    """
    Test the initialization of the grouper class
    """
    input_path = './tests/unit/input/input.json'
    countrate = Grouper(input_path)
    assert countrate.input_path == input_path, f"Expected input path {input_path}, but got {countrate.input_path}"
    assert countrate.output_dir == './tests/unit/output', f"Expected output directory './tests/unit/output', but got {countrate.output_dir}"
    assert countrate.energy_MeV == 1.0, f"Expected energy 1.0, but got {countrate.energy_MeV}"
    assert countrate.fissiles == {'U235': 0.8, 'U238': 0.2}, f"Expected fissile targets {{'U235': 0.8, 'U238': 0.2}}, but got {countrate.fissiles}"

    return
