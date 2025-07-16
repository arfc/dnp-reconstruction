from mosden.concentrations import Concentrations

def test_concentrations_init():
    """
    Test the initialization of the Concentrations class.
    """

    input_path = './tests/unit/input/input.json'
    conc = Concentrations(input_path)
    assert conc.input_path == input_path, f"Expected input path {input_path}, but got {conc.input_path}"
    assert conc.output_dir == './tests/unit/output', f"Expected output directory './tests/unit/output', but got {conc.output_dir}"
    assert conc.energy_MeV == 1.0, f"Expected energy 1.0, but got {conc.energy_MeV}"
    assert conc.fissiles == {'U235': 0.8, 'U238': 0.2}, f"Expected fissile targets {{'U235': 0.8, 'U238': 0.2}}, but got {conc.fissiles}"

    return