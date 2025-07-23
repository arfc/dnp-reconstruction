from mosden.groupfit import Grouper
import numpy as np
import pytest

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

@pytest.mark.slow
def test_grouper_pulse_fitting():
    """
    Test the nonlinear least squares fitting for a pulse case
    """
    input_path = './tests/unit/input/input.json'
    grouper = Grouper(input_path)
    count_data: dict[str: np.ndarray[float]] = dict()

    half_lives = [60, 50, 40, 30, 20, 10]
    yields = [0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    times = np.linspace(0, 600, 100)
    counts = np.zeros(len(times))
    for group in range(len(half_lives)):
        lam = np.log(2) / half_lives[group]
        a = yields[group]
        counts += a * lam * np.exp(-lam * times)
    
    parameters = yields + half_lives
    func_counts = grouper._pulse_fit_function(times, parameters)
    assert np.isclose(func_counts, counts).all(), 'Pulse counts mismatch between hand calculation and function evaluation'

    grouper.irrad_type = 'pulse'
    grouper.num_groups = 6
    count_data['times'] = times
    count_data['counts'] = counts
    count_data['sigma counts'] = np.zeros(len(counts))
    data = grouper._nonlinear_least_squares(count_data=count_data)
    test_yields = [data[key]['yield'] for key in data.keys()]
    test_half_lives = [data[key]['half_life'] for key in data.keys()]
    test_yields.sort(reverse=True)
    test_half_lives.sort(reverse=True)
    for group in range(grouper.num_groups):
        assert np.isclose(test_yields[group], yields[group]), f'Group {group+1} yields - {test_yields[group]=} != {yields[group]=}'
        assert np.isclose(test_half_lives[group], half_lives[group]), f'Group {group+1} half lives mismatch - {test_yields[group]=} != {half_lives[group]=}'
    
    
@pytest.mark.slow
def test_grouper_saturation_fitting():
    """
    Test the nonlinear least squares fitting for a saturation case
    """
    input_path = './tests/unit/input/input.json'
    grouper = Grouper(input_path)
    count_data: dict[str: np.ndarray[float]] = dict()
    grouper.t_ex = 0

    half_lives = [60, 50, 40, 30, 20, 10]
    yields = [0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    times = np.linspace(0, 600, 100)
    counts = np.zeros(len(times))
    for group in range(len(half_lives)):
        lam = np.log(2) / half_lives[group]
        a = yields[group]
        counts += a * np.exp(-lam * times)
    
    parameters = yields + half_lives
    func_counts = grouper._saturation_fit_function(times, parameters)
    assert np.isclose(func_counts, counts).all(), 'Pulse counts mismatch between hand calculation and function evaluation'

    grouper.irrad_type = 'saturation'
    grouper.num_groups = 6
    count_data['times'] = times
    count_data['counts'] = counts
    count_data['sigma counts'] = np.zeros(len(counts))
    data = grouper._nonlinear_least_squares(count_data=count_data)
    test_yields = [data[key]['yield'] for key in data.keys()]
    test_half_lives = [data[key]['half_life'] for key in data.keys()]
    test_yields.sort(reverse=True)
    test_half_lives.sort(reverse=True)
    for group in range(grouper.num_groups):
        assert np.isclose(test_yields[group], yields[group]), f'Group {group+1} yields - {test_yields[group]=} != {yields[group]=}'
        assert np.isclose(test_half_lives[group], half_lives[group]), f'Group {group+1} half lives mismatch - {test_yields[group]=} != {half_lives[group]=}'