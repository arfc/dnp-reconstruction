import os
from pathlib import Path
import json
import itertools
import shutil
from copy import deepcopy
import subprocess
from mosden.utils.chemical_schemes import MSBR_scheme

base_input_file = './input.json'
analysis_list = list()

name = 'tintex'
residence_time_analysis = {
    'meta': {
        'name': name,
        'run_full': False,
        'run_post': False,
        'overwrite': True,
    },
    'incore_s': [0.1, 1, 10, 100],
    'excore_s': [1],
    'multi_id': [name]
}
analysis_list.append(residence_time_analysis)

name = 'chem_long'
chemical_long_analysis = {
    'meta': {
        'name': name,
        'run_full': False,
        'run_post': False,
        'overwrite': True
    },
    'reprocessing': [MSBR_scheme(),
                     MSBR_scheme(include_long=False)],
    'incore_s': [10],
    'excore_s': [10],
    'multi_id': [name]
}
analysis_list.append(chemical_long_analysis)

name = 'chem_bool'
chemical_bool_analysis = {
    'meta': {
        'name': name,
        'run_full': False,
        'run_post': False,
        'overwrite': True
    },
    'reprocessing': [MSBR_scheme(),
                     MSBR_scheme(rate_scaling=0.0)],
    'incore_s': [10],
    'excore_s': [10],
    'multi_id': [name]
}
analysis_list.append(chemical_bool_analysis)

name = 'spacing_times'
spacing_times_analysis = {
    'meta': {
        'name': name,
        'run_full': False,
        'run_post': False,
        'overwrite': True
    },
    'decay_time_spacing': ['log', 'linear'],
    'multi_id': [name]
}
analysis_list.append(spacing_times_analysis)

name = 'decay_time_nodes'
decay_times_analysis = {
    'meta': {
        'name': name,
        'run_full': False,
        'run_post': False,
        'overwrite': True
    },
    'num_decay_times': [50, 100, 200, 400, 800],
    'multi_id': [name]
}
analysis_list.append(decay_times_analysis)

name = 'total_decay_time'
total_decay_analysis = {
    'meta': {
        'name': name,
        'run_full': False,
        'run_post': False,
        'overwrite': True
    },
    'decay_time': [150, 300, 600, 1200, 2400, 4800],
    'multi_id': [name]
}
analysis_list.append(total_decay_analysis)

name = 'detailed_decay'
detailed_decay_analysis = {
    'meta': {
        'name': name,
        'run_full': True,
        'run_post': False,
        'overwrite': True
    },
    'decay_time': [1200, 2400],
    'num_decay_times': [800, 1600],
    'multi_id': [name]
}
analysis_list.append(detailed_decay_analysis)



def replace_value(input_data: dict, key: str, new_val) -> bool:
    """
    Recursively search through dict d to find key and replace its value.
    Returns True if replacement was made, False otherwise.
    """
    if key in input_data.keys():
        input_data[key] = new_val
        return True

    for val in input_data.values():
        if isinstance(val, dict):
            if replace_value(val, key, new_val):
                return True
    return False

def create_directory(dir_name) -> None:
    if os.path.isdir(dir_name) and analysis['meta']['overwrite']:
        shutil.rmtree(dir_name)
    os.makedirs(dir_name, exist_ok=analysis['meta']['overwrite'])
    return None

def populate_inputs(analysis: dict, dir_path: str) -> list[str]:
    paths = []
    dir_path = Path(dir_path)

    with open("input.json", "r") as f:
        base_data = json.load(f)

    component_keys = [k for k in analysis.keys() if k != "meta"]
    value_combinations = itertools.product(*(analysis[k] for k in component_keys))
    for idx, combination in enumerate(value_combinations, start=1):
        new_data = deepcopy(base_data)

        for key, val in zip(component_keys, combination):
            replaced = replace_value(new_data, key, val)
            if not replaced:
                raise KeyError(f'{key} not found in input file')

        filename = 'input.json'
        file_dir = dir_path / str(idx)
        file_path = file_dir / filename
        new_data['file_options']['processed_data_dir'] = str(file_dir)
        new_data['file_options']['output_dir'] = str(file_dir)
        new_data['file_options']['log_file'] = str(file_dir) + '/log.log'
        new_data['name'] = str(combination)
        if analysis['meta']['run_full']:
            create_directory(file_dir)

        with open(file_path, "w") as f:
            json.dump(new_data, f, indent=2)

        paths.append(str(file_path))
    return paths

def run_mosden(analysis: dict, input_paths: list[str]):
    if analysis['meta']['run_full']:
        command = ['mosden', '-a'] + input_paths
    elif analysis['meta']['run_post']:
        command = ['mosden', '-post'] + input_paths
    else:
        return None
    subprocess.run(command)
    return None

if __name__ == '__main__':
    for analysis in analysis_list:
        dir_name = f'./{analysis["meta"]["name"]}'
        if analysis['meta']['run_full'] or analysis['meta']['run_post']:
            input_paths = populate_inputs(analysis, dir_name)
            run_mosden(analysis, input_paths)
    pass