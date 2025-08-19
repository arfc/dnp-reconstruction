import os
from pathlib import Path
import json
import itertools
import shutil
from copy import deepcopy
import subprocess

base_input_file = './input.json'
analysis_list = list()

residence_time_analysis = {
    'meta': {
        'name': 'residence_time',
        'run': True,
        'overwrite': True
    },
    'incore_s': [1, 2],
    'excore_s': [1, 2],
}
analysis_list.append(residence_time_analysis)












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

def create_directory(analysis: dict) -> str:
    dir_name = f'./{analysis["meta"]["name"]}'
    if os.path.isdir(dir_name) and analysis['meta']['overwrite']:
        shutil.rmtree(dir_name)
    os.makedirs(dir_name, exist_ok=analysis['meta']['overwrite'])
    return dir_name

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

        filename = f"{analysis['meta']['name']}_{idx}.json"
        file_path = dir_path / filename

        with open(file_path, "w") as f:
            json.dump(new_data, f, indent=2)

        paths.append(str(file_path))

    return paths

def run_mosden(input_paths: list[str]):
    command = ['mosden', '-a'] + input_paths
    subprocess.run(command)

if __name__ == '__main__':
    for analysis in analysis_list:
        # Create directory for each analysis
        dir_path = create_directory(analysis)
        # Populate input files
        input_paths = populate_inputs(analysis, dir_path)
        # Run MoSDeN and MultiPostProc
        run_mosden(input_paths)
    pass