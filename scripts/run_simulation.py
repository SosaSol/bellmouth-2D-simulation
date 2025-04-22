import yaml
import itertools
import os
from copy import deepcopy
from your_module import run_one_simulation  # You should define this

def load_config(path):
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def generate_output_dir(params):
    fname = f"ELL-{params['Mw']}-{params['Mb']}-{int(params['Kx'] * 100)}-{int(params['Ky'] * 100)}-{int(params['r'] * 1e3)}-{int(params['t'] * 1e3)}"
    return os.path.join("outputs", fname)

def is_sweep(param):
    return isinstance(param, dict) and 'sweep' in param

def expand_sweep_params(config):
    sweep_keys = []
    sweep_values = []

    geom = config['geometry']
    for key, value in geom.items():
        if is_sweep(value):
            sweep_keys.append(key)
            sweep_values.append(value['sweep'])

    for combo in itertools.product(*sweep_values):
        new_config = deepcopy(config)
        for i, key in enumerate(sweep_keys):
            new_config['geometry'][key] = combo[i]

        new_config['simulation']['output_dir'] = generate_output_dir(new_config['geometry'])
        yield new_config

def main():
    config = load_config("sweep_config.yaml")  # or sys.argv[1] for CLI input

    if any(is_sweep(v) for v in config['geometry'].values()):
        print("Running parameter sweep...")
        for sim_config in expand_sweep_params(config):
            os.makedirs(sim_config['simulation']['output_dir'], exist_ok=True)
            run_one_simulation(sim_config)
    else:
        print("Running single simulation...")
        config['simulation']['output_dir'] = generate_output_dir(config['geometry'])
        os.makedirs(config['simulation']['output_dir'], exist_ok=True)
        run_one_simulation(config)

if __name__ == "__main__":
    main()
