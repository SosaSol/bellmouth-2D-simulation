import yaml
import json
import os

def load_config(config_file):
    """
    Load configuration from a YAML or JSON file and return it as a Python dictionary.
    """
    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Config file {config_file} not found")
    
    with open(config_file, 'r') as file:
        if config_file.endswith('.yaml') or config_file.endswith('.yml'):
            config = yaml.safe_load(file)
        elif config_file.endswith('.json'):
            config = json.load(file)
        else:
            raise ValueError("Unsupported config file format. Use YAML or JSON.")
    
    return config

def save_config(config, config_file):
    """
    Save a Python dictionary as a YAML or JSON config file.
    """
    with open(config_file, 'w') as file:
        if config_file.endswith('.yaml') or config_file.endswith('.yml'):
            yaml.dump(config, file)
        elif config_file.endswith('.json'):
            json.dump(config, file, indent=4)
        else:
            raise ValueError("Unsupported config file format. Use YAML or JSON.")
    print(f"Config saved to {config_file}")

def merge_configs(base_config, override_config):
    """
    Merge two configuration dictionaries. Override values from the second config.
    """
    merged_config = base_config.copy()
    merged_config.update(override_config)
    return merged_config

def validate_config(config, required_keys):
    """
    Validate the config by checking if all required keys are present.
    """
    missing_keys = [key for key in required_keys if key not in config]
    if missing_keys:
        raise ValueError(f"Missing required keys in config: {', '.join(missing_keys)}")
    print("Config validated successfully")
