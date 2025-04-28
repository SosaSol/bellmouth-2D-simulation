import os
import glob
import numpy as np
import pandas as pd

### Functions

def read_last_value(filepath: str) -> float:
    """Read the last value from a file, ignoring comments and empty lines."""
    with open(filepath, 'r') as f:
        lines = [line for line in f if not line.startswith('#') and line.strip()]
    if not lines:
        raise ValueError(f"No valid data found in {filepath}")
    return float(lines[-1].split()[-1])

def read_yplus(filepath: str) -> float:
    """Read the last Max yPlus value from a file, ignoring comments and empty lines."""
    with open(filepath, 'r') as f:
        lines = [line for line in f if not line.startswith('#') and line.strip()]
    if not lines:
        raise ValueError(f"No valid data found in {filepath}")
    return float(lines[-1].split()[-2])

def extract_values(folder_path: str, folder_name: str, geometry_type: str) -> tuple[float, float, float] | None:
    """Extract head loss, mass flow, and max yPlus from post-processing files."""
    try:
        inlet_file  = glob.glob(os.path.join(folder_path, 'postProcessing/averageTotalPressureInlet1/0/surfaceFieldValue.dat'))[0]
        outlet_file = glob.glob(os.path.join(folder_path, 'postProcessing/averageTotalPressureOutlet1/0/surfaceFieldValue.dat'))[0]
        massflow_file = glob.glob(os.path.join(folder_path, 'postProcessing/massFlowRateOutlet1/0/surfaceFieldValue.dat'))[0]
        yplus_file = glob.glob(os.path.join(folder_path, 'postProcessing/yPlus1/0/yPlus.dat'))[0]

        p_inlet = read_last_value(inlet_file)
        p_outlet = read_last_value(outlet_file)
        massflow = read_last_value(massflow_file)
        max_yplus = read_yplus(yplus_file)

        delta_p_total = p_inlet - p_outlet

        # Reporting
        print(f"Geometry: {folder_name} ({geometry_type})")
        print(f"  - Mass flow rate (kg/s)   : {massflow:.4f}")
        print(f"  - Head loss (Pa)          : {delta_p_total:.2f}")
        print(f"  - Max yPlus               : {max_yplus:.2f}")

        if max_yplus > 5.0:
            print("  WARNING: yPlus > 5 detected!")
        elif max_yplus >= 1.0:
            print("  WARNING: yPlus not fully below 1!")
        else:
            print("  OK: yPlus all below 1.")
        print("-" * 60)

        return delta_p_total, massflow, max_yplus

    except (IndexError, FileNotFoundError) as e:
        print(f"Missing or incomplete post-processing data in {folder_name}: {e}")
    except Exception as e:
        print(f"Error processing {folder_name}: {e}")
    return None

def ensure_dir(path: str) -> None:
    """Ensure that a directory exists, creating it if necessary."""
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Created directory: {path}")

### Main script

outputs_dir = 'outputs'
postprocessing_dir = 'postProcessingOutputs'
ensure_dir(postprocessing_dir)

if not os.path.exists(outputs_dir):
    print(f"Error: '{outputs_dir}' folder not found.")
    exit(1)

folders = [f for f in os.listdir(outputs_dir) if os.path.isdir(os.path.join(outputs_dir, f))]
print(f"Found {len(folders)} geometries.\n")

results = {}

for folder in folders:
    folder_path = os.path.join(outputs_dir, folder)

    if folder.startswith('ELL-'):
        geometry_type = 'ELL'
        parts = folder.split('-')
        if len(parts) != 8:
            print(f"Unexpected ELL folder format: {folder}. Skipping.")
            continue
        _, Mw, Mb, Kx, Ky, r, L, t = parts
        Mb, Mw = int(Mb), int(Mw)
        params_key = f"{Kx}_{Ky}_{r}_{L}"

    elif folder.startswith('IN-'):
        geometry_type = 'IN'
        parts = folder.split('-')
        if len(parts) != 5:
            print(f"Unexpected IN folder format: {folder}. Skipping.")
            continue
        _, Mw, L, t = parts
        Mb = 1
        params_key = f"{L}"

    elif folder.startswith('RAD-'):
        geometry_type = 'RAD'
        parts = folder.split('-')
        if len(parts) != 6:
            print(f"Unexpected RAD folder format: {folder}. Skipping.")
            continue
        _, Mw, Mb, K, L, t = parts
        Mb, Mw = int(Mb), int(Mw)
        params_key = f"{K}_{L}"

    else:
        print(f"Unknown geometry type in folder {folder}. Skipping.")
        continue

    result = extract_values(folder_path, folder, geometry_type)
    if result is None:
        print(f"Skipping folder {folder} due to data issues.\n")
        continue

    delta_p_total, massflow, max_yplus = result
    key = (geometry_type, params_key)
    results.setdefault(key, {})[(Mb, Mw)] = {
        'head_loss': delta_p_total,
        'massflow': massflow,
        'max_yplus': max_yplus
    }

# Save results into labeled DataFrames
for (geom, param_key), data in results.items():
    mbmw_pairs = list(data.keys())
    Mb_values = sorted(set(mb for mb, mw in mbmw_pairs))
    Mw_values = sorted(set(mw for mb, mw in mbmw_pairs))

    # Create empty DataFrames
    head_loss_df = pd.DataFrame(index=Mb_values, columns=Mw_values, dtype=float)
    massflow_df  = pd.DataFrame(index=Mb_values, columns=Mw_values, dtype=float)
    yplus_df     = pd.DataFrame(index=Mb_values, columns=Mw_values, dtype=float)

    for (mb, mw), vals in data.items():
        head_loss_df.loc[mb, mw] = vals['head_loss']
        massflow_df.loc[mb, mw] = vals['massflow']
        yplus_df.loc[mb, mw] = vals['max_yplus']

    base_filename = f"{geom}_{param_key}"
    head_loss_df.to_csv(os.path.join(postprocessing_dir, f"head_losses_{base_filename}.csv"))
    massflow_df.to_csv(os.path.join(postprocessing_dir, f"massflow_{base_filename}.csv"))
    yplus_df.to_csv(os.path.join(postprocessing_dir, f"max_yplus_{base_filename}.csv"))

print("\nPost-processing complete. CSV files saved in 'postProcessingOutputs/'.")
