import os
import glob
import numpy as np

### Functions

def read_last_value(filepath) -> float:
    """Read the last value from a file, ignoring comments and empty lines."""
    
    # Open the file and read lines, filtering out comments and empty lines
    with open(filepath, 'r') as f: 
        lines = [line for line in f.readlines() if not line.startswith('#') and line.strip()] 
    
    # Check if lines are empty
    if not lines: 
        raise ValueError(f"No valid data in {filepath}")
    
    last_line = lines[-1] # Get the last line
    return float(last_line.split()[-1]) # Extract the last value

def read_yplus(filepath) -> float:
    """Read the last Max yPlus value from a file, ignoring comments and empty lines."""

    # Open the file and read lines, filtering out comments and empty lines
    with open(filepath, 'r') as f: # Read the file
        lines = [line for line in f.readlines() if not line.startswith('#') and line.strip()]
    
    # Check if lines are empty
    if not lines: 
        raise ValueError(f"No valid data in {filepath}")
    
    last_line = lines[-1] # Get the last line
    max_yplus_value = float(last_line.split()[-2]) # Extract the second last value (Max yPlus)
    return max_yplus_value

def extract_values(folder_path, folder_name, geometry_type) -> tuple[float, float, float]:
    """Extract values from the specified folder path."""
    try:
        inlet_pressure_file  = glob.glob(os.path.join(folder_path, 'postProcessing/averageTotalPressureInlet1/0/surfaceFieldValue.dat'))[0]
        outlet_pressure_file = glob.glob(os.path.join(folder_path, 'postProcessing/averageTotalPressureOutlet1/0/surfaceFieldValue.dat'))[0]
        massflow_file        = glob.glob(os.path.join(folder_path, 'postProcessing/massFlowRateOutlet1/0/surfaceFieldValue.dat'))[0]
        yplus_file           = glob.glob(os.path.join(folder_path, 'postProcessing/yPlus1/0/yPlus.dat'))[0]
        
        p_inlet   = read_last_value(inlet_pressure_file)
        p_outlet  = read_last_value(outlet_pressure_file)
        massflow  = read_last_value(massflow_file)
        max_yplus = read_yplus(yplus_file)

        delta_p_total = p_inlet - p_outlet

        # Output
        print(f"Geometry: {folder_name} ({geometry_type})")
        print(f"  - Mass flow rate (kg/s)   : {massflow:.4f}")
        print(f"  - Head loss (Pa)          : {delta_p_total:.2f}")
        print(f"  - Max yPlus               : {max_yplus:.2f}")

        if max_yplus > 5.0:
            print("  WARNING: yPlus > 5 detected!")
        if max_yplus >= 1.0:
            print("  WARNING: yPlus not fully below 1!")
        else:
            print("  OK: yPlus all below 1.")
        print("-" * 60)
        return delta_p_total, massflow, max_yplus
    except IndexError:
        print(f"Missing files in {folder}. Skipping...")
    except Exception as e:
        print(f"Error processing {folder}: {e}")

def ensure_dir(path) -> None:
    """Ensure that a directory exists. If not, create it."""
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Directory {path} created.")

### Main script

outputs_dir = 'outputs'
postprocessing_dir = 'postProcessingOutputs'
ensure_dir(postprocessing_dir)

# Check if the outputs directory exists
if not os.path.exists(outputs_dir):
    print(f"Folder {outputs_dir} not found.")
    exit()

# List all folders in the outputs directory
folders = [f for f in os.listdir(outputs_dir) if os.path.isdir(os.path.join(outputs_dir, f))] 

print(f"Found {len(folders)} geometries.\n")

# Store results
results = {}

for folder in folders:
    folder_path = os.path.join(outputs_dir, folder)
    geometry_type = None

    # Parse the folder name
    if folder.startswith('ELL-'):
        geometry_type = 'ELL'
        # Extract parameters from the folder name
        _, Mw, Mb, Kx, Ky, r, L, t = folder.split('-')
        params_key = f"{Kx}_{Ky}_{r}_{L}"
        Mb, Mw = int(Mb), int(Mw)

    elif folder.startswith('IN-'):
        geometry_type = 'IN'
        # Extract parameters from the folder name
        _, Mw, L, t = folder.split('-')
        Mb = 1  # Set Mb to 1 for IN geometry
        params_key = f"{L}"

    elif folder.startswith('RAD-'):
        geometry_type = 'RAD'
        # Extract parameters from the folder name
        _, Mw, Mb, K, L, t = folder.split('-')
        params_key = f"{K}_{L}"

    else:
        print(f"Unknown geometry type in folder {folder}. Skipping.")
        continue
    

    # Extract values from the folder
    result = extract_values(folder_path=folder_path, folder_name=folder, geometry_type=geometry_type)
    if result is None:
        print(f"Skipping folder {folder} due to missing data.")
        continue
    delta_p_total, massflow, max_yplus = result

    
    key = (geometry_type, params_key)

    if key not in results:
        results[key] = {}

    results[key][(Mb, Mw)] = {
        'head_loss':    delta_p_total,
        'massflow':     massflow,
        'max_yplus':    max_yplus
    }

    
# Save the matrices
# Loop over each geometry and parameter combination.
for (geom, param_key), data in results.items():
    mbmw_pairs = list(data.keys()) # Extract all (Mb, Mw) pairs for this configuration.
    Mb_values = sorted(set(mb for mb, mw in mbmw_pairs)) # Get all unique Mb values.
    Mw_values = sorted(set(mw for mb, mw in mbmw_pairs)) # Get all unique Mw values.

    # Create a mapping from Mb and Mw to their row/column indices.
    Mb_index = {mb: i for i, mb in enumerate(Mb_values)}
    Mw_index = {mw: i for i, mw in enumerate(Mw_values)}

    # Initialize matrices with NaN values
    head_loss_matrix    = np.full((len(Mb_values), len(Mw_values)), np.nan)
    massflow_matrix     = np.full((len(Mb_values), len(Mw_values)), np.nan)
    yplus_matrix        = np.full((len(Mb_values), len(Mw_values)), np.nan)

    # Loop through all (Mb, Mw) points.
    for (mb, mw), vals in data.items():
        i = Mb_index[mb]
        j = Mw_index[mw]
        head_loss_matrix[i, j] = vals['head_loss']
        massflow_matrix[i, j] = vals['massflow']
        yplus_matrix[i, j] = vals['max_yplus']

    # Save CSVs
    base_filename = f"{geom}_{param_key}"

    np.savetxt(os.path.join(postprocessing_dir, f"head_losses_{base_filename}.csv"), head_loss_matrix, delimiter=",", fmt="%.6f")
    np.savetxt(os.path.join(postprocessing_dir, f"massflow_{base_filename}.csv"),    massflow_matrix,  delimiter=",", fmt="%.6f")
    np.savetxt(os.path.join(postprocessing_dir, f"max_yplus_{base_filename}.csv"),   yplus_matrix,     delimiter=",", fmt="%.6f")

print("\nPost-processing finished. CSV files saved in postProcessingOutputs/")
