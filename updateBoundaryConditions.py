from pathlib import Path
# import re

# Define the path to your OpenFOAM "boundary" file
BOUNDARY_FILE = Path("constant/polyMesh/boundary")

# Map of patch names to desired boundary types
boundary_updates = {
    "inlet": {
        "type": "patch",
        "physicalType": "inlet"
    },
    "outlet": {
        "type": "patch",
        "physicalType": "outlet"
    },
    "wall": {
        "type": "wall",
        "physicalType": "wall"
    },
    "symmetry": {
        "type": "symmetryPlane",
        "physicalType": "symmetry"
    },
    "frontAndBack": {
        "type": "empty",
        "physicalType": "empty"
    }
}

def update_boundary_file(path):
    with open(path, "r") as file:
        lines = file.readlines()

    new_lines = []
    in_patch = False
    current_patch = ""
    
    for line in lines:
        stripped = line.strip()
        
        # Detect start of patch
        if stripped in boundary_updates:
            in_patch = True
            current_patch = stripped
            new_lines.append(line)
            continue

        # If inside a patch block, look for lines to update
        if in_patch and "type" in stripped and current_patch in boundary_updates:
            new_type = boundary_updates[current_patch]["type"]
            new_lines.append(f"        type            {new_type};\n")
            continue

        if in_patch and "physicalType" in stripped and current_patch in boundary_updates:
            new_ptype = boundary_updates[current_patch]["physicalType"]
            new_lines.append(f"        physicalType    {new_ptype};\n")
            in_patch = False  # Done updating this patch
            continue

        new_lines.append(line)

    # Write the updated file back
    with open(path, "w") as file:
        file.writelines(new_lines)
    print(f"[INFO] Updated boundary file at: {path}")

if __name__ == "__main__":
    update_boundary_file(BOUNDARY_FILE)
