# canonical_orbital_area_occupancy

This script estimates the occupancy of canonical orbitals in a specified region by determining how much each orbital is located inside or outside that region.

## Requirements

Ensure that the following files and folder are present in your working directory:
- `input.yaml`
- `posinp.xyz`
- `localization.input`
- Folder: `data`

Also, make sure to set the following parameters in your configuration:
- `nvirt: 0`
- `norbv` to your desired number

## Input File: localization.input

An example `localization.input` file for 109_2ASP is provided in the working directory. In this file, the **Index and Radii** section contains several lines, each with two numbers:
- The **first number** is the index of the atom that serves as the center of the specified region.
- The **second number** is the radius around that center.

> **Note:** The atom indices in the `posinp.xyz` file start at 1.

## Usage

Simply compile (make) and run the script without any additional commands.

## Output

After running the script, two output files will be generated:
- `area_occupancy_occ.out`
- `area_occupancy_virt.out`

Each line in these files corresponds to an orbital index. The value (ranging from **-1** to **+1**) on each line indicates how much the corresponding orbital is outside the specified region:
- A value of **-1** means the orbital is completely inside the region.
- A value of **+1** means the orbital is completely outside the region.

## Summary

- **Input Files/Folders:** `input.yaml`, `posinp.xyz`, `localization.input`, `data/`
- **Configuration Parameters:** Set `nvirt: 0` and `norbv` as required.
- **Usage:** Make and run the script.
- **Output Files:** `area_occupancy_occ.out` and `area_occupancy_virt.out` with orbital occupancy values.

Happy computing!
