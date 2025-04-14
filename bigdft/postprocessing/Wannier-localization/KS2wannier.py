import subprocess
import yaml
import sys
import os

wann_executable = "/data/wannier90-3.1.0/wannier90.x"
b2w_executable = "/data/bigdft/bigdft-suite/build/bigdft/src/BigDFT2Wannier"

def execute_command_with_realtime_output(command):
    """
    Execute a terminal command and show its output in real-time.

    Args:
        command (list or str): The command to execute.
    """
    # Start the process
    process = subprocess.Popen(
        command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    # Read and print output line by line
    while True:
        output = process.stdout.readline()
        if output == "" and process.poll() is not None:  # Check if the process is done
            break
        if output:
            print(output.strip())

    # Capture and print any remaining error messages
    error_output = process.stderr.read()
    if error_output:
        print("Error:", error_output.strip())

        # Terminate the subprocess
        process.terminate()
        process.wait()

        # Exit the script with an error code
        sys.exit(1)

    # Ensure the process has completed
    process.wait()

# Load settings from `b2w.yaml`
def load_b2w_yaml():
    yaml_file_path = "b2w.yaml"
    if os.path.exists(yaml_file_path):
        with open(yaml_file_path, "r") as file:
            settings = yaml.safe_load(file)
        return settings if settings else {}
    else:
        raise FileNotFoundError(f"YAML file '{yaml_file_path}' not found.")

# Load settings from `b2w.yaml`
settings = load_b2w_yaml()

# Extract settings with default values if not provided in `b2w.yaml`
n_wann = settings.get("n_wann", None)  # Required field
if not n_wann:
    raise ValueError("The 'n_wann' value is missing or empty in 'b2w.yaml'.")

n_iter = settings.get("n_iter", 300)
n_occ = settings.get("n_occ", "all")
n_virt_tot = settings.get("n_virt_tot", 0)
n_virt = settings.get("n_virt", 0)
name = settings.get("name", "default_sys")
xyz = settings.get("xyz", "posinp")
input_yaml_name = settings.get("input_yaml_name", "input")
write_UNK = settings.get("write_UNK", True)
additional_settings = settings.get("additional_settings", "")

if write_UNK:
   write_UNK = "T"
else:
   write_UNK = "F"

# Extract `projections` as a long string from YAML
projections = settings.get("projections", "")
# if not projections:
#     print("The 'projections' value is missing or empty in 'b2w.yaml'.\n \
#                       Using random projections.")

with open ("./data/wavefunction-k001-NR.b000001", "r") as f:
  f.readline() #skip the first line
  grid_lengths = f.readline().split()
  grid_nums = f.readline().split()

unit_cell_carts = [f"{float(grid_lengths[i])*float(grid_nums[i]):.6f}" for i in range(3)]

if isinstance(n_occ, str) and n_occ.lower()=="all":
  with open("log.yaml", "r") as f:
      log_data = yaml.safe_load(f)

  if "Total Number of Orbitals" in log_data:
      n_occ = log_data["Total Number of Orbitals"]
  else:
      raise KeyError("Key 'Total Number of Orbitals' not found in log.yaml")

num_bands = n_occ+n_virt
slwf_num = n_wann

if num_bands<n_wann:
   raise ValueError("Not enough orbs to build wannier functions. (nocc+nvirt<nw)")

win1 = f""" postproc_setup = .true.
 write_u_matrices = true
 num_wann        =  {num_bands}
 num_iter        = {n_iter}
 slwf_num       = {n_wann}
 write_xyz = true

 begin atoms_cart
 ang
"""

win2 = f""" end atoms_cart

 begin unit_cell_cart
 bohr
    {unit_cell_carts[0]}      0.000000      0.000000
     0.000000     {unit_cell_carts[1]}      0.000000
     0.000000      0.000000     {unit_cell_carts[2]}
 end unit_cell_cart
 
 begin projections
random
"""

win3 = """
 end projections

 mp_grid    : 1 1 1
 gamma_only : true

 begin kpoints
 0.0 0.0 0.0
 end kpoints"""

inter = f"""{name}                  # Name of the .win file
form                  # Format : cube or etsf
F  F   {n_occ}             #  Use resolution of the identity (ROI), write ROI states, No. of occupied orbitals
F    {n_virt_tot}    {n_virt}          # Pre-check, n_virt_tot, n_virt
{write_UNK}    F    F    F     # Write_UNKp.s, write_spherical_harmonics, write_angular_parts,  write_radial_parts
  {int(grid_nums[0])*2}  {int(grid_nums[1])*2}  {int(grid_nums[2])*2}       # Number of points for each axis in the cubic  BigDFT representation (information needed by Wannier90)
 data"""

# Modify the YAML file------------------------------------ Begin
with open(f"{input_yaml_name}.yaml", "r") as file:
    data = yaml.safe_load(file)

data['dft']['norbv'] = int(n_virt)
data['dft']['nvirt'] = 0

with open(f"{input_yaml_name}.yaml", "w") as file:
    yaml.dump(data, file, default_flow_style=False)
# Modify the YAML file------------------------------------ End

with open(f"./{xyz}.xyz", "r") as xyzfile:
  n_atom = int(xyzfile.readline().split()[0])
  xyzfile.readline()# skip the second line
  
  atoms = xyzfile.readlines()

with open(f"{name}.win", "w") as winfile:
  winfile.write(win1)

  atom_count = 0
  for atom in atoms:
    winfile.write(atom)
    atom_count += 1

  if n_atom!=atom_count:
    raise ValueError(f"The number of atoms doesn't match the header.(In ./{xyz}.xyz)")
  
  winfile.write(win2)

  projections_lines = projections.splitlines()
  for projection in projections_lines:
      if projection.strip() and projection[0].isdigit():  # Check if the line starts with a digit
          atom_ind = int(projection.split(":")[0].strip())  # Extract the atom index
          _, x, y, z = atoms[atom_ind - 1].split()  # Get coordinates from `atoms`
          proj_to_write = f"c={x},{y},{z}:{projection.split(':', 1)[1].strip()}"  # Format projection
          winfile.write(proj_to_write + '\n')  # Write to winfile with newline
      else:
          winfile.write(projection + '\n')  # Write the line as-is if it doesn't start with a digit

  
  winfile.write(win3)

  if additional_settings:
     winfile.write("\n")
     winfile.write("\n")
     winfile.write(additional_settings)


execute_command_with_realtime_output(f"{wann_executable} {name}")


with open(f"{name}.win", "r") as winfile:
  lines = winfile.readlines()
  lines[0] = " postproc_setup = .false.\n"

with open(f"{name}.win", "w") as winfile:
  winfile.writelines(lines)

with open("input.inter", "w") as interfile:
   interfile.write(inter)

execute_command_with_realtime_output(b2w_executable)

execute_command_with_realtime_output(f"{wann_executable} {name}")
