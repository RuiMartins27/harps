import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time

input_template_path = "input/test.in"
output_path = "Outputs/E_amplitude.txt"

# Electron densities to simulate
n_e_values = [2e17, 5e17, 1e18, 2e18, 5e18, 1e19, 2e19]

# Read the template once
with open(input_template_path, "r") as f:
    input_template = f.read()

# Function to update the ElectronDensity in the input string
def update_input_density(template_str, new_density):
    lines = template_str.splitlines()
    new_lines = []
    for line in lines:
        if "ElectronDensity" in line:
            new_line = f"    ElectronDensity           =   {new_density:.2E}          # [m^-3]"
        else:
            new_line = line
        new_lines.append(new_line)
    return "\n".join(new_lines)

# Simulation domain setup (must match input file dimensions)
n_points = 512
z_length = 0.03
x_sim = np.linspace(0, z_length, n_points)

plt.figure(figsize=(10, 6))

# Loop over densities
for n_e in n_e_values:
    # Update input and write it back
    updated_input = update_input_density(input_template, n_e)
    with open(input_template_path, "w") as f:
        f.write(updated_input)

    # Run simulation
    result = subprocess.run('./harps.exe input/test.in', shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Simulation failed for n_e = {n_e:.1e}")
        print(result.stderr)
        continue

    # Wait a moment if needed for output files to be written
    time.sleep(0.5)

    # Load output
    with open(output_path, "r") as f:
        lines = f.readlines()
        sim_values = list(map(float, lines[1].split()[1:]))  # skip first number

    sim_values_norm = np.array(sim_values) / sim_values[0]
    plt.plot(x_sim * 100, sim_values_norm, label=f"Sim $n_e = {n_e:.1e}$")

plt.xlabel('x [cm]')
plt.ylabel('Amplitude (normalized)')
plt.title('Simulated Exponential Decay for Various $n_e$')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
plt.savefig("time_1D", dpi=300)
plt.close()
