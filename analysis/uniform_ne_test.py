import numpy as np
import matplotlib.pyplot as plt
import subprocess
import time

# Physical constants
e_charge = 1.602e-19        # Elementary charge [C]
m_e = 9.109e-31             # Electron mass [kg]
epsilon_0 = 8.854e-12       # Vacuum permittivity [F/m]
mu_0 = 4 * np.pi * 1e-7     # Vacuum permeability [H/m]

# Common parameters
nu = 1e10                   # Collision frequency [Hz]
omega = 2 * np.pi * 2.45e9  # Angular frequency [rad/s]
epsilon_r_real = 1.0        # Relative real permittivity
z_length = 0.03             # Diameter of reactor
n_points = 512
n_pml = 24
x_sim = np.linspace(0, z_length, n_points)
x = np.linspace(0, z_length, n_points)
z_pml = z_length*(n_points - n_pml) / n_points

input_template_path = "input/test_uniform_ne.in"
output_path = "Outputs/E_amplitude.txt"

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

# Electron densities to compare
n_e_values = [2e17, 5e17, 1e18, 2e18, 5e18, 1e19, 2e19]  # m⁻³

plt.figure(figsize=(10, 6))

with open(input_template_path, "r") as f:
    input_template = f.read()

# To store relative error data
relative_errors = []


for n_e in n_e_values:
    updated_input = update_input_density(input_template, n_e)
    with open(input_template_path, "w") as f:
        f.write(updated_input)

    # Run simulation
    result = subprocess.run('mpirun -np 4 ./harps.exe input/test_uniform_ne.in', shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"Simulation failed for n_e = {n_e:.1e}")
        print(result.stderr)
        continue

    sigma_complex = (e_charge**2 * n_e) / (m_e * (nu + 1j * omega))
    epsilon_r_complex = epsilon_r_real - 1j * sigma_complex / (omega * epsilon_0)
    k_complex = omega * np.sqrt(mu_0 * epsilon_0 * epsilon_r_complex)
    k_imag = np.imag(k_complex)
    wave = np.exp(k_imag * x)

    plt.plot(x * 100, wave, label=f"Theory $n_e = {n_e:.1e}$ m⁻³")

    with open("Outputs/E_amplitude.txt", "r") as f:
        lines = f.readlines()
        if len(lines) >= 2 and lines[1].strip() == "UNIFORM_GRID":
            sim_values = list(map(float, lines[2].split()))  # Assuming values are in the third line
        else:
            print("Unexpected output format in E_amplitude.txt. The grid type is not UNIFORM_GRID.")
            
    plt.plot(x_sim * 100, sim_values, 'k--', linewidth=2, label=f"Simulation $n_e = {n_e:.1e}$ m⁻³")

    # Calculate relative error: set to 0 if x_sim >= z_pml
    rel_error = np.zeros_like(sim_values)
    for i, x_val in enumerate(x_sim):
        if x_val < z_pml:
            rel_error[i] = np.where(np.abs(wave[i]) > 1e-10, np.abs(sim_values[i] - abs(wave[i])) / np.abs(wave[i]), 0)

    relative_errors.append(rel_error)

# Add the vertical line to the plot
plt.axvline(x=z_pml*100, color='r', linestyle='--', label='PML Boundary')

plt.xlabel('x [cm]')
plt.ylabel('Amplitude (normalized)')
plt.title('Exponential Decay: Theoretical vs Simulation')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("analysis/analysis_output/electron_density_test.png", dpi=300)
plt.close()

# Plot relative error
plt.figure(figsize=(10, 6))

for idx, n_e in enumerate(n_e_values):
    plt.plot(x_sim * 100, relative_errors[idx], label=f"$n_e = {n_e:.1e}$")

plt.xlabel('x [cm]')
plt.ylabel('Relative Error')
plt.title('Relative Error between Simulation and Theory')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("analysis/analysis_output/relative_error.png", dpi=300)
plt.close()
