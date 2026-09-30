import numpy as np
from numpy import cumsum
import matplotlib.pyplot as plt
import subprocess

# Physical constants
e_charge    = 1.602e-19         # Elementary charge [C]
m_e         = 9.109e-31         # Electron mass [kg]
epsilon_0   = 8.854e-12         # Vacuum permittivity [F/m]
mu_0        = 4*np.pi*1e-7      # Vacuum permeability [H/m]

# Common parameters
nu              = 1e10                  # Collision frequency [Hz]
omega           = 2*np.pi*2.45e9        # Angular frequency [rad/s]
epsilon_r_real  = 1.0                   # Relative real permittivity
z_length        = 0.15                  # Diameter of reactor
n_points        = 1024
n_pml           = 24                    # Number of points in PML

plot_reflector = False

x_sim   = np.linspace(0, z_length, n_points)
x       = np.linspace(0, z_length, n_points)
z_pml   = z_length * (n_points - n_pml) / n_points

density_file_path   = "input/test_non_uniform_ne_elec_dens_values.dat"
locations_file_path = "input/test_non_uniform_ne_elec_dens_points.dat"
reflector_file_path = "input/test_non_uniform_ne_reflector_points.dat"

def create_gaussian_density(n_points, peak_density, width, x_center=None):
    if(x_center==None): x_center = z_length / 2 
    x_vals = np.linspace(0, z_length, n_points)
    gaussian = peak_density * np.exp(-0.5 * ((x_vals - x_center) / width) ** 2)
    return gaussian

def create_MOD_gaussian_density(n_points, peak_density, width, x_center=None):
    if(x_center==None): x_center = z_length / 2 
    x_vals = np.linspace(0, z_length, n_points)
    gaussian = peak_density * np.exp(-0.5 * (abs(x_vals - x_center) / width) ** 3)
    return gaussian

def create_constant_density(n_points, peak_density):
    constant_density = np.full(n_points, peak_density)
    return constant_density

def create_stepped_density(n_points, first_density, second_density, third_density, pos_step_1, pos_step_2):
    step1 = int(pos_step_1*n_points)
    step2 = int(pos_step_2* n_points)
    density = np.zeros(n_points)

    density[:step1] = first_density
    density[step1:step2] =second_density
    density[step2:] = third_density

    return density

# Write the Gaussian density values to the .dat file
def write_density_file(density_file_path, densities):
    with open(density_file_path, "w") as f:
        f.write("{")
        for density in densities:
            f.write(f"({density:.5e})\n")
        f.write("}")

def write_location_file(locations_file_path, n_points):
    with open(locations_file_path, "w") as f:
        f.write("{")
        for i in range(n_points):
            f.write(f"(0,0,{i})\n")
        f.write("}")

def write_reflector_file(locations_file_path, point):
    with open(locations_file_path, "w") as f:
        f.write("{"); f.write(f"(0,0,{point})\n"); f.write("}")

# Set the parameters for the Gaussian density distribution
peak_density = 3e18
width = 0.003

electron_density_values = create_gaussian_density(n_points, peak_density, width)
#electron_density_values = create_MOD_gaussian_density(n_points, peak_density, width, z_length/3)
#electron_density_values = create_constant_density(n_points, 2e18)
#electron_density_values = create_stepped_density(n_points, 1e15, 8e18, 1e15, 0.30, 0.45)
write_density_file(density_file_path, electron_density_values)
write_location_file(locations_file_path, n_points)
write_reflector_file(reflector_file_path, 614)

if plot_reflector:
    average_absorbed_power_reflector = []
    for i in range(n_points):
        if i < n_points*0.1:
            average_absorbed_power_reflector.append(0)
            continue

        write_reflector_file(reflector_file_path, i)

        # Simulation
        result = subprocess.run('mpirun -np 1 ./harps.exe input/test_non_uniform_ne.in', shell=True, capture_output=True, text=True)
        #print(result.stdout)
        if result.returncode != 0:
            print(f"Simulation failed for non-uniform electron density")
            print(result.stderr)
            exit()

        with open("Outputs/E_amplitude.txt", "r") as f:
            lines = f.readlines()
            sim_values = list(map(float, lines[1].split())) 
        sim_values_normalized = (sim_values - np.min(sim_values)) / (np.max(sim_values) - np.min(sim_values))

        with open("Outputs/absorbed_power_density.txt", "r") as f:
            lines = f.readlines()
            sim_values_pabs = list(map(float, lines[1].split())) 

        #print sum of absorbed power density
        pabs_average = np.sum(sim_values_pabs)/n_points
        print(f"Average absorbed power density: {pabs_average:.5e} W/m^3")

        average_absorbed_power_reflector.append(pabs_average)

    plt.plot(x_sim * 100, average_absorbed_power_reflector, label=f"Simulation $P_{{abs}}$")
    plt.xlabel('Reflector Position [cm]')
    plt.ylabel('Average P$_{abs}$ [W/m$^3$]')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("analysis/analysis_output/electron_density_non_uniform_reflector.png", dpi=300)
    plt.close()

# Simulation
result = subprocess.run('mpirun -np 1 ./harps.exe input/test_non_uniform_ne.in', shell=True, capture_output=True, text=True)
#print(result.stdout)
if result.returncode != 0:
    print(f"Simulation failed for non-uniform electron density")
    print(result.stderr)
    exit()

with open("Outputs/E_amplitude.txt", "r") as f:
    lines = f.readlines()

    if len(lines) >= 2 and lines[1].strip() == "UNIFORM_GRID":
        sim_values = list(map(float, lines[2].split()))  # Assuming values are in the third line
    else:
        print("Unexpected output format in E_amplitude.txt. The grid type is not UNIFORM_GRID.")

sim_values_normalized = (sim_values - np.min(sim_values)) / (np.max(sim_values) - np.min(sim_values))
sim_values_normalized = sim_values

with open("Outputs/absorbed_power_density.txt", "r") as f:
    lines = f.readlines()

    if len(lines) >= 2 and lines[1].strip() == "UNIFORM_GRID":
        sim_values_pabs = list(map(float, lines[2].split()))  # Assuming values are in the third line
    else:
        print("Unexpected output format in E_amplitude.txt. The grid type is not UNIFORM_GRID.")

#print sum of absorbed power density
pabs_average = np.sum(sim_values_pabs)/n_points
print(f"Average absorbed power density: {pabs_average:.5e} W/m^3")
sim_values_normalized_pabs = (sim_values_pabs - np.min(sim_values_pabs)) / (np.max(sim_values_pabs) - np.min(sim_values_pabs))
sim_values_normalized_pabs = np.max(sim_values) * np.array(sim_values_pabs) / np.max(sim_values_pabs)


# Create figure with primary axis for wave amplitudes
plt.figure(figsize=(12, 6))
ax1 = plt.gca()

#ax1.plot(x * 100, wave, label="Theory") # there is no analytic solution for non-uniform density
ax1.plot(x_sim * 100, sim_values_normalized, 'k--', linewidth=2, label="Simulation Field")
ax1.plot(x_sim * 100, sim_values_normalized_pabs, color='tab:orange', linewidth=2, label="Simulation $P_{abs}$")
#ax1.axvline(x=z_pml * 100, color='r', linestyle='--', label='PML Boundary')
#ax1.axvline(x=(z_length - z_pml) * 100, color='r', linestyle='--')

# Primary axis labels
ax1.set_xlabel('x [cm]')
ax1.set_ylabel('Amplitude [V/m]', color='black')
ax1.tick_params(axis='y', labelcolor='black')

# Create secondary axis for electron density
ax2 = ax1.twinx()
ax2.plot(x_sim * 100, electron_density_values, 'g-.', linewidth=1.5, label="Electron Density")
ax2.set_ylabel('Electron Density [m$^{-3}$]', color='green')
ax2.tick_params(axis='y', labelcolor='green')
ax2.set_ylim(bottom=0)  # Start from 0 for density

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.title('Wave Amplitudes and Electron Density Profile')
plt.grid(True)
plt.tight_layout()
plt.savefig("analysis/analysis_output/electron_density_non_uniform_test.png", dpi=300)
plt.close()