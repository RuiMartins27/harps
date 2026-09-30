import numpy as np
import matplotlib.pyplot as plt
import os

harps_dir = os.path.dirname(os.path.abspath(__file__)) + "/../"

def load_3d_data(filename):
    with open(harps_dir+filename) as f:
        # Read dimensions from header
        n_x, n_y, n_z, length_x, length_y, length_z = map(float, f.readline().split())        
        n_x, n_y, n_z = int(n_x), int(n_y), int(n_z)

        header = f.readline().strip()

        if header == "NON_UNIFORM_GRID":
            x_coords = np.array(list(map(float, f.readline().split())))
            y_coords = np.array(list(map(float, f.readline().split())))
            z_coords = np.array(list(map(float, f.readline().split())))
            
            data = np.loadtxt(f)
            data = data.reshape((n_x, n_y, n_z))
            return data, x_coords, y_coords, z_coords
        
        else:
            # Uniform grid mode
            x_coords = np.linspace(0, length_x, n_x)
            y_coords = np.linspace(0, length_y, n_y)
            z_coords = np.linspace(0, length_z, n_z)

            data = np.loadtxt(f)
            data = data.reshape((n_x, n_y, n_z))
            return data, x_coords, y_coords, z_coords

data, grid_x, grid_y, grid_z = load_3d_data("Outputs/absorbed_power_density.txt")
n_x, n_y, n_z = data.shape

length_x = grid_x[-1] - grid_x[0]
length_y = grid_y[-1] - grid_y[0]
length_z = grid_z[-1] - grid_z[0]

omega = 2*np.pi*2.45e9 
mu_0 = 4*np.pi*1e-7
epsilon_0 = 8.854e-12
c_light = 1/np.sqrt(mu_0*epsilon_0)

a = 0.08636
b = 0.04318
R_inner = 0.014

P_abs = 770  # Target absorbed power in W

E_0 = 2.6e-5 # Peak electric field in V/m (wihtout plasma)

center_x = 0 # symmetric
#center_x = a/2

center_y = 0.046 
center_z = b/2


dr = 0.0004  # Radial bin width in meters


if n_z == 1: length_z = b


print_file_radius = True
graph_radius = True
print_file_z = False
graph_z = False
theoretical_quantities = False


# --- Define radial bins ---
X, Y = np.meshgrid(grid_x, grid_y, indexing='ij')
R = np.sqrt((X - center_x)*(X - center_x) + (Y - center_y)*(Y - center_y))

# Flatten the 2D slices
r_vals = R.flatten()
r_bins = np.arange(0, R_inner + dr, dr)
r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])

# get the voxel volume and absorbed power density
voxel_volume = np.zeros((n_x, n_y, n_z))
pabs_dens = np.zeros((n_x, n_y, n_z))
for i in range(n_x):
    for j in range(n_y):
        for k in range(n_z):
            if n_x == 1: dz = length_x
            elif (i < n_x-1): dx = grid_x[i+1] - grid_x[i]
            else: grid_x[i] - grid_x[i-1]

            if n_y == 1: dy = length_y
            elif (j < n_y-1): dy = grid_y[j+1] - grid_y[j]
            else: grid_y[j] - grid_y[j-1]

            if (n_z == 1): dz = length_z
            elif(k < n_z-1): dz = grid_z[k+1] - grid_z[k]
            else: grid_z[k] - grid_z[k-1]

            voxel_volume[i, j, k] = dx * dy * dz
            pabs_dens[i, j, k] = data[i, j, k]

# --- 1. Radial power density profile over all z ---
radial_power = np.zeros_like(r_centers)
radial_volume = np.zeros_like(r_centers)

for k in range(n_z):
    for i in range(n_x):
        for j in range(n_y):
            r = R[i, j]
            bin_idx = np.searchsorted(r_bins, r) - 1
            if 0 <= bin_idx < len(radial_power):
                radial_power[bin_idx] += pabs_dens[i, j, k] * voxel_volume[i, j, k]
                radial_volume[bin_idx] += voxel_volume[i, j, k]

# Normalize by volume to get average absorbed power density
radial_power_density = np.divide(
    radial_power,
    radial_volume,
    out=np.zeros_like(radial_power),
    where=radial_volume > 0
)

# --- Save and print to file ---
if(print_file_radius):
    with open(harps_dir + "analysis/analysis_output/p_abs_radius.txt", "w") as f:
        f.write("Radius (m) | Absorbed Power Density (W/m³)\n")
        for r, p in zip(r_centers, radial_power_density):
            f.write(f"{r:.4f}     {p:.6e}\n")

if(graph_radius):
    plt.figure()
    plt.plot([rad*100 for rad in r_centers], radial_power_density)
    plt.xlabel("Radius [cm]", fontsize=16)
    plt.ylabel("$P_{dens}$ [W/m³]", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #plt.title("Radial Profile Averaged Over z")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(harps_dir + "analysis/analysis_output/p_abs_radius.png", dpi=300)
    plt.close()

# --- 2. Power density profile over the z slices ---
axial_power = np.zeros(n_z)
axial_volume = np.zeros(n_z)

for k in range(n_z):
    for i in range(n_x):
        for j in range(n_y):
            r = R[i, j]

            if(r < R_inner):
                axial_power[k]  += pabs_dens[i, j, k] * voxel_volume[i, j, k]
                axial_volume[k] += voxel_volume[i, j, k]

# Normalize to get average absorbed power density in each z-slice
axial_power_density = np.divide(
    axial_power,
    axial_volume,
    out=np.zeros_like(axial_power),
    where=axial_volume > 0
)

z_centers = grid_z

if(print_file_z):
    print("Z Position (m) | Absorbed Power Density (W/m³)")
    for z, p in zip(z_centers, axial_power_density):
        print(f"{z:.4f}        {p:.6e}")

if(graph_z):
    plt.figure()
    plt.plot(z_centers, axial_power_density)
    plt.xlabel("Axial Position z (m)", fontsize=16)
    plt.ylabel("Average Absorbed Power Density (W/m³)", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title("Axial Profile Averaged Over Radius")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(harps_dir + "analysis/analysis_output/p_abs_z.png", dpi=300)
    plt.close()


# --- Theoretical Quantities ---
if(theoretical_quantities):
    pabs_total = np.sum(pabs_dens * voxel_volume)
    print(f"Absorbed power: {pabs_total:.6e} W/m^3")

    total_volume = np.sum(voxel_volume)
    print(f"Total volume: {total_volume:.6e} m³")

    average_absorbed_power_density = pabs_total / total_volume
    print(f"Average absorbed power density: {average_absorbed_power_density:.6e} W/m³")

    omega_c = c_light * np.pi / a
    impedance = np.sqrt(mu_0 / epsilon_0)
    impedance_TE = impedance / np.sqrt(1 - (omega_c / omega)**2)
    input_power = a * b * E_0**2 / (4 * impedance_TE)

    print(f"\nTE10 mode impedance: {impedance_TE:.6e} Ohm")
    print(f"Theoretical input power through the waveguide: {input_power:.6e} W")
    print(f"To scale this to {P_abs:.2f} W, multiply the excitation wave by {(P_abs / pabs_total):.6e}")