import numpy as np
from numpy import cumsum
import matplotlib.pyplot as plt
import subprocess

def reflector(x, y, z, i, j, k, points, param):
    center_z = param[0]
    b = param[1]
    y_reflect = param[2]

    if((abs(y - y_reflect) < 4e-4) and z > (center_z - b/2) and z < (center_z + b/2)):
        points.append((int(i), int(j), int(k)))

def write_metal_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set_metal, parameters, name_files):
    """Generate and write permittivity data files based on custom condition function."""
    points = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                function_set_metal(x, y, z, i, j, k, points, parameters)
    
    write_points("input/" + name_files + "_metal_points.dat", points)

def write_points(filename, points):
    """Write points to file."""
    with open(filename, 'w') as f:
        f.write("{")
        for i, (x, y, z) in enumerate(points):
            f.write(f"({int(x)}, {int(y)}, {int(z)})")
        f.write("}")

def load_3d_data(filename):
    with open(filename) as f:
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
        
def generate_grid(L_x, L_y, L_z, N_x, N_y, N_z, uniform_grid=True):
    if uniform_grid:
        x_grid = np.linspace(0, L_x, N_x)
        y_grid = np.linspace(0, L_y, N_y)
        z_grid = np.linspace(0, L_z, N_z)
    else:
        x_grid = np.loadtxt("input/grid_x.in", ndmin=1)
        y_grid = np.loadtxt("input/grid_y.in", ndmin=1)
        z_grid = np.loadtxt("input/grid_z.in", ndmin=1)
    return x_grid, y_grid, z_grid


N_x = 256; N_y = N_x; N_z = 1
L_x = 0.08636; L_y = 0.1400; L_z = 0.04318
uniform_grid = False

x_grid, y_grid, z_grid = generate_grid(L_x, L_y, L_z, N_x, N_y, N_z, uniform_grid)

n_points = 15

absorbed_power_reflector = []
y_reflector_array = []
for i in range(n_points):

    y_reflector = 0.115 + i*(0.13-0.115)/n_points
    y_reflector_array.append(100*y_reflector)
    print(f"Reflector position: {100*y_reflector:.2f} cm")

    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,     reflector, [L_z/2, 1e5*0.04518, y_reflector], "example6")

    # Simulation
    result = subprocess.run('mpirun -np 32 ./harps.exe input/example6.in ', shell=True, capture_output=True, text=True)
    #print(result.stdout)
    
    if result.returncode != 0:
        print(f"Simulation failed for 2D reflector")
        print(result.stderr)
        exit()

    sim_values_pabs, x_grid, y_grid, z_grid = load_3d_data("Outputs/absorbed_power_density.txt")

    voxel_volume = np.zeros((N_x, N_y, N_z))
    for i in range(N_x):
        for j in range(N_y):
            for k in range(N_z):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]

                if N_x == 1: dz = L_x
                else: dx = x_grid[i+1] - x_grid[i] if i < N_x-1 else x_grid[i] - x_grid[i-1]

                if N_y == 1: dy = L_y
                else: dy = y_grid[j+1] - y_grid[j] if j < N_y-1 else y_grid[j] - y_grid[j-1]

                if N_z == 1: dz = L_z
                else: dz = z_grid[k+1] - z_grid[k] if k < N_z-1 else z_grid[k] - z_grid[k-1]

                voxel_volume[i, j, k] = dx * dy * dz

    pabs_total = np.sum(sim_values_pabs * voxel_volume)
    print(f"Absorbed power: {pabs_total:.6e} W/m^3")

    absorbed_power_reflector.append(pabs_total)

plt.plot(y_reflector_array,absorbed_power_reflector, label=f"Simulation $P_{{abs}}$")
plt.xlabel('Reflector Position [cm]')
plt.ylabel('Average P$_{abs}$ [W/m$^3$]')
plt.grid(True)
plt.tight_layout()
plt.savefig("analysis/analysis_output/reflector_2D.png", dpi=300)
plt.close()