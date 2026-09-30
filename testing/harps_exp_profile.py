import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

# Define the datasets from your output tables
datasets = {
    "CO2_100mbar_Steeg": pd.DataFrame({
        "r_m": [0.0000, 0.0005, 0.0010, 0.0015, 0.0020, 0.0030, 0.0100],
        "ne_m3": [2.388643e+18, 2.955209e+18, 2.631286e+18, 2.435312e+18, 2.342868e+18, 2.086064e+18, 1.000000e+11],
        "mu_real_m2_Vs": [5.354722, 5.394787, 5.403212, 5.374862, 5.367857, 5.064549, 0.583787],
        "mu_imag_m2_Vs": [-5.784804, -5.399507, -5.528224, -5.451809, -5.297783, -3.398110, -0.033578]
    }),
    "CO2_250mbar_Steeg": pd.DataFrame({
        "r_m": [0.0000, 0.0005, 0.0010, 0.0015, 0.0020, 0.0030, 0.0100],
        "ne_m3": [3.383855e+19, 2.735110e+19, 8.574091e+18, 9.436043e+18, 3.800423e+18, 1.167419e+18, 1.000000e+11],
        "mu_real_m2_Vs": [5.052097, 4.761586, 4.919960, 4.312253, 4.332875, 2.997263, 0.235582],
        "mu_imag_m2_Vs": [-3.960007, -3.356084, -3.686174, -2.318860, -2.400659, -0.967332, -0.005414]
    }),
    "Ar_20mbar_Huber": pd.DataFrame({
        "r_m": [0.000157, 0.000351, 0.000843, 0.001351, 0.001851, 0.010000],
        "ne_m3": [4.765957e+19, 4.714894e+19, 4.340426e+19, 3.574468e+19, 2.331916e+19, 1.000000e+11],
        "mu_real_m2_Vs": [3.726352, 3.732205, 3.842609, 4.346500, 4.712823, 3.663855],
        "mu_imag_m2_Vs": [-8.983942, -8.967484, -8.778149, -7.688560, -6.165136, -1.907750]
    }),
    "Ar_88mbar_Huber": pd.DataFrame({
        "r_m": [0.000090, 0.000403, 0.000597, 0.001090, 0.001910, 0.010000],
        "ne_m3": [8.817021e+19, 7.914893e+19, 7.302127e+19, 5.174468e+19, 1.429787e+19, 1.000000e+11],
        "mu_real_m2_Vs": [4.360392, 4.359367, 4.152418, 2.728806, 1.920415, 1.074177],
        "mu_imag_m2_Vs": [-3.582568, -3.564287, -3.196182, -1.068092, -0.512936, -0.139166]
    }),
    "Ar_1000mbar_Khazem": pd.DataFrame({
        "r_m": [0.000055, 0.000203, 0.000399, 0.000605, 0.000801, 0.001000, 0.010000],
        "ne_m3": [3.069622e+20, 2.941945e+20, 1.529332e+20, 1.226199e+19, 1.226199e+17, 1.226199e+13, 1.000000e+11],
        "mu_real_m2_Vs": [0.422754, 0.414578, 0.367255, 0.300680, 0.217180, 0.167316, 0.096661],
        "mu_imag_m2_Vs": [-0.022099, -0.021371, -0.017194, -0.011199, -0.005505, -0.003371, -0.001129]
    })
}


import numpy as np
import pandas as pd
import os
import re
import math
import matplotlib.pyplot as plt
from scipy import stats
import subprocess
from scipy.interpolate import interp1d

HARPS_DIR = "./"
f = 2.45e9                  # frequency in Hz
k_B = 1.380649e-23          # Boltzmann constant in J/K
e_charge = 1.602176634e-19  # elementary charge in C
MU_0 = 4e-7 * np.pi        # vacuum permeability in H/m
EPSILON_0 = 8.854187817e-12 # vacuum permittivity in F/m
m_e = 9.10938356e-31        # electron mass in kg
omega = 2*np.pi*f           # angular frequency in rad/s


E_CHARGE = 1.602176634e-19
C_LIGHT = 1/np.sqrt(MU_0*EPSILON_0)

P_in = 770                   # input power in W
E_0 = 21000

flag_center_coordinates = False

def write_values(filename, values):
    """Write float values to file."""
    with open(filename, 'w') as f:
        f.write("{")
        for i, value in enumerate(values):
            f.write(f"({value:.4e})")
            if i < len(values) - 1:
                f.write(" \n")
        f.write("}")

def write_points(filename, points):
    """Write points to file."""
    with open(filename, 'w') as f:
        f.write("{")
        for i, (x, y, z) in enumerate(points):
            f.write(f"({int(x)}, {int(y)}, {int(z)})\n")
        f.write("}")

def interpolated_radial_profile(x, y, z, i, j, k, points, values, param):
    interp_func = param[0]
    center_x    = param[1]
    center_y    = param[2]
    r_max       = param[3]

    r = np.sqrt((x - center_x)**2 + (y - center_y)**2)

    if r > r_max:
        return

    val = interp_func(r)

    points.append((int(i), int(j), int(k)))
    values.append(float(val))

def write_real_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set, parameters, name_files, name_var):
    """Generate and write real-valued data files based on custom condition function."""
    points = []
    values = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                function_set(x, y, z, i, j, k, points, values, parameters)
    
    write_values(HARPS_DIR + "/input/" + name_files + "_" + name_var + "_values.dat", values)
    write_points(HARPS_DIR + "/input/" + name_files + "_" + name_var + "_points.dat", points)

def generate_grid(n_x, n_y, n_z, length_x, length_y, length_z, refinementFactor_x=-1, refinementFactor_y=-1, refinementFactor_z=-1,
    x_BL=0, x_RL=0, x_RR=0, x_BR=0, y_BL=0, y_RL=0, y_RR=0, y_BR=0, z_BL=0, z_RL=0, z_RR=0, z_BR=0, x_grid=None, y_grid=None, z_grid=None, x_file=None, y_file=None, z_file=None):

    def process_direction(n, length, refinementFactor, BL, BR, RL, RR, grid, file):
        if grid is not None and len(grid) > 0:
            arr = np.array(grid)
            if not np.all(np.diff(arr) >= 0):
                raise ValueError("Input grid is not monotonically increasing.")
            d = np.diff(arr)
            return arr, d, True
        
        elif file is not None and os.path.exists(file):
            arr = np.loadtxt(file, ndmin=1)
            if not np.all(np.diff(arr) >= 0):
                raise ValueError("Input file grid is not monotonically increasing.")
            d = np.diff(arr)
            return arr, d, True
        
        elif refinementFactor != -1:
            if n < 2:
                raise ValueError("n must be at least 2 for refinement.")
            if (abs(BL) < 1e-6) and (abs(RL) < 1e-6) and (BR > RR):
                symmetric = True
            elif (BL > RL or RL > RR or RR > BR) and not(abs(BL) < 1e-6) and (abs(RL) < 1e-6) and (BR > RR):
                raise ValueError("Key points must be monotonically increasing: BL <= RL <= RR <= BR. It can also be only one sided (0, 0, RR, BR)")
            
            arr = np.zeros(n)
            d = np.zeros(n-1)

            normal = (BL+length - BR+(RR-RL)*refinementFactor+refinementFactor/(refinementFactor-1)*(RL-BL+BR-RR)*np.log(refinementFactor))/(n-1)
            refined = normal / refinementFactor

            arr[0] = 0.0
            for i in range(1, n):
                prev = arr[i-1]
                if prev < BL:
                    arr[i] = prev + normal
                elif prev < RL:
                    frac = (prev - BL) / (RL - BL)
                    arr[i] = prev + normal*(1 - frac) + refined*frac
                elif prev < RR:
                    arr[i] = prev + refined
                elif prev < BR:
                    frac = (prev - RR) / (BR - RR)
                    arr[i] = prev + refined*(1 - frac) + normal*frac
                else:
                    arr[i] = prev + normal

            alpha = length/arr[n-1]
            for i in range(1, n):
                if(arr[i]>RR): arr[i] = alpha*arr[i]

            d = np.diff(arr)

            return arr, d, True
        else:
            # Uniform grid
            arr = np.linspace(0, length, n)
            d = np.diff(arr)

            return arr, d, False

    x, dx, non_uniform_x = process_direction(n_x, length_x, refinementFactor_x, x_BL, x_BR, x_RL, x_RR, x_grid, x_file)
    y, dy, non_uniform_y = process_direction(n_y, length_y, refinementFactor_y, y_BL, y_BR, y_RL, y_RR, y_grid, y_file)
    z, dz, non_uniform_z = process_direction(n_z, length_z, refinementFactor_z, z_BL, z_BR, z_RL, z_RR, z_grid, z_file)

    non_uniform_grid = non_uniform_x or non_uniform_y or non_uniform_z
    return x, dx, y, dy, z, dz, non_uniform_grid

def load_3d_data(filename):
    try:
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

                # Selecting a line from the data
                #data = data[n_x // 2, :, n_z // 2]; data = data[np.newaxis, :, np.newaxis]

                if(flag_center_coordinates == True):
                    x_coords = x_coords - center_x
                    y_coords = y_coords - center_y
                    z_coords = z_coords - center_z

                return data, x_coords, y_coords, z_coords
            
            else:
                # Uniform grid mode
                x_coords = np.linspace(0, length_x, n_x)
                y_coords = np.linspace(0, length_y, n_y)
                z_coords = np.linspace(0, length_z, n_z)

                data = np.loadtxt(f)
                data = data.reshape((n_x, n_y, n_z))

                if(flag_center_coordinates == True):
                    x_coords = x_coords - center_x
                    y_coords = y_coords - center_y
                    z_coords = z_coords - center_z

                return data, x_coords, y_coords, z_coords
        
    except FileNotFoundError:
        print(f"Warning: File '{filename}' not found. Skipping.")
        return None, None, None, None

def get_radial_profile(data, grid_x, grid_y, grid_z, dr = 0.0003, plot=False, title="power_density"):
    if data.all() == None:
        return None, None
    n_x, n_y, n_z = data.shape
    
    # Precompute Δ arrays
    dx = np.diff(grid_x)
    dy = np.diff(grid_y)
    dz = np.diff(grid_z)

    if(len(dz)<1): dz = [L_z]
    
    # Radial bins
    r_bins = np.arange(0, R_in + dr, dr)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    radial_quantity = np.zeros_like(r_centers)
    radial_volume = np.zeros_like(r_centers)

    # Compute distance grid in XY plane
    X, Y = np.meshgrid(grid_x, grid_y, indexing='ij')
    R = np.sqrt((X - center_x)**2 + (Y - center_y)**2)
    
    # Loop over all voxels
    for i in range(n_x):
        for j in range(n_y):
            for k in range(n_z):
                # Compute local voxel volume (same as C++)
                vol = (cell_size(dx, i, n_x) * cell_size(dy, j, n_y) * cell_size(dz, k, n_z))
                r = R[i, j]
                bin_idx = np.searchsorted(r_bins, r) - 1
                if 0 <= bin_idx < len(r_centers):
                    radial_quantity[bin_idx] += data[i, j, k] * vol
                    radial_volume[bin_idx] += vol

    # Compute radial average
    with np.errstate(divide='ignore', invalid='ignore'):
        radial_quantity_density = np.divide(radial_quantity, radial_volume, out=np.zeros_like(radial_quantity), where=radial_volume > 0)

    return r_centers, radial_quantity_density

def create_1d_profiles(r_centers, E_field_1D ,ne_0_exp, sig2_exp, p_exp, idx_cond, Tg_0_exp, P_in_exp):
    if r_centers is None:
        r_centers = np.linspace(0, 0.014, 20)
        E_reduced = np.full_like(r_centers, 100)

        n_e = ne_0_exp*np.exp(-r_centers**2/(2*(0.25*sig2_exp*sig2_exp)))

        T_g = 700+(Tg_0_exp-700)*np.exp(-r_centers**2/(2*(sig2_exp*sig2_exp)))

        N = p_exp/(k_B*T_g)
    else:
        n_e = ne_0_exp*np.exp(-r_centers**2/(2*(0.25*sig2_exp*sig2_exp)))

        T_g = 700+(Tg_0_exp-700)*np.exp(-r_centers**2/(2*(sig2_exp*sig2_exp)))

        N = p_exp/(k_B*T_g)
        E_reduced = (P_in_exp/P_in)*E_field_1D*1e21/N/np.sqrt(2) # Td



    # depending which gas is used
    match idx_cond:
        case 0: # N2
            mu_interp = get_mu_interpolator(HARPS_DIR + "testing/mu_n2.csv")
        case 1: # CO2
            mu_interp = get_mu_interpolator(HARPS_DIR + "testing/mu_co2.csv")

    mu_DC = mu_interp(E_reduced)/N
    mu_hat = mu_DC/(1 - (1j*omega*mu_DC*m_e)/e_charge)

    #print(T_g)
    #print(N)
    #print(E_reduced)
    #print(mu_DC)
    #print(mu_hat)

    n_e_interp = interp1d(r_centers, n_e, kind="linear", bounds_error=False, fill_value=0 )
    mu_re_interp = interp1d(r_centers, np.real(mu_hat), kind="linear", bounds_error=False, fill_value=0)
    mu_im_interp = interp1d(r_centers, -np.imag(mu_hat), kind="linear", bounds_error=False, fill_value=0)

    return n_e_interp, mu_re_interp, mu_im_interp


def get_mu_interpolator(filename):
    # skipinitialspace=True handles the "comma + space" issue automatically
    df = pd.read_csv(filename, skipinitialspace=True)
    
    # Just in case there are other weird spaces, strip them all
    df.columns = df.columns.str.strip()
    
    # SORTING IS CRITICAL: Your new data is in descending order. 
    # interp1d requires the x-axis (EoN_Td) to be increasing.
    df = df.sort_values(by='EoN_Td')
    
    return interp1d(df['EoN_Td'], df['muN'], kind='linear', fill_value="extrapolate")

def cell_size(d, idx, n):
        if n == 1:
            return d[0]
        if idx == 0:
            return 0.5 * d[0]
        if idx == n - 1:
            return 0.5 * d[-1]
        return 0.5 * (d[idx] + d[idx - 1])

def harps_exp_profile(ne_0_exp, sig2_exp, p_exp, idx_cond, Tg_0_exp, P_in_exp):
        data, x_coords, y_coords, z_coords = load_3d_data(HARPS_DIR + "Outputs/E_amplitude.txt")
        r_centers, E_field_1D = get_radial_profile(data, x_coords, y_coords, z_coords)
        ne_interp, mu_r_interp, mu_i_interp = create_1d_profiles(r_centers, E_field_1D ,ne_0_exp, sig2_exp, p_exp, idx_cond, Tg_0_exp, P_in_exp)
        write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [ne_interp, center_x, center_y, 0.014], "2D_run", "elec_dens")
        write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [mu_r_interp, center_x, center_y, 0.014], "2D_run", "real_mobi")
        write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [mu_i_interp, center_x, center_y, 0.014], "2D_run", "imag_mobi")

def calculate_skin_depth(ne, mu_real, mu_imag):
    complex_conductivity = e_charge * ne * (mu_real + 1j * mu_imag)
    #complex_permittivity = real_permittivity - 1j * complex_conductivity / (angular_frequency * EPSILON_0)

    sigma_real = np.real(complex_conductivity)
    sigma_imag = np.imag(complex_conductivity)

    # Components of k² = ω² μ0 ε0 - i ω μ0 σ
    k2_re = omega*omega * MU_0 * EPSILON_0 - omega * MU_0 * sigma_imag
    k2_im = - omega * MU_0 * sigma_real

    # Compute Im(k)
    k_imag = np.sqrt((np.sqrt(k2_re*k2_re + k2_im*k2_im) - k2_re) / 2.0)
    k_imag = np.where(k2_im < 0, -k_imag, k_imag)

    # Skin depth δ = 1 / Im(k)
    skin_depth = np.where(k_imag != 0, 1.0 / np.abs(k_imag), np.inf)

    return skin_depth

def calibrate_empty_waveguide_E(R_in, y_r):
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, reflector, [center_z, 1e5*b, y_r], "2D_run")
    
    #result = subprocess.run('mpirun --bind-to socket -np 1 ./harps.exe input/2D_run_no_plasma.in', shell=True, capture_output=True, text=True)
    result = subprocess.run('./harps.exe input/2D_run_no_plasma.in', shell=True, capture_output=True, text=True)
    
    data_E_amplitude, grid_x, grid_y, grid_z = load_3d_data("Outputs/E_amplitude.txt")

    radii_samples = [0, R_in / 4, R_in / 2, R_in, 2 * R_in]
    integrated_values = [0]

    for r_sample in radii_samples[1:]:
        _, integrated_E_field = calculate_first_moment( data_E_amplitude, grid_x, grid_y, grid_z, r_max=r_sample)
        integrated_values.append(integrated_E_field)

    return interp1d(radii_samples, integrated_values, kind="cubic", fill_value="extrapolate")

def calculate_radial_average(r_centers, radial_data, r_plasma=0.014):
    dr = r_centers[1] - r_centers[0]
    r_max = r_centers[-1] + dr/2

    radial_total = 0
    total_area = 0

    for i in range(len(r_centers)):
        r = r_centers[i]
        if(r > r_plasma): continue

        area = 2 * r * dr
        radial_total  += radial_data[i] * area
        total_area += area

    radial_average = radial_total/total_area

    return radial_average

def calculate_first_moment(data, grid_x, grid_y, grid_z, r_max):
    n_x, n_y, n_z = data.shape
    
    # Precompute Δ arrays
    dx = np.diff(grid_x)
    dy = np.diff(grid_y)
    dz = np.diff(grid_z)
    
    if(len(dz)<1): dz = [L_z]

    first_moment = 0.0
    total_value = 0.0

    # Loop over all voxels
    for i in range(n_x):
        for j in range(n_y):
            for k in range(n_z):
                vol = (cell_size(dx, i, n_x) * cell_size(dy, j, n_y) * cell_size(dz, k, n_z))
                
                r = np.sqrt((grid_x[i] - center_x)**2 + (grid_y[j] - center_y)**2)
                if(r > r_max): continue

                value = data[i, j, k]

                first_moment += r * value * vol
                total_value += value * vol

    return first_moment, total_value

def calculate_second_moment(data, grid_x, grid_y, grid_z, first_moment, total_value):
    n_x, n_y, n_z = data.shape
    
    # Precompute Δ arrays
    dx = np.diff(grid_x)
    dy = np.diff(grid_y)
    dz = np.diff(grid_z)
    
    if(len(dz)<1): dz = [L_z]

    second_moment = 0.0

    # Loop over all voxels
    for i in range(n_x):
        for j in range(n_y):
            for k in range(n_z):
                vol = (cell_size(dx, i, n_x) * cell_size(dy, j, n_y) * cell_size(dz, k, n_z))
                
                r = np.sqrt((grid_x[i] - center_x)**2 + (grid_y[j] - center_y)**2)
                if(r > R_in): continue

                value = data[i, j, k]

                second_moment += ((r - first_moment)*(r - first_moment)) * value * vol

    second_moment = np.sqrt(second_moment / total_value)

    return second_moment

def calculate_optical_depth(r_centers, skin_depth_radial, area_factor=True):
    R_in = r_centers[-1]
    optical_depth = 0.0

    for i in range(1,len(r_centers)):
        r = r_centers[i]
        delta = skin_depth_radial[i]
        if(delta <= 0): continue

        if(area_factor): area_coefficient = r/(0.5*R_in)
        else: area_coefficient = 1.0

        optical_depth += (r - r_centers[i-1])*area_coefficient/delta

    return optical_depth

def write_metal_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set_metal, parameters, name_files):
    """Generate and write permittivity data files based on custom condition function."""
    points = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                function_set_metal(x, y, z, i, j, k, points, parameters)
    
    write_points(HARPS_DIR + "/input/" + name_files + "_metal_points.dat", points)

def reflector(x, y, z, i, j, k, points, param):
    center_z = param[0]
    b = param[1]
    y_reflect = param[2]

    if(y >= y_reflect and z > (center_z - b/2) and z < (center_z + b/2)):
        points.append((int(i), int(j), int(k)))


def read_total_power(path="Outputs/p_abs.txt"):
    try:
        with open(path, "r") as f:
            return float(f.read().strip())
    except FileNotFoundError:
        print(f"Warning: {path} not found, returning 0.")
        return 0.0
    except ValueError:
        print(f"Warning: invalid content in {path}, returning 0.")
        return 0.0
    
def run_harps(y_r, verbose=False):
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, reflector, [center_z, 1e5*b, y_r], "2D_run")


    #result = subprocess.run('mpirun --bind-to socket -np 1 ./harps.exe input/2D_run.in', shell=True, capture_output=True, text=True)
    result = subprocess.run('./harps.exe input/2D_run.in', shell=True, capture_output=True, text=True)
    if(verbose): 
        print("\nHARPS Output:"); print(result.stdout)
        print("\nHARPS Errors:"); print(result.stderr); print(f"Return code: {result.returncode}")

    try:
        with open("Outputs/p_abs.txt", "r") as f:
            pabs = 2 * float(f.read().strip())
    except FileNotFoundError:
        print("Warning: p_abs.txt not found, appending 0.")
        pabs = 0

    return pabs

def optimize_reflector(f, a, b, tol=3e-3, max_iter=15, flag_verbose=True):
    # Brent's method for 1D minimization (not done by me)
    # Cache to avoid re-evaluating expensive f
    cache = {}
    def eval_cached(x):
        # normalize x to reduce duplicate floating errors (optional)
        key = float(x)
        if key not in cache:
            if(flag_verbose): print(f"Running simulation for y_r = {key:.3f}")
            cache[key] = f(key)
            if(flag_verbose): print(f"p_abs for y_r = {cache[key]:.3f}")
            # store point in lists
        return cache[key]

    K = (3.0 - np.sqrt(5.0)) / 2.0   # ~0.38196601125 (used for golden steps fallback)

    # Initial points: pick x at middle
    x = w = v = a + 0.5 * (b - a)
    fx = fw = fv = eval_cached(x)

    d = e = b - a  # movement amounts

    for iteration in range(1, max_iter + 1):
        if(iteration == max_iter):
            if(flag_verbose): print("MAXIMUM ITERATIONS REACHED IN REFLECTOR OPTIMIZATION")

        m = 0.5 * (a + b)
        tol1 = tol * abs(x) + 1e-12
        tol2 = 2.0 * tol1

        # check stopping criterion
        if abs(x - m) <= (tol2 - 0.5 * (b - a)):
            break

        p = q = r = 0.0
        parabolic_accepted = False

        if abs(e) > tol1:
            # Attempt a parabolic interpolation step
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            if q != 0.0:
                if q > 0.0:
                    p = -p
                q = abs(q)
                if (abs(p) < abs(0.5 * q * e)) and (p > q * (a - x)) and (p < q * (b - x)):
                    d = p / q
                    u = x + d
                    if (u - a) < tol2 or (b - u) < tol2:
                        d = math.copysign(tol1, m - x)
                    parabolic_accepted = True

        if not parabolic_accepted:
            # Golden-section step
            if x < m:
                e = b - x
            else:
                e = a - x
            d = K * e

        # make sure step is at least tol1
        if abs(d) >= tol1:
            u = x + d
        else:
            u = x + math.copysign(tol1, d)

        fu = eval_cached(u)

        # Update a,b, and bookkeeping points
        if fu > fx:
            # new best
            if u < x:
                b = x
            else:
                a = x
            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if (fu >= fw) or (w == x):
                v, fv = w, fw
                w, fw = u, fu
            elif (fu > fv) or (v == x) or (v == w):
                v, fv = u, fu

    fx = eval_cached(x)

    return x, fx

def perform_run(flag_print_quantities = False):
    y_r, _ = optimize_reflector(run_harps, 0.07, 0.15, 1e-2, max_iter=12, flag_verbose=False)
    #get_integrated_Efield_empty = calibrate_empty_waveguide_E(R_in, y_r)
    pabs = run_harps(y_r, True)

    data_pabs, grid_x, grid_y, grid_z = load_3d_data("Outputs/absorbed_power_density.txt")
    r_centers, radial_pabs = get_radial_profile(data_pabs, grid_x, grid_y, grid_z, dr=0.0003, plot=True, title="absorbed_power_density")

    data_E_amplitude, grid_x, grid_y, grid_z = load_3d_data("Outputs/E_amplitude.txt")
    r_centers, radial_E_amplitude = get_radial_profile(data_E_amplitude, grid_x, grid_y, grid_z, dr=0.0003, plot=True, title="E_amplitude")

    data_ne, grid_x, grid_y, grid_z = load_3d_data("Outputs/electron_density.txt")
    r_centers, radial_ne = get_radial_profile(data_ne, grid_x, grid_y, grid_z, dr=0.0003, plot=False, title="electron_density")

    data_mu_real, grid_x, grid_y, grid_z = load_3d_data("Outputs/real_mobility.txt")
    r_centers, radial_mu_real = get_radial_profile(data_mu_real, grid_x, grid_y, grid_z, dr=0.0003, plot=False, title="mu_real")

    data_mu_imag, grid_x, grid_y, grid_z = load_3d_data("Outputs/imag_mobility.txt")
    r_centers, radial_mu_imag = get_radial_profile(data_mu_imag, grid_x, grid_y, grid_z, dr=0.0003, plot=False, title="mu_imag")

    # store in .txt file radius, E_field, p_abs
    with open(HARPS_DIR + "analysis/analysis_output/radial_profiles.txt", "w") as f:
        f.write("# Radius(m)    E_amplitude(V/m)    p_abs(W/m^3)\n")
        for i in range(len(r_centers)):
            f.write(f"{r_centers[i]:.6e},    {radial_E_amplitude[i]:.6e},    {radial_pabs[i]:.6e}\n")

    # Relevant quantities
    real_cond_first_moment, integrated_real_cond = calculate_first_moment(data_mu_real*data_ne*E_CHARGE, grid_x, grid_y, grid_z, r_max=R_in)
    skin_depth = calculate_skin_depth(radial_ne, radial_mu_real, radial_mu_imag)
    first_moment_pabs, total_pabs = calculate_first_moment(data_pabs, grid_x, grid_y, grid_z, r_max=R_in)
    total_pabs = read_total_power()
    second_moment = calculate_second_moment(data_pabs, grid_x, grid_y, grid_z, first_moment_pabs/total_pabs, total_pabs)
    optical_depth = calculate_optical_depth(r_centers, skin_depth, True)
    #optical_depth_no_area = calculate_optical_depth(r_centers, skin_depth, False)
    E_field_first_moment, integrated_E_field = calculate_first_moment(data_E_amplitude, grid_x, grid_y, grid_z, r_max=2*real_cond_first_moment/integrated_real_cond)


    plasma_size = real_cond_first_moment/integrated_real_cond
    integrated_real_cond *= 2;
    total_pabs = 2*total_pabs; 
    #total_pabs = pabs; # using the value from HARPS directly
    pabs_ratio = total_pabs / P_in; 
    max_pabs = np.max(data_pabs);
    max_pabs_radius = r_centers[np.argmax(radial_pabs)]
    avg_pabs_radius = 2*first_moment_pabs/total_pabs; 
    field_penetration_coeff = integrated_E_field/get_integrated_Efield_empty(2*plasma_size)
    avg_E_field_radius = E_field_first_moment/integrated_E_field
    E_field_radius_ratio = avg_E_field_radius / plasma_size

    avg_skin_depth = calculate_radial_average(r_centers, skin_depth, 1.25*plasma_size); 
    min_skin_depth = np.min(skin_depth[np.where(skin_depth>0)]);

    avg_real_cond_plasma = integrated_real_cond / (np.pi * plasma_size * plasma_size * L_z);
    penetration_coeff = avg_pabs_radius / plasma_size
    deposition_size_coefficient = second_moment / plasma_size

    if(flag_print_quantities):
        print("=== Simulation Results ===")
        print(f"  Integrated real conductivity: {integrated_real_cond:.6e} S m^2")
        print(f"  Total absorbed power: {total_pabs:.3f} W")
        print(f"  Reflector position y_r: {y_r:.6f} m")
        print(f"  Absorbed power ratio: {pabs_ratio:.6f}")
        print(f"  Average power density radius: {avg_pabs_radius*1e3:.6f} mm")
        print(f"  Power density second moment: {second_moment*1e3:.6f} mm")
        print(f"  Radius of maximum power density: {max_pabs_radius*1e3:.6f} mm")
        print(f"  Average skin depth: {avg_skin_depth*1e3:.6f} mm")
        print(f"  Minimum skin depth: {min_skin_depth*1e3:.6f} mm")
        print(f"  Plasma size: {plasma_size*1e3:.6f} mm")
        print(f"  Average real conductivity in plasma: {avg_real_cond_plasma:.6e} S/m")
        print(f"  Penetration coefficient: {penetration_coeff:.6f}")
        print(f"  Optical depth: {optical_depth:.6f}")
        print(f"  Field penetration coefficient: {field_penetration_coeff:.6f}")
        print(f"  E-field radius / plasma size: {E_field_radius_ratio:.6f}")
    
    return integrated_real_cond, pabs_ratio, penetration_coeff, deposition_size_coefficient, avg_skin_depth, plasma_size, avg_real_cond_plasma, optical_depth, y_r, field_penetration_coeff, E_field_radius_ratio, avg_E_field_radius



def harps_table_profile(df, r_max=0.014):
    """
    Interpolates radial profiles directly from tabular data and writes 3D input files.
    """
    r_data = df['r_m'].values
    ne_data = df['ne_m3'].values
    mu_r_data = df['mu_real_m2_Vs'].values
    mu_i_data = np.abs(df['mu_imag_m2_Vs'].values)  # Convert to positive magnitude for HARPS

    # Ensure profile starts at r = 0 by duplicating innermost point if missing
    if r_data[0] > 0.0:
        r_data = np.insert(r_data, 0, 0.0)
        ne_data = np.insert(ne_data, 0, ne_data[0])
        mu_r_data = np.insert(mu_r_data, 0, mu_r_data[0])
        mu_i_data = np.insert(mu_i_data, 0, mu_i_data[0])

    ne_interp = interp1d(r_data, ne_data, kind="linear", bounds_error=False, fill_value=0.0)
    mu_r_interp = interp1d(r_data, mu_r_data, kind="linear", bounds_error=False, fill_value=0.0)
    mu_i_interp = interp1d(r_data, mu_i_data, kind="linear", bounds_error=False, fill_value=0.0)

    write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [ne_interp, center_x, center_y, r_max], "2D_run", "elec_dens")
    write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [mu_r_interp, center_x, center_y, r_max], "2D_run", "real_mobi")
    write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [mu_i_interp, center_x, center_y, r_max], "2D_run", "imag_mobi")

def process_all_datasets():
    results = []
    
    for name, df in datasets.items():
        print(f"\n==========================================")
        print(f"Processing Dataset: {name}")
        print(f"==========================================")
        
        # 1. Generate grid files directly from dataset table
        harps_table_profile(df, r_max=R_in)
        
        # 2. Run HARPS simulation & optimize reflector position
        metrics = perform_run(flag_print_quantities=True)
        
        results.append({
            "Dataset": name,
            "Integrated_Real_Cond_S_m2": metrics[0],
            "Pabs_Ratio": metrics[1],
            "Penetration_Coeff": metrics[2],
            "Deposition_Size_Coeff": metrics[3],
            "Avg_Skin_Depth_mm": metrics[4] * 1e3,
            "Plasma_Size_mm": metrics[5] * 1e3,
            "Avg_Real_Cond_S_m": metrics[6],
            "Optical_Depth": metrics[7],
            "Reflector_Yr_m": metrics[8],
            "Field_Penetration_Coeff": metrics[9],
            "E_Field_Radius_Ratio": metrics[10]
        })

    summary_df = pd.DataFrame(results)
    summary_df.to_csv(HARPS_DIR + "analysis/analysis_output/datasets_summary.csv", index=False)
    return summary_df


N_x = 256; N_y = 256; N_z = 1
RefinementFactorX = 4; RefinementFactorY = 4; RefinementFactorZ = -1
x_BL = 0.023; x_RL = 0.031; x_RR = 0.057; x_BR = 0.063
y_BL = 0.026; y_RL = 0.034; y_RR = 0.058; y_BR = 0.066

# Define geometry parameters
a = 0.08636; b = 0.04318; R_in = 0.014; R_out = 0.015;
L_x = a; L_y = 0.1450; L_z = b; center_x = a/2; center_y = 0.046; center_z = b/2
symmetric = True

if(symmetric):
    x_offset = L_x/2; N_x = N_x//2; L_x = L_x/2; center_x = 0; x_BL = 0; x_RL = 0; x_RR = x_RR - L_x; x_BR = x_BR - L_x
    x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, RefinementFactorX, RefinementFactorY, RefinementFactorZ,
        x_BL, x_RL, x_RR, x_BR, y_BL, y_RL, y_RR, y_BR, x_grid=None, y_grid=None, z_grid=None)
else:
    x_offset = 0
    x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, RefinementFactorX, RefinementFactorY, RefinementFactorZ,
        x_BL, x_RL, x_RR, x_BR, y_BL, y_RL, y_RR, y_BR, x_grid=None, y_grid=None, z_grid=None)


get_integrated_Efield_empty = calibrate_empty_waveguide_E(R_in, 0.081)

# Execute the processing pipeline across all 5 datasets
summary = process_all_datasets()