import pandas as pd
import numpy as np
import subprocess
import os
import glob
import shutil
from datetime import datetime
import re
import math
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.lines as mlines
from scipy.interpolate import interp1d

import_1d_profiles = False

# Configuration
HARPS_EXE = "./harps.exe"  # Path to HARPS executable
HARPS_DIR = "./"
INPUT_TEMPLATE = "input/2D_run.in"  # Base configuration file
OUTPUT_BASE = "Outputs/2D_runs"  # Base directory for results

free = -0.001234567  # Arbitrary small value for free space excitation
EPSILON_0 = 8.854187817e-12
MU_0 = 4 * np.pi * 1e-7
E_CHARGE = 1.602176634e-19
C_LIGHT = 1/np.sqrt(MU_0*EPSILON_0)

angular_frequency = 2 * np.pi * 2.45e9  # Example frequency: 2.45 GHz
real_permittivity = 1.0  # Free space

flag_print_quantities = True
flag_center_coordinates = False

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

def harps_table_profile(df, r_max=0.014):
    r_data = df['r_m'].values
    ne_data = df['ne_m3'].values
    mu_r_data = df['mu_real_m2_Vs'].values
    mu_i_data = np.abs(df['mu_imag_m2_Vs'].values)  # Convert to positive magnitude for HARPS

    if r_data[0] > 0.0:
        r_data = np.insert(r_data, 0, 0.0)
        ne_data = np.insert(ne_data, 0, ne_data[0])
        mu_r_data = np.insert(mu_r_data, 0, mu_r_data[0])
        mu_i_data = np.insert(mu_i_data, 0, mu_i_data[0])

    ne_interp = interp1d(r_data, ne_data, kind="linear", bounds_error=False, fill_value=0.0)
    mu_r_interp = interp1d(r_data, mu_r_data, kind="linear", bounds_error=False, fill_value=mu_r_data[-1])
    mu_i_interp = interp1d(r_data, mu_i_data, kind="linear", bounds_error=False, fill_value=mu_i_data[-1])

    write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [ne_interp, center_x, center_y, r_max], "2D_run", "elec_dens")
    write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [mu_r_interp, center_x, center_y, r_max], "2D_run", "real_mobi")
    write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [mu_i_interp, center_x, center_y, r_max], "2D_run", "imag_mobi")

def process_all_datasets():
    dataset_results = []
    for name, df in datasets.items():
        print(f"\n==========================================")
        print(f"Processing Dataset: {name}")
        print(f"==========================================")
        
        harps_table_profile(df, r_max=R_in)
        res = perform_run(flag_print_quantities=True, custom_profile=True)
        dataset_results.append(res)

    return np.array(dataset_results)

def process_1d_directory(directory_path="testing/self_consistent/harps"):
    """Reads all 1D profile txt files in the directory (excluding those with 'mm' in the filename), runs HARPS, and returns array of results."""
    pattern = os.path.join(directory_path, "*.txt")
    files = sorted(glob.glob(pattern))
    
    if not files:
        # Fallback in case files don't have .txt extension
        files = sorted(glob.glob(os.path.join(directory_path, "*")))

    # Filter out files that contain 'mm' in their filename
    files = [f for f in files if "mm" not in os.path.basename(f)]
    files = [f for f in files if "5_slm" not in os.path.basename(f)]
    files = [f for f in files if "600mbar_600W" not in os.path.basename(f)]
    files = [f for f in files if "100mbar_600W" not in os.path.basename(f)]
    files = [f for f in files if "1200W" not in os.path.basename(f)]
    files = [f for f in files if "2000W" not in os.path.basename(f)]

    results = []
    for filepath in files:
        if os.path.isdir(filepath):
            continue

        print(f"\n==========================================")
        print(f"Processing 1D Self-Consistent Profile: {filepath}")
        print(f"==========================================")

        try:
            ne_interp, mu_r_interp, mu_i_interp = load_1d_profiles(filepath)

            write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [ne_interp, center_x, center_y, R_in], "2D_run", "elec_dens")
            write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [mu_r_interp, center_x, center_y, R_in], "2D_run", "real_mobi")
            write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [mu_i_interp, center_x, center_y, R_in], "2D_run", "imag_mobi")

            res = perform_run(flag_print_quantities=True, custom_profile=True)
            results.append(res)
        except Exception as e:
            print(f"Error processing {filepath}: {e}")

    return np.array(results)

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "lines.linewidth": 1.5
})

def plot_paper_ready_scatter(x_label, x_data, quantities, filename, 
                             highlights_x=None, highlights_y=None, 
                             highlights_green_x=None, highlights_green_y=None):
    n = len(quantities)
    if n == 0:
        return

    if n == 1:
        nrows, ncols = 1, 1
    elif n == 2:
        nrows, ncols = 1, 2
    else:
        nrows, ncols = 2, 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows), sharex=True)
    axes = np.atleast_1d(axes).flatten()

    sublabels = ['(a)', '(b)', '(c)', '(d)']

    x_data = np.asarray(x_data)
    mask = x_data >= 4e-8
    x_filtered = x_data[mask]

    for i, (y_label, y_data) in enumerate(quantities.items()):
        y_data = np.asarray(y_data)

        if i == 0:
            y_filtered = y_data[mask]
        else:
            y_smooth = y_data.copy()
            upper = 1.03
            lower = 0.95
            transition = (y_smooth >= lower) & (y_smooth <= upper)
            y_smooth[transition] = lower + (y_smooth[transition] - lower) * (1 - lower) / (upper - lower)
            y_smooth[y_smooth > upper] = 1.0
            y_filtered = y_smooth[mask]

        ax = axes[i]
        
        # Blue: Randomly Generated Plasma Profiles
        ax.scatter(x_filtered, y_filtered, color="tab:blue", alpha=0.25, s=5, rasterized=False, zorder=1)

        bins = np.logspace(np.log10(x_filtered.min()), np.log10(x_filtered.max()), 25)

        bin_median, bin_edges, _ = stats.binned_statistic(
            x_filtered, y_filtered, statistic='median', bins=bins
        )
        bin_75, _, _ = stats.binned_statistic(x_filtered, y_filtered, statistic=lambda x: np.percentile(x, 75), bins=bins)
        bin_25, _, _ = stats.binned_statistic(x_filtered, y_filtered, statistic=lambda x: np.percentile(x, 25), bins=bins)

        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        ax.fill_between(bin_centers, bin_25, bin_75, color='black', alpha=0.2, zorder=2)
        ax.plot(bin_centers, bin_median, color='black', lw=2, zorder=3)

        # Red: Self-Consistently Calculated Plasma Profiles
        if highlights_x and highlights_y and y_label in highlights_y:
            hx = highlights_x[x_label] if x_label in highlights_x else highlights_x[list(highlights_x.keys())[0]]
            hy = highlights_y[y_label]
            ax.scatter(hx, hy, color="firebrick", s=40, marker='o', edgecolors='white', linewidths=0.6, zorder=6)
        
        # Green: Experimentally Obtained Plasma Profiles
        if highlights_green_x and highlights_green_y and y_label in highlights_green_y:
            hx = highlights_green_x[x_label] if x_label in highlights_green_x else highlights_green_x[list(highlights_green_x.keys())[0]]
            hy = highlights_green_y[y_label]
            ax.scatter(hx, hy, color="forestgreen", s=40, marker='s', edgecolors='black', linewidths=0.6, zorder=6)

        ax.set_xscale('log')
        ax.set_ylabel(y_label)
        ax.grid(True, which="both", ls=":", alpha=0.5)
        ax.text(0.05, 0.88, sublabels[i], transform=ax.transAxes, fontweight='bold')

        if i >= (n - ncols):
            ax.set_xlabel(x_label)

    # Clean up empty subplots if any
    for j in range(n, len(axes)):
        fig.delaxes(axes[j])

    # Custom Legend Proxy Handles for top legend
    legend_elements = [
        mlines.Line2D([0], [0], marker='o', color='w', label='Randomly Generated',
                      markerfacecolor='tab:blue', markersize=6, alpha=0.6),
        mlines.Line2D([0], [0], marker='o', color='w', label='Self-Consistently Calculated',
                      markerfacecolor='firebrick', markeredgecolor='black', markeredgewidth=0.5, markersize=8),
        mlines.Line2D([0], [0], marker='s', color='w', label='Experimentally Obtained',
                      markerfacecolor='forestgreen', markeredgecolor='black', markeredgewidth=0.5, markersize=8)
    ]

    # Place Legend Box at the top above subplots
    fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 1.0),
               ncol=3, frameon=True, framealpha=0.95, edgecolor='gray', fontsize=11)

    # Adjust layout to accommodate top legend without clipping
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    plt.savefig(filename + '.pdf', bbox_inches='tight')
    plt.close(fig)


def load_1d_profiles(filename):
    data = np.loadtxt(filename)

    radius     = data[:, 2]   # radius
    elec_dens  = data[:, 4]   # ne
    mu_real    = data[:, 3]   # real mobility
    mu_imag    = data[:, 8]   # imaginary mobility

    T_g        = data[:, 7]   # gas temperature

    print("Max ne:", np.max(elec_dens))
    print("Max T_g:", np.max(T_g))

    idx = np.argsort(radius)
    radius    = radius[idx]
    elec_dens = elec_dens[idx]
    mu_real   = mu_real[idx]
    mu_imag   = mu_imag[idx]

    ne_interp = interp1d(radius, elec_dens, kind="linear", bounds_error=False, fill_value=0.0)
    mu_r_interp = interp1d(radius, mu_real, kind="linear", bounds_error=False, fill_value=mu_real[-1])
    mu_i_interp = interp1d(radius, mu_imag, kind="linear", bounds_error=False, fill_value=mu_imag[-1])

    return ne_interp, mu_r_interp, mu_i_interp

def interpolated_radial_profile(x, y, z, i, j, k, points, values, param):
    interp_func = param[0]
    center_x    = param[1]
    center_y    = param[2]
    r_max       = param[3]

    r = np.sqrt((x - center_x)**2 + (y - center_y)**2)

    if r > r_max:
        return

    val = float(interp_func(r))

    if val <= 0:
        return

    points.append((int(i), int(j), int(k)))
    values.append(val)


def write_complex_field(filename, field_values, component):
    with open(filename, 'w') as f:
        f.write("{")
        for i, field in enumerate(field_values):
            if component == 0: f.write(f"[Complex({field.real},{field.imag});F;F]\n")
            if component == 1: f.write(f"[F;Complex({field.real},{field.imag});F]\n")
            if component == 2: f.write(f"[F;F;Complex({field.real},{field.imag})]\n")

            if i < len(field_values) - 1:
                f.write("  ")
        f.write("}")

def write_complex_field_full(filename, field_values):
    with open(filename, 'w') as f:
        f.write("{")
        for i, field in enumerate(field_values):
            f.write("[")
            for j in range(3):
                if j < len(field) and np.abs(field[j].imag - free) >= 1e-9:
                    f.write(f"Complex({field[j].real},{field[j].imag})")
                else:
                    f.write("F")
                
                if j < 2:
                    f.write(";")
            f.write("]\n")

            if i < len(field_values) - 1:
                f.write("  \n")
        f.write("}")
        
def write_values(filename, values):
    with open(filename, 'w') as f:
        f.write("{")
        for i, value in enumerate(values):
            f.write(f"({value:.4e})")
            if i < len(values) - 1:
                f.write(" \n")
        f.write("}")

def write_points(filename, points):
    with open(filename, 'w') as f:
        f.write("{")
        for i, (x, y, z) in enumerate(points):
            f.write(f"({int(x)}, {int(y)}, {int(z)})\n")
        f.write("}")

def write_real_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set, parameters, name_files, name_var):
    points = []
    values = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                function_set(x, y, z, i, j, k, points, values, parameters)
    
    write_values(HARPS_DIR + "/input/" + name_files + "_" + name_var + "_values.dat", values)
    write_points(HARPS_DIR + "/input/" + name_files + "_" + name_var + "_points.dat", points)
    
def write_excitation_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set_excitation, parameters, name_files):
    points = []
    field_values = []

    component = parameters[0]
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                fields = function_set_excitation(x, y, z, parameters)
                if fields == None: continue
                
                field_values.append(fields)
                points.append((i, j, k))

    if(component < 3): write_complex_field(HARPS_DIR + "/input/"  + name_files + "_excitation_values.dat", field_values, component)
    else: write_complex_field_full(HARPS_DIR + "/input/"  + name_files + "_excitation_values.dat", field_values)
    write_points(HARPS_DIR + "/input/" + name_files + "_excitation_points.dat", points)

def write_metal_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set_metal, parameters, name_files):
    points = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                function_set_metal(x, y, z, i, j, k, points, parameters)
    
    write_points(HARPS_DIR + "/input/" + name_files + "_metal_points.dat", points)


def mod_gaussian_profile(x, y, z, i, j, k, points, values, param):
    ne_0 = param[0]
    center_x = param[1]
    center_y = param[2]
    center_z = param[3]
    sigma_r = param[4]
    sigma_z = param[5]
    shift_r = param[6]

    radius = np.sqrt((x-center_x)*(x-center_x) + (y-center_y)*(y-center_y))

    if(radius<=0.014):
        if(sigma_z > 0):
            ne = ne_0*np.exp(-0.5*((radius-shift_r)/sigma_r)**2)*np.exp(-0.5*(np.abs(z-center_z)/sigma_z)**2)
        else:
            ne = ne_0*np.exp(-0.5*((radius-shift_r)/sigma_r)**2)
    else:
        ne = 1.01e11

    if(ne >= 1e11):
        points.append((int(i), int(j), int(k)))
        values.append(ne)

def mu_real_profile(x, y, z, i, j, k, points, values, param):
    mu_real_0 = param[0]
    center_x = param[1]
    center_y = param[2]
    p = param[3]
    Tg_0 = param[4]
    sigma_Tg = param[5]

    radius = np.sqrt((x-center_x)*(x-center_x) + (y-center_y)*(y-center_y))

    if(radius<=0.014):
        temperature = Tg_0*np.exp(-0.5*(radius/sigma_Tg)**2)
        n_number = p/(1.380649e-23*temperature)
        mu_real = 5.5*np.exp(-np.log(n_number/mu_real_0)*np.log(n_number/mu_real_0)/2/2/2)
    else:
        return
    
    points.append((int(i), int(j), int(k)))
    values.append(mu_real)

def mu_imag_profile(x, y, z, i, j, k, points, values, param):
    mu_imag_0 = param[0]
    center_x = param[1]
    center_y = param[2]
    p = param[3]
    Tg_0 = param[4]
    sigma_Tg = param[5]

    radius = np.sqrt((x-center_x)*(x-center_x) + (y-center_y)*(y-center_y))

    if(radius<=0.014):
        temperature = Tg_0*np.exp(-0.5*(radius/sigma_Tg)**2)
        n_number = p/(1.380649e-23*temperature)
        mu_imag = -17+17/(1+np.exp(-(n_number-1.2e23)/mu_imag_0))
    else:
        return
    
    points.append((int(i), int(j), int(k)))
    values.append(mu_imag)

def ring_interface_3D(x, y, z, i, j, k, points, values, param):
    epsilon_r = param[0]
    R_inner = param[1]
    R_outer = param[2]
    center_x = param[3]
    center_y = param[4]

    radius = np.sqrt((x-center_x)*(x-center_x) + (y-center_y)*(y-center_y))

    if(radius > R_inner and radius < R_outer):
        points.append((int(i), int(j), int(k)))
        values.append(epsilon_r)

def reflector(x, y, z, i, j, k, points, param):
    center_z = param[0]
    b = param[1]
    y_reflect = param[2]

    if(y >= y_reflect and z > (center_z - b/2) and z < (center_z + b/2)):
        points.append((int(i), int(j), int(k)))

def square_waveguide_excitation_3D_BC(x, y, z, parameters):
    alpha = parameters[1]
    a = parameters[2]
    x_0 = parameters[3]

    if (y > 1e-6): return None

    return [complex(0, 0), complex(0, 0), complex(-alpha*np.sin(np.pi*(x+x_0)/a), 0)]


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
                raise ValueError("Key points must be monotonically increasing: BL <= RL <= RR <= BR.")
            
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
            arr = np.linspace(0, length, n)
            d = np.diff(arr)

            return arr, d, False

    x, dx, non_uniform_x = process_direction(n_x, length_x, refinementFactor_x, x_BL, x_BR, x_RL, x_RR, x_grid, x_file)
    y, dy, non_uniform_y = process_direction(n_y, length_y, refinementFactor_y, y_BL, y_BR, y_RL, y_RR, y_grid, y_file)
    z, dz, non_uniform_z = process_direction(n_z, length_z, refinementFactor_z, z_BL, z_BR, z_RL, z_RR, z_grid, z_file)

    non_uniform_grid = non_uniform_x or non_uniform_y or non_uniform_z
    return x, dx, y, dy, z, dz, non_uniform_grid

def cell_size(d, idx, n):
        if n == 1:
            return d[0]
        if idx == 0:
            return 0.5 * d[0]
        if idx == n - 1:
            return 0.5 * d[-1]
        return 0.5 * (d[idx] + d[idx - 1])

def update_input_file(filename, params):
    with open(filename, 'r') as f:
        text = f.read()

    for key, value in params.items():
        text = re.sub(
            rf'({key}\s*=\s*)([0-9Ee\.\+-]+)',
            rf'\g<1>{value}',
            text
        )

    with open(filename, 'w') as f:
        f.write(text)


def run_harps(y_r, verbose=False, custom_profile=False):
    params = { "Nx": N_x, "Ny": N_y, "Nz": N_z, "LengthX": L_x, "LengthY": L_y, "LengthZ": L_z,
            "RefinementFactorX": RefinementFactorX, "x_BL": x_BL, "x_RL": x_RL, "x_RR": x_RR, "x_BR": x_BR,
            "RefinementFactorY": RefinementFactorY, "y_BL": y_BL, "y_RL": y_RL, "y_RR": y_RR, "y_BR": y_BR,
            "RealMobility": RealMobility, "ImagMobility": ImagMobility}
    update_input_file("input/2D_run.in", params)

    if not custom_profile:
        if import_1d_profiles:
            ne_interp, mu_r_interp, mu_i_interp = load_1d_profiles("Outputs/1d_run_output.txt")
            write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [ne_interp, center_x, center_y, 0.014], "2D_run", "elec_dens")
            write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [mu_r_interp, center_x, center_y, 0.014], "2D_run", "real_mobi")
            write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, interpolated_radial_profile, [mu_i_interp, center_x, center_y, 0.014], "2D_run", "imag_mobi")
        else:
            write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile, [ne_0, center_x, center_y, center_z, ne_sigma_r, -1, ne_shift_r], "2D_run", "elec_dens")
            write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mu_real_profile, [mu_R_0, center_x, center_y, pressure, Tg_0, ne_sigma_r], "2D_run", "real_mobi")
            write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mu_imag_profile, [mu_I_N, center_x, center_y, pressure, Tg_0, ne_sigma_r], "2D_run", "imag_mobi")

    write_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, ring_interface_3D, [tube_perm, R_in, R_out, center_x, center_y], "2D_run", "permittivity")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,E_0, a, x_offset], "2D_run")
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,     reflector, [center_z, 1e5*b, y_r], "2D_run")

    result = subprocess.run('./harps.exe input/2D_run.in', shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print("\nHARPS Error Output:")
        print(result.stderr)
        raise RuntimeError(f"HARPS simulation failed with exit code {result.returncode}")

    if verbose: 
        print("\nHARPS Output:"); print(result.stdout)

    try:
        with open("Outputs/p_abs.txt", "r") as f:
            pabs = 2 * float(f.read().strip())
    except FileNotFoundError:
        print("Warning: p_abs.txt not found, appending 0.")
        pabs = 0

    return pabs

def optimize_reflector(f, a, b, tol=3e-3, max_iter=15, flag_verbose=False, custom_profile=False):
    cache = {}
    def eval_cached(x):
        key = float(x)
        if key not in cache:
            if(flag_verbose): print(f"Running simulation for y_r = {key:.3f}")
            cache[key] = f(key, custom_profile=custom_profile)
            if(flag_verbose): print(f"p_abs for y_r = {cache[key]:.3f}")
        return cache[key]

    K = (3.0 - np.sqrt(5.0)) / 2.0

    x = w = v = a + 0.5 * (b - a)
    fx = fw = fv = eval_cached(x)

    d = e = b - a

    for iteration in range(1, max_iter + 1):
        if(iteration == max_iter):
            if(flag_verbose): print("MAXIMUM ITERATIONS REACHED IN REFLECTOR OPTIMIZATION")

        m = 0.5 * (a + b)
        tol1 = tol * abs(x) + 1e-12
        tol2 = 2.0 * tol1

        if abs(x - m) <= (tol2 - 0.5 * (b - a)):
            break

        p = q = r = 0.0
        parabolic_accepted = False

        if abs(e) > tol1:
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
            if x < m:
                e = b - x
            else:
                e = a - x
            d = K * e

        if abs(d) >= tol1:
            u = x + d
        else:
            u = x + math.copysign(tol1, d)

        fu = eval_cached(u)

        if fu > fx:
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


def load_3d_data(filename):
    try:
        with open(filename) as f:
            n_x, n_y, n_z, length_x, length_y, length_z = map(float, f.readline().split())        
            n_x, n_y, n_z = int(n_x), int(n_y), int(n_z)

            header = f.readline().strip()

            if header == "NON_UNIFORM_GRID":
                x_coords = np.array(list(map(float, f.readline().split())))
                y_coords = np.array(list(map(float, f.readline().split())))
                z_coords = np.array(list(map(float, f.readline().split())))
                
                data = np.loadtxt(f)
                data = data.reshape((n_x, n_y, n_z))

                if(flag_center_coordinates == True):
                    x_coords = x_coords - center_x
                    y_coords = y_coords - center_y
                    z_coords = z_coords - center_z

                return data, x_coords, y_coords, z_coords
            
            else:
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
    
def calculate_skin_depth(ne, mu_real, mu_imag):
    complex_conductivity = E_CHARGE * ne * (mu_real + 1j * mu_imag)

    sigma_real = np.real(complex_conductivity)
    sigma_imag = np.imag(complex_conductivity)

    k2_re = angular_frequency*angular_frequency * MU_0 * EPSILON_0 - angular_frequency * MU_0 * sigma_imag
    k2_im = - angular_frequency * MU_0 * sigma_real

    k_imag = np.sqrt((np.sqrt(k2_re*k2_re + k2_im*k2_im) - k2_re) / 2.0)
    k_imag = np.where(k2_im < 0, -k_imag, k_imag)

    skin_depth = np.where(k_imag != 0, 1.0 / np.abs(k_imag), np.inf)

    return skin_depth

def get_radial_profile(data, grid_x, grid_y, grid_z, dr = 0.0003, plot=False, title="power_density"):
    n_x, n_y, n_z = data.shape
    
    dx = np.diff(grid_x)
    dy = np.diff(grid_y)
    dz = np.diff(grid_z)

    if(len(dz)<1): dz = [L_z]
    
    r_bins = np.arange(0, R_in + dr, dr)
    r_centers = 0.5 * (r_bins[:-1] + r_bins[1:])
    radial_quantity = np.zeros_like(r_centers)
    radial_volume = np.zeros_like(r_centers)

    X, Y = np.meshgrid(grid_x, grid_y, indexing='ij')
    R = np.sqrt((X - center_x)**2 + (Y - center_y)**2)
    
    for i in range(n_x):
        for j in range(n_y):
            for k in range(n_z):
                vol = (cell_size(dx, i, n_x) * cell_size(dy, j, n_y) * cell_size(dz, k, n_z))
                r = R[i, j]
                bin_idx = np.searchsorted(r_bins, r) - 1
                if 0 <= bin_idx < len(r_centers):
                    radial_quantity[bin_idx] += data[i, j, k] * vol
                    radial_volume[bin_idx] += vol

    with np.errstate(divide='ignore', invalid='ignore'):
        radial_quantity_density = np.divide(radial_quantity, radial_volume, out=np.zeros_like(radial_quantity), where=radial_volume > 0)

    if plot:
        plt.figure()
        plt.plot(r_centers*1e3, radial_quantity_density, label=title)
        plt.xlabel('Radius (mm)')
        plt.ylabel(title)
        plt.grid()
        plt.savefig(HARPS_DIR + "analysis/analysis_output/" + title + "_radial.png", dpi=300)

    return r_centers, radial_quantity_density

def read_total_power(path="Outputs/power_poynting_total.txt"):
    try:
        with open(path, "r") as f:
            return float(f.read().strip())
    except FileNotFoundError:
        print(f"Warning: {path} not found, returning 0.")
        return 0.0
    except ValueError:
        print(f"Warning: invalid content in {path}, returning 0.")
        return 0.0

def calculate_radial_average(r_centers, radial_data, r_plasma=0.014):
    dr = r_centers[1] - r_centers[0]

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
    
    dx = np.diff(grid_x)
    dy = np.diff(grid_y)
    dz = np.diff(grid_z)
    
    if(len(dz)<1): dz = [L_z]

    first_moment = 0.0
    total_value = 0.0

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
    
    dx = np.diff(grid_x)
    dy = np.diff(grid_y)
    dz = np.diff(grid_z)
    
    if(len(dz)<1): dz = [L_z]

    second_moment = 0.0

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

def calibrate_empty_waveguide_E(R_in, y_r):
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, reflector, [center_z, 1e5*b, y_r], "2D_run")
    
    result = subprocess.run('./harps.exe input/2D_run_no_plasma.in', shell=True, capture_output=True, text=True)
    
    data_E_amplitude, grid_x, grid_y, grid_z = load_3d_data("Outputs/E_amplitude.txt")

    radii_samples = [0, R_in / 4, R_in / 2, R_in, 2 * R_in]
    integrated_values = [0]

    for r_sample in radii_samples[1:]:
        _, integrated_E_field = calculate_first_moment( data_E_amplitude, grid_x, grid_y, grid_z, r_max=r_sample)
        integrated_values.append(integrated_E_field)

    return interp1d(radii_samples, integrated_values, kind="cubic", fill_value="extrapolate")

def randomize_params(flag_print=False):
    globals().update({
        name: (np.random.uniform(low, high) if flag == 0 else
         10 ** np.random.uniform(np.log10(low), np.log10(high)))
        for name, (low, high, flag) in param_ranges.items()
    })
    if flag_print:
        print("Randomized parameters:")
        for name in param_ranges.keys():
            print(f"  {name} = {globals()[name]}")

def perform_run(flag_print_quantities = False, custom_profile=False):
        y_r, _ = optimize_reflector(run_harps, 0.07, 0.15, 1e-2, max_iter=12, flag_verbose=False, custom_profile=custom_profile)
        pabs = run_harps(y_r, True, custom_profile=custom_profile)

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

        with open(HARPS_DIR + "analysis/analysis_output/radial_profiles.txt", "w") as f:
            f.write("# Radius(m)    E_amplitude(V/m)    p_abs(W/m^3)\n")
            for i in range(len(r_centers)):
                f.write(f"{r_centers[i]:.6e},    {radial_E_amplitude[i]:.6e},    {radial_pabs[i]:.6e}\n")

        real_cond_first_moment, integrated_real_cond = calculate_first_moment(data_mu_real*data_ne*E_CHARGE, grid_x, grid_y, grid_z, r_max=R_in)
        imag_cond_first_moment, integrated_imag_cond = calculate_first_moment(data_mu_imag*data_ne*E_CHARGE, grid_x, grid_y, grid_z, r_max=R_in)
        skin_depth = calculate_skin_depth(radial_ne, radial_mu_real, radial_mu_imag)
        first_moment_pabs, total_pabs = calculate_first_moment(data_pabs, grid_x, grid_y, grid_z, r_max=R_in)
        total_pabs = pabs
        second_moment = calculate_second_moment(data_pabs, grid_x, grid_y, grid_z, first_moment_pabs/total_pabs, total_pabs)
        optical_depth = calculate_optical_depth(r_centers, skin_depth, True)
        E_field_first_moment, integrated_E_field = calculate_first_moment(data_E_amplitude, grid_x, grid_y, grid_z, r_max=2*real_cond_first_moment/integrated_real_cond)

        plasma_size = real_cond_first_moment/integrated_real_cond
        integrated_real_cond *= 2; integrated_imag_cond *= 2
        integrated_conductivity = np.sqrt(integrated_real_cond**2 + integrated_imag_cond**2)
        total_pabs = total_pabs; 
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
            print(f"  Integrated imaginary conductivity: {integrated_imag_cond:.6e} S m^2")
            print(f"  Integrated conductivity: {integrated_conductivity:.6e} S m^2")
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
        
        return integrated_real_cond, abs(integrated_imag_cond), integrated_conductivity, min(pabs_ratio, 0.999), penetration_coeff, deposition_size_coefficient, avg_skin_depth, plasma_size, avg_real_cond_plasma, optical_depth, y_r, min(field_penetration_coeff, 0.999), E_field_radius_ratio, avg_E_field_radius


if __name__ == "__main__":
    # Define grid parameters
    N_x = 256; N_y = 256; N_z = 1
    RefinementFactorX = 4; RefinementFactorY = 4; RefinementFactorZ = -1
    x_BL = 0.023; x_RL = 0.031; x_RR = 0.057; x_BR = 0.063
    y_BL = 0.026; y_RL = 0.034; y_RR = 0.058; y_BR = 0.066

    # Define geometry parameters
    a = 0.08636; b = 0.04318; R_in = 0.014; R_out = 0.015;
    L_x = a; L_y = 0.1450; L_z = b; center_x = a/2; center_y = 0.046; center_z = b/2
    symmetric = True

    omega_c = C_LIGHT * np.pi / a
    impedance = np.sqrt(MU_0/ EPSILON_0)
    impedance_TE = impedance / np.sqrt(1 - (omega_c / angular_frequency)**2)

    # Generate grid with the correct geometry
    if(symmetric):
        x_offset = L_x/2; N_x = N_x//2; L_x = L_x/2; center_x = 0; x_BL = 0; x_RL = 0; x_RR = x_RR - L_x; x_BR = x_BR - L_x
        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, RefinementFactorX, RefinementFactorY, RefinementFactorZ,
            x_BL, x_RL, x_RR, x_BR, y_BL, y_RL, y_RR, y_BR, x_grid=None, y_grid=None, z_grid=None)
    else:
        x_offset = 0
        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, RefinementFactorX, RefinementFactorY, RefinementFactorZ,
            x_BL, x_RL, x_RR, x_BR, y_BL, y_RL, y_RR, y_BR, x_grid=None, y_grid=None, z_grid=None)
            

    # Define parameter arrays
    P_in = 770.1905  # Input power in Watts
    tube_perm = 4.5;    E_0 = np.sqrt(P_in*4*impedance_TE/(a*b))
    RealMobility = 0.8; ImagMobility = -2.5
    y_r = 0.103
    
    ne_0 = 4.5e18;        ne_sigma_r = 0.0035;      ne_shift_r = 0.000
    mu_R_0 = 0.3;         mu_I_N = -4.5e23;         pressure = 30000;       Tg_0 = 5000  

    # Run without plasma to calibrate E-coefficient
    get_integrated_Efield_empty = calibrate_empty_waveguide_E(R_in, 0.081)

    #y_r, _ = optimize_reflector(run_harps, 0.08, L_y, 1e-2, max_iter=10, flag_verbose=True)
    #pabs = run_harps(y_r, True)
    #perform_run(True)



    # ====================
    # === Perform runs ===
    # ====================
    N_runs = 7002

    param_ranges = {
        "ne_0": (0.4e18, 100e18,1), "ne_sigma_r": (0.0005, 0.0065,0), "ne_shift_r": (5e-5, 6e-3,1),
        "mu_R_0": (5e22, 4e23,0), "mu_I_N": (5e22, 4e23,0), "pressure": (2000, 100000,1), "Tg_0": (800, 9000,0)}

    results = []
    for i in range(N_runs):
        print(f"\n=== RUN {i+1}/{N_runs} ===")
        randomize_params(True)

        # Run the main function
        try:
            res = perform_run(True)
            results.append(res)
        except Exception as e:
            print(f"Run {i+1} failed: {e}")

    results = np.array(results)

    np.save(HARPS_DIR + f"analysis/analysis_output/results7002.npy", results)
    #results = np.load(HARPS_DIR + f"analysis/analysis_output/results7001.npy")

    if results.shape[0] == 0:
        raise RuntimeError("No successful runs!")

    x_axes = {
        "Integrated real conductivity [S·m²]": results[:, 0],
        "Integrated imaginary conductivity [S·m²]": results[:, 1],
        "Integrated conductivity [S·m²]": results[:, 2],
        "Optical Depth": results[:, 9],
    }
    quantities = {
        "Absorbed power ratio": results[:, 3],
        "Field Penetration Coefficient": results[:, 11],
    }

    # === Red Highlights (Process 1D Self-Consistent Profiles) ===
    run_1d_self_consistent = False  # Set to True to re-run simulations
    red_results_path = HARPS_DIR + f"analysis/analysis_output/red1.npy"

    if run_1d_self_consistent or not os.path.exists(red_results_path):
        red_results = process_1d_directory("testing/self_consistent/harps")
        np.save(red_results_path, red_results)
    else:
        red_results = np.load(red_results_path)

    highlights_x = {
        "Integrated real conductivity [S·m²]": red_results[:, 0],
        "Integrated imaginary conductivity [S·m²]": red_results[:, 1],
        "Integrated conductivity [S·m²]": red_results[:, 2],
        "Optical Depth": red_results[:, 9],
    }
    highlights_y = {
        "Absorbed power ratio": red_results[:, 3],
        "Field Penetration Coefficient": red_results[:, 11],
    }

    # === Green Highlights (Process Experimental Datasets) ===
    run_green = False  # Set to True to re-run simulations
    green_results_path = HARPS_DIR + f"analysis/analysis_output/green1.npy"
    if run_green or not os.path.exists(green_results_path):
        dataset_results = process_all_datasets()
        np.save(green_results_path, dataset_results)
    else:
        dataset_results = np.load(green_results_path)

    highlights_green_x = {
        "Integrated real conductivity [S·m²]": dataset_results[:, 0],
        "Integrated imaginary conductivity [S·m²]": dataset_results[:, 1],
        "Integrated conductivity [S·m²]": dataset_results[:, 2],
        "Optical Depth": dataset_results[:, 9],
    }
    highlights_green_y = {
        "Absorbed power ratio": dataset_results[:, 3],
        "Field Penetration Coefficient": dataset_results[:, 11],
    }

    for x_label, x_data in x_axes.items():
        safe_name = x_label.replace(' ', '_').replace('(', '').replace(')', '').replace('·', '')
        fname = HARPS_DIR + f"analysis/analysis_output/paper_fig_{safe_name}"
        plot_paper_ready_scatter(x_label, x_data, quantities, fname,
            highlights_x=highlights_x, highlights_y=highlights_y, highlights_green_x=highlights_green_x, highlights_green_y=highlights_green_y)