import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.colors import LogNorm
import argparse
from mpl_toolkits.mplot3d import Axes3D

#Defining center of glass tube
flag_center_coordinates = False
contour = False

center_x = 0
center_y = 0.046
center_z = 0.0

folder = "Outputs/"
#folder = "Outputs/temp/"


plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "stix",         # Makes math look like Times New Roman
})


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

                # Remove points where z > 0.5
                z_mask = z_coords <= 0.5
                z_coords = z_coords[z_mask]
                data = data[:, :, z_mask]

                return data, x_coords, y_coords, z_coords
        
    except FileNotFoundError:
        print(f"Warning: File '{filename}' not found. Skipping.")
        return None, None, None, None


def plot_3d_scatter(data, output_image, title, x_vals, y_vals, z_vals, xlabel="X [cm]", ylabel="Y [cm]", zlabel="Z [cm]"):
    n_x, n_y, n_z = data.shape

    if n_x == 1 or n_y == 1 or n_z == 1:
        print("Data is 2D -> 3D graph not done, please use plot_3d_slices instead.", title)
        return

    # Ensure inputs are numpy arrays
    x_vals = np.array(100*x_vals)
    y_vals = np.array(100*y_vals)
    z_vals = np.array(100*z_vals)

    # Validate dimensions
    if not (len(x_vals) == n_x and len(y_vals) == n_y and len(z_vals) == n_z):
        raise ValueError("Grid dimensions do not match data shape.")

    X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing="ij")

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    X_flat = X.flatten()
    Y_flat = Y.flatten()
    Z_flat = Z.flatten()
    data_flat = data.flatten()

    min_data = np.min(data_flat)
    max_data = np.max(data_flat)

    # Avoid division by zero if data is constant
    if max_data - min_data > 0:
        normalized_data = (data_flat - min_data) / (max_data - min_data)
    else:
        normalized_data = np.zeros_like(data_flat)

    colors = plt.cm.inferno(normalized_data)
    alpha_values = np.clip(normalized_data, 1e-5, 1)
    point_sizes = 100 * alpha_values
    colors[..., 3] = alpha_values

    ax.scatter(X_flat, Y_flat, Z_flat, c=colors, s=point_sizes, marker='o')

    norm = mcolors.Normalize(vmin=min_data, vmax=max_data)
    cbar = plt.colorbar(mappable=plt.cm.ScalarMappable(cmap="inferno", norm=norm), ax=ax)
    cbar.set_label(title, fontsize=15)
    cbar.ax.tick_params(labelsize=15)

    ax.tick_params(labelsize=14)
    ax.set_xlabel(xlabel, fontsize=16)
    ax.set_ylabel(ylabel, fontsize=16)
    ax.set_zlabel(zlabel, fontsize=16)

    plt.tight_layout()
    plt.savefig(output_image, dpi=300)
    plt.close()

def plot_3d_slices(data, output_prefix, title, x_vals, y_vals, z_vals, axis="z", color_map="inferno"):
    n_x, n_y, n_z = data.shape

    x_cm = 100 * x_vals
    y_cm = 100 * y_vals
    z_cm = 100 * z_vals

    if axis == "x":
        for i in range(n_x):
            if n_x > 10 and i % 8 and (i != n_x // 2):
                continue
            plt.figure(figsize=(8, 6))
            Y, Z = np.meshgrid(y_cm, z_cm, indexing="ij")
            
            
            #data_slice = np.ma.masked_where(data[i, :, :] <= 0, data[i, :, :])
            #if data_slice.count() == 0:
            #    plt.close()
            #    continue
            #norm = LogNorm(vmin=np.nanmin(data[data > 0.03]),vmax=np.nanmax(data))

            if(contour): mesh = plt.contourf(Y, Z, data[i, :, :], levels=100, cmap=color_map) # ,norm=norm)
            else: mesh = plt.pcolormesh(Y, Z, data[i, :, :], cmap=color_map, shading="auto")    #, vmin=0, vmax=0.5)  # ,norm=norm)
            cbar = plt.colorbar(mesh)
            cbar.set_label(f"{title} (Slice x={x_cm[i]:.2f} cm)" if n_x > 1 else title, fontsize=15)
            cbar.ax.tick_params(labelsize=16)
            plt.xlabel("Y [cm]", fontsize=18)
            plt.ylabel("Z [cm]", fontsize=18)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.savefig(f"{output_prefix}_slice_x_{i}.png", dpi=300, bbox_inches="tight")
            plt.savefig(f"{output_prefix}_slice_x_{i}.pdf", bbox_inches="tight")
            plt.close()

    elif axis == "y":
        for j in range(n_y):
            if n_y > 10 and j % 8 and (j != n_y // 2):
                continue
            plt.figure(figsize=(4, 6))
            X, Z = np.meshgrid(x_cm, z_cm, indexing="ij")
            if(contour): mesh = plt.contourf(X, Z, data[:, j, :], levels=100, cmap=color_map)
            else: mesh = plt.pcolormesh(X, Z, data[:, j, :], cmap=color_map, shading="auto")
            cbar = plt.colorbar(mesh)
            cbar.set_label(f"{title} (Slice y={y_cm[j]:.2f} cm)" if n_y > 1 else title, fontsize=18)
            cbar.ax.tick_params(labelsize=16)
            plt.xlabel("$x$ [cm]", fontsize=18)
            plt.ylabel("$z$ [cm]", fontsize=18)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.savefig(f"{output_prefix}_slice_y_{j}.pdf", dpi=300, bbox_inches="tight")
            plt.close()

    else:  # axis == "z"
        for k in range(n_z):
            if n_z > 10 and k % 8 and (k != n_z // 2):
                continue
            plt.figure(figsize=(4, 6))
            X, Y = np.meshgrid(x_cm, y_cm, indexing="ij")
            if(contour): mesh = plt.contourf(X, Y, data[:, :, k], levels=100, cmap=color_map)
            else: mesh = plt.pcolormesh(X, Y, data[:, :, k], cmap=color_map, shading="auto")
            cbar = plt.colorbar(mesh)
            cbar.set_label(f"{title} (Slice z={z_cm[k]:.2f} cm)" if n_z > 1 else title, fontsize=18)
            cbar.ax.tick_params(labelsize=16)
            plt.xlabel("$x$ [cm]", fontsize=18)
            plt.ylabel("$y$ [cm]", fontsize=18)
            plt.xticks(fontsize=16)
            plt.yticks(fontsize=16)
            plt.savefig(f"{output_prefix}_slice_z_{k}.png", dpi=300, bbox_inches="tight")
            plt.close()

def plot_3d_voxel(data, x_vals, y_vals, z_vals, output_image, title):
    n_x, n_y, n_z = data.shape

    if not (len(x_vals) == n_x and len(y_vals) == n_y and len(z_vals) == n_z):
        raise ValueError("Grid dimensions do not match data shape.")

    # Convert to numpy arrays and get voxel edges
    x_edges = 0.5 * (np.pad(x_vals, (1, 0), mode='edge') + np.pad(x_vals, (0, 1), mode='edge'))
    y_edges = 0.5 * (np.pad(y_vals, (1, 0), mode='edge') + np.pad(y_vals, (0, 1), mode='edge'))
    z_edges = 0.5 * (np.pad(z_vals, (1, 0), mode='edge') + np.pad(z_vals, (0, 1), mode='edge'))

    X, Y, Z = np.meshgrid(x_edges, y_edges, z_edges, indexing="ij")

    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Normalize data for colormap
    min_data, max_data = np.min(data), np.max(data)
    if max_data - min_data > 0:
        normalized_data = (data - min_data) / (max_data - min_data)
    else:
        normalized_data = np.zeros_like(data)

    colors = plt.cm.inferno(normalized_data)
    alpha_values = np.clip(normalized_data, 2e-1, 1)
    colors[..., 3] = alpha_values

    ax.voxels(X, Y, Z, data > (0.001 * np.max(data)), facecolors=colors, edgecolors='k', linewidth=0.1)

    ax.set_xlabel("X [cm]")
    ax.set_ylabel("Y [cm]")
    ax.set_zlabel("Z [cm]")
    ax.set_title(title)

    norm = mcolors.Normalize(vmin=min_data, vmax=max_data)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.inferno, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label='Value')

    plt.tight_layout()
    plt.savefig(output_image, dpi=300, bbox_inches="tight")
    plt.close()

def plot_1d_line(data, direction, idx_1, idx_2, output_prefix, title, grid_x, grid_y, grid_z):
    plt.figure(figsize=(10,6))        

    if direction == 'x':
        plt.xlabel("x [cm]", fontsize=16)
        plt.plot(grid_x, data[:, idx_1, idx_2])
    elif direction == 'y':
        plt.xlabel("y [cm]", fontsize=16)
        plt.plot(grid_y, data[idx_1, :, idx_2])
        print(max(data[idx_1, :40, idx_2]))
        print(min(data[idx_1, :40, idx_2]))
    elif direction == 'z':
        plt.xlabel("z [cm]", fontsize=16)
        plt.plot(grid_z, data[idx_1, idx_2, :])
    else:
        raise ValueError("Direction must be 'x', 'y', or 'z'.")

    plt.ylabel(title)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_1d.png", dpi=300)
    plt.close()
    print(f"Plot saved for {title}")

def smart_plot(data, output_prefix, title, grid_x, grid_y, grid_z, color_map="inferno"):
    if data is None:
        print("Data is None -> skipping plot.")
        return
    
    shape = data.shape
    dims = np.count_nonzero(np.array(shape) > 1)

    if dims == 3:
        #plot_3d_scatter(data, f"{output_prefix}_3d.png", title, grid_x, grid_y, grid_z)
        plot_3d_slices(data, output_prefix, title, grid_x, grid_y, grid_z, axis=axis, color_map=color_map)
        #plot_3d_voxel(data, grid_x, grid_y, grid_z, f"{output_prefix}_voxel_3d.png", title)
    elif dims == 2:
        plot_3d_slices(data, output_prefix, title, grid_x, grid_y, grid_z, axis=axis, color_map=color_map)
    elif dims == 1:
        print("Data is 1D -> using line plot instead.")
        nx, ny, nz = data.shape
        data = np.squeeze(data)  # Remove singleton dimensions
    
        grid_x = 100*grid_x; grid_y = 100*grid_y; grid_z = 100*grid_z

        plt.figure(figsize=(10,6))        
        
        if nx != 1:
            plt.xlabel("x [cm]", fontsize=16)
            plt.plot(grid_x, data)
        elif ny != 1:
            plt.xlabel("y [cm]", fontsize=16)
            plt.plot(grid_y, data)
        elif nz != 1:
            plt.xlabel("z [cm]", fontsize=16)
            plt.plot(grid_z, data)

        plt.ylabel(title)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_1d.png", dpi=300)
        plt.close()
    else:
        print("Data has invalid shape:", shape)

    print(f"Plot saved for {title}")


parser = argparse.ArgumentParser(description="3D Graph Plotter")
parser.add_argument('axis', nargs='?', default='z', choices=['x', 'y', 'z'], help="Axis to plot along (default: 'z')")
args = parser.parse_args()
axis = args.axis

# Usage
graph_outputs = True
graph_P_abs_i_ratios = False
graph_f_grad = False
graph_differences = False
graph_conductivity = True
graph_components = True

data, grid_x, grid_y, grid_z = load_3d_data(folder + "absorbed_power_density.txt")
smart_plot(data/1e6, "analysis/analysis_output/absorbed_power_density", "$P_{abs}$ [W / cm$^{3}$]", grid_x, grid_y, grid_z, color_map="inferno")

if(graph_conductivity):
    ne, grid_x, grid_y, grid_z = load_3d_data(folder + "electron_density.txt")
    mu_re, _, _, _ = load_3d_data(folder + "real_mobility.txt")

    e_charge = 1.60217663e-19
    cond = ne*mu_re*e_charge

    smart_plot(cond, "analysis/analysis_output/Conductivity", "$\\sigma_{Re}$", grid_x, grid_y, grid_z)

    data, grid_x, grid_y, grid_z = load_3d_data(folder + "real_mobility.txt")
    smart_plot(data, "analysis/analysis_output/real_mobility", "$\\mu_{Re}$ [m$^{2}$/(V⋅s)]", grid_x, grid_y, grid_z, color_map="viridis")

    data, grid_x, grid_y, grid_z = load_3d_data(folder + "imag_mobility.txt")
    smart_plot(data, "analysis/analysis_output/imag_mobility", "$\\mu_{Im}$ [m$^{2}$/(V⋅s)]", grid_x, grid_y, grid_z, color_map="viridis")

    data, grid_x, grid_y, grid_z = load_3d_data(folder + "electron_density.txt")
    smart_plot(data, "analysis/analysis_output/electron_density", "n$_e$ [m$^{-3}$]", grid_x, grid_y, grid_z, color_map="viridis")

if(graph_outputs):
    if(graph_components):
        data, grid_x, grid_y, grid_z = load_3d_data(folder + "Ez_amp.txt")
        smart_plot(data/1e3, "analysis/analysis_output/Ez_amp", "$E_{z,amp}$ [kV/m]", grid_x, grid_y, grid_z, color_map="plasma")

        data, grid_x, grid_y, grid_z = load_3d_data(folder + "Ex_amp.txt")
        smart_plot(data/1e3, "analysis/analysis_output/Ex_amp", "$E_{x,amp}$ [kV/m]", grid_x, grid_y, grid_z, color_map="plasma")

        data, grid_x, grid_y, grid_z = load_3d_data(folder + "Ey_amp.txt")
        smart_plot(data/1e3, "analysis/analysis_output/Ey_amp", "$E_{y,amp}$ [kV/m]", grid_x, grid_y, grid_z, color_map="plasma")

    data, grid_x, grid_y, grid_z = load_3d_data(folder + "E_amplitude.txt")
    smart_plot(data/1e3, "analysis/analysis_output/E_amplitude", "$E_{amp}$ [kV/m]", grid_x, grid_y, grid_z, color_map="plasma")

    data, grid_x, grid_y, grid_z = load_3d_data(folder + "power_flux_amp.txt")
    smart_plot(data, "analysis/analysis_output/power_flux_amp", "$J_{power}$ [W /m$^{2}$]", grid_x, grid_y, grid_z, color_map="inferno")

    data, grid_x, grid_y, grid_z = load_3d_data(folder + "power_poynting.txt")
    smart_plot(data, "analysis/analysis_output/power_poynting", "$\\nabla S$ [W /m$^{3}$]", grid_x, grid_y, grid_z, color_map="inferno")

    data, grid_x, grid_y, grid_z = load_3d_data(folder + "power_flux_x.txt")
    smart_plot(data, "analysis/analysis_output/power_flux_x", "$J_{power,x}$ [W /m$^{2}$]", grid_x, grid_y, grid_z, color_map="inferno")
    data, grid_x, grid_y, grid_z = load_3d_data(folder + "power_flux_y.txt")
    smart_plot(data, "analysis/analysis_output/power_flux_y", "$J_{power,y}$ [W /m$^{2}$]", grid_x, grid_y, grid_z, color_map="inferno")
    data, grid_x, grid_y, grid_z = load_3d_data(folder + "power_flux_z.txt")
    smart_plot(data, "analysis/analysis_output/power_flux_z", "$J_{power,z}$ [W /m$^{2}$]", grid_x, grid_y, grid_z, color_map="inferno")


    data, grid_x, grid_y, grid_z = load_3d_data(folder + "E_real.txt")
    if data: smart_plot(data/1e3, "analysis/analysis_output/E_real", "$E_{real}$ [kV/m]", grid_x, grid_y, grid_z, color_map="plasma")

    data, grid_x, grid_y, grid_z = load_3d_data(folder + "E_imag.txt")
    if data: smart_plot(data/1e3, "analysis/analysis_output/E_imag", "$E_{imag}$ [kV/m]", grid_x, grid_y, grid_z, color_map="plasma")


if(graph_P_abs_i_ratios):
    E_amp, grid_x, grid_y, grid_z = load_3d_data(folder + "E_amplitude.txt")
    Ex_amp, _, _, _ = load_3d_data(folder + "Ex_amp.txt")
    Ey_amp, _, _, _ = load_3d_data(folder + "Ey_amp.txt")
    Ez_amp, _, _, _  = load_3d_data(folder + "Ez_amp.txt")

    # Avoid divide-by-zero by adding a small epsilon
    epsilon = 1e-18
    E_amp_sq = E_amp**2 + epsilon

    Ex_ratio = (Ex_amp**2) / E_amp_sq
    Ey_ratio = (Ey_amp**2) / E_amp_sq
    Ez_ratio = (Ez_amp**2) / E_amp_sq

    # Plot the ratios
    smart_plot(Ex_ratio, "analysis/analysis_output/Ex_ratio", "$|E_{x}|^2 / |E|^2$", grid_x, grid_y, grid_z)
    smart_plot(Ey_ratio, "analysis/analysis_output/Ey_ratio", "$|E_{y}|^2 / |E|^2$", grid_x, grid_y, grid_z)
    smart_plot(Ez_ratio, "analysis/analysis_output/Ez_ratio", "$|E_{z}|^2 / |E|^2$", grid_x, grid_y, grid_z)

if(graph_f_grad):
    data, grid_x, grid_y, grid_z = load_3d_data(folder + "f_grad_cond_x.txt")
    smart_plot(data, "analysis/analysis_output/f_grad_cond_x", "$\\nabla f_{cond} {x,amp}$ [$m^{-1}$]", grid_x, grid_y, grid_z)

    data, grid_x, grid_y, grid_z = load_3d_data(folder + "f_grad_cond_y.txt")
    smart_plot(data, "analysis/analysis_output/f_grad_cond_y", "$\\nabla f_{cond} {y,amp}$ [$m^{-1}$]", grid_x, grid_y, grid_z)

    data, grid_x, grid_y, grid_z = load_3d_data(folder + "f_grad_cond_z.txt")
    smart_plot(data, "analysis/analysis_output/f_grad_cond_z", "$\\nabla f_{cond} {z,amp}$ [$m^{-1}$]", grid_x, grid_y, grid_z)

if(graph_differences):
    P_abs_no_grad, grid_x, grid_y, grid_z = load_3d_data(folder + "absorbed_power_density.txt")
    P_abs_with_grad, _, _, _ = load_3d_data(folder + "absorbed_power_density_with_grad.txt")

    P_abs_difference = P_abs_with_grad - P_abs_no_grad
    P_abs_relative_difference = P_abs_difference / (P_abs_no_grad + 1)

    smart_plot(P_abs_difference, "analysis/analysis_output/P_abs_difference", "$P_{abs, with grad} - P_{abs, no grad}$ [W/m$^{3}$]", grid_x, grid_y, grid_z)
    smart_plot(P_abs_relative_difference, "analysis/analysis_output/P_abs_relative_difference", "$\\frac{P_{abs, with grad} - P_{abs, no grad}}{P_{abs, no grad}}$", grid_x, grid_y, grid_z)



print("3D Visualizations saved successfully in the analysis_output folder.")
