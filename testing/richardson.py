import numpy as np
from scipy.interpolate import RegularGridInterpolator

contour = True

def plot_3d_slices(data, output_prefix, title, x_vals, y_vals, z_vals, axis="z"):
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
            if(contour): mesh = plt.contourf(Y, Z, data[i, :, :], levels=100, cmap="inferno")
            else: mesh = plt.pcolormesh(Y, Z, data[i, :, :], cmap="inferno", shading="auto")
            cbar = plt.colorbar(mesh)
            cbar.set_label(f"{title} (Slice x={x_cm[i]:.2f} cm)" if n_x > 1 else title, fontsize=15)
            cbar.ax.tick_params(labelsize=15)
            plt.xlabel("Y [cm]", fontsize=16)
            plt.ylabel("Z [cm]", fontsize=16)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.savefig(f"{output_prefix}_slice_x_{i}.png", dpi=300, bbox_inches="tight")
            plt.close()

    elif axis == "y":
        for j in range(n_y):
            if n_y > 10 and j % 8 and (j != n_y // 2):
                continue
            plt.figure(figsize=(8, 6))
            X, Z = np.meshgrid(x_cm, z_cm, indexing="ij")
            if(contour): mesh = plt.contourf(X, Z, data[:, j, :], levels=100, cmap="inferno")
            else: mesh = plt.pcolormesh(X, Z, data[:, j, :], cmap="inferno", shading="auto")
            cbar = plt.colorbar(mesh)
            cbar.set_label(f"{title} (Slice y={y_cm[j]:.2f} cm)" if n_y > 1 else title, fontsize=15)
            cbar.ax.tick_params(labelsize=15)
            plt.xlabel("X [cm]", fontsize=16)
            plt.ylabel("Z [cm]", fontsize=16)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.savefig(f"{output_prefix}_slice_y_{j}.png", dpi=300, bbox_inches="tight")
            plt.close()

    else:  # axis == "z"
        for k in range(n_z):
            if n_z > 10 and k % 8 and (k != n_z // 2):
                continue
            plt.figure(figsize=(8, 6))
            X, Y = np.meshgrid(x_cm, y_cm, indexing="ij")
            if(contour): mesh = plt.contourf(X, Y, data[:, :, k], levels=100, cmap="inferno")
            else: mesh = plt.pcolormesh(X, Y, data[:, :, k], cmap="inferno", shading="auto")
            cbar = plt.colorbar(mesh)
            cbar.set_label(f"{title} (Slice z={z_cm[k]:.2f} cm)" if n_z > 1 else title, fontsize=15)
            cbar.ax.tick_params(labelsize=15)
            plt.xlabel("X [cm]", fontsize=16)
            plt.ylabel("Y [cm]", fontsize=16)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.savefig(f"{output_prefix}_slice_z_{k}.png", dpi=300, bbox_inches="tight")
            plt.close()


def richardson_extrapolation(fine, coarse, order, scale):
    factor = scale ** order
    return (factor*fine - coarse)/(factor - 1)

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

                return data, x_coords, y_coords, z_coords
            
            else:
                # Uniform grid mode
                x_coords = np.linspace(0, length_x, n_x)
                y_coords = np.linspace(0, length_y, n_y)
                z_coords = np.linspace(0, length_z, n_z)

                data = np.loadtxt(f)
                data = data.reshape((n_x, n_y, n_z))

                return data, x_coords, y_coords, z_coords
        
    except FileNotFoundError:
        print(f"Warning: File '{filename}' not found. Skipping.")
        return None, None, None, None

def interp_to_target(source_data, source_x, source_y, source_z, target_x, target_y, target_z):
    interp = RegularGridInterpolator((source_x, source_y, source_z), source_data, bounds_error=False, fill_value=None)

    X, Y, Z = np.meshgrid(target_x, target_y, target_z, indexing='ij')
    points = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=-1)

    return  interp(points).reshape(X.shape)

data_coarse, grid_course_x, grid_course_y, grid_course_z = load_3d_data("Outputs/E_amplitude_512.txt")
data_fine, grid_fine_x, grid_fine_y, grid_fine_z = load_3d_data("Outputs/E_amplitude_1024.txt")

order = 2
scale = 2

data_coarse_on_fine = interp_to_target(data_coarse, grid_course_x, grid_course_y, grid_course_z, grid_fine_x, grid_fine_y, grid_fine_z)
data_richardson = richardson_extrapolation(data_fine, data_coarse_on_fine, order, scale)


# Comparing
data_coarse, grid_course_x, grid_course_y, grid_course_z = load_3d_data("Outputs/E_amplitude_512.txt")
data_coarse_on_fine = interp_to_target(data_coarse, grid_course_x, grid_course_y, grid_course_z, grid_fine_x, grid_fine_y, grid_fine_z)

print("Max values:")
max_coarse = np.max(data_coarse)
max_fine = np.max(data_fine)
max_richardson = np.max(data_richardson)

print(max_coarse, max_fine, max_richardson)

sum_richardson = np.sum(data_richardson)
err_coarse = np.linalg.norm(data_richardson - data_coarse_on_fine)/np.linalg.norm(data_richardson)
err_fine   = np.linalg.norm(data_richardson - data_fine)/np.linalg.norm(data_richardson)


print("Error estimates (sum):")
print("Coarse:", err_coarse)
print("Fine:", err_fine)

# plot diferences
difference_coarse = (data_richardson - data_coarse_on_fine)/np.max(data_richardson)
difference_fine   = (data_richardson - data_fine)/np.max(data_richardson)

import matplotlib.pyplot as plt

plot_3d_slices(difference_coarse, "testing/richardson_difference_coarse.png", "Differences to Richardson", grid_fine_x, grid_fine_y, grid_fine_z, axis="z")
plot_3d_slices(difference_fine, "testing/richardson_difference_fine.png", "Differences to Richardson", grid_fine_x, grid_fine_y, grid_fine_z, axis="z")

