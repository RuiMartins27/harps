
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

FLAG_CENTER_COORDINATES = True
CENTER_X, CENTER_Y, CENTER_Z = 0.0000, 0.056, 0.0

plt.rcParams.update({"font.family": "serif", "font.size": 11, "mathtext.fontset": "stix"})


def load_3d_data(filepath: Path, center_coords: bool = FLAG_CENTER_COORDINATES, center_xyz: tuple = (CENTER_X, CENTER_Y, CENTER_Z)):
    filepath = Path(filepath)
    if not filepath.is_file():
        print(f"Error: File '{filepath}' not found.")
        return None, None, None, None
    with open(filepath, "r") as f:
        n_x, n_y, n_z, len_x, len_y, len_z = map(float, f.readline().split())
        n_x, n_y, n_z = int(n_x), int(n_y), int(n_z)
        header = f.readline().strip()
        if header == "NON_UNIFORM_GRID":
            x_coords = np.fromstring(f.readline(), sep=" ")
            y_coords = np.fromstring(f.readline(), sep=" ")
            z_coords = np.fromstring(f.readline(), sep=" ")
        else:
            x_coords, y_coords, z_coords = np.linspace(0, len_x, n_x), np.linspace(0, len_y, n_y), np.linspace(0, len_z, n_z)
        data = np.loadtxt(f).reshape((n_x, n_y, n_z))
    if center_coords:
        x_coords, y_coords, z_coords = x_coords - center_xyz[0], y_coords - center_xyz[1], z_coords - center_xyz[2]
    if header != "NON_UNIFORM_GRID":
        z_mask = z_coords <= 0.5
        z_coords, data = z_coords[z_mask], data[:, :, z_mask]
    return data, x_coords, y_coords, z_coords


def azimuthal_average(data, x_coords, y_coords, z_coords, axis_x=CENTER_X, axis_y=CENTER_Y, n_r=100, r_max=0.0135, verbose=False):
    n_x, n_y, n_z = data.shape
    X, Y = np.meshgrid(x_coords, y_coords, indexing="ij")
    R = np.hypot(X - axis_x, Y - axis_y)
    if verbose:
        r_full = min(x_coords.max() - axis_x, axis_x - x_coords.min(), y_coords.max() - axis_y, axis_y - y_coords.min())
        print(f"Full 360° coverage up to r = {r_full:.4g} m.")
    r_max = r_max if r_max is not None else R.max()
    n_r = n_r if n_r is not None else max(n_x, n_y) * 2
    r_edges = np.linspace(0.0, r_max, n_r + 1)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    bin_idx = np.digitize(R.ravel(), r_edges) - 1
    valid = (bin_idx >= 0) & (bin_idx < n_r)
    bin_idx = bin_idx[valid]
    counts = np.bincount(bin_idx, minlength=n_r).astype(float)
    counts[counts == 0] = np.nan
    data_flat, data_rz = data.reshape(n_x * n_y, n_z)[valid, :], np.empty((n_r, n_z))
    for k in range(n_z):
        sums = np.bincount(bin_idx, weights=data_flat[:, k], minlength=n_r)
        data_rz[:, k] = sums / counts
    return r_centers, z_coords, data_rz


def plot_rz(r, z, data_rz, output_path: Path, colorbar_label: str, unit: str = "mm", color_map: str = "inferno", log_scale: bool = False, mirror: bool = False, contour_lines: bool = False, aspect: str = "auto", levels: int = 100):
    unit_scale = 1e3 if unit == "mm" else (1e2 if unit == "cm" else 1.0)
    unit_label = f"[{unit}]" if unit != "m" else "[m]"
    r_scaled, z_scaled = r * unit_scale, z * (1e2 if unit in ["mm", "cm"] else 1.0)
    z_unit_label = "[cm]" if unit in ["mm", "cm"] else "[m]"
    if mirror:
        r_plot, data_plot = np.concatenate([-r_scaled[::-1], r_scaled]), np.concatenate([data_rz[::-1, :], data_rz], axis=0)
    else:
        r_plot, data_plot = r_scaled, data_rz
    Z, R = np.meshgrid(z_scaled, r_plot)
    fig, ax = plt.subplots(figsize=(8, 3), dpi=300)
    norm = LogNorm(vmin=np.nanmin(data_plot[data_plot > 0]), vmax=np.nanmax(data_plot[data_plot > 0])) if log_scale else None

    # Filled contour plot for smooth shading
    cntr = ax.contourf(Z, R, data_plot, levels=levels, cmap=color_map, norm=norm)
    if contour_lines:
        ax.contour(Z, R, data_plot, levels=10, colors="white", linewidths=0.5, alpha=0.6)

    ax.set_xlabel(f"Axial position $z$ {z_unit_label}")
    ax.set_ylabel(f"Radial position $r$ {unit_label}")
    ax.set_aspect(aspect)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.12)
    cbar = fig.colorbar(cntr, cax=cax)
    cbar.set_label(colorbar_label)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    output_file = output_path.with_suffix(".png")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Plot saved successfully to: {output_file}")
    plt.close(fig)


def integrate_total_power_rz(r, z, data_rz):
    radial_integral = np.trapezoid(data_rz * (2.0 * np.pi * r[:, None]), r, axis=0)
    return np.trapezoid(radial_integral, z, axis=0)


def main():
    folder, output = Path("Outputs"), Path("analysis/analysis_output/absorbed_power_density_rz")
    n_r, r_max, unit = 20, 0.013, "mm"
    log_scale, mirror, contour_lines, aspect = False, False, False, "auto"
    input_file = folder / "absorbed_power_density_2D_600W.txt"

    data, grid_x, grid_y, grid_z = load_3d_data(input_file)
    if data is None:
        return

    axis_x, axis_y = (0.0, 0.0) if FLAG_CENTER_COORDINATES else (CENTER_X, CENTER_Y)
    r, z, data_rz = azimuthal_average(data, grid_x, grid_y, grid_z, axis_x=axis_x, axis_y=axis_y, n_r=n_r, r_max=r_max)
    total_power = integrate_total_power_rz(r, z, data_rz)
    print(f"Total volume-integrated absorbed power: {total_power:.6g} W")

    plot_rz(r, z, data_rz / 1e6, output_path=output, colorbar_label=r"$P_{\mathrm{abs}}\quad [\mathrm{W / cm}^3]$", unit=unit, color_map="inferno", log_scale=log_scale, mirror=mirror, contour_lines=contour_lines, aspect=aspect)


if __name__ == "__main__":
    main()