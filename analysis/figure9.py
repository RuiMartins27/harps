from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

FLAG_CENTER_COORDINATES = True
CENTER_X, CENTER_Y, CENTER_Z = 0.0000, 0.056, 0.0

plt.rcParams.update({"font.family": "serif", "font.size": 12, "mathtext.fontset": "stix"})


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
    r_max = r_max if r_max is not None else R.max()
    n_r = n_r if n_r is not None else max(n_x, n_y) * 2
    r_edges = np.linspace(0.0, r_max, n_r + 1)
    bin_idx = np.digitize(R.ravel(), r_edges) - 1
    valid = (bin_idx >= 0) & (bin_idx < n_r)
    bin_idx = bin_idx[valid]
    counts = np.bincount(bin_idx, minlength=n_r).astype(float)
    counts[counts == 0] = np.nan
    data_flat, data_rz = data.reshape(n_x * n_y, n_z)[valid, :], np.empty((n_r, n_z))
    for k in range(n_z):
        sums = np.bincount(bin_idx, weights=data_flat[:, k], minlength=n_r)
        data_rz[:, k] = sums / counts
    return 0.5 * (r_edges[:-1] + r_edges[1:]), z_coords, data_rz


def plot_comparison_rz(r1, z1, data_rz1, r2, z2, data_rz2, output_path: Path, colorbar_label: str, unit: str = "mm", color_map: str = "inferno", log_scale: bool = False, aspect: str = "auto", levels: int = 100):
    unit_scale = 1e3 if unit == "mm" else (1e2 if unit == "cm" else 1.0)
    unit_label = f"[{unit}]" if unit != "m" else "[m]"
    z_unit_label = "[cm]" if unit in ["mm", "cm"] else "[m]"
    r_s1, z_s1 = r1 * unit_scale, z1 * (1e2 if unit in ["mm", "cm"] else 1.0)
    r_s2, z_s2 = r2 * unit_scale, z2 * (1e2 if unit in ["mm", "cm"] else 1.0)
    Z1, R1 = np.meshgrid(z_s1, r_s1)
    Z2, R2 = np.meshgrid(z_s2, r_s2)
    fig, axes = plt.subplots(2, 1, figsize=(8, 5), dpi=300, sharex=True)
    all_data = np.concatenate([data_rz1, data_rz2], axis=None)
    norm = LogNorm(vmin=np.nanmin(all_data[all_data > 0]), vmax=np.nanmax(all_data[all_data > 0])) if log_scale else None
    vmin, vmax = np.nanmin(all_data), np.nanmax(all_data)
    
    cntr0 = axes[0].contourf(Z1, R1, data_rz1, levels=levels, cmap=color_map, norm=norm, vmin=None if log_scale else vmin, vmax=None if log_scale else vmax)
    cntr1 = axes[1].contourf(Z2, R2, data_rz2, levels=levels, cmap=color_map, norm=norm, vmin=None if log_scale else vmin, vmax=None if log_scale else vmax)
    
    axes[0].set_ylabel(f"Radial $r$ {unit_label}")
    axes[1].set_ylabel(f"Radial $r$ {unit_label}")
    axes[1].set_xlabel(f"Axial position $z$ {z_unit_label}")
    axes[0].tick_params(labelbottom=False)

    # Add panel labels (a) and (b) on the top left of each plot
    for ax, label in zip(axes, ["(a)", "(b)"]):
        ax.text(0.02, 0.95, label, transform=ax.transAxes, fontsize=14, fontweight="bold",
                verticalalignment="top", bbox=dict(boxstyle="square,pad=0.2", facecolor="white", alpha=0.95, edgecolor="none"))

    for ax in axes:
        ax.set_aspect(aspect)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        
    fig.subplots_adjust(right=0.82)
    cax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    cbar = fig.colorbar(cntr1, cax=cax)
    cbar.set_label(colorbar_label)
    
    output_file = output_path.with_suffix(".png")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    print(f"Plot saved successfully to: {output_file}")
    plt.close(fig)


def main():
    folder, output = Path("Outputs"), Path("analysis/analysis_output/figure9_grad_power")
    n_r, r_max, unit = 12, 0.013, "mm"
    file_no_grad = folder / "absorbed_power_density_2D_no_grad_600W.txt"
    file_grad = folder / "absorbed_power_density_2D_600W.txt"
    data1, gx1, gy1, gz1 = load_3d_data(file_no_grad)
    data2, gx2, gy2, gz2 = load_3d_data(file_grad)
    if data1 is None or data2 is None:
        return
    axis_x, axis_y = (0.0, 0.0) if FLAG_CENTER_COORDINATES else (CENTER_X, CENTER_Y)
    r1, z1, data_rz1 = azimuthal_average(data1, gx1, gy1, gz1, axis_x=axis_x, axis_y=axis_y, n_r=n_r, r_max=r_max)
    r2, z2, data_rz2 = azimuthal_average(data2, gx2, gy2, gz2, axis_x=axis_x, axis_y=axis_y, n_r=n_r, r_max=r_max)
    plot_comparison_rz(r1, z1, data_rz1 / 1e6, r2, z2, data_rz2 / 1e6, output_path=output, colorbar_label=r"$P_{\mathrm{abs}}\quad [\mathrm{W / cm}^3]$", unit=unit)


if __name__ == "__main__":
    main()