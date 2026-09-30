from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "mathtext.fontset": "stix",
    "axes.labelsize": 14, "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13,
    "legend.title_fontsize": 14, "axes.linewidth": 1.0, "grid.linewidth": 0.5, "lines.linewidth": 2.0,
})

CENTER_Y = 0.056
CENTER_X_3D = 0.04318
CENTER_X_2D = 0.0
FOLDER = Path("Outputs")

FILE_PABS_3D, FILE_PABS_2D = FOLDER / "absorbed_power_density_3D_600W.txt", FOLDER / "absorbed_power_density_2D_600W.txt"
FILE_E_3D, FILE_E_2D = FOLDER / "E_amplitude_3D_600W.txt", FOLDER / "E_amplitude_2D_600W.txt"

def load_3d_data(filepath: Path):
    if not (filepath := Path(filepath)).is_file(): return None, None, None, None
    with open(filepath, "r") as f:
        n_x, n_y, n_z, len_x, len_y, len_z = map(float, f.readline().split())
        header = f.readline().strip()
        if header == "NON_UNIFORM_GRID":
            x_coords = np.fromstring(f.readline(), sep=" ")
            y_coords = np.fromstring(f.readline(), sep=" ")
            z_coords = np.fromstring(f.readline(), sep=" ")
        else:
            x_coords, y_coords, z_coords = np.linspace(0, len_x, int(n_x)), np.linspace(0, len_y, int(n_y)), np.linspace(0, len_z, int(n_z))
        data = np.loadtxt(f).reshape((int(n_x), int(n_y), int(n_z)))
    if header != "NON_UNIFORM_GRID":
        mask = z_coords <= 0.5
        z_coords, data = z_coords[mask], data[:, :, mask]
    return data, x_coords, y_coords, z_coords

def azimuthal_average_3d(data, x_coords, y_coords, z_coords, axis_x, axis_y, n_r=20, r_max=0.0135):
    R = np.hypot(*np.meshgrid(x_coords - axis_x, y_coords - axis_y, indexing="ij"))
    r_edges = np.linspace(0.0, r_max, n_r + 1)
    bin_idx = np.digitize(R.ravel(), r_edges) - 1
    valid = (bin_idx >= 0) & (bin_idx < n_r)
    counts = np.bincount(bin_idx[valid], minlength=n_r).astype(float)
    counts[counts == 0] = np.nan
    data_flat = data.reshape(data.shape[0] * data.shape[1], data.shape[2])[valid, :]
    data_rz = np.empty((n_r, data.shape[2]))
    for k in range(data.shape[2]): data_rz[:, k] = np.bincount(bin_idx[valid], weights=data_flat[:, k], minlength=n_r) / counts
    return 0.5 * (r_edges[:-1] + r_edges[1:]), z_coords, data_rz

def azimuthal_average_2d(data, y_coords, z_coords, n_r=20, r_max=0.0135):
    r_raw = np.abs(y_coords - CENTER_Y)
    r_edges = np.linspace(0.0, r_max, n_r + 1)
    bin_idx = np.digitize(r_raw, r_edges) - 1
    data_rz, p_yz = np.empty((n_r, data.shape[2])), data[0, :, :]
    for k in range(n_r):
        mask = bin_idx == k
        data_rz[k, :] = np.nanmean(p_yz[mask, :], axis=0) if np.any(mask) else np.nan
    return 0.5 * (r_edges[:-1] + r_edges[1:]), z_coords, data_rz

def compute_total_power_3d(data, x, y, z):
    return np.trapezoid(np.trapezoid(np.trapezoid(data, x=z, axis=2), x=y, axis=1), x=x, axis=0)

def compute_total_power_2d(data, x, y, z):
    if len(x) != 1: raise ValueError(f"Expected a single x-plane for 2D data, got len(x)={len(x)}")
    p_yz, r = data[0, :, :], np.abs(y - CENTER_Y)
    order = np.argsort(r)
    r, p_yz = r[order], p_yz[order, :]
    return np.trapezoid(np.trapezoid(2.0 * np.pi * r[:, None] * p_yz, x=z, axis=1), x=r, axis=0)

def process_and_plot_data(datasets, y_label, y_scale_factor, filename, z_slices, rescale_power=False, calc_global_l2=False):
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    z_colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(z_slices)))
    rz_results = {}

    for filepath, label, ls, center_x in datasets:
        data, x, y, z = load_3d_data(filepath)
        if data is None:
            print(f"Warning: {filepath} not found — skipping {label}.")
            continue
        if rescale_power:
            p_pabs_3d, x_p, y_p, z_p = load_3d_data(FILE_PABS_3D if label == "3D" else FILE_PABS_2D)
            total = compute_total_power_3d(p_pabs_3d, x_p, y_p, z_p) if label == "3D" else compute_total_power_2d(p_pabs_3d, x_p, y_p, z_p)
            print(f"{label} (Power Normalization): integrated power = {total:.6f} W")
            if not np.isfinite(total) or total <= 0: continue
            data = data * (600.0 / total)

        if label == "3D":
            r, z_coords, data_rz = azimuthal_average_3d(data, x, y, z, axis_x=center_x, axis_y=CENTER_Y, n_r=20, r_max=0.0135)
        else:
            r, z_coords, data_rz = azimuthal_average_2d(data, y, z, n_r=20, r_max=0.0135)
        rz_results[label] = (r, z_coords, data_rz)

        for z_cm, color in zip(z_slices, z_colors):
            idx = np.argmin(np.abs(z_coords - z_cm / 100.0))
            ax.plot(r * 1e2, data_rz[:, idx] / y_scale_factor, ls, color=color, label="_nolegend_")

    if not rz_results: return

    if "3D" in rz_results and "2D" in rz_results:
        r_3d, z_3d, rz_3d = rz_results["3D"]
        r_2d, z_2d, rz_2d = rz_results["2D"]

        if calc_global_l2:
            z_min, z_max = 0.005, 0.095
            mask_3d = (z_3d >= z_min) & (z_3d <= z_max)
            z_eval = z_3d[mask_3d]
            rz_3d_domain = rz_3d[:, mask_3d]
            rz_2d_interp = np.zeros_like(rz_3d_domain)
            for i in range(rz_3d.shape[0]): rz_2d_interp[i, :] = np.interp(z_eval, z_2d, rz_2d[i, :])
            l2_abs_domain = np.linalg.norm(rz_3d_domain - rz_2d_interp)
            l2_rel_domain = l2_abs_domain / np.linalg.norm(rz_3d_domain)
            print(f"\nGeneral Domain (r < 1.35 cm, 0.5 cm <= z <= 9.5 cm) L2 Error:\n  Absolute:  {l2_abs_domain:.6e}\n  Relative:  {l2_rel_domain:.6e} ({l2_rel_domain*100:.2f}%)")

        print("\nSlice-wise L2 Errors:")
        for z_cm in z_slices:
            idx_3d = np.argmin(np.abs(z_3d - z_cm / 100.0))
            idx_2d = np.argmin(np.abs(z_2d - z_cm / 100.0))
            l2_abs = np.linalg.norm(rz_3d[:, idx_3d] - rz_2d[:, idx_2d])
            l2_rel = l2_abs / np.linalg.norm(rz_3d[:, idx_3d])
            print(f"  z = {z_cm} cm | Absolute: {l2_abs:.6e} | Relative: {l2_rel:.6e} ({l2_rel*100:.2f}%)")

    ax.set_xlabel(r"Radial position $r$ [cm]")
    ax.set_ylabel(y_label)
    ax.set_xlim(0, 1.35)
    ax.set_ylim(bottom=0)
    ax.tick_params(direction="in", top=True, right=True, which="both")
    ax.grid(True, alpha=0.5, linestyle=":")

    h1 = [Line2D([0], [0], color="k", linestyle=ls, label=l) for _, l, ls, _ in datasets]
    legend1 = ax.legend(handles=h1, loc="upper right", frameon=True, framealpha=0.9, facecolor="white", edgecolor="none", title="EM Solver")
    ax.add_artist(legend1)

    h2 = [Line2D([0], [0], color=c, linestyle="-", label=f"$z = {z_cm}$ cm") for z_cm, c in zip(z_slices, z_colors)]
    ax.legend(handles=h2, loc="right", bbox_to_anchor=[1, 0.59], frameon=True, framealpha=0.9, facecolor="white", edgecolor="none", title="Axial slice")

    output_dir = Path("analysis/analysis_output")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / filename
    fig.savefig(output_file.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_file.with_suffix(".pdf"), dpi=300, bbox_inches="tight")
    plt.close(fig)

def main():
    z_slices = [3.5, 5]
    
    print("=== Processing Absorbed Power Density ===")
    pabs_datasets = [
        (FILE_PABS_3D, "3D", "-", CENTER_X_3D),
        (FILE_PABS_2D, "2D", "--", CENTER_X_2D),
    ]
    process_and_plot_data(
        datasets=pabs_datasets,
        y_label=r"$P_{\mathrm{abs}}\quad [\mathrm{W/cm}^3]$",
        y_scale_factor=1e6,
        filename="figure4_power_radial_3D_2D",
        z_slices=z_slices,
        rescale_power=True,
        calc_global_l2=True
    )

    print("\n=== Processing Electric Field Amplitude ===")
    e_datasets = [
        (FILE_E_3D, "3D", "-", CENTER_X_3D),
        (FILE_E_2D, "2D", "--", CENTER_X_2D),
    ]
    process_and_plot_data(
        datasets=e_datasets,
        y_label=r"$|E|\quad [\mathrm{V/m}]$",
        y_scale_factor=1.0,
        filename="figure4_E_field_radial_3D_2D",
        z_slices=z_slices,
        rescale_power=False,
        calc_global_l2=False
    )

if __name__ == "__main__":
    main()