import numpy as np
import matplotlib.pyplot as plt

# --- 1. Paper-Ready Styling ---
plt.rcParams.update({
    "text.usetex": False, # Set to True if you have LaTeX (e.g., MiKTeX/TeX Live)
    "font.family": "serif",
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 13.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 2.4,
    "legend.frameon": True,
    "legend.edgecolor": "black",
    "legend.fancybox": False
})

# --- 2. Data Organization ---
# Replace these with your actual file paths
regimes = {
    "Weak": ["output_harps_weak.txt", "output_WKB_weak.txt", "output_constant_E.txt"],
    "Average": ["output_harps_avg.txt", "output_WKB_avg.txt", "output_constant_E.txt"],
    "Strong": ["output_harps_strong.txt", "output_WKB_strong.txt", "output_constant_E.txt"]
}

model_labels = ["2D Maxwell", "WKB", "E Constant"]
regime_keys = ["Weak", "Average", "Strong"]
L_z = 0.04
target_power = 770 * 1e-6 # MW

# --- 3. Plotting Function ---
def plot_comparison(metric_type='E'):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.1), sharex=True, constrained_layout=True)
    
    for idx, regime in enumerate(regime_keys):
        ax = axes[idx]
        files = regimes[regime]
        
        # We need a reference for E-field normalization per regime
        # Usually HARPS is the first file in the list
        ref_data = np.loadtxt(files[0], delimiter=",")
        E_ref_max = np.max(ref_data[:, 1])
        
        for i, f in enumerate(files):
            # Load: r=0, E=1, p=2
            arr = np.loadtxt(f, delimiter=",")
            r = arr[:, 0] * 100 # Convert to cm
            
            if metric_type == 'E':
                val = arr[:, 1] * 1e-3 # Convert to kV/m
                # Optional: Scale to HARPS max as per your original logic
                val *= (E_ref_max * 1e-3) / np.max(arr[:, 1] * 1e-3)
                ylabel = "Electric Field Amplitude [kV/m]"
                filename = "Figure_E_Field_Comparison"
            else:
                p = arr[:, 2]
                # Cylindrical integral normalization
                integral = np.trapezoid(L_z * 2 * np.pi * p * (arr[:, 0]), arr[:, 0])
                if integral != 0:
                    val = (p / integral) * target_power
                ylabel = "Absorbed Power Density [MW/m³]"
                filename = "Figure_Power_Density_Comparison"

            ax.plot(r, val, label=model_labels[i])

        # After plotting all curves in this subplot
        ymax = 0.0
        for line in ax.lines:
            ymax = max(ymax, np.max(line.get_ydata()))

        ax.set_ylim(0, ymax * 1.05)

        # Formatting subplots
        ax.set_title(f"({chr(97+idx)}) {regime} Plasma", loc='left', fontweight='bold')
        ax.set_xlabel("Radial position r [cm]")
        ax.grid(True, linestyle='--', alpha=0.4)
        
        if idx == 0:
            ax.set_ylabel(ylabel)
        if idx == 1:
            # Place legend only on the middle plot or top to save space
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.16), ncol=3)

    # Save outputs
    plt.savefig(f"{filename}.pdf", bbox_inches='tight')
    plt.savefig(f"{filename}.png", dpi=600, bbox_inches='tight')
    print(f"Saved {filename}")

# --- 4. Generate Figures ---
plot_comparison(metric_type='E')
plot_comparison(metric_type='p')