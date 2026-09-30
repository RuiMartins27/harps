import numpy as np
import matplotlib.pyplot as plt

# --- 1. Paper-Ready Styling ---
plt.rcParams.update({
    "text.usetex": False, # Set to True if you have LaTeX (e.g., MiKTeX/TeX Live)
    "font.family": "serif",
    "axes.labelsize": 14.2,
    "axes.titlesize": 14.2,
    "xtick.labelsize": 13.2,
    "ytick.labelsize": 13.2,
    "legend.fontsize": 14.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 2.6,
    "legend.frameon": True,
    "legend.edgecolor": "black",
    "legend.fancybox": False
})

# --- 2. Data & Settings ---
regimes = {
    "Weak": ["output_harps_weak.txt", "output_WKB_weak.txt", "output_constant_E.txt"],
    "Average": ["output_harps_avg.txt", "output_WKB_avg.txt", "output_constant_E.txt"],
    "Strong": ["output_harps_strong.txt", "output_WKB_strong.txt", "output_constant_E.txt"]
}

model_labels = ["2D Maxwell", "WKB", "E Constant"]
regime_keys = ["Weak", "Average", "Strong"]
L_z = 0.04
target_power = 770 * 1e-6 # MW

# --- 3. Plotting ---
# Create a 2 row x 3 column figure
fig, axes = plt.subplots(2, 3, figsize=(14, 8), sharex=True, constrained_layout=True)

for col, regime in enumerate(regime_keys):
    files = regimes[regime]
    
    # Load reference for E-field normalization
    ref_data = np.loadtxt(files[0], delimiter=",")
    E_ref_max = np.max(ref_data[:, 1])

    for i, f in enumerate(files):
        arr = np.loadtxt(f, delimiter=",")
        r = arr[:, 0] * 100 # cm
        
        # --- Row 0: Electric Field ---
        ax_e = axes[0, col]
        val_e = arr[:, 1] * 1e-3 # kV/m
        val_e *= (E_ref_max * 1e-3) / np.max(arr[:, 1] * 1e-3) # Normalize per your logic
        line = ax_e.plot(r, val_e, label=model_labels[i])
        
        # --- Row 1: Power Density ---
        ax_p = axes[1, col]
        p_raw = arr[:, 2]
        integral = np.trapezoid(L_z * 2 * np.pi * p_raw * (arr[:, 0]), arr[:, 0])
        val_p = (p_raw / integral) * target_power if integral != 0 else p_raw
        ax_p.plot(r, val_p, label=model_labels[i])

    # Subplot Formatting
    for row in range(2):
        ax = axes[row, col]
        ax.grid(True, linestyle='--', alpha=0.4)
        
        # After plotting all curves in this subplot
        ymax = 0.0
        for line in ax.lines:
            ymax = max(ymax, np.max(line.get_ydata()))

        ax.set_ylim(0, ymax * 1.05)

        # Set Titles only on the top row
        if row == 0:
            ax.set_title(f"({chr(97+col)}) {regime} Plasma", loc='left', fontweight='bold')
        
        # Set X-label only on the bottom row
        if row == 1:
            ax.set_xlabel("Radial position r [cm]")

# Axis Labels (Only on the leftmost plots)
axes[0, 0].set_ylabel("Electric Field Amplitude [kV/m]")
axes[1, 0].set_ylabel("Absorbed Power Density [MW/m³]")

# Single Legend for the whole figure
# Positioning it below the subplots
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.06))

# Save output
plt.savefig("Combined_Plasma_Comparison.pdf", bbox_inches='tight')