import os
import numpy as np
import matplotlib.pyplot as plt

# --- Global Journal Style Settings ---
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
    'figure.titlesize': 14
})

# Define file directory and labels
base_dir = "self_consistent/harps/v_profile"
files = {
    "v_ratio_1.txt": "V-Ratio = 1",
    "v_ratio_4.txt": "V-Ratio = 4",
    "v_ratio_10.txt": "V-Ratio = 10",
    "v_ratio_40.txt": "V-Ratio = 40"
}

# Set up a single unified figure
plt.figure(figsize=(8, 6))

for filename, label in files.items():
    file_path = os.path.join(base_dir, filename)
    
    # Check if file exists before processing
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found. Skipping.")
        continue
        
    # Load raw text data
    data = np.loadtxt(file_path)
    
    # Read all variables explicitly as requested
    times     = data[:, 0]
    # data[:, 1] is the index/counter column
    radius    = data[:, 2]
    mu_real   = data[:, 3]
    elec_dens = data[:, 4]
    elec_temp = data[:, 5]
    p_abs     = data[:, 6]
    gas_temp  = data[:, 7]  # This is the target 'tg' variable
    mu_im     = data[:, 8]
    T_v       = data[:, 9]
    
    # Plot tg vs Radius for this profile on the same set of axes
    plt.plot(radius, elec_dens, label=label, linewidth=2, marker='o', markersize=3)

# Label axes and title
plt.xlabel('Radius (m)')
plt.ylabel('Electron Density (m⁻³)')

# Add legend with a clear box border matching typical journal presentation
plt.legend(frameon=True, edgecolor='gainsboro', facecolor='white')

# Auto-adjust padding
plt.tight_layout()

# Save complete high-res visualization
output_plot = os.path.join(base_dir, "elec_dens_combined_comparison.png")
plt.savefig(output_plot, dpi=300)
print(f"Combined plot successfully saved to: {output_plot}")