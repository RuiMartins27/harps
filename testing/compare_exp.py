import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# --- Global Journal Style Settings ---
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linestyle': '--',
    'figure.titlesize': 14
})

# 1. Define file paths
sim_files = [
    "self_consistent/harps/300mbar_300W.txt",
    "self_consistent/harps/10slm_50mbar_130W.txt",
    "self_consistent/harps/5slm_300mbar_200W.txt",
    "self_consistent/harps/20slm_300mbar_800W.txt"
]

exp_files = [
    "exp_data/radial_profile_sasha_300W_200mbar.csv",
    "exp_data/radial_profile_gatti_120W.csv",
    "exp_data/radial_profile_peter_200W.csv",
    "exp_data/radial_profile_peter_800W.csv"
]

# Cleaned up titles (Labels a, b, c, d will be prepended dynamically)
titles = [
    "Sasha: 300W / 10 slm / 300 mbar",
    "Gatti: 120W / 10 slm / 50 mbar",
    "Peter: 200W / 5 slm / 300 mbar",
    "Peter: 800W / 20 slm / 300 mbar"
]

# 2. Initialize a 2x2 subplot figure (optimized aspect ratio for papers)
fig, axs = plt.subplots(2, 2, figsize=(12, 9.5))
axs = axs.ravel()  # Flatten for simple looping

# Defining a cohesive, high-contrast academic color palette
color_sim = '#004488'  # Elegant Deep Blue
color_exp = '#D55E00'  # Clear Crimson/Orange

for i in range(4):
    # --- Load and plot Simulation Data ---
    if os.path.exists(sim_files[i]):
        try:
            sim_data = pd.read_csv(sim_files[i], sep=r'\s+', header=None)
            sim_r = sim_data.iloc[:, 2] * 1000.0  # Convert m to mm
            sim_T = sim_data.iloc[:, 7]           # Gas Temperature (K)
            
            axs[i].plot(sim_r, sim_T, label='Simulation', color=color_sim, linewidth=2, zorder=2)
        except Exception as e:
            print(f"Error reading simulation file {sim_files[i]}: {e}")
    else:
        print(f"Simulation file not found: {sim_files[i]}")

    # --- Load and plot Experimental Data ---
    if os.path.exists(exp_files[i]):
        try:
            exp_data = pd.read_csv(exp_files[i])
            exp_r = exp_data['r_pos_mm']
            exp_T = exp_data['T_g_K']
            
            # Dynamically adjust marker size based on data density to prevent big blobs
            m_size = 14 if len(exp_r) > 50 else 28
            
            if 'std_Tg_K' in exp_data.columns:
                exp_std = exp_data['std_Tg_K']
                # Plot error bars lighter/thinner so they don't choke out dense plots
                axs[i].errorbar(exp_r, exp_T, yerr=exp_std, fmt='none', ecolor=color_exp,
                                elinewidth=0.8, capsize=1.5, alpha=0.4, zorder=3)
                
                # Overlay markers cleanly on top
                axs[i].scatter(exp_r, exp_T, label='Experiment', color=color_exp, 
                            s=m_size, edgecolors='k', linewidths=0.6, zorder=4)
            else:
                axs[i].scatter(exp_r, exp_T, label='Experiment', color=color_exp, 
                            s=m_size, edgecolors='k', linewidths=0.6, zorder=4)
        except Exception as e:
            print(f"Error reading experimental file {exp_files[i]}: {e}")
    else:
        print(f"Experimental file not found: {exp_files[i]}")

    # --- Formatting each Subplot ---
    # Add standard journal panel identifiers (a), (b), (c), (d)
    panel_label = f"({chr(97 + i)})"
    axs[i].set_title(f"{panel_label} {titles[i]}", loc='left', pad=10, fontweight='semibold')
    
    axs[i].set_xlabel(r'$r$ (mm)')
    axs[i].set_ylabel(r'$T_g$ (K)')
    
    # Borderless, clean legend layout
    axs[i].legend(loc='upper right', frameon=False, fontsize=10)
    
    # Smart tick allocation
    axs[i].xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
    axs[i].yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))

# Adjust spacing cleanly to ensure no overlapping text elements
plt.tight_layout()

# Save options
output_filename = 'sim_vs_exp_comparison.png'
plt.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Plot successfully saved as '{output_filename}'")

plt.savefig('sim_vs_exp_comparison.pdf', dpi=300, bbox_inches='tight')
print(f"Plot successfully saved as 'sim_vs_exp_comparison.pdf'")

plt.close()