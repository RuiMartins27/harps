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

# 1. Define explicit single file paths and title
sim_file = "self_consistent/harps/10slm_300mbar_400W.txt"
exp_file = "exp_data/radial_profile_omid.csv"
title = "Omid: 400W / 10 slm / 300 mbar"

# 2. Initialize a single figure (optimized for single-column journal layout)
fig, ax = plt.subplots(figsize=(6.5, 4.8))

# Cohesive, high-contrast academic color palette
color_sim = '#004488'  # Elegant Deep Blue
color_exp = '#D55E00'  # Clear Crimson/Orange

# --- Load and plot Simulation Data ---
if os.path.exists(sim_file):
    try:
        sim_data = pd.read_csv(sim_file, sep=r'\s+', header=None)
        sim_r = sim_data.iloc[:, 2] * 1000.0  # Convert m to mm
        sim_T = sim_data.iloc[:, 7]           # Gas Temperature (K)
        
        ax.plot(sim_r, sim_T, label='Simulation', color=color_sim, linewidth=2, zorder=2)
    except Exception as e:
        print(f"Error reading simulation file {sim_file}: {e}")
else:
    print(f"Simulation file not found: {sim_file}")

# --- Load and plot Experimental Data ---
if os.path.exists(exp_file):
    try:
        exp_data = pd.read_csv(exp_file)
        exp_r = exp_data['r_pos_mm']
        exp_T = exp_data['T_g_K']
        
        # Dynamically adjust marker size based on data density
        m_size = 14 if len(exp_r) > 50 else 28
        
        if 'std_Tg_K' in exp_data.columns:
            exp_std = exp_data['std_Tg_K']
            ax.errorbar(exp_r, exp_T, yerr=exp_std, fmt='none', ecolor=color_exp,
                        elinewidth=0.8, capsize=1.5, alpha=0.4, zorder=3)
            
            ax.scatter(exp_r, exp_T, label='Experiment', color=color_exp, 
                       s=m_size, edgecolors='k', linewidths=0.6, zorder=4)
        else:
            ax.scatter(exp_r, exp_T, label='Experiment', color=color_exp, 
                       s=m_size, edgecolors='k', linewidths=0.6, zorder=4)
    except Exception as e:
        print(f"Error reading experimental file {exp_file}: {e}")
else:
    print(f"Experimental file not found: {exp_file}")

# --- Formatting ---
ax.set_title(title, loc='left', pad=10, fontweight='semibold')
ax.set_xlabel(r'$r$ (mm)')
ax.set_ylabel(r'$T_g$ (K)')

# Borderless, clean legend layout
ax.legend(loc='upper right', frameon=False, fontsize=10)

# Smart tick allocation
ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))

# Adjust layout cleanly
plt.tight_layout()

# Save options
output_filename = 'sim_vs_exp_comparison'
plt.savefig(f'{output_filename}.png', dpi=300, bbox_inches='tight')
plt.savefig(f'{output_filename}.pdf', dpi=300, bbox_inches='tight')

print(f"Plots successfully saved as {output_filename}.png and .pdf")
plt.close()