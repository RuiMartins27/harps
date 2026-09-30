import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os
from matplotlib.lines import Line2D

electron_charge = 1.602e-19  # C

# --- 1. Scientific Styling ---
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "lines.linewidth": 2.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "figure.dpi": 300
})

# --- 2. Define Data Structures ---
harps_files = [ 'self_consistent/harps/100mbar_600W.txt', 'self_consistent/harps/200mbar_600W.txt', 
                'self_consistent/harps/300mbar_600W.txt', 'self_consistent/harps/600mbar_600W.txt']
static_files = [ 'self_consistent/static/100mbar_600W.txt', 'self_consistent/static/200mbar_600W.txt', 
                 'self_consistent/static/300mbar_600W.txt', 'self_consistent/static/600mbar_600W.txt']
pressures = [100, 200, 300, 600]

columns = ['time', 'idx', 'radius', 'mu_real', 'n_e', 'T_e', 'p_abs', 'T_g', 'mu_imag', 'T_v']

# --- 3. Plotting Setup (2 Rows, 3 Columns) ---
fig, axs = plt.subplots(2, 3, figsize=(14, 7.5))
axes = axs.flatten()

# Explicitly map panels for organization
ax_pabs = axes[0]  # Top Left
ax_E    = axes[1]  # Top Middle
ax_mu   = axes[2]  # Top Right
ax_ne   = axes[3]  # Bottom Left
ax_Tg   = axes[4]  # Bottom Middle
ax_leg  = axes[5]  # Bottom Right (Reserved for Unified Legend Zone)

# Viridis colormap matching your pressure variants
colors = cm.viridis(np.linspace(0.1, 0.8, len(pressures)))

# Baseline static plot for Power Density (Using the 100 mbar reference configuration)
if os.path.exists(static_files[0]):
    df_s_base = pd.read_csv(static_files[0], sep='\s+', names=columns, engine='python')
    ax_pabs.plot(df_s_base['radius']*100, df_s_base['p_abs']*1e-6, color='black', linestyle='--')

# E-field calculation lambda function
calc_E = lambda p_abs, n_e, mu_real: np.sqrt(p_abs / (n_e * electron_charge * mu_real))

for i, p in enumerate(pressures):
    # Process HARPS Framework
    if os.path.exists(harps_files[i]):
        df_h = pd.read_csv(harps_files[i], sep='\s+', names=columns, engine='python')
        
        # Calculate E (V/m) and scale by 1e-2 to convert to V/cm
        E_h = calc_E(df_h['p_abs'], df_h['n_e'], df_h['mu_real']) * 1e-2
        
        # Top Row Plots
        ax_pabs.plot(df_h['radius']*100, df_h['p_abs']*1e-6, color=colors[i], linestyle='-')
        ax_E.plot(df_h['radius'][:-1]*100, E_h[:-1], color=colors[i], linestyle='-')
        ax_mu.plot(df_h['radius']*100, df_h['mu_real'], color=colors[i], linestyle='-')
        
        # Bottom Row Plots (Applying 1e-18 scaling for Electron Density as requested)
        ax_ne.plot(df_h['radius']*100, df_h['n_e']*1e-18, color=colors[i], linestyle='-')
        ax_Tg.plot(df_h['radius']*100, df_h['T_g'], color=colors[i], linestyle='-')
    
    # Process Static Framework
    if os.path.exists(static_files[i]):
        df_s = pd.read_csv(static_files[i], sep='\s+', names=columns, engine='python')
        
        # Calculate E for Static files
        E_s = calc_E(df_s['p_abs'], df_s['n_e'], df_s['mu_real']) * 1e-2
        
        # Top Row Plots (Static)
        ax_E.plot(df_s['radius']*100, E_s, color=colors[i], linestyle='--')
        ax_mu.plot(df_s['radius']*100, df_s['mu_real'], color=colors[i], linestyle='--')
        
        # Bottom Row Plots (Static)
        ax_ne.plot(df_s['radius']*100, df_s['n_e']*1e-18, color=colors[i], linestyle='--')
        ax_Tg.plot(df_s['radius']*100, df_s['T_g'], color=colors[i], linestyle='--')

# --- 4. Labels and Cosmetics ---

# Panel (a) - Power Density
ax_pabs.set_title(r'(a) Power Density')
ax_pabs.set_ylabel(r'$p_{\mathrm{abs}}$ [W/cm$^3$]')

# Panel (b) - Electric Field
ax_E.set_title(r'(b) Electric Field')
ax_E.set_ylabel(r'$E$ [V/cm]')

# Panel (c) - Real Mobility
ax_mu.set_title(r'(c) Real Mobility')
ax_mu.set_ylabel(r'$\mu_{\mathrm{real}}$ [m$^2$/(V$\cdot$s)]')

# Panel (d) - Electron Density (Linear Scale with scientific multiplier)
ax_ne.set_title(r'(d) Electron Density')
ax_ne.set_ylabel(r'$n_e$ [$10^{18}$ m$^{-3}$]')

# Panel (e) - Gas Temperature
ax_Tg.set_title(r'(e) Gas Temperature')
ax_Tg.set_ylabel(r'$T_g$ [K]')

for ax in [ax_pabs, ax_E, ax_mu, ax_ne, ax_Tg]:
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.set_xlabel('Radius $r$ [cm]')

# --- 5. Unified Custom Legend Zone ---
ax_leg.axis('off')  # Clear axis background entirely

# Configuration lists (Ensure 'colors' vector matches your script's setup)
powers = ['100 mbar', '200 mbar', '300 mbar', '600 mbar']

# Part A: Absorbed Power Configurations (Clean 2-column layout, no baseline)
power_lines = [Line2D([0], [0], color=colors[j], lw=2.5, linestyle='-') for j in range(len(powers))]
power_labels = powers

leg_power = ax_leg.legend(
    power_lines, power_labels, 
    loc='upper center', 
    bbox_to_anchor=(0.5, 0.98),
    frameon=True, 
    title="Pressure Configuration",
    title_fontsize=16,
    fontsize=14,
    ncol=2,
    columnspacing=1.5,
    handletextpad=0.6
)
# Style frame edge for a polished, academic look
leg_power.get_frame().set_edgecolor('#cccccc')
leg_power.get_frame().set_linewidth(0.8)
ax_leg.add_artist(leg_power)  # Lock the power legend down

# Part B: Simulation Framework Styles (Now hosting the Static Baseline)
framework_lines = [
    Line2D([0], [0], color='gray', lw=2.5, linestyle='-'),
    Line2D([0], [0], color='gray', lw=2.5, linestyle='--'),
    Line2D([0], [0], color='black', lw=2.5, linestyle='--')  # Moved Static Baseline here
]
framework_labels = ['2D EM solver', 'Static']

leg_framework = ax_leg.legend(
    framework_lines, framework_labels, 
    loc='lower center', 
    bbox_to_anchor=(0.5, 0.02),
    frameon=True, 
    title="Simulation Framework",
    title_fontsize=16,
    fontsize=14,
    ncol=1,
    labelspacing=0.5,
    handletextpad=0.6
)
# Match the frame styling of the top box
leg_framework.get_frame().set_edgecolor('#cccccc')
leg_framework.get_frame().set_linewidth(0.8)

plt.tight_layout()
plt.savefig('Martins2026_f8_grid.pdf')
plt.savefig('Martins2026_f8_grid.png', dpi=300)