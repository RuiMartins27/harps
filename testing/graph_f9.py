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
    "legend.fontsize": 11,
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
harps_files = [ 'self_consistent/harps/200mbar_150W.txt', 'self_consistent/harps/200mbar_300W.txt',
                'self_consistent/harps/200mbar_600W.txt', 'self_consistent/harps/200mbar_1200W.txt']
static_files = [ 'self_consistent/static/200mbar_150W.txt', 'self_consistent/static/200mbar_300W.txt',
                 'self_consistent/static/200mbar_600W.txt', 'self_consistent/static/200mbar_1200W.txt']
powers = [150, 300, 600, 1200]

columns = ['time', 'idx', 'radius', 'mu_real', 'n_e', 'T_e', 'p_abs', 'T_g', 'mu_imag', 'T_v']

# --- 3. Plotting Setup (2 Rows, 4 Columns) ---
fig, axs = plt.subplots(2, 4, figsize=(18, 7.5))
axes = axs.flatten()

# Explicitly map panels for organization
ax_pabs  = axes[0]  # Top Row, Col 1
ax_E     = axes[1]  # Top Row, Col 2
ax_mu    = axes[2]  # Top Row, Col 3
ax_sigma = axes[3]  # Top Row, Col 4 (New Conductivity Plot)
ax_ne    = axes[4]  # Bottom Row, Col 1
ax_Tg    = axes[5]  # Bottom Row, Col 2
ax_empty = axes[6]  # Bottom Row, Col 3 (Empty spacer zone)
ax_leg   = axes[7]  # Bottom Row, Col 4 (Reserved for Unified Legend Zone)

ax_empty.axis('off')  # Hide the unused spacer panel

colors = cm.plasma(np.linspace(0.1, 0.8, len(powers)))

# Baseline static plot for Power Density (Using the 300 W reference configuration)
if os.path.exists(static_files[1]):
    df_s_base = pd.read_csv(static_files[1], sep='\s+', names=columns, engine='python')
    ax_pabs.plot(df_s_base['radius']*100, df_s_base['p_abs']*1e-6, color='black', linestyle='--')

# Lambda functions for calculation
calc_E = lambda p_abs, n_e, mu_real: np.sqrt(p_abs / (n_e * electron_charge * mu_real))
calc_sigma = lambda n_e, mu_real: n_e * electron_charge * mu_real

for i, p in enumerate(powers):
    # Process HARPS
    if os.path.exists(harps_files[i]):
        df_h = pd.read_csv(harps_files[i], sep='\s+', names=columns, engine='python')
        
        # Calculate E (V/m) and scale by 1e-2 to convert to V/cm
        E_h = calc_E(df_h['p_abs'], df_h['n_e'], df_h['mu_real']) * 1e-2
        sigma_h = calc_sigma(df_h['n_e'], df_h['mu_real'])
        
        # Plot HARPS Framework
        ax_pabs.plot(df_h['radius']*100, df_h['p_abs']*1e-6, color=colors[i], linestyle='-')
        ax_E.plot(df_h['radius'][:-1]*100, E_h[:-1], color=colors[i], linestyle='-')
        ax_mu.plot(df_h['radius']*100, df_h['mu_real'], color=colors[i], linestyle='-')
        ax_sigma.plot(df_h['radius']*100, sigma_h, color=colors[i], linestyle='-')
        ax_ne.plot(df_h['radius']*100, df_h['n_e'], color=colors[i], linestyle='-')
        ax_Tg.plot(df_h['radius']*100, df_h['T_g'], color=colors[i], linestyle='-')
    
    # Process Static
    if os.path.exists(static_files[i]):
        df_s = pd.read_csv(static_files[i], sep='\s+', names=columns, engine='python')
        
        # Calculate E and Conductivity for Static files
        E_s = calc_E(df_s['p_abs'], df_s['n_e'], df_s['mu_real']) * 1e-2
        sigma_s = calc_sigma(df_s['n_e'], df_s['mu_real'])
        
        # Plot Static Framework
        ax_E.plot(df_s['radius']*100, E_s, color=colors[i], linestyle='--')
        ax_mu.plot(df_s['radius']*100, df_s['mu_real'], color=colors[i], linestyle='--')
        ax_sigma.plot(df_s['radius']*100, sigma_s, color=colors[i], linestyle='--')
        ax_ne.plot(df_s['radius']*100, df_s['n_e'], color=colors[i], linestyle='--')
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

# Panel (d) - Electrical Conductivity
ax_sigma.set_title(r'(d) Conductivity')
ax_sigma.set_ylabel(r'$\sigma$ [S/m]')

# Panel (e) - Electron Density (Log Scale)
ax_ne.set_title(r'(e) Electron Density')
ax_ne.set_ylabel(r'$n_e$ [m$^{-3}$]')
ax_ne.set_yscale('log')

# Panel (f) - Gas Temperature
ax_Tg.set_title(r'(f) Gas Temperature')
ax_Tg.set_ylabel(r'$T_g$ [K]')

for ax in [ax_pabs, ax_E, ax_mu, ax_sigma, ax_ne, ax_Tg]:
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.set_xlabel('Radius $r$ [cm]')

# --- 5. Unified Custom Legend Zone ---
ax_leg.axis('off')  # Clear axis background entirely

power_labels_text = ['150 W', '300 W', '600 W', '1200 W']

# Part A: Absorbed Power Configurations
power_lines = [Line2D([0], [0], color=colors[j], lw=2.5, linestyle='-') for j in range(len(powers))]
power_labels = power_labels_text

leg_power = ax_leg.legend(
    power_lines, power_labels, 
    loc='upper center', 
    bbox_to_anchor=(0.5, 0.98),
    frameon=True, 
    title="Total Power Input",
    title_fontsize=13,
    fontsize=11,
    ncol=2,
    columnspacing=1.2,
    handletextpad=0.6
)
leg_power.get_frame().set_edgecolor('#cccccc')
leg_power.get_frame().set_linewidth(0.8)
ax_leg.add_artist(leg_power)

# Part B: Simulation Framework Styles
framework_lines = [
    Line2D([0], [0], color='gray', lw=2.5, linestyle='-'),
    Line2D([0], [0], color='gray', lw=2.5, linestyle='--'),
    Line2D([0], [0], color='black', lw=2.5, linestyle='--')
]
framework_labels = ['2D EM solver', 'Static', 'Static baseline']

leg_framework = ax_leg.legend(
    framework_lines, framework_labels, 
    loc='lower center', 
    bbox_to_anchor=(0.5, 0.02),
    frameon=True, 
    title="Simulation Framework",
    title_fontsize=13,
    fontsize=11,
    ncol=1,
    labelspacing=0.5,
    handletextpad=0.6
)
leg_framework.get_frame().set_edgecolor('#cccccc')
leg_framework.get_frame().set_linewidth(0.8)

plt.tight_layout()
plt.savefig('Martins2026_f9_grid.pdf')
plt.savefig('Martins2026_f9_grid.png', dpi=300)