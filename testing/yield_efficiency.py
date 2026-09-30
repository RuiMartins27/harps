import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from matplotlib.ticker import ScalarFormatter, NullFormatter, LogLocator

# --- 1. Data Setup (Kept from your snippet) ---
pressures = [50, 100, 200, 300, 600]
powers = [150, 300, 600, 1200, 2000]

yield_p_harps = [0.072, 0.032, 0.095, 0.113, 0.118]
eff_p_harps = [0.840, 0.366, 1.100, 1.309, 1.366]
yield_p_static = [1.039, 0.943, 0.396, 0.166, 0.064]
eff_p_static = [12.079, 10.958, 4.606, 1.929, 0.742]

yield_w_harps = [0.003, 0.029, 0.095, 0.149, 0.200]
eff_w_harps = [0.139, 0.674, 1.100, 0.866, 0.697]
yield_w_static = [0.005, 0.047, 0.396, 1.004, 1.406]
eff_w_static = [0.228, 1.084, 4.606, 5.838, 4.904]


# --- 2. Styling (Updated with larger fonts and pro-ticks) ---
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 16,
    "axes.labelsize": 17,
    "axes.titlesize": 16,
    "legend.fontsize": 16, # Adjusted for a balanced "large" look
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 7,
    "xtick.minor.size": 4,
    "ytick.major.size": 7,
    "ytick.minor.size": 4,
    "lines.linewidth": 3,
    "figure.dpi": 300
})

fig, (ax_p, ax_w) = plt.subplots(1, 2, figsize=(14, 6))

def plot_performance(ax, x_data, y_yield_h, y_eff_h, y_yield_s, y_eff_s, x_label, title):
    ax_eff = ax.twinx()
    
    # Plotting lines
    l1, = ax.plot(x_data, y_yield_h, 'o-', color='#1f77b4', markersize=9)
    l2, = ax_eff.plot(x_data, y_eff_h, 's--', color='#1f77b4', markersize=9, markerfacecolor='none')
    l3, = ax.plot(x_data, y_yield_s, 'o-', color='#d62728', markersize=9)
    l4, = ax_eff.plot(x_data, y_eff_s, 's--', color='#d62728', markersize=9, markerfacecolor='none')
    
    ax.set_xlabel(x_label)
    ax.set_xscale('log')
    ax.set_ylabel('Yield [%]')
    ax_eff.set_ylabel('Efficiency [%]') # or Efficiency [%] per your preference
    ax.set_title(title)
    
    # Professional Ticking for Log Axis
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xticks(x_data) # Ensure your specific points have major ticks
    
    ax.grid(True, which='major', linestyle=':', alpha=0.6)
    ax.grid(True, which='minor', linestyle=':', alpha=0.3)
    
    return [l1, l2, l3, l4]

# --- 3. Execute Plotting ---
plot_performance(ax_p, pressures, yield_p_harps, eff_p_harps, yield_p_static, eff_p_static, 
                 'Pressure [mbar]', '(a)                 ')

plot_performance(ax_w, powers, yield_w_harps, eff_w_harps, yield_w_static, eff_w_static, 
                 'Power [W]', '              (b)')

# --- 2. THE UNIFIED LEGEND (Bottom Centered) ---
legend_elements = [
    Line2D([0], [0], color='#1f77b4', lw=3, label='Self-Consistent'),
    Line2D([0], [0], color='#d62728', lw=3, label='Static'),
    Line2D([0], [0], marker='o', color='gray', label='Yield', markerfacecolor='gray', markersize=10, linestyle='none'),
    Line2D([0], [0], marker='s', color='gray', label='Efficiency', markerfacecolor='none', markersize=10, linestyle='none')
]

# ncol=4 puts them all in one horizontal line
fig.legend(handles=legend_elements, loc='lower center', ncol=4, 
           frameon=True, edgecolor='black', bbox_to_anchor=(0.5, 0.92), fontsize=14.5)

# Add space at the bottom for the legend
plt.subplots_adjust(bottom=0.22, wspace=0.35)

plt.tight_layout()
plt.savefig('Martins2026_f10.pdf', bbox_inches='tight')