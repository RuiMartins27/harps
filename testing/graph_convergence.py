import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogLocator, NullFormatter

# 1. Data Points
h_values = [1e-3, 5e-4, 2.5e-4, 1.25e-4, 6.25e-5]
error_values = [3.426e-3, 9.099e-4, 2.95e-4, 1.44e-4, 1.06e-4]
h_vals, err_vals = zip(*sorted(zip(h_values, error_values)))

# 2. Scientific Styling Update (Bigger Fonts)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 16,             # General base font
    "axes.labelsize": 17,        # X and Y labels
    "axes.titlesize": 17,
    "legend.fontsize": 16.3,       # Legend
    "xtick.labelsize": 15.8,       # Tick numbers
    "ytick.labelsize": 15.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 8,
    "xtick.minor.size": 4,
    "ytick.major.size": 8,
    "ytick.minor.size": 4,
    "xtick.top": True,
    "ytick.right": True
})

fig, ax = plt.subplots(figsize=(7, 6))

# 3. Plotting
ax.loglog(h_vals, err_vals, marker='s', linestyle='-', color='black', 
           linewidth=2.5, markersize=8, label='Numerical Error')

h_ref = np.array([min(h_vals), max(h_vals)])
ax.loglog(h_ref, (h_ref**2)*5e3, 'k--', alpha=0.5, label='Order 2 Reference', linewidth=2)

# 4. ENHANCED TICK DENSITY (Log Scale)
# Major ticks every power of 10
ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=12))
ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=12))

# Minor ticks at 2, 3, 4, 5, 6, 7, 8, 9 within each decade
ax.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=12))
ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=12))

# Ensure minor ticks don't have labels (standard for log plots)
ax.xaxis.set_minor_formatter(NullFormatter())
ax.yaxis.set_minor_formatter(NullFormatter())

# 5. Labels and Grid
ax.set_xlabel(r'Step Size $h$ (m)', labelpad=10)
ax.set_ylabel(r'Normalized Error', labelpad=10)
ax.grid(True, which="both", ls=":", linewidth=0.6, alpha=0.6)
ax.legend(loc='best', frameon=True, edgecolor='black')

plt.tight_layout()
plt.savefig('convergence_analytical.pdf', bbox_inches='tight')
plt.show()