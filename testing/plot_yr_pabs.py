import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# --- Data ---
y_r = np.array([6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 9.8, 10, 10.1, 10.2, 10.3, 10.5, 10.7, 11, 11.5, 12, 12.5, 13, 13.2, 13.5, 13.8, 14])
p_abs = np.array([57.25, 84.57, 107.05, 134.12, 160.12, 180.62, 203.00, 220.95, 227.06, 229.56, 460.47/2, 460.76/2, 228.57, 224.07, 220.04, 198.73, 154.42, 86.78, 31.02, 8.31, 5.37, 7.76, 14.91, 25.04])

p_abs *= 2  # Scaling as per your requirement
p_in = 770.0

# --- Styling ---
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "legend.fontsize": 14,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "lines.linewidth": 2.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
})

plt.figure(figsize=(9, 4))

# 1. Create a smooth spline for the trend
X_smooth = np.linspace(y_r.min(), y_r.max(), 300)
spl = make_interp_spline(y_r, p_abs, k=3)
plt.plot(X_smooth, spl(X_smooth), color='#1f77b4', linewidth=3, alpha=0.7, zorder=1)

# 2. Plot actual data points (Open markers for academic look)
plt.scatter(y_r, p_abs, color='#1f77b4', edgecolor='black', marker='o', s=20, 
            label='Simulated Data', zorder=3, facecolors='white', linewidths=2)

# 3. Vertical dashed line at Maximum
plt.axhline(y = p_in , color='red', linestyle='--', linewidth=2.5, 
            label=f'Input Power: {p_in:.1f} W', zorder=2)

# 4. Gentle Dashed lines for the peak (Vertical and Horizontal)
# These don't go all the way across the plot, just to the axes for a cleaner look
# Calculate exact peak from data
max_idx = np.argmax(p_abs)
peak_x = y_r[max_idx]
peak_y = p_abs[max_idx]
plt.vlines(x=peak_x, ymin=0, ymax=peak_y, color='gray', linestyle='--', linewidth=1, alpha=0.6, zorder=1)
plt.hlines(y=peak_y, xmin=y_r.min(), xmax=peak_x, color='gray', linestyle='--', linewidth=1, alpha=0.6, zorder=1)

# Annotate the peak value slightly above the point
plt.text(peak_x, peak_y + 20, f'Optimal Position: {peak_y:.1f} W', ha='center', fontsize=14, fontweight='bold')

# --- Labels and Clean Up ---
plt.xlabel(r"Reflector Position $y_r$ [cm]")
plt.ylabel(r"Total Absorbed Power [W]")

plt.grid(True, linestyle=':', alpha=0.6)
plt.legend(
    loc='upper left', 
    bbox_to_anchor=(0.02, 0.9),  # 2% from left, 5% below the top
    frameon=True, 
    edgecolor='black'
)

# Adjust limits for breathing room
plt.xlim(y_r.min() - 0.2, y_r.max() + 0.2)
plt.ylim(0, p_in * 1.06) 

plt.tight_layout()
plt.savefig("absorbed_power_vs_position.pdf", dpi=300)