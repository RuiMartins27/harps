import numpy as np
import matplotlib.pyplot as plt

# --- Global Configuration for Presentations ---
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 17,
    "axes.labelsize": 19,
    "legend.fontsize": 17,
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": '--'
})

data = {
    "Static Profile": {
        "x": [300, 600, 1200],
        "y": [5661, 8069, 9000],
        "color": "#2c3e50", # Dark Navy
        "marker": "o"
    },
    "Coupled at 2.45 GHz": {
        "x": [300, 600, 1200, 2000],
        "y": [5378, 6315, 6365, 6199],
        "color": "#e74c3c", # Soft Red
        "marker": "s"
    },
}

fig, ax = plt.subplots(figsize=(10, 6.5))

for label, values in data.items():
    ax.plot(values["x"], values["y"], 
            label=label, 
            color=values["color"], 
            marker=values["marker"], 
            markersize=10, 
            linewidth=2.5, 
            markerfacecolor='white', # Hollow effect
            markeredgewidth=2)

# Customizing Axes
ax.set_xlabel("Power [W]") # Replace with your actual unit
ax.set_ylabel("Peak Temperature [K]")
ax.legend(frameon=True, facecolor='white', framealpha=1, loc='lower right')

ax.set_ylim(4900, 9200)
ax.set_xlim(200, 3050)

plt.tight_layout()
plt.savefig("Tg_Saturation_Comparison.png", dpi=300, bbox_inches="tight")