import matplotlib.pyplot as plt

# Data for the lines
data_const_tau = {
    "$R_{in}$ = 8.0 mm": {
        "x": [300, 600, 1200, 2000, 3000],
        "y": [1.83e18, 4.10e18, 8.36e18, 1.17e19, 2.12e19]
    },
    "$R_{in}$ = 10.0 mm": {
        "x": [600, 1200, 2000],
        "y": [2.42e18, 3.10e18, 6.68e18]
    },
    "$R_{in}$ = 13.5 mm": {
        "x": [300, 600, 1200, 2000],
        "y": [1.26e18, 2.29e18, 1.91e18, 3.34e18]
    },
    "$R_{in}$ = 20.0 mm": {
        "x": [600, 1200, 2000, 3000],
        "y": [2.13e18, 2.32e18, 1.68e18, 2.46e18]
    },
    "$R_{in}$ = 30.0 mm": {
        "x": [300, 600, 1200, 2000, 3000],
        "y": [8.63e17,1.96e18, 1.81e18, 1.43e18, 1.24e18]
    }
}

data_const_flow = {
    "$R_{in}$ = 8.0 mm": {
        "x": [600, 1200, 2000, 3000],
        "y": [3.53e18, 7.65e18, 1.13e19, 1.85e19]
    },
    "$R_{in}$ = 10.0 mm": {
        "x": [600, 1200, 2000],
        "y": [2.62e18, 2.72e18, 6.13e18]
    },
    "$R_{in}$ = 13.5 mm": {
        "x": [300, 600, 1200, 2000],
        "y": [1.26e18, 2.29e18, 1.91e18, 3.34e18]
    },
    "$R_{in}$ = 20.0 mm": {
        "x": [600, 1200, 2000, 3000],
        "y": [2.13e18, 2.37e18, 1.66e18, 2.46e18]
    },
    "$R_{in}$ = 30.0 mm": {
        "x": [600, 1200, 2000, 3000],
        "y": [1.76e18, 2.32e18, 1.68e18, 2.46e18]
    }
}


data = data_const_tau  # data_const_flow or data_const_tau

# --- Paper Ready Settings ---
plt.rcParams.update({
    "font.family": "serif",     # Standard for publications
    "font.size": 16,
    "axes.labelsize": 16,
    "axes.titlesize": 16,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "axes.grid": True,
    "grid.linestyle": "--",
    "grid.alpha": 0.6
})

# Professional markers and colors
markers = ['o', 's', '^', 'd', 'P']  # Circle, square, triangle, diamond, pentagon
colors = ['#1f77b4', '#2ca02c',  '#ff7f0e', '#d62728', "#602d91"]  # Distinct colors for clarity

plt.figure(figsize=(8, 6), dpi=300)

# Plot the data
for i, (label, values) in enumerate(data.items()):
    plt.plot(values["x"], values["y"], 
             marker=markers[i % len(markers)], 
             color=colors[i % len(colors)],
             label=label, 
             linewidth=1.5, 
             markersize=8, 
             markeredgewidth=1.5,
             markerfacecolor='white') # Hollow-style markers for clarity

plt.xlim(250, 3050)
plt.ylim(1e18, 2e19) 
plt.yscale('log')

plt.xlabel('Power Input [W]', fontweight='bold')
plt.ylabel('Peak Electron Density [m⁻³]', fontweight='bold')
plt.legend(frameon=True, loc='best', edgecolor='black')

plt.tight_layout()
plt.savefig('radius_ne_sat_pro.png', bbox_inches='tight')