import matplotlib.pyplot as plt

# Data for the lines
data_const_tau = {
    "$R_{in}$ = 8.0 mm": {
        "x": [300, 600, 1200, 2000, 3000],
        "y": [5433, 6606, 7592, 7721, 7180]
    },
    "$R_{in}$ = 10.0 mm": {
        "x": [600, 1200, 2000],
        "y": [6314, 6346, 6445]
    },
    "$R_{in}$ = 13.5 mm": {
        "x": [300, 600, 1200, 2000],
        "y": [5378, 6315, 6365, 6199]
    },
    "$R_{in}$ = 20.0 mm": {
        "x": [600, 1200, 2000, 3000],
        "y": [6342, 6455, 6333, 6141]
    },
    "$R_{in}$ = 30.0 mm": {
        "x": [300, 600, 1200, 2000, 3000],
        "y": [4808, 6314, 6386, 6296, 6261]
    }
}

data_const_flow = {
    "$R_{in}$ = 8.0 mm": {
        "x": [600, 1200, 2000, 3000],
        "y": [6426, 8575, 9000, 9000]
    },
    "$R_{in}$ = 10.0 mm": {
        "x": [600, 1200, 2000],
        "y": [6240, 6260, 6405]
    },
    "$R_{in}$ = 13.5 mm": {
        "x": [300, 600, 1200, 2000],
        "y": [5378, 6315, 6365, 6199]
    },
    "$R_{in}$ = 20.0 mm": {
        "x": [600, 1200, 2000, 3000],
        "y": [6355, 6380, 6258, 6108]
    },
    "$R_{in}$ = 30.0 mm": {
        "x": [600, 1200, 2000, 3000],
        "y": [6231, 6213, 6175, 6141]
    }
}

data_other_freq = {
    "$R_{in}$ = 13.5 mm, Static Profile": {
        "x": [300, 600, 1200],
        "y": [5661, 8069, 9000]
    },
    "$R_{in}$ = 13.5 mm, 2.45 GHz, 10 slm": {
        "x": [300, 600, 1200, 2000],
        "y": [5378, 6315, 6365, 6199]
    },
    "$R_{in}$ = 13.5 mm, 915 MHz, 10 slm": {
        "x": [300, 600, 1200, 2000, 3000],
        "y": [5388, 6493, 8051, 8444, 8468]
    },
    "$R_{in}$ = 30 mm, 915 MHz, 10 slm": {
        "x": [600, 1200, 2000, 3000],
        "y": [6344, 6380, 6350, 6277]
    },
    "$R_{in}$ = 30 mm, 915 MHz, 50 slm": {
        "x": [300, 600, 1200, 2000, 3000],
        "y": [4713,6427, 7206, 8126, 8412]
    }
}


data = data_other_freq  # data_const_flow or data_const_tau

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
plt.ylim(4500, 9100) 

plt.xlabel('Power Input [W]', fontweight='bold')
plt.ylabel('Peak Temperature [K]', fontweight='bold')
plt.legend(frameon=True, loc='best', edgecolor='black')

plt.tight_layout()
plt.savefig('radius_tg_sat_pro.png', bbox_inches='tight')