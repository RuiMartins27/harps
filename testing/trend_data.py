import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker

data = [
    # Inflow = 5
    {"Inflow [slm]": 5, "p [mbar]": 200, "P_in [W]": 600,  "max(Tg) [K]": 6366, "R_50,Tg [cm]": 1.205, "max(Pabs)/P_in": 0.105, "Sigma_R [S/m^2]": 3.208e-5, "max_ne [x10^18 m^-3]": 2.30, "max_Te [K]": 10017},
    {"Inflow [slm]": 5, "p [mbar]": 200, "P_in [W]": 1200, "max(Tg) [K]": 6395, "R_50,Tg [cm]": 1.272, "max(Pabs)/P_in": 0.126, "Sigma_R [S/m^2]": 6.636e-5, "max_ne [x10^18 m^-3]": 2.47, "max_Te [K]": 12298},
    {"Inflow [slm]": 5, "p [mbar]": 200, "P_in [W]": 2000, "max(Tg) [K]": 6224, "R_50,Tg [cm]": 1.306, "max(Pabs)/P_in": 0.164, "Sigma_R [S/m^2]": 12.76e-5,  "max_ne [x10^18 m^-3]": 3.72, "max_Te [K]": 58030},
    {"Inflow [slm]": 5, "p [mbar]": 300, "P_in [W]": 200,  "max(Tg) [K]": 4531, "R_50,Tg [cm]": 0.969, "max(Pabs)/P_in": 0.454, "Sigma_R [S/m^2]": 0.234e-5,  "max_ne [x10^18 m^-3]": 0.67, "max_Te [K]": 9984},
    
    # Inflow = 20
    {"Inflow [slm]": 20, "p [mbar]": 200, "P_in [W]": 600,  "max(Tg) [K]": 6364, "R_50,Tg [cm]": 1.070, "max(Pabs)/P_in": 0.183, "Sigma_R [S/m^2]": 2.171e-5, "max_ne [x10^18 m^-3]": 2.53, "max_Te [K]": 11302},
    {"Inflow [slm]": 20, "p [mbar]": 200, "P_in [W]": 1200, "max(Tg) [K]": 6372, "R_50,Tg [cm]": 1.272, "max(Pabs)/P_in": 0.123, "Sigma_R [S/m^2]": 4.567e-5, "max_ne [x10^18 m^-3]": 1.99, "max_Te [K]": 10933},
    {"Inflow [slm]": 20, "p [mbar]": 200, "P_in [W]": 2000, "max(Tg) [K]": 5098, "R_50,Tg [cm]": 1.340, "max(Pabs)/P_in": 0.188, "Sigma_R [S/m^2]": 9.473e-5, "max_ne [x10^18 m^-3]": 3.23, "max_Te [K]": 37831},
    {"Inflow [slm]": 20, "p [mbar]": 300, "P_in [W]": 800,  "max(Tg) [K]": 6768, "R_50,Tg [cm]": 1.036, "max(Pabs)/P_in": 0.245, "Sigma_R [S/m^2]": 1.896e-5, "max_ne [x10^18 m^-3]": 4.51, "max_Te [K]": 10761},

    # Inflow = 10
    {"Inflow [slm]": 10, "p [mbar]": 50,  "P_in [W]": 600,  "max(Tg) [K]": 3196, "R_50,Tg [cm]": 1.306, "max(Pabs)/P_in": 0.150, "Sigma_R [S/m^2]": 17.70e-5, "max_ne [x10^18 m^-3]": 1.23, "max_Te [K]": 13373},
    {"Inflow [slm]": 10, "p [mbar]": 100, "P_in [W]": 600,  "max(Tg) [K]": 4163, "R_50,Tg [cm]": 1.306, "max(Pabs)/P_in": 0.144, "Sigma_R [S/m^2]": 5.774e-5, "max_ne [x10^18 m^-3]": 0.775, "max_Te [K]": 11796}, 
    {"Inflow [slm]": 10, "p [mbar]": 200, "P_in [W]": 150,  "max(Tg) [K]": 3728, "R_50,Tg [cm]": 1.070, "max(Pabs)/P_in": 0.398, "Sigma_R [S/m^2]": 0.303e-5, "max_ne [x10^18 m^-3]": 0.436, "max_Te [K]": 9913},
    {"Inflow [slm]": 10, "p [mbar]": 200, "P_in [W]": 300,  "max(Tg) [K]": 5378, "R_50,Tg [cm]": 1.104, "max(Pabs)/P_in": 0.265, "Sigma_R [S/m^2]": 0.956e-5, "max_ne [x10^18 m^-3]": 1.26, "max_Te [K]": 10442},
    {"Inflow [slm]": 10, "p [mbar]": 200, "P_in [W]": 600,  "max(Tg) [K]": 6315, "R_50,Tg [cm]": 1.171, "max(Pabs)/P_in": 0.142, "Sigma_R [S/m^2]": 2.553e-5, "max_ne [x10^18 m^-3]": 2.29, "max_Te [K]": 10715},
    {"Inflow [slm]": 10, "p [mbar]": 200, "P_in [W]": 1200, "max(Tg) [K]": 6364, "R_50,Tg [cm]": 1.272, "max(Pabs)/P_in": 0.130, "Sigma_R [S/m^2]": 5.610e-5, "max_ne [x10^18 m^-3]": 1.906, "max_Te [K]": 10819},
    {"Inflow [slm]": 10, "p [mbar]": 200, "P_in [W]": 2000, "max(Tg) [K]": 6199, "R_50,Tg [cm]": 1.306, "max(Pabs)/P_in": 0.168, "Sigma_R [S/m^2]": 9.230e-5, "max_ne [x10^18 m^-3]": 3.34, "max_Te [K]": 10283},
    {"Inflow [slm]": 10, "p [mbar]": 300, "P_in [W]": 300,  "max(Tg) [K]": 5444, "R_50,Tg [cm]": 1.036, "max(Pabs)/P_in": 0.349, "Sigma_R [S/m^2]": 0.489e-5, "max_ne [x10^18 m^-3]": 1.27, "max_Te [K]": 10093},
    {"Inflow [slm]": 10, "p [mbar]": 300, "P_in [W]": 600,  "max(Tg) [K]": 6432, "R_50,Tg [cm]": 1.070, "max(Pabs)/P_in": 0.239, "Sigma_R [S/m^2]": 1.511e-5, "max_ne [x10^18 m^-3]": 3.35, "max_Te [K]": 9975},
    {"Inflow [slm]": 10, "p [mbar]": 600, "P_in [W]": 300,  "max(Tg) [K]": 4869, "R_50,Tg [cm]": 0.935, "max(Pabs)/P_in": 0.562, "Sigma_R [S/m^2]": 0.102e-5, "max_ne [x10^18 m^-3]": 0.766, "max_Te [K]": 9834.5},
    {"Inflow [slm]": 10, "p [mbar]": 600, "P_in [W]": 600,  "max(Tg) [K]": 6466, "R_50,Tg [cm]": 0.867, "max(Pabs)/P_in": 0.665, "Sigma_R [S/m^2]": 0.483e-5, "max_ne [x10^18 m^-3]": 5.78, "max_Te [K]": 9547.5},
    {"Inflow [slm]": 10, "p [mbar]": 300, "P_in [W]": 400,  "max(Tg) [K]": 5644, "R_50,Tg [cm]": 1.070, "max(Pabs)/P_in": 0.297, "Sigma_R [S/m^2]": 0.734e-5, "max_ne [x10^18 m^-3]": 1.72, "max_Te [K]": 9796},
    {"Inflow [slm]": 10, "p [mbar]": 50, "P_in [W]": 130,  "max(Tg) [K]": 2863, "R_50,Tg [cm]": 1.070, "max(Pabs)/P_in": 0.38, "Sigma_R [S/m^2]": 3.009e-5, "max_ne [x10^18 m^-3]": 0.963, "max_Te [K]": 12055},
    {"Inflow [slm]": 10, "p [mbar]": 300, "P_in [W]": 800,  "max(Tg) [K]": 6507, "R_50,Tg [cm]": 1.137, "max(Pabs)/P_in": 0.185, "Sigma_R [S/m^2]": 2.264e-5, "max_ne [x10^18 m^-3]": 4.204, "max_Te [K]": 9677},
    {"Inflow [slm]": 10, "p [mbar]": 300, "P_in [W]": 1200,  "max(Tg) [K]": 6879, "R_50,Tg [cm]": 1.171, "max(Pabs)/P_in": 0.159, "Sigma_R [S/m^2]": 3.559e-5, "max_ne [x10^18 m^-3]": 6.561, "max_Te [K]": 10043},

    # Inflow = 30
    {"Inflow [slm]": 30, "p [mbar]": 300,  "P_in [W]": 1200,  "max(Tg) [K]": 8016, "R_50,Tg [cm]": 0.969, "max(Pabs)/P_in": 0.199, "Sigma_R [S/m^2]": 3.255e-5, "max_ne [x10^18 m^-3]": 8.23, "max_Te [K]": 11864}
]

# Create DataFrame
df = pd.DataFrame(data)

sns.set_theme(style="whitegrid")
# High contrast color palette
hc_palette = sns.color_palette("bright")

# ADDED max_Te to this list (now exactly 6 variables)
vars_to_plot_lines = ["max(Tg) [K]", "R_50,Tg [cm]", "max(Pabs)/P_in", "Sigma_R [S/m^2]", "max_ne [x10^18 m^-3]", "max_Te [K]"]

def create_6_panel_plot(data, x_col, hue_col, title, filename):
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle(title, fontsize=18, fontweight='bold')
    axes = axes.flatten()
    
    for i, var in enumerate(vars_to_plot_lines):
        sns.lineplot(
            data=data, x=x_col, y=var, hue=hue_col, 
            marker='o', markersize=8, linewidth=2.5, ax=axes[i], palette=hc_palette[:data[hue_col].nunique()]
        )
        axes[i].set_title(f'{var} vs {x_col}', fontsize=12)
        axes[i].set_xlabel(x_col, fontsize=11)
        axes[i].set_ylabel(var, fontsize=11)
        axes[i].set_ylim(bottom=0)
        if var == "max_Te [K]":
            axes[i].set_ylim(0, 15000)
        
    # We no longer need to remove an empty subplot here since all 6 are used
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

# 1. Plot at 200 mbar
df_200mbar = df[df["p [mbar]"] == 200].sort_values(by="P_in [W]")
create_6_panel_plot(
    data=df_200mbar, x_col="P_in [W]", hue_col="Inflow [slm]", 
    title="Effect of Input Power at Constant Pressure (200 mbar)", 
    filename="plot_200mbar.png"
)

# 1.2 Plot at 300 mbar
df_300mbar = df[df["p [mbar]"] == 300].sort_values(by="P_in [W]")
create_6_panel_plot(
    data=df_300mbar, x_col="P_in [W]", hue_col="Inflow [slm]", 
    title="Effect of Input Power at Constant Pressure (300 mbar)", 
    filename="plot_300mbar.png"
)

# 2. Plot at 10 slm
df_10slm = df[df["Inflow [slm]"] == 10].sort_values(by="p [mbar]")
create_6_panel_plot(
    data=df_10slm, x_col="p [mbar]", hue_col="P_in [W]", 
    title="Effect of Pressure at Constant Inflow (10 slm)", 
    filename="plot_10slm.png"
)

# 3. Plot at 600 W
df_600W = df[df["P_in [W]"] == 600].sort_values(by="p [mbar]")
create_6_panel_plot(
    data=df_600W, x_col="p [mbar]", hue_col="Inflow [slm]", 
    title="Effect of Pressure and Inflow at Constant Input Power (600 W)", 
    filename="plot_600W.png"
)

# 4. Plot with Sigma as X axis
df["Sigma_R [10 ^4 S/m^2]"] = df["Sigma_R [S/m^2]"] * 1e4

sns.set_theme(style="whitegrid")

# ADDED max_Te to this list (now 5 variables)
vars_to_plot_scatter = ["max(Tg) [K]", "R_50,Tg [cm]", "max(Pabs)/P_in", "max_ne [x10^18 m^-3]", "max_Te [K]"]

# Changed grid from 2x2 to 2x3 to accommodate the 5th plot
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Relationship of Plasma Parameters vs. Electrical Conductivity ($\Sigma_R$)', fontsize=16, fontweight='bold')

axes = axes.flatten()

for i, var in enumerate(vars_to_plot_scatter):
    sns.scatterplot(
        data=df, 
        x="Sigma_R [10 ^4 S/m^2]", 
        y=var, 
        hue="p [mbar]", 
        size="P_in [W]", 
        sizes=(50, 250), 
        palette="viridis", 
        alpha=0.8,
        ax=axes[i]
    )
    axes[i].set_title(f'{var} vs $\Sigma_R$', fontsize=13)
    axes[i].set_xlabel('$\Sigma_R$ $[10^4 S/m^2]$', fontsize=11)
    axes[i].set_ylabel(var, fontsize=11)

    # 1. Make x-axis logarithmic
    axes[i].set_xscale('log')
    
    # 2. THE FIX: Use LogLocator with numticks instead of locator_params(nbins)
    axes[i].xaxis.set_major_locator(ticker.LogLocator(numticks=5))
    
    # Use scientific notation formatting nicely for the x-axis
    axes[i].xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    axes[i].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    # Rotate the ticks slightly for better readability
    for label in axes[i].get_xticklabels():
        label.set_rotation(30)
        label.set_horizontalalignment('right')
    
    # Anchor the legend to the top-right plot (index 2 in a 3-column grid)
    if i == 2:
        axes[i].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    else:
        axes[i].get_legend().remove()

# Remove the empty 6th subplot
fig.delaxes(axes[5])

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.savefig('plot_sigma_x_cm2.png', dpi=300)
plt.close()

print("All 4 plots generated successfully.")