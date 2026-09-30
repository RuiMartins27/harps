import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# --- 1. Paper-Ready Styling ---
plt.rcParams.update({
    "text.usetex": False, # Set to True if you have LaTeX (e.g., MiKTeX/TeX Live)
    "font.family": "serif",
    "axes.labelsize": 14.2,
    "axes.titlesize": 14.2,
    "xtick.labelsize": 13.2,
    "ytick.labelsize": 13.2,
    "legend.fontsize": 13.5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "lines.linewidth": 2.7,
    "legend.frameon": True,
    "legend.edgecolor": "black",
    "legend.fancybox": False
})

# --- 2. Data Configuration ---
# Update these paths with your actual filenames for each regime
regimes = {
    "Weak": {
        "2D Maxwell": "self_consistent/harps/200mbar_300W.txt",
        "WKB": "self_consistent/WKB/200mbar_300W.txt",
        "Constant E": "self_consistent/E_constant/200mbar_300W.txt"
    },
    "Average": {
        "2D Maxwell": "self_consistent/harps/200mbar_600W.txt",
        "WKB": "self_consistent/WKB/200mbar_600W.txt",
        "Constant E": "self_consistent/E_constant/200mbar_600W.txt"
    },
    "Strong": {
        "2D Maxwell": "self_consistent/harps/200mbar_1200W.txt",
        "WKB": "self_consistent/WKB/200mbar_1200W.txt",
        "Constant E": "self_consistent/E_constant/200mbar_1200W.txt"
    }
}

model_order = ["2D Maxwell", "WKB", "Constant E"]
regime_keys = ["Weak", "Average", "Strong"]
columns = ['time', 'idx', 'radius', 'mu_real', 'n_e', 'T_e', 'p_abs', 'T_g', 'mu_imag', 'T_v']

def get_last_profile(file_path):
    """Reads the file and extracts only the last available timestep profile."""
    try:
        df = pd.read_csv(file_path, sep=r'\s+', names=columns, engine='python')
        # Get the maximum time value to extract the steady-state profile
        last_time = df['time'].max()
        profile = df[df['time'] == last_time].copy()
        # Ensure it is sorted by radius
        profile = profile.sort_values('radius')
        return profile
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return None

# --- 3. Plotting ---
def plot_electron_density():
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.1), sharey=False, constrained_layout=True)
    
    filename = "Electron_Density_Comparison"
    
    for idx, regime_name in enumerate(regime_keys):
        ax = axes[idx]
        regime_files = regimes[regime_name]
        
        for model in model_order:
            file_path = regime_files[model]
            data = get_last_profile(file_path)
            
            if data is not None:
                # Plot n_e vs radius (convert radius to cm if it's in meters)
                # Assuming radius in file is in meters, multiply by 100
                r_cm = data['radius'] * 100 
                ax.plot(r_cm, data['n_e']*1e-18, label=model)

        # Subplot formatting
        ax.set_title(f"({chr(97+idx)}) {regime_name} Plasma", loc='left', fontweight='bold')
        ax.set_xlabel("Radial position $r$ [cm]")
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlim(0, 1.38) # Adjust based on your tube radius
        
        # Only set Y label for the first plot
        if idx == 0:
            ax.set_ylabel(r"Electron Density [$10^{18}$ m$^{-3}$]")
        
        # Place legend under the middle plot
        if idx == 1:
            ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3)

    # Save outputs
    plt.savefig(f"{filename}.pdf", bbox_inches='tight')
    plt.savefig(f"{filename}.png", dpi=300, bbox_inches='tight')
    print(f"Successfully saved {filename}.pdf and .png")

if __name__ == "__main__":
    plot_electron_density()