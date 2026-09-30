import pandas as pd
import matplotlib.pyplot as plt
import os

# --- 1. Scientific Styling Configuration ---
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "lines.linewidth": 3,
    "figure.dpi": 300,
    "text.usetex": False  # Set to True if you have LaTeX installed
})

# --- 2. Load Data ---
file_path = 'self_consistent/harps/200mbar_600W.txt'

if not os.path.exists(file_path):
    print(f"File {file_path} not found. Ensure the directory structure is correct.")
else:
    columns = ['time', 'idx', 'radius', 'mu_real', 'n_e', 'T_e', 'p_abs', 'T_g', 'mu_imag', 'T_v']
    df = pd.read_csv(file_path, sep='\s+', names=columns, engine='python')

    df['radius'] *= 100
    df['p_abs']  *= 1e-6

    # --- 3. Plotting ---
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12, 5))

    # --- Plot (a): Electron Density and Temperature vs Radius ---
    ax_a_twin = ax_a.twinx()
    
    # Normalized profiles to compare shapes (n_e ~10^18 vs T_e ~10^4)
    ne_norm = df['n_e'] / df['n_e'].max()
    te_norm = df['T_e'] / df['T_e'].max()

    l1, = ax_a.plot(df['radius'], ne_norm, color='#d62728', label='$n_e / n_{e,\mathrm{max}}$')
    l2, = ax_a.plot(df['radius'], te_norm, color='#1f77b4', label='$T_e / T_{e,\mathrm{max}}$')
    l3, = ax_a_twin.plot(df['radius'], df['p_abs'], color='black', linestyle='--', label='$p_{\mathrm{abs}}$')

    ax_a.set_xlabel('Radius $r$ [cm]')
    ax_a.set_ylabel('Normalized $n_e, T_e$')
    ax_a_twin.set_ylabel('$p_{\mathrm{abs}}$ [W/cm$^3$]')
    ax_a.set_title('(a)')
    ax_a.grid(True, linestyle=':', alpha=0.6)
    ax_a.legend([l1, l2, l3], [l.get_label() for l in [l1, l2, l3]], loc='lower left', frameon=True)

    # --- Plot (b): Heavy Species Temperatures vs Radius ---
    ax_b_twin = ax_b.twinx()

    l4, = ax_b.plot(df['radius'], df['T_g'], color='#1f77b4', label='$T_g$')
    l5, = ax_b.plot(df['radius'], df['T_v'], color='#2ca02c', label='$T_v$')
    l6, = ax_b_twin.plot(df['radius'], df['p_abs'], color='black', linestyle='--', label='$p_{\mathrm{abs}}$')

    ax_b.set_xlabel('Radius $r$ [cm]')
    ax_b.set_ylabel('Temperature [K]')
    ax_b_twin.set_ylabel('$p_{\mathrm{abs}}$ [W/cm$^3$]')
    ax_b.set_title('(b)')
    ax_b.grid(True, linestyle=':', alpha=0.6)
    ax_b.legend([l4, l5, l6], [l.get_label() for l in [l4, l5, l6]], loc='lower left', frameon=True)

    plt.tight_layout()
    plt.savefig('Martins2026_f7.pdf', bbox_inches='tight')