import numpy as np
import matplotlib.pyplot as plt

# --- Physical Constants ---
e = 1.60217663e-19    # Elementary charge [C]
eps0 = 8.854187e-12   # Vacuum permittivity [F/m]
c = 3e8               # Speed of light [m/s]

# --- Parameters ---
ne = 3e18             # Electron density [m^-3]
eps_r = 1.0           # Background relative permittivity
mu_complex = 0.8 - 2.5j # Complex mobility [m^2/(V·s)]

# Define frequency range
freq_ghz = np.linspace(0.1, 3.0, 500) 
omega = 2 * np.pi * freq_ghz * 1e9

# --- Calculations ---
# 1. Complex conductivity: sigma = e * ne * mu
sigma_complex = e * ne * mu_complex

# 2. Complex relative permittivity: eps_hat = eps_r - i*sigma / (omega * eps0)
epsilon_hat = eps_r - (1j * sigma_complex) / (omega * eps0)

# 3. Propagation constant: k = (omega / c) * sqrt(epsilon_hat)
k_complex = (omega / c) * np.sqrt(epsilon_hat)

# 4. Skin depth is the reciprocal of the imaginary part of k
alpha = np.imag(k_complex)
delta = 1 / np.abs(alpha)

# --- Plotting ---
plt.style.use('seaborn-v0_8-paper')
fig, ax = plt.subplots(figsize=(8, 6), dpi=300) # Slightly larger figure for better scaling

ax.plot(freq_ghz, delta*100, color='#2e7d32', linewidth=3, label=r'Skin Depth $\delta(\omega)$')

# --- Formatting Labels ---
ax.set_xlabel(r'Microwave Frequency [GHz]', fontsize=19, labelpad=10)
ax.set_ylabel(r'Skin Depth $\delta$ [cm]', fontsize=19, labelpad=10)

# --- Increasing Tick Label Size ---
ax.tick_params(axis='both', which='major', labelsize=15)

# Grid and limits
ax.grid(True, linestyle='--', alpha=0.6)
ax.set_xlim(0, 3)
ax.set_ylim(0, 2.5) 

# Parameter Box - Adjusted font size for readability
textstr = (f'$n_e = {ne:.1e}$ m$^{{-3}}$\n'
           f'$\mu = {mu_complex.real} - {abs(mu_complex.imag)}j$ m$^2$/Vs')
props = dict(boxstyle='round', facecolor='wheat', alpha=0.3)
ax.text(0.50, 0.95, textstr, transform=ax.transAxes, fontsize=16,
        verticalalignment='top', bbox=props)

plt.tight_layout()
plt.savefig('skin_depth_plot.png', bbox_inches='tight')

print("Skin depth at 3GHz:", delta[-1]*100, "centimeters")