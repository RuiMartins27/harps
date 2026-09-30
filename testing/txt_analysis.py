import numpy as np
import matplotlib.pyplot as plt


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from collections import defaultdict

electron_charge = 1.602e-19  # C


print_values = False


# --- Plotting Configuration for "Paper Ready" look ---
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 16,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'legend.fontsize': 16,
    'figure.titlesize': 18,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
})

# Load data
filename = "self_consistent/harps/10slm_300mbar_1200W.txt"


data = np.genfromtxt(filename, delimiter=None)
power_string = filename.split('_')[-1].replace('.txt', '').replace('W', '')
P_in = int(power_string)

times = data[:, 0]
radius = data[:, 2]
mu_real = data[:, 3]
elec_dens = data[:, 4]
elec_temp = data[:, 5]
p_abs = data[:, 6]
gas_temp = data[:, 7]
mu_im = data[:, 8]
T_v = data[:, 9]

print("Data loaded from:", filename, "\n")

if(print_values):
    print("Times:", times)
    print("Radius:", radius)
    print("Real part of mobility:", mu_real)
    print("Electron density:", elec_dens)
    print("Electron temperature:", elec_temp)
    print("Absorbed power:", p_abs)
    print("Gas temperature:", gas_temp)
    print("Imaginary part of mobility:", mu_im)
    print("Vibrational temperature:", T_v)

# Max Temperature
print(f"Max(Tg)         = {np.max(gas_temp):.1f} K")

# Radius in which temperature is 50% the Tg in the center
half_Tg = np.max(gas_temp) / 2
radius_half_Tg = radius[np.argmin(np.abs(gas_temp - half_Tg))]
print(f"R_50_Tg         = {radius_half_Tg*100:.3f} cm")

# Max power absorbed distribution
print(f"Max(p_abs/P_in) = {np.max(p_abs)*1e-6/P_in:.3f} cm^-3")


# Integrated Real Conductivity (considering it's cilindrical symmetry and integrating over the radius)
conductivity = elec_dens * electron_charge * mu_real
integrated_conductivity = 2 * np.pi * np.trapezoid(conductivity * radius, radius)
print(f"\\Sigma_R        = {integrated_conductivity:.3e} S m^2")

# Max electron density
print(f"Max(ne)         = {np.max(elec_dens):.3e} m^-3")

# Max electron temperature
print(f"Max(Te)         = {np.max(elec_temp):.1f} K")