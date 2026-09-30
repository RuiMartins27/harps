import numpy as np
import matplotlib.pyplot as plt

# Load data from output file
radius, E_amp, n_e, skin_depth, P_abs = [], [], [], [], []
with open("output.txt") as f:
    for line in f:
        parts = line.replace(",", "").split()
        if len(parts) >= 2:
            r = float(parts[0])
            E = float(parts[1])
            P = float(parts[2])
            radius.append(r)
            E_amp.append(E)
            P_abs.append(P)

radius = np.array(radius)
E_amp = np.array(E_amp)
P_abs = np.array(P_abs)


# Electric field amplitude
plt.figure(figsize=(7,5))
plt.plot(radius*1e3, E_amp, lw=2)
plt.xlabel("y (mm)")
plt.ylabel("|E| (V/m)")
plt.title("Electric Field Amplitude vs Radius")
plt.grid(True)
plt.tight_layout()
plt.savefig("electric_field.png", dpi=300)
plt.close()

# Absorbed power density
plt.figure(figsize=(7,5))
plt.plot(radius*1e3, P_abs, lw=2, color="red")
plt.xlabel("y (mm)")
plt.ylabel("Absorbed Power Density (W/m³)")
plt.title("Power Absorption vs Radius")
plt.grid(True)
plt.tight_layout()
plt.savefig("absorbed_power.png", dpi=300)
plt.close()



print("Saved plots")
