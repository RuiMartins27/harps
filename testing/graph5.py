import numpy as np
import matplotlib.pyplot as plt


N = np.array([32 ,64, 128, 256, 512, 1024])  # must be a NumPy array
n_pts = N*N

time = np.array([0.3,0.26, 0.59, 1.52, 4.75, 19.3])  # in seconds

error_amp = [0.011323, 0.004669, 0.003078,  0.002718, 0.002626, 0.002605]  # from previous runs

# ---------------------------------------------------------
# Fit errors to power law: error = A + C * (1/N)^p  (log-log fit)
# ---------------------------------------------------------
invN = 1.0 / N
print(invN)



# ---------------------------------------------------------
# Fit time to power-law scaling: time = C * N^p
# ---------------------------------------------------------
logN = np.log(n_pts)
logT = np.log(time)

p_time, logC_time = np.polyfit(logN, logT, 1)
C_time = np.exp(logC_time)

print(f"Time fit: time ≈ {C_time:.4e} * N^{p_time:.3f}")

# Plot time scaling (log-log)
plt.figure(figsize=(6,5))
plt.loglog(n_pts, time, 'o-', label='Data')
plt.loglog(n_pts, C_time * n_pts**p_time, '--', label=f'Fit: p={p_time:.2f}')
plt.xlabel("n_pts")
plt.ylabel("Time (s)")
plt.title("Runtime scaling (power law)")
plt.grid(True, which='both', ls='--')
plt.legend()
plt.savefig("testing/time_scaling.png", dpi=300)