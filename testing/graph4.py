import numpy as np
import matplotlib.pyplot as plt


N = np.array([32, 48, 64, 96])  # must be a NumPy array
n_pts = N*N*N
time = np.array([13.2, 45.7, 151.2, 1080.0])  # in seconds

error_x = [0.1621, 0.0870, 0.1003, 0.0446] 
error_y = [0.0614, 0.0565, 0.0236, 0.0105]  
error_amp = [0.01475, 0.01015, 0.00810, 0.00360] 

# fit to quadratic scaling and plot
log_n = np.log(n_pts)
log_time = np.log(time)
p, logA = np.polyfit(log_n, log_time, 1)
A = np.exp(logA)

print(f"Fitted power law: time = {A:.3e} * n_pts^{p:.2f}")
time_fit = A * n_pts**p

# Plot
plt.figure(figsize=(8,5))
plt.loglog(n_pts, time, 'o', label='Data')
plt.loglog(n_pts, time_fit, '-', label=f'Fit: time ~ n_pts^{p:.2f}')
plt.xlabel('Number of points (n^3)')
plt.ylabel('Computation time [s]')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.savefig("testing/t_n.png", dpi=300)

# Fit errors to error ~ (1/N)^p using log-log
invN = 1.0 / N

# Fit errors to error ~ (1/N)^p
def fit_error(error):
    log_error = np.log(error)
    log_invN = np.log(invN)
    p, logC = np.polyfit(log_invN, log_error, 1)  # linear fit in log-log
    C = np.exp(logC)
    return C, p

C_x, p_x = fit_error(error_x)
C_y, p_y = fit_error(error_y)
C_amp, p_amp = fit_error(error_amp)

print(f"Error_x fit: C={C_x:.4f}, p={p_x:.2f}")
print(f"Error_y fit: C={C_y:.4f}, p={p_y:.2f}")
print(f"Error_amp fit: C={C_amp:.4f}, p={p_amp:.2f}")

# Plot
plt.figure(figsize=(8,5))
plt.loglog(invN, error_x, 'o', label='Error X')
plt.loglog(invN, C_x*invN**p_x, '-', label=f'Fit X: p={p_x:.2f}')
plt.loglog(invN, error_y, 's', label='Error Y')
plt.loglog(invN, C_y*invN**p_y, '-', label=f'Fit Y: p={p_y:.2f}')
plt.loglog(invN, error_amp, '^', label='Error Amp')
plt.loglog(invN, C_amp*invN**p_amp, '-', label=f'Fit Amp: p={p_amp:.2f}')

plt.xlabel('1/N')
plt.ylabel('Error')
plt.legend()
plt.grid(True, which='both')
plt.savefig("testing/errors_h.png", dpi=300)

