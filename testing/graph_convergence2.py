import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

# 1. Data Setup
N = np.array([64, 128, 256, 512, 768])
P = np.array([196.35, 195.27, 192.66, 191.75, 191.45])
P_ref = 191.21
errors = np.abs((P - P_ref) / P_ref)

# 2. Styling (Using your specified large fonts)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "font.size": 16,             # General base font
    "axes.labelsize": 17,        # X and Y labels
    "axes.titlesize": 17,
    "legend.fontsize": 16.3,       # Legend
    "xtick.labelsize": 15.8,       # Tick numbers
    "ytick.labelsize": 15.8,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 8,
    "xtick.minor.size": 4,
    "ytick.major.size": 8,
    "ytick.minor.size": 4,
    "xtick.top": True,
    "ytick.right": True
})

fig, ax = plt.subplots(figsize=(7, 6))

# 3. Plotting
ax.loglog(N, errors, 'ko-', label='Numerical Result', linewidth=2.5, markersize=8)
ref_n = np.array([128, 768])
ref_err = errors[-1] * (ref_n / N[-1])**(-2) 
ax.loglog(ref_n, ref_err, 'k--', alpha=0.5, label=r'Order 2 Reference ($N^{-2}$)', linewidth=2)

# 1. Use the official Mathtext formatter
# 'minor_thresholds' forced to infinity ensures every tick gets a label
y_formatter = ticker.LogFormatterMathtext(labelOnlyBase=False, minor_thresholds=(np.inf, np.inf))

# 2. Apply to both major and minor
ax.yaxis.set_major_formatter(y_formatter)
ax.yaxis.set_minor_formatter(y_formatter)

# 3. Use your existing locator logic
ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=12))
ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=[4], numticks=12))


# 5. Finishing
ax.set_xlabel('Number of Points ($N$)', labelpad=12)
ax.set_ylabel(r'Relative Error $\frac{|P - P_{ref}|}{P_{ref}}$', labelpad=12)
ax.set_xlim(50, 1000)
ax.grid(True, which="both", ls=":", alpha=0.6)
ax.legend(loc='best', frameon=True, edgecolor='black')

plt.tight_layout()
plt.savefig('absorbed_power_convergence.pdf', bbox_inches='tight')