import matplotlib.pyplot as plt
import numpy as np

# Load data
data = np.loadtxt("Outputs/results_yr_pabs.txt")

# Split into yr and Pabs
yr = data[:, 0]
pabs = data[:, 1]

# Sort by yr for a smooth curve
sort_idx = np.argsort(yr)
yr_sorted = yr[sort_idx]
pabs_sorted = pabs[sort_idx]

# Create the plot (wider, less tall)
plt.figure(figsize=(8, 3))
plt.scatter(100*yr, 2*pabs, c="blue", s=30, label="Samples")
plt.plot(100*yr_sorted, 2*pabs_sorted, alpha=0.7, label="Trend")

plt.xlabel("Reflector Position [cm]", fontsize=16)
plt.ylabel("Absorbed Power [W]", fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)

plt.tight_layout()
plt.savefig("analysis/analysis_output/reflector.png", dpi=300)
