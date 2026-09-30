import numpy as np

# ── 1.  Load data produced by extract_colormap_data.py ───────────────────────
Tg     = np.load("Tg.npy")       # shape (n_r, n_z)
r_axis = np.load("r_axis.npy")   # shape (n_r,)   [mm]

# ── 2.  Average and std over the axial (z) dimension  ────────────────────────
T_mean = Tg.mean(axis=1)    # shape (n_r,)
T_std  = Tg.std(axis=1)     # shape (n_r,)

# ── 3.  Pack into a single (n_r, 3) array ────────────────────────────────────
#   profile[i] = [r_pos (mm), T_g (K), std_Tg (K)]
profile = np.column_stack([r_axis, T_mean, T_std])

# ── 4.  Save ──────────────────────────────────────────────────────────────────
np.save("radial_profile.npy", profile)
np.savetxt(
    "radial_profile.csv",
    profile,
    delimiter=",",
    header="r_pos_mm,T_g_K,std_Tg_K",
    comments="",
)

# ── 5.  Summary ───────────────────────────────────────────────────────────────
print(f"profile shape : {profile.shape}  →  (n_r, 3)")
print(f"Columns       : r_pos [mm] | T_g [K] | std_Tg [K]")
print(f"\nSample (every 50th row):")
print(f"{'r (mm)':>10}  {'T_g (K)':>10}  {'std (K)':>10}")
for row in profile[::50]:
    print(f"{row[0]:>10.3f}  {row[1]:>10.1f}  {row[2]:>10.1f}")
print("\nSaved  radial_profile.npy  and  radial_profile.csv")

# ── 6.  Optional quick plot ───────────────────────────────────────────────────
try:
    import matplotlib.pyplot as plt

    r   = profile[:, 0]
    Tm  = profile[:, 1]
    Ts  = profile[:, 2]

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.plot(Tm, r, "k-", linewidth=1.5, label="mean $T_g$")
    ax.fill_betweenx(r, Tm - Ts, Tm + Ts, alpha=0.25, label="±1 std")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("Transverse position (mm)")
    ax.set_title("Axially averaged temperature profile")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig("radial_profile.png", dpi=150)
    print("Verification plot saved to radial_profile.png")
except ImportError:
    print("matplotlib not available – skipping plot")