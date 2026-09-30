import numpy as np
from PIL import Image
from scipy.spatial import cKDTree

# ── 1.  Load image ────────────────────────────────────────────────────────────
IMG_PATH = "exp_data/Peter_N2_800W.png"

img = Image.open(IMG_PATH).convert("RGB")
arr = np.array(img)          # shape (H, W, 3), dtype uint8

# ── 2.  Hard-coded pixel coordinates (determined from the 895×563 px image) ──
#
#   Plot frame  (inside the tick marks, excluding axis lines themselves)
PLOT_TOP   = 13    # first data row
PLOT_BOT   = 464   # last  data row
PLOT_LEFT  = 101   # first data col
PLOT_RIGHT = 704   # last  data col

#   Colour-bar strip  (single-pixel-wide column sampled through its centre)
CB_COL = 733       # x-coordinate of colour-bar centre column
CB_TOP = 12        # top    pixel of colour bar  → T_MAX
CB_BOT = 465       # bottom pixel of colour bar  → T_MIN

# ── 3.  Physical axis ranges (read from axis labels in the image) ─────────────
T_MAX = 7000.0   # K  (top    of colour bar)
T_MIN = 2000.0   # K  (bottom of colour bar)

Z_MIN =  0.0     # mm, axial position (left  edge of plot)
Z_MAX = 17.5     # mm, axial position (right edge of plot)
R_MAX =  6.0     # mm, transverse position (top    of plot, positive)
R_MIN = -6.0     # mm, transverse position (bottom of plot, negative)

# ── 4.  Build colour-bar look-up table ────────────────────────────────────────
cb_rgb = arr[CB_TOP : CB_BOT + 1, CB_COL, :3].astype(np.float32)   # (N_cb, 3)
n_cb   = len(cb_rgb)

# Temperature assigned to each colour-bar row (linear mapping top→T_MAX)
T_lut  = np.linspace(T_MAX, T_MIN, n_cb)                            # (N_cb,)

# KD-tree on RGB space for fast nearest-neighbour colour matching
tree = cKDTree(cb_rgb)

# ── 5.  Map every plot pixel to a temperature ─────────────────────────────────
plot_rgb = arr[PLOT_TOP : PLOT_BOT + 1,
               PLOT_LEFT : PLOT_RIGHT + 1, :3].astype(np.float32)

n_r, n_z = plot_rgb.shape[:2]

_, nn_idx = tree.query(plot_rgb.reshape(-1, 3), k=1, workers=-1)
Tg = T_lut[nn_idx].reshape(n_r, n_z)    # shape (n_r, n_z), float32

# ── 6.  Physical coordinate arrays ───────────────────────────────────────────
r_axis = np.linspace(R_MAX, R_MIN, n_r)   # mm, top → bottom
z_axis = np.linspace(Z_MIN, Z_MAX, n_z)   # mm, left → right

# ── 7.  Quick sanity check ───────────────────────────────────────────────────
print(f"Tg shape : {Tg.shape}  →  (n_r={n_r}, n_z={n_z})")
print(f"r_axis   : {r_axis[0]:.2f} … {r_axis[-1]:.2f} mm")
print(f"z_axis   : {z_axis[0]:.2f} … {z_axis[-1]:.2f} mm")
print(f"T range  : {Tg.min():.0f} … {Tg.max():.0f} K  (mean {Tg.mean():.0f} K)")

# ── 8.  Optional: save to disk ───────────────────────────────────────────────
np.save("Tg.npy",     Tg)
np.save("r_axis.npy", r_axis)
np.save("z_axis.npy", z_axis)
print("Saved Tg.npy, r_axis.npy, z_axis.npy")

# ── 9.  Optional: quick re-plot to verify the reconstruction ─────────────────
try:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4))
    pcm = ax.pcolormesh(z_axis, r_axis, Tg, cmap="jet",
                        vmin=T_MIN, vmax=T_MAX, shading="auto")
    fig.colorbar(pcm, ax=ax, label="Temperature (K)")
    ax.set_xlabel("Axial position (mm)")
    ax.set_ylabel("Transverse position (mm)")
    ax.set_title("Reconstructed Tg from Peter_N2_800W.png")
    plt.tight_layout()
    plt.savefig("Tg_reconstructed.png", dpi=150)
    print("Verification plot saved to Tg_reconstructed.png")
except ImportError:
    print("matplotlib not available – skipping verification plot")