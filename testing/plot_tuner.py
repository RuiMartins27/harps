import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

#Defining center of glass tube
flag_center_coordinates = False
contour = False

center_x = 0
center_y = 0.046
center_z = 0.0

def load_3d_data(filename):
    try:
        with open(filename) as f:
            # Read dimensions from header
            n_x, n_y, n_z, length_x, length_y, length_z = map(float, f.readline().split())        
            n_x, n_y, n_z = int(n_x), int(n_y), int(n_z)

            header = f.readline().strip()

            if header == "NON_UNIFORM_GRID":
                x_coords = np.array(list(map(float, f.readline().split())))
                y_coords = np.array(list(map(float, f.readline().split())))
                z_coords = np.array(list(map(float, f.readline().split())))
                
                data = np.loadtxt(f)
                data = data.reshape((n_x, n_y, n_z))

                # Selecting a line from the data
                #data = data[n_x // 2, :, n_z // 2]; data = data[np.newaxis, :, np.newaxis]

                if(flag_center_coordinates == True):
                    x_coords = x_coords - center_x
                    y_coords = y_coords - center_y
                    z_coords = z_coords - center_z

                return data, x_coords, y_coords, z_coords
            
            else:
                # Uniform grid mode
                x_coords = np.linspace(0, length_x, n_x)
                y_coords = np.linspace(0, length_y, n_y)
                z_coords = np.linspace(0, length_z, n_z)

                data = np.loadtxt(f)
                data = data.reshape((n_x, n_y, n_z))

                if(flag_center_coordinates == True):
                    x_coords = x_coords - center_x
                    y_coords = y_coords - center_y
                    z_coords = z_coords - center_z

                # Remove points where z > 0.5
                z_mask = z_coords <= 0.5
                z_coords = z_coords[z_mask]
                data = data[:, :, z_mask]

                return data, x_coords, y_coords, z_coords
        
    except FileNotFoundError:
        print(f"Warning: File '{filename}' not found. Skipping.")
        return None, None, None, None
    

# --- Global Style Configuration ---
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 16,
    "axes.labelsize": 17,
    "axes.titlesize": 17,
    "legend.fontsize": 16,
    "xtick.labelsize": 14.5,
    "ytick.labelsize": 14.5,
    "xtick.direction": "in",
    "ytick.direction": "in"
})


def plot_yz_field_comparison(y_coords, z_coords, data_no_tuner, data_tuner):
    # Assuming coordinates are in meters, converting to cm for the paper
    Y_mesh, Z_mesh = np.meshgrid(y_coords * 100, z_coords * 100, indexing="ij")

    fig = plt.figure(figsize=(14, 3))
    # width_ratios: two plots + one wider colorbar
    # Using a 0.15 width ratio for the colorbar to make it nice and wide
    gs = GridSpec(1, 3, width_ratios=[1, 1, 0.06], wspace=0.07)
    
    # Identify the X-middle index
    mid_x = data_no_tuner.shape[0] // 2
    
    # Extract the YZ slices [mid_x, :, :]
    # We transpose these (.T) if needed to match the Y-Z orientation in your files
    slice_no_tuner = data_no_tuner[mid_x, :, :]
    slice_tuner = data_tuner[mid_x, :, :]
    
    # Consistent scale for comparison
    vmax = max(np.max(slice_no_tuner), np.max(slice_tuner))
    
    # --- (a) No Tuner ---
    ax0 = fig.add_subplot(gs[0, 0])
    mesh0 = ax0.pcolormesh(Y_mesh, Z_mesh, slice_no_tuner, cmap='plasma', 
                           vmin=0, vmax=vmax, shading='auto')
    ax0.set_title("(a) No Tuner")
    ax0.set_xlabel("$y$ [cm]")
    ax0.set_ylabel("$z$ [cm]")
    ax0.tick_params(top=True, right=True)

    # --- (b) Tuner Optimized ---
    ax1 = fig.add_subplot(gs[0, 1])
    mesh1 = ax1.pcolormesh(Y_mesh, Z_mesh, slice_tuner, cmap='plasma', 
                           vmin=0, vmax=vmax, shading='auto')
    ax1.set_title("(b) Tuner Optimized")
    ax1.set_xlabel("$y$ [cm]")
    ax1.set_yticklabels([]) # Hide redundant Y-labels
    ax1.tick_params(top=True, right=True)

    # --- Shared Wide Colorbar ---
    cax = fig.add_subplot(gs[0, 2])
    cbar = fig.colorbar(mesh1, cax=cax)
    cbar.set_label(r"$E_{\mathrm{amp}}$ [kV/m]", fontsize=17, labelpad=7)
    
    # Instead of tight_layout, we use subplots_adjust to prevent overlap
    fig.subplots_adjust(left=0.08, right=0.92, top=0.90, bottom=0.15)
    
    plt.savefig("YZ_Field_Comparison.png", dpi=300, bbox_inches="tight")
    print("Figure saved as YZ_Field_Comparison.pdf")

# --- Load and Execute ---
# Using your load_3d_data function
data_no, x_c, y_c, z_c = load_3d_data("no_tuner/E_amplitude.txt")
data_tu, _, _, _ = load_3d_data("tuner_optimized/E_amplitude.txt")

if data_no is not None and data_tu is not None:
    plot_yz_field_comparison(y_c, z_c, data_no/1e3, data_tu/1e3)