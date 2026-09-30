import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec

# --- Global Configuration ---
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "stix",
    "font.size": 16,
    "axes.labelsize": 19.5,
    "legend.fontsize": 18,
    "xtick.labelsize": 17,
    "ytick.labelsize": 17,
})

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


def add_overlays(ax, line_y, semicircle_center=(center_x*100, center_y*100), semicircle_radius=1.45):
    """Adds the physical geometry overlays to the plots."""
    # 1. Horizontal line representing the reflector
    ax.axhline(y=line_y+0.1, color='lightgray', linestyle='-', linewidth=5, alpha=0.8)
    
    # 2. Semicircle representing the chamber boundary/plasma region
    # Based on your image, it looks like a half-circle on the left edge
    circ = patches.Arc(semicircle_center, width=semicircle_radius*2, 
                       height=semicircle_radius*2, theta1=270, theta2=90, 
                       color='lightblue', linewidth=3.5, alpha=0.7)
    ax.add_patch(circ)

def plot_field_and_absorption(x_coords, y_coords, field_data, absorption_list, line_positions):
    X, Y = np.meshgrid(x_coords, y_coords, indexing="ij")

    # --- FIGURE 1: Electric Field Magnitude (Single Panel) ---
    fig1, ax1 = plt.subplots(figsize=(4.5, 6))
    mesh1 = ax1.pcolormesh(X, Y, field_data, cmap='plasma', shading='auto')
    add_overlays(ax1, line_y=10.0) # Specific position from your 'optimal' data
    
    cbar1 = fig1.colorbar(mesh1)
    cbar1.set_label(r'$E_{amp}$ [kV/m]', fontsize=18)
    ax1.set_xlabel('$x$ [cm]')
    ax1.set_ylabel('$y$ [cm]')
    fig1.tight_layout()
    fig1.savefig("field_distribution_paper.png", dpi=300)
    fig1.savefig("field_distribution_paper.pdf")

    # --- FIGURE 2: Absorption Comparison (3 Panels, Shared Colorbar) ---
    fig2 = plt.figure(figsize=(12, 6))
    gs = GridSpec(1, 4, width_ratios=[1, 1, 1, 0.1], wspace=0.1)
    
    # Calculate global vmax for consistent comparison
    vmax = max(np.max(d) for d in absorption_list)
    titles = ["$y_r=8.0$ cm", "$y_r=10.2$ cm (Optimal)", "$y_r=12.0$ cm"]
    
    axes = [fig2.add_subplot(gs[0, i]) for i in range(3)]
    
    for i, ax in enumerate(axes):
        mesh = ax.pcolormesh(X, Y, absorption_list[i], cmap='inferno', vmin=0, vmax=vmax, shading='auto')
        add_overlays(ax, line_y=line_positions[i])
        ax.set_title(titles[i])
        ax.set_xlabel('$x$ [cm]')
        if i == 0:
            ax.set_ylabel('$y$ [cm]')
        else:
            ax.set_yticklabels([]) # Hide y-labels for shared look

    # Create the shared colorbar in the 4th slot
    cax = fig2.add_subplot(gs[0, 3])
    cbar2 = fig2.colorbar(mesh, cax=cax)
    cbar2.set_label(r'$P_{abs}$ [W/cm$^3$]', fontsize=18)

    fig2.savefig("absorption_comparison_shared.png", dpi=300, bbox_inches="tight")
    fig2.savefig("absorption_comparison_shared.pdf", bbox_inches="tight")
    plt.show()


data_E_3d, x_raw, y_raw, z_raw = load_3d_data("reflector_position/optimal/E_amplitude.txt")

if data_E_3d is not None:
    # Convert meters to centimeters for the axes
    x_cm = x_raw * 100
    y_cm = y_raw * 100
    
    # Select the first slice (index 0) of the Z-axis to make it 2D
    e_field_2d = data_E_3d[:, :, 0]/1e3

    # 2. Load Absorption Density for all three cases
    def get_abs_slice(path):
        data_3d, _, _, _ = load_3d_data(path)
        return data_3d[:, :, 0] if data_3d is not None else None

    p_abs_pre = get_abs_slice("reflector_position/pre_optimal/absorbed_power_density.txt")
    p_abs_opt = get_abs_slice("reflector_position/optimal/absorbed_power_density.txt")
    p_abs_post = get_abs_slice("reflector_position/post_optimal/absorbed_power_density.txt")

    # 3. Create the list and run the plot
    p_list = [p_abs_pre/1e6, p_abs_opt/1e6, p_abs_post/1e6]
    
    # Ensure all data loaded correctly before plotting
    if all(d is not None for d in p_list):
        # Adjusted positions to match the Y-scale of your data (e.g., ~10 cm)
        pos_list = [8.0, 10.2, 12] 
        
        plot_field_and_absorption(x_cm, y_cm, e_field_2d, p_list, pos_list)
    else:
        print("Error: One or more absorption files could not be loaded.")
else:
    print("Error: Field amplitude file could not be loaded.")
