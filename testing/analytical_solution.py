import numpy as np
import matplotlib.pyplot as plt


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

                return data, x_coords, y_coords, z_coords
            
            else:
                # Uniform grid mode
                x_coords = np.linspace(0, length_x, n_x)
                y_coords = np.linspace(0, length_y, n_y)
                z_coords = np.linspace(0, length_z, n_z)

                data = np.loadtxt(f)
                data = data.reshape((n_x, n_y, n_z))

                return data, x_coords, y_coords, z_coords
        
    except FileNotFoundError:
        print(f"Warning: File '{filename}' not found. Skipping.")
        return None, None, None, None
    

plotting_enabled = True

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "axes.labelsize": 14,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "text.usetex": False  # Set to True if you have LaTeX installed on your system
})

# Compute in 3D space the analytical solution for microwave inside homogeneous lossy medium inside rectangular waveguide in mode TE10
E_0 = 21000

omega = 2*np.pi*2.45e9      # angular frequency [rad/s]
mu_0 = 1.2566370614e-6           # vacuum permeability [H/m]
epsilon_0 = 8.854187817e-12    # vacuum permittivity [F/m]
electron_charge = 1.602176634e-19   # electron charge [C]

mobility_real = 0.8     # m2/(V⋅s)
mobility_imag = -2.5    # m2/(V⋅s)
electron_density = 6e16

a = 0.08638
b = 0.04318

sigma = electron_density * electron_charge * (mobility_real + 1j*mobility_imag)  # conductivity [S/m]
epsilon_r = 1 - 1j * sigma / (omega * epsilon_0)

print(f"Relative permittivity: {epsilon_r}")


k0 = omega * np.sqrt(mu_0*epsilon_0)
kc = np.pi/a
beta_c = np.sqrt(k0*k0*epsilon_r - kc*kc )

print(f"k0: {k0}")
print(f"kc: {kc}")
print(f"Beta_c: {beta_c}")

nx = ny = nz = 100

x = np.linspace(0, a,   nx)
y = np.linspace(0, b,   ny)
z = np.linspace(0, 0.5, nz)

X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

Ex = np.zeros_like(X)
Ey = (1j*omega*mu_0/kc) * np.sin(kc*X) * np.exp(-1j*beta_c*Z)
Ez = np.zeros_like(X)

Hx = (beta_c / kc) * np.sin(kc*X) * np.exp(-1j * beta_c * Z)
Hy = np.zeros_like(X)
Hz = (1j*omega*epsilon_0*epsilon_r/kc) * np.cos(kc*X) * np.exp(-1j*beta_c*Z)


alpha = E_0/np.max(Ey)

Ex = Ex*alpha; Ey = Ey*alpha; Ez = Ez*alpha
Hx = Hx*alpha; Hy = Hy*alpha; Hz = Hz*alpha

# Plotting the magnitude of the electric field component Ey at y = b/2
if(plotting_enabled):
    plt.figure(figsize=(8, 6))
    plt.contourf(X[:, nx//2, :], Z[:, nx//2, :], np.abs(np.real(Ey[:, nx//2, :])), levels=100, cmap='inferno')
    plt.colorbar(label='|Real(Ey)| (V/m)')
    plt.title('Electric Field Ey at y = b/2')
    plt.xlabel('x (m)')
    plt.ylabel('z (m)')
    plt.savefig("testing/analytical_Ey.png", dpi=300)


# Comparing my solution with analytical solution

data_amp, grid_x, grid_y, grid_z = load_3d_data("Outputs/E_amplitude.txt")
data_real, grid_x, grid_y, grid_z = load_3d_data("Outputs/E_real.txt")
data_imag, grid_x, grid_y, grid_z = load_3d_data("Outputs/E_imag.txt")

n_pml = 800

# plot over axis x = a/2 and y = b/2 both solutions
x_idx = np.size(grid_x)//2
y_idx = np.size(grid_y)//2

E_num_slice = data_real[x_idx, y_idx, :-n_pml]
E_anal_slice = np.abs(np.real(alpha*(1j*omega*mu_0/kc) * np.sin(kc*grid_x[x_idx]) * np.exp(1j*np.pi/1.838-1j*beta_c*grid_z[:-n_pml])))

E_anal_slice = E_anal_slice/np.max(E_anal_slice) #Normalize analytical anad numerical for better comparison
E_num_slice = E_num_slice/np.max(E_num_slice)

plt.figure(figsize=(8, 5.2))

plt.plot(grid_z[:-n_pml], E_anal_slice,  label="Analytical ($TE_{10}$)", color='black',linewidth=2.1, zorder=1)
plt.plot(grid_z[:-n_pml], E_num_slice, label="Numerical Solution", color='#1f77b4', linestyle='None',marker='o', markerfacecolor='none', # Open circle
         markeredgewidth=1.7,markersize=7,markevery=10,zorder=2)

plt.tick_params(direction='in', top=True, right=True, labelsize= 14)
plt.xlabel(r"Position $z$ (m)", fontsize=16)
plt.ylabel(r"Normalized $Real(|E_y|)$", fontsize=16)

plt.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.7)
plt.legend(loc='upper right', frameon=True, fancybox=False, fontsize=15, edgecolor='black')

plt.tight_layout()
plt.savefig("testing/analytical_vs_numerical.png", dpi=300)
plt.savefig("testing/analytical_vs_numerical.pdf")


# plot error between numerical and analytical solution over axis x = a/2 and y = b/2
error = np.abs(E_num_slice - E_anal_slice)
L2_error = np.linalg.norm(error) / np.sqrt(len(error))
print(f"Norma Error: {np.linalg.norm(L2_error)}")

plt.figure(figsize=(10,6))
plt.plot(grid_z[0:(np.size(grid_z)-n_pml)], error[0:(np.size(grid_z)-n_pml)], label="Relative Error", color='red', linewidth=2)
plt.xlabel("z (m)")
plt.ylabel("Error")
plt.title("Error between Numerical and Analytical Solution at x = a/2, y = b/2")
plt.legend()
plt.grid(True)
plt.savefig("testing/error_analytical_numerical.png", dpi=300)


# plot over axis y = b/2 and z = 0.25 both solutions
y_idx = np.size(grid_y)//2
z_idx = np.size(grid_z)//2

E_num_slice = data_real[:, y_idx, z_idx]
E_anal_slice = np.abs(np.real(alpha*(1j*omega*mu_0/kc) * np.sin(kc*grid_x) * np.exp(1j*np.pi/1.8215-1j*beta_c*grid_z[z_idx])))

E_anal_slice = E_anal_slice/np.max(E_anal_slice) #Normalize analytical anad numerical for better comparison
E_num_slice = E_num_slice/np.max(E_num_slice)

plt.figure(figsize=(10,6))
plt.plot(grid_x, E_num_slice, label="Numerical HARPS", linewidth=2)
plt.plot(grid_x, E_anal_slice, "--", label="Analytical TE10", linewidth=2)

plt.xlabel("x (m)")
plt.ylabel("|E_y| (V/m)")
plt.title("Electric Field at y = b/2 z = 0.25 m")
plt.legend()
plt.grid(True)
plt.show()
