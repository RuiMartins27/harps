import numpy as np
import subprocess
import os
import argparse
import matplotlib.pyplot as plt

harps_dir = os.path.dirname(os.path.abspath(__file__)) + "/../"
script_dir = os.path.dirname(os.path.abspath(__file__))

import numpy as np

# Constants
I = 1j
e_charge = 1.602e-19  # C
eps0 = 8.854e-12      # F/m
c0 = 3.0e8            # m/s
free = -0.00123456


def E_field_gauss_trig(x, y, z, alpha):
    Ex = np.sin(alpha*x) * np.sin(alpha*y) * np.sin(alpha*z)
    Ey = np.sin(alpha*x) * np.sin(alpha*y) * np.sin(alpha*z)
    Ez = np.sin(alpha*x) * np.sin(alpha*y) * np.sin(alpha*z)

    return Ex, Ey, Ez

def E_field_gauss_trig_2D(x, y, z, alpha):
    Ex = 0
    Ey = 0
    Ez = np.sin(alpha*x)*np.sin(alpha*y)

    return Ex, Ey, Ez

def ne_profile_gauss(x, y, z, ne0, x0, y0, z0, L):
    rx = x - x0
    ry = y - y0
    rz = z - z0
    r2 = rx*rx + ry*ry + rz*rz

    ne = ne0*np.exp(-r2/(L*L))
    return ne


def mms_source_all_equalE(x, y, z, sigma0, x0, y0, z0, L, omega, eps0, k0, alpha):
    """
    Manufactured-source S(x,y,z) for:
       ∇²E + ∇[(∇σ·E)/(i ω ε0 ε̂_r)] + k0² ε̂_r E = S.

    Uses E(x,y,z) = sin(alpha x) sin(alpha y) sin(alpha z) on all components and σ = σ0 exp(-|r|² / L²).
    """

    # --- E0 and partials ---
    E0 = np.sin(alpha*x) * np.sin(alpha*y) * np.sin(alpha*z)

    dE0_dx =  alpha * np.cos(alpha*x) * np.sin(alpha*y) * np.sin(alpha*z)
    dE0_dy =  alpha * np.sin(alpha*x) * np.cos(alpha*y) * np.sin(alpha*z)
    dE0_dz =  alpha * np.sin(alpha*x) * np.sin(alpha*y) * np.cos(alpha*z)

    # Laplacian term factor (componentwise)
    lap_factor = -3.0 * alpha**2

    # --- Gaussian sigma and its derivatives ---
    rx = x - x0
    ry = y - y0
    rz = z - z0
    r2 = rx*rx + ry*ry + rz*rz

    sigma = sigma0 * np.exp(- r2 / (L*L))

    # first derivatives
    sig_x = -2.0 * rx / (L*L) * sigma
    sig_y = -2.0 * ry / (L*L) * sigma
    sig_z = -2.0 * rz / (L*L) * sigma

    # second derivatives needed for grad S1
    # sigma_xx = sigma*(4*rx^2/L^4 - 2/L^2), etc.
    sig_xx = sigma * (4.0*rx*rx / (L**4) - 2.0 / (L**2))
    sig_yy = sigma * (4.0*ry*ry / (L**4) - 2.0 / (L**2))
    sig_zz = sigma * (4.0*rz*rz / (L**4) - 2.0 / (L**2))

    # mixed second derivatives (sigma_xy, sigma_xz, sigma_yz)
    sig_xy = sigma * (4.0*rx*ry / (L**4))
    sig_xz = sigma * (4.0*rx*rz / (L**4))
    sig_yz = sigma * (4.0*ry*rz / (L**4))

    # --- S1 and grad S1 ---
    S1 = sig_x + sig_y + sig_z

    # grad S1 = (partial_x S1, partial_y S1, partial_z S1)
    # partial_x S1 = sigma_xx + sigma_xy + sigma_xz, etc.
    dS1_dx = sig_xx + sig_xy + sig_xz
    dS1_dy = sig_xy + sig_yy + sig_yz
    dS1_dz = sig_xz + sig_yz + sig_zz

    # --- complex permittivity and A ---
    eps_hat = 1.0 - I * sigma / (omega * eps0)
    A = I * omega * eps0 * eps_hat   # A = i ω ε0 ε̂_r

    # --- assemble S components ---
    # common "mass + laplacian" term (complex)
    common = lap_factor + (k0*k0) * eps_hat   # scalar complex

    Sx = common * E0
    Sy = common * E0
    Sz = common * E0

    # gradient of phi:
    # grad phi = ( grad E0 * S1 + E0 * grad S1 ) / A  -  (E0*S1 * grad sigma) / A^2

    # first part: (grad E0 * S1) / A
    part1_x = dE0_dx * S1 / A
    part1_y = dE0_dy * S1 / A
    part1_z = dE0_dz * S1 / A

    # second part: (E0 * grad S1) / A
    part2_x = E0 * dS1_dx / A
    part2_y = E0 * dS1_dy / A
    part2_z = E0 * dS1_dz / A

    # third part: - (E0*S1 * grad sigma) / A^2
    A_sq = A * A
    part3_x = - E0 * S1 * sig_x / A_sq
    part3_y = - E0 * S1 * sig_y / A_sq
    part3_z = - E0 * S1 * sig_z / A_sq

    # add grad phi to the S components
    Sx = Sx + part1_x + part2_x + part3_x
    Sy = Sy + part1_y + part2_y + part3_y
    Sz = Sz + part1_z + part2_z + part3_z

    return Sx, Sy, Sz

def mms_source_all_equalE_numerical(x_grid, y_grid, z_grid, sigma0, x0, y0, z0, L, omega, eps0, k0, alpha):
    """
    Numerical Manufactured-source S(x,y,z) using non uniform grid for:
       ∇²E + ∇[(∇σ·E)/(i ω ε0 ε̂_r)] + k0² ε̂_r E = S.

    Uses E(x,y,z) = sin(alpha x) sin(alpha y) sin(alpha z) on all components and σ = σ0 exp(-|r|² / L²).
    """
    Nx = len(x_grid)
    Ny = len(y_grid)
    Nz = len(z_grid)

    dx = np.diff(x_grid)
    dy = np.diff(y_grid)
    dz = np.diff(z_grid)

    Ex = np.zeros((Nx, Ny, Nz), dtype=complex)
    Ey = np.zeros((Nx, Ny, Nz), dtype=complex)
    Ez = np.zeros((Nx, Ny, Nz), dtype=complex)

    sig = np.zeros((Nx, Ny, Nz), dtype=complex)

    phi = np.zeros((Nx, Ny, Nz), dtype=complex)
    
    Sx = np.zeros((Nx, Ny, Nz), dtype=complex)
    Sy = np.zeros((Nx, Ny, Nz), dtype=complex)
    Sz = np.zeros((Nx, Ny, Nz), dtype=complex)

    for i in range(1,Nx):
        for j in range(1,Ny):
            for k in range(1,Nz):
                Ex[i,j,k] = np.sin(alpha*x_grid[i]) * np.sin(alpha*y_grid[j]) * np.sin(alpha*z_grid[k])
                Ey[i,j,k] = Ex[i,j,k]
                Ez[i,j,k] = Ex[i,j,k]

                sig[i,j,k] = sigma0 * np.exp(- ((x_grid[i]-x0)**2 + (y_grid[j]-y0)**2 + (z_grid[k]-z0)**2) / (L*L))

    # Compute grad sigma
    for i in range(1,Nx):
        for j in range(1,Ny):
            for k in range(1,Nz):
                if i == 0 or i == Nx-1 or j == 0 or j == Ny-1 or k == 0 or k == Nz-1: 
                    continue

                x = x_grid[i]
                y = y_grid[j]
                z = z_grid[k]

                rx = x - x0
                ry = y - y0
                rz = z - z0
                r2 = rx*rx + ry*ry + rz*rz

                sigma = sigma0 * np.exp(- r2 / (L*L))

                eps_hat = 1.0 - I * sigma / (omega * eps0)
                A = I * omega * eps0 * eps_hat   # A = i ω ε0 ε̂_r

                grad_sig_x = -sig[i-1,j,k]*dx[i]/((dx[i-1])*(dx[i]+dx[i-1])) + sig[i+1,j,k]*dx[i-1]/((dx[i])*(dx[i]+dx[i-1])) + \
                          sig[i,j,k]*(dx[i]*dx[i]-dx[i-1]*dx[i-1])/((dx[i]*dx[i-1])*(dx[i]+dx[i-1]))
                grad_sig_y = -sig[i,j-1,k]*dy[j]/((dy[j-1])*(dy[j]+dy[j-1])) + sig[i,j+1,k]*dy[j-1]/((dy[j])*(dy[j]+dy[j-1])) + \
                          sig[i,j,k]*(dy[j]*dy[j]-dy[j-1]*dy[j-1])/((dy[j]*dy[j-1])*(dy[j]+dy[j-1]))
                grad_sig_z = -sig[i,j,k-1]*dz[k]/((dz[k-1])*(dz[k]+dz[k-1])) + sig[i,j,k+1]*dz[k-1]/((dz[k])*(dz[k]+dz[k-1])) + \
                          sig[i,j,k]*(dz[k]*dz[k]-dz[k-1]*dz[k-1])/((dz[k]*dz[k-1])*(dz[k]+dz[k-1]))

                phi[i,j,k] = (grad_sig_x*Ex[i,j,k] + grad_sig_y*Ey[i,j,k] + grad_sig_z*Ez[i,j,k])/A

    # Compute S field
    for i in range(1,Nx):
        for j in range(1,Ny):
            for k in range(1,Nz):
                if i == 0 or i == Nx-1 or j == 0 or j == Ny-1 or k == 0 or k == Nz-1: 
                    continue

                x = x_grid[i]
                y = y_grid[j]
                z = z_grid[k]

                rx = x - x0
                ry = y - y0
                rz = z - z0
                r2 = rx*rx + ry*ry + rz*rz

                sigma = sigma0 * np.exp(- r2 / (L*L))

                eps_hat = 1.0 - I * sigma / (omega * eps0)
                
                # Laplacian term factor (componentwise) + k0² ε̂_r E term
                Sx[i,j,k] = -2*Ex[i,j,k]*(1/(dx[i]*dx[i-1]) + 1/(dy[j]*dy[j-1]) + 1/(dz[k]*dz[k-1])) + \
                            2*Ex[i+1,j,k]/(dx[i]* (dx[i]+dx[i-1])) + 2*Ex[i-1,j,k]/(dx[i-1]* (dx[i]+dx[i-1])) + \
                            2*Ex[i,j+1,k]/(dy[j]* (dy[j]+dy[j-1])) + 2*Ex[i,j-1,k]/(dy[j-1]* (dy[j]+dy[j-1])) + \
                            2*Ex[i,j,k+1]/(dz[k]* (dz[k]+dz[k-1])) + 2*Ex[i,j,k-1]/(dz[k-1]* (dz[k]+dz[k-1])) + \
                            (k0*k0) * eps_hat * Ex[i,j,k]
                Sy[i,j,k] = -2*Ey[i,j,k]*(1/(dx[i]*dx[i-1]) + 1/(dy[j]*dy[j-1]) + 1/(dz[k]*dz[k-1])) + \
                            2*Ey[i+1,j,k]/(dx[i]* (dx[i]+dx[i-1])) + 2*Ey[i-1,j,k]/(dx[i-1]* (dx[i]+dx[i-1])) + \
                            2*Ey[i,j+1,k]/(dy[j]* (dy[j]+dy[j-1])) + 2*Ey[i,j-1,k]/(dy[j-1]* (dy[j]+dy[j-1])) + \
                            2*Ey[i,j,k+1]/(dz[k]* (dz[k]+dz[k-1])) + 2*Ey[i,j,k-1]/(dz[k-1]* (dz[k]+dz[k-1])) + \
                            (k0*k0) * eps_hat * Ey[i,j,k]
                Sz[i,j,k] = -2*Ez[i,j,k]*(1/(dx[i]*dx[i-1]) + 1/(dy[j]*dy[j-1]) + 1/(dz[k]*dz[k-1])) + \
                            2*Ez[i+1,j,k]/(dx[i]* (dx[i]+dx[i-1])) + 2*Ez[i-1,j,k]/(dx[i-1]* (dx[i]+dx[i-1])) + \
                            2*Ez[i,j+1,k]/(dy[j]* (dy[j]+dy[j-1])) + 2*Ez[i,j-1,k]/(dy[j-1]* (dy[j]+dy[j-1])) + \
                            2*Ez[i,j,k+1]/(dz[k]* (dz[k]+dz[k-1])) + 2*Ez[i,j,k-1]/(dz[k-1]* (dz[k]+dz[k-1])) + \
                            (k0*k0) * eps_hat * Ez[i,j,k]

                # Gradient Term (numerical)
                Sx[i,j,k] = Sx[i,j,k] - phi[i-1,j,k]*dx[i]/((dx[i-1])*(dx[i]+dx[i-1])) + \
                                        phi[i+1,j,k]*dx[i-1]/((dx[i])*(dx[i]+dx[i-1])) + \
                                        phi[i,j,k]*(dx[i]*dx[i]-dx[i-1]*dx[i-1])/((dx[i]*dx[i-1])*(dx[i]+dx[i-1]))
            
                Sy[i,j,k] = Sy[i,j,k] - phi[i,j-1,k]*dy[j]/((dy[j-1])*(dy[j]+dy[j-1])) + \
                                        phi[i,j+1,k]*dy[j-1]/((dy[j])*(dy[j]+dy[j-1])) + \
                                        phi[i,j,k]*(dy[j]*dy[j]-dy[j-1]*dy[j-1])/((dy[j]*dy[j-1])*(dy[j]+dy[j-1]))

                Sz[i,j,k] = Sz[i,j,k] - phi[i,j,k-1]*dz[k]/((dz[k-1])*(dz[k]+dz[k-1])) + \
                                        phi[i,j,k+1]*dz[k-1]/((dz[k])*(dz[k]+dz[k-1])) + \
                                        phi[i,j,k]*(dz[k]*dz[k]-dz[k-1]*dz[k-1])/((dz[k]*dz[k-1])*(dz[k]+dz[k-1]))


    return Sx, Sy, Sz


def mms_source_2D(x, y, z, sigma0, x0, y0, z0, L, omega, eps0, k0, alpha):
    """
    Manufactured-source S(x,y,z) for:
       ∇²E + k0² ε̂_r E = S.

    Uses E(x,y,z) = sin(alpha x) sin(alpha y) on z direction only  and σ = σ0 exp(-|r|² / L²).
    """
    Ex = 0
    Ey = 0
    Ez = np.sin(alpha*x) * np.sin(alpha*y)

    lap_factor = -2*alpha**2

    sigma = ne_profile_gauss(x, y, z, sigma0, x0, y0, z0, L)

    eps_hat = 1.0 - I * sigma / (omega*eps0)
    common = lap_factor + (k0*k0)*eps_hat   # scalar complex

    Sx = 0
    Sy = 0
    Sz = common * Ez

    return Sx, Sy, Sz

def generate_grid(n_x, n_y, n_z, length_x, length_y, length_z, refinementFactor_x=-1, refinementFactor_y=-1, refinementFactor_z=-1,
    x_BL=0, x_RL=0, x_RR=0, x_BR=0, y_BL=0, y_RL=0, y_RR=0, y_BR=0, z_BL=0, z_RL=0, z_RR=0, z_BR=0, x_grid=None, y_grid=None, z_grid=None, x_file=None, y_file=None, z_file=None):

    def process_direction(n, length, refinementFactor, BL, BR, RL, RR, grid, file):
        if grid is not None and len(grid) > 0:
            arr = np.array(grid)
            if not np.all(np.diff(arr) >= 0):
                raise ValueError("Input grid is not monotonically increasing.")
            d = np.diff(arr)
            return arr, d, True
        
        elif file is not None and os.path.exists(file):
            arr = np.loadtxt(file, ndmin=1)
            if not np.all(np.diff(arr) >= 0):
                raise ValueError("Input file grid is not monotonically increasing.")
            d = np.diff(arr)
            return arr, d, True
        
        elif refinementFactor > 1:
            if n < 2:
                raise ValueError("n must be at least 2 for refinement.")
            if (abs(BL) < 1e-6) and (abs(RL) < 1e-6) and (BR > RR):
                symmetric = True
            elif (BL > RL or RL > RR or RR > BR) and not(abs(BL) < 1e-6) and (abs(RL) < 1e-6) and (BR > RR):
                raise ValueError("Key points must be monotonically increasing: BL <= RL <= RR <= BR. It can also be only one sided (0, 0, RR, BR)")
            
            arr = np.zeros(n)
            d = np.zeros(n-1)

            normal = (BL+length - BR+(RR-RL)*refinementFactor+refinementFactor/(refinementFactor-1)*(RL-BL+BR-RR)*np.log(refinementFactor))/(n-1)
            refined = normal / refinementFactor

            arr[0] = 0.0
            for i in range(1, n):
                prev = arr[i-1]
                if prev < BL:
                    arr[i] = prev + normal
                elif prev < RL:
                    frac = (prev - BL) / (RL - BL)
                    arr[i] = prev + normal*(1 - frac) + refined*frac
                elif prev < RR:
                    arr[i] = prev + refined
                elif prev < BR:
                    frac = (prev - RR) / (BR - RR)
                    arr[i] = prev + refined*(1 - frac) + normal*frac
                else:
                    arr[i] = prev + normal

            alpha = length/arr[n-1]
            for i in range(1, n):
                if(arr[i]>RR): arr[i] = alpha*arr[i]

            d = np.diff(arr)

            return arr, d, True
        else:
            # Uniform grid
            arr = np.linspace(0, length, n)
            d = np.diff(arr)

            return arr, d, False

    x, dx, non_uniform_x = process_direction(n_x, length_x, refinementFactor_x, x_BL, x_BR, x_RL, x_RR, x_grid, x_file)
    y, dy, non_uniform_y = process_direction(n_y, length_y, refinementFactor_y, y_BL, y_BR, y_RL, y_RR, y_grid, y_file)
    z, dz, non_uniform_z = process_direction(n_z, length_z, refinementFactor_z, z_BL, z_BR, z_RL, z_RR, z_grid, z_file)

    non_uniform_grid = non_uniform_x or non_uniform_y or non_uniform_z
    return x, dx, y, dy, z, dz, non_uniform_grid


def write_points(filename, points):
    """Write points to file."""
    with open(filename, 'w') as f:
        f.write("{")
        for i, (x, y, z) in enumerate(points):
            f.write(f"({int(x)}, {int(y)}, {int(z)})\n")
        f.write("}")

def write_values(filename, values):
    """Write float values to file."""
    with open(filename, 'w') as f:
        f.write("{")
        for i, value in enumerate(values):
            f.write(f"({value:.4g})")
            if i < len(values) - 1:
                f.write("  \n")
        f.write("}")

def write_complex_field_full(filename, field_values):
    with open(filename, 'w') as f:
        f.write("{")
        for i, field in enumerate(field_values):
            f.write("[")
            for j in range(3):
                if j < len(field) and np.abs(field[j].imag - free) >= 1e-9:
                    f.write(f"Complex({field[j].real:.4g},{field[j].imag:.4g})")
                else:
                    f.write("F")
                
                if j < 2:  # Add semicolon between elements except after the last one
                    f.write(";")
            f.write("]")

            if i < len(field_values) - 1:
                f.write("  \n")
        f.write("}")

def write_excitation_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, name_files):
    points = []
    field_values = []
    
    if flag_numerical:
        Sx_array, Sy_array, Sz_array = mms_source_all_equalE_numerical(x_grid, y_grid, z_grid, sigma0, x0, y0, z0, s_gaussian, omega, eps0, k0, alpha)
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    if i == 0 or i == Nx-1 or j == 0 or j == Ny-1: 
                        field_values.append((0, 0, 0))
                        points.append((i, j, k))
                        continue

                    if not flag_2D and (k == 0 or k == Nz-1):
                        field_values.append((0, 0, 0))
                        points.append((i, j, k))
                        continue

                    Sx = Sx_array[i,j,k]
                    Sy = Sy_array[i,j,k]
                    Sz = Sz_array[i,j,k]

                    field_values.append((Sx, Sy, Sz))
                    points.append((i, j, k))
    else: 
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    if i == 0 or i == Nx-1 or j == 0 or j == Ny-1: 
                        field_values.append((0, 0, 0))
                        points.append((i, j, k))
                        continue

                    if not flag_2D and (k == 0 or k == Nz-1):
                        field_values.append((0, 0, 0))
                        points.append((i, j, k))
                        continue

                    x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                    if flag_2D:
                        Sx, Sy, Sz = mms_source_2D(x, y, z, sigma0, x0, y0, z0, s_gaussian, omega, eps0, k0, alpha)
                    else:
                        Sx, Sy, Sz = mms_source_all_equalE(x, y, z, sigma0, x0, y0, z0, s_gaussian, omega, eps0, k0, alpha)
                    
                    field_values.append((Sx, Sy, Sz))
                    points.append((i, j, k))

            
    write_complex_field_full(script_dir + "/../input/" + name_files + "_excitation_values.dat", field_values)
    write_points(script_dir + "/../input/" + name_files + "_excitation_points.dat", points)  

def write_elec_dens_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, name_files):
    """Generate and write electron density data files based on custom condition function."""
    points = []
    values = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]

                ne = ne_profile_gauss(x, y, z, ne_0, x0, y0, z0, s_gaussian)
                if ne > 1e11:
                    values.append(ne_profile_gauss(x, y, z, ne_0, x0, y0, z0, s_gaussian))
                    points.append((i, j, k))
    
    write_values(script_dir + "/../input/" + name_files + "_elec_dens_values.dat", values)
    write_points(script_dir + "/../input/" + name_files + "_elec_dens_points.dat", points)


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
    

flag_2D = 1
flag_numerical      = 0
flag_write_input    = 1
flag_run_harps      = 1
flag_plot_results   = True

if __name__ == "__main__":     
    # parameters
    f = 2.45e9             # Hz
    ne_0 = 0*8e17          # m^-3
    x0, y0, z0 = 0.2, 0.2, 0.2
    L_BOX = 0.4
    s_gaussian = 0.1  

    sigma0 = e_charge*ne_0*(1-2.0j)
    omega = 2.0*np.pi*f
    k0 = omega / c0
    alpha = 2*np.pi*1.5/L_BOX
    

    # grid points

    if flag_2D:
        N = 256;        Nx = N;    Ny = N
        Nz = 1
        z_grid = np.array([L_BOX/2])
        z0 = L_BOX/2
    else:
        N = 64;        Nx = N;    Ny = N;    Nz = N
        z_grid = None

    x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(Nx, Ny, Nz, L_BOX, L_BOX, L_BOX, -2, -2, -2, 0.13, 0.17, 0.23, 0.27,
                                                        0.13, 0.17, 0.23, 0.27, 0.13, 0.17, 0.23, 0.27, x_grid=None, y_grid=None, z_grid=z_grid)

    if flag_write_input:
        print("Writing input files...")
        write_excitation_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, "MMS")   # source term for MMS
        write_elec_dens_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, "MMS") # Conductivity profile used

    if flag_run_harps:
        print("Running HARPS...")
        if flag_2D:
            subprocess.run('mpirun --bind-to socket -np 1 ./harps.exe input/MMS_2D.in', shell=True, capture_output=True, text=True)
        else:
            subprocess.run('mpirun --bind-to socket -np 64 ./harps.exe input/MMS.in', shell=True, capture_output=True, text=True)


    if flag_plot_results:
        X, Y = np.meshgrid(x_grid, y_grid)

        if flag_2D:
            Z = np.full_like(X, L_BOX/2)

            Ex, Ey, Ez = E_field_gauss_trig_2D(X,Y,Z, alpha)

            if not flag_numerical:
                Sx, Sy, Sz = mms_source_2D(X, Y, Z, sigma0, x0, y0, z0, s_gaussian, omega, eps0, k0, alpha)

            # plot E field
            plt.figure(figsize=(10,4))
            plt.subplot(1,2,1)
            plt.contourf(X, Y, np.abs(Ez), 40, cmap='RdBu_r')
            #plt.scatter(X, Y, c=np.abs(Ex), s=6, cmap='RdBu_r')
            plt.colorbar()
            plt.title('E_x field (V/m)')
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')

            # plot S field
            plt.subplot(1,2,2)
            plt.contourf(X, Y, np.abs(Sz), 40, cmap='RdBu_r')
            #plt.scatter(X, Y, c=np.abs(Sx), s=6, cmap='RdBu_r')
            plt.colorbar()
            plt.title('S_z field (V/m³)')
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')

            plt.tight_layout()
            plt.savefig('analysis/analysis_output/MMS_E_S_fields.png', dpi=300)

            if flag_run_harps:
                data_x, grid_x, grid_y, grid_z = load_3d_data("Outputs/E_real.txt")

                X, Y = np.meshgrid(x_grid, y_grid)
                E_z_slice = data_x[:, :, 0]

                plt.figure(figsize=(5,4))
                plt.contourf(X, Y, np.abs(E_z_slice.T), 40, cmap='RdBu_r')
                plt.colorbar()
                plt.title('HARPS E_z field (V/m)')
                plt.xlabel('x (m)')
                plt.ylabel('y (m)')
                plt.tight_layout()
                plt.savefig('analysis/analysis_output/MMS_HARPS_Ez_field.png', dpi=300)
                
                
                error = np.abs(np.abs(Ez) - E_z_slice.T)
                plt.figure(figsize=(5,4))
                plt.contourf(X, Y, error, 40, cmap='RdBu_r')
                plt.colorbar()
                plt.title('Error in E_z field (V/m)')
                plt.xlabel('x (m)')
                plt.ylabel('y (m)')
                plt.tight_layout()
                plt.savefig('analysis/analysis_output/MMS_HARPS_Ez_field_error.png', dpi=300)

                print(f"norm of error in Ez field: {np.linalg.norm(error)/np.linalg.norm(Ez)}")
        else:
            Z = np.full_like(X, L_BOX*5/12)  # create 2D array filled with z_slice value
            z_slice = L_BOX*5/12
            z_index = (np.abs(z_grid - z_slice)).argmin()

            Ex, Ey, Ez = E_field_gauss_trig(X,Y,Z, alpha)
            if flag_numerical:
                Sx, Sy, Sz = mms_source_all_equalE_numerical(x_grid, y_grid, z_grid, sigma0, x0, y0, z0, s_gaussian, omega, eps0, k0, alpha)
            else:
                Sx, Sy, Sz = mms_source_all_equalE(X, Y, Z, sigma0, x0, y0, z0, s_gaussian, omega, eps0, k0, alpha)
            
            # plot E field
            plt.figure(figsize=(10,4))
            plt.subplot(1,2,1)
            plt.contourf(X, Y, np.abs(Ex), 40, cmap='RdBu_r')
            #plt.scatter(X, Y, c=np.abs(Ex), s=6, cmap='RdBu_r')
            plt.colorbar()
            plt.title('E_x field (V/m)')
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')

            # plot S field
            #if flag_numerical: np.save(f"./analysis/analysis_output/Sx.npy", Sx)
            #else: Sx_old = np.load(f"./analysis/analysis_output/Sx.npy")

            plt.subplot(1,2,2)
            if flag_numerical: plt.contourf(X, Y, np.abs(Sx[:, :, z_index]), 40, cmap='RdBu_r')
            else: plt.contourf(X, Y, np.abs(Sx), 40, cmap='RdBu_r')
            #plt.scatter(X, Y, c=np.abs(Sx), s=6, cmap='RdBu_r')
            plt.colorbar()
            plt.title('Error Analytical - Numerical (V/m³)')
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')

            plt.tight_layout()
            plt.savefig('analysis/analysis_output/MMS_E_S_fields.png', dpi=300)

            if flag_run_harps:
                #plotting Ex on z =L_BOX/4 slice from HARPS output
                data_x, grid_x, grid_y, grid_z = load_3d_data("Outputs/Ex_amp.txt")
                
                X, Y = np.meshgrid(x_grid, y_grid)
                z_slice = L_BOX*5/12
                z_index = (np.abs(grid_z - z_slice)).argmin()
                E_x_slice = data_x[:, :, z_index]

                plt.figure(figsize=(5,4))
                plt.contourf(X, Y, np.abs(E_x_slice.T), 40, cmap='RdBu_r')
                plt.colorbar()
                plt.title('HARPS E_x field (V/m)')
                plt.xlabel('x (m)')
                plt.ylabel('y (m)')
                plt.tight_layout()
                plt.savefig('analysis/analysis_output/MMS_HARPS_Ex_field.png', dpi=300)

                # plotting Ex on y =L_BOX/4 slice from HARPS output
                y_slice = L_BOX*5/12
                y_index = (np.abs(grid_y - y_slice)).argmin()
                E_x_slice_y = data_x[:, y_index, :]
                X_z, Z_z = np.meshgrid(grid_x, grid_z)

                plt.figure(figsize=(5,4))
                plt.contourf(X_z, Z_z, np.abs(E_x_slice_y.T), 40, cmap='RdBu_r')
                plt.colorbar()  
                plt.title('HARPS E_x field (V/m)')
                plt.xlabel('x (m)')
                plt.ylabel('z (m)')
                plt.tight_layout()
                plt.savefig('analysis/analysis_output/MMS_HARPS_Ex_field_y_slice.png', dpi=300)

                # compute and plot error
                Z = np.full_like(X, z_grid[z_index])
                Ex, Ey, Ez = E_field_gauss_trig(X, Y, Z, alpha)
                error = np.abs(np.abs(Ex) - E_x_slice.T)
                plt.figure(figsize=(5,4))
                #plt.contourf(X, Y, error, 40, cmap='RdBu_r')
                plt.scatter(X, Y,  error, s=6, cmap='RdBu_r')
                plt.colorbar()  
                plt.title('Error in E_x field (V/m)')
                plt.xlabel('x (m)')
                plt.ylabel('y (m)')
                plt.tight_layout()
                plt.savefig('analysis/analysis_output/MMS_HARPS_Ex_field_error.png', dpi=300)