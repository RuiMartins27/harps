import numpy as np
import os
import argparse
import matplotlib.pyplot as plt

harps_dir = os.path.dirname(os.path.abspath(__file__)) + "/../"
script_dir = os.path.dirname(os.path.abspath(__file__))

free = -0.00123456

def write_complex_field(filename, field_values, component):
    """Write complex field values to file."""
    with open(filename, 'w') as f:
        f.write("{")
        for i, field in enumerate(field_values):
            if component == 0: f.write(f"[Complex({field.real:.4g},{field.imag:.4g});F;F]\n")
            if component == 1: f.write(f"[F;Complex({field.real:.4g},{field.imag:.4g});F]\n")
            if component == 2: f.write(f"[F;F;Complex({field.real:.4g},{field.imag:.4g})]\n")

            if i < len(field_values) - 1:
                f.write("  ")
        f.write("}")

def write_complex_field_full(filename, field_values):
    """Write complex field values to file."""
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
            f.write("]\n")

            if i < len(field_values) - 1:
                f.write("  \n")
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

def write_points(filename, points):
    """Write points to file."""
    with open(filename, 'w') as f:
        f.write("{")
        for i, (x, y, z) in enumerate(points):
            f.write(f"({int(x)}, {int(y)}, {int(z)})\n")
        f.write("}")



def write_permittivity_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set_permittivity, parameters, name_files):
    """Generate and write permittivity data files based on custom condition function."""
    points = []
    values = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                function_set_permittivity(x, y, z, i, j, k, points, values, parameters)
    
    write_values(script_dir + "/../input/" + name_files + "_permittivity_values.dat", values)
    write_points(script_dir + "/../input/" + name_files + "_permittivity_points.dat", points)

def write_elec_dens_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set_elec_dens, parameters, name_files):
    """Generate and write electron density data files based on custom condition function."""
    points = []
    values = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                function_set_elec_dens(x, y, z, i, j, k, points, values, parameters)
    
    write_values(script_dir + "/../input/" + name_files + "_elec_dens_values.dat", values)
    write_points(script_dir + "/../input/" + name_files + "_elec_dens_points.dat", points)

def write_mu_real_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set_elec_dens, parameters, name_files):
    """Generate and write electron density data files based on custom condition function."""
    points = []
    values = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                function_set_elec_dens(x, y, z, i, j, k, points, values, parameters)
    
    write_values(script_dir + "/../input/" + name_files + "_mu_real_values.dat", values)
    write_points(script_dir + "/../input/" + name_files + "_mu_real_points.dat", points)

def write_mu_imag_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set_elec_dens, parameters, name_files):
    """Generate and write electron density data files based on custom condition function."""
    points = []
    values = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                function_set_elec_dens(x, y, z, i, j, k, points, values, parameters)
    
    write_values(script_dir + "/../input/" + name_files + "_mu_imag_values.dat", values)
    write_points(script_dir + "/../input/" + name_files + "_mu_imag_points.dat", points)

def write_excitation_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set_excitation, parameters, name_files):
    """Generate grid points and field values."""
    points = []
    field_values = []

    component = parameters[0]
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                fields = function_set_excitation(x, y, z, parameters)
                if fields == None: continue
                if j > 0: continue
                
                field_values.append(fields)

                points.append((i, j, k))

            
    if(component < 3): write_complex_field(script_dir + "/../input/" + name_files + "_excitation_values.dat", field_values, component)
    else: write_complex_field_full(script_dir + "/../input/" + name_files + "_excitation_values.dat", field_values)
    write_points(script_dir + "/../input/" + name_files + "_excitation_points.dat", points)

def write_metal_data(Nx, Ny, Nz, x_grid, y_grid, z_grid, function_set_metal, parameters, name_files):
    """Generate and write permittivity data files based on custom condition function."""
    points = []
    
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                x = x_grid[i]; y = y_grid[j]; z = z_grid[k]
                function_set_metal(x, y, z, i, j, k, points, parameters)
    
    write_points(script_dir + "/../input/" + name_files + "_metal_points.dat", points)


def gaussian_field(x, y, z, parameters):
    b_x = parameters[1]
    b_y = parameters[2]
    b_z = parameters[3]
    sigma = parameters[4]

    return complex(np.exp(-((x-b_x)*(x-b_x)+(y-b_y)*(y-b_y)+(z-b_z)*(z-b_z))/(2*sigma*sigma)),0)

def square_waveguide_excitation(x, y, z, parameters):
    a = parameters[1]
    alpha = parameters[2]
    y_exc = parameters[3]
    x_min = parameters[4]
    x_max = parameters[5]
    x_offset = parameters[6]

    if (abs(y-y_exc) > 5e-4 or x<x_min or x>x_max): return None

    #return [complex(-alpha*np.sin(np.pi*x/a), 0), complex(0, free), complex(0, free)]
    #return [complex(0, free), complex(-alpha*np.sin(np.pi*x/a), 0), complex(0, free)]

    return [complex(0, free), complex(0, free), complex(-alpha*np.sin(np.pi*(x+x_offset)/a), 0)]

def constant_excitation(x, y, z, parameters):
    if (abs(y-parameters[2])<4e-3 and x>parameters[3] and x<parameters[4]): alpha = parameters[1]
    else: return

    return complex(alpha,0)

def square_waveguide_excitation_3D(x, y, z, parameters):
    alpha = parameters[1]
    a = parameters[2]
    b = parameters[3]
    center_z = parameters[4]
    y_exc = parameters[5]
    x_min = parameters[6]
    x_max = parameters[7]

    if (abs(y-y_exc) > 2e-3 or z < (center_z - b/2) or z > (center_z + b/2) or x<x_min or x>x_max): return None

    return [complex(0, free), complex(0, free), complex(-alpha*np.sin(np.pi*x/a), 0)]

def square_waveguide_excitation_3D_BC(x, y, z, parameters):
    alpha = parameters[1]
    a = parameters[2]
    x_0 = parameters[3]
    center_z = parameters[4]
    b = parameters[5]


    if z < (center_z - b/2) or z > (center_z + b/2): return [complex(0, 0), complex(0, 0), complex(0, 0) ]
    
    if (x_0 < 0):
        if (z > 1e-6): return None
        return [complex(0, 0), complex(-alpha*np.sin(np.pi*(x+x_0)/a), 0), complex(0, 0) ]
    
    if (y > 1e-6): return None

    return [complex(0, 0), complex(0, 0), complex(-alpha*np.sin(np.pi*(x+x_0)/a), 0)]


def diagonal_interface(x, y, z, i, j, k, points, values, param):
    epsilon_r = param[0]
    diag = param[1]

    if(x + y > diag):
        points.append((int(i), int(j), int(k)))
        values.append(epsilon_r)

def ring_interface_3D(x, y, z, i, j, k, points, values, param):
    epsilon_r = param[0]
    R_inner = param[1]
    R_outer = param[2]
    center_x = param[3]
    center_y = param[4]
    rotated = param[5]

    if rotated: radius = np.sqrt((z-center_x)*(z-center_x) + (y-center_y)*(y-center_y))
    else: radius = np.sqrt((x-center_x)*(x-center_x) + (y-center_y)*(y-center_y))


    if(radius > R_inner and radius < R_outer):
        points.append((int(i), int(j), int(k)))
        values.append(epsilon_r)

def gaussian_profile(x, y, z, i, j, k, points, values, param):
    ne_0 = param[0]
    center_x = param[1]
    center_y = param[2]
    center_z = param[3]
    sigma_r = param[4]
    sigma_z = param[5]

    radius = np.sqrt((x-center_x)*(x-center_x) + (y-center_y)*(y-center_y))

    if(radius<0.014):
        if(sigma_z > 0):
            ne = ne_0*np.exp(-0.5*(radius/sigma_r)**2)*np.exp(-0.5*(np.abs(z-center_z)/sigma_z)**2)
        else:
            ne = ne_0*np.exp(-0.5*(radius/sigma_r)**2)
    else:
        ne = 0.9e11

    if(ne >= 1e11):
        points.append((int(i), int(j), int(k)))
        values.append(ne)

def mod_gaussian_profile(x, y, z, i, j, k, points, values, param):
    ne_0 = param[0]
    center_x = param[1]
    center_y = param[2]
    center_z = param[3]
    sigma_x = param[4]
    sigma_y = param[5]
    sigma_z = param[6]
    power = param[7]

    if(abs(x-center_x) < sigma_x*3 and abs(y-center_y) < sigma_y*3 and abs(z-center_z) < sigma_z*power):
        ne = ne_0*np.exp(-0.5*(abs(x-center_x)/sigma_x)**power-0.5*(abs(y-center_y)/sigma_y)**power-0.5*(abs(z-center_z)/sigma_z)**power)
    else:
        ne = 0.9e11
    

    if(ne >= 1e11):
        points.append((int(i), int(j), int(k)))
        values.append(ne)

def turn_1D_into_grid(x, y, z, i, j, k, points, values, param):
    radius = param[0]
    values_r = param[1]
    center_x = param[2]
    center_y = param[3]
    R_in = param[4]

    # linearly interpolate the values to the grid
    r = np.sqrt((x - center_x)*(x - center_x) + (y - center_y)*(y - center_y))
    value = np.interp(r, radius, values_r)

    if(r < R_in):
        points.append((int(i), int(j), int(k)))
        values.append(value)

def quartz_tube(x, y, z, i, j, k, points, values, parameters):
    epsilon_r = parameters[0]
    R_inner = parameters[1]
    R_outer = parameters[2]
    center_y = parameters[3]

    radius = np.abs(y-center_y)

    if(radius > R_inner and radius < R_outer):
        points.append((int(i), int(j), int(k)))
        values.append(epsilon_r)

def excitation_TE10_z(x, y, z, parameters):
    alpha = parameters[1]
    b = parameters[2]
    center_z = parameters[3]
    
    if (y > 1e-6): return None

    if (abs(z - center_z) < b/2): return [complex(0, 0), complex(0, 0), complex(-alpha, 0)]
    else: return [complex(0, 0), complex(0, 0), complex(0, 0)]

def waveguide_metal(x, y, z, i, j, k, points, param):
    center_y = param[0]
    center_z = param[1]
    w = param[2]
    t_wg = param[3]
    R_outer = param[4]
    y_reflector = param[5]

    radius_z = np.abs(center_z - z)
    radius_y = np.abs(center_y - y)

    if(radius_y > R_outer and (radius_z > w/2 and radius_z < (w/2 + t_wg)) or (radius_z < w/2 and y > y_reflector)):
        points.append((int(i), int(j), int(k)))


def real_waveguide(x, y, z, i, j, k, points, param):
    center_x = param[0]
    center_y = param[1]
    center_z = param[2]
    R_outer = param[3]
    b = param[4]


    radius = np.sqrt((x-center_x)*(x-center_x) + (y-center_y)*(y-center_y))

    #if(radius > R_outer and ((abs(z - (center_z - b/2)) < 1e-3)  or (abs(z - (center_z + b/2)) < 1e-3))):

    if(radius > R_outer and (z < (center_z - b/2)  or z > (center_z + b/2) )):
        points.append((int(i), int(j), int(k)))

def reflector(x, y, z, i, j, k, points, param):
    center_z = param[0]
    b = param[1]
    y_reflect = param[2]


    if(y >= y_reflect and z > (center_z - b/2) and z < (center_z + b/2)):
        points.append((int(i), int(j), int(k)))

def reflector_waveguide(x, y, z, i, j, k, points, param):
    center_x = param[0]
    center_y = param[1]
    center_z = param[2]
    b = param[3]
    y_reflect = param[4]
    t_WG = param[5]
    R_out = param[6]

    radius = np.sqrt((x-center_x)*(x-center_x) + (y-center_y)*(y-center_y))

    if(y >= y_reflect):
        points.append((int(i), int(j), int(k)))
        return
    elif radius > R_out and ((z > (center_z - b/2 - t_WG) and z < (center_z - b/2)) or ( z > (center_z + b/2) and z < (center_z + b/2 + t_WG))):
        points.append((int(i), int(j), int(k)))


def three_stub_tuner(x, y, z, i, j, k, points, param):
    d = param[0]
    center_x = param[1]
    y_1 = param[2]
    y_2 = param[3]
    y_3 = param[4]
    L_1 = param[5]
    L_2 = param[6]
    L_3 = param[7]

    d_2_squared = d*d/4

    radius_1_squared = (x-center_x)*(x-center_x) + (y-y_1)*(y-y_1)
    radius_2_squared = (x-center_x)*(x-center_x) + (y-y_2)*(y-y_2)
    radius_3_squared = (x-center_x)*(x-center_x) + (y-y_3)*(y-y_3)

    if(radius_1_squared < d_2_squared and z < L_1):
        points.append((int(i), int(j), int(k)))
    elif(radius_2_squared < d_2_squared and z < L_2):
        points.append((int(i), int(j), int(k)))
    elif(radius_3_squared < d_2_squared and z < L_3):
        points.append((int(i), int(j), int(k)))


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
        
        elif refinementFactor != -1:
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

def read_radial_profile(filename):
    # Skip the header and read all but the last line
    with open(filename) as f:
        total_lines = sum(1 for _ in f)
    data = np.loadtxt(filename, skiprows=1, max_rows=total_lines - 2)
    radius = data[:, 0]      # First column: Radius (m)
    values = data[:, 1]      # Second column: Absorbed Power Density (W/m³)
    return radius, values

from scipy.interpolate import RegularGridInterpolator

def load_harps_2d_field(filename, fill_val=0.0):
    with open(filename, 'r') as f:
        # Remove empty lines
        lines = [line.strip() for line in f.readlines() if line.strip()]
    
    # Line 0: Nx Ny Nz Lx Ly Lz
    dims = lines[0].split()
    ny, nz = int(dims[1]), int(dims[2])
    
    # Line 1 is GRID_TYPE (e.g., NON_UNIFORM_GRID)
    # Line 2 is X grid (Nx=1 point, so we skip using it for interpolation)
    
    # Line 3: Y grid (Radial/Transverse axis)
    y_grid = np.array(lines[3].split(), dtype=float)
    
    # Line 4: Z grid (Axial axis)
    z_grid = np.array(lines[4].split(), dtype=float)
    
    # Line 5+: The flattened block of field data
    data_str = " ".join(lines[5:])
    data = np.array(data_str.split(), dtype=float).reshape((ny, nz))
    
    # Create the 2D interpolator mapping (y, z) -> value.
    # We use fill_val for 3D corners or lengths that extend beyond the 2D domain.
    interp = RegularGridInterpolator((y_grid, z_grid), data, 
                                     bounds_error=False, fill_value=fill_val)
    return interp

def map_2d_to_3d_cylindrical(x, y, z, i, j, k, points, values, parameters):
    interp = parameters[0]
    cx = parameters[1]
    cy = parameters[2]
    z_offset = parameters[3]  # New parameter: shift in Z between 3D and 2D
    
    # Calculate radial distance
    r = np.sqrt((x - cx)**2 + (y - cy)**2)
    
    # Map back to 2D coordinates
    y_2d = cy - r 
    z_2d = z - z_offset       # Shift the Z coordinate to match the 2D center
    
    # Evaluate interpolator
    val = interp((y_2d, z_2d)).item()
    
    points.append((i, j, k))
    values.append(val)

def map_2d_to_2d_direct(x, y, z, i, j, k, points, values, parameters):
    """
    Directly maps a 2D coordinate to a 2D interpolator without cylindrical radial math.
    """
    interp = parameters[0]
    
    # In 2D, y and z map directly.
    val = interp((y, z)).item()
    
    # Append grid indices as a tuple
    points.append((i, j, k))
    values.append(val)
    
#Usage
parser = argparse.ArgumentParser(description="Writer for input files")
parser.add_argument('problem', nargs='?', default=8, type=int, help="Problem index (0-10), default is 8")

# Next three arguments: z1, z2, z3 (real numbers)
parser.add_argument( 'z1', nargs='?', default=0.119, type=float, help="First real number (tuner param 1 or replector position in meters)" )
parser.add_argument( 'z2', nargs='?', default=0.25, type=float, help="Second real number (tuner param 2)" )
parser.add_argument( 'z3', nargs='?', default=0.25, type=float, help="Third real number (tuner param 3)" )

args = parser.parse_args()
problem = args.problem
z1 = args.z1; z2 = args.z2; z3 = args.z3

if(problem == 1):
    N_x = 100; N_y = 1; N_z = 1600
    a = 0.08636; b = 0.04318
    L_x = a; L_y = b; L_z = 1
     
    symmetric = False
    if(symmetric):
        x_offset = L_x/2; N_x = N_x//2; L_x = L_x/2; center_x = 0
        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, -1, -1, x_grid=None, y_grid=None, z_grid=None)

    else:
        x_offset = 0
        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, -1, -1,  x_grid=None, y_grid=None, z_grid=None)

    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,21000, a, -1e-15, L_z/2, 1e5*b],"example1")

if(problem == 3):
    N_x = 256; N_y = N_x; N_z = 1
    L_x = 0.9; L_y = L_x; L_z = 0.1
    uniform_grid = True

    x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, -1, -1)

    write_excitation_data(N_x,1,1, x_grid, y_grid, z_grid, gaussian_field,[2,L_x/2,0,0,0.03],"example3")
    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, diagonal_interface, [4.5,L_x], "example3")

if(problem == 4):
    N_x = 256; N_y = N_x; N_z = 1
    L_x = 0.8; L_y = L_x; L_z = 0.1

    x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, -1, -1)

    write_elec_dens_data(N_x, N_y, 1, x_grid, y_grid, z_grid, gaussian_profile, [5e18,80,L_x/2,80,L_y/2,0,0], "example4")
    write_excitation_data(N_x, N_y, 1, x_grid, y_grid, z_grid, constant_excitation,[1,1,0.05,0.02,0.78],"example4")

if(problem == 6):
    N_x = 256; N_y = 256; N_z = 1
    a = 0.08636; b = 0.04318; R_in = 0.0135; R_out = 0.015; y_r = z1
    L_x = a; L_y = 0.1400; L_z = b; center_x = a/2; center_y = 0.046; center_z = b/2


    symmetric = True
    if(symmetric):
        x_offset = L_x/2; N_x = N_x//2; L_x = L_x/2; center_x = 0

        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 4, 4, -1,
            0.000, 0.00, 0.008, 0.012, 0.034, 0.038, 0.054, 0.058, x_grid=None, y_grid=None, z_grid=None)

    else:
        x_offset = 0
        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, -1, -1,
            0.000, 0.00, 0.009, 0.015, 0.031, 0.037, 0.055, 0.061, x_grid=None, y_grid=None, z_grid=None)

    write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile, [4.5e18, center_x, center_y, center_z, 0.0035, 0.0035, 1e9, 3], "example6")
    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,ring_interface_3D, [4.5, R_in, R_out, center_x, center_y, 0], "example6")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,21000, a, x_offset, L_z/2, 1e5*b],"example6")

    if(z2 > 0.8): write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,     reflector, [center_z, 1e5*b, y_r], "example6" + str(int(z2)))
    else:  write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,     reflector, [center_z, 1e5*b, y_r], "example6")
    
    #write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation,[7, a, 21000, 0.01, L_x/50, L_x*49/50, x_offset],"example6")

if(problem == 7):
    N_x = 64; N_y = N_x; N_z = N_x
    a = 0.08636; b = 0.04318; R_in = 0.014; R_out = 0.015; 
    L_x = a; L_y = 0.122; L_z = 0.075; center_x = a/2; center_y = L_y/2; center_z = 0.075/2

    symmetric = False
    if(symmetric):
        x_offset = L_x/2; N_x = N_x//2; L_x = L_x/2; center_x = 0

        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 6, 6, -1,
            0.0, 0.0, 0.013, 0.019, 0.046, 0.076, 0.086, 0.015, 0.022, 0.043, 0.050, x_grid=None, y_grid=None, z_grid=None)

    else:
        x_offset = 0
        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 6, 6, -1,
            0.023, 0.031, 0.057, 0.063, 0.046, 0.076, 0.086, 0.015, 0.022, 0.043, 0.050, x_grid=None, y_grid=None, z_grid=None)

    write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile,   [0.5e18, center_x, center_y, center_z, 0.005, 0.005, 0.011, 3], "example7")
    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,ring_interface_3D,      [4.5, R_in, R_out, center_x, center_y, 0], "example7")
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,       real_waveguide,         [center_x, center_y, center_z, R_out, b], "example7")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,21000, a, x_offset, L_z/2, 1e5*b],"example7")
     
if(problem == 8):
    N_x = 64; N_y = 160; N_z = 64
    a = 0.08636; b = 0.04318; R_in = 0.014; R_out = 0.015
    L_x = a; L_y = 0.36; L_z = b;  center_x = a/2; center_y = 0.12*2.5; center_z = b/2
     
    symmetric = False
    if(symmetric):
        x_offset = L_x/2; N_x = N_x//2; L_x = L_x/2; center_x = 0

        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 6, 3, -1,
            0.0, 0.0, 0.013, 0.019, 0.12, 0.286, 0.315, 0.325, x_grid=None, y_grid=None, z_grid=None)

    else:
        x_offset = 0
        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 6, 3, -1,
            0.023, 0.031, 0.057, 0.063, 0.12, 0.286, 0.315, 0.325, x_grid=None, y_grid=None, z_grid=None)

    
    write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile,   [2.5e18, center_x, center_y, center_z, 0.0043, 0.0043, 0.011, 3], "example8")
    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,ring_interface_3D,      [4.5, R_in, R_out, center_x, center_y, 0], "example8")
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,       three_stub_tuner,        [a/10, center_x, 0.12*1.1, 0.12*1.55, 0.12*2, L_z*z1, L_z*z2, L_z*z3], "example8")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,21000, a, x_offset, center_z, b],"example8")

if(problem == 9):
    N_x = 64; N_y = N_x; N_z = 2*N_x
    a = 0.08636; b = 0.04318; R_in = 0.014; R_out = 0.015; 
    L_x = a; L_y = 0.0875; L_z = b; center_x = a/2; center_y = 0.046; center_z = L_z/2

    symmetric = False
    if(symmetric):
        x_offset = L_x/2; N_x = N_x//2; L_x = L_x/2; center_x = 0

        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 8, 8, -1,
            0.0, 0.0, 0.014, 0.024, 0.022, 0.032, 0.060, 0.070, x_grid=None, y_grid=None, z_grid=None)

    else:
        x_offset = 0
        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 8, 8, -1,
            0.023, 0.031, 0.055, 0.063, 0.026, 0.034, 0.058, 0.066, 0.0016, 0.0096, 0.0336, 0.0416, x_grid=None, y_grid=None, z_grid=None)

    write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile,   [6.25e18, center_x, center_y, center_z, 0.004, 0.004, 0.03, 3], "example9")
    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,ring_interface_3D,      [4.5, R_in, R_out, center_x, center_y, 0], "example9")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,21000, a, x_offset, center_z, b],"example9")

if(problem == 10):
    N_x = 256; N_y = 256; N_z = 1
    a = 0.08636; b = 0.04318; R_in = 0.014; R_out = 0.015; y_r = z1
    L_x = a; L_y = 0.1400; L_z = b; center_x = a/2; center_y = 0.046; center_z = b/2

    first_run = True
    symmetric = True

    if(symmetric):
        x_offset = L_x/2; N_x = N_x//2; L_x = L_x/2; center_x = 0

        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 4, 4, -1,
            0.000, 0.00, 0.013, 0.019, 0.026, 0.034, 0.058, 0.066, x_grid=None, y_grid=None, z_grid=None)

    else:
        x_offset = 0
        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 4, 4, -1,
            0.023, 0.031, 0.057, 0.063, 0.026, 0.034, 0.058, 0.066, x_grid=None, y_grid=None, z_grid=None)

    if(first_run):
        #write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile, [8e17, center_x, center_y, center_z, 0.0042, -1], "example6")
        write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,ring_interface_3D, [4.5, R_in, R_out, center_x, center_y], "example6")
        write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,21000, a, x_offset],"example6")
        write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,     reflector, [center_z, 1e5*b, y_r], "example6")

        #write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation,[7, a, 21000, 0.01, L_x/50, L_x*49/50, x_offset],"example6")
    else:
        radius, values = read_radial_profile("../1d_codes/results/elec_dens_radius.txt")
        write_elec_dens_data(N_x, N_y, 1, x_grid, y_grid, z_grid,   turn_1D_into_grid, [radius, values, center_x, center_y, R_in], "radial_1D")
        radius, values = read_radial_profile("../1d_codes/results/mu_real_radius.txt")
        write_mu_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, turn_1D_into_grid, [radius, values, center_x, center_y, R_in], "radial_1D")
        radius, values = read_radial_profile("../1d_codes/results/mu_imag_radius.txt")
        write_mu_imag_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, turn_1D_into_grid, [radius, values, center_x, center_y, R_in], "radial_1D")

if(problem == 11):
    N_x = 1; N_y = 384; N_z = 512
    a = 0.08636; b = 0.04318; R_in = 0.014; R_out = 0.015; y_r = z1
    L_x = 0.025; L_y = 0.185; L_z = b; center_y = 0.146; center_z = b/2; center_x = 0; x_offset = b

    x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, 4, -1,
            0.01, 0.02, 0.03, 0.04, 0.13, 0.136, 0.156, 0.162, x_grid=None, y_grid=None, z_grid=None)

    write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile,   [1.0e18, center_x, center_y, center_z, 0.0043, 0.009], "example11")
    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, ring_interface_3D, [4.5, R_in, R_out, center_x, center_y], "example11")
    #write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,     reflector, [center_z, 1e5*b, y_r], "example11")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,21000, a, x_offset, L_z/2, 1e5*b],"example11")

if(problem == 12):
    N_x = 1; N_y = 512; N_z = 512
    a = 0.08636; b = 0.04318; R_in = 0.014; R_out = 0.015; y_r = z1; thickness_wg = 0.02
    L_x = 0.02; L_y = 0.2400; L_z = 0.1200; center_y = 0.146; center_z = L_z/2; center_x = 0; x_offset = b

    x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, 4, -1,
            0.01, 0.02, 0.03, 0.04, 0.13, 0.136, 0.156, 0.162, x_grid=None, y_grid=None, z_grid=None)

    write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile,   [6.25e18, center_x, center_y, center_z, 0.002, 0.002, 0.01, 3], "example12")
    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, ring_interface_3D, [4.5, R_in, R_out, center_x, center_y, 0], "example12")
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,     reflector_waveguide, [center_x, center_y, center_z, b, y_r, thickness_wg, R_out], "example12")
    write_excitation_data(N_x, 1, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,21000, a, x_offset, L_z/2, b],"example12")


if(problem == 13):
    # Rotated waveguide in 3D (no tuner)
    N_x = 64*2; N_y = 64; N_z = 64
    a = 0.08636; b = 0.04318; R_in = 0.014; R_out = 0.015; 
    L_x = a; L_y = 0.0875; L_z = b; center_x = a/2; center_y = 0.046; center_z = L_z/2

    symmetric = False
    x_offset = 0
    x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, 8, 8,
        0.023, 0.031, 0.055, 0.063, 0.026, 0.034, 0.058, 0.066, 0.0016, 0.0096, 0.0336, 0.0416, x_grid=None, y_grid=None, z_grid=None)

    write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile,   [0.7e18, center_x, center_y, center_z, 0.0043, 0.0043, 0.0043, 3], "example13")
    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,ring_interface_3D,      [4.5, R_in, R_out, center_z, center_y, 1], "example13")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,21000, a, x_offset, center_z, b],"example13")

if(problem == 14):
    N_x = 256; N_y = 256; N_z = 1
    a = 0.123825*2; b = 0.123825; R_in = 0.0300; R_out = 0.0302; y_r = z1
    L_x = a; L_y = 0.3600; L_z = b; center_x = a/2; center_y = 0.115; center_z = b/2


    symmetric = True
    if(symmetric):
        x_offset = L_x/2; N_x = N_x//2; L_x = L_x/2; center_x = 0

        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 4, 4, -1,
            0.000, 0.00, 0.026, 0.036, 0.079, 0.089, 0.141, 0.151, x_grid=None, y_grid=None, z_grid=None)

    else:
        x_offset = 0
        x_grid, dx, y_grid, dy, z_grid, dz, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, -1, -1,
            0.000, 0.00, 0.013, 0.019, 0.121, 0.128, 0.153, 0.159, x_grid=None, y_grid=None, z_grid=None)

    write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile, [4.5e18, center_x, center_y, center_z, 0.0035, 0.0035, 1e9, 3], "example14")
    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,ring_interface_3D, [4.5, R_in, R_out, center_x, center_y, 0], "example14")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,7000, a, x_offset, L_z/2, 1e5*b],"example14")

    if(z2 > 0.8): write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,     reflector, [center_z, 1e5*b, y_r], "example6" + str(int(z2)))
    else:  write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,     reflector, [center_z, 1e5*b, y_r], "example14")


if(problem == 15):
    N_x = 1; N_y = 256; N_z = 384
    a = 0.08636; b = 0.04318; thickness_wg = 0.02; R_in = 0.0135; R_out = 0.015
    L_x = a; L_y = 0.16; L_z = 0.120
    center_y = 0.056; center_z = L_z/2 - 0.01

    x_grid, _, y_grid, _, z_grid, _, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, 4, -1,
                                                            0.01, 0.02, 0.03, 0.04, 0.04, 0.046, 0.066, 0.072, x_grid=None, y_grid=None, z_grid=None)
    
    load_2d_harps = True
    if(load_2d_harps):
        path_2d = "Outputs/" 
        
        # Load interpolators
        interp_elec = load_harps_2d_field(path_2d + "electron_density_2D_Pele.txt", fill_val=1e11)
        interp_mu_r = load_harps_2d_field(path_2d + "real_mobility_2D_Pele.txt", fill_val=0.0)
        interp_mu_i = load_harps_2d_field(path_2d + "imag_mobility_2D_Pele.txt", fill_val=0.0)

        # Write data using direct 2D mapping
        write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, map_2d_to_2d_direct, [interp_elec], "2D_RZ")
        write_mu_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, map_2d_to_2d_direct, [interp_mu_r], "2D_RZ")
        write_mu_imag_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, map_2d_to_2d_direct, [interp_mu_i], "2D_RZ")
    else:
        center_x = 0
        write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile,   [6.25e18, center_x, center_y, center_z, 0.002, 0.002, 0.01, 3], "2D_RZ")

    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, quartz_tube, [4.5, R_in, R_out, center_y], "2D_RZ")
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, waveguide_metal, [center_y, center_z, b, thickness_wg, R_out, 1e6], "2D_RZ")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  excitation_TE10_z,[7,21000, b, center_z],"2D_RZ")

if(problem == 16):
    N_x = 64; N_y = 64; N_z = 192
    a = 0.08636; b = 0.04318; thickness_wg = 0.02; R_in = 0.0135; R_out = 0.015
    L_x = 0.08636; L_y = 0.100; L_z = 0.100
    center_x = L_x/2; center_y = 0.056; center_z = L_z/2

    x_grid, _, y_grid, _, z_grid, _, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, 8, 8, -1,
                                                            0.02318, 0.03318, 0.05318, 0.06318, 0.036, 0.046, 0.066, 0.076, x_grid=None, y_grid=None, z_grid=None)

    load_2d_harps = True
    if(load_2d_harps):
        center_z_2d = 0.100 / 2
        z_offset = center_z - center_z_2d
        path_2d = "Outputs/" 
    
        interp_elec = load_harps_2d_field(path_2d + "electron_density_2D_Pele.txt", fill_val=1e11)
        interp_mu_r = load_harps_2d_field(path_2d + "real_mobility_2D_Pele.txt", fill_val=0.0)
        interp_mu_i = load_harps_2d_field(path_2d + "imag_mobility_2D_Pele.txt", fill_val=0.0)

        write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  map_2d_to_3d_cylindrical, [interp_elec, center_x, center_y, z_offset], "example16")                       
        write_mu_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, map_2d_to_3d_cylindrical, [interp_mu_r, center_x, center_y, z_offset], "example16")                  
        write_mu_imag_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, map_2d_to_3d_cylindrical, [interp_mu_i, center_x, center_y, z_offset], "example16")
    else:
        write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile,   [6.25e18, center_x, center_y, center_z, 0.004, 0.004, 0.03, 3], "example16")
    
    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, ring_interface_3D,      [4.5, R_in, R_out, center_x, center_y, 0], "example16")
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, waveguide_metal, [center_y, center_z, b, thickness_wg, R_out, 1], "example16")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  square_waveguide_excitation_3D_BC,[7,21000, a, 0, L_z/2, b],"example16")


if(problem == 17):
    N_x = 1; N_y = 256; N_z = 384
    a = 0.08636; b = 0.04318; thickness_wg = 0.02; R_in = 0.0130; R_out = 0.015
    L_x = a; L_y = 0.160; L_z = 0.120
    center_y = 0.056; center_z = L_z/2 - 0.01

    x_grid, _, y_grid, _, z_grid, _, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, 4, -1,
                                                            0.01, 0.02, 0.03, 0.04, 0.04, 0.046, 0.066, 0.072, x_grid=None, y_grid=None, z_grid=None)
    
    load_2d_harps = False
    if(load_2d_harps):
        path_2d = "Outputs/" 
        
        # Load interpolators
        interp_elec = load_harps_2d_field(path_2d + "electron_density_2D.txt", fill_val=1e11)
        interp_mu_r = load_harps_2d_field(path_2d + "real_mobility_2D.txt", fill_val=0.0)
        interp_mu_i = load_harps_2d_field(path_2d + "imag_mobility_2D.txt", fill_val=0.0)

        # Write data using direct 2D mapping
        write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, map_2d_to_2d_direct, [interp_elec], "example17")
        write_mu_real_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, map_2d_to_2d_direct, [interp_mu_r], "example17")
        write_mu_imag_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, map_2d_to_2d_direct, [interp_mu_i], "example17")
    else:
        center_x = 0
        write_elec_dens_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,   mod_gaussian_profile,   [6.25e18, center_x, center_y, center_z, 0.002, 0.002, 0.01, 3], "example17")

    write_permittivity_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, quartz_tube, [4.5, R_in, R_out, center_y], "example17")
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, waveguide_metal, [center_y, center_z, b, thickness_wg, R_out, z1], "example17")
    write_excitation_data(N_x, N_y, N_z, x_grid, y_grid, z_grid,  excitation_TE10_z,[7,21000, b, center_z],"example17")


if(problem == -17):
    N_x = 1; N_y = 256; N_z = 384
    a = 0.08636; b = 0.04318; thickness_wg = 0.02; R_in = 0.0130; R_out = 0.015
    L_x = a; L_y = 0.160; L_z = 0.120
    center_y = 0.056; center_z = L_z/2 - 0.01

    x_grid, _, y_grid, _, z_grid, _, is_non_uniform = generate_grid(N_x, N_y, N_z, L_x, L_y, L_z, -1, 4, -1,
                                                            0.01, 0.02, 0.03, 0.04, 0.04, 0.046, 0.066, 0.072, x_grid=None, y_grid=None, z_grid=None)
    
    write_metal_data(N_x, N_y, N_z, x_grid, y_grid, z_grid, waveguide_metal, [center_y, center_z, b, thickness_wg, R_out, z1], "example17")