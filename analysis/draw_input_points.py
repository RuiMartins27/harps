import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re

def read_points_from_file(filename):
    with open(filename, 'r') as file:
        data = file.read()
    
    # Extract tuples using regex
    points = re.findall(r'\((\d+),\s*(\d+),\s*(\d+)\)', data)
    
    # Convert extracted strings to integers
    return [(int(x), int(y), int(z)) for x, y, z in points]

def plot_3d_points(points, output_filename="output_plot.png"):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    x_vals = [p[0] for p in points]
    y_vals = [p[1] for p in points]
    z_vals = [p[2] for p in points]
    
    ax.scatter(x_vals, y_vals, z_vals, c='b', marker='o')
    
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')

    ax.set_xlim(0, 64)
    ax.set_ylim(0, 64)
    ax.set_zlim(0, 192)
    
    # Save the plot to the specified file
    plt.savefig(output_filename)
    print(f"Plot saved to {output_filename}")

if __name__ == "__main__":
    #filename = "input/example16_excitation_points.dat"
    #filename = "input/example16_metal_points.dat"
    #filename = "input/example16_permittivity_points.dat"
    filename = "input/example16_elec_dens_points.dat"
    points = read_points_from_file(filename)
    
    # Call the function with a desired output filename
    plot_3d_points(points, output_filename="analysis/analysis_output/graph_input_points.png")
