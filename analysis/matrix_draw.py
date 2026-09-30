import numpy as np
import matplotlib.pyplot as plt
import re

def visualize_large_matrix_sparsity(filename):
    try:
        with open(filename, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return
    
    # Extract matrix size from the first line
    size_pattern = re.compile(r'Matrix size: (\d+) x (\d+)')
    size_match = size_pattern.search(lines[0])
    if not size_match:
        print("Error: Could not find matrix size in the file.")
        return
    
    rows_size = int(size_match.group(1))
    cols_size = int(size_match.group(2))
    
    # Extract nonzero elements
    row_indices = []
    col_indices = []
    
    element_pattern = re.compile(r'^(\d+),\s*(\d+),\s*\((-?[\d\.]+),\s*(-?[\d\.]+)\)')
    
    for line in lines[2:]:  # Skip header lines
        match = element_pattern.match(line)
        if match:
            row_indices.append(int(match.group(1)))
            col_indices.append(int(match.group(2)))
    
    # Plot sparsity pattern
    plt.figure(figsize=(100, 100))
    plt.plot(col_indices, row_indices, marker=',', linestyle='None', color='blue', markersize=0.01)
    plt.gca().invert_yaxis()
    plt.xlim(-1, cols_size)
    plt.ylim(rows_size, -1)
    plt.title(f'Sparsity Pattern of {rows_size}x{cols_size} Matrix')
    plt.xlabel('Column Index')
    plt.ylabel('Row Index')
    
    actual_nonzeros = len(row_indices)
    info_text = f'Non-zeros found: {actual_nonzeros}'
    plt.text(0.5, -0.05, info_text, 
             horizontalalignment='center', verticalalignment='center', 
             transform=plt.gca().transAxes)
    
    plt.tight_layout()
    
    # Save as vector graphics for infinite resolution
    plt.savefig('analysis/analysis_output/large_matrix_sparsity.svg', format='svg')
    plt.savefig('analysis/analysis_output/large_matrix_sparsity.pdf', format='pdf')
    
    plt.show()

# Example usage
visualize_large_matrix_sparsity('Outputs/matrix.txt')
