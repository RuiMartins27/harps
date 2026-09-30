#ifndef UTILS_H
#define UTILS_H

#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
#include <cassert>

using Complex = std::complex<double>;


void printVecToFile(Vec v, std::string filename){
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename.c_str(), &viewer);
    VecView(v, viewer);
    PetscViewerDestroy(&viewer);
}

void printVecToFileBinary(Vec v, std::string filename){
    PetscViewer viewer;
    PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &viewer);
    VecView(v, viewer);
    PetscViewerDestroy(&viewer);
}

void printMatrixToFile(Mat A, std::string filename) {
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename.c_str(), &viewer);
    MatView(A, viewer);
    PetscViewerDestroy(&viewer);
}

void exportField2D(double* field, const std::string& filename, int n_x, int n_y) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Could not open file for writing field\n";
        return;
    }

    for (int i = 0; i < n_x; ++i) {
        for (int j = 0; j < n_y; ++j) {
            out << field[i * n_y + j] << " ";
        }
        out << "\n"; 
    }

    out.close();
}

void exportField3D(double* field, const std::string& filename, int n_x, int n_y, int n_z, double lengthX, double lengthY, double lengthZ,
                     bool non_uniform_grid, std::vector<double> x_grid, std::vector<double> y_grid, std::vector<double> z_grid) {

    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Could not open file for writing field\n";
        return;
    }
    
    // Write dimensions as header
    out << n_x << " " << n_y << " " << n_z << " " 
        << lengthX << " " << lengthY << " " << lengthZ << "\n";

    if(non_uniform_grid) {
        // Write grid points
        out << "NON_UNIFORM_GRID\n";
        for (int i = 0; i < n_x; ++i) out << x_grid[i] << " ";
        out << "\n";
        for (int j = 0; j < n_y; ++j) out << y_grid[j] << " ";
        out << "\n";
        for (int k = 0; k < n_z; ++k) out << z_grid[k] << " ";
        out << "\n";
    }else{
        out << "UNIFORM_GRID\n";
    }
    
    for (int i = 0; i < n_x; ++i) {
        for (int j = 0; j < n_y; ++j) {
            for (int k = 0; k < n_z; ++k) {
                out << field[i * n_y * n_z + j * n_z + k] << " ";
            }
            out << "\n";
        }
        out << "\n";
    }
    out.close();
}

void printGrid(const std::string& filename, int n_x, int n_y, int n_z, double lengthX, double lengthY, double lengthZ,
    bool non_uniform_grid, std::vector<double> x_grid, std::vector<double> y_grid, std::vector<double> z_grid) {

    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Could not open file for writing grid\n";
        return;
    }

    // Write dimensions as header
    out << n_x << " " << n_y << " " << n_z << " " 
    << lengthX << " " << lengthY << " " << lengthZ << "\n";

    if(non_uniform_grid) {
        // Write grid points
        out << "NON_UNIFORM_GRID\n";
            for (int i = 0; i < n_x; ++i) out << x_grid[i] << " ";
            out << "\n";
            for (int j = 0; j < n_y; ++j) out << y_grid[j] << " ";
            out << "\n";
            for (int k = 0; k < n_z; ++k) out << z_grid[k] << " ";
            out << "\n";
        }else{
        out << "UNIFORM_GRID\n";
    }

    out.close();
}

std::vector<double> computeError(const double* reference, const double* computed, int size) {
    std::vector<double> error(size);
    for (int i = 0; i < size; ++i) error[i] = computed[i] - reference[i];
    return error;
}

std::vector<double> computeRelativeError(const double* reference, const double* computed, int size) {
    std::vector<double> relativeError(size);
    for (int i = 0; i < size; ++i) {
        if (std::abs(reference[i]) > 1e-14)
            relativeError[i] = computed[i] - reference[i] / std::abs(reference[i]);
        else if (std::abs(computed[i]) > 1e-8)
            relativeError[i] = computed[i];
        else
            relativeError[i] = 0;
    }
    return relativeError;
}

void applyReLU(double* input, int size) {
    for (int i = 0; i < size; ++i) {
        input[i] = std::max(0.0, input[i]);
    }
}

double linear_interp(double x, const std::vector<double>& xs, const std::vector<double>& ys) {
    if (x <= xs.front()) return ys.front();
    if (x >= xs.back()) return ys.back();
    for (std::size_t i = 0; i < xs.size() - 1; ++i) {
        if (x >= xs[i] && x < xs[i+1]) {
            double t = (x - xs[i]) / (xs[i+1] - xs[i]);
            return ys[i] * (1.0 - t) + ys[i+1] * t;
        }
    }
    return ys.back();
}

double linear_interpolation_binary(double x0, double y0, double x1, double y1, double x) {
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

template <typename T>
class CubicHermiteSpline {
public:
    CubicHermiteSpline(const T * x_ptr, const T * y_ptr, const T * m_ptr, const size_t size):
        x_(x_ptr, x_ptr + size),  y_(y_ptr, y_ptr + size), m_(m_ptr, m_ptr + size), size_(size) {}

    T get_interpolated_value(const T x) const {
        const size_t idx = binary_search_(x);
        if(idx == size_ - 1)
            return y_[size_-1];
        const T t = (x - x_[idx]) / (x_[idx+1] - x_[idx]);
        return interp_func_(t, idx);
    }

private:
    bool is_x_in_boundary(const size_t idx, const T x) const {
        return (x_[idx] <= x) && (x < x_[idx+1]);
    }

    size_t binary_search_(const T x) const {
        assert((x_[0] <= x) && (x <= x_[size_-1]));
        size_t idx_l = 0, idx_r = size_ - 1, idx = size_ / 2;
        while (1) {
            if(idx_r - idx_l == 1)  {
                if(is_x_in_boundary(idx, x))
                    return idx;
                else
                    return (idx + 1);
            }
            if(is_x_in_boundary(idx, x))
                return idx;
            else if(x_[idx+1] <= x) {
                idx_l = idx;
                idx = (idx_r - idx_l) / 2 + idx_l;
            }
            else {
                idx_r = idx;
                idx = (idx_r - idx_l) / 2 + idx_l;
            }
        }
    }

    T interp_func_(const T t, const size_t idx) const {
        return (2 * std::pow(t, 3) - 3 * std::pow(t, 2) + 1) * y_[idx] + (std::pow(t, 3) - 2 * std::pow(t, 2) + t) * (x_[idx+1] - x_[idx]) * m_[idx] +
               (-2 * std::pow(t, 3) + 3 * std::pow(t, 2)) * y_[idx+1] + (std::pow(t, 3) - std::pow(t, 2))*(x_[idx+1] - x_[idx]) * m_[idx+1];
    }

    std::vector<T> x_, y_, m_;
    size_t size_;
};


template <typename T>
class MonotoneCubicInterpolation {
public:
    MonotoneCubicInterpolation(const T * x_ptr, const T * y_ptr, const size_t size) {
        std::vector<T> delta(size, 0);
        std::vector<T> m(size, 0);
        for(unsigned int i=0; i<size-1; ++i) delta[i] = (y_ptr[i+1] - y_ptr[i]) / (x_ptr[i+1] - x_ptr[i]);
        for(unsigned int i=1; i<size-1; ++i) m[i] = (delta[i-1] + delta[i]) / 2;
        m[0] = delta[0];
        m[size-1] = delta[size-2];
        for(unsigned int i=0; i<(size-1); ++i) {
            if(std::abs(delta[i]) < keps)
                m[i] = m[i+1] = 0;
        }
        spliner_ptr_ = std::make_unique<CubicHermiteSpline<T>>(x_ptr, y_ptr, m.data(), size);
    }

    T operator()(const T x) const {
        return spliner_ptr_ -> get_interpolated_value(x);
    }

private:
    static constexpr const T keps = 1e-10;
    std::unique_ptr<CubicHermiteSpline<T>> spliner_ptr_;
};


std::pair<int,int> find_indices(const std::vector<double> &arr, double r) {
    int n = arr.size();
    
    // Handle outside range
    if (r <= arr.front()) return {0, 1};
    if (r >= arr.back())  return {n-2, n-1};

    // Binary search
    int low = 0, high = n - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (arr[mid] <= r && r < arr[mid+1]) {
            return {mid, mid+1};
        } else if (arr[mid] < r) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return {-1, -1}; // should never happen if r is within bounds
}


std::vector<double> complexToAbsArray(const std::vector<std::vector<std::complex<double>>>& complexVec, int component_index) {
    std::vector<double> absArray(complexVec.size());
    for (size_t i = 0; i < complexVec.size(); ++i)
        absArray[i] = std::abs(complexVec[i][component_index]);
    return absArray;
}


#endif // UTILS_H
