#pragma once

#include <Eigen/Dense>
#include <unsupported/Eigen/FFT>
#include <complex>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace Eigen {
namespace Ext {

/**
 * Solve a circulant linear system A*x = b using FFT.
 * 
 * A circulant matrix is defined by its first row vector 'a', where each row
 * is a cyclic shift of the previous row. This function solves the system
 * efficiently using the Fast Fourier Transform (O(n log n) instead of O(n^3)).
 * 
 * Algorithm:
 * 1. Compute eigenvalues of circulant matrix via FFT of first row
 * 2. Compute FFT of right-hand side vector b
 * 3. Divide in Fourier space (element-wise)
 * 4. Inverse FFT to get solution
 * 
 * @param a First row of the circulant matrix (defines the entire matrix)
 * @param b Right-hand side vector
 * @return Solution vector x such that A*x = b
 * @throws std::runtime_error if matrix is singular (zero eigenvalue detected)
 */
inline Eigen::VectorXd solveCirculantEigenFFT(const Eigen::VectorXd &a, const Eigen::VectorXd &b) {
    Eigen::Index n = a.size();
    Eigen::FFT<double> fft;

    double tol = 1e-12;
    
    // Convert Eigen::VectorXd to std::vector for FFT
    std::vector<double> a_vec(a.data(), a.data() + n);
    std::vector<double> b_vec(b.data(), b.data() + n);
    
    // 1. Compute eigenvalues of circulant matrix A via FFT of first row
    std::vector<std::complex<double>> psi;
    fft.fwd(psi, a_vec);

    // 2. FFT of right-hand side b
    std::vector<std::complex<double>> b_hat;
    fft.fwd(b_hat, b_vec);

    // 3. Element-wise division in Fourier space
    std::vector<std::complex<double>> y(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        if (std::abs(psi[i]) < tol) {
            throw std::runtime_error("Circulant matrix is singular or nearly singular (zero eigenvalue).");
        }
        y[i] = b_hat[i] / psi[i];
    }

    // 4. Inverse FFT to get solution x
    std::vector<std::complex<double>> x_complex;
    fft.inv(x_complex, y);

    // Extract real part (Eigen's inv() already normalizes, so no division by n needed)
    Eigen::VectorXd x(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        x[i] = x_complex[i].real();  // discard small imaginary parts from numerical errors
    }

    return x;
}

} // namespace Ext
} // namespace Eigen