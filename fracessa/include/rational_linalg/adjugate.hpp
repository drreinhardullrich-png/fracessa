#ifndef RATIONAL_LINALG_ADJUGATE_HPP
#define RATIONAL_LINALG_ADJUGATE_HPP

#include <rational_linalg/matrix.hpp>
#include <rational_linalg/lu.hpp>

namespace rational_linalg {

// Compute adjugate using cofactor expansion with LUFactor
template<typename T>
inline Matrix<T> adjugate(const Matrix<T>& A) {
    const size_t n = A.rows();
    
    // Handle edge cases
    if (n == 0) {
        return Matrix<T>(0, 0);
    }
    
    if (n == 1) {
        Matrix<T> result(1, 1);
        result(0, 0) = T(1);
        return result;
    }
    
    // Build cofactor matrix
    Matrix<T> cofactor_matrix(n, n);
    
    // Reuse a single minor matrix to avoid n² allocations
    Matrix<T> minor(n - 1, n - 1);
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            // Extract minor M_ij (remove row i, column j) into reused matrix
            size_t minor_row = 0;
            for (size_t row = 0; row < n; ++row) {
                if (row == i) continue;
                size_t minor_col = 0;
                for (size_t col = 0; col < n; ++col) {
                    if (col == j) continue;
                    minor(minor_row, minor_col) = A(row, col);
                    ++minor_col;
                }
                ++minor_row;
            }
            
            // Compute determinant of minor using LUFactor
            LUFactor<T> lu(minor);
            T det_minor = lu.determinant();
            
            // Compute cofactor: (-1)^(i+j) * det(M_ij)
            // Use bitwise AND for sign: (i+j) & 1 is faster than % 2
            T sign = ((i + j) & 1) == 0 ? T(1) : T(-1);
            cofactor_matrix(i, j) = sign * det_minor;
        }
    }
    
    // Adjugate is the transpose of the cofactor matrix
    return cofactor_matrix.transpose();
}

} // namespace rational_linalg

#endif // RATIONAL_LINALG_ADJUGATE_HPP

