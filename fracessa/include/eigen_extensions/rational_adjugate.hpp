#pragma once

#include <Eigen/Dense>
#include <eigen_extensions/bareiss_lu.hpp>
#include <fracessa/matrix_old.hpp>

namespace Eigen {
namespace Ext {

// Compute adjugate using cofactor expansion with BareissLUFactor
inline RationalMatrix adjugate(const RationalMatrix& A) {
    const Eigen::Index n = A.rows();
    
    // Handle edge cases
    if (n == 0) {
        return RationalMatrix(0, 0);
    }
    
    if (n == 1) {
        RationalMatrix result(1, 1);
        result(0, 0) = rational(1);
        return result;
    }
    
    // Build cofactor matrix
    RationalMatrix cofactor_matrix(n, n);
    
    // Reuse a single minor matrix to avoid n² allocations
    RationalMatrix minor(n - 1, n - 1);
    
    for (Eigen::Index i = 0; i < n; ++i) {
        for (Eigen::Index j = 0; j < n; ++j) {
            // Extract minor M_ij (remove row i, column j) into reused matrix
            Eigen::Index minor_row = 0;
            for (Eigen::Index row = 0; row < n; ++row) {
                if (row == i) continue;
                Eigen::Index minor_col = 0;
                for (Eigen::Index col = 0; col < n; ++col) {
                    if (col == j) continue;
                    minor(minor_row, minor_col) = A(row, col);
                    ++minor_col;
                }
                ++minor_row;
            }
            
            // Compute determinant of minor using BareissLUFactor
            BareissLUFactor<RationalMatrix> lu(minor);
            rational det_minor = lu.determinant();
            
            // Compute cofactor: (-1)^(i+j) * det(M_ij)
            // Use bitwise AND for sign: (i+j) & 1 is faster than % 2
            rational sign = ((i + j) & 1) == 0 ? rational(1) : rational(-1);
            cofactor_matrix(i, j) = sign * det_minor;
        }
    }
    
    // Adjugate is the transpose of the cofactor matrix
    return cofactor_matrix.transpose();
}

} // namespace Ext
} // namespace Eigen

