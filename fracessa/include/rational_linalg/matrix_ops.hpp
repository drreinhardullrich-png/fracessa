#ifndef RATIONAL_LINALG_MATRIX_OPS_HPP
#define RATIONAL_LINALG_MATRIX_OPS_HPP

#include <rational_linalg/matrix.hpp>
#include <fracessa/bitset64.hpp>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Type aliases for Eigen double matrices (still used for double operations)
using DoubleMatrix = Eigen::MatrixXd;
using DoubleVector = Eigen::VectorXd;

namespace rational_linalg {

// Convert Matrix<T> to DoubleMatrix
template<typename T>
inline DoubleMatrix to_double(const Matrix<T>& A)
{
    DoubleMatrix result = DoubleMatrix::Zero(static_cast<Eigen::Index>(A.rows()), static_cast<Eigen::Index>(A.cols()));
    for (size_t i = 0; i < A.rows(); ++i) {
        for (size_t j = 0; j < A.cols(); ++j) {
            result(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = static_cast<double>(A(i, j));
        }
    }
    return result;
}

// String representation
template<typename T>
inline std::string to_string(const Matrix<T>& A)
{
    std::ostringstream oss;
    for (size_t i = 0; i < A.rows(); ++i) {
        oss << "\t\t\t";  // Add tabs before each line
        for (size_t j = 0; j < A.cols(); ++j) {
            oss << A(i, j) << ",";
        }
        if (i < A.rows() - 1) {
            oss << std::endl;
        }
    }
    return oss.str();
}

// Extract principal submatrix from matrix using bitset64 mask
// Optimized: only iterates over SET bits using for_each_set_bit_no_exit for better performance
// NOTE: submatrix must already be correctly sized (support_size x support_size) before calling!
template<typename T>
inline void principal_submatrix(const Matrix<T>& A, size_t /*dimension*/, const bitset64& support, size_t /*support_size*/, Matrix<T>& submatrix)
{
    // Only iterate over SET bits for efficiency
    size_t row = 0;
    support.for_each_set_bit_no_exit([&](unsigned i) {
        size_t col = 0;
        support.for_each_set_bit_no_exit([&](unsigned j) {
            submatrix(row, col) = A(static_cast<size_t>(i), static_cast<size_t>(j));
            ++col;
        });
        ++row;
    });
}

// Get bordered matrix A from support for Matrix<T> (split version for optimization)
template<typename T>
inline void get_kkt_bordering(const Matrix<T>& game_matrix, const bitset64& support, size_t support_size, Matrix<T>& A)
{
    size_t n = support_size + 1;
    // Only resize if needed, and matrix is square, so check for rows is sufficient
    if (A.rows() != n) {
        A = Matrix<T>(n, n);
    } else {
        // Set all elements to zero
        for (size_t i = 0; i < n * n; ++i) {
            A.data()[i] = T(0);
        }
    }
    
    size_t row = 0;
    size_t rows = game_matrix.rows();
    size_t cols = game_matrix.cols();

    // Fill rows 0 to support_size-1: submatrix from game_matrix, then -1 in last column
    for (size_t i = 0; i < rows; ++i) {
        if (support.test(i)) {
            size_t column = 0;
            for (size_t j = 0; j < cols; ++j) {
                if (support.test(j)) {
                    A(row, column) = game_matrix(i, j);
                    column++;
                }
            }
            A(row, support_size) = T(-1);
            row++;
        }
    }
    
    // Fill last row: all 1s, then 0 in last column
    for (size_t i = 0; i < support_size; ++i) {
        A(support_size, i) = T(1);
    }
    A(support_size, support_size) = T(0);
}

// Get RHS vector b for Matrix<T> (split version for optimization)
// Returns Matrix<T>(n, 1) as column vector
template<typename T>
inline void get_kkt_rhs(size_t support_size, Matrix<T>& b)
{
    size_t n = support_size + 1;
    // Only resize if needed
    if (b.rows() != n || b.cols() != 1) {
        b = Matrix<T>(n, 1);
    } else {
        // Set all elements to zero
        for (size_t i = 0; i < n; ++i) {
            b(i, 0) = T(0);
        }
    }
    // Last element is 1
    b(support_size, 0) = T(1);
}

// Check if matrix is positive definite (LDLT decomposition for rational)
template<typename T>
inline bool is_positive_definite_rational(const Matrix<T>& A)
{
    // LDLT decomposition for rational numbers
    size_t n = A.rows();
    Matrix<T> D(n, n);
    Matrix<T> L(n, n);
    
    // Initialize L as identity
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            L(i, j) = (i == j) ? T(1) : T(0);
        }
    }

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            T aSum = T(0);
            for (size_t k = 0; k < j; ++k) {
                aSum += L(i, k) * L(j, k) * D(k, k);
            }
            L(i, j) = (T(1) / D(j, j)) * (A(i, j) - aSum);
        }
        T bSum = T(0);
        for (size_t k = 0; k < i; ++k) {
            bSum += L(i, k) * L(i, k) * D(k, k);
        }
        D(i, i) = A(i, i) - bSum;
        if (D(i, i) <= T(0)) {
            return false;
        }
    }
    return true;
}

// Check if all entries of a matrix are greater than zero
template<typename T>
inline bool all_entries_greater_zero(const Matrix<T>& A) {
    for (size_t i = 0; i < A.rows(); ++i) {
        for (size_t j = 0; j < A.cols(); ++j) {
            if (A(i, j) <= T(0)) {
                return false;
            }
        }
    }
    return true;
}

} // namespace rational_linalg

// Double matrix operations (still use Eigen types, keep in global namespace for compatibility)
namespace matrix_ops {

// Get bordered matrix A from support (split version for optimization) - Double version
inline void get_kkt_bordering(const DoubleMatrix& game_matrix, const bitset64& support, size_t support_size, DoubleMatrix& A)
{
    size_t n = support_size + 1;
    // Only resize if needed, and matrix is square, so check for rows is sufficient
    if (A.rows() != static_cast<Eigen::Index>(n)) {
        A = DoubleMatrix::Zero(n, n);
    } else {
        A.setZero();
    }
    
    size_t row = 0;
    Eigen::Index rows = static_cast<Eigen::Index>(game_matrix.rows());
    Eigen::Index cols = static_cast<Eigen::Index>(game_matrix.cols());

    // Fill rows 0 to support_size-1: submatrix from game_matrix, then -1 in last column
    for (Eigen::Index i = 0; i < rows; ++i) {
        if (support.test(static_cast<size_t>(i))) {
            size_t column = 0;
            for (Eigen::Index j = 0; j < cols; ++j) {
                if (support.test(static_cast<size_t>(j))) {
                    A(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) = game_matrix(i, j);
                    column++;
                }
            }
            A(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(support_size)) = -1.0;
            row++;
        }
    }
    
    // Fill last row: all 1s, then 0 in last column
    for (size_t i = 0; i < support_size; ++i) {
        A(static_cast<Eigen::Index>(support_size), static_cast<Eigen::Index>(i)) = 1.0;
    }
    A(static_cast<Eigen::Index>(support_size), static_cast<Eigen::Index>(support_size)) = 0.0;
}

// Get RHS vector b (split version for optimization) - Double version
inline void get_kkt_rhs(size_t support_size, DoubleVector& b)
{
    size_t n = support_size + 1;
    // Only resize if needed
    if (b.size() != static_cast<Eigen::Index>(n)) {
        b = DoubleVector::Zero(n);
    }
    // Set all elements to 0 first
    b.setZero();
    // Last element is 1
    b(static_cast<Eigen::Index>(support_size)) = 1.0;
}

} // namespace matrix_ops

#endif // RATIONAL_LINALG_MATRIX_OPS_HPP

