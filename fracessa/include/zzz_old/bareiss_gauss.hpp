#ifndef RATIONAL_LINALG_BAREISS_GAUSS_HPP
#define RATIONAL_LINALG_BAREISS_GAUSS_HPP

#include <rational_linalg/matrix.hpp>
#include <cmath>
#include <limits>

namespace rational_linalg {

/*
 * BareissGauss (fraction-free Gaussian elimination) solver for Matrix<T>
 *
 * Speed:
 *   - Complexity: O(n^3), same as standard Gaussian elimination.
 *   - Much slower in practice with rational or big integer types
 *     because intermediate numbers grow rapidly.
 *   - For floating-point, Bareiss is slower than standard methods
 *     because it avoids floating-point division tricks and uses exact arithmetic.
 *
 * Stability / Accuracy:
 *   - Perfect for rational or symbolic types; no rounding errors.
 *   - Detects singular matrices exactly (pivot = 0).
 *   - Great for ill-conditioned matrices in exact arithmetic because
 *     no division by small numbers until necessary.
 *   - Cannot suffer from numerical instability in exact arithmetic.
 *
 * Use case:
 *   - Exact solutions, guaranteed singularity detection,
 *     symbolic or rational arithmetic.
 */

template<typename T>
class BareissGauss {
public:
    using Scalar = T;
    using VectorType = Matrix<T>;

    explicit BareissGauss(const Matrix<T>& A) : m_A(A) {}

    bool solve(const Matrix<T>& b, Matrix<T>& x) const {
        // b must be a column vector
        if (b.cols() != 1) return false;
        
        const size_t n = m_A.rows();
        if (m_A.cols() != n || b.rows() != n) return false;

        // Build augmented matrix [A | b]
        Matrix<T> M(n, n + 1);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                M(i, j) = m_A(i, j);
            }
            M(i, n) = b(i, 0);
        }

        T divPrev = T(1);

        // Bareiss elimination with partial pivoting
        for (size_t k = 0; k < n - 1; ++k) {
            // Partial pivoting: find row with maximum absolute value in column k
            size_t max_row = k;
            T max_val = M(k, k);
            if (max_val < T(0)) max_val = -max_val;
            
            for (size_t i = k + 1; i < n; ++i) {
                T val = M(i, k);
                if (val < T(0)) val = -val;
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }
            
            // Swap rows if necessary
            if (max_row != k) {
                M.swap_rows(k, max_row);
            }
            
            const T pivot = M(k, k);
            if (pivot == T(0)) {
                return false;
            }

            // Bareiss elimination step
            for (size_t i = k + 1; i < n; ++i) {
                for (size_t j = k + 1; j <= n; ++j) {
                    M(i, j) = (M(i, j) * pivot - M(i, k) * M(k, j)) / divPrev;
                }
                M(i, k) = T(0);
            }

            divPrev = pivot;
        }

        if (M(n - 1, n - 1) == T(0)) return false;

        // Back substitution
        x = Matrix<T>(n, 1);
        for (size_t i = n; i-- > 0; ) {
            T sum = M(i, n);
            for (size_t j = i + 1; j < n; ++j) {
                sum -= M(i, j) * x(j, 0);
            }

            const T pivot = M(i, i);
            if (pivot == T(0)) return false;

            x(i, 0) = sum / pivot;
        }

        return true;
    }

private:
    Matrix<T> m_A;
};

/*
 * GaussDouble - Optimized standard Gaussian elimination for Matrix<double>
 *
 * Speed:
 *   - Complexity: O(n^3), same as BareissGauss
 *   - Much faster than BareissGauss for double because:
 *     - Uses direct division (no fraction-free arithmetic overhead)
 *     - Optimized for floating-point operations
 *     - Better cache locality
 *
 * Stability / Accuracy:
 *   - Uses partial pivoting for numerical stability
 *   - Tolerance-based singularity detection
 *   - Standard floating-point precision limitations apply
 *
 * Use case:
 *   - Fast solving of double-precision linear systems
 *   - When exact arithmetic is not required
 */
class GaussDouble {
public:
    using Scalar = double;
    using VectorType = Matrix<double>;

    explicit GaussDouble(const Matrix<double>& A) : m_A(A) {}

    bool solve(const Matrix<double>& b, Matrix<double>& x) const {
        // b must be a column vector
        if (b.cols() != 1) return false;
        
        const size_t n = m_A.rows();
        if (m_A.cols() != n || b.rows() != n) return false;

        // Build augmented matrix [A | b]
        Matrix<double> M(n, n + 1);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                M(i, j) = m_A(i, j);
            }
            M(i, n) = b(i, 0);
        }

        const double epsilon = std::numeric_limits<double>::epsilon() * n;

        // Standard Gaussian elimination with partial pivoting
        for (size_t k = 0; k < n - 1; ++k) {
            // Partial pivoting: find row with maximum absolute value in column k
            size_t max_row = k;
            double max_val = std::abs(M(k, k));
            
            for (size_t i = k + 1; i < n; ++i) {
                double val = std::abs(M(i, k));
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }
            
            // Check for singularity
            if (max_val < epsilon) {
                return false;
            }
            
            // Swap rows if necessary
            if (max_row != k) {
                M.swap_rows(k, max_row);
            }
            
            const double pivot = M(k, k);

            // Standard elimination step (direct division, no fraction-free)
            for (size_t i = k + 1; i < n; ++i) {
                const double factor = M(i, k) / pivot;
                for (size_t j = k + 1; j <= n; ++j) {
                    M(i, j) -= factor * M(k, j);
                }
                M(i, k) = 0.0; // Explicitly zero (for clarity, though not strictly necessary)
            }
        }

        // Check last pivot
        if (std::abs(M(n - 1, n - 1)) < epsilon) {
            return false;
        }

        // Back substitution
        x = Matrix<double>(n, 1);
        for (size_t i = n; i-- > 0; ) {
            double sum = M(i, n);
            for (size_t j = i + 1; j < n; ++j) {
                sum -= M(i, j) * x(j, 0);
            }

            const double pivot = M(i, i);
            if (std::abs(pivot) < epsilon) {
                return false;
            }

            x(i, 0) = sum / pivot;
        }

        return true;
    }

private:
    Matrix<double> m_A;
};

// Helper struct for automatic solver selection
// Automatically selects GaussDouble for double, BareissGauss<T> for other types
template<typename T>
struct SolverSelector {
    using type = BareissGauss<T>;
};

template<>
struct SolverSelector<double> {
    using type = GaussDouble;
};

// Type alias for automatic solver selection
template<typename T>
using LinearSolver = typename SolverSelector<T>::type;

} // namespace rational_linalg

#endif // RATIONAL_LINALG_BAREISS_GAUSS_HPP

