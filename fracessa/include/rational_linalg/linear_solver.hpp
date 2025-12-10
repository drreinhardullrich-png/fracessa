#ifndef RATIONAL_LINALG_LINEAR_SOLVER_HPP
#define RATIONAL_LINALG_LINEAR_SOLVER_HPP

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

    explicit BareissGauss(const Matrix<T>& Ab) : M(Ab) {}

    bool solve(Matrix<T>& x) {
        const size_t n = M.rows();

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
                M(i, k) = T(0); // Zero out eliminated element (not updated in inner loop)
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

            if (pivot == T(0)) 
                return false;            

            T temp_x = sum / pivot;
            if (temp_x <= T(0))
                return false; //early exit if sum is not >zero! then it cannot be a solution (since x on the simplex is always positive!)
            else
                x(i, 0) = temp_x;
        }

        return true;
    }

private:
    Matrix<T> M;
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

    explicit GaussDouble(const Matrix<double>& Ab) : M(Ab) {}

    bool solve(Matrix<double>& x) {
        const size_t n = M.rows();

        const double errorbound = 1e-3 * n; // huge margin, false positives eliminated by rational check

        const double epsilon = std::numeric_limits<double>::epsilon(); // * n; keep epsilon super small to avoid false negative. false positive will be detected by rational afterwards!!!

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
                M(i, k) = 0.0; // Zero out eliminated element (not updated in inner loop)
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

            double temp_x = sum / pivot;
            // Check for not >zero in solution component with huge margin
            // False positives will be eliminated by rational checks later
            if (temp_x < -errorbound) 
                return false;
            else
                x(i, 0) = temp_x;
        }

        return true;
    }

private:
    Matrix<double> M;
};

/*
 * GaussRational<T> - Generic standard Gaussian elimination for Matrix<T>
 *
 * Speed:
 *   - Complexity: O(n^3), same as BareissGauss
 *   - Faster than BareissGauss for types that support efficient division:
 *     - Uses direct division (no fraction-free arithmetic overhead)
 *     - Better for rational types when intermediate values don't grow too large
 *     - More efficient for floating-point types
 *
 * Stability / Accuracy:
 *   - Uses partial pivoting for numerical stability
 *   - Exact zero detection for singularity (no tolerance needed for exact types)
 *   - For floating-point types, may suffer from numerical instability
 *     compared to BareissGauss for ill-conditioned matrices
 *
 * Use case:
 *   - Fast solving of linear systems when exact arithmetic is available
 *   - Alternative to BareissGauss when division is efficient
 *   - For rational types when matrix is well-conditioned
 */
template<typename T>
class GaussRational {
public:
    using Scalar = T;
    using VectorType = Matrix<T>;

    explicit GaussRational(const Matrix<T>& Ab) : M(Ab) {}

    bool solve(Matrix<T>& x) {
        const size_t n = M.rows();

        // Standard Gaussian elimination with partial pivoting
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
            
            // Check for singularity (exact zero check for generic types)
            if (max_val == T(0)) {
                return false;
            }
            
            // Swap rows if necessary
            if (max_row != k) {
                M.swap_rows(k, max_row);
            }
            
            const T pivot = M(k, k);

            // Standard elimination step (direct division, no fraction-free)
            for (size_t i = k + 1; i < n; ++i) {
                const T factor = M(i, k) / pivot;
                for (size_t j = k + 1; j <= n; ++j) {
                    M(i, j) -= factor * M(k, j);
                }
                M(i, k) = T(0); // Zero out eliminated element (not updated in inner loop)
            }
        }

        // Check last pivot
        if (M(n - 1, n - 1) == T(0)) {
            return false;
        }

        // Back substitution
        x = Matrix<T>(n, 1);
        for (size_t i = n; i-- > 0; ) {
            T sum = M(i, n);
            for (size_t j = i + 1; j < n; ++j) {
                sum -= M(i, j) * x(j, 0);
            }

            const T pivot = M(i, i);

            if (pivot == T(0)) 
                return false;

            T temp_x = sum / pivot;
            if (temp_x <= T(0))
                return false; //early exit if sum is not >zero! then it cannot be a solution (since x on the simplex is always positive!)
            else
                x(i, 0) = temp_x;               
        }

        return true;
    }

private:
    Matrix<T> M;
};

// Helper struct for automatic solver selection
// Automatically selects GaussDouble for double, GaussRational<T> for other types
template<typename T>
struct SolverSelector {
    //using type = GaussRational<T>;
    using type = BareissGauss<T>; //leave this here!!! for comparison!!!
};

template<>
struct SolverSelector<double> {
    using type = GaussDouble;
};

// Type alias for automatic solver selection
template<typename T>
using LinearSolver = typename SolverSelector<T>::type;

} // namespace rational_linalg

#endif // RATIONAL_LINALG_LINEAR_SOLVER_HPP

