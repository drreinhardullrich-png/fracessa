#ifndef RATIONAL_LINALG_BAREISS_GAUSS_HPP
#define RATIONAL_LINALG_BAREISS_GAUSS_HPP

#include <rational_linalg/matrix.hpp>
#include <stdexcept>

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

} // namespace rational_linalg

#endif // RATIONAL_LINALG_BAREISS_GAUSS_HPP

