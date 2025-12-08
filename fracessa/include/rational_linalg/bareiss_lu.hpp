#ifndef RATIONAL_LINALG_BAREISS_LU_HPP
#define RATIONAL_LINALG_BAREISS_LU_HPP

#include <rational_linalg/matrix.hpp>
#include <stdexcept>

namespace rational_linalg {

/*
 * BareissLUFactor (fraction-free LU factorization) for Matrix<T>
 *
 * Speed:
 *   - Complexity: O(n^3), same as standard LU.
 *   - Much slower in practice with rational or big integer types
 *     because intermediate numbers grow rapidly.
 *   - For floating-point, Bareiss is slower than PartialPivLU
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
class BareissLUFactor {
public:
    using Scalar = T;
    using VectorType = Matrix<T>;

    BareissLUFactor(const Matrix<T>& A) {
        compute(A);
    }

    // ------------------------------------------------------------------------
    // Perform fraction-free LU factorization using Bareiss algorithm
    // ------------------------------------------------------------------------
    void compute(const Matrix<T>& A) {
        const size_t n = A.rows();
        m_n = n;
        m_L = Matrix<T>::identity(n);
        m_U = A;
        m_P = Matrix<T>::identity(n);
        m_swap_count = 0;

        T divPrev = T(1);

        for (size_t k = 0; k < n - 1; ++k) {

            // ----- Partial Pivoting -----
            size_t max_row = k;
            T max_val = m_U(k, k);
            if (max_val < T(0)) max_val = -max_val;
            
            for (size_t i = k + 1; i < n; ++i) {
                T val = m_U(i, k);
                if (val < T(0)) val = -val;
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }
            
            if (max_row != k) {
                m_U.swap_rows(k, max_row);
                m_P.swap_rows(k, max_row);
                // L rows up to k-1 also swap to maintain LU decomposition
                if (k > 0) {
                    for (size_t j = 0; j < k; ++j) {
                        T tmp = m_L(k, j);
                        m_L(k, j) = m_L(max_row, j);
                        m_L(max_row, j) = tmp;
                    }
                }
                m_swap_count++;
            }
            
            const T pivot = m_U(k, k);
            if (pivot == T(0)) {
                m_is_singular = true;
                return;
            }

            // ----- Bareiss Fraction-Free Update -----
            for (size_t i = k + 1; i < n; ++i) {
                m_L(i, k) = m_U(i, k) / pivot;   // store multiplier

                for (size_t j = k + 1; j < n; ++j) {
                    m_U(i, j) =
                        (m_U(i, j) * pivot - m_U(i, k) * m_U(k, j)) / divPrev;
                }
                m_U(i, k) = T(0);
            }

            divPrev = pivot;
        }

        m_is_singular = (m_U(n - 1, n - 1) == T(0));
    }

    // ------------------------------------------------------------------------
    // Compute exact determinant
    // ------------------------------------------------------------------------
    T determinant() const {
        if (m_is_singular) return T(0);

        T det = T(1);

        // determinant(P) = sign of permutation = (-1)^(swap_count)
        if (m_swap_count % 2 == 1) det = T(-1);

        for (size_t i = 0; i < m_n; ++i)
            det *= m_U(i, i);

        return det;
    }

    // ------------------------------------------------------------------------
    // Compute inverse matrix via LU solve
    // ------------------------------------------------------------------------
    Matrix<T> inverse() const {
        Matrix<T> Inv(m_n, m_n);

        if (m_is_singular)
            throw std::runtime_error("Matrix is singular");

        for (size_t col = 0; col < m_n; ++col) {
            Matrix<T> e = Matrix<T>::Zero(m_n);
            e(col, 0) = T(1);
            Matrix<T> col_result = solve(e);
            // Copy column result into inverse
            for (size_t i = 0; i < m_n; ++i) {
                Inv(i, col) = col_result(i, 0);
            }
        }

        return Inv;
    }

    // ------------------------------------------------------------------------
    // Solve Ax = b using computed LU: A = P^T * L * U
    // ------------------------------------------------------------------------
    Matrix<T> solve(const Matrix<T>& b) const {
        // b must be a column vector
        if (b.cols() != 1 || b.rows() != m_n) {
            throw std::runtime_error("solve: b must be a column vector of size n");
        }
        
        // Apply permutation: bp = P * b
        Matrix<T> bp(m_n, 1);
        for (size_t i = 0; i < m_n; ++i) {
            T sum = T(0);
            for (size_t j = 0; j < m_n; ++j) {
                sum += m_P(i, j) * b(j, 0);
            }
            bp(i, 0) = sum;
        }
        
        Matrix<T> y(m_n, 1);
        Matrix<T> x(m_n, 1);

        // Forward substitution: L y = bp
        for (size_t i = 0; i < m_n; ++i) {
            T sum = bp(i, 0);
            for (size_t j = 0; j < i; ++j)
                sum -= m_L(i, j) * y(j, 0);
            y(i, 0) = sum;
        }

        // Back substitution: U x = y
        for (size_t i = m_n; i-- > 0; ) {
            T sum = y(i, 0);
            for (size_t j = i + 1; j < m_n; ++j)
                sum -= m_U(i, j) * x(j, 0);
            x(i, 0) = sum / m_U(i, i);
        }

        return x;
    }

    bool isSingular() const { return m_is_singular; }

private:
    size_t m_n;
    Matrix<T> m_L, m_U, m_P;
    bool m_is_singular = false;
    int m_swap_count = 0;
};

} // namespace rational_linalg

#endif // RATIONAL_LINALG_BAREISS_LU_HPP

