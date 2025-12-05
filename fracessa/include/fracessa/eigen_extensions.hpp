#pragma once

/*
 * BareissLU (fraction-free LU) solver
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

 
// -----------------------------------------------------------------------------
// Eigen Extensions: BareissLU
// Fraction-free Gaussian elimination (Bareiss algorithm)
// API modeled after Eigen's PartialPivLU / FullPivLU
// Correct version: compute(A) only stores matrix, solve(b) creates augmented matrix internally
// -----------------------------------------------------------------------------

#include <Eigen/Dense>
#include <boost/multiprecision/cpp_int.hpp>
#include <cmath>

namespace Eigen {
namespace Ext {

template <typename MatrixType>
class BareissLU {
public:
    using Scalar      = typename MatrixType::Scalar;
    using Index       = typename MatrixType::Index;
    using VectorType  = Eigen::Matrix<Scalar, MatrixType::RowsAtCompileTime, 1>;

    BareissLU() : m_isInitialized(false), m_isSingular(false) {}

    explicit BareissLU(const MatrixType& A) {
        compute(A);
    }

    BareissLU& compute(const MatrixType& A) {
        m_isInitialized = true;
        m_isSingular = false;
        m_A = A;
        return *this;
    }

    bool isInitialized() const { return m_isInitialized; }
    bool isInvertible() const { return !m_isSingular; }
    bool isSingular() const { return m_isSingular; }

    bool solve(const VectorType& b, VectorType& x) const {
        if (!m_isInitialized) return false;

        const Index n = m_A.rows();
        if (m_A.cols() != n || b.size() != n) return false;

        // Build augmented matrix [A | b]
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> M(n, n + 1);
        M.block(0, 0, n, n) = m_A;
        M.col(n) = b;

        Scalar divPrev = Scalar(1);

        // Bareiss elimination with partial pivoting
        for (Index k = 0; k < n - 1; ++k) {
            // Partial pivoting: find row with maximum absolute value in column k
            Index max_row = k;
            Scalar max_val = (M(k, k) < Scalar(0)) ? -M(k, k) : M(k, k);  // abs for any Scalar type
            for (Index i = k + 1; i < n; ++i) {
                Scalar val = M(i, k);
                Scalar abs_val = (val < Scalar(0)) ? -val : val;  // abs for any Scalar type
                if (abs_val > max_val) {
                    max_val = abs_val;
                    max_row = i;
                }
            }
            
            // Swap rows if necessary
            if (max_row != k) {
                for (Index j = k; j <= n; ++j) {
                    Scalar tmp = M(k, j);
                    M(k, j) = M(max_row, j);
                    M(max_row, j) = tmp;
                }
            }
            
            Scalar pivot = M(k, k);
            if (pivot == Scalar(0)) {
                return false;
            }

            for (Index i = k + 1; i < n; ++i) {
                for (Index j = k + 1; j <= n; ++j) { // include RHS
                    M(i, j) = (M(i, j) * pivot - M(i, k) * M(k, j)) / divPrev;
                }
                M(i, k) = Scalar(0);
            }

            divPrev = pivot;
        }

        if (M(n - 1, n - 1) == Scalar(0)) return false;

        // Back substitution
        x.resize(n);
        for (Index i = n - 1; i >= 0; --i) {
            Scalar sum = M(i, n);
            for (Index j = i + 1; j < n; ++j) sum -= M(i, j) * x(j);

            Scalar pivot = M(i, i);
            if (pivot == Scalar(0)) return false;

            x(i) = sum / pivot;
        }

        return true;
    }

    VectorType solve(const VectorType& b) const {
        VectorType x;
        solve(b, x);
        return x;
    }

private:
    MatrixType m_A;
    bool m_isInitialized;
    bool m_isSingular;
};

template <typename MatrixType>
inline BareissLU<MatrixType> bareiss(const MatrixType& A) {
    return BareissLU<MatrixType>(A);
}

} // namespace Ext
} // namespace Eigen
