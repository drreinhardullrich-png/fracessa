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

#include <Eigen/Dense>

namespace Eigen {
namespace Ext {

template <typename MatrixType>
class BareissGauss {
public:
    using Scalar      = typename MatrixType::Scalar;
    using Index       = typename MatrixType::Index;
    using VectorType  = Eigen::Matrix<Scalar, MatrixType::RowsAtCompileTime, 1>;

    explicit BareissGauss(const MatrixType& A) : m_A(A) {}

    bool solve(const VectorType& b, VectorType& x) const {
        const Index n = m_A.rows();
        if (m_A.cols() != n || b.size() != n) return false;

        // Build augmented matrix [A | b]
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> M(n, n + 1);
        M.leftCols(n) = m_A;
        M.col(n) = b;

        Scalar divPrev = Scalar(1);

        // Bareiss elimination with partial pivoting
        for (Index k = 0; k < n - 1; ++k) {
            // Partial pivoting: find row with maximum absolute value in column k
            Index max_row = k;
            Scalar max_val = M(k, k);
            if (max_val < Scalar(0)) max_val = -max_val;
            
            for (Index i = k + 1; i < n; ++i) {
                Scalar val = M(i, k);
                if (val < Scalar(0)) val = -val;
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }
            
            // Swap rows if necessary
            if (max_row != k) {
                M.row(k).swap(M.row(max_row));
            }
            
            const Scalar pivot = M(k, k);
            if (pivot == Scalar(0)) {
                return false;
            }

            // Bareiss elimination step
            for (Index i = k + 1; i < n; ++i) {
                for (Index j = k + 1; j <= n; ++j) {
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
            for (Index j = i + 1; j < n; ++j) {
                sum -= M(i, j) * x(j);
            }

            const Scalar pivot = M(i, i);
            if (pivot == Scalar(0)) return false;

            x(i) = sum / pivot;
        }

        return true;
    }

private:
    MatrixType m_A;
};

} // namespace Ext
} // namespace Eigen
