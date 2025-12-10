#pragma once
#include <Eigen/Dense>
#include <cmath>

namespace Eigen {
namespace Ext {

template <typename MatrixType>
class BareissLUFactor {
public:
    using Scalar     = typename MatrixType::Scalar;
    using Index      = typename MatrixType::Index;
    using VectorType = Eigen::Matrix<Scalar, MatrixType::RowsAtCompileTime, 1>;

    BareissLUFactor(const MatrixType& A) {
        compute(A);
    }

    // ------------------------------------------------------------------------
    // Perform fraction-free LU factorization using Bareiss algorithm
    // ------------------------------------------------------------------------
    void compute(const MatrixType& A) {
        const Index n = A.rows();
        m_n = n;
        m_L = MatrixType::Identity(n, n);
        m_U = A;
        m_P.setIdentity(n, n);
        m_swap_count = 0;

        Scalar divPrev = Scalar(1);

        for (Index k = 0; k < n - 1; ++k) {

            // ----- Partial Pivoting -----
            Index max_row = k;
            Scalar max_val = m_U(k, k);
            if (max_val < Scalar(0)) max_val = -max_val;
            
            for (Index i = k + 1; i < n; ++i) {
                Scalar val = m_U(i, k);
                if (val < Scalar(0)) val = -val;
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }
            
            if (max_row != k) {
                m_U.row(k).swap(m_U.row(max_row));
                m_P.row(k).swap(m_P.row(max_row));
                // L rows up to k-1 also swap to maintain LU decomposition
                if (k > 0)
                    m_L.block(k, 0, 1, k).swap(m_L.block(max_row, 0, 1, k));
                m_swap_count++;
            }
            
            const Scalar pivot = m_U(k, k);
            if (pivot == Scalar(0)) {
                m_is_singular = true;
                return;
            }

            // ----- Bareiss Fraction-Free Update -----
            for (Index i = k + 1; i < n; ++i) {
                m_L(i, k) = m_U(i, k) / pivot;   // store multiplier

                for (Index j = k + 1; j < n; ++j) {
                    m_U(i, j) =
                        (m_U(i, j) * pivot - m_U(i, k) * m_U(k, j)) / divPrev;
                }
                m_U(i, k) = Scalar(0);
            }

            divPrev = pivot;
        }

        m_is_singular = (m_U(n - 1, n - 1) == Scalar(0));
    }

    // ------------------------------------------------------------------------
    // Compute exact determinant
    // ------------------------------------------------------------------------
    Scalar determinant() const {
        if (m_is_singular) return Scalar(0);

        Scalar det = Scalar(1);

        // determinant(P) = sign of permutation = (-1)^(swap_count)
        if (m_swap_count % 2 == 1) det = Scalar(-1);

        for (Index i = 0; i < m_n; ++i)
            det *= m_U(i, i);

        return det;
    }

    // ------------------------------------------------------------------------
    // Compute inverse matrix via LU solve
    // ------------------------------------------------------------------------
    MatrixType inverse() const {
        MatrixType Inv(m_n, m_n);

        if (m_is_singular)
            throw std::runtime_error("Matrix is singular");

        for (Index col = 0; col < m_n; ++col) {
            VectorType e = VectorType::Zero(m_n);
            e(col) = Scalar(1);
            Inv.col(col) = solve(e);
        }

        return Inv;
    }

    // ------------------------------------------------------------------------
    // Solve Ax = b using computed LU: A = P^T * L * U
    // ------------------------------------------------------------------------
    VectorType solve(const VectorType& b) const {
        VectorType bp = m_P * b;
        VectorType y(m_n), x(m_n);

        // Forward substitution: L y = bp
        for (Index i = 0; i < m_n; ++i) {
            Scalar sum = bp(i);
            for (Index j = 0; j < i; ++j)
                sum -= m_L(i, j) * y(j);
            y(i) = sum;
        }

        // Back substitution: U x = y
        for (Index i = m_n - 1; i >= 0; --i) {
            Scalar sum = y(i);
            for (Index j = i + 1; j < m_n; ++j)
                sum -= m_U(i, j) * x(j);
            x(i) = sum / m_U(i, i);
        }

        return x;
    }

    bool isSingular() const { return m_is_singular; }

private:
    Index m_n;
    MatrixType m_L, m_U;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> m_P;
    bool m_is_singular = false;
    int m_swap_count = 0;
};

} // namespace Ext
} // namespace Eigen
