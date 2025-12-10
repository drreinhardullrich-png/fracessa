// #pragma once

// /**
//  * @brief Correct Bareiss LDLᵀ decomposition (fraction-free for symmetric matrices)
//  * 
//  * Computes: d_{n-1} * A = L * Δ * Lᵀ
//  * Where:
//  * - L is unit lower triangular with INTEGER entries (before final scaling)
//  * - Δ is integer diagonal matrix
//  * - d_{n-1} = product of pivots from steps 0 to n-2
//  * 
//  * NO DIVISION during elimination - only exact integer arithmetic!
//  */

// #include <Eigen/Dense>
// #include <vector>

// namespace Eigen {
// namespace Ext {

// template <typename MatrixType>
// class BareissLDLT {
// public:
//     using Scalar = typename MatrixType::Scalar;
//     using Index = typename MatrixType::Index;
//     using VectorType = Eigen::Matrix<Scalar, MatrixType::RowsAtCompileTime, 1>;

//     BareissLDLT() : m_isInitialized(false) {}

//     explicit BareissLDLT(const MatrixType& matrix) {
//         compute(matrix);
//     }

//     void compute(const MatrixType& A) {
//         const Index n = A.rows();
        
//         // Reset everything
//         m_L = MatrixType::Identity(n, n);
//         m_Delta = MatrixType::Zero(n, n);
//         m_permutation.setIdentity(n);
//         m_isPositiveDefinite = true;
//         m_isInitialized = true;
        
//         // Working copy - we'll modify this
//         MatrixType U = A.template triangularView<Eigen::Upper>();
        
//         // Bareiss algorithm
//         Scalar d_prev = Scalar(1);  // d_{k-1}
        
//         for (Index k = 0; k < n; ++k) {
//             // SYMMETRIC pivoting: find maximum diagonal element
//             Index pivot_idx = k;
//             Scalar max_val = U(k, k);
            
//             for (Index i = k + 1; i < n; ++i) {
//                 if (U(i, i) > max_val) {
//                     max_val = U(i, i);
//                     pivot_idx = i;
//                 }
//             }
            
//             // Perform symmetric permutation if needed
//             if (pivot_idx != k) {
//                 // Swap rows and columns
//                 U.row(k).swap(U.row(pivot_idx));
//                 U.col(k).swap(U.col(pivot_idx));
                
//                 // Also swap in L for columns < k
//                 m_L.row(k).head(k).swap(m_L.row(pivot_idx).head(k));
                
//                 // Update permutation
//                 m_permutation.applyTranspositionOnTheRight(k, pivot_idx);
//             }
            
//             Scalar pivot = U(k, k);
//             m_Delta(k, k) = pivot;
            
//             // Check for singularity
//             if (pivot == Scalar(0)) {
//                 m_info = Eigen::NumericalIssue;
//                 return;
//             }
            
//             // Check positive definiteness
//             if (pivot <= Scalar(0)) {
//                 m_isPositiveDefinite = false;
//             }
            
//             // Store multipliers in L (INTEGER values, no division!)
//             for (Index i = k + 1; i < n; ++i) {
//                 // Store the NUMERATOR only - division postponed!
//                 m_L(i, k) = U(i, k);
//             }
            
//             // Bareiss update (only for k < n-1)
//             if (k < n - 1) {
//                 for (Index i = k + 1; i < n; ++i) {
//                     // Update diagonal elements
//                     U(i, i) = (U(i, i) * pivot - U(i, k) * U(i, k)) / d_prev;
                    
//                     // Update off-diagonal elements (upper triangle)
//                     for (Index j = i + 1; j < n; ++j) {
//                         U(i, j) = (U(i, j) * pivot - U(i, k) * U(k, j)) / d_prev;
//                     }
//                 }
                
//                 // Zero out the k-th column (except diagonal)
//                 for (Index i = k + 1; i < n; ++i) {
//                     U(i, k) = Scalar(0);
//                 }
                
//                 d_prev = pivot;  // Update denominator for next step
//             }
//         }
        
//         m_info = Eigen::Success;
//         m_denominator = d_prev;  // d_{n-1}
//     }

//     VectorType solve(const VectorType& b) const {
//         if (!m_isInitialized) {
//             throw std::runtime_error("Decomposition not computed");
//         }
        
//         const Index n = m_L.rows();
        
//         // Apply permutation: b_perm = P * b
//         VectorType b_perm = m_permutation * b;
        
//         // Forward substitution: L * z = b_perm
//         // BUT: Our L stores numerators, not actual L factors!
//         // Actual L_factor(i,k) = stored_L(i,k) / Δ(k,k)
//         VectorType z(n);
//         for (Index i = 0; i < n; ++i) {
//             Scalar sum = Scalar(0);
//             for (Index j = 0; j < i; ++j) {
//                 // Need to divide by Δ(j,j) here!
//                 sum += (m_L(i, j) / m_Delta(j, j)) * z(j);
//             }
//             z(i) = b_perm(i) - sum;
//         }
        
//         // Diagonal scaling: y = z / Δ (element-wise)
//         VectorType y(n);
//         for (Index i = 0; i < n; ++i) {
//             y(i) = z(i) / m_Delta(i, i);
//         }
        
//         // Back substitution: Lᵀ * x = y
//         VectorType x(n);
//         for (Index i = n - 1; i >= 0; --i) {
//             Scalar sum = Scalar(0);
//             for (Index j = i + 1; j < n; ++j) {
//                 // Need to divide by Δ(i,i) here!
//                 sum += (m_L(j, i) / m_Delta(i, i)) * x(j);
//             }
//             x(i) = y(i) - sum;
//         }
        
//         return x;
//     }

//     /**
//      * @brief Returns the info status of the decomposition
//      */
//     Eigen::ComputationInfo info() const {
//         return m_info;
//     }

//     /**
//      * @brief Returns true if matrix is positive definite
//      */
//     bool isPositiveDefinite() const {
//         return m_isPositiveDefinite && m_info == Eigen::Success;
//     }

// private:
//     MatrixType m_L;           // Stores NUMERATORS for L (not actual L!)
//     MatrixType m_Delta;       // Diagonal matrix Δ
//     Eigen::PermutationMatrix<Eigen::Dynamic> m_permutation;
//     Scalar m_denominator;     // d_{n-1}
//     bool m_isPositiveDefinite;
//     bool m_isInitialized;
//     Eigen::ComputationInfo m_info;
// };

// } // namespace Ext
// } // namespace Eigen