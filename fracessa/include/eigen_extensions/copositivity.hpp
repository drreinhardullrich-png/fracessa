#include <cmath>
#include <unordered_map>
#include <Eigen/Dense>
#include <fracessa/bitset64.hpp>
#include <fracessa/matrix.hpp>

// Hash function for bitset64 to use with std::unordered_map
struct bitset64_hash {
    std::size_t operator()(const bitset64& bs) const noexcept {
        return bs.hash();
    }
};

// 0 = False, 1 = True, -1 = Unknown
static std::unordered_map<bitset64, int8_t, bitset64_hash> memo;

inline bool checkRecursive(const Eigen::MatrixXd& A, const bitset64& mask, int n) {
    // 1. Check Cache
    auto it = memo.find(mask);
    if (it != memo.end() && it->second != -1) {
        return it->second == 1;
    }

    size_t current_dim = mask.count();

    // 2. Base Case: 1x1 Matrix
    if (current_dim == 1) {
        unsigned idx = mask.find_first();
        bool result = A(static_cast<Eigen::Index>(idx), static_cast<Eigen::Index>(idx)) > 0;
        memo[mask] = result ? 1 : 0;
        return result;
    }

    // 3. Recursive Step: Check all submatrices of size (current_dim - 1)
    // We generate them by turning off one bit at a time from the current mask.
    for (int i = 0; i < n; ++i) {
        // If index i is in the current set
        if (mask.test(i)) {
            bitset64 sub_mask = mask;
            sub_mask.reset(i); // Turn off bit i
            if (!checkRecursive(A, sub_mask, n)) {
                memo[mask] = 0;
                return false;
            }
        }
    }

    // 4. If all submatrices are valid, check Determinant condition
    // Construct the actual submatrix for this mask using principal_submatrix
    Eigen::MatrixXd subMat;
    matrix_ops::principal_submatrix(A, static_cast<size_t>(n), mask, current_dim, subMat);
    
    double det = subMat.determinant();
    const double EPSILON = 1e-9;

    if (det < -EPSILON) {
        // Check Adjugate
        Eigen::MatrixXd adj = subMat.inverse() * det; 
        bool all_positive = true;
        for (int r = 0; r < static_cast<int>(current_dim); ++r) {
            for (int c = 0; c < static_cast<int>(current_dim); ++c) {
                if (adj(r, c) <= EPSILON) {
                    all_positive = false; 
                    break;
                }
            }
            if (!all_positive) break;
        }

        if (all_positive) {
            memo[mask] = 0; // Not strict copositive
            return false;
        }
    }

    // Passed all checks
    memo[mask] = 1;
    return true;
}

inline bool isStrictlyCopositiveMemoized(const Eigen::MatrixXd& A) {
    int n = A.rows();
    // Clear memo table for new computation
    memo.clear();
    // Start recursion with all bits set (111...1) representing the full matrix
    bitset64 full_mask;
    full_mask.set_all(n);
    return checkRecursive(A, full_mask, n);
}
