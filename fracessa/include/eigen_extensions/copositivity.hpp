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

// Memoization cache: maps bitset64 mask to bool result
static std::unordered_map<bitset64, bool, bitset64_hash> memo;

// Recursive Check (Hadeler Criterion) for Rationals
inline bool checkRecursive(const RationalMatrix& A, const bitset64& mask) {
    // 1. Check Cache
    auto it = memo.find(mask);
    if (it != memo.end()) {
        return it->second;
    }

    int current_dim = static_cast<int>(mask.count());

    // 2. Base Case: 1x1 Matrix
    if (current_dim == 1) {
        unsigned idx = mask.find_first();
        // Check if diagonal element > 0 (Rational comparison)
        bool result = A(static_cast<Eigen::Index>(idx), static_cast<Eigen::Index>(idx)) > rational(0);
        memo[mask] = result;
        return result;
    }

    // 3. Recursive Step: Check all submatrices of size (current_dim - 1)
    // We iterate ONLY over the bits that are currently set to turn them off one by one
    if (!mask.for_each_set_bit([&](unsigned i) {
        bitset64 sub_mask = mask;
        sub_mask.reset(i); // Turn off bit i representing row/col i
        return checkRecursive(A, sub_mask); // Return false to stop early if check fails
    })) {
        memo[mask] = false; // Fail early
        return false;
    }

    // 4. Determinant / Adjugate Check
    // If all proper principal submatrices are strictly copositive,
    // A is strictly copositive UNLESS (det(A) <= 0 AND adj(A) > 0)
    
    RationalMatrix subMat;
    // Use principal_submatrix which already uses optimized find_first()/find_next() iteration
    matrix_ops::principal_submatrix(A, static_cast<size_t>(A.rows()), mask, static_cast<size_t>(current_dim), subMat);

    // Use matrix_ops::determinant
    rational det = matrix_ops::determinant(subMat);

    if (det <= rational(0)) {
        // Compute Adjugate
        RationalMatrix adj = matrix_ops::adjugate(subMat);

        // Check if Adjugate is Strictly Positive (> 0)
        // Use matrix_ops::all_entries_greater_zero
        if (matrix_ops::all_entries_greater_zero(adj)) {
            memo[mask] = false; // Violates Hadeler condition
            return false;
        }
    }

    // Passed all checks
    memo[mask] = true;
    return true;
}

// Main Entry Point
inline bool isStrictlyCopositiveMemoized(const RationalMatrix& A) {
    int n = static_cast<int>(A.rows());
    // Clear cache for new computation
    memo.clear();
    
    // Reserve map size to avoid reallocations (heuristic: worst case 2^n subsets)
    memo.reserve(1 << n);
    
    // Create mask with lower n bits set
    bitset64 full_mask;
    full_mask.set_all(static_cast<unsigned>(n));

    return checkRecursive(A, full_mask);
}
