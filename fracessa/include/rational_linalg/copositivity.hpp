#ifndef RATIONAL_LINALG_COPOSITIVITY_HPP
#define RATIONAL_LINALG_COPOSITIVITY_HPP

#include <rational_linalg/matrix.hpp>
#include <rational_linalg/bareiss_lu.hpp>
#include <rational_linalg/adjugate.hpp>
#include <fracessa/bitset64.hpp>
#include <unordered_map>

// Hash function for bitset64 to use with std::unordered_map
struct bitset64_hash {
    std::size_t operator()(const bitset64& bs) const noexcept {
        return bs64::hash(bs);
    }
};

namespace rational_linalg {

// Memoization cache: maps bitset64 mask to bool result
template<typename T>
class CopositivityChecker {
private:
    std::unordered_map<bitset64, bool, bitset64_hash> memo;

    // Recursive Check (Hadeler Criterion) for Rationals
    bool checkRecursive(const Matrix<T>& A, const bitset64& mask) {
        // 1. Check Cache
        auto it = memo.find(mask);
        if (it != memo.end()) {
            return it->second;
        }

        size_t current_dim = bs64::count(mask);

        // 2. Base Case: 1x1 Matrix
        if (current_dim == 1) {
            unsigned idx = bs64::find_first(mask);
            // Check if diagonal element > 0 (Rational comparison)
            bool result = A(static_cast<size_t>(idx), static_cast<size_t>(idx)) > T(0);
            memo[mask] = result;
            return result;
        }

        // 3. Recursive Step: Check all submatrices of size (current_dim - 1)
        // We iterate ONLY over the bits that are currently set to turn them off one by one
        if (!bs64::for_each_set_bit(mask, [&](unsigned i) {
            bitset64 sub_mask = mask;
            bs64::reset(sub_mask, i); // Turn off bit i representing row/col i
            return checkRecursive(A, sub_mask); // Return false to stop early if check fails
        })) {
            memo[mask] = false; // Fail early
            return false;
        }

        // 4. Determinant / Adjugate Check
        // If all proper principal submatrices are strictly copositive,
        // A is strictly copositive UNLESS (det(A) <= 0 AND adj(A) > 0)
        
        Matrix<T> subMat(current_dim, current_dim);
        // Use principal_submatrix which already uses optimized find_first()/find_next() iteration
        principal_submatrix(A, A.rows(), mask, current_dim, subMat);

        // Use BareissLUFactor for determinant and inverse
        BareissLUFactor<T> lu(subMat);
        T det = lu.determinant();

        if (det <= T(0)) {
            // Compute Adjugate: adj(A) = det(A) * A^(-1)
            // For singular matrices (det = 0), we can't use the inverse formula
            // because A^(-1) doesn't exist. We need to compute adjugate differently.
            Matrix<T> adj;
            if (lu.isSingular()) {
                // Fall back to cofactor method for singular matrices
                adj = adjugate(subMat);
            } else {
                // For det < 0, we can use the inverse formula, is cheaper!!!
                adj = det * lu.inverse();
            }

            // Check if Adjugate is Strictly Positive (> 0)
            if (all_entries_greater_zero(adj)) {
                memo[mask] = false; // Violates Hadeler condition
                return false;
            }
        }

        // Passed all checks
        memo[mask] = true;
        return true;
    }

public:
    // Main Entry Point
    bool isStrictlyCopositiveMemoized(const Matrix<T>& A) {
        size_t n = A.rows();
        // Clear cache for new computation
        memo.clear();
        
        // Reserve map size to avoid reallocations (heuristic: worst case 2^n subsets)
        memo.reserve(1ULL << n);
        
        // Create mask with lower n bits set
        bitset64 full_mask = 0ULL;
        bs64::set_all(full_mask, static_cast<unsigned>(n));

        return checkRecursive(A, full_mask);
    }
};

// Convenience function
template<typename T>
inline bool isStrictlyCopositiveMemoized(const Matrix<T>& A) {
    CopositivityChecker<T> checker;
    return checker.isStrictlyCopositiveMemoized(A);
}

} // namespace rational_linalg

#endif // RATIONAL_LINALG_COPOSITIVITY_HPP

