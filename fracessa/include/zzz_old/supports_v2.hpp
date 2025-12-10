// supports_v2.hpp
#pragma once
#include <cstddef>
#include <vector>
#include <algorithm>
#include <fracessa/bitset64.hpp>
#include <boost/math/special_functions/binomial.hpp>

// Force-inline hint (same as bitset64.hpp)
#if defined(_MSC_VER)
#  define FORCE_INLINE __forceinline
#else
#  define FORCE_INLINE __attribute__((always_inline)) inline
#endif

/// High-performance Supports class for managing support sets
/// All hot-path methods are FORCE_INLINE for maximum performance
class Supports {
private:
    std::vector<std::vector<bitset64>> supports_;
    size_t dimension_;
    bool is_cs_;
    
    /// Calculate combinadic rank of a bitset x with popcount k in space of n bits
    /// Returns the zero-based position in the lexicographically sorted list
    FORCE_INLINE size_t combinadic_rank(const bitset64& x, size_t n, size_t k) const noexcept {
        size_t rank = 0;
        size_t remaining_k = k;
        
        // Scan bit positions from n-1 down to 0
        for (size_t i = n; i > 0 && remaining_k > 0; ) {
            --i;
            if (x.test(static_cast<unsigned>(i))) {
                // Add count of all combinations smaller than x at this position
                // This is binom(i, remaining_k) - all numbers with '0' at position i
                // that use remaining_k ones from the remaining i lower positions
                if (i >= remaining_k) {
                    rank += static_cast<size_t>(
                        boost::math::binomial_coefficient<double>(i, remaining_k)
                    );
                }
                --remaining_k;
            }
        }
        return rank;
    }
    
    /// Find the smallest bitset with popcount k that is greater than target (integer comparison)
    /// This finds the lexicographically smallest combination with k bits that has integer value > target
    FORCE_INLINE bitset64 find_smallest_bitset_gt(size_t k, const bitset64& target) const noexcept {
        // Start with the smallest possible combination with k bits: (1 << k) - 1
        uint64_t min_combination = (1ULL << k) - 1ULL;
        bitset64 candidate(min_combination);
        
        // If the smallest combination is already > target, return it
        if (candidate > target) {
            return candidate;
        }
        
        // Otherwise, iterate through combinations with k bits to find the first one > target
        // We use next_bitset_with_same_popcount to iterate efficiently
        bitset64 current = candidate;
        uint64_t max_valid = (1ULL << dimension_) - 1ULL;
        uint64_t max_combination = 0;
        for (size_t j = 0; j < k; ++j) {
            max_combination |= (1ULL << (dimension_ - 1 - j));
        }
        
        // Iterate through combinations until we find one > target or exhaust all possibilities
        while (current <= target) {
            bitset64 next = current.next_bitset_with_same_popcount();
            
            // Check if we've gone beyond valid range or lost the popcount
            if (next.count() != k || (next > bitset64(max_valid))) {
                // Return maximum combination as fallback
                return bitset64(max_combination);
            }
            
            current = next;
            
            // Safety check: if we've reached the maximum, return it
            if (current == bitset64(max_combination)) {
                return current;
            }
        }
        
        return current;
    }
    
public:
    /// Constructor - stores parameters, does not initialize supports
    FORCE_INLINE Supports(size_t dimension, bool is_cs) noexcept
        : supports_(dimension)
        , dimension_(dimension)
        , is_cs_(is_cs)
    {}
    
    /// Initialize all supports - must be called after construction
    inline void initialize() {
        // Reserve space for each support size using binomial coefficients and set coprime flags
        std::vector<bool> is_coprime(dimension_);
        for (size_t i = 0; i < dimension_; ++i) {
            supports_[i].reserve(static_cast<uint64_t>(
                boost::math::binomial_coefficient<double>(dimension_, i + 1)
            ));            
            if (is_cs_) 
                    is_coprime[i] = (boost::integer::gcd(i+1, dimension_) == 1);
        }
        
        // Populate supports based on is_cs_ flag
        if (is_cs_) {
            bitset64::iterate_all_supports(dimension_, [&](const bitset64& support) {
                size_t current_index = support.count() - 1;
                if (is_coprime[current_index]) {
                    // Only add if it's the smallest representation (canonical form)
                    if (support.is_smallest_representation(dimension_)) {
                        supports_[current_index].push_back(support);
                    }
                } else {
                    supports_[current_index].push_back(support);
                }
            });
        } else {
            bitset64::iterate_all_supports(dimension_, [&](const bitset64& support) {
                supports_[support.count() - 1].push_back(support);
            });
        }
    }
    
    /// Get const reference to supports for a given support size (1-indexed)
    /// CRITICAL hot path - must be FORCE_INLINE
    FORCE_INLINE const std::vector<bitset64>& get_supports(size_t support_size) const noexcept {
        return supports_[support_size-1];
    }
    
    /// Remove all supersets of the given subset, starting from from_size
    /// Hot path - FORCE_INLINE for maximum performance
    /// Uses Combinadic Rank Calculation as upper bound, then binary search
    FORCE_INLINE void remove_supersets(const bitset64& subset, uint64_t support_size = 0) noexcept {
        if (support_size == 0) {
            support_size = subset.count();
        }
        for (size_t i = support_size; i < dimension_; ++i) { //index support_size means erase from support_size+1 on!!!
            auto& vec = supports_[i];
            if (vec.empty()) continue;
            
            // Elements in supports_[i] have popcount (i+1)
            size_t target_popcount = i + 1;
            
            // Strategy: Use combinadic rank to get an upper bound, then binary search within that bound
            // 1. Find the smallest bitset with popcount (i+1) that is > subset
            // 2. Calculate its combinadic rank - this gives us an upper bound (since elements may have been removed)
            // 3. Use binary search from beginning to that rank to find the actual first element > subset
            
            // Find target bitset: smallest with popcount (i+1) that is > subset
            bitset64 target = find_smallest_bitset_gt(target_popcount, subset);
            
            // Calculate combinadic rank of target - this is an upper bound
            size_t upper_bound_rank = combinadic_rank(target, dimension_, target_popcount);
            
            // Clamp upper bound to vector size
            size_t search_end = std::min(upper_bound_rank, vec.size());
            
            // Binary search from beginning to upper_bound_rank to find first element > subset
            auto start_it = std::upper_bound(
                vec.begin(),
                vec.begin() + search_end,
                subset
            );
            
            // Only check elements from start_it onwards
            if (start_it != vec.end()) {
                vec.erase(
                    std::remove_if(
                        start_it,
                        vec.end(),
                        [=](const bitset64& x) { return subset.is_subset_of(x); }
                    ),
                    vec.end()
                );
            }
        }
    }
    
    // /// Get size of supports for a given support size
    // FORCE_INLINE size_t size(size_t support_size) const noexcept {
    //     return supports_[support_size].size();
    // }
    
    // /// Check if supports for a given size are empty
    // FORCE_INLINE bool empty(size_t support_size) const noexcept {
    //     return supports_[support_size].empty();
    // }
};

