// supports.hpp
#pragma once
#include <sys/types.h>
#include <vector>
#include <algorithm>
#include <cstddef>
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
    FORCE_INLINE void remove_supersets(const bitset64& subset, uint64_t support_size = 0) noexcept {
        if (support_size == 0) {
            support_size = subset.count();
        }
        for (size_t i = support_size; i < dimension_; ++i) { //index support_size means erase from support_size+1 on!!!
            supports_[i].erase(
                std::remove_if(
                    supports_[i].begin(),
                    supports_[i].end(),
                    [=](const bitset64& x) { return subset.is_subset_of(x); }
                ),
                supports_[i].end()
            );
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

