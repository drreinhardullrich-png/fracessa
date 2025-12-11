// supports.hpp
#pragma once
#include <cstddef>
#include <vector>
#include <algorithm>
#include <fracessa/bitset64.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/integer/common_factor.hpp>

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
            for (uint64_t bits = 1ULL; bits < (1ULL << dimension_); ++bits) {
                bitset64 support = bits;
                size_t current_index = bs64::count(support) - 1;
                if (is_coprime[current_index]) {
                    // Only add if it's the smallest representation (canonical form)
                    if (bs64::is_smallest_representation(support, dimension_)) {
                        supports_[current_index].push_back(support);
                    }
                } else {
                    supports_[current_index].push_back(support);
                }
            }
        } else {
            for (uint64_t bits = 1ULL; bits < (1ULL << dimension_); ++bits) {
                bitset64 support = bits;
                supports_[bs64::count(support) - 1].push_back(support);
            }
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
            support_size = bs64::count(subset);
        }
        for (size_t i = support_size; i < dimension_; ++i) { //index support_size means erase from support_size+1 on!!!
            auto& vec = supports_[i];
            // Binary search: find first element where x > subset (only for dimension >= 10)
            // this is an educated guess. for dimensions the overhead of the binary search is bigger than the time saved
            auto start_it = (dimension_ >= 10) 
                ? std::upper_bound(vec.begin(), vec.end(), subset)
                : vec.begin();
            // Only check elements from start_it onwards
            if (start_it != vec.end()) {
                vec.erase(
                std::remove_if(
                        start_it,
                        vec.end(),
                    [=](const bitset64& x) { return bs64::is_subset_of(subset, x); }
                ),
                    vec.end()
            );
            }
        }
    }
    
    /// Remove all supersets for a batch of subsets
    /// More efficient than calling remove_supersets multiple times NOT TRUE!!!!!!
    FORCE_INLINE void remove_supersets_batch(const std::vector<bitset64>& subset_list, size_t support_size) noexcept {
        for (size_t i = support_size; i < dimension_; ++i) { //index support_size means erase from support_size+1 on!!!
            auto& vec = supports_[i];
            vec.erase(
                std::remove_if(
                    vec.begin(),
                    vec.end(),
                    [&](const bitset64& x) { 
                        return std::any_of(
                            subset_list.begin(), 
                            subset_list.end(), 
                            [&](const bitset64& s) { 
                                return bs64::is_subset_of(s, x); 
                            }
                        ); 
                    }
                ),
                vec.end()
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

