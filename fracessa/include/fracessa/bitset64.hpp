// bitset64.hpp
#pragma once
#include <cstdint>
#include <cstddef>
#include <string>

// Force-inline hint
#if defined(_MSC_VER)
#  define FORCE_INLINE __forceinline
#else
#  define FORCE_INLINE __attribute__((always_inline)) inline
#endif

// Platform-specific intrinsics
#ifdef _MSC_VER
#include <intrin.h>
#endif

// Portable popcount wrapper
FORCE_INLINE unsigned popcount64(uint64_t x) noexcept {
#ifdef _MSC_VER
    return (unsigned)_mm_popcnt_u64(x);
#else
    return (unsigned)__builtin_popcountll(x);
#endif
}

// Portable count trailing zeros wrapper
FORCE_INLINE unsigned ctz64(uint64_t x) noexcept {
#ifdef _MSC_VER
    unsigned long index;
    if (_BitScanForward64(&index, x)) {
        return (unsigned)index;
    }
    return 64; // undefined behavior case, but we check for 0 before calling
#else
    return (unsigned)__builtin_ctzll(x);
#endif
}

/// Ultra-optimized bitset for n <= 64.
/// - stores only uint64_t bits (8 bytes, same as uint64_t)
/// - nbits must be provided by caller when needed for masking
/// - all operations are inlined and branch-light
/// - NO bounds checks for maximum performance

// Type alias: bitset64 is just uint64_t
typedef uint64_t bitset64;

// Namespace for bitset64 operations
namespace bs64 {

    // --------------------------
    // Modifying operations (take bitset64&)
    // --------------------------
    
    FORCE_INLINE void set(bitset64& bits, unsigned pos) noexcept {
        bits |= (1ULL << pos);
    }
    
    FORCE_INLINE void reset(bitset64& bits, unsigned pos) noexcept {
        bits &= ~(1ULL << pos);
    }
    
    FORCE_INLINE void set_all(bitset64& bits, unsigned nbits) noexcept {
        bits = (1ULL << nbits) - 1ULL; //careful with nbits == 0
    }
    
    // circular rotate right by exactly 1 bit (in place) for a bitmask of size nbits!
    FORCE_INLINE void rot_one_right(bitset64& bits, unsigned nbits) noexcept {
        uint64_t mask = (1ULL << nbits) - 1ULL;
        uint64_t low = bits & mask;
        uint64_t lo = low << (nbits - 1);
        uint64_t hi = low >> 1;
        bits = (hi | lo) & mask;
    }

    // --------------------------
    // Non-modifying operations (take bitset64 by value)
    // --------------------------
    
    FORCE_INLINE bool test(bitset64 bits, unsigned pos) noexcept {
        return (bits >> pos) & 1ULL;
    }
    
    FORCE_INLINE unsigned count(bitset64 bits) noexcept {
        return popcount64(bits);
    }
    
    FORCE_INLINE unsigned find_first(bitset64 bits) noexcept {
        return ctz64(bits);
    }
    
    // find next bit after pos
    FORCE_INLINE unsigned find_next(bitset64 bits, unsigned pos) noexcept {
        unsigned p = pos + 1;
        uint64_t w = bits & (~0ULL << p); // zero below p
        if (w) return ctz64(w);
        return 64; // no more bits found
    }
    
    // this ⊆ o  <=> (this & ~o) == 0
    FORCE_INLINE bool is_subset_of(bitset64 bits, bitset64 o) noexcept {
        return (bits & ~o) == 0ULL;
    }
    
    // Return this \ o (set difference)
    FORCE_INLINE bitset64 subtract(bitset64 bits, bitset64 o) noexcept {
        return bits & ~o;
    }
    
    // Get the lowest set bit as a bitset64 (only that bit set)
    FORCE_INLINE bitset64 lowest_set_bit(bitset64 bits) noexcept {
        if (bits == 0ULL) return 0ULL;
        unsigned pos = find_first(bits);
        bitset64 result = 0ULL;
        set(result, pos);
        return result;
    }
    
    //careful! no check for the biggest element here! done outside!
    FORCE_INLINE bitset64 next_bitset_with_same_popcount(bitset64 bits) noexcept {
        uint64_t t = bits | (bits - 1);
        return (t + 1) | (((~t & -~t) - 1) >> (ctz64(bits) + 1));
    }
    
    // Get smallest representation by circular rotation (for circular symmetric matrices) with bitmask of size nbits
    FORCE_INLINE bitset64 smallest_representation(bitset64 bits, unsigned nbits) noexcept {
        if (nbits == 0 || bits == 0ULL) return bits;
        uint64_t mask = (1ULL << nbits) - 1ULL;
        uint64_t min_val = bits & mask;
        
        // Early exit: 0 is absolute minimum
        if (min_val == 0ULL) return 0ULL;
        
        uint64_t current = min_val;
        for (unsigned i = 1; i < nbits; i++) {
            // Circular rotate right by 1
            uint64_t lo = current << (nbits - 1);
            uint64_t hi = current >> 1;
            current = (hi | lo) & mask;
            
            // Early exit: found absolute minimum (0 is smallest possible)
            if (current == 0ULL) return 0ULL;
            
            if (current < min_val) {
                min_val = current;
                // Early exit optimization: if we found a new minimum that's very small,
                // continue but the 0 check above will catch the absolute minimum
            }
        }
        return min_val;
    }
    
    // Check if this bitset is in its smallest representation (canonical form)
    // Uses early exit: returns false immediately if any rotation is smaller than original
    // Much faster than smallest_representation() == *this for canonical checks
    // Optimized version: reduces masking operations for better performance
    FORCE_INLINE bool is_smallest_representation(bitset64 bits, unsigned nbits) noexcept {
        uint64_t mask = (1ULL << nbits) - 1ULL;
        uint64_t original = bits & mask;
        
        uint64_t current = original;
        unsigned shift_left = nbits - 1;
        
        for (unsigned i = 1; i < nbits; i++) {
            current = ((current >> 1) | (current << shift_left)) & mask;
            if (current < original) {
                return false; //there exists a smaller representation, ie. this one cannot be canonical!
            }
        }
        // All rotations are >= original, so it's canonical
        return true;
    }
    
    // --------------------------
    // Template functions
    // --------------------------
    
    // Helper function to iterate through all non-empty support sets
    // Calls callback for each bitset64 from 1 to (1<<nbits)-1
    template<typename F>
    FORCE_INLINE void iterate_all_supports(unsigned nbits, F&& callback) {
        for (uint64_t bits = 1ull; bits < (1ull << nbits); bits++) {
            callback(bits);
        }
    }
    
    // Iterate over all set bits, calling callback with each bit position
    // Returns false if callback returns false (early exit), true if all iterations complete
    // Use this for copositivity checks where early exit is needed
    template<typename F>
    FORCE_INLINE bool for_each_set_bit(bitset64 bits, F&& callback) noexcept {
        for (unsigned i = find_first(bits); i < 64; i = find_next(bits, i)) {
            if (!callback(i)) {
                return false; // Early exit
            }
        }
        return true; // All iterations completed
    }
    
    // Iterate over all set bits, calling callback with each bit position
    // Always completes all iterations (no early exit)
    // Use this for checkstab.cpp and matrix operations where all bits must be processed
    template<typename F>
    FORCE_INLINE void for_each_set_bit_no_exit(bitset64 bits, F&& callback) noexcept {
        for (unsigned i = find_first(bits); i < 64; i = find_next(bits, i)) {
            callback(i);
        }
    }
    
    // --------------------------
    // String functions
    // --------------------------
    
    // Convert to string representation (decimal representation of uint64)
    inline std::string to_string(bitset64 bits) noexcept {
        return std::to_string(bits);
    }
    
    // Convert to bitstring representation (MSB first, like std::bitset::to_string())
    // Only outputs bits 0 to dimension-1 (rightmost dimension bits)
    // Example: dimension=5, bits 0,3,4 set -> "10011"
    inline std::string to_bitstring(bitset64 bits, unsigned dimension) noexcept {
        if (dimension == 0) return "";
        std::string result;
        result.reserve(dimension);
        // Output from highest bit (dimension-1) to lowest bit (0)
        for (int i = static_cast<int>(dimension) - 1; i >= 0; i--) {
            result += (test(bits, static_cast<unsigned>(i)) ? '1' : '0');
        }
        return result;
    }
    
    // --------------------------
    // Hash function
    // --------------------------
    
    // portable hash (fast)
    FORCE_INLINE std::size_t hash(bitset64 bits) noexcept {
        // 64-bit mix (FNV-like)
        uint64_t x = bits;
        x ^= x >> 33;
        x *= 0xff51afd7ed558ccdULL;
        x ^= x >> 33;
        x *= 0xc4ceb9fe1a85ec53ULL;
        x ^= x >> 33;
        return (std::size_t)x;
    }
    
    // --------------------------
    // Operators
    // --------------------------
    // Note: Operators (&, |, ^, ~, ==, !=, <, <=, >, >=) work natively on uint64_t,
    // so no need to define them here since bitset64 is just a typedef for uint64_t
    
} // namespace bs64
