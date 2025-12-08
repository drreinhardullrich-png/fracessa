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
class bitset64 {
public:
    // --------------------------
    // ctor / copy / move
    // --------------------------
    FORCE_INLINE bitset64() noexcept : bits_(0ULL) {}
    FORCE_INLINE bitset64(uint64_t bits) noexcept : bits_(bits) {}
    FORCE_INLINE bitset64(const bitset64& o) noexcept = default;
    FORCE_INLINE bitset64(bitset64&& o) noexcept = default;
    FORCE_INLINE bitset64& operator=(const bitset64& o) noexcept = default;
    FORCE_INLINE bitset64& operator=(bitset64&& o) noexcept = default;

    //careful! no check for the biggest element here! done outside!
    FORCE_INLINE bitset64 next_bitset_with_same_popcount() const noexcept {
        uint64_t t = bits_ | (bits_ - 1);
        return bitset64((t + 1) | (((~t & -~t) - 1) >> (ctz64(bits_) + 1)));
    }


    // Helper function to iterate through all non-empty support sets
    // Calls callback for each bitset64 from 1 to (1<<nbits)-1
    template<typename F>
    static FORCE_INLINE void iterate_all_supports(unsigned nbits, F&& callback) {
        for (uint64_t bits = 1ull; bits < (1ull << nbits); bits++) {
            callback(bitset64(bits));
        }
    }

    // Iterate over all set bits, calling callback with each bit position
    // Returns false if callback returns false (early exit), true if all iterations complete
    // Use this for copositivity checks where early exit is needed
    template<typename F>
    FORCE_INLINE bool for_each_set_bit(F&& callback) const noexcept {
        for (unsigned i = find_first(); i < 64; i = find_next(i)) {
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
    FORCE_INLINE void for_each_set_bit_no_exit(F&& callback) const noexcept {
        for (unsigned i = find_first(); i < 64; i = find_next(i)) {
            callback(i);
        }
    }

    // --------------------------
    // clear / set all
    // --------------------------
    //FORCE_INLINE void clear() noexcept { bits_ = 0ULL; }
    FORCE_INLINE void set_all(unsigned nbits) noexcept { bits_ = (1ULL << nbits) - 1ULL; } //careful with nbits == 0

    // --------------------------
    // single-bit ops (NO bounds checks, NO nbits needed)
    // --------------------------
    FORCE_INLINE void set(unsigned pos) noexcept {
        bits_ |= (1ULL << pos);
    }
    FORCE_INLINE void reset(unsigned pos) noexcept {
        bits_ &= ~(1ULL << pos);
    }
    FORCE_INLINE bool test(unsigned pos) const noexcept {
        return (bits_ >> pos) & 1ULL;
    }

    // --------------------------
    // Return this \ o (set difference)
    FORCE_INLINE bitset64 subtract(const bitset64 &o) const noexcept {
        return bitset64(bits_ & ~o.bits_);
    }

    // Get the lowest set bit as a bitset64 (only that bit set)
    FORCE_INLINE bitset64 lowest_set_bit() const noexcept {
        if (bits_ == 0ULL) return bitset64(0ULL);
        unsigned pos = find_first();
        bitset64 result(0ULL);
        result.set(pos);
        return result;
    }

    // Get smallest representation by circular rotation (for circular symmetric matrices) with bitmask of size nbits
    FORCE_INLINE bitset64 smallest_representation(unsigned nbits) const noexcept {
        if (nbits == 0 || bits_ == 0ULL) return *this;
        uint64_t mask = (1ULL << nbits) - 1ULL;
        uint64_t min_val = bits_ & mask;
        
        // Early exit: 0 is absolute minimum
        if (min_val == 0ULL) return bitset64(0ULL);
        
        uint64_t current = min_val;
        for (unsigned i = 1; i < nbits; i++) {
            // Circular rotate right by 1
            uint64_t lo = current << (nbits - 1);
            uint64_t hi = current >> 1;
            current = (hi | lo) & mask;
            
            // Early exit: found absolute minimum (0 is smallest possible)
            if (current == 0ULL) return bitset64(0ULL);
            
            if (current < min_val) {
                min_val = current;
                // Early exit optimization: if we found a new minimum that's very small,
                // continue but the 0 check above will catch the absolute minimum
            }
        }
        return bitset64(min_val);
    }

    // Check if this bitset is in its smallest representation (canonical form)
    // Uses early exit: returns false immediately if any rotation is smaller than original
    // Much faster than smallest_representation() == *this for canonical checks
    // Optimized version: reduces masking operations for better performance
    FORCE_INLINE bool is_smallest_representation(unsigned nbits) const noexcept {
        uint64_t mask = (1ULL << nbits) - 1ULL;
        uint64_t original = bits_ & mask;
        
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
    // non-modifying bitwise ops (NO nbits needed for &, |, ^)
    // --------------------------
    FORCE_INLINE bitset64 operator&(const bitset64 &o) const noexcept { return bitset64(bits_ & o.bits_); }
    FORCE_INLINE bitset64 operator|(const bitset64 &o) const noexcept { return bitset64(bits_ | o.bits_); }
    FORCE_INLINE bitset64 operator^(const bitset64 &o) const noexcept { return bitset64(bits_ ^ o.bits_); }
    // Note: operator~() without nbits returns all bits set - use with care
    FORCE_INLINE bitset64 operator~() const noexcept { return bitset64(~bits_); }

    // --------------------------
    // tests (NO nbits needed)
    // --------------------------
    // this ⊆ o  <=> (this & ~o) == 0
    FORCE_INLINE bool is_subset_of(const bitset64 &o) const noexcept { return (bits_ & ~o.bits_) == 0ULL; }
    FORCE_INLINE bool operator==(const bitset64 &o) const noexcept { return bits_ == o.bits_; }
    FORCE_INLINE bool operator!=(const bitset64 &o) const noexcept { return bits_ != o.bits_; }
    FORCE_INLINE bool operator<=(const bitset64 &o) const noexcept { return bits_ <= o.bits_; }
    FORCE_INLINE bool operator<(const bitset64 &o) const noexcept { return bits_ < o.bits_; }
    FORCE_INLINE bool operator>=(const bitset64 &o) const noexcept { return bits_ >= o.bits_; }
    FORCE_INLINE bool operator>(const bitset64 &o) const noexcept { return bits_ > o.bits_; }

    // --------------------------
    // popcount / size of support (NO nbits needed)
    // --------------------------
    FORCE_INLINE unsigned count() const noexcept { return popcount64(bits_); }

    // --------------------------
    // find-first and next (NO nbits needed for find_first)
    // --------------------------
    // caller must check against nbits, be careful with nbits == 0
    FORCE_INLINE unsigned find_first() const noexcept {
        return ctz64(bits_);
    }

    // find next bit after pos
    FORCE_INLINE unsigned find_next(unsigned pos) const noexcept {
        unsigned p = pos + 1;
        uint64_t w = bits_ & (~0ULL << p); // zero below p
        if (w) return ctz64(w);
        return 64; // no more bits found
    }

    // circular rotate right by exactly 1 bit (in place) for a bitmask of size nbits!
    FORCE_INLINE void rot_one_right(unsigned nbits) noexcept {
        uint64_t mask = (1ULL << nbits) - 1ULL;
        uint64_t low = bits_ & mask;
        uint64_t lo = low << (nbits - 1);
        uint64_t hi = low >> 1;
        bits_ = (hi | lo) & mask;
    }
    
    // Convert to string representation (decimal representation of uint64)
    inline std::string to_string() const noexcept {
        return std::to_string(bits_);
    }

    // Convert to bitstring representation (MSB first, like std::bitset::to_string())
    // Only outputs bits 0 to dimension-1 (rightmost dimension bits)
    // Example: dimension=5, bits 0,3,4 set -> "10011"
    inline std::string to_bitstring(unsigned dimension) const noexcept {
        if (dimension == 0) return "";
        std::string result;
        result.reserve(dimension);
        // Output from highest bit (dimension-1) to lowest bit (0)
        for (int i = static_cast<int>(dimension) - 1; i >= 0; i--) {
            result += (test(static_cast<unsigned>(i)) ? '1' : '0');
        }
        return result;
    }

    // portable hash (fast)
    FORCE_INLINE std::size_t hash() const noexcept {
        // 64-bit mix (FNV-like)
        uint64_t x = bits_;
        x ^= x >> 33;
        x *= 0xff51afd7ed558ccdULL;
        x ^= x >> 33;
        x *= 0xc4ceb9fe1a85ec53ULL;
        x ^= x >> 33;
        return (std::size_t)x;
    }

private:
    uint64_t bits_;
};
