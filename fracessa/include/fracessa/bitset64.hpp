// bitset64.hpp
#pragma once
#include <cstdint>
#include <cstddef>
#include <functional>
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

// Very small utility for safe shift: defined behavior for shift >= 64 -> returns 0
FORCE_INLINE constexpr uint64_t shl64_safe(uint64_t x, unsigned s) noexcept {
    return (s >= 64) ? 0ULL : (x << s);
}
FORCE_INLINE constexpr uint64_t shr64_safe(uint64_t x, unsigned s) noexcept {
    return (s >= 64) ? 0ULL : (x >> s);
}

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

// Helper to compute mask for nbits (low nbits set)
FORCE_INLINE constexpr uint64_t compute_mask(unsigned nbits) noexcept {
    if (nbits >= 64) return ~0ULL;
    if (nbits == 0) return 0ULL;
    return (1ULL << nbits) - 1ULL;
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

    // construct from raw mask (low bits used). nbits must be <=64.
    static FORCE_INLINE bitset64 from_mask(unsigned nbits, uint64_t mask) noexcept {
        return bitset64(mask & compute_mask(nbits));
    }

    // Helper function to iterate through all non-empty support sets
    // Calls callback for each bitset64 from 1 to (1<<nbits)-1
    // Optimized: reuses a single bitset64 object, updating it with set_mask()
    template<typename F>
    static FORCE_INLINE void iterate_all_supports(unsigned nbits, F&& callback) {
        bitset64 support(0ULL);
        uint64_t mask = compute_mask(nbits);
        for (uint64_t bits = 1ull; bits < (1ull << nbits); bits++) {
            support.bits_ = bits & mask;
            callback(support, nbits);  // Pass nbits to callback
        }
    }

    // --------------------------
    // clear / set all
    // --------------------------
    FORCE_INLINE void clear() noexcept { bits_ = 0ULL; }
    FORCE_INLINE void set_all(unsigned nbits) noexcept { bits_ = compute_mask(nbits); }

    // --------------------------
    // single-bit ops (NO bounds checks, NO nbits needed)
    // --------------------------
    FORCE_INLINE void set(unsigned pos) noexcept {
        bits_ |= (1ULL << pos);
    }
    FORCE_INLINE void reset(unsigned pos) noexcept {
        bits_ &= ~(1ULL << pos);
    }
    FORCE_INLINE void flip(unsigned pos) noexcept {
        bits_ ^= (1ULL << pos);
    }
    FORCE_INLINE bool test(unsigned pos) const noexcept {
        return (bits_ >> pos) & 1ULL;
    }

    // --------------------------
    // bulk ops (in-place, NO nbits needed)
    // --------------------------
    FORCE_INLINE void inplace_and(const bitset64 &o) noexcept { bits_ &= o.bits_; }
    FORCE_INLINE void inplace_or (const bitset64 &o) noexcept { bits_ |= o.bits_; }
    FORCE_INLINE void inplace_xor(const bitset64 &o) noexcept { bits_ ^= o.bits_; }
    // this := this \ o
    FORCE_INLINE void inplace_subtract(const bitset64 &o) noexcept { bits_ &= ~o.bits_; }

    // Return this \ o (set difference)
    FORCE_INLINE bitset64 subtract(const bitset64 &o, unsigned nbits) const noexcept {
        return bitset64((bits_ & ~o.bits_) & compute_mask(nbits));
    }

    // Get the lowest set bit as a bitset64 (only that bit set)
    FORCE_INLINE bitset64 lowest_set_bit() const noexcept {
        if (bits_ == 0ULL) return bitset64(0ULL);
        unsigned pos = find_first();
        bitset64 result(0ULL);
        result.set(pos);
        return result;
    }

    // Get smallest representation by circular rotation (for circular symmetric matrices)
    FORCE_INLINE bitset64 smallest_representation(unsigned nbits) const noexcept {
        if (nbits == 0 || bits_ == 0ULL) return *this;
        uint64_t mask = compute_mask(nbits);
        uint64_t min_val = bits_ & mask;
        uint64_t current = bits_ & mask;
        for (unsigned i = 1; i < nbits; i++) {
            // Circular rotate right by 1
            uint64_t low = current & mask;
            uint64_t lo = shl64_safe(low, nbits - 1);
            uint64_t hi = shr64_safe(low, 1);
            current = (hi | lo) & mask;
            if (current < min_val) {
                min_val = current;
            }
        }
        return bitset64(min_val);
    }

    // --------------------------
    // non-modifying bitwise ops (NO nbits needed for &, |, ^)
    // --------------------------
    FORCE_INLINE bitset64 operator&(const bitset64 &o) const noexcept { return bitset64(bits_ & o.bits_); }
    FORCE_INLINE bitset64 operator|(const bitset64 &o) const noexcept { return bitset64(bits_ | o.bits_); }
    FORCE_INLINE bitset64 operator^(const bitset64 &o) const noexcept { return bitset64(bits_ ^ o.bits_); }
    // Note: operator~() without nbits returns all bits set - use with care
    FORCE_INLINE bitset64 operator~() const noexcept { return bitset64(~bits_); }
    // Complement with masking
    FORCE_INLINE bitset64 complement(unsigned nbits) const noexcept { return bitset64((~bits_) & compute_mask(nbits)); }

    // --------------------------
    // tests (NO nbits needed)
    // --------------------------
    FORCE_INLINE bool any() const noexcept { return bits_ != 0ULL; }
    FORCE_INLINE bool none() const noexcept { return bits_ == 0ULL; }
    // this ⊆ o  <=> (this & ~o) == 0
    FORCE_INLINE bool is_subset_of(const bitset64 &o) const noexcept { return (bits_ & ~o.bits_) == 0ULL; }
    FORCE_INLINE bool operator==(const bitset64 &o) const noexcept { return bits_ == o.bits_; }
    FORCE_INLINE bool operator!=(const bitset64 &o) const noexcept { return bits_ != o.bits_; }

    // --------------------------
    // popcount / size of support (NO nbits needed)
    // --------------------------
    FORCE_INLINE unsigned count() const noexcept { return popcount64(bits_); }

    // --------------------------
    // find-first and next (NO nbits needed for find_first)
    // --------------------------
    // returns 64 if none (caller must check against nbits)
    FORCE_INLINE unsigned find_first() const noexcept {
        if (bits_ == 0ULL) return 64;
        return ctz64(bits_);
    }

    // find next bit after pos; returns 64 if none
    FORCE_INLINE unsigned find_next(unsigned pos) const noexcept {
        unsigned p = pos + 1;
        if (p >= 64) return 64;
        uint64_t w = bits_ & (~0ULL << p); // zero below p
        if (w) return ctz64(w);
        return 64;
    }

    // iterate set bits: f(unsigned idx)
    template<typename F>
    FORCE_INLINE void for_each_set(F &&f) const {
        uint64_t v = bits_;
        while (v) {
            unsigned tz = ctz64(v);
            f(tz);
            v &= v - 1; // clear lowest set
        }
    }

    // --------------------------
    // shifts (need nbits for masking)
    // --------------------------
    // logical left shift by s (fills low bits with 0)
    FORCE_INLINE bitset64 shl(unsigned s, unsigned nbits) const noexcept {
        if (s == 0 || bits_ == 0ULL) return *this;
        if (s >= nbits) return bitset64(0ULL);
        uint64_t mask = compute_mask(nbits);
        uint64_t r = (shl64_safe(bits_, s) & mask);
        return bitset64(r);
    }

    // logical right shift by s (fills high bits with 0)
    FORCE_INLINE bitset64 shr(unsigned s, unsigned nbits) const noexcept {
        if (s == 0 || bits_ == 0ULL) return *this;
        if (s >= nbits) return bitset64(0ULL);
        uint64_t mask = compute_mask(nbits);
        uint64_t r = (shr64_safe(bits_, s) & mask);
        return bitset64(r);
    }

    // circular rotate right by s (0 <= s < nbits)
    FORCE_INLINE bitset64 rot_r(unsigned s, unsigned nbits) const noexcept {
        if (nbits == 0) return *this;
        s %= nbits;
        if (s == 0) return *this;
        uint64_t mask = compute_mask(nbits);
        uint64_t low = bits_ & mask;
        uint64_t lo = shl64_safe(low, nbits - s);
        uint64_t hi = shr64_safe(low, s);
        return bitset64((hi | lo) & mask);
    }

    // circular rotate left by s
    FORCE_INLINE bitset64 rot_l(unsigned s, unsigned nbits) const noexcept {
        if (nbits == 0) return *this;
        s %= nbits;
        if (s == 0) return *this;
        return rot_r(nbits - s, nbits);
    }

    // convenience: return the underlying mask (low bits contain data)
    FORCE_INLINE uint64_t to_uint64() const noexcept { return bits_; }
    FORCE_INLINE uint64_t to_uint64(unsigned nbits) const noexcept { return bits_ & compute_mask(nbits); }
    
    // Convert to string representation (decimal representation of uint64)
    inline std::string to_string() const noexcept {
        return std::to_string(bits_);
    }

    // set from a raw mask (low bits used)
    FORCE_INLINE void set_mask(uint64_t mask, unsigned nbits) noexcept { bits_ = mask & compute_mask(nbits); }

    // Convert to bitstring representation (MSB first, like std::bitset::to_string())
    // Only outputs bits 0 to dimension-1 (rightmost dimension bits)
    // Example: dimension=5, bits 0,3,4 set -> "10011"
    inline std::string to_bitstring(unsigned dimension) const noexcept {
        if (dimension == 0) return "";
        if (dimension > 64) dimension = 64;
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

    // Make bits_ accessible for iterate_all_supports optimization
    uint64_t bits_;
};
