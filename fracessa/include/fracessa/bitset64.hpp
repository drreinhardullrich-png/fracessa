// bitset64.hpp
#pragma once
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <functional>

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

/// Ultra-optimized runtime-sized bitset for n <= 64.
/// - nbits is fixed at construction time (0..64)
/// - stores bits in a single uint64_t (low bits used)
/// - all operations are inlined and branch-light
class bitset64 {
public:
    // --------------------------
    // ctor / copy / move
    // --------------------------
    explicit FORCE_INLINE bitset64(unsigned nbits = 64) noexcept
        : nbits_((nbits > 64) ? 64u : nbits), bits_(0ULL)
    {
        mask_last_word();
    }

    FORCE_INLINE bitset64(const bitset64& o) noexcept = default;
    FORCE_INLINE bitset64(bitset64&& o) noexcept = default;
    FORCE_INLINE bitset64& operator=(const bitset64& o) noexcept = default;
    FORCE_INLINE bitset64& operator=(bitset64&& o) noexcept = default;

    // construct from raw mask (low bits used). nbits must be <=64.
    static FORCE_INLINE bitset64 from_mask(unsigned nbits, uint64_t mask) noexcept {
        bitset64 b(nbits);
        b.bits_ = mask & b.last_mask_;
        return b;
    }

    // Helper function to iterate through all non-empty support sets
    // Calls callback for each bitset64 from 1 to (1<<nbits)-1
    // Optimized: reuses a single bitset64 object, updating it with set_mask()
    template<typename F>
    static FORCE_INLINE void iterate_all_supports(unsigned nbits, F&& callback) {
        bitset64 support(nbits);  // Create once, reuse
        for (uint64_t mask = 1ull; mask < (1ull << nbits); mask++) {
            support.set_mask(mask);  // Fast: just bits_ = mask & last_mask_
            callback(support);       // Pass by const reference
        }
    }

    // --------------------------
    // info
    // --------------------------
    FORCE_INLINE unsigned size() const noexcept { return nbits_; }
    FORCE_INLINE bool empty() const noexcept { return nbits_ == 0; }
    FORCE_INLINE uint64_t mask() const noexcept { return last_mask_; }

    // --------------------------
    // clear / set all
    // --------------------------
    FORCE_INLINE void clear() noexcept { bits_ = 0ULL; }
    FORCE_INLINE void set_all() noexcept { bits_ = last_mask_; }

    // --------------------------
    // single-bit ops
    // --------------------------
    FORCE_INLINE void set(unsigned pos) noexcept {
        assert(pos < nbits_);
        bits_ |= (1ULL << pos);
    }
    FORCE_INLINE void reset(unsigned pos) noexcept {
        assert(pos < nbits_);
        bits_ &= ~(1ULL << pos);
    }
    FORCE_INLINE void flip(unsigned pos) noexcept {
        assert(pos < nbits_);
        bits_ ^= (1ULL << pos);
    }
    FORCE_INLINE bool test(unsigned pos) const noexcept {
        assert(pos < nbits_);
        return (bits_ >> pos) & 1ULL;
    }

    // --------------------------
    // bulk ops (in-place)
    // --------------------------
    FORCE_INLINE void inplace_and(const bitset64 &o) noexcept { bits_ &= o.bits_; }
    FORCE_INLINE void inplace_or (const bitset64 &o) noexcept { bits_ |= o.bits_; }
    FORCE_INLINE void inplace_xor(const bitset64 &o) noexcept { bits_ ^= o.bits_; }
    // this := this \ o
    FORCE_INLINE void inplace_subtract(const bitset64 &o) noexcept { bits_ &= ~o.bits_; bits_ &= last_mask_; }

    // Return this \ o (set difference)
    FORCE_INLINE bitset64 subtract(const bitset64 &o) const noexcept { return from_mask(nbits_, bits_ & ~o.bits_); }

    // Get the lowest set bit as a bitset64 (only that bit set)
    FORCE_INLINE bitset64 lowest_set_bit() const noexcept {
        if (bits_ == 0ULL) return bitset64(nbits_);
        unsigned pos = find_first();
        bitset64 result(nbits_);
        result.set(pos);
        return result;
    }

    // Get smallest representation by circular rotation (for circular symmetric matrices)
    FORCE_INLINE bitset64 smallest_representation() const noexcept {
        if (nbits_ == 0 || bits_ == 0ULL) return *this;
        bitset64 min_val = *this;
        bitset64 current = *this;
        for (unsigned i = 1; i < nbits_; i++) {
            current = current.rot_r(1);
            if (current.to_uint64() < min_val.to_uint64()) {
                min_val = current;
            }
        }
        return min_val;
    }

    // --------------------------
    // non-modifying bitwise ops
    // --------------------------
    FORCE_INLINE bitset64 operator&(const bitset64 &o) const noexcept { return from_mask(nbits_, bits_ & o.bits_); }
    FORCE_INLINE bitset64 operator|(const bitset64 &o) const noexcept { return from_mask(nbits_, bits_ | o.bits_); }
    FORCE_INLINE bitset64 operator^(const bitset64 &o) const noexcept { return from_mask(nbits_, bits_ ^ o.bits_); }
    FORCE_INLINE bitset64 operator~() const noexcept { return from_mask(nbits_, (~bits_) & last_mask_); }

    // --------------------------
    // tests
    // --------------------------
    FORCE_INLINE bool any() const noexcept { return bits_ != 0ULL; }
    FORCE_INLINE bool none() const noexcept { return bits_ == 0ULL; }
    // this ⊆ o  <=> (this & ~o) == 0
    FORCE_INLINE bool is_subset_of(const bitset64 &o) const noexcept { return (bits_ & ~o.bits_) == 0ULL; }
    FORCE_INLINE bool operator==(const bitset64 &o) const noexcept { return nbits_ == o.nbits_ && bits_ == o.bits_; }
    FORCE_INLINE bool operator!=(const bitset64 &o) const noexcept { return !(*this == o); }

    // --------------------------
    // popcount / size of support
    // --------------------------
    FORCE_INLINE unsigned count() const noexcept { return popcount64(bits_); }

    // --------------------------
    // find-first and next
    // --------------------------
    // returns nbits_ if none
    FORCE_INLINE unsigned find_first() const noexcept {
        if (bits_ == 0ULL) return nbits_;
        return ctz64(bits_);
    }

    // find next bit after pos; returns nbits_ if none
    FORCE_INLINE unsigned find_next(unsigned pos) const noexcept {
        assert(pos < nbits_);
        unsigned p = pos + 1;
        if (p >= nbits_) return nbits_;
        uint64_t w = bits_ & (~0ULL << p); // zero below p
        if (w) return ctz64(w);
        return nbits_;
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
    // shifts
    // --------------------------
    // logical left shift by s (fills low bits with 0)
    FORCE_INLINE bitset64 shl(unsigned s) const noexcept {
        if (s == 0 || bits_ == 0ULL) return *this;
        if (s >= nbits_) return bitset64(nbits_);
        uint64_t r = shl64_safe(bits_, s) & last_mask_;
        return from_mask(nbits_, r);
    }

    // logical right shift by s (fills high bits with 0)
    FORCE_INLINE bitset64 shr(unsigned s) const noexcept {
        if (s == 0 || bits_ == 0ULL) return *this;
        if (s >= nbits_) return bitset64(nbits_);
        uint64_t r = shr64_safe(bits_, s) & last_mask_;
        return from_mask(nbits_, r);
    }

    // circular rotate right by s (0 <= s < nbits_)
    FORCE_INLINE bitset64 rot_r(unsigned s) const noexcept {
        if (nbits_ == 0) return *this;
        s %= nbits_;
        if (s == 0) return *this;
        // rotate within nbits_ domain: do rotation on masked bits
        uint64_t low = bits_ & last_mask_;
        uint64_t lo = shl64_safe(low, nbits_ - s);
        uint64_t hi = shr64_safe(low, s);
        return from_mask(nbits_, (hi | lo) & last_mask_);
    }

    // circular rotate left by s
    FORCE_INLINE bitset64 rot_l(unsigned s) const noexcept {
        if (nbits_ == 0) return *this;
        s %= nbits_;
        if (s == 0) return *this;
        return rot_r(nbits_ - s);
    }

    // convenience: return the underlying mask (low bits contain data)
    FORCE_INLINE uint64_t to_uint64() const noexcept { return bits_ & last_mask_; }

    // set from a raw mask (low bits used)
    FORCE_INLINE void set_mask(uint64_t mask) noexcept { bits_ = mask & last_mask_; }

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
    unsigned nbits_;     // 0..64
    uint64_t bits_;      // low nbits_ used
    uint64_t last_mask_; // mask with nbits_ low bits set

    FORCE_INLINE void mask_last_word() noexcept {
        if (nbits_ >= 64) last_mask_ = ~0ULL;
        else if (nbits_ == 0) last_mask_ = 0ULL;
        else last_mask_ = (1ULL << nbits_) - 1ULL;
        bits_ &= last_mask_;
    }
};
