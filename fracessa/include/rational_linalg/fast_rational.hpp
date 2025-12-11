#pragma once
#include <cstdint>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <utility>
#include <string>
#include <algorithm>

class rational_overflow : public std::runtime_error {
public:
    explicit rational_overflow(const char* msg)
        : std::runtime_error(msg) {}
};

class rational64 {
private:
    int64_t n_;   // numerator
    int64_t d_;   // denominator (> 0 always)

    // -------------------------------------------------------------
    // GCD: extremely fast binary-GCD for int64
    // Uses compiler intrinsics when available, fallback otherwise
    // -------------------------------------------------------------
    static inline int64_t gcd_int64(int64_t a, int64_t b) noexcept {
        // Handle negatives by taking absolute value
        // Note: Safe because we check for INT64_MIN separately if needed
        if (a < 0) a = -a;
        if (b < 0) b = -b;
        if (a == 0) return b;
        if (b == 0) return a;
        
        // Use compiler intrinsics for counting trailing zeros
        #if defined(__GNUC__) || defined(__clang__)
            int shift = __builtin_ctzll(static_cast<unsigned long long>(a | b));
            a >>= __builtin_ctzll(static_cast<unsigned long long>(a));
        do {
                b >>= __builtin_ctzll(static_cast<unsigned long long>(b));
            if (a > b) std::swap(a, b);
            b -= a;
        } while (b != 0);
        return a << shift;
        #else
            // Fallback: standard Euclidean algorithm (slower but portable)
            while (b != 0) {
                int64_t temp = b;
                b = a % b;
                a = temp;
            }
            return a;
        #endif
    }

    // -------------------------------------------------------------
    // Reduce fraction & normalize sign
    // -------------------------------------------------------------
    inline void normalize() {
        // Note: d_ is guaranteed non-zero by constructor check
        if (n_ == 0) {
            d_ = 1;
            return;
        }

        int64_t g = gcd_int64(n_, d_);
        n_ /= g;
        d_ /= g;

        // keep denominator positive
        // CRITICAL FIX: Handle INT64_MIN negation safely
        if (d_ < 0) {
            if (n_ == INT64_MIN) {
                // INT64_MIN cannot be negated (undefined behavior)
                // This case should be extremely rare, but we must handle it
                // If d_ is negative and n_ is INT64_MIN, we can't represent it
                // The only mathematically valid case is INT64_MIN / -1, which
                // would be INT64_MAX+1, but that overflows. So we throw.
                throw rational_overflow("INT64_MIN with negative denominator cannot be normalized");
            }
            n_ = -n_;
            d_ = -d_;
        }
    }

    // -------------------------------------------------------------
    // Builtin overflow-checked arithmetic helper
    // Uses compiler intrinsics when available, fallback otherwise
    // -------------------------------------------------------------
    static inline int64_t checked_mul(int64_t a, int64_t b) {
        #if defined(__GNUC__) || defined(__clang__)
        int64_t out;
        if (__builtin_mul_overflow(a, b, &out))
            throw rational_overflow("int64 overflow in multiplication");
        return out;
        #else
            // Fallback: manual overflow check
            // Check if a * b would overflow
            if (a > 0 && b > 0) {
                if (a > INT64_MAX / b) throw rational_overflow("int64 overflow in multiplication");
            } else if (a < 0 && b < 0) {
                if (a < INT64_MAX / b) throw rational_overflow("int64 overflow in multiplication");
            } else if (a > 0 && b < 0) {
                if (b < INT64_MIN / a) throw rational_overflow("int64 overflow in multiplication");
            } else if (a < 0 && b > 0) {
                if (a < INT64_MIN / b) throw rational_overflow("int64 overflow in multiplication");
            }
            return a * b;
        #endif
    }

    static inline int64_t checked_add(int64_t a, int64_t b) {
        #if defined(__GNUC__) || defined(__clang__)
        int64_t out;
        if (__builtin_add_overflow(a, b, &out))
            throw rational_overflow("int64 overflow in addition");
        return out;
        #else
            // Fallback: manual overflow check
            if ((b > 0 && a > INT64_MAX - b) || (b < 0 && a < INT64_MIN - b))
                throw rational_overflow("int64 overflow in addition");
            return a + b;
        #endif
    }

    static inline int64_t checked_sub(int64_t a, int64_t b) {
        #if defined(__GNUC__) || defined(__clang__)
        int64_t out;
        if (__builtin_sub_overflow(a, b, &out))
            throw rational_overflow("int64 overflow in subtraction");
        return out;
        #else
            // Fallback: manual overflow check
            if ((b > 0 && a < INT64_MIN + b) || (b < 0 && a > INT64_MAX + b))
                throw rational_overflow("int64 overflow in subtraction");
            return a - b;
        #endif
    }

public:
    // -------------------------------------------------------------
    // Constructors
    // -------------------------------------------------------------
    rational64() noexcept : n_(0), d_(1) {}

    rational64(int64_t numerator) noexcept
        : n_(numerator), d_(1) {}

    rational64(int64_t num, int64_t den)
        : n_(num), d_(den) {
        if (den == 0) {
            throw rational_overflow("denominator cannot be zero");
        }
        normalize();
    }

    // -------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------
    inline int64_t numerator() const noexcept { return n_; }
    inline int64_t denominator() const noexcept { return d_; }
    
    // Legacy accessors for backward compatibility
    inline int64_t num() const noexcept { return n_; }
    inline int64_t den() const noexcept { return d_; }

    // -------------------------------------------------------------
    // Arithmetic operators
    // -------------------------------------------------------------
    inline rational64 operator+(const rational64& r) const {
        // (a/b + c/d) = (ad + bc) / bd
        int64_t ad = checked_mul(n_, r.d_);
        int64_t bc = checked_mul(r.n_, d_);
        int64_t num = checked_add(ad, bc);
        int64_t den = checked_mul(d_, r.d_);
        return rational64(num, den);
    }

    inline rational64 operator-(const rational64& r) const {
        int64_t ad = checked_mul(n_, r.d_);
        int64_t bc = checked_mul(r.n_, d_);
        int64_t num = checked_sub(ad, bc);
        int64_t den = checked_mul(d_, r.d_);
        return rational64(num, den);
    }

    inline rational64 operator*(const rational64& r) const {
        // multiply first — gcd-reduce early to avoid overflow
        int64_t g1 = gcd_int64(n_, r.d_);
        int64_t g2 = gcd_int64(r.n_, d_);

        int64_t a = n_ / g1;
        int64_t d1 = r.d_ / g1;
        int64_t c = r.n_ / g2;
        int64_t d2 = d_ / g2;

        int64_t num = checked_mul(a, c);
        int64_t den = checked_mul(d1, d2);

        return rational64(num, den);
    }

    inline rational64 operator/(const rational64& r) const {
        if (r.n_ == 0)
            throw rational_overflow("division by zero");

        int64_t g1 = gcd_int64(n_, r.n_);
        int64_t g2 = gcd_int64(d_, r.d_);

        int64_t a = n_ / g1;
        int64_t d1 = r.n_ / g1;
        int64_t c = r.d_ / g2;
        int64_t d2 = d_ / g2;

        int64_t num = checked_mul(a, c);
        int64_t den = checked_mul(d2, d1);

        return rational64(num, den);
    }

    // -------------------------------------------------------------
    // Comparisons
    // -------------------------------------------------------------
    inline bool operator==(const rational64& r) const noexcept {
        return (n_ == r.n_) && (d_ == r.d_);
    }

    inline bool operator<(const rational64& r) const {
        // compare a/b < c/d  =>  ad < cb
        // Quick check: if signs differ, result is obvious
        bool this_neg = (n_ < 0);
        bool r_neg = (r.n_ < 0);
        if (this_neg != r_neg) {
            return this_neg;  // Negative < positive
        }
        
        // Same sign: compare ad < cb
        // Try 64-bit first if safe (optimization)
        if (!this_neg && n_ > 0 && r.n_ > 0) {
            // Both positive - check for overflow risk
            if (n_ <= INT64_MAX / r.d_ && r.n_ <= INT64_MAX / d_) {
                // Safe to use 64-bit
                int64_t lhs = n_ * r.d_;
                int64_t rhs = r.n_ * d_;
                return lhs < rhs;
            }
        }
        
        // Need extended precision or fallback
        #if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
            // Use 128-bit arithmetic (GCC/Clang)
        __int128 lhs = (__int128)n_ * r.d_;
        __int128 rhs = (__int128)r.n_ * d_;
        return lhs < rhs;
        #else
            // Fallback: use double (may lose precision for very large numbers)
            // This is acceptable for comparison purposes
            return to_double() < r.to_double();
        #endif
    }

    inline bool operator!=(const rational64& r) const noexcept { return !(*this == r); }
    inline bool operator>(const rational64& r) const { return r < *this; }
    inline bool operator<=(const rational64& r) const { return !(*this > r); }
    inline bool operator>=(const rational64& r) const { return !(*this < r); }

    // -------------------------------------------------------------
    // Unary operators
    // -------------------------------------------------------------
    inline rational64 operator+() const noexcept {
        return *this;
    }

    inline rational64 operator-() const {
        if (n_ == INT64_MIN && d_ == 1) {
            throw rational_overflow("cannot negate INT64_MIN");
        }
        return rational64(-n_, d_);
    }

    // -------------------------------------------------------------
    // Compound assignment operators
    // -------------------------------------------------------------
    inline rational64& operator+=(const rational64& r) {
        *this = *this + r;
        return *this;
    }

    inline rational64& operator-=(const rational64& r) {
        *this = *this - r;
        return *this;
    }

    inline rational64& operator*=(const rational64& r) {
        *this = *this * r;
        return *this;
    }

    inline rational64& operator/=(const rational64& r) {
        *this = *this / r;
        return *this;
    }

    // -------------------------------------------------------------
    // Arithmetic operators with integers
    // -------------------------------------------------------------
    inline rational64 operator+(int64_t i) const {
        return *this + rational64(i);
    }

    friend inline rational64 operator+(int64_t i, const rational64& r) {
        return rational64(i) + r;
    }

    inline rational64 operator-(int64_t i) const {
        return *this - rational64(i);
    }

    friend inline rational64 operator-(int64_t i, const rational64& r) {
        return rational64(i) - r;
    }

    inline rational64 operator*(int64_t i) const {
        return *this * rational64(i);
    }

    friend inline rational64 operator*(int64_t i, const rational64& r) {
        return rational64(i) * r;
    }

    inline rational64 operator/(int64_t i) const {
        return *this / rational64(i);
    }

    friend inline rational64 operator/(int64_t i, const rational64& r) {
        return rational64(i) / r;
    }

    inline rational64& operator+=(int64_t i) {
        return *this += rational64(i);
    }

    inline rational64& operator-=(int64_t i) {
        return *this -= rational64(i);
    }

    inline rational64& operator*=(int64_t i) {
        return *this *= rational64(i);
    }

    inline rational64& operator/=(int64_t i) {
        return *this /= rational64(i);
    }

    // -------------------------------------------------------------
    // Comparison operators with integers
    // -------------------------------------------------------------
    inline bool operator==(int64_t i) const noexcept {
        return (d_ == 1) && (n_ == i);
    }

    friend inline bool operator==(int64_t i, const rational64& r) noexcept {
        return r == i;
    }

    inline bool operator!=(int64_t i) const noexcept {
        return !(*this == i);
    }

    friend inline bool operator!=(int64_t i, const rational64& r) noexcept {
        return r != i;
    }

    inline bool operator<(int64_t i) const {
        return *this < rational64(i);
    }

    friend inline bool operator<(int64_t i, const rational64& r) {
        return rational64(i) < r;
    }

    inline bool operator>(int64_t i) const {
        return *this > rational64(i);
    }

    friend inline bool operator>(int64_t i, const rational64& r) {
        return rational64(i) > r;
    }

    inline bool operator<=(int64_t i) const {
        return *this <= rational64(i);
    }

    friend inline bool operator<=(int64_t i, const rational64& r) {
        return rational64(i) <= r;
    }

    inline bool operator>=(int64_t i) const {
        return *this >= rational64(i);
    }

    friend inline bool operator>=(int64_t i, const rational64& r) {
        return rational64(i) >= r;
    }

    // -------------------------------------------------------------
    // Utility functions
    // -------------------------------------------------------------
    inline rational64 abs() const {
        if (n_ < 0) {
            if (n_ == INT64_MIN && d_ == 1) {
                throw rational_overflow("cannot compute abs of INT64_MIN");
            }
            return rational64(-n_, d_);
        }
        return *this;
    }

    inline rational64 inverse() const {
        if (n_ == 0) {
            throw rational_overflow("cannot invert zero");
        }
        return rational64(d_, n_);
    }

    inline bool is_zero() const noexcept {
        return n_ == 0;
    }

    inline bool is_positive() const noexcept {
        return n_ > 0;
    }

    inline bool is_negative() const noexcept {
        return n_ < 0;
    }

    inline bool is_integer() const noexcept {
        return d_ == 1;
    }

    // -------------------------------------------------------------
    // Convert to double
    // -------------------------------------------------------------
    inline double to_double() const noexcept {
        // d_ is guaranteed non-zero by constructor
        return static_cast<double>(n_) / static_cast<double>(d_);
    }

    // -------------------------------------------------------------
    // Static constants
    // -------------------------------------------------------------
    static const rational64 ZERO;
    static const rational64 ONE;
    static const rational64 MINUS_ONE;
};

// -------------------------------------------------------------
// Static constant definitions
// -------------------------------------------------------------
inline const rational64 rational64::ZERO(0, 1);
inline const rational64 rational64::ONE(1, 1);
inline const rational64 rational64::MINUS_ONE(-1, 1);

// -------------------------------------------------------------
// Hash function for std::unordered_map/std::unordered_set
// -------------------------------------------------------------
namespace std {
    template<>
    struct hash<rational64> {
        size_t operator()(const rational64& r) const noexcept {
            // Combine hash of numerator and denominator
            size_t h1 = std::hash<int64_t>{}(r.numerator());
            size_t h2 = std::hash<int64_t>{}(r.denominator());
            // Combine hashes (boost::hash_combine style)
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };
}

// -------------------------------------------------------------
// Stream operators
// -------------------------------------------------------------
inline std::ostream& operator<<(std::ostream& os, const rational64& r) {
    if (r.denominator() == 1) {
        return os << r.numerator();
    } else {
        return os << r.numerator() << "/" << r.denominator();
    }
}

inline std::istream& operator>>(std::istream& is, rational64& r) {
    int64_t num, den = 1;
    char c;
    
    if (is >> num) {
        if (is.get(c) && c == '/') {
            if (!(is >> den)) {
                is.setstate(std::ios::failbit);
                return is;
            }
            if (den == 0) {
                is.setstate(std::ios::failbit);
                return is;
            }
        } else {
            is.putback(c);
        }
        r = rational64(num, den);
    }
    return is;
}

// -------------------------------------------------------------
// String conversion
// -------------------------------------------------------------
inline std::string to_string(const rational64& r) {
    if (r.denominator() == 1) {
        return std::to_string(r.numerator());
    } else {
        return std::to_string(r.numerator()) + "/" + std::to_string(r.denominator());
    }
}
