#ifndef TYPES_RATIONAL_H
#define TYPES_RATIONAL_H

#include <cstdint>
#include <type_traits>
#include <boost/rational.hpp>
#include <boost/safe_numerics/safe_integer.hpp>
#include <boost/safe_numerics/exception_policies.hpp>
#include <string>
#include <sstream>

//*************************** boost ***********************************************************
// Use Boost multiprecision cpp_rational (pure C++, no external dependencies)
#include <boost/multiprecision/cpp_int.hpp>
typedef boost::multiprecision::cpp_rational rational;
//************************************************************************************************

//*************************** GMP ***********************************************************
// Use GMP mpq_class for arbitrary precision rational numbers
// COMMENTED OUT - now using boost::multiprecision::cpp_rational
// #include <gmpxx.h>
// typedef mpq_class rational;
//************************************************************************************************

// Small rational type using safe int64 (faster for small values, limited range)
// Explicitly define safe int64 with throw_exception policy to guarantee exceptions are thrown on overflow
// The exception_policy template takes 4 parameters:
// 1. Arithmetic error policy (throw_exception)
// 2. Undefined behavior policy (throw_exception)
// 3. Implementation-defined behavior policy (throw_exception)
// 4. Uninitialized value policy (ignore_exception - less critical, can be ignored)
using safe_int64_t = boost::safe_numerics::safe<
    int64_t,
    boost::safe_numerics::native,
    boost::safe_numerics::exception_policy<
        boost::safe_numerics::throw_exception,
        boost::safe_numerics::throw_exception,
        boost::safe_numerics::throw_exception,
        boost::safe_numerics::ignore_exception
    >
>;
typedef boost::rational<safe_int64_t> small_rational;

// // Constants for zero-overhead performance
// const rational ZERO = rational(0);
// const rational ONE = rational(1);
// const rational MINUS_ONE = rational(-1);

// Overload to_string for boost::rational types
// Format: if denominator == 1 then output numerator, else numerator/denominator
template<typename IntType>
inline std::string to_string(const boost::rational<IntType>& r) {
    IntType num = r.numerator();
    IntType den = r.denominator();
    
    if (den == IntType(1)) {
        std::ostringstream oss;
        oss << num;
        return oss.str();
    } else {
        std::ostringstream oss;
        oss << num << "/" << den;
        return oss.str();
    }
}

// Overload to_string for boost::multiprecision::cpp_rational
// Format: if denominator == 1 then output numerator, else numerator/denominator
inline std::string to_string(const boost::multiprecision::cpp_rational& r) {
    boost::multiprecision::cpp_int num = boost::multiprecision::numerator(r);
    boost::multiprecision::cpp_int den = boost::multiprecision::denominator(r);
    
    if (den == 1) {
        return num.str();
    } else {
        std::ostringstream oss;
        oss << num << "/" << den;
        return oss.str();
    }
}

// Convert boost::rational<T> to double
// Handles safe<int64_t> by extracting to int64_t first
template<class T>
double rational_to_double(const boost::rational<T>& r) {
    // For safe<int64_t>, extract to int64_t first, then to double
    if constexpr (std::is_same_v<T, safe_int64_t>) {
        int64_t num = static_cast<int64_t>(r.numerator());
        int64_t den = static_cast<int64_t>(r.denominator());
        return static_cast<double>(num) / static_cast<double>(den);
    } else {
        // For other integer types, direct conversion
        return static_cast<double>(r.numerator()) /
               static_cast<double>(r.denominator());
    }
}

// Overload for boost::multiprecision::cpp_rational
inline double rational_to_double(const boost::multiprecision::cpp_rational& r) {
    return static_cast<double>(r);
}

// Convert small_rational to rational
inline rational small_to_rational(const small_rational& val) {
    // Extract underlying int64_t values from safe<int64_t>
    int64_t num = static_cast<int64_t>(val.numerator());
    int64_t den = static_cast<int64_t>(val.denominator());
    // Construct rational from numerator and denominator
    return rational(num) / rational(den);
}

// Overload to_string for GMP mpq_class
// Format: if denominator == 1 then output numerator, else numerator/denominator
// COMMENTED OUT - now using boost::multiprecision::cpp_rational
// inline std::string to_string(const mpq_class& r) {
//     mpz_class num = r.get_num();
//     mpz_class den = r.get_den();
//     
//     if (den == 1) {
//         return num.get_str();
//     } else {
//         std::ostringstream oss;
//         oss << num.get_str() << "/" << den.get_str();
//         return oss.str();
//     }
// }

#endif // TYPES_RATIONAL_H

