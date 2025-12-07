#ifndef RATIONAL_EIGEN_H
#define RATIONAL_EIGEN_H

#include <fracessa/types.hpp>
#include <Eigen/Core>
#include <limits>
#include <cmath>

#ifdef USE_GMP_RATIONAL
    #include <boost/multiprecision/gmp.hpp>
#else
    #include <boost/multiprecision/cpp_int.hpp>
#endif

namespace Eigen {

// Specialization of NumTraits for rational (GMP or cpp_rational backend)
template<>
struct NumTraits<rational> : GenericNumTraits<rational>
{
    typedef rational Real;
    typedef rational NonInteger;
    typedef rational Nested;
    typedef rational Literal;
    
    enum {
        IsInteger = 0,
        IsSigned = 1,
        IsComplex = 0,
        RequireInitialization = 1,
        ReadCost = 1,
        AddCost = 3,
        MulCost = 3
    };
    
    static inline rational epsilon() {
        // Rational numbers have exact precision, so epsilon is effectively 0
        // Return a very small rational (1/2^100)
        #ifdef USE_GMP_RATIONAL
            // Use GMP mpz_int for large powers
            boost::multiprecision::mpz_int base(2);
            boost::multiprecision::mpz_int power = boost::multiprecision::pow(base, 100);
        #else
            // Use cpp_int for large powers
            boost::multiprecision::cpp_int base(2);
            boost::multiprecision::cpp_int power = boost::multiprecision::pow(base, 100);
        #endif
        return rational(1) / rational(power);
    }
    
    static inline rational dummy_precision() {
        return epsilon();
    }
    
    static inline int digits10() {
        // Rational numbers have arbitrary precision
        // Return a large value to indicate high precision
        return std::numeric_limits<double>::digits10 * 2;
    }
    
    static inline rational highest() {
        // Return a very large rational number
        // Using a large power of 2
        #ifdef USE_GMP_RATIONAL
            // Use GMP mpz_int for large powers
            boost::multiprecision::mpz_int base(2);
            boost::multiprecision::mpz_int power = boost::multiprecision::pow(base, 1000);
        #else
            // Use cpp_int for large powers
            boost::multiprecision::cpp_int base(2);
            boost::multiprecision::cpp_int power = boost::multiprecision::pow(base, 1000);
        #endif
        return rational(power);
    }
    
    static inline rational lowest() {
        return -highest();
    }
    
    static inline rational min() {
        return lowest();
    }
    
    static inline rational max() {
        return highest();
    }
    
    static inline rational infinity() {
        // Rational numbers don't have infinity
        // Return a very large number as approximation
        return highest();
    }
    
    static inline rational quiet_NaN() {
        // Rational numbers don't have NaN
        // Return 0 as a safe default
        return rational(0);
    }
};

} // namespace Eigen

#endif // RATIONAL_EIGEN_H

