#ifndef RATIONAL_H
#define RATIONAL_H

#include <cstdint>
#include <boost/rational.hpp>
#include <boost/safe_numerics/safe_integer.hpp>

//*************************** boost ***********************************************************
// Use Boost multiprecision cpp_rational (pure C++, no external dependencies)
#include <boost/multiprecision/cpp_int.hpp>
typedef boost::multiprecision::cpp_rational rational;
//************************************************************************************************

// Small rational type using safe int64 (faster for small values, limited range)
typedef boost::rational<boost::safe_numerics::safe<int64_t>> small_rational;

// Constants for zero-overhead performance
const rational ZERO = rational(0);
const rational ONE = rational(1);
const rational MINUS_ONE = rational(-1);

#endif // RATIONAL_H

