#ifndef TYPES_H
#define TYPES_H

//*************************** boost ***********************************************************
// Choose backend based on CMake option USE_GMP_RATIONAL
#ifdef USE_GMP_RATIONAL
    // GMP backend: faster, requires libgmp-dev
    #include <boost/multiprecision/gmp.hpp>
    typedef boost::multiprecision::mpq_rational rational;
#else
    // Pure C++ backend: no dependencies, but slower
    #include <boost/multiprecision/cpp_int.hpp>
    typedef boost::multiprecision::cpp_rational rational;
#endif
//************************************************************************************************

#endif // TYPES_H

