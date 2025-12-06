#ifndef HELPER_H
#define HELPER_H

#include <sstream>
#include <string>
#include <vector>
#include <bitset>

//*************************** boost ***********************************************************
#include <boost/math/special_functions/binomial.hpp>
#include <boost/integer/common_factor.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/algorithm/string.hpp>
//#include "boost/multiprecision/gmp.hpp"

typedef boost::multiprecision::cpp_rational rational;
//typedef boost::multiprecision::gmp_rational rational;
//typedef mpq_rational rational;
//typedef number< rational_adaptor< cpp_int_backend< 255, 256, signed_magnitude, checked, void> > > rational;
//************************************************************************************************

// Named constants to replace magic numbers
static constexpr size_t CANDIDATE_RESERVE_MULTIPLIER = 10;
static constexpr uint64_t UINT64_MAX_VALUE = 18446744073709551615ULL;


#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif




#ifdef _MSC_VER
#  include <intrin.h>
#  define  __builtin_popcount  __popcnt64
#endif

static inline size_t support_size_from_int(uint64_t support)
{
    return __builtin_popcount((unsigned long long)support);
}

static inline size_t postion_of_lowest_setbit(uint64_t bitsetm) //zero based!!!!
{
    size_t m = 0;
    while (true) {
        if ((bitsetm & (1ull << m)) != 0)
            break;
        m++;
    }
    return m;
}

static inline uint64_t shift_right(uint64_t x, size_t dimension)
{
    return ((x>>1) | (x << (dimension-1))) & (UINT64_MAX_VALUE>>(64-dimension));

}

static inline uint64_t smallest_representation(uint64_t bitsetm, size_t dimension)
{
    uint64_t min_value = UINT64_MAX_VALUE;
    for (size_t i=0; i<dimension; i++) {
        bitsetm = shift_right(bitsetm,dimension);
        if (bitsetm < min_value)
            min_value = bitsetm;
    }
    return min_value;
}


#endif // HELPER_H
