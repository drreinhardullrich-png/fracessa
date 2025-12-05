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



#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif


enum class ReasonEss : int
{
    true_pure_ess = 1,
    true_posdef_double = 2,
    true_posdef_rational = 3,
    true_copositive = 4,
    false_not_posdef_and_kay_0_1 = 5,
    false_not_partial_copositive = 6,
    false_not_copositive = 7
};


#ifdef _MSC_VER
#  include <intrin.h>
#  define  __builtin_popcount  __popcnt64
#endif

static inline size_t support_size_from_int(uint64_t support)
{
	return __builtin_popcount((unsigned long long)support);

	//size_t m = 0;
	//for (size_t i = 0; i < 60; i++) {
	//	if ((support & (1ull << i)) != 0)
	//		m++;
	//}
	//return m;
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
    //see: http://blog.regehr.org/archives/1063 with n=1 and 32=dimension
    //but then the left shift has to be corrected, shift 11111... to the right for 64-dimension, then &
    return ((x>>1) | (x << (dimension-1))) & (18446744073709551615ull>>(64-dimension));

}

static inline uint64_t smallest_representation(uint64_t bitsetm, size_t dimension)
{
    uint64_t min_value=18446744073709551615ull;
    for (size_t i=0; i<dimension; i++) {
        bitsetm = shift_right(bitsetm,dimension);
        if (bitsetm < min_value)
            min_value = bitsetm;
    }
    return min_value;
}


#endif // HELPER_H
