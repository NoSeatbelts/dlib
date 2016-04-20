#ifndef MATH_UTIL_H
#define MATH_UTIL_H
#include <stdint.h>
#include <stdlib.h>

namespace dlib {
    double hellinger(double* arr1, double* arr2, size_t length);
    static inline int64_t ipow(int base, int exp)
    {
        int64_t ret = 1;
        while(exp) {
            if(exp & 1) ret *= base;
            exp >>= 1;
            base *= base;
        }
        return ret;
    }
}



#endif /* MATH_UTIL_H */
