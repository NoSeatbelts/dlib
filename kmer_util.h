#ifndef KMER_UTIL_H
#define KMER_UTIL_H
#include <assert.h>
#include "logging_util.h"
// Largest odd kmer that can be held in a uint64_t
#define MAX_KMER 31
#define DEFAULT_KMER 21

#ifndef num2nuc
# ifndef NUM2NUC_STR
#  define NUM2NUC_STR "ACGTN"
# endif
# define num2nuc(x) NUM2NUC_STR[(uint8_t)x]
#endif

#ifdef __cplusplus
namespace dlib {
#endif

    inline void kmer2cstr(uint64_t kmer, int k, char *buf)
    {
        buf[k] = '\0';
        while(k) *buf++ = num2nuc((kmer >> (2 * --k)) & 0x3u);
        //LOG_DEBUG("kmer %lu has now become string '%s'.\n", kmer, start);
    }

#ifdef __cplusplus
}
#endif

#endif
