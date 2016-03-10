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

    // Used to determine the direction in which to encode a kmer
    inline int cstr_rc_lt(char *seq, int k, int cpos) {
        char *_seq1 = cpos + seq, *_seq2 = _seq1 + k - 1;
        for(;k;--k) {
            if(*_seq1 != nuc_cmpl(*_seq2)) return *_seq1 < nuc_cmpl(*_seq2);
            ++_seq1, --_seq2;
        }
        return 0; // This is reverse-complementarily palindromic. Doesn't matter: it's the same string.
    }


#ifdef __cplusplus
}
#endif

#endif
