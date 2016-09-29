#ifndef KMER_UTIL_H
#define KMER_UTIL_H
#include "logging_util.h"
#include "bam_util.h"
// Largest odd kmer that can be held in a uint64_t
#define MAX_KMER 31
#define DEFAULT_KMER 21

#ifndef num2nuc
# ifndef NUM2NUC_STR
#  define NUM2NUC_STR "ACGTN"
# endif
# define num2nuc(x) NUM2NUC_STR[(uint8_t)x]
#endif

#ifndef BINFINITY
#define BINFINITY -1ull
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

    static inline int bam_is_lt(uint8_t *seq, int cpos, int8_t k) {
        int _cpos = cpos + k - 1;
        uint8_t ki1, ki2;
        for(int i = 0; i < k; ++i, --_cpos) {
            ki1 = bam_seqi(seq, i + cpos);
            ki2 = bam_seqi_cmpl(seq, _cpos);
            if(ki1 != ki2)
                return ki1 < ki2;
        }
        return 0;
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
    inline void bam_set_kmer(uint64_t &ret, uint8_t *seq, int cpos, int is_lt,
                             uint64_t kmer_mask, int k) {
        ret = 0;
        if(is_lt) {
            for(int i = cpos; i < cpos + k; ++i) {
                ret <<= 2;
                switch(bam_seqi(seq, i)) {
                    case htseq::HTS_A:
                        break;
                    case htseq::HTS_C:
                        ret |= 1; break;
                    case htseq::HTS_G:
                        ret |= 2; break;
                    case htseq::HTS_T:
                        ret |= 3; break;
                    default:
                        ret = BINFINITY; return;
                }
                ret &= kmer_mask;
            }
        } else {
            for(int i = cpos + k - 1; i >= cpos;--i) {
                ret <<= 2;
                switch (bam_seqi(seq, i)) {
                    case htseq::HTS_A:
                        ret |= 3; break;
                    case htseq::HTS_C:
                        ret |= 2; break;
                    case htseq::HTS_G:
                        ret |= 1; break;
                    case htseq::HTS_T:
                        break;
                    default:
                        ret = BINFINITY; return;
                }
                ret &= kmer_mask;
            }
        }
    } /*bam_set_kmer*/

    inline void cstr_set_kmer(uint64_t &ret, char *seq, int cpos, int is_lt,
                              uint64_t kmer_mask, int k) {
        ret = 0;
        //LOG_DEBUG("Seq: %s at pointer %p.", seq, (void *)seq);
        if(is_lt) {
            seq += cpos;
            for(;k; --k) {
                ret <<= 2;
                switch (*seq++) {
                  case 'A': case 'a':
                    break;
                  case 'C': case 'c':
                    ret |= 1; break;
                  case 'G': case 'g':
                    ret |= 2; break;
                  case 'T': case 't':
                    ret |= 3; break;
                  default:
                    ret = BINFINITY; return;
                }
                ret &= kmer_mask;
            }
        } else {
            seq += cpos + k;
            for(;k; --k) {
                ret <<= 2;
                switch (*--seq) {
                  case 'A': case 'a':
                    ret |= 3; break;
                  case 'C': case 'c':
                    ret |= 2; break;
                  case 'G': case 'g':
                    ret |= 1; break;
                  case 'T': case 't':
                    break;
                  default:
                    ret = BINFINITY; return;
                }
                ret &= kmer_mask;
            }
        }
    } /*cstr_set_kmer*/
#endif /* ifdef __cplusplus */

#ifdef __cplusplus
}
#endif

#endif
