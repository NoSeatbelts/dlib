#include "kmer_util.h"

namespace dlib {
//Set a return kmer (overlap_kmer) for the matched sequence, and len_pos, which contains the length of the kmer
//and the index of the kmer match in ka as the first 16 and last 16 bits of a uint32_t, respectively.

unsigned cstr_longest_overlap(char *s1, char *s2, unsigned n) {
    // There has to be a better than cubic algorithm for this.
    unsigned max_runlen = 0u, runlen, max_index = 0u;
    for(unsigned i = 0u; i < n; ++i) {
        for(unsigned j = 0u; j < n; ++j) {
            if(n - i < max_runlen) break;
            runlen = 0u;
            for(unsigned k = 0u; k + j < n; ++k) {
                if(s1[k + i] != s2[k + j]) break;
                ++runlen;
            }
            if(runlen > max_runlen) max_runlen = runlen, max_index = i;
        }
    }
    LOG_DEBUG("max_index: %u.\n", max_index);
    return max_runlen;
}

inline uint32_t get_runlen(uint64_t ka, uint64_t kb, int i, int j, const int k) {
    uint32_t ret(0);
    while(encoded_base(ka, k, i++) == encoded_base(kb, k, j++) && i < k && j < k) ++ret;
    fprintf(stderr, "runlen: %u.\n", ret);
    return ret;
}

void kmer_overlap(uint64_t ka, uint64_t kb, const uint32_t k,
                  uint64_t *overlap_kmer, uint32_t *len_pos) {
#if !NDEBUG
    char kmer1[32], kmer2[32];
    kmer2cstr(ka, 31u, kmer1);
    kmer2cstr(kb, 31u, kmer2);
    fprintf(stderr, "Maximum run should be %i.\n", cstr_longest_overlap(kmer1, kmer2, k));
#endif
    uint_fast8_t max_runlen = 0, max_runindex = 0, runlen;
    uint_fast8_t i, j, m;
    for(i = 0; i < k; ++i) {
        for(j = 0; j < k; ++j) {
            runlen = 0;
            for(m = 0; m + j < k && m + i < k; ++m) {
                LOG_DEBUG("encoded bases: %i, %i. kmer string bases: %c, %c.\n", encoded_base(ka, k, i + m), encoded_base(kb, k, j + m), kmer1[i + m], kmer2[j + m]);
                if(encoded_base(ka, k, i + m) != encoded_base(kb, k, j + m)) break;
                if(kmer1[i + m] != kmer2[j + m]) {
                    LOG_DEBUG("at indices %i, %i. encoded bases: %i, %i. kmer string bases: %c, %c.\n", i + m, j + m, encoded_base(ka, k, i + m), encoded_base(kb, k, j + m), kmer1[i + m], kmer2[j + m]);
                    assert(false);
                }
                ++runlen;
            }
            if(runlen > max_runlen) max_runlen = runlen, max_runindex = i;
        }
        // Should we be short-circuing?
        if(k - i < max_runlen) break;
    }
    // Set return values.
    if(max_runlen < MIN_RUNLEN) {
        *overlap_kmer = BINFINITY;
        *len_pos = 0;
    } else {
        LOG_DEBUG("len: %u. index: %u.\n", max_runlen, max_runindex);
        *overlap_kmer = (ka >> (2 * max_runindex)) & (BF >> (2 * (32 - max_runlen)));
        *len_pos = ((uint32_t)max_runlen << 16) | max_runindex;
    }
}

} // namespace dlib
