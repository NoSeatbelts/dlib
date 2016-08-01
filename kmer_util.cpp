#include "kmer_util.h"

namespace dlib {
//Set a return kmer (overlap_kmer) for the matched sequence, and len_pos, which contains the length of the kmer
//and the index of the kmer match in ka as the first 16 and last 16 bits of a uint32_t, respectively.
void kmer_overlap(uint64_t ka, uint64_t kb, const uint32_t k,
                  uint64_t *overlap_kmer, uint32_t *len_pos) {
    uint16_t runlen(0), max_runlen(0), max_runindex(0);
    for(uint_fast8_t i(0); i < k; ++i) {
        for(uint_fast8_t j(0); j < k; ++j) {
            if(encoded_base(ka, k, i) == encoded_base(kb, k, j)) ++runlen;
            else {
                if(runlen > max_runlen) {
                    max_runlen = runlen;
                    max_runindex = i;
                    if(k - j < max_runlen) break;
                    // Should we be short-circuing?
                }
                runlen = 0;
            }
        }
        // Should we be short-circuing?
        if(k - i < max_runlen) break;
    }
    // Set return values.
    if(max_runlen < MIN_RUNLEN) {
        *overlap_kmer = BINFINITY;
        *len_pos = 0;
    } else {
        *overlap_kmer = (ka >> (2 * max_runindex)) & (BF >> (2 * (32 - max_runlen)));
        *len_pos = ((uint32_t)runlen << 16) | max_runindex;
    }
}

} // namespace dlib
