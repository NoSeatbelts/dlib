#ifndef SORT_UTIL_H
#define SORT_UTIL_H

#include "bam_util.h"

// bam flag macros
#define bam_is_r1(b) (!!((b)->core.flag&BAM_FREAD1))
#define bam_is_r2(b) (!!((b)->core.flag&BAM_FREAD2))
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE != 0)
#define IS_MATE_REVERSE(bam) (((bam)->core.flag&BAM_FMREVERSE) != 0)
#define IS_READ2(bam) (((bam)->core.flag&BAM_FREAD2) != 0)
#define IS_READ1(bam) (((bam)->core.flag&BAM_FREAD1) != 0)


// Unclipped start sort key macros
#define ucs_sort_mate_key(a) ((uint64_t)(a->core.mtid + 1) << 32 | bam_itag(a, "MU") << 1 | bam_is_mrev(a))
#define ucs_sort_core_key(a) (((uint64_t)(a->core.tid + 1) << 32) | ((bam_itag(a, "SU")+1) << 2) | (bam_is_rev(a) << 1) |bam_is_r1(a))

#define ucs_se_sort_key(a) (uint64_t)((uint64_t)a->core.tid<<32|(bam_aux2i(bam_aux_get(b, "SU"))+1)<<2|bam_is_rev(a))

// Pos sort key macros
#define bmfsort_core_key(a) (uint64_t)((uint64_t)a->core.tid<<32|(a->core.pos+1)<<2|bam_is_rev(a)<<1|bam_is_r1(a))
#define bmfsort_mate_key(a) (uint64_t)((uint64_t)a->core.mtid<<32|(a->core.mpos+1)<<1|bam_is_mrev(a))

#define bmfsort_se_key(a) (uint64_t)((uint64_t)a->core.tid<<32|(bam_aux2i(bam_aux_get(b, "SU"))+1)<<2|bam_is_rev(a))


#endif /* SORT_UTIL_H */
