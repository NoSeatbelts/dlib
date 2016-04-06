#ifndef BED_UTIL_H
#define BED_UTIL_H

#ifndef __STDC_LIMIT_MACROS
#    define __STDC_LIMIT_MACROS
#endif

#include <stddef.h>
#include <stdio.h>
#include <inttypes.h>
#ifdef __cplusplus
extern "C" {
#endif
    #include <stdint.h>
    #include "htslib/khash.h"
    #include "htslib/sam.h"
    #include "htslib/vcf.h"
    #include "htslib/faidx.h"
#ifdef __cplusplus
}
#endif
#include "dlib/logging_util.h"
#include "dlib/misc_util.h"
#include "dlib/io_util.h"
#ifdef __cplusplus
#include <vector>
#include <algorithm>
#endif

#define DEFAULT_PADDING 0u
#define NO_ID_STR ((char *)"MissingContigName")


// Bed interval query utility macros.

/* @func get_start
 * @abstract Returns the start of a bed interval in a uint64_t
 * :param: ivl [uint64_t] Interval to extract start from.
 * :returns: [uint32_t] (start)
 */
#define get_start(ivl) (int32_t)((ivl) >> 32)

/* @func get_stop
 * @abstract Returns the stop of a bed interval in a uint64_t
 * :param: ivl [uint64_t] Interval to extract start from.
 * :returns: [uint32_t] (stop)
 */
#define get_stop(ivl) ((int32_t)(ivl))

/* @func to_ivl
 * @abstract Converts a start and stop combination and converts it into an interval.
 * :param: start [int32_t/uint32_t] start position of the interval
 * :param: stop [int32_t/uint32_t] stop position of the interval
 * :returns: [uint64_t] Interval encoded in start/stop format.
 */
#define to_ivl(start, stop) (start > 0 ? (((uint64_t)start) << 32 | (stop)): stop)

// Like bam_endpos, but doesn't check that the read is mapped, as that's already been checked.
#ifndef bam_getend
#define bam_getend(b) ((b)->core.pos + bam_cigar2rlen((b)->core.n_cigar, bam_get_cigar(b)))
#endif


/*
 * struct region_set, aka region_set_t
 * Struct holding an array of intervals and a count for the number of intervals.
 * Currently not expecting tons of bed intervals, so sizing is handled directly,
 * not with a used/max protocol like a vector or khash.
 * Essentially a vector.
 */
typedef struct region_set {
    uint64_t *intervals;
    char *contig_name;
    uint64_t n;
} region_set_t;
/*
 * khash_t(bed) is now a bed file: a region_set_t as a value for the key which is the contig
 * of the interval.
 */
KHASH_MAP_INIT_INT(bed, region_set_t)


#ifdef __cplusplus
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
namespace dlib {
#endif
    void sort_bed_hash(khash_t(bed) *hash);
    khash_t(bed) *parse_bed_hash(const char *path, bam_hdr_t *header, uint32_t padding);
#ifdef __cplusplus
    // Overloading
    khash_t(bed) *parse_bed_hash(const char *path, faidx_t *fai, uint32_t padding);
#endif
    static int intcmp(const void *a, const void *b); // Compare intervals for sorting by start
    void sort_bed(khash_t(bed) *bed);
    khash_t(bed) *build_ref_hash(bam_hdr_t *header);
    void *bed_read(const char *fn);
    void bed_destroy_hash(void *);
    size_t get_nregions(khash_t(bed) *h);

    static inline int bed_test(bam1_t *b, khash_t(bed) *h)
    {
        khint_t k;
        if((b->core.flag & BAM_FUNMAP) == 0)
            if((k = kh_get(bed, h, b->core.tid)) != kh_end(h))
                for(uint64_t i = 0; i < kh_val(h, k).n; ++i)
                    if(get_start(kh_val(h, k).intervals[i]) <= bam_getend(b))
                        if(b->core.pos <= get_stop(kh_val(h, k).intervals[i]))
                            return 1;
        return 0;
        // Returns 0 if b is unmapped or k == kh_end(h)
    }

#ifdef __cplusplus
    std::vector<khiter_t> make_sorted_keys(khash_t(bed) *h);
    class ParsedBed {
        khash_t(bed) *contig_hash;
        std::vector<khiter_t> sorted_keys;
    public:
        int test(int tid, int pos) {
            khiter_t k;
            if((k = kh_get(bed, contig_hash, tid)) != kh_end(contig_hash))
                for(uint64_t i = 0; i < kh_val(contig_hash, k).n; ++i)
                    if(get_start(kh_val(contig_hash, k).intervals[i]) <= pos)
                        if(pos < get_stop(kh_val(contig_hash, k).intervals[i]))
                            return 1;
            return 0;
        }
        int bcf1_test(bcf1_t *vrec) {
            return test(vrec->rid, vrec->pos);
        }
        int bam1_test(bam1_t *b) {
            return bed_test(b, contig_hash);
        }
        ~ParsedBed() {
            bed_destroy_hash(contig_hash);
        }
        ParsedBed(const char *path, faidx_t *fai, uint32_t padding=DEFAULT_PADDING) :
            contig_hash(parse_bed_hash(path, fai, padding)),
            sorted_keys(make_sorted_keys(contig_hash))
        {
            sort_bed_hash(contig_hash);
        }
        ParsedBed(const char *path, bam_hdr_t *header, uint32_t padding=DEFAULT_PADDING) :
            contig_hash(parse_bed_hash(path, header, padding)),
            sorted_keys(make_sorted_keys(contig_hash))
        {
            sort_bed_hash(contig_hash);
        }
    };
    #endif



    static inline int vcf_bed_test(bcf1_t *b, khash_t(bed) *h)
    {
        khint_t k;
        if((k = kh_get(bed, h, b->rid)) != kh_end(h))
            for(uint64_t i = 0; i < kh_val(h, k).n; ++i)
                if(b->pos >= get_start(kh_val(h, k).intervals[i]))
                    if(b->pos < get_stop(kh_val(h, k).intervals[i]))
                        return 1;
        return 0;
    }

#ifdef __cplusplus

} /* namespace dlib */
#else
static int intcmp(const void *a, const void *b) {
    return get_start(*((uint64_t *)a)) == get_start(*((uint64_t *)b)) ? get_stop(*((uint64_t *)a)) < get_stop(*((uint64_t *)a))
                                        : get_start(*((uint64_t *)a)) < get_start(*((uint64_t *)a));
}
#endif
#endif /* BED_UTIL_H */
