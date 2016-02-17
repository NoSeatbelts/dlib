#ifndef BED_UTIL_H
#define BED_UTIL_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>
#include "htslib/khash.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "logging_util.h"
#include "mem_util.h"

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
#define to_ivl(start, stop) (start > 0 ? ((start) << 32 | (stop)): stop)

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
extern "C" {
#endif
int intcmp(const void *a, const void *b); // Compare intervals for sorting by start
void sort_bed(khash_t(bed) *bed);
khash_t(bed) *parse_bed_hash(char *path, bam_hdr_t *header, uint32_t padding);
khash_t(bed) *build_ref_hash(bam_hdr_t *header);
void *bed_read(const char *fn);
void bed_destroy_hash(void *);
size_t get_nregions(khash_t(bed) *h);

// Like bam_endpos, but doesn't check that the read is mapped, as that's already been checked.
#define bam_getend(b) ((b)->core.pos + bam_cigar2rlen((b)->core.n_cigar, bam_get_cigar(b)))

static inline int bed_test(bam1_t *b, khash_t(bed) *h)
{
	khint_t k;
	if((k = kh_get(bed, h, b->core.tid)) == kh_end(h)) return 0;
	for(uint64_t i = 0; i < kh_val(h, k).n; ++i) {
		if(get_start(kh_val(h, k).intervals[i]) <= bam_getend(b) &&
				b->core.pos <= get_stop(kh_val(h, k).intervals[i])) {
			return 1;
		}
	}
	return 0;
}

static inline int vcf_bed_test(bcf1_t *b, khash_t(bed) *h)
{
	khint_t k;
	if((k = kh_get(bed, h, b->rid)) == kh_end(h))
		return 0;
	for(uint64_t i = 0; i < kh_val(h, k).n; ++i) {
		if(b->pos >= get_start(kh_val(h, k).intervals[i]) && b->pos <= get_stop(kh_val(h, k).intervals[i]))
			return 1;
	}
	return 0;
}
#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
#	ifdef __GNUC__
#		include <parallel/algorithm>
#	else
#		include <algorithm>
#	endif
#	include <vector>
std::vector<khiter_t> make_sorted_keys(khash_t(bed) *h);
#endif


#endif /* BED_UTIL_H */
