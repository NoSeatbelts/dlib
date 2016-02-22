#ifndef BED_UTIL_H
#define BED_UTIL_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <inttypes.h>
#include "htslib/khash.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "dlib/logging_util.h"
#include "dlib/mem_util.h"
#include "dlib/bam_util.h"

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

#ifdef __cplusplus
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

class ParsedBed;
class RegionSet {
	friend ParsedBed;
	std::vector<uint64_t> intervals;
	std::string contig_name;
	std::vector<std::string> region_names;
	RegionSet(): intervals(std::vector<uint64_t>()),
			contig_name(""){
	}
public:
	RegionSet(int start, int stop, char *refname, char *region_name): intervals(1, to_ivl(start, stop)),
			contig_name(refname ? refname: NO_ID_STR),
			region_names(1, region_name){
	}
	void add_region(int start, int stop, char *region_name) {
		intervals.push_back(to_ivl(start, stop));
		region_names.push_back(region_name);
		assert(region_names.size() == intervals.size());
	}
};
class ParsedBed {
	std::vector<int> sorted_keys;
	std::unordered_map<int, RegionSet> contig_hash;
	int test(int tid, int pos) {
		auto match = contig_hash.find(tid);
		if(match == contig_hash.end())
			return 0;
		for(auto& ivl: match->second.intervals)
			if(pos >= get_start(ivl) && pos < get_stop(ivl)) return 1;
		return 0;
	}
	int bcf1_test(bcf1_t *vrec) {
		return test(vrec->rid, vrec->pos);
	}
public:
	int bam1_test(bam1_t *b) {
		if(b->core.flag & BAM_FUNMAP) return 0;
		auto match = contig_hash.find(b->core.tid);
		if(match == contig_hash.end())
			return 0;
		for(auto& ivl: match->second.intervals)
			if(bam_getend(b) >= get_start(ivl) && b->core.pos < get_stop(ivl)) return 1;
		return 0;
	}
	ParsedBed(const char *path, bam_hdr_t *header, uint32_t padding=DEFAULT_PADDING) {
		contig_hash = std::unordered_map<int, RegionSet>();
		FILE *ifp = fopen(path, "r");
		char *line = NULL;
		char *tok = NULL;
		size_t len = 0;
		ssize_t read;
		int tid;
		int64_t start, stop;
		std::unordered_map<int, RegionSet>::iterator it;
		std::unordered_set<int> keyset;
		while ((read = getline(&line, &len, ifp)) != -1) {
			switch(*line) {
				case '\0': case '#': continue;
			}
			tok = strtok(line, "\t");
			tid = bam_name2id(header, tok);
			keyset.insert(tid);
			tok = strtok(NULL, "\t");
			start = strtoll(tok, NULL, 10) - padding;
			tok = strtok(NULL, "\t");
			stop = strtoll(tok, NULL, 10) + padding;
			tok = strtok(NULL, "\t");
			if((it = contig_hash.find(tid)) == contig_hash.end())
				it->second = RegionSet(start, stop, header->target_name[tid], tok ? tok: NO_ID_STR);
			else it->second.add_region(start, stop, tok ? tok: NO_ID_STR);
		}
		sorted_keys = std::vector<int>(keyset.begin(), keyset.end());
		keyset.clear();
		std::sort(sorted_keys.begin(), sorted_keys.end());
		fclose(ifp);
	}
};
#endif

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
khash_t(bed) *parse_bed_hash(const char *path, bam_hdr_t *header, uint32_t padding);
khash_t(bed) *build_ref_hash(bam_hdr_t *header);
void *bed_read(const char *fn);
void bed_destroy_hash(void *);
size_t get_nregions(khash_t(bed) *h);

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
