#ifndef BAM_UTIL_H
#define BAM_UTIL_H
#include "io_util.h"
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include <unistd.h>
#include "sam_opts.h"
#include "sam.h"
#include "bam.h"
#include "char_util.h"
#include "compiler_util.h"

typedef void (*pair_fn)(bam1_t *b,bam1_t *b1);
typedef void (*single_fn)(bam1_t *b);
typedef void (*single_aux)(bam1_t *b, void *data);
typedef int (*single_aux_check)(bam1_t *b, void *data);
static inline void add_unclipped_mate_starts(bam1_t *b1, bam1_t *b2);
void abstract_pair_iter(samFile *in, bam_hdr_t *hdr, samFile *ofp, pair_fn function);
void abstract_single_filter(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux_check function, void *data);
void abstract_single_data(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux function, void *data);
void abstract_single_iter(samFile *in, bam_hdr_t *hdr, samFile *out, single_fn function);

enum htseq {
	HTS_A = 1,
	HTS_C = 2,
	HTS_G = 4,
	HTS_T = 8,
	HTS_N = 15
};

// bam utility macros.

/* @func inc_tag increments a numeric tag with a given type
 * :param: p [bam1_t *] One bam record
 * :param: b [bam1_t *] Second bam record
 * :param: key [const char *] Bam aux key
 */
#define inc_tag(p, b, key, type) *(type *)(bam_aux_get(p, key) + 1) += *(type *)(bam_aux_get(b, key) + 1);

/* @func inc_tag_int increments an integer tag with another integer tag.
 * :param: p [bam1_t *] One bam record
 * :param: b [bam1_t *] Second bam record
 * :param: key [const char *] Bam aux key
 */

#define inc_tag_int(p, b, key) inc_tag(p, b, key, int)
/* @func inc_tag_float increments a float tag
 * :param: p [bam1_t *] One bam record
 * :param: b [bam1_t *] Second bam record
 * :param: key [const char *] Bam aux key
 */
#define inc_tag_float(p, b, key) inc_tag(p, b, key, float)

/* @func set_base sets the nucleotide at index i in read p to be set to base at index i in read b.
 * :param: p [bam1_t *] One bam record
 * :param: b [bam1_t *] Second bam record
 * :param: i [index] Base position in read
 */
#define set_base(pSeq, bSeq, i) (pSeq)[(i)>>1] = ((bam_seqi(bSeq, i) << (((~i) & 1) << 2)) | (((pSeq)[(i)>>1]) & (0xf0U >> (((~i) & 1) << 2))))

/* @func n_base sets the nucleotide at index i in read p to N.
 * :param: p [bam1_t *] One bam record
 * :param: i [index] Base position in read
 */
#define n_base(pSeq, i) pSeq[(i)>>1] |= (0xf << ((~(i) & 1) << 2));

/* Just an array-checking utility for debugging. I don't see much use for this. */
#define check_fa(arr, fm, len) \
		do {\
		for(int i##arr = 0; i##arr < len; ++i##arr) {\
			if(arr[i##arr] > fm){\
				fprintf(stderr, "[E:%s] %u arr value greater than FM %u.\n", arr[i##arr], fm);\
				exit(EXIT_FAILURE);\
			}\
		}\
		} while(0)

/* @func process_mei_tag
 * This is deprecated. Keeping it around for no real reason.
 */
inline void process_mei_tag(bam1_t *b) {
	uint8_t *const tag_ptr = bam_aux_get(b, "ME");
	if(UNLIKELY(!tag_ptr)) {
		fprintf(stderr, "[E:%s] Expected ME tag not present. Abort mission! Qname: %s.", __func__,
				(char *)bam_get_qname(b));
		exit(EXIT_FAILURE);
	}
	if(bam_aux2i(tag_ptr)) {
		b->core.pos = b->core.mpos;
		b->core.tid = b->core.mtid;
		b->core.flag |= BAM_FUNMAP;
	} else {
		b->core.mpos = b->core.pos;
		b->core.mtid = b->core.tid;
		b->core.flag |= BAM_FMUNMAP;
	}
}

static inline int get_unclipped_start(bam1_t *b)
{
	if(b->core.flag & BAM_FUNMAP)
		return -1;
	int offset = 0;
	uint32_t i;
	const uint32_t *const cigar = bam_get_cigar(b);
	for(i = 0; i < b->core.n_cigar; ++i) {
		if(bam_cigar_op(cigar[i]) == 0) break; // 'M' in cigar
		else offset += bam_cigar_oplen(cigar[i]);
	}
	return b->core.pos + ((b->core.flag & BAM_FREVERSE) ? offset: -offset);
}


/*  @func add_unclipped_mate_starts
 *  @abstract Adds the unclipped start positions for each read and its mate
 */
static inline void add_unclipped_mate_starts(bam1_t *b1, bam1_t *b2) {
	const int32_t ucs1 = get_unclipped_start(b1); const int32_t ucs2 = get_unclipped_start(b2);
	bam_aux_append(b2, "MU", 'i', sizeof(int32_t), (uint8_t *)&ucs1);
	bam_aux_append(b1, "MU", 'i', sizeof(int32_t), (uint8_t *)&ucs2);
	bam_aux_append(b2, "SU", 'i', sizeof(int32_t), (uint8_t *)&ucs2);
	bam_aux_append(b1, "SU", 'i', sizeof(int32_t), (uint8_t *)&ucs1);
}

#define check_bam_tag(bamrec, tag) \
	do {\
		if(!bam_aux_get(bamrec, tag)) {\
		fprintf(stderr, "[E:%s] Required bam tag '%s' not found. Abort mission!\n", __func__, tag),\
		exit(EXIT_FAILURE);\
		}\
	} while(0)

CONST static inline void *array_tag(bam1_t *b, const char *tag) {
	const uint8_t *const data = bam_aux_get(b, tag);
	return data ? (void *)(data + sizeof(int) + 2): NULL;
}

void bam_plp_set_maxcnt(bam_plp_t, int);
#endif // BAM_UTIL_H
