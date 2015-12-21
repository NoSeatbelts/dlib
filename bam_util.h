#ifndef PAIR_UTIL_H
#define PAIR_UTIL_H
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
/* @func INC_TAG increments an int teag
 * :param: p [bam1_t *] One bam record
 * :param: b [bam1_t *] Second bam record
 * :param: key [const char *] Bam aux key
 */
#define inc_tag(p, b, key, type) *(type *)(bam_aux_get(p, key) + 1) += *(type *)(bam_aux_get(b, key) + 1);
#define inc_tag_int(p, b, key) inc_tag(p, b, key, int)
#define inc_tag_float(p, b, key) inc_tag(p, b, key, float)
#define set_base(pSeq, bSeq, i) (pSeq)[(i)>>1] = ((bam_seqi(bSeq, i) << (((~i) & 1) << 2)) | (((pSeq)[(i)>>1]) & (0xf0U >> (((~i) & 1) << 2))))
#define n_base(pSeq, i) pSeq[(i)>>1] |= (0xf << ((~(i) & 1) << 2));

#define check_fa(arr, fm, len) \
		do {\
		for(int i##arr = 0; i##arr < len; ++i##arr) {\
			if(arr[i##arr] > fm){\
				fprintf(stderr, "[E:%s] %" PRIu32 " arr value greater than FM %" PRIu32 ".\n",\
                        __func__, arr[i##arr], fm);\
				exit(EXIT_FAILURE);\
			}\
		}\
		} while(0)


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
	}
	else {
		b->core.mpos = b->core.pos;
		b->core.mtid = b->core.tid;
		b->core.flag |= BAM_FMUNMAP;
	}
	return;
}

static inline void add_unclipped_mate_starts(bam1_t *b1, bam1_t *b2) {
    uint32_t i, offset1 = 0u, offset2 = 0u;
	int32_t ucs1 = b1->core.pos, ucs2 = b2->core.pos;
	const uint32_t *cigar1 = bam_get_cigar(b1);
	const uint32_t *cigar2 = bam_get_cigar(b2);
	for(i = 0; i < b1->core.n_cigar; ++i) {
		if(!(cigar1[i]&0xf)) // 'M' in cigar.
			break;
		else
			offset1 += cigar1[i] >> BAM_CIGAR_SHIFT;
	}
	for(i = 0; i < b2->core.n_cigar; ++i) {
		if(!(cigar2[i]&0xf))
			break;
		else
			offset2 += cigar2[i] >> BAM_CIGAR_SHIFT;
	}
	if(b1->core.flag & BAM_FREVERSE)
		ucs1 += offset1;
	else
		ucs1 -= offset1;
	if(b2->core.flag & BAM_FREVERSE)
		ucs2 += offset2;
	else
		ucs2 -= offset2;
	bam_aux_append(b2, "MU", 'I', sizeof(uint32_t), (uint8_t *)&ucs1);
	bam_aux_append(b1, "MU", 'I', sizeof(uint32_t), (uint8_t *)&ucs2);
	bam_aux_append(b2, "SU", 'I', sizeof(uint32_t), (uint8_t *)&ucs2);
	bam_aux_append(b1, "SU", 'I', sizeof(uint32_t), (uint8_t *)&ucs1);
	return;
}

#define check_bam_tag(bamrec, tag) \
	do {\
		if(!bam_aux_get(bamrec, tag)) {\
		fprintf(stderr, "[E:%s] Required bam tag '%s' not found. Abort mission!\n", __func__, tag),\
		exit(EXIT_FAILURE);\
		}\
	} while(0)

CONST static inline void *array_tag(bam1_t *b, const char *tag) {
	const uint8_t *data = bam_aux_get(b, tag);
	return data ? (void *)(data + sizeof(int) + 2): NULL;
}

extern void bam_plp_set_maxcnt(bam_plp_t, int);

#endif
