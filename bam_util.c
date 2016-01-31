#include "bam_util.h"

bam_plp_t bam_plp_maxcnt_init(bam_plp_auto_f func, void *data, int maxcnt)
{
	bam_plp_t iter = bam_plp_init(func, data);
	bam_plp_set_maxcnt(iter, maxcnt);
	return iter;
}

void abstract_single_data(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux function, void *data)
{
	bam1_t *b = bam_init1();
	while (LIKELY(sam_read1(in, hdr, b) >= 0))
		function(b, data), sam_write1(out, hdr, b);
	bam_destroy1(b);
}

void abstract_single_iter(samFile *in, bam_hdr_t *hdr, samFile *out, single_fn function)
{
	bam1_t *b = bam_init1();
	while (LIKELY(sam_read1(in, hdr, b) >= 0))
		function(b), sam_write1(out, hdr, b);
	bam_destroy1(b);
}

void abstract_single_filter(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux_check function, void *data)
{
	bam1_t *b;
	b = bam_init1();
	while (LIKELY(sam_read1(in, hdr, b) >= 0)) {
		if(function(b, data))
			continue;
		sam_write1(out, hdr, b);
	}
	bam_destroy1(b);
}

#ifdef __cplusplus

void abstract_pair_set(samFile *in, bam_hdr_t *hdr, samFile *ofp, std::set<pair_fn> functions)
{
	bam1_t *b = bam_init1(), *b1 = bam_init1();
	while (LIKELY(sam_read1(in, hdr, b) >= 0)) {
		if(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
		if(b->core.flag & BAM_FREAD1) {
			bam_copy1(b1, b); continue;
		}
		for(auto f: functions) f(b1, b);
		sam_write1(ofp, hdr, b), sam_write1(ofp, hdr, b1);
	}
	bam_destroy1(b), bam_destroy1(b1);
}

#endif


void abstract_pair_iter(samFile *in, bam_hdr_t *hdr, samFile *ofp, pair_fn function)
{
	bam1_t *b = bam_init1(), *b1 = bam_init1();
	while (LIKELY(sam_read1(in, hdr, b) >= 0)) {
		if(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
			continue;
		if(b->core.flag & BAM_FREAD1) {
			bam_copy1(b1, b); continue;
		}
		function(b1, b);
		sam_write1(ofp, hdr, b), sam_write1(ofp, hdr, b1);
	}
	bam_destroy1(b), bam_destroy1(b1);
}

int bampath_has_tag(char *bampath, const char *tag)
{
	samFile *fp = sam_open(bampath, "r");
	bam_hdr_t *header = sam_hdr_read(fp);
	if(!header || !fp) {
		LOG_ERROR("Could not open bam file at '%s' for reading. Abort!\n", bampath);
	}
	bam1_t *b = bam_init1();
	if(sam_read1(fp, header, b) < 0) {
		LOG_ERROR("Empty bam file at '%s'. Abort!\n", bampath);
	}
	int ret = !!bam_aux_get(b, tag);
	bam_destroy1(b);
	bam_hdr_destroy(header);
	sam_close(fp);
	return ret;
}

void check_bam_tag_exit(char *bampath, const char *tag)
{
	if(!(strcmp(bampath, "-") && strcmp(bampath, "stdin"))) {
			LOG_WARNING("Could not check for bam tag without exhausting a pipe. "
						"Tag '%s' has not been verified.\n", tag);
			return;
	}
	if(!bampath_has_tag(bampath, tag)) {
		LOG_ERROR("Required bam tag '%s' missing from bam file at path '%s'. Abort!\n", tag, bampath);
	}
}

