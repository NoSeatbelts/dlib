#include "bed_util.h"


khash_t(bed) *parse_bed_hash(char *path, bam_hdr_t *header, uint32_t padding)
{
	khash_t(bed) *ret = kh_init(bed);
	FILE *ifp = fopen(path, "r");
	char *line = NULL;
	char *tok = NULL;
	size_t len = 0;
	ssize_t read;
	uint32_t tid;
	uint64_t start, stop;
	int khr;
	khint_t k;
	while ((read = getline(&line, &len, ifp)) != -1) {
		if(line[0] == '\0' || line[0] == '#') // Empty line or comment line
			continue;
		tok = strtok(line, "\t");
		tid = (uint32_t)bam_name2id(header, tok);
		tok = strtok(NULL, "\t");
		start = strtoull(tok, NULL, 10);
		tok = strtok(NULL, "\t");
		stop = strtoull(tok, NULL, 10);
		k = kh_get(bed, ret, tid);
		if(k == kh_end(ret)) {
			k = kh_put(bed, ret, tid, &khr);
			kh_val(ret, k).intervals = (uint64_t *)calloc(1, sizeof(uint64_t));
			kh_val(ret, k).intervals[0] = to_ivl(start - padding, stop + padding);
			kh_val(ret, k).n = 1;
		}
		else {
			kh_val(ret, k).intervals = (uint64_t *)realloc(kh_val(ret, k).intervals, ++kh_val(ret, k).n * sizeof(uint64_t));
			if(!kh_val(ret, k).intervals) {
				fprintf(stderr, "[E:%s] Could not allocate memory. Abort mission!\n", __func__);
				exit(EXIT_FAILURE);
			}
			kh_val(ret, k).intervals[kh_val(ret, k).n - 1] = to_ivl(start - padding, stop + padding);
#if !NDEBUG
			fprintf(stderr, "[D:%s] Number of intervals in bed file "
					"for contig %u, ('%s'): %lu\n", __func__, tid, header->target_name[tid], kh_val(ret, k).n);
#endif
		}
	}
	sort_bed(ret);
	return ret;
}

size_t get_nregions(khash_t(bed) *h)
{
	size_t ret = 0uL;
	for(khiter_t ki = kh_begin(h); ki != kh_end(h); ++ki){
		if(!kh_exist(h, ki))
			continue;
		ret += kh_val(h, ki).n;
	}
	return ret;
}


/*
 * Compare two interval objects
 */
int intcmp(const void *a, const void *b)
{
	return *(uint64_t *)a - *(uint64_t *)b;
}


void bed_destroy_hash(void *arg)
{
	khash_t(bed) *b = (khash_t(bed) *)arg;
	khint_t ki;
	for(ki = kh_begin(b); ki != kh_end(b); ++ki) {
		if(!kh_exist(b, ki))
			continue;
		cond_free(kh_val(b, ki).intervals);
		kh_val(b, ki).n = 0;
	}
	cond_free(b);
}

void sort_bed(khash_t(bed) *bed)
{
	for(khint_t k = kh_begin(bed); k != kh_end(bed); ++k) {
		if(!kh_exist(bed, k))
			continue;
		qsort(kh_val(bed, k).intervals, kh_val(bed, k).n, sizeof(uint64_t), &intcmp);
	}
	return;
}

