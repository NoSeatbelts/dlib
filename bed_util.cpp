#include "bed_util.h"

#define N_IVL_BINS 5

#ifdef __cplusplus
std::vector<std::pair<khint_t, khiter_t>> make_sorted_keys(khash_t(bed) *h) {
	std::vector<std::pair<khint_t, khiter_t>> keyset;
	for(khiter_t ki = kh_begin(aux->bed); ki != kh_end(h); ++ki) {
		if(kh_exist(h, ki)) keyset.push_back(std::pair<khint_t, khiter_t>(kh_key(h, ki), ki));
	}
#	ifdef __GNUC__
			__gnu_parallel::sort(keyset.begin(), keyset.end(), [](std::pair<khint_t, khiter_t> p1, std::pair<khint_t, khiter_t> p2) {
#	else
			std::sort(keyset.begin(), keyset.end(), [](std::pair<khint_t, khiter_t> p1, std::pair<khint_t, khiter_t> p2) {
#	endif
				return p1.first < p2.first;
	});
	return keyset;
}
#endif

void print_bed_hash(FILE *fp, khash_t(bed) *bed)
{
	for(khiter_t k = kh_begin(bed); k != kh_end(bed); ++k)
		if(kh_exist(bed, k))
			for(unsigned j = 0; j < kh_val(bed, k).n; ++j)
				fprintf(fp, "Contig: %s. Start: %i. Stop: %i\n",
						kh_val(bed, k).contig_name, get_start(kh_val(bed, k).intervals[j]),
						get_stop(kh_val(bed, k).intervals[j]));
}

khash_t(bed) *build_ref_hash(bam_hdr_t *header) {
	khiter_t k;
	int khr;
	khash_t(bed) *ret = kh_init(bed);
	uint64_t binsize;
	for(int i = 0; i < header->n_targets; ++i) {
		k = kh_put(bed, ret, i, &khr);
		kh_val(ret, k).n = N_IVL_BINS;
		kh_val(ret, k).intervals = (uint64_t *)calloc(N_IVL_BINS, sizeof(uint64_t));
		binsize = (header->target_len[i] - 1) / N_IVL_BINS;
		for(int j = 0; j < N_IVL_BINS - 1; ++j) {
			LOG_INFO("ivl sizes: %lu, %lu.\n", (uint64_t)(j * binsize), (uint64_t)((j + 1) * binsize) - 1);
			if((j + 1) * binsize >= header->target_len[i] - 1)
				kh_val(ret, k).intervals[j] = to_ivl((uint64_t)(j * binsize), header->target_len[i] - 1);
			else
				kh_val(ret, k).intervals[j] = to_ivl((uint64_t)(j * binsize), (uint64_t)((j + 1) * binsize) - 1);
			kh_val(ret, k).intervals[j] = to_ivl((uint64_t)(j * binsize), ((j + 1) * binsize >= header->target_len[i] - 1) ? (uint64_t)header->target_len[i] - 1:(uint64_t)((j + 1) * binsize) - 1);
		}
		kh_val(ret, k).contig_name = strdup(header->target_name[i]);
	}
#if !NDEBUG
	print_bed_hash(stderr, ret);
#endif
	return ret;
}


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
	size_t region_num = 0;
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
			kstring_t ks = {0, 0, NULL};
			ksprintf(&ks, "|%s|tid:%u|region_num:%lu|", ((tok = strtok(NULL, "\t")) != NULL) ? tok: NO_ID_STR, kh_key(ret, k), ++region_num);
			kh_val(ret, k).contig_name = ks_release(&ks);
			fputs(kh_val(ret, k).contig_name, stderr);
		} else {
			kh_val(ret, k).intervals = (uint64_t *)realloc(kh_val(ret, k).intervals, ++kh_val(ret, k).n * sizeof(uint64_t));
			if(!kh_val(ret, k).intervals) {
                LOG_ERROR("Could not allocate memory. Abort mission!\n");
			}
			kh_val(ret, k).intervals[kh_val(ret, k).n - 1] = to_ivl(start - padding, stop + padding);
		}
	}
	sort_bed(ret);
	return ret;
}

size_t get_nregions(khash_t(bed) *h)
{
	size_t ret = 0uL;
	for(khiter_t ki = kh_begin(h); ki != kh_end(h); ++ki)
		if(kh_exist(h, ki))
			ret += kh_val(h, ki).n;
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
		if(!kh_exist(b, ki)) continue;
		cond_free(kh_val(b, ki).intervals);
		cond_free(kh_val(b, ki).contig_name);
		kh_val(b, ki).n = 0;
	}
	kh_destroy(bed, b);
}

void sort_bed(khash_t(bed) *bed)
{
	for(khint_t k = kh_begin(bed); k != kh_end(bed); ++k)
		if(kh_exist(bed, k))
			qsort(kh_val(bed, k).intervals, kh_val(bed, k).n, sizeof(uint64_t), &intcmp);
	return;
}

