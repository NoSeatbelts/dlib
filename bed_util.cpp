#include "bed_util.h"
#include <zlib.h>

#ifdef __cplusplus
    #include <algorithm>

    namespace dlib {
    std::vector<khiter_t> make_sorted_keys(khash_t(bed) *h) {
        std::vector<std::pair<khint_t, khiter_t>> keyset;
        for(khiter_t ki = kh_begin(aux->bed); ki != kh_end(h); ++ki)
            if(kh_exist(h, ki))
                keyset.push_back(std::pair<khint_t, khiter_t>(kh_key(h, ki), ki));
        std::sort(keyset.begin(), keyset.end(), [](std::pair<khint_t, khiter_t> p1, std::pair<khint_t, khiter_t> p2) {
                return p1.first < p2.first;
        });
        std::vector<khiter_t> ret;
        for(auto tup: keyset) ret.push_back(tup.second);
        return ret;
    }

    template<typename K, typename V>
    std::vector<K> make_sorted_keys(std::unordered_map<K,V> hashmap) {
        std::vector<std::pair<K,V>> keyset;
        keyset.reserve(hashmap.size());
        for(auto& pair: hashmap) keyset.push_back(pair);
        std::sort(keyset.begin(), keyset.end(), [](std::pair<K,V> p1, std::pair<K,V> p2){
            return p1.second < p2.second;
        });
        std::vector<K> ret();
        ret.reserve(keyset.size());
        for(auto& pair: keyset) ret.push_back(pair.first);
        return ret;
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
        for(int i = 0; i < header->n_targets; ++i) {
            k = kh_put(bed, ret, i, &khr);
            kh_val(ret, k).n = 1;
            kh_val(ret, k).intervals = (uint64_t *)calloc(1, sizeof(uint64_t));
            *kh_val(ret, k).intervals = to_ivl(0, header->target_len[i] - 1);
            kh_val(ret, k).contig_name = strdup(header->target_name[i]);
        }
    #if !NDEBUG
        print_bed_hash(stderr, ret);
    #endif
        sort_bed(ret);
        return ret;
    }

    void sort_bed_hash(khash_t(bed) *hash) {
        for(khiter_t k = kh_begin(hash); k != kh_end(hash); ++k)
            if(kh_exist(hash, k))
#ifdef __cplusplus
                std::sort(kh_val(hash, k).intervals, kh_val(hash, k).intervals + kh_val(hash, k).n, [](uint64_t a, uint64_t b){
                    return get_start(a) == get_start(b) ? get_stop(a) < get_stop(b)
                                                        : get_start(a) < get_start(b);
                });
#else
                qsort(kh_val(hash, k).intervals, kh_val(hash, k).n, &intcmp);
#endif
    }

    int faidx_seq_num(faidx_t *fai, const char *name) {
        if(!faidx_has_seq(fai, name)) return -1; // Not found. Abort!
        int n = faidx_nseq(fai);
        for(int i = 0; i < n; ++i)
            if(strcmp(faidx_iseq(fai, i), name) == 0)
                return i;
        return -1; // This never happens
    }


    khash_t(bed) *parse_bed_hash(const char *path, faidx_t *fai, uint32_t padding)
    {
        if(!path) return NULL;
        khash_t(bed) *ret = kh_init(bed);
        gzFile ifp = gzopen(path, "rb");
        if(ifp == Z_NULL) {
            LOG_INFO("Could not open file at %s for reading. Returning NULL.\n");
            return NULL;
        }
        const size_t bufsize = 15000;
        char *line = (char *)malloc(bufsize);
        char *tmp = NULL;
        char *tok = NULL;
        uint32_t tid;
        uint64_t start, stop;
        int khr;
        size_t region_num = 0;
        khint_t k;
        while ((tmp = gzgets(ifp, line, bufsize)) != NULL) {
            switch(*line) {
                case '\0': case '#': continue;
            }
            tok = strtok(line, "\t");
            if((tid = (uint32_t)faidx_seq_num(fai, tok)) == (uint32_t)-1) LOG_EXIT("tid %s not found in fasta reference. Abort!\n", tok);
            tok = strtok(NULL, "\t");
            start = strtoull(tok, NULL, 10);
            tok = strtok(NULL, "\t");
            stop = strtoull(tok, NULL, 10);
            if((k = kh_get(bed, ret, tid)) == kh_end(ret)) {
                k = kh_put(bed, ret, tid, &khr);
                kh_val(ret, k).intervals = (uint64_t *)calloc(1, sizeof(uint64_t));
                *kh_val(ret, k).intervals = to_ivl(start - padding, stop + padding);
                kh_val(ret, k).n = 1;
                kstring_t ks = {0, 0, NULL};
                ksprintf(&ks, "|%s|tid:%u|region_num:%lu|", ((tok = strtok(NULL, "\t")) != NULL) ? tok: NO_ID_STR, kh_key(ret, k), ++region_num);
                kh_val(ret, k).contig_name = ks_release(&ks);
            } else {
                kh_val(ret, k).intervals = (uint64_t *)realloc(kh_val(ret, k).intervals,
                                                               (kh_val(ret, k).n + 1)* sizeof(uint64_t));
                if(!kh_val(ret, k).intervals) LOG_EXIT("Could not allocate memory. Abort mission!\n");
                kh_val(ret, k).intervals[kh_val(ret, k).n++] = to_ivl(start - padding, stop + padding);
            }
        }
        free(line);
        gzclose(ifp);
        sort_bed(ret);
        return ret;
    }

    khash_t(bed) *parse_bed_hash(const char *path, bam_hdr_t *header, uint32_t padding)
    {
        LOG_DEBUG("Path to bed: %s.\n", path);
        if(!path) return NULL;
        khash_t(bed) *ret = kh_init(bed);
        gzFile ifp = gzopen(path, "rb");
        if(ifp == Z_NULL) {
            LOG_INFO("Could not open file at %s for reading. Returning NULL.\n");
            return NULL;
        }
        const size_t bufsize = 15000;
        char *line = (char *)malloc(bufsize);
        char *tmp = NULL;
        char *tok = NULL;
        uint32_t tid;
        uint64_t start, stop;
        int khr;
        size_t region_num = 0;
        khint_t k;
        while ((tmp = gzgets(ifp, line, bufsize)) != NULL) {
            switch(*line) {
                case '\0': case '#': continue;
            }
            tok = strtok(line, "\t");
            tid = (uint32_t)bam_name2id(header, tok);
            tok = strtok(NULL, "\t");
            start = strtoull(tok, NULL, 10);
            tok = strtok(NULL, "\t");
            stop = strtoull(tok, NULL, 10);
            if((k = kh_get(bed, ret, tid)) == kh_end(ret)) {
                k = kh_put(bed, ret, tid, &khr);
                kh_val(ret, k).intervals = (uint64_t *)calloc(1, sizeof(uint64_t));
                *kh_val(ret, k).intervals = to_ivl(start - padding, stop + padding);
                kh_val(ret, k).n = 1;
                kstring_t ks = {0, 0, NULL};
                ksprintf(&ks, "|%s|tid:%u|region_num:%lu|",
                         ((tok = strtok(NULL, "\t")) != NULL) ? tok
                                                              : NO_ID_STR,
                         kh_key(ret, k), ++region_num);
                kh_val(ret, k).contig_name = ks_release(&ks);
            } else {
                kh_val(ret, k).intervals = (uint64_t *)realloc(kh_val(ret, k).intervals,
                                                               (kh_val(ret, k).n + 1)* sizeof(uint64_t));
                if(!kh_val(ret, k).intervals) LOG_EXIT("Could not allocate memory. Abort mission!\n");
                kh_val(ret, k).intervals[kh_val(ret, k).n++] = to_ivl(start - padding, stop + padding);
            }
        }
        free(line);
        gzclose(ifp);
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
                std::sort(kh_val(bed, k).intervals, kh_val(bed, k).intervals + kh_val(bed, k).n);
        return;
    }

#ifdef __cplusplus
} /* namespace dlib */
#endif
