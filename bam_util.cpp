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

void abstract_pair_set(samFile *in, bam_hdr_t *hdr, samFile *ofp, std::unordered_set<pair_fn> functions)
{
    bam1_t *b = bam_init1(), *b1 = bam_init1();
    while (LIKELY(sam_read1(in, hdr, b) >= 0)) {
        if(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        if(b->core.flag & BAM_FREAD1) {
            bam_copy1(b1, b); continue;
        }
        for(auto f: functions) f(b1, b);
        sam_write1(ofp, hdr, b1), sam_write1(ofp, hdr, b);
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
        if(strcmp(bam_get_qname(b1), bam_get_qname(b)))
            LOG_EXIT("Is the bam name sorted? Reads in 'pair' don't have the same name (%s, %s). Abort!\n", bam_get_qname(b1), bam_get_qname(b));
        function(b1, b);
        sam_write1(ofp, hdr, b1), sam_write1(ofp, hdr, b);
    }
    bam_destroy1(b), bam_destroy1(b1);
}

int bampath_has_tag(char *bampath, const char *tag)
{
    samFile *fp = sam_open(bampath, "r");
    bam_hdr_t *header = sam_hdr_read(fp);
    if(!header || !fp) {
        LOG_EXIT("Could not open bam file at '%s' for reading. Abort!\n", bampath);
    }
    bam1_t *b = bam_init1();
    if(sam_read1(fp, header, b) < 0) {
        LOG_EXIT("Empty bam file at '%s'. Abort!\n", bampath);
    }
    int ret = !!bam_aux_get(b, tag);
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(fp);
    return ret;
}

void check_bam_tag_exit(char *bampath, const char *tag)
{
    LOG_DEBUG("Bam at %s has the tag %s.\n", bampath, tag);
    if(!(strcmp(bampath, "-") && strcmp(bampath, "stdin"))) {
            LOG_WARNING("Could not check for bam tag without exhausting a pipe. "
                        "Tag '%s' has not been verified.\n", tag);
            return;
    }
    if(!bampath_has_tag(bampath, tag)) LOG_EXIT("Required bam tag '%s' missing from bam file at path '%s'. Abort!\n", tag, bampath);
}

void check_bam_tag(bam1_t *b, const char *tag)
{
    if(!bam_aux_get(b, tag)) LOG_EXIT((char *)"Required bam tag '%s' not found. Abort mission!\n",tag);
}

#ifdef __cplusplus
namespace dlib {
    /*
     * Applies fn
     */
    int BamHandle::for_each(std::function<int (bam1_t *, void *)> fn, BamHandle& ofp, void *data) {
        int ret;
        while(next() >= 0) {
            if((ret = fn(rec, data)) != 0) continue;
            ofp.write(rec);
        }
        return 0;
    }
    int BamHandle::for_each_pair(std::function<int (bam1_t *, bam1_t *, void *)> fn, BamHandle& ofp, void *data) {
        int ret;
        bam1_t *r1 = bam_init1();
        while(next() >= 0) {
            if(rec->core.flag & BAM_FREAD1) {
                bam_copy1(r1, rec);
                continue;
            }
            assert(strcmp(bam_get_qname(r1), bam_get_qname(rec)) == 0);
            if((ret = fn(r1, rec, data)) != 0) continue;
            ofp.write(r1), ofp.write(rec);
        }
        bam_destroy1(r1);
        return 0;
    }
    int BamHandle::write() {return write(rec);}
    int bam_apply_function(char *infname, char *outfname,
                           std::function<int (bam1_t *, void *)> func, void *data, const char *mode) {
        BamHandle in(infname);
        BamHandle out(outfname, in.header, mode);
        return in.for_each(func, out, data);
    }
    int bam_pair_apply_function(char *infname, char *outfname,
            pair_aux_fn fn, void *data, const char *mode) {
        BamHandle in(infname);
        BamHandle out(outfname, in.header, mode);
        return in.for_each_pair(fn, out, data);
    }
    int BamHandle::bed_plp_auto(khash_t(bed) *bed, std::function<int (const bam_pileup1_t *, int, void *)> fn,
                                BedPlpAuxBase *auxen) {
        plp = bam_plp_init((bam_plp_auto_f)&read_bam, (void *)auxen);
        auxen->handle = this;
        for(khiter_t ki = kh_begin(bed); ki != kh_end(bed); ++ki) {
            if(!kh_exist(bed, ki)) continue;
            for(unsigned i = 0; i < kh_val(bed, ki).n; ++i) {
                const int start = get_start(kh_val(bed, ki).intervals[i]);
                const int stop = get_stop(kh_val(bed, ki).intervals[i]);
                const int bed_tid = (int)kh_key(bed, ki);
                int tid, n_plp, pos;
                if(iter) hts_itr_destroy(iter);
                iter = bam_itr_queryi(idx, bed_tid, start, stop);
                bam_plp_reset(plp);
                while(bam_plp_auto(plp, &tid, &pos, &n_plp) != 0) {
                    assert(tid == bed_tid);
                    if(pos < start) continue;
                    if(pos >= stop) break;
                    fn(pileups, n_plp, auxen);
                }
            }
        }
        return 0;
    }
    /*
     * Does not own aux.
     */
    class BedPlpAux {
        dlib::BamHandle in;
        dlib::BamHandle out;
        BedPlpAuxBase *aux;
        khash_t(bed) *bed;
        BedPlpAux(char *inpath, char *outpath, char *bedpath, BedPlpAuxBase *data):
            in(inpath),
            out(outpath, in.header),
            aux(data),
            bed(parse_bed_hash(bedpath, in.header, data->padding))
        {
            if(!bed) LOG_EXIT("Failed to open bedfile %s. Abort!\n", bedpath);
        }
        ~BedPlpAux() {
            if(bed) bed_destroy_hash(bed);
        }
        int process(std::function<int (const bam_pileup1_t *, int, void *)> fn) {
            return in.bed_plp_auto(bed,fn, aux);
        }
    };
}

#endif /* ifdef __cplusplus */
