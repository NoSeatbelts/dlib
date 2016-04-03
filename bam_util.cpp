#include "bam_util.h"


#ifdef __cplusplus
namespace dlib {
#endif


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

    std::string bam2cppstr(bam1_t *b, std::string& qname)
    {
        char *qual, *seqbuf;
        int i;
        uint8_t *seq, *rvdata;
        uint32_t *pv, *fa;
        int8_t t;
        kstring_t ks = {0, 0, nullptr};
        ksprintf(&ks, "@%s PV:B:I", qname.c_str());
        pv = (uint32_t *)dlib::array_tag(b, "PV");
        fa = (uint32_t *)dlib::array_tag(b, "FA");
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", pv[i]);
        kputs("\tFA:B:I", &ks);
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", fa[i]);
        ksprintf(&ks, "\tFM:i:%i\tFP:i:%i\tNC:i:%i",
                bam_itag(b, "FM"), bam_itag(b, "FP"), bam_itag(b, "NC"));
        if((rvdata = bam_aux_get(b, "RV")) != nullptr)
            ksprintf(&ks, "\tRV:i:%i", bam_aux2i(rvdata));
        kputc('\n', &ks);
        seq = bam_get_seq(b);
        seqbuf = (char *)malloc(b->core.l_qseq + 1);
        for (i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = seq_nt16_str[bam_seqi(seq, i)];
        seqbuf[i] = '\0';
        if (b->core.flag & BAM_FREVERSE) { // reverse complement
            for(i = 0; i < b->core.l_qseq>>1; ++i) {
                t = seqbuf[b->core.l_qseq - i - 1];
                seqbuf[b->core.l_qseq - i - 1] = nuc_cmpl(seqbuf[i]);
                seqbuf[i] = nuc_cmpl(t);
            }
            if(b->core.l_qseq&1) seqbuf[i] = nuc_cmpl(seqbuf[i]);
        }
        seqbuf[b->core.l_qseq] = '\0';
        assert(strlen(seqbuf) == (uint64_t)b->core.l_qseq);
        kputs(seqbuf, &ks);
        kputs("\n+\n", &ks);
        qual = (char *)bam_get_qual(b);
        for(i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = 33 + qual[i];
        if (b->core.flag & BAM_FREVERSE) { // reverse
            for (i = 0; i < b->core.l_qseq>>1; ++i) {
                t = seqbuf[b->core.l_qseq - 1 - i];
                seqbuf[b->core.l_qseq - 1 - i] = seqbuf[i];
                seqbuf[i] = t;
            }
        }
        assert(strlen(seqbuf) == (uint64_t)b->core.l_qseq);
        kputs(seqbuf, &ks), free(seqbuf);
        kputc('\n', &ks);
        std::string ret(ks.s), free(ks.s);
        return ret;
    }

    std::string bam2cppstr(bam1_t *b)
    {
        char *qual, *seqbuf;
        int i;
        uint8_t *seq, *rvdata;
        uint32_t *pv, *fa;
        int8_t t;
        kstring_t ks{0, 0, nullptr};
        ksprintf(&ks, "@%s PV:B:I", bam_get_qname(b));
        pv = (uint32_t *)dlib::array_tag(b, "PV");
        fa = (uint32_t *)dlib::array_tag(b, "FA");
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", pv[i]);
        kputs("\tFA:B:I", &ks);
        for(i = 0; i < b->core.l_qseq; ++i) ksprintf(&ks, ",%u", fa[i]);
        ksprintf(&ks, "\tFM:i:%i\tFP:i:%i", bam_itag(b, "FM"), bam_itag(b, "FP"));
        if((rvdata = bam_aux_get(b, "NC")) != nullptr)
            ksprintf(&ks, "\tNC:i:%i", bam_aux2i(rvdata));
        if((rvdata = bam_aux_get(b, "RV")) != nullptr)
            ksprintf(&ks, "\tRV:i:%i", bam_aux2i(rvdata));
        kputc('\n', &ks);
        seq = bam_get_seq(b);
        seqbuf = (char *)malloc(b->core.l_qseq + 1);
        for (i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = seq_nt16_str[bam_seqi(seq, i)];
        seqbuf[i] = '\0';
        if (b->core.flag & BAM_FREVERSE) { // reverse complement
            for(i = 0; i < b->core.l_qseq>>1; ++i) {
                t = seqbuf[b->core.l_qseq - i - 1];
                seqbuf[b->core.l_qseq - i - 1] = nuc_cmpl(seqbuf[i]);
                seqbuf[i] = nuc_cmpl(t);
            }
            if(b->core.l_qseq&1) seqbuf[i] = nuc_cmpl(seqbuf[i]);
        }
        seqbuf[b->core.l_qseq] = '\0';
        assert(strlen(seqbuf) == (uint64_t)b->core.l_qseq);
        kputs(seqbuf, &ks);
        kputs("\n+\n", &ks);
        qual = (char *)bam_get_qual(b);
        for(i = 0; i < b->core.l_qseq; ++i) seqbuf[i] = 33 + qual[i];
        if (b->core.flag & BAM_FREVERSE) { // reverse
            for (i = 0; i < b->core.l_qseq>>1; ++i) {
                t = seqbuf[b->core.l_qseq - 1 - i];
                seqbuf[b->core.l_qseq - 1 - i] = seqbuf[i];
                seqbuf[i] = t;
            }
        }
        assert(strlen(seqbuf) == (uint64_t)b->core.l_qseq);
        kputs(seqbuf, &ks), free(seqbuf);
        kputc('\n', &ks);
        std::string ret(ks.s), free(ks.s);
        return ret;
    }

    std::string get_SO(bam_hdr_t *hdr) {
        char *end, *so_start;
        std::string ret;
        if (strncmp(hdr->text, "@HD", 3) != 0) goto NA;
        if ((end = strchr(hdr->text, '\n')) == 0) goto NA;
        *end = '\0';

        if((so_start = strstr(hdr->text, "SO:")) == nullptr) goto NA;
        ret = std::string(so_start + strlen("SO:"));
        *end = '\n';
        return ret;

        NA:
        LOG_WARNING("Sort order not found. Returning N/A.\n");
        return std::string("N/A");
    }

    void resize_stack(tmp_stack_t *stack, size_t n) {
        if(n > stack->max) {
            stack->max = n;
            stack->a = (bam1_t **)realloc(stack->a, sizeof(bam1_t *) * n);
            if(!stack->a) LOG_EXIT("Failed to reallocate memory for %lu bam1_t * objects. Abort!\n", stack->max);
        } else if(n < stack->n){
            for(uint64_t i = stack->n;i > n;) free(stack->a[--i]->data);
            stack->max = n;
            stack->a = (bam1_t **)realloc(stack->a, sizeof(bam1_t *) * n);
        }
    }

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


void abstract_pair_iter(samFile *in, bam_hdr_t *hdr, samFile *ofp, pair_aux_fn function, void *aux)
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
        if(function(b1, b, aux) == 0)
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

#ifdef __cplusplus
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
        LOG_EXIT("Broken function.\n");
        int ret;
        bam1_t *r1 = bam_init1();
        while(next() >= 0) {
            if(rec->core.flag & BAM_FREAD1) {
                bam_copy1(r1, rec);
                continue;
            }
            if((rec->core.flag & BAM_FREAD2) == 0) continue;
            assert(strcmp(bam_get_qname(r1), bam_get_qname(rec)) == 0);
            if((ret = fn(r1, rec, data)) != 0) continue;
            if(strcmp(bam_get_qname(r1), bam_get_qname(rec)))
                LOG_EXIT("WTF\n");
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
        BamHandle in;
        BamHandle out;
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
} /* namespace dlib */
#endif /* ifdef __cplusplus */
