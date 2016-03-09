#ifndef BAM_UTIL_H
#define BAM_UTIL_H
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include "htslib/sam.h"
#include "char_util.h"
#include "compiler_util.h"
#include "io_util.h"
#include "logging_util.h"
#include "misc_util.h"
#include "bed_util.h"
#ifdef __cplusplus
#include <functional>
#include <unordered_set>
#endif

typedef void (*pair_fn)(bam1_t *b, bam1_t *b1);
typedef int (*pair_aux_fn)(bam1_t *b, bam1_t *b1, void *data);
typedef void (*single_fn)(bam1_t *b);
typedef void (*single_aux)(bam1_t *b, void *data);
typedef int (*single_aux_check)(bam1_t *b, void *data);
typedef int (*plp_fn)(const bam_pileup1_t *plp, int n_plp, void *data);

#ifdef __cplusplus
namespace dlib {
    void abstract_pair_set(samFile *in, bam_hdr_t *hdr, samFile *ofp, std::unordered_set<pair_fn> functions);
    std::string get_SO(bam_hdr_t *hdr);
    std::string bam2cppstr(bam1_t *b);
    class BamRec {
    public:
        bam1_t *b;
        BamRec(): b(bam_init1()){
            LOG_DEBUG("Initialized rec at pointer %p.\n", (void *)b);
        }
        // Copy
        BamRec(bam1_t *b) : b(bam_dup1(b))
        {
        }
        BamRec(BamRec& other) :
        b(bam_dup1(other.b))
        {
        }
        ~BamRec() {
            if(b) bam_destroy1(b);
        }
    };
    class BamHandle;
    class BedPlpAuxBase {
    public:
        unsigned minMQ;
        int padding;
        BamHandle *handle;
        BedPlpAuxBase(unsigned _minMQ, BamHandle *_handle=NULL, int _padding=DEFAULT_PADDING):
            minMQ(_minMQ),
            padding(_padding),
            handle(_handle)
        {
        }
    };
    class BamHandle {
    public:
        int is_write;
        samFile *fp;
        hts_itr_t *iter;
        bam_hdr_t *header;
        hts_idx_t *idx;
        const bam_pileup1_t *pileups;
        bam_plp_t plp;
        bam1_t *rec;
        // Read constructor
        BamHandle(const char *path):
            is_write(0),
            fp(sam_open(path, "r")),
            iter(NULL),
            header(sam_hdr_read(fp)),
            idx(bam_index_load(path)),
            pileups(NULL),
            plp(NULL),
            rec(bam_init1())
        {
            if(!fp) LOG_EXIT("Could not open input bam %s for reading. Abort!\n", path);
            if(!idx) LOG_WARNING("Could not load index file for input bam, just FYI.\n");
        }
        // Write constructor
        BamHandle(const char *path, bam_hdr_t *hdr, const char *mode = "wb"):
            is_write(1),
            fp(sam_open(path, mode)),
            iter(NULL),
            header(bam_hdr_dup(hdr)),
            idx(NULL),
            pileups(NULL),
            plp(NULL),
            rec(NULL)
        {
            sam_hdr_write(fp, header);
            if(fp == NULL) LOG_EXIT("Could not open output bam %s for reading. Abort!\n", path);
        }
        ~BamHandle() {
            if(fp) sam_close(fp), fp = NULL;
            if(iter) hts_itr_destroy(iter), iter = NULL;
            if(header) bam_hdr_destroy(header);
            if(idx) hts_idx_destroy(idx);
            if(plp) bam_plp_destroy(plp);
            if(rec) bam_destroy1(rec);
        }
        int for_each_pair(std::function<int (bam1_t *, bam1_t *, void *)> fn, BamHandle& ofp, void *data=NULL);
        int for_each(std::function<int (bam1_t *, void *)> fn, BamHandle& ofp, void *data=NULL);
        int write();
        int read(BamRec b);
        int write(bam1_t *b) {
            return sam_write1(fp, header, b);
        }
        int write(BamRec b) {
            return write(b.b);
        }
        int read(bam1_t *b) {
            return iter ? bam_itr_next(fp, iter, b)
                        : sam_read1(fp, header, b);
        }
        int next() {
            if(read(rec) < 0) {
                LOG_INFO("StopIteration: Finished iterating through bam %s.\n", fp->fn);
                return -1;
            }
            return 1;
        }
        int bed_plp_auto(khash_t(bed) *bed, std::function<int (const bam_pileup1_t *, int, void *)> fn,
                         BedPlpAuxBase *auxen);
    };
    static inline int bam_readrec(BedPlpAuxBase *data, bam1_t *b) {
        return data->handle->iter ? bam_itr_next(data->handle->fp, data->handle->iter, data->handle->rec)
                                  : sam_read1(data->handle->fp, data->handle->header, data->handle->rec);
    }
    static int read_bam(BedPlpAuxBase *data, bam1_t *b) {
        int ret;
        for(;;) {
            if((ret = bam_readrec(data, b)) >= 0) break;
            if(b->core.flag & (BAM_FSECONDARY | BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP))
                continue;
        }
        return ret;
    }

    int bam_apply_function(char *infname, char *outfname,
                           std::function<int (bam1_t *, void *)> func, void *data=NULL, const char *mode="wb");
    int bam_pair_apply_function(char *infname, char *outfname,
            pair_aux_fn fn, void *data=NULL, const char *mode="wb");
    /*
     * Finds the index in an array tag to use for a bam_pileup1_t struct.
     */
    static inline int arr_qpos(const bam_pileup1_t *plp)
    {
        return (plp->b->core.flag & BAM_FREVERSE) ? plp->b->core.l_qseq - 1 - plp->qpos
                                                  : plp->qpos;
    }

    static const char *dlib_tags[] = {
            "MU",
            "SU",
            "LM",
            "LR",
            "SC",
            "ML",
            "AF",
            "FM"
    };

    static inline void nuke_dlib_tags(bam1_t *b) {
        uint8_t *data;
        for(auto tag: dlib_tags)
            if((data = bam_aux_get(b, tag)) != NULL)
                bam_aux_del(b, data);
    }

#endif /* ifdef __cplusplus */


/*
 *bam_seqi_cmpl: returns the complement of bam_seqi.
 **/
#ifdef __cplusplus
#    define bam_seqi_cmpl(seq, index) dlib::seq_nt16_rc[bam_seqi(seq, index)]
#else
#    define bam_seqi_cmpl(seq, index) seq_nt16_rc[bam_seqi(seq, index)]
#endif




#ifndef SEQ_TABLE_DEFS
#define SEQ_TABLE_DEFS
static const int8_t seq_comp_table[] = {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};
static const int8_t nt16_num_table[] = {0, 1, 2, 0, 3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 5};
static const uint8_t seq_nt16_rc[] = {15, 8, 4, 15, 2, 15, 15, 15, 1, 15, 15, 15, 15, 15, 15, 15};
#endif

#define seqnt2num(character) dlib::nt16_num_table[(int8_t)character]

// Like bam_endpos, but doesn't check that the read is mapped, as that's already been checked.
#ifndef bam_getend
#define bam_getend(b) ((b)->core.pos + bam_cigar2rlen((b)->core.n_cigar, bam_get_cigar(b)))
#endif


typedef struct tmp_stack {
    size_t n, max;
    bam1_t **a;
} tmp_stack_t;

static inline void stack_insert(tmp_stack_t *stack, bam1_t *b)
{
    if (stack->n == stack->max) {
        stack->max = stack->max? stack->max<<1 : 0x10000;
        stack->a = (bam1_t**)realloc(stack->a, sizeof(bam1_t*) * stack->max);
    }
    stack->a[stack->n++] = bam_dup1(b);
}

static inline void add_unclipped_mate_starts(bam1_t *b1, bam1_t *b2);
void abstract_pair_iter(samFile *in, bam_hdr_t *hdr, samFile *ofp, pair_fn function);
void abstract_single_filter(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux_check function, void *data);
void abstract_single_data(samFile *in, bam_hdr_t *hdr, samFile *out, single_aux function, void *data);
void abstract_single_iter(samFile *in, bam_hdr_t *hdr, samFile *out, single_fn function);
void resize_stack(tmp_stack_t *stack, size_t n);

static inline void seq_nt16_cpy(char *read_str, uint8_t *seq, int len, int is_rev) {
    if(is_rev) {
        for(--len;len != -1; --len) *read_str++ = seq_nt16_str[bam_seqi_cmpl(seq, len)];
        *read_str++ = '\0';
    } else {
        read_str += len;
        *read_str-- = '\0';
        for(--len;len != -1; --len) *read_str-- = seq_nt16_str[bam_seqi(seq, len)];
    }
}

static inline void bam_seq_cpy(char *read_str, bam1_t *b) {
    seq_nt16_cpy(read_str, (uint8_t *)bam_get_seq(b), b->core.l_qseq, b->core.flag & BAM_FREVERSE);
}


static inline int32_t int_tag_zero(uint8_t *data)
{
    return data ? bam_aux2i(data): 0;
}

CONST static inline int32_t get_unclipped_start(bam1_t *b)
{
    if(b->core.flag & BAM_FUNMAP) return -1;
    const uint32_t *cigar = bam_get_cigar(b);
    int32_t ret = b->core.pos;
    for(int i = 0; i < b->core.n_cigar; ++i) {
        switch(bam_cigar_op(cigar[i])) {
            case BAM_CSOFT_CLIP:
            case BAM_CDEL:
            case BAM_CREF_SKIP:
            case BAM_CPAD:
                ret -= bam_cigar_oplen(cigar[i]); break;
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                return ret;
            /*
            case BAM_CINS:
            case BAM_CHARD_CLIP:
                // DO nothing
            */
        }
    }
    return ret;
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

int bampath_has_tag(char *bampath, const char *tag);
void check_bam_tag_exit(char *bampath, const char *tag);


CONST static inline void *array_tag(bam1_t *b, const char *tag) {
    uint8_t *data = bam_aux_get(b, tag);
    if(!data) {
        LOG_EXIT("Missing tag %s. Abort!\n", tag);
    }
    char tagtype = *data++;
    if(UNLIKELY(tagtype != 'B')) LOG_EXIT("Incorrect byte %c where B expected in array tag for key %s. Abort!\n", tagtype, tag);
    switch(*data++) {
        case 'i': case 'I': case 's': case 'S': case 'f': case 'c': case 'C':
            break;
        default:
            LOG_EXIT("Unrecognized tag type %c.\n", *(data - 1));
    }
    return (void *)(data + sizeof(int));
}

#define cigarop_sc_len(cigar) ((((cigar) & 0xfU) == BAM_CSOFT_CLIP) ? bam_cigar_oplen(cigar): 0)

CONST static inline int bam_sc_len_cigar(bam1_t *b, uint32_t *cigar, int n_cigar)
{
    const int clen1 = cigarop_sc_len(cigar[0]);
    const int clen2 = cigarop_sc_len(cigar[n_cigar - 1]);
    return MAX2(clen1, clen2);
}

CONST static inline int bam_sc_len(bam1_t *b)
{
    return (b->core.flag & BAM_FUNMAP) ? 0: bam_sc_len_cigar(b, bam_get_cigar(b), b->core.n_cigar);
}

CONST static inline float bam_frac_align(bam1_t *b)
{
    if(b->core.flag & BAM_FUNMAP) return 0.;
    int sum = 0;
    uint32_t *cigar = bam_get_cigar(b);
    for(unsigned i = 0; i < b->core.n_cigar; ++i)
        if(bam_cigar_op(cigar[i]) & (BAM_CMATCH | BAM_CEQUAL | BAM_CDIFF))
            sum += bam_cigar_oplen(cigar[i]);
    return (float)sum / b->core.l_qseq;
}

static inline void add_sc_lens(bam1_t *b1, bam1_t *b2) {
       const int sc1 = bam_sc_len(b1); const int sc2 = bam_sc_len(b2);
       bam_aux_append(b2, "SC", 'i', sizeof(int), (uint8_t *)&sc2);
       bam_aux_append(b2, "ML", 'i', sizeof(int), (uint8_t *)&sc1);
       bam_aux_append(b1, "SC", 'i', sizeof(int), (uint8_t *)&sc1);
       bam_aux_append(b1, "ML", 'i', sizeof(int), (uint8_t *)&sc2);
}

static inline void add_qseq_len(bam1_t *b1, bam1_t *b2) {
       bam_aux_append(b2, "LM", 'i', sizeof(int), (uint8_t *)&b1->core.l_qseq);
       bam_aux_append(b1, "LM", 'i', sizeof(int), (uint8_t *)&b2->core.l_qseq);
}

/*  @func add_unclipped_mate_starts
 *  @abstract Adds the unclipped start positions for each read and its mate
 */
static inline void add_fraction_aligned(bam1_t *b1, bam1_t *b2) {
       const float frac1 = bam_frac_align(b1);
       const float frac2 = bam_frac_align(b2);
       bam_aux_append(b2, "AF", 'f', sizeof(float), (uint8_t *)&frac2);
       bam_aux_append(b2, "MF", 'f', sizeof(float), (uint8_t *)&frac1);
       bam_aux_append(b1, "AF", 'f', sizeof(float), (uint8_t *)&frac1);
       bam_aux_append(b1, "MF", 'f', sizeof(float), (uint8_t *)&frac2);
}

void bam_plp_set_maxcnt(bam_plp_t, int);

enum htseq {
    HTS_A = 1,
    HTS_C = 2,
    HTS_G = 4,
    HTS_T = 8,
    HTS_N = 15
};

#ifdef __cplusplus
} // namespace dlib
#endif

// bam utility macros.

/* true if both in pair are unmapped, false if only one.
 */
#define bam_pair_unmapped(b2) (((b2)->core.flag & BAM_FPAIR_UNMAPPED) == BAM_FPAIR_UNMAPPED)

/* @func bam_set_base sets the nucleotide at index i in read p to be set to base at index i in read b.
 * :param: p [bam1_t *] One bam record
 * :param: b [char] Nucleotide to set
 * :param: i [index] Base position in read
 */
#define set_base(pSeq, base, i) (pSeq)[(i)>>1] = ((seq_nt16_table[(int8_t)base] << (((~i) & 1) << 2)) | (((pSeq)[(i)>>1]) & (0xf0U >> (((~i) & 1) << 2))))

/* @func bam_set_base sets the nucleotide at index i in read p to be set to base at index i in read b.
 * :param: pSeq [uint8_t *] One bam record's seq ptr
 * :param: bSeq [uint8_t *] Second bam record's seq ptr
 * :param: i [index] Base position in read
 */
#define bam_set_base(pSeq, bSeq, i) (pSeq)[(i)>>1] = ((bam_seqi(bSeq, i) << (((~i) & 1) << 2)) | (((pSeq)[(i)>>1]) & (0xf0U >> (((~i) & 1) << 2))))

/* @func n_base sets the nucleotide at index i in read p to N.
 * :param: p [uint8_t *] One bam record's seq ptr
 * :param: i [index] Base position in read
 */
#define n_base(pSeq, i) pSeq[(i)>>1] |= (0xf << ((~(i) & 1) << 2));


/* @func bam_itag gets an integer tag for a given key.
 * Warning: will segfault if not present! Use check_bam_tag_exit to check first.
 * :param: b [bam1_t *] Bam record
 * :param: key [const char *] Key for tag
 */
#define bam_itag(b, key) bam_aux2i(bam_aux_get(b, key))

/* Just an array-checking utility for debugging. I don't see much use for this. */
#define check_fa(arr, fm, len) \
    do {\
        for(int i##arr = 0; i##arr < len; ++i##arr) {\
            if(arr[i##arr] > fm){\
                LOG_EXIT((char *)"%u arr value greater than FM %u.\n", arr[i##arr], fm);\
            }\
        }\
    } while(0)


#endif // BAM_UTIL_H
