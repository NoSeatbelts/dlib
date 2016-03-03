#include "vcf_util.h"

namespace dlib {

    void bcf_add_bam_contigs(bcf_hdr_t *hdr, bam_hdr_t *src) {
        for(int i = 0; i < src->n_targets; ++i)
            if(bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%u>", src->target_name[i], src->target_len[i]))
                LOG_EXIT("Could not add header line '##contig=<ID=%s,length=%u>'", src->target_name[i], src->target_len[i]);
    }

}
