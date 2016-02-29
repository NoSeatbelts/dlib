#include "vcf_util.h"

void bcf_add_bam_contigs(bcf_hdr_t *hdr, bam_hdr_t *src) {
    for(int i = 0; i < src->n_targets; ++i) {
        bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_CTG, "ID", (char*)src->target_name[i], NULL);
        if ( hrec ) continue;
        hrec = (bcf_hrec_t*) calloc(1,sizeof(bcf_hrec_t));
        hrec->key = strdup("contig");
        bcf_hrec_add_key(hrec, "ID", strlen("ID"));
        bcf_hrec_set_val(hrec, hrec->nkeys-1, (char*)src->target_name[i], strlen(src->target_name[i]), 0);
        bcf_hdr_add_hrec(hdr, hrec);
    }
}
