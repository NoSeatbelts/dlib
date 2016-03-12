#include "dlib/bed_util.h"
#include <time.h>

int main(int argc, char *argv[]) {
    samFile *fp = sam_open("test.bam", "r");
    bam_hdr_t *hdr = sam_hdr_read(fp);
    khash_t(bed) *cbed = dlib::parse_bed_hash("test.bed", hdr, 0);
    dlib::ParsedBed cppbed = dlib::ParsedBed("test.bed", hdr, 0);
    clock_t startc, afterc, startcpp, aftercpp;
    startc = clock();
    int c;
    bam1_t *b = bam_init1();
    sam_read1(fp, hdr, b);
    const size_t n_iter = 1uL << 26;
    for(uint64_t i = 0; i < n_iter; ++i)
        c = dlib::bed_test(b, cbed);
    afterc = clock();
    fprintf(stderr, "#s for c style: %f\n", ((double)afterc - startc) / CLOCKS_PER_SEC);
    startcpp = clock();
    for(uint64_t i = 0; i < n_iter; ++i)
        c = cppbed.bam1_test(b);
    aftercpp = clock();
    fprintf(stderr, "#s for cpp style: %f\n", ((double)aftercpp - startcpp) / CLOCKS_PER_SEC);
    dlib::bed_destroy_hash(cbed);
    bam_destroy1(b);
    return 0;
}
