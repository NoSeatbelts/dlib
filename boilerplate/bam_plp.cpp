/*
 * This is untested, though I am attempting to write a generic boilerplate for
 * pileup-based tools.
 */

#include "dlib/bam_util.h"
#include "dlib/bed_util.h"
#include <getopt.h>
int usage() {
    fprintf(stderr, "<-l output_compression_level> in.bam out.bam\n"
                    "Use - for stdin or stdout.\n");
    return EXIT_FAILURE;
}

class BedCovAux : dlib::BedPlpAuxBase {

};



int bed_plp_core(char *infname, char *outfname) {

}

int main(int argc, char *argv[]) {
    if(argc < 3) {
        return usage();
    }
    if(strcmp(argv[1], "--help") == 0) {
        usage();
        return EXIT_SUCCESS;
    }
    int c;
    char out_mode[4] = "wb";
    while((c = get_opt(argc, argv, "l:h?")) > -1) {
        switch(c) {
        case 'l':
            out_mode[2] = atoi(optarg) % 10 + '0'; break;
        case 'h': case '?': usage(); return EXIT_SUCCESS;
        default: LOG_INFO("Unrecognized option %c\n", optopt); return usage();
        }
    }
    if(argc - 2 != optind) {
        LOG_EXIT("Required: precisely two positional arguments (in bam, out bam).\n");
    }

    return EXIT_SUCCESS;
}
