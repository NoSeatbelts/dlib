#include "dlib/bam_util.h"
#include <getopt.h>

int usage() {
    fprintf(stderr, "<-l output_compression_level> in.bam out.bam\n"
                    "Use - for stdin or stdout.\n");
    return EXIT_FAILURE;
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
    // Actually this function. You can't really apply a null function....
    single_aux_check fn = NULL;
    std::function<int (bam1_t *, void *)> fn = NULL;
    // Actually create your type for data and then provide it if needed.
    void *data = NULL;
    dlib::bam_apply_function(argv[optind], argv[optind + 1], fn, data, out_mode);
    return EXIT_SUCCESS;
}
