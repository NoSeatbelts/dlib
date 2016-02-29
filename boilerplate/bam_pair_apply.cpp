#include "dlib/bam_util.h"
#include <getopt.h>

int usage() {
    fprintf(stderr, "<-l output_compression_level> in.namesrt.bam out.bam\n"
                    "Use - for stdin or stdout.\n");
    return EXIT_FAILURE;
}

class opts {
public:
	int level;
	opts(): level(-1){}
};

int main(int argc, char *argv[]) {
    if(argc < 3) {
        return usage();
    }
    if(strcmp(argv[1], "--help") == 0) {
        usage();
        return EXIT_SUCCESS;
    }
    int c;
    opts opts();
    opts.level = -1;
    char out_mode[4] = "wb";
    // Add more options here to take parameters.
    while((c = get_opt(argc, argv, "l:h?")) > -1) {
        switch(c) {
        case 'l':
            opts.level = atoi(optarg) % 10; out_mode[2] = opts.level + '0'; break;
        case 'h': case '?': usage(); return EXIT_SUCCESS;
        default: LOG_INFO("Unrecognized option %c\n", optopt); return usage();
        }
    }
    if(argc - 2 != optind) {
        LOG_EXIT("Required: precisely two positional arguments (in bam, out bam).\n");
    }
    // Actually this function. You can't really apply a null function....
    pair_aux_fn fn = NULL;
    // Actually create your type for data and then provide it if needed.
    dlib::bam_pair_apply_function(argv[optind], argv[optind + 1], fn, (void *)opts, out_mode);
    return EXIT_SUCCESS;
}
