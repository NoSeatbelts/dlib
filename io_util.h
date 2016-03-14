#ifndef O_IO_UTIL_H
#define O_IO_UTIL_H
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <zlib.h>
#include "logging_util.h"


// FUNCTIONS
// Function Declarations
#ifdef __cplusplus
namespace dlib {
    int my_system (const char *command, const char *executable="/bin/bash");
#else
    int my_system (const char *command, const char *executable);
#endif
    int isfile(char *fname); // True if file exists.
    int bash_system(const char *command); // Call command with bash instead of sh
    FILE *open_ofp(const char *outfname);
    FILE *open_ifp(const char *infname);
    int file_has_ext(char *fn, const char *ext);
    int is_bgzipped_vcf(char *fn);
    void check_popen(const char *cmd);
    void check_call(const char *cmd);
    gzFile open_gzfile(char *infname);
    int count_lines(const char *fname);
#ifdef __cplusplus
}
#endif

#endif
