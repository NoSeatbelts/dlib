#ifndef MISC_UTIL_H
#define MISC_UTIL_H
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "compiler_util.h"
#include "logging_util.h"
#include "htslib/kstring.h"

#ifdef __cplusplus
#    include <string>
namespace dlib {
    void string_fmt_time(std::string& ret);
}
#endif

// Cribbed from nlopt
#define MAX2(a, b) ((a) > (b) ? (a): (b))

// From Chromium
#define COUNT_OF(x) ((sizeof(x)/sizeof(0[x])) / ((size_t)(!(sizeof(x) % sizeof(0[x])))))


#endif /* MISC_UTIL_H */
