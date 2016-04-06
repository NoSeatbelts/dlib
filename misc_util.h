#ifndef MISC_UTIL_H
#define MISC_UTIL_H
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "logging_util.h"
#include "compiler_util.h"
#include "logging_util.h"
#include "htslib/kstring.h"

#ifdef __cplusplus
#    include <string>
namespace dlib {
    void string_fmt_time(std::string& ret);
}
#endif

#if !NDEBUG
#   define cond_free(var)\
    do {\
        if(var) {\
            /*LOG_DEBUG("About to free variable at %p (%s).\n", (void *)var, #var); */\
            free(var);\
            var = NULL;\
        }\
    } while(0)
#else
#   define cond_free(var) do {if(var) {free(var); var = NULL;}} while(0)
#endif

// Cribbed from nlopt
#define MAX2(a, b) ((a) > (b) ? (a): (b))

// From Chromium
#define COUNT_OF(x) ((sizeof(x)/sizeof(0[x])) / ((size_t)(!(sizeof(x) % sizeof(0[x])))))

#ifdef LOG_MEMCPY
#define memcpy(dest, src, size)\
    do {\
        LOG_DEBUG("Copying to %p from %p a number of bytes %lu.\n", (void *)(dest), (void *)(src), size);\
        memcpy(dest, src, size);\
    } while(0)
#endif



#endif /* MISC_UTIL_H */
