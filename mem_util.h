#ifndef MEM_UTIL_H
#define MEM_UTIL_H

#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include "logging_util.h"
#include "compiler_util.h"

// Leave semicolon out so that it looks like a normal function.
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

#endif
