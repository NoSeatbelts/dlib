#ifndef NIX_UTIL_H
#define NIX_UTIL_H
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include "logging_util.h"

typedef struct rlimit rlimit_t;

#ifdef __cplusplus
extern "C" {
#endif

void increase_nofile_limit(int soft_limit);
int get_fileno_limit();

#ifdef __cplusplus
}
#endif

#endif /* NIX_UTIL_H */
