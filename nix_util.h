#ifndef NIX_UTIL_H
#define NIX_UTIL_H
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
typedef struct rlimit rlimit_t;

void increase_nofile_limit(int soft_limit);
int get_fileno_limit();

#endif /* NIX_UTIL_H */
