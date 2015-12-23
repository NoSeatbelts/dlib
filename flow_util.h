#ifndef FLOW_UTIL_H
#define FLOW_UTIL_H

#include "stdlib.h"
#include "stdio.h"

void log_abort(char *str) {
	fprintf(stderr, "[E:%s] ", __func__);
	fprintf(stderr, str);
	exit(EXIT_FAILURE);
}

#ifndef ABORT
#define ABORT(message) do {fprintf(stderr, "[E:%s] ", __func__); fprintf(stderr, message); exit(EXIT_FAILURE);} while(0)
#endif
#endif
