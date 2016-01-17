#ifndef KMER_UTIL_H
#include "logging_util.h"
// Largest odd kmer that can be held in a uint64_t
#define MAX_KMER 31
#define DEFAULT_KMER 21

#ifdef __cplusplus
extern "C" {
#endif

	inline void kmer2cstr(uint64_t kmer, int k, char *buf)
	{
		buf += k;
		while(k) *--buf = ((kmer >> (2 * --k)) & 0x3);
		LOG_DEBUG("kmer %lu has now become string '%s'.\n", kmer, buf);
	}

#ifdef __cplusplus
}
#endif

#endif
