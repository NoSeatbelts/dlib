#ifndef VCF_UTIL_H
#define VCF_UTIL_H

#include "htslib/vcf.h"

#define check_vcf_open(path, fp, header) \
	do {\
		if((fp = vcf_open(path, "r")) == NULL) {\
			LOG_ERROR("Could not open input [bv]cf %s. Abort!\n", path);\
		}\
		if((header = vcf_hdr_read(fp)) == NULL) {\
			LOG_ERROR("Could not read header from input [bv]cf %s. Abort!\n", path);\
		}\
	} while(0)

#endif /* VCF_UTIL_H */
