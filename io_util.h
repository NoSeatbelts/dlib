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


#ifdef __cplusplus
extern "C" {
#endif
// FUNCTIONS
// Function Declarations
int isfile(char *fname); // True if file exists.
int bash_system(const char *command); // Call command with bash instead of sh
#ifdef __cplusplus
int my_system (const char *command, const char *executable="/bin/bash");
#else
int my_system (const char *command, const char *executable);
#endif

FILE *open_ofp(char *infname);
FILE *open_ifp(char *infname);
int file_has_ext(char *fn, const char *ext);
int is_bgzipped_vcf(char *fn);

// Inline Function Definitions
static inline gzFile open_gzfile(char *infname) {
	if(strcmp(infname, "-") == 0 || strcmp(infname, "stdin") == 0) {
		LOG_DEBUG("Reading from standard in because infname is %s.\n", infname);
		return gzdopen(STDIN_FILENO, "r"); // Opens stdin.
	} else {
		LOG_DEBUG("Reading from %s.\n", infname);
		return gzopen(infname, "r");
	}
}

// Function Macros
/* CHECK_POPEN
 * Executes cmd with popen and exits if it returns a non-zero exit status.
 * cmd [arg/char *] Command to execute via popen
 */
#define CHECK_POPEN(cmd) \
	do {\
		LOG_DEBUG("About to call '%s' via popen.\n", cmd);\
		if(pclose(popen(cmd, "w"))) {\
			LOG_ERROR("Command '%s' failed. Abort!\n", cmd);\
		}\
	} while(0)

/* CHECK_CALL
 * Executes cmd with system and exits if it returns a non-zero exit status.
 * cmd [arg/char *] Command to execute via popen
 */
#ifndef CHECK_CALL
#	if !NDEBUG
#		define CHECK_CALL(buff)\
	do {\
		LOG_DEBUG("Now check calling command '%s'.\n", buff); \
		if(system(buff) < 0) {\
			LOG_ERROR("System call failed. Command: '%s'.\n", buff);\
		}\
	} while(0)

#	else
#		define CHECK_CALL(buff) \
	if(system(buff) < 0) {\
		LOG_ERROR("System call failed. Command: '%s'.\n", buff);\
	}

#	endif
#endif

static int count_lines(char *fname) {
	int ret = 0;
	FILE *fp = fopen(fname, "r");
	if(!fp) {
		fprintf(stderr, "[E:%s] Could not open file %s. Abort mission!\n", __func__, fname);
		exit(EXIT_FAILURE);
	}
	start:
	switch(getc(fp)) {
		case EOF: fclose(fp); return ret;
		case '\n': ++ret;
	}
	goto start;
}
#ifdef __cplusplus
} /* extern "C" */
#endif


#endif
