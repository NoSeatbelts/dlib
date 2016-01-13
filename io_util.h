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
int isfile(char *fname); // True if file exists.
int bash_system(const char *command); // Call command with bash instead of sh
#ifdef __cplusplus
int my_system (const char *command, const char *executable="/bin/bash");
#else
int my_system (const char *command, const char *executable);
#endif

inline gzFile open_gzfile(char *infname); // Opens gzFile from file or stdin ('-', 'stdin')

// Inline Function Definitions
inline gzFile open_gzfile(char *infname) {
	if(strcmp(infname, "-") == 0 || strcmp(infname, "stdin") == 0) {
		LOG_DEBUG("Reading from standard in because infname is %s.\n", infname);
		return gzdopen(STDIN_FILENO, "r"); // Opens stdin.
	} else {
		LOG_DEBUG("Reading from %s.\n", infname);
		return gzopen(infname, "r");
	}
}

static inline FILE *open_ofp(char *infname) {
	if(strcmp(infname, "-") == 0 || strcmp(infname, "stdout") == 0) {
        LOG_DEBUG("Writing to standard out because infname is %s.\n", infname);
		return stdout; // Opens stdin.
	} else {
        LOG_DEBUG("Writing to %s.\n", infname);
		return fopen(infname, "w");
	}
}

// Function Macros
/* CHECK_POPEN
 * Executes cmd with popen and exits if it returns a non-zero exit status.
 * cmd [arg/char *] Command to execute via popen
 */
#define CHECK_POPEN(cmd) \
	do {\
		if(pclose(popen(cmd, "w"))) {\
			LOG_ERROR("Command '%s' failed. Abort!\n", __func__, cmd);\
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





#endif
