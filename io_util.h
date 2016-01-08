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

#define CHECK_POPEN(cmd) \
	do {\
		FILE *tmp_popen_##cmd = popen(cmd, "w");\
		if(pclose(tmp_popen_##cmd)) {\
			fprintf(stderr, "[E:%s] Command '%s' failed. Abort!\n", __func__, cmd);\
			exit(EXIT_FAILURE);\
		}\
	} while(0)

#ifndef CHECK_CALL
#if !NDEBUG
#define CHECK_CALL(buff) \
	fprintf(stderr, "[D:%s] Now check calling command '%s'.\n", __func__, buff); \
	if(system(buff) < 0) fprintf(stderr, "[D:%s] System call failed. Command: '%s'.\n", __func__, buff)
#else
#define CHECK_CALL(buff) \
	if(system(buff) < 0) fprintf(stderr, "[D:%s] System call failed. Command: '%s'.\n", __func__, buff)
#endif
#endif

int isfile(char *fname);

inline gzFile open_gzfile(char *infname) {
	if(strcmp(infname, "-") == 0 || strcmp(infname, "stdin") == 0) {
		fprintf(stderr, "Reading from standard in because infname is %s.\n", infname);
		return gzdopen(STDIN_FILENO, "r"); // Opens stdin.
	}
	else {
		fprintf(stderr, "Reading from %s.\n", infname);
		return gzopen(infname, "r");
	}
}

inline FILE *open_ofp(char *infname) {
	if(strcmp(infname, "-") == 0 || strcmp(infname, "stdout") == 0) {
		fprintf(stderr, "Reading from standard in because infname is %s.\n", infname);
		return stdout; // Opens stdin.
	}
	else {
		fprintf(stderr, "Reading from %s.\n", infname);
		return fopen(infname, "r");
	}
}

int bash_system (const char *command);

#endif
