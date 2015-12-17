#include "nix_util.h"

void increase_nofile_limit(int soft_limit)
{
	rlimit_t rl;
	getrlimit(RLIMIT_NOFILE, &rl);
	rl.rlim_cur = (soft_limit > rl.rlim_cur) ? soft_limit : rl.rlim_cur;
	if(setrlimit(RLIMIT_NOFILE, &rl)) {
		fprintf(stderr, "Could not increase the soft limit for number "
						"of open files in a directory to %i. The hard "
						"limit needs to be changed, which may require sudo privileges."
						"Abort mission!\n", soft_limit);
		exit(EXIT_FAILURE);
	}
	return;
}

int get_fileno_limit() {
	rlimit_t rl;
	getrlimit(RLIMIT_NOFILE, &rl);
	return rl.rlim_max;
}
