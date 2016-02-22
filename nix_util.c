#include "nix_util.h"

/*
 * @func increase_nofile_limit
 * Gets the current limit for the number of open files to input argument.
 * Errors on failure.
 * :param: soft_limit [int] New soft limit
 */
void increase_nofile_limit(int soft_limit)
{
	rlimit_t rl;
	getrlimit(RLIMIT_NOFILE, &rl);
	rl.rlim_cur = (soft_limit > rl.rlim_cur) ? soft_limit : rl.rlim_cur;
	if(setrlimit(RLIMIT_NOFILE, &rl)) {
		LOG_EXIT("Could not increase the soft limit for number "
					"of open files in a directory to %i. The hard "
					"limit needs to be changed, which may require sudo privileges."
					"Abort mission!\n", soft_limit);
	}
}

/*
 * @func get_fileno_limit
 * Gets the hard limit for the number of open files on this system.
 * :returns: [int] Hard limit for number of open files.
 */
int get_fileno_limit() {
	rlimit_t rl;
	getrlimit(RLIMIT_NOFILE, &rl);
	return rl.rlim_max;
}
