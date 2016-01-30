#include "io_util.h"

int isfile(char *fname)
{
	return access(fname, F_OK) != -1;
}


int my_system (const char *command, const char *executable)
{
	int status;
	pid_t pid = fork ();
	if (pid == 0) {
		/* This is the child process.  Execute the shell command. */
		execl(executable, executable, "-c", command, NULL);
		_exit(EXIT_FAILURE);
	} else if(pid < 0) status = -1;
	/* The fork failed.  Report failure.  */
	else if (waitpid (pid, &status, 0) != pid) {
	/* This is the parent process.  Wait for the child to complete.  */
		status = -1;
		LOG_WARNING("Called process '%s' failed. Check return status (%i).\n", command, status);
	}
	return status;
}

int bash_system (const char *command)
{
  return my_system(command, "/bin/bash");
}

FILE *open_ifp(char *path) {
	return (!path || !*path || path[0] == '-') ? stdin: fopen(path, "r");
}

FILE *open_ofp(char *path) {
	if(!path || !*path || strcmp(path, "-") == 0 || strcmp(path, "stdout") == 0) {
		LOG_DEBUG("Writing to standard out because path is %s.\n", path);
		return stdout; // Opens stdin.
	}
	LOG_DEBUG("Writing to %s.\n", path);
	return fopen(path, "w");
}

int is_bgzipped_vcf(char *fn)
{
	return strcmp(strrchr(fn, '.') - 4, ".vcf.gz") == 0;
}

int file_has_ext(char *fn, const char *ext)
{
	return strcmp(strrchr(fn, '.') + 1, ext) == 0;
}
