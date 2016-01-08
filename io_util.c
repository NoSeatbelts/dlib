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
		fprintf(stderr, "[W:%s] called process '%s' failed. Check return status.\n", __func__, command);
	}
	return status;
}

int bash_system (const char *command)
{
  return my_system(command, "/bin/bash");
}
