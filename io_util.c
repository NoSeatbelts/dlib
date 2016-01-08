#include "io_util.h"

int isfile(char *fname)
{
	return access(fname, F_OK) != -1;
}

/* Execute the command using this shell program.  */
#define SHELL "/bin/bash"

int bash_system (const char *command)
{
  int status;
  pid_t pid;
  pid = fork ();
  if (pid == 0)
    {
      /* This is the child process.  Execute the shell command. */
      execl (SHELL, SHELL, "-c", command, NULL);
      _exit (EXIT_FAILURE);
    }
  else if (pid < 0)
    /* The fork failed.  Report failure.  */
    status = -1;
  else
    /* This is the parent process.  Wait for the child to complete.  */
    if (waitpid (pid, &status, 0) != pid)
      status = -1;
  return status;
}
