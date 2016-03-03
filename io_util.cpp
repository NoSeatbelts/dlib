#include "io_util.h"

namespace dlib {

    int isfile(char *fname)
    {
        return access(fname, F_OK) != -1;
    }

    /* check_popen
     * Executes cmd with popen and exits if it returns a non-zero exit status.
     * cmd [arg/char *] Command to execute via popen
     */
    void check_popen(const char *cmd) {
        LOG_DEBUG("Now check popen-ing command '%s'.\n", cmd);
        int ret;
        if((ret = pclose(popen(cmd, "w"))) != 0) {
            LOG_EXIT("System call failed. Command: '%s'. Return code: %i.\n", cmd, ret);
        }
    }

    /* check_call
     * Executes cmd with system and exits if it returns a non-zero exit status.
     * cmd [arg/char *] Command to execute via popen
     */
    void check_call(const char *cmd) {
        LOG_DEBUG("Now check calling command '%s'.\n", cmd);
        int ret;
        if((ret = system(cmd)) != 0) {
            LOG_EXIT("System call failed. Command: '%s'. Return code: %i.\n", cmd, ret);
        }
    }


    gzFile open_gzfile(char *infname) {
        if(strcmp(infname, "-") == 0 || strcmp(infname, "stdin") == 0) {
            LOG_DEBUG("Reading from standard in because infname is %s.\n", infname);
            return gzdopen(STDIN_FILENO, "r"); // Opens stdin.
        } else {
            LOG_DEBUG("Reading from %s.\n", infname);
            return gzopen(infname, "r");
        }
    }


    int count_lines(char *fname) {
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
        char *tmp = strrchr(fn, '.') - 4;
        return tmp < fn ? 0: strcmp(tmp, ".vcf.gz") == 0;
    }

    int file_has_ext(char *fn, const char *ext)
    {
        return strcmp(strrchr(fn, '.') + 1, ext) == 0;
    }

}
