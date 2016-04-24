#include "cstr_util.h"

#ifdef __cplusplus
namespace dlib {
    std::vector<std::string> tokenize(const char *str, char c) {
        std::vector<std::string> ret;
        do {
            const char *begin = str;
            while(*str && *str != c) ++str;
            ret.emplace_back(begin, str);
        } while(*str++);
        return ret;
    }

    /*
     * @func rand_string
     * Stolen from stackoverflow.
     * Fills the str with size random characters from the charset string
     * and appends a terminal null character.
     */
    char *rand_string(char *str, size_t size)
    {
        if(!size) {
            *str = '\0';
            return str;
        }
        srand(time(NULL)); // Pick a seed!
        const char charset[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKSTFUOMGZWTF";
        str[size--] = '\0';
        while(size--)
            str[size] = charset[(rand() % (sizeof charset - 1))];
        return str;
    }

    char *trim_ext(char *fname)
    {
        char *ret, *found_pos;
        LOG_DEBUG("Now trimming char * %s.\n", fname);
        ret = (char *)malloc((strlen(fname) + 1));

        if((found_pos = strrchr(fname, '.')) == nullptr) LOG_EXIT("Filename '%s' has no '.'.).\n", fname);
        memcpy(ret, fname, (found_pos - fname));
        ret[found_pos - fname] = '\0';
        return ret;
    }

    char *make_default_outfname(char *fname, const char *suffix)
    {
        char buf[200];
        char *prefix = trim_ext(fname);
        strcpy(buf, prefix);
        strcat(buf, suffix);
        char *ret = strdup(buf);
        free(prefix);
        return ret;
    }
}

#endif /* ifdef __cplusplus */
