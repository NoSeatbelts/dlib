#include "cstr_util.h"

#ifdef __cplusplus
namespace dlib {
    std::vector<std::string> tokenize(const char *str, char c) {
        std::vector<std::string> ret;
        do {
            const char *begin = str;
            while(*str && *str != c)
                str++;
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
        srand(time(NULL)); // Pick a seed!
        const char charset[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKSTFUOMGZWTF";
        if (size) {
            --size;
            for (size_t n = 0; n < size; n++) {
                str[n] = charset[(int)(rand() % (int) (sizeof charset - 1))];
            }
            str[size] = '\0';
        }
        return str;
    }
    size_t safe_stringprintf(std::string& str, const char *fmt, ...) {
        va_list ap;
        size_t l;
        va_start(ap, fmt);
        if((l = snprintf(const_cast<char *>(str.data()), str.size(), fmt, ap)) > str.size()) {
            str.resize(l);
            snprintf(const_cast<char *>(str.data()), str.size(), fmt, ap);
        }
        va_end(ap);
        return l;
    }
}

#endif
