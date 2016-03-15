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
}

#endif
