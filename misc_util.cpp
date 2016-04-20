#include "misc_util.h"
#include <time.h>

namespace dlib {

    void string_fmt_time(std::string& ret) {
        time_t now = time(NULL);
        struct tm *t = localtime(&now);
        while(!strftime((char *)ret.data(), ret.capacity(), "%Y:%m:%d %H:%M", t))
            ret.reserve(ret.capacity() ? ret.capacity() << 1: 16);
    }

}
