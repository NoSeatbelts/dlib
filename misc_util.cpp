#include "misc_util.h"

void string_fmt_time(std::string& ret) {
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    while(strftime((char *)ret.data(), ret.capacity(), "%Y:%m:%d %H:%M", t) == 0)
        ret.reserve(ret.size() << 1);
    return;
}
