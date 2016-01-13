#ifndef LOGGING_UTIL_H
#define LOGGING_UTIL_H

#if USE_FPRINTF_MACRO
#define fprintf(fp, ...) \
    do {\
        fprintf(fp, "[%s] ", __func__);\
        fprintf(fp, ##__VA_ARGS__);\
        } while(0)\

#endif


#define LOG_DEBUG(str, ...) \
    do {\
        fprintf(stderr, "[%s:%d] " str , __func__, __LINE__, ##__VA_ARGS__);\
        } while(0)\

#endif /* LOGGING_UTIL_H */
