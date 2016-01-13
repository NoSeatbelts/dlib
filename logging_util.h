#ifndef LOGGING_UTIL_H
#define LOGGING_UTIL_H

#if USE_FPRINTF_MACRO
#define fprintf(fp, ...) \
    do {\
        fprintf(fp, "[%s] ", __func__);\
        fprintf(fp, ##__VA_ARGS__);\
        } while(0)\

#endif


#if !NDEBUG
#define LOG_DEBUG(str, ...)\
    do {\
        fprintf(stderr, "[D:%s:%d] " str, __func__, __LINE__, ##__VA_ARGS__);\
    } while(0)
#else
#define LOG_DEBUG(str, ...)
#endif

#define LOG_INFO(str, ...) \
    do {\
        fprintf(stderr, "[%s] " str, __func__, ##__VA_ARGS__);\
    } while(0)

#define LOG_ERROR(str, ...) \
    do {\
        fprintf(stderr, "[E:%s:%d] " str, __func__, __LINE__, ##__VA_ARGS__);\
        exit(EXIT_FAILURE);\
    } while(0)

#endif /* LOGGING_UTIL_H */
