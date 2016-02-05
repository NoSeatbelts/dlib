#ifndef LOGGING_UTIL_H
#define LOGGING_UTIL_H

#include <stdarg.h>

#if USE_FPRINTF_MACRO
#	define fprintf(fp, ...) \
    do {\
        fprintf(fp, "[%s] ", __func__);\
        fprintf(fp, ##__VA_ARGS__);\
        } while(0)

#endif

#define LOG_INFO(...) log_info((char *)__func__, ##__VA_ARGS__);
#define LOG_WARNING(...) log_warning((char *)__func__, ##__VA_ARGS__);
#define LOG_ERROR(...) log_error((char *)__func__, __LINE__, ##__VA_ARGS__);
#if !NDEBUG
#	define LOG_DEBUG(...) log_debug((char *)__func__, __LINE__, ##__VA_ARGS__);
#else
#	define LOG_DEBUG(...)
#endif

static inline void log_debug(const char *func, int line, char *fmt, ...) {
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, (char *)"[D:%s:%d] ", func, line);
	vfprintf(stderr, fmt, args);
	va_end(args);
}

static inline void log_warning(const char *func, char *fmt, ...) {
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, (char *)"[W:%s] ", func);
	vfprintf(stderr, fmt, args);
	va_end(args);
}

static inline void log_error(const char *func, int line, char *fmt, ...) {
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, (char *)"[E:%s:%d] ", func, line);
	vfprintf(stderr, fmt, args);
	va_end(args);
	exit(EXIT_FAILURE);
}

static inline void log_info(const char *func, char *fmt, ...) {
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, (char *)"[%s] ", func);
	vfprintf(stderr, fmt, args);
	va_end(args);
}

#define LOG_ASSERT(condition) log_assert(__func__, __LINE__, condition, (const char *)(#condition))

static inline void log_assert(const char *func, int line, int assertion, const char *assert_str) {
	if(assertion) return;
	LOG_ERROR((char *)"Assertion '%s' failed.", assert_str);
}

#endif /* LOGGING_UTIL_H */
