#ifndef COMPILER_UTIL_H
#define COMPILER_UTIL_H

#define __CSTD_LIMIT_MACROS

#ifdef __GNUC__
#    define LIKELY(x) __builtin_expect((x),1)
#    define UNLIKELY(x) __builtin_expect((x),0)
#    define PURE __attribute__((pure))
#else
#    define LIKELY(x) (x)
#    define UNLIKELY(x) (x)
#    define PURE
#endif /* #ifdef __GNUC__ */

#ifdef __GNUC__
#    ifndef UNUSED
#        define UNUSED(x) x __attribute__((unused))
#    endif
#    ifndef UNUSED_FUNC
#        define UNUSED_FUNC(x) __attribute__((__unused__)) x
#    endif
#    define ALWAYS_INLINE __attribute__((always_inline)) inline
#    define CONST __attribute__((const))
#else
#    define CONST
#    define UNUSED(x) x
#    define UNUSED_FUNC(x) x
#endif /* #ifdef __GNUC__ */

#ifndef FOREVER
#    define FOREVER for(;;)
#endif

#endif /* COMPILER_UTIL_H */
