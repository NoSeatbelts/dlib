#ifndef MISC_UTIL_H
#define MISC_UTIL_H

// Cribbed from nlopt
#define MAX2(a, b) ((a) > (b) ? (a): (b))

// From Chromium
#define COUNT_OF(x) ((sizeof(x)/sizeof(0[x])) / ((size_t)(!(sizeof(x) % sizeof(0[x])))))

#endif /* MISC_UTIL_H */
