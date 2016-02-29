#ifndef SHEN_H
#define SHEN_H
#include <math.h>
#include "compiler_util.h"
#include "khash.h"
#include "stdint.h"

#ifndef SHEN_KHASH_INIT
#define SHEN_KHASH_INIT
KHASH_MAP_INIT_INT(shen, uint64_t)
#endif

/*
 * @func shannon_entroy
 * :param: str [char *] Input string
 * :returns: [double] Calculated shannon entropy.
 */
static inline double UNUSED_FUNC(shannon_entropy) (char *str)
{
    khash_t(shen) *h = kh_init(shen);
    double shen = 0.0;
    khiter_t k;
    char *const start = str;
    for(;*str;) {
        k = kh_get(shen, h, *str);
        if(k == kh_end(h)) {
            k = kh_put(shen, h, *str, &ret);
            kh_val(h, k) = 1uL;
        } else ++kh_val(h, k);
        ++str;
    }
    double f;
    const double n_elem = (double)(str - start);
    for(k = 0; k != kh_end(h); ++k) {
        if(!kh_exist(h, k)) continue;
        f = kh_val(h, k) / n_elem;
        shen -= f * log(f);
    }
    kh_destroy(shen, h);
    return shen;
}

/*
 * @func shannon_entropy_acgtn
 * Caculates shannon entropy assuming only A, C, G, T, and N as possible characters.
 * Ignores all ^[acgtACGT]
 * :param: str [char *] Input string
 * :returns: [double] Calculated shannon entropy.
 */
static inline double UNUSED_FUNC(shannon_entropy_acgtn) (char *str)
{
    uint64_t counts[4] = {0uL, 0uL, 0uL, 0uL};
    double shen = 0.0;
    uint64_t i;
    for(i = 0; str[i]; ++i) {
        switch(str[i]) {
            case 'a':
            case 'A': ++counts[0]; break;
            case 'c':
            case 'C': ++counts[1]; break;
            case 'g':
            case 'G': ++counts[2]; break;
            case 't':
            case 'T': ++counts[3]; break;
        }
    }
    double f;
    for(int j = 0; j < 4; ++j) {
        f = (double)counts[j] / i;
        shen -= f * log(f);
    }
    return shen;
}


#endif
