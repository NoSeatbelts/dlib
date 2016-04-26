#ifndef CSTR_UTIL_H
#define CSTR_UTIL_H

#include <zlib.h>
#include <time.h>
#include "logging_util.h"
#include "char_util.h"
#include "kseq.h"

#ifndef KSEQ_DEC_GZ
#define KSEQ_DEC_GZ
KSEQ_INIT(gzFile, gzread)
#endif

#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 30
#endif

#ifndef SEQBUF_SIZE
#define SEQBUF_SIZE 300
#endif

#define kputsnl(literal, ks) kputsn(literal, sizeof(literal) - 1, ks)

#ifdef __cplusplus
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

namespace dlib {

/* sprintf's to the buffer in string.
 * If not long enough, everything breaks. Be careful!
*/
//#define stringprintf(str, ...) str.resize(sprintf(const_cast<char *>(str.data()), ##__VA_ARGS__))

    std::vector<std::string> tokenize(const char *str, char c='\t');

    extern "C" char *rand_string(char *str, size_t size);

    static inline int strhd_thresh(std::string const& str1, std::string const& str2, size_t mmlim) {
        return strhd_thresh(str1.c_str(), str2.c_str(), mmlim);
    }

    static inline int strhd_thresh(const char *str1, const char *str2, size_t mmlim) {
        while(*str1) {
            if(*str1 != *str2)
                if(mmlim-- == 0)
                    return 0;
            ++str1, ++str2;
        }
        return 1;
    }

    /*
     * Fast, signed atoi.
     */
    CONST static inline int fast_delim_atoi(char *str, char delim=',')
    {
        int ret = 0;
        int sign = 1;
        switch(*str) {
            // While retuurn the wrong value if *str is ',': don't start a number with ,!
            case '-': sign = -1; break;
            case '+': break;
            default: ret = *str - '0';
        }
        ++str;
        while(*str) {
            if(*str != delim)
                ret = ret*10 + (*str++ - '0');
            else ++str; // Ignore it
        }
        return ret * sign;
    }


#endif

    /*
     * @func kfill_rc
     * :param: seq [kseq_t *] Input kseq object
     * :param: buffer [char *] String to copy into.
     * Fills a buffer with reverse-complemented fastq lines
     * (seq line through the end of the quality line,
     *  with a terminal newline and a null character)
     */
    static inline void kfill_rc(kseq_t *seq, char *buffer) {
        uint64_t i = seq->seq.l;
        char *seqp = seq->seq.s + i;
        // Add seq field.
        for(; i ; --i) *buffer++ = nuc_cmpl(*--seqp);
        // Add "\n+\n"
        *buffer++ = '\n'; *buffer++ = '+'; *buffer++ = '\n';
        // Add reversed quality
        i = seq->qual.l;
        seqp = seq->qual.s + i;
        for(; i ; --i) *buffer++ = *--seqp;
        // Terminate with newline and null character.
        *buffer++ = '\n'; *buffer++ = '\0';
    }

    /*
     * @func fill_rc
     * :param: str [char *] Input string
     * :param: buffer [char *] String to copy into.
     * :param: len [size_t] Length of input string.
     * Fills a buffer with reverse-complemented characters.
     */
    static inline void fill_rc(char *str, char *buffer, size_t len) {
        str += len; // Skip to the end of the string.
        for(; len; --len)
            *buffer++ = nuc_cmpl(*--str);
        *buffer++ = '\0';
    }

    /*
     * @func fill_rv
     * :param: str [char *] Input string
     * :param: buffer [char *] String to copy into.
     * :param: len [size_t] Length of input string.
     * Fills a buffer with reversed characters.
     */
    static inline void fill_rv(char *str, char *buffer, size_t len) {
        str += len; // Skip to the end of the string.
        for(; len; --len)
            *buffer++ = *--str;
        *buffer++ = '\0';
    }

    /*
     * Returns a null-terminated string with the extension and terminal period removed.
     * Warning: Must be freed!
     */
    char *trim_ext(char *fname);

    /*
     * Fast positive atoi
     */
    CONST static inline int fp_atoi(char *str)
    {
        int ret = *str++ - '0';
        while(*str)
            ret = ret*10 + (*str++ - '0');
        return ret;
    }

    /*
     * Fast, signed atoi.
     */
    CONST static inline int fast_atoi(char *str)
    {
        int ret = 0;
        int sign = 1;
        switch(*str) {
            case '-': sign = -1; break;
            case '+': break;
            default: ret = *str - '0';
        }
        ++str;
        while(*str)
            ret = ret*10 + (*str++ - '0');
        return ret * sign;
    }

    static inline char *revcmp(char *dest, char *src, uint64_t l)
    {
        src += l;
        for(; l; --l)
            *dest++ = nuc_cmpl(*--src);
        *dest++ = '\0';
        return dest;
    }


    CONST static inline int lex_memcmp(char *s1, char *s2, size_t l)
    {
        for(; l; --l) {
            if(*s1 != *s2) return *s1 < *s2;
            ++s1, ++s2;
        }
        return -1;
    }

    CONST static inline int lex_strlt(char *s1, char *s2)
    {
        while(*s1) {
            if(*s1 != *s2) return *s2 < *s1;
            ++s1, ++s2;
        }
        return -1;
    }

    CONST static inline int lex_lt(char *s, size_t l)
    {
        char *s2 = s + l - 1;
        while(*s) {
            if(*s != *s2) return *s < *s2;
            ++s; --s2;
        }
        return -1; // Palindromic
    }


    /*
     * Returns a null-terminated string with the default outfname.
     * Warning: Must be freed!
     */

    char *make_default_outfname(char *fname, const char *suffix);

    static inline char *kstrdup(kstring_t *ks)
    {
        char *ret = (char *)malloc((ks->l + 1) * sizeof(char));
        memcpy(ret, ks->s, ks->l + 1);
        return ret;
    }

#ifdef __cplusplus
}
#endif

#endif
