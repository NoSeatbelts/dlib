#ifndef CSTR_UTIL_H
#define CSTR_UTIL_H
#include <inttypes.h>
#include <time.h>
#include <zlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <math.h>
#include "kseq.h"
#include "kstring.h"
#include "compiler_util.h"
#include "logging_util.h"
#include "char_util.h"

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

namespace dlib {

	/*
	 * @func rand_string
	 * Stolen from stackoverflow.
	 * Fills the str with size random characters from the charset string
	 * and appends a terminal null character.
	 */
	static inline char *rand_string(char *str, size_t size)
	{
	    srand(time(NULL)); // Pick a seed!
	    const char charset[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKSTFUOMGZWTF";
	    if (size) {
	        --size;
	        for (size_t n = 0; n < size; n++) {
	            str[n] = charset[(int)(rand() % (int) (sizeof charset - 1))];
	        }
	        str[size] = '\0';
	    }
	    return str;
	}


	/*
	 * @func fill_csv_buffer
	 * Used to write out uint32_t arrays as textual bam tags.
	 * :param: readlen [int] Length of read
	 * :param: arr [uint32_t *] Array of values to put into the buffer.
	 * :param: buffer [char *] Buffer for the values.
	 * :param: prefix_typecode [const char *] typecode and prefix to append for a bam tag.
	 */
	static void fill_csv_buffer(int readlen, uint32_t *arr, char *buffer, const char *prefix_typecode)
	{
	    kstring_t ks = {0, 0, NULL};
	    kputs(prefix_typecode, &ks);
	    for(int i = 0; i < readlen; i++)
	        ksprintf(&ks, ",%u", arr[i]);
	    strcpy(buffer, ks.s);
	    free(ks.s);
	}

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
	static inline char *trim_ext(char *fname)
	{
	    LOG_DEBUG("Now trimming char * %s.\n", fname);
	    char *ret = (char *)malloc((strlen(fname) + 1) * sizeof(char ));
	    char *found_pos = strrchr(fname, '.');
	    if(!found_pos) {
	        LOG_EXIT("Could not trim file name's extension. Looks like it's missing a '.' (name: '%s').\n", fname);
	    }
	    memcpy(ret, fname, (found_pos - fname) * sizeof(char));
	    ret[found_pos - fname] = '\0';
	    return ret;
	}

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
	static inline char *make_default_outfname(char *fname, const char *suffix)
	{
	    char buf[200];
	    char *prefix = trim_ext(fname);
	    strcpy(buf, prefix);
	    strcat(buf, suffix);
	    char *ret = strdup(buf);
	    free(prefix);
	    return ret;
	}

	static inline char *kstrdup(kstring_t *ks)
	{
	    char *ret = (char *)malloc((ks->l + 1) * sizeof(char));
	    memcpy(ret, ks->s, ks->l + 1);
	    return ret;
	}

}
#endif
